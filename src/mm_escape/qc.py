"""
Stage 04 — per-cell quality control and doublet detection.

WHAT THIS STAGE DECIDES
-----------------------
Which cells enter the analysis. That matters more here than in a typical scRNA
project because the headline metric, `frac_double_negative`, is a *fraction of
zeros*: any filter that preferentially keeps or drops shallow cells moves it
directly. So every threshold in this module is derived from this cohort's own
distributions and written down, and nothing is inherited.

THREE THINGS THAT ARE NOT DEFAULTS
----------------------------------
1. **Thresholds are derived per cohort, never pooled.** The four cohorts differ by
   ~1.9x in genes detected per cell (MMRF 1916, WU2 1210, Donor 1103, WU1 1023 —
   sample-level medians, pre-QC). A pooled MAD would flag much of WashU cohort 1 as
   low-quality for a reason that is batch, not biology, and WashU cohort 1 is 23 of
   the 54 myeloma samples. See `mad_thresholds`.

2. **`sc-best-practices`' numbers are not copied.** Its 5-MAD counts and 8%
   mitochondrial cap are defaults derived from healthy PBMC/BMMC. This is myeloma
   marrow, and malignant plasma cells are professional secretors with unusual
   library composition. The MAD *procedure* is adopted; the *values* are recomputed
   and reported by `qc_report`.

3. **The deposit is already filtered, differently in each cohort — and that is
   measured here rather than inherited.** An earlier reading of GSE223060's stated
   Seurat filter concluded it "was not applied to what is deposited", reasoning from
   a cohort-wide average UMI count. That average pooled MMRF with WashU and hid a
   per-cohort truth. What the files show:

       WU1, WU2   UMI < 10,000   UMI >= 1,000   pct_mt < 20%   genes >= 200
       MMRF       uncensored     UMI >= 1,000   pct_mt < 10%   genes >= 200
       Donor      uncensored     uncensored     pct_mt < 20%   genes >= 200

   The 10,000-UMI ceiling on the two WashU cohorts is the consequential one. Malignant
   plasma cells are the highest-RNA-content cells in marrow, so that ceiling did not
   remove a random slice — measured in the uncensored cohorts, the band above it is
   enriched 3-21x for `TNFRSF17` and **20-70x for `GPRC5D`**. 36 of the 54 myeloma
   samples therefore had the antigen-positive tail of their own tumours removed before
   deposit, which inflates `frac_double_negative` for those cohorts.

   Nothing is done about it *here*: truncating MMRF and Donor to match would discard
   42% of MMRF's cells to make every cohort equally damaged. It is carried forward as
   a quantified confounder, and stage 08 runs the truncate-everything-at-10k version
   as a sensitivity analysis. See `notebooks/04_qc.ipynb` and
   `results/04_qc/umi_censoring_effect.csv`.

`pct_counts_in_top_20_genes` IS TRACKED BUT DOES NOT FILTER
-----------------------------------------------------------
It is tracked because it is one of the few handles on ambient RNA available at all
here — SoupX/DecontX need unfiltered matrices, which this deposit does not have.

It does **not** delete cells, and that is a deliberate, measured departure from the
generic recipe rather than an omission. In most tissues a library dominated by a
handful of transcripts means an empty-ish droplet full of soup. In myeloma marrow it
also means a plasma cell: a professional secretor whose normal state is to have its
library dominated by immunoglobulin. On this cohort the high-`top20` tail is 21x
enriched for `TNFRSF17` (BCMA), so filtering on it would preferentially delete
antigen-POSITIVE malignant cells and inflate `frac_double_negative` — biased in the
project's own direction of interest. See `DEFAULT_FILTERS` for the numbers.

DIRECTION OF THE REMAINING BIAS
-------------------------------
QC cannot remove the two errors that matter downstream, and it is worth being explicit
that it *shifts* rather than resolves them: dropping shallow cells reduces dropout
(which inflates the escape fraction) while leaving ambient contamination (which
deflates it) untouched. Stage 08 bounds both; this stage only has to avoid making
either worse in a cohort-correlated way, which is what per-cohort thresholds buy.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from . import config, io

__all__ = [
    "add_qc_metrics",
    "mad_outlier",
    "mad_thresholds",
    "flag_outliers",
    "DEFAULT_FILTERS",
    "ALL_FLAGS",
    "MAD_METRICS",
    "detect_doublets",
    "run_sample_qc",
    "run_cohort_qc",
    "collect_obs",
    "cohort_thresholds",
    "qc_report",
    "load_checkpoints",
    "DoubletDetectionError",
]


class DoubletDetectionError(RuntimeError):
    """Raised when the scDblFinder bridge cannot run.

    Never silently downgraded to "no doublets" — a sample whose doublets were never
    looked for must not be indistinguishable from one where none were found.
    """


# ---------------------------------------------------------------------------
# QC metrics
# ---------------------------------------------------------------------------

def add_qc_metrics(adata: AnnData, *, inplace: bool = True) -> AnnData:
    """Attach the QC metrics this stage filters on.

    Adds scanpy's standard set plus `pct_counts_mt`, `pct_counts_ribo`,
    `pct_counts_hb` and `pct_counts_in_top_20_genes`, and the two log1p columns the
    MAD rule is applied to.

    Gene classes are matched on the **deposited symbols**, because this runs before
    stage 05's gene-space harmonization — the object still carries whichever HGNC
    vintage its Cell Ranger reference used. The prefixes below are stable across both
    builds (`MT-`, `RPS`/`RPL`, `HB[ABDEGQZ]`), which is why prefix matching is safe
    here where a curated symbol list would not be.
    """
    if not inplace:
        adata = adata.copy()

    names = adata.var_names.str.upper()
    adata.var["mt"] = names.str.startswith("MT-")
    adata.var["ribo"] = names.str.startswith(("RPS", "RPL"))
    # HBB/HBA1/HBA2/... but not HBEGF or HBP1 — the trailing character class is what
    # keeps this from silently annotating unrelated genes as haemoglobin.
    adata.var["hb"] = names.str.contains(r"^HB[ABDEGQZ]\d?$", regex=True)

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        percent_top=[20],
        log1p=True,
        inplace=True,
    )
    return adata


# ---------------------------------------------------------------------------
# MAD-based outlier calling
# ---------------------------------------------------------------------------

def mad_outlier(values: Sequence[float] | np.ndarray, n_mads: float = 5.0) -> np.ndarray:
    """Boolean mask of points more than `n_mads` MADs from the median.

    MAD is the raw median absolute deviation, NOT scaled by 1.4826. That is the
    convention `sc-best-practices` uses and the one the `n_mads` counts here are
    calibrated to; mixing the two silently changes every threshold by a third.

    A zero MAD (a constant or near-constant metric) would make every non-median point
    infinitely deviant, so it returns all-False instead — a degenerate metric should
    filter nothing, not everything.
    """
    values = np.asarray(values, dtype=float)
    median = np.median(values)
    mad = np.median(np.abs(values - median))
    if not np.isfinite(mad) or mad == 0:
        return np.zeros(values.shape, dtype=bool)
    return np.abs(values - median) > n_mads * mad


def mad_thresholds(
    frame: pd.DataFrame,
    metrics: Iterable[str],
    *,
    n_mads: float = 5.0,
) -> pd.DataFrame:
    """The numeric interval each metric's MAD rule implies, as a reportable table.

    `flag_outliers` could compute the masks without materialising this, but then the
    thresholds this cohort actually produces would exist only inside a boolean array.
    CLAUDE.md requires them written down, so they are returned as data: one row per
    metric with `median`, `mad`, `lower`, `upper` and how many cells each side drops.
    """
    rows = []
    for metric in metrics:
        values = frame[metric].to_numpy(dtype=float)
        median = float(np.median(values))
        mad = float(np.median(np.abs(values - median)))
        lower, upper = median - n_mads * mad, median + n_mads * mad
        rows.append(
            {
                "metric": metric,
                "n_cells": int(values.size),
                "median": median,
                "mad": mad,
                "n_mads": n_mads,
                "lower": lower,
                "upper": upper,
                "n_below": int((values < lower).sum()) if mad > 0 else 0,
                "n_above": int((values > upper).sum()) if mad > 0 else 0,
            }
        )
    return pd.DataFrame(rows)


#: The metrics the MAD rule is *computed* on. Computed is not the same as filtered
#: on — see `DEFAULT_FILTERS`, which deliberately excludes one of these.
#: `pct_counts_mt` is handled separately because its rule is one-sided.
MAD_METRICS = (
    "log1p_total_counts",
    "log1p_n_genes_by_counts",
    "pct_counts_in_top_20_genes",
)

#: Default MAD counts: the standard 5-MAD procedure on all three metrics, which is
#: what `sc-best-practices` uses and is deliberately NOT tightened here. An earlier
#: draft put `pct_counts_in_top_20_genes` at 3 MADs on the argument that it is the
#: soup-detector and deserves a narrower band. That argument may well be right, but
#: making it *before looking* is the same error as copying the tutorial's numbers —
#: it just fails in the stricter direction. Notebook 04 computes what 3 vs 5 costs
#: per cohort and records the decision; this is the starting point it starts from.
#:
#: Note the band is two-sided. High `pct_counts_in_top_20_genes` is a low-complexity
#: droplet (in this tissue, usually ambient Ig); low is an unusually complex library,
#: which correlates with doublets. `mad_thresholds` reports the two sides separately
#: as `n_below`/`n_above` so the notebook can see which one is doing the work.
DEFAULT_N_MADS = {
    "log1p_total_counts": 5.0,
    "log1p_n_genes_by_counts": 5.0,
    "pct_counts_in_top_20_genes": 5.0,
}


#: Which flags compose `obs["outlier"]`. **`outlier_top20` is computed but excluded**,
#: and this is the stage's one real departure from the generic recipe. The evidence,
#: measured on this cohort rather than assumed:
#:
#: A 5-MAD upper bound on `pct_counts_in_top_20_genes` flags 17% of MMRF cells and
#: 15% of WashU 1's. Inspecting what those cells are (MMRF_1695, top decile) shows
#: they are dominated by `IGKC` (25% of counts) and haemoglobin (`HBB`+`HBA1/2`, 32%)
#: — i.e. two populations, plasma cells and erythroid debris, not one class of bad
#: droplet. And the plasma-cell half is the project's entire subject:
#:
#:     TNFRSF17 (BCMA) detected in 21.8% of the high-top20 decile vs 0.8% elsewhere
#:     SDC1     (CD138)              18.8%                    vs 0.0%
#:
#: A plasma cell is a professional secretor; a library dominated by a handful of
#: immunoglobulin transcripts is its *normal* state, not a defect. Filtering on this
#: metric would therefore preferentially delete antigen-POSITIVE malignant plasma
#: cells, which inflates `frac_double_negative` — biased in the direction the project
#: is looking, which is the worst kind of artifact to leave in.
#:
#: The metric is still computed, reported per cohort, and available: it is a genuine
#: ambient-Ig handle for stage 08, which needs one because SoupX/DecontX cannot run
#: on this deposit. It is just not allowed to delete cells here. The erythroid half is
#: `pct_counts_hb`'s business and, more properly, stage 06's.
DEFAULT_FILTERS = (
    "outlier_counts",
    "outlier_genes",
    "outlier_mt",
    "outlier_min_genes",
)

#: Every flag this module writes, filtering or not.
ALL_FLAGS = (*DEFAULT_FILTERS[:2], "outlier_top20", *DEFAULT_FILTERS[2:])


def flag_outliers(
    adata: AnnData,
    *,
    group_key: str = "cohort",
    n_mads: dict[str, float] | None = None,
    mt_n_mads: float = 3.0,
    mt_max_pct: float = 20.0,
    min_genes: int = 200,
    filters: Sequence[str] = DEFAULT_FILTERS,
) -> pd.DataFrame:
    """Flag low-quality cells, deriving every threshold **within `group_key`**.

    Returns the per-group threshold table and writes these boolean `obs` columns:

        outlier_counts, outlier_genes, outlier_top20, outlier_mt,
        outlier_min_genes, outlier      (the OR of `filters`, NOT of all five)

    `outlier_top20` is computed and reported but is **not** in `DEFAULT_FILTERS` —
    read that constant before changing it, the reason is specific to this tissue and
    is measured rather than assumed.

    `group_key="cohort"` is the whole point and is not a tunable preference — see the
    module docstring. Passing `None` pools, which exists only so the notebook can show
    the pooled result side by side and demonstrate what it would have cost.

    **Mitochondrial fraction is one-sided.** A cell with unusually FEW mitochondrial
    reads is not low quality, and a symmetric MAD band would discard it. The rule is
    `pct_counts_mt > median + mt_n_mads*MAD`, additionally capped at `mt_max_pct` so a
    cohort whose median is already high cannot license an absurd ceiling. Both the
    derived and the capped value are reported so it is visible which one bound.

    `min_genes` is an absolute floor applied on top of the MAD rule. It is deliberately
    low (200): its job is to remove cells that are unambiguously not cells, not to do
    the filtering, which the MAD rule does. A high floor here would filter WashU
    cohort 1 harder than MMRF for a depth reason, which is exactly what per-cohort
    thresholds exist to avoid.
    """
    if "log1p_total_counts" not in adata.obs:
        raise ValueError("Run add_qc_metrics() before flag_outliers().")

    n_mads = {**DEFAULT_N_MADS, **(n_mads or {})}
    obs = adata.obs
    if group_key is None:
        groups = pd.Series("all", index=obs.index)
    else:
        if group_key not in obs:
            raise ValueError(
                f"obs has no {group_key!r} column. Cohort reaches the object via "
                f"io.load_manifest(); did this object come from somewhere else?"
            )
        groups = obs[group_key].astype(str)

    flags = {
        name: pd.Series(False, index=obs.index)
        for name in ("outlier_counts", "outlier_genes", "outlier_top20", "outlier_mt")
    }
    metric_flag = dict(zip(MAD_METRICS, ("outlier_counts", "outlier_genes", "outlier_top20")))

    tables = []
    for group, index in obs.groupby(groups, observed=True).groups.items():
        block = obs.loc[index]
        # Each metric carries its own n_mads (see DEFAULT_N_MADS), so the table is
        # built one metric at a time rather than with a single shared value.
        table = pd.concat(
            [mad_thresholds(block, [m], n_mads=n_mads[m]) for m in MAD_METRICS],
            ignore_index=True,
        )
        for metric in MAD_METRICS:
            flags[metric_flag[metric]].loc[index] = mad_outlier(
                block[metric].to_numpy(), n_mads[metric]
            )

        mt = block["pct_counts_mt"].to_numpy(dtype=float)
        mt_median = float(np.median(mt))
        mt_mad = float(np.median(np.abs(mt - mt_median)))
        mt_derived = mt_median + mt_n_mads * mt_mad
        mt_cut = min(mt_derived, mt_max_pct)
        flags["outlier_mt"].loc[index] = mt > mt_cut

        table = pd.concat(
            [
                table,
                pd.DataFrame([{
                    "metric": "pct_counts_mt",
                    "n_cells": len(block),
                    "median": mt_median,
                    "mad": mt_mad,
                    "n_mads": mt_n_mads,
                    "lower": -np.inf,          # one-sided: low mt is not a defect
                    "upper": mt_cut,
                    "n_below": 0,
                    "n_above": int((mt > mt_cut).sum()),
                }]),
            ],
            ignore_index=True,
        )
        table.insert(0, "group", group)
        table["mt_derived_upper"] = mt_derived
        table["mt_capped_at"] = mt_max_pct
        table["mt_bound_by"] = "cap" if mt_derived > mt_max_pct else "mad"
        tables.append(table)

    for name, mask in flags.items():
        adata.obs[name] = mask.to_numpy()
    adata.obs["outlier_min_genes"] = (
        adata.obs["n_genes_by_counts"].to_numpy() < min_genes
    )
    unknown = sorted(set(filters) - set(ALL_FLAGS))
    if unknown:
        raise ValueError(f"unknown filter(s) {unknown}; known: {list(ALL_FLAGS)}")
    adata.obs["outlier"] = np.logical_or.reduce(
        [adata.obs[name].to_numpy() for name in filters]
    )
    adata.uns["qc_filters"] = list(filters)

    thresholds = pd.concat(tables, ignore_index=True)
    thresholds["group_key"] = "pooled" if group_key is None else group_key
    adata.uns["qc_thresholds"] = thresholds
    return thresholds


# ---------------------------------------------------------------------------
# Doublets — scDblFinder over the rpy2 bridge
# ---------------------------------------------------------------------------

_SCDBLFINDER_R = """
function(counts, seed) {
    set.seed(seed)
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    sce <- scDblFinder::scDblFinder(sce)
    data.frame(
        score = SingleCellExperiment::colData(sce)$scDblFinder.score,
        class = as.character(SingleCellExperiment::colData(sce)$scDblFinder.class),
        stringsAsFactors = FALSE
    )
}
"""


def detect_doublets(
    adata: AnnData,
    *,
    seed: int = 0,
    inplace: bool = True,
) -> pd.DataFrame:
    """Run scDblFinder on one sample and write `doublet_score` / `doublet_class`.

    Requires the `mm-qc` environment — R, rpy2 and bioconductor-scdblfinder live
    there and nowhere else, deliberately (see CLAUDE.md's env split). Raises
    `DoubletDetectionError` rather than returning empty results, because a sample
    that was never screened must not look like a clean one.

    **Run per sample, never on a concatenated object.** scDblFinder simulates
    doublets by combining cells from the same droplet suspension; pooling samples
    lets it manufacture cross-sample "doublets" that no droplet could contain, and
    it also loses each sample's own doublet rate, which scales with loading density
    and therefore with chemistry — the exact axis this cohort varies on.

    Counts are passed as an integer genes x cells matrix, R's orientation, which is
    the transpose of AnnData's. Getting that backwards produces a run that completes
    and returns nonsense.
    """
    try:
        import anndata2ri  # noqa: F401
        from rpy2.robjects import r
        from rpy2.robjects.conversion import localconverter
        import rpy2.robjects as ro
    except ImportError as error:
        raise DoubletDetectionError(
            f"scDblFinder needs the mm-qc environment (rpy2 + anndata2ri + "
            f"bioconductor-scdblfinder): {error}. Run this notebook on the mm-qc "
            f"kernel; mm-core carries no R."
        ) from error

    import scipy.sparse as sp

    counts = adata.X
    matrix = (counts if sp.issparse(counts) else sp.csr_matrix(counts)).T.tocsc()
    matrix = matrix.astype(np.int32)

    try:
        with localconverter(ro.default_converter + anndata2ri.scipy2ri.converter):
            r_counts = ro.conversion.py2rpy(matrix)
        result = r(_SCDBLFINDER_R)(r_counts, seed)
        with localconverter(ro.default_converter + ro.pandas2ri.converter):
            frame = ro.conversion.rpy2py(result)
    except Exception as error:  # rpy2 raises RRuntimeError and friends
        raise DoubletDetectionError(
            f"scDblFinder failed on {adata.obs['sample_name'].iat[0]}: {error}"
        ) from error

    if len(frame) != adata.n_obs:
        raise DoubletDetectionError(
            f"scDblFinder returned {len(frame)} rows for {adata.n_obs} cells — the "
            f"matrix was probably passed in the wrong orientation."
        )

    frame.index = adata.obs_names
    if inplace:
        adata.obs["doublet_score"] = frame["score"].to_numpy(dtype=float)
        adata.obs["doublet_class"] = pd.Categorical(frame["class"].astype(str))
    return frame


# ---------------------------------------------------------------------------
# Driving the stage
# ---------------------------------------------------------------------------

def _checkpoint_path(directory: Path, sample_name: str) -> Path:
    return Path(directory) / f"{sample_name}.h5ad"


def run_sample_qc(
    sample: str,
    manifest: pd.DataFrame | None = None,
    *,
    checkpoint_dir: Path | None = None,
    doublets: bool = True,
    overwrite: bool = False,
    seed: int = 0,
    **flag_kwargs,
) -> AnnData:
    """QC one sample and checkpoint it, resuming if the checkpoint already exists.

    Per-sample checkpointing is the resumability mechanism: 62 samples with an R
    call each is long enough that a crash 50 samples in must not cost the first 50.
    `overwrite=True` forces a recompute.

    **Cells are annotated, not deleted.** The object written to disk holds every
    barcode with `obs["outlier"]`, `obs["doublet_class"]` and `obs["keep"]` set. That
    costs a little disk and buys the ability to ask, at stage 08, whether a result
    depends on the QC — which for a metric that is a fraction of zeros is a question
    that will be asked. `run_cohort_qc(..., filter=True)` applies the mask.

    Note the MAD thresholds written here are **within-sample**, which is not what the
    stage ships: `run_cohort_qc` re-derives them per cohort across all samples, which
    needs every sample's metrics at once. The per-sample pass exists to compute the
    metrics and the doublet calls, both of which are genuinely per sample.
    """
    checkpoint_dir = Path(checkpoint_dir or config.RESULTS_DIR / "04_qc" / "samples")
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    path = _checkpoint_path(checkpoint_dir, sample)

    if path.exists() and not overwrite:
        import anndata
        return anndata.read_h5ad(path)

    adata = io.read_sample(sample, manifest=manifest)
    add_qc_metrics(adata)
    flag_outliers(adata, group_key=None, **flag_kwargs)

    if doublets:
        detect_doublets(adata, seed=seed)
        adata.obs["is_doublet"] = (
            adata.obs["doublet_class"].astype(str) == "doublet"
        )
    else:
        adata.obs["doublet_score"] = np.nan
        adata.obs["doublet_class"] = pd.Categorical(["not_run"] * adata.n_obs)
        adata.obs["is_doublet"] = False

    adata.obs["keep"] = ~(adata.obs["outlier"] | adata.obs["is_doublet"])
    adata.write_h5ad(path)
    return adata


def collect_obs(
    checkpoint_dir: Path | None = None,
    samples: Iterable[str] | None = None,
) -> pd.DataFrame:
    """Concatenate just the `obs` tables of the per-sample checkpoints.

    The cohort threshold pass needs every cell's QC metrics and nothing else, so it
    reads `obs` rather than the matrices. That keeps this stage memory-flat — the
    alternative, concatenating 62 count matrices to compute four medians, would cost
    several GB for no gain, and stage 05 reads the checkpoints directly anyway.
    """
    import anndata

    checkpoint_dir = Path(checkpoint_dir or config.RESULTS_DIR / "04_qc" / "samples")
    paths = (
        [_checkpoint_path(checkpoint_dir, name) for name in samples]
        if samples is not None
        else sorted(checkpoint_dir.glob("*.h5ad"))
    )
    missing = [str(p) for p in paths if not p.exists()]
    if missing:
        raise FileNotFoundError(
            f"{len(missing)} checkpoint(s) absent, e.g. {missing[:3]}. "
            f"Run run_sample_qc() for them first."
        )
    frames = [anndata.read_h5ad(p, backed="r").obs for p in paths]
    obs = pd.concat(frames, axis=0)
    if not obs.index.is_unique:
        raise ValueError(
            "cell ids collided across checkpoints — read_sample prefixes them with "
            "sample_name, so this means two samples share a name."
        )
    return obs


def cohort_thresholds(
    obs: pd.DataFrame,
    *,
    group_key: str = "cohort",
    n_mads: dict[str, float] | None = None,
    mt_n_mads: float = 3.0,
    mt_max_pct: float = 20.0,
    min_genes: int = 200,
    filters: Sequence[str] = DEFAULT_FILTERS,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Derive the outlier thresholds per `group_key` from a pooled `obs`.

    Returns `(flags, thresholds)`: a boolean frame indexed like `obs`, and the
    reportable per-group threshold table.

    Wraps `flag_outliers` on a metadata-only AnnData so there is exactly one
    implementation of the rule — the shell object holds `obs` and a zero-width `X`,
    which is legal and costs nothing.

    `group_key="cohort"` is the stage's central decision, not a tunable preference:
    the cohorts differ ~1.9x in genes/cell, so a pooled MAD would filter WashU
    cohort 1 (23 of the 54 myeloma samples) harder than MMRF for a batch reason.
    Pass `None` to pool, which exists so the notebook can show what that would cost.
    """
    shell = AnnData(
        X=np.zeros((len(obs), 0), dtype=np.float32),
        obs=obs.copy(),
        var=pd.DataFrame(index=pd.Index([], name="deposited_symbol")),
    )
    thresholds = flag_outliers(
        shell,
        group_key=group_key,
        n_mads=n_mads,
        mt_n_mads=mt_n_mads,
        mt_max_pct=mt_max_pct,
        min_genes=min_genes,
        filters=filters,
    )
    columns = [*ALL_FLAGS, "outlier"]
    return shell.obs[columns].copy(), thresholds


def run_cohort_qc(
    manifest: pd.DataFrame | None = None,
    samples: Iterable[str] | None = None,
    *,
    checkpoint_dir: Path | None = None,
    doublets: bool = True,
    overwrite: bool = False,
    group_key: str = "cohort",
    write_back: bool = True,
    **flag_kwargs,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run the whole stage: per-sample QC, then per-cohort thresholds.

    Returns `(obs, thresholds)` — the pooled cell-level QC table and the per-cohort
    threshold table. It deliberately does **not** return a concatenated matrix; see
    `collect_obs` for why, and `load_checkpoints` if one is actually wanted.

    Two passes, and the split is the point:

      1. **per sample** — QC metrics and doublet calls, both genuinely per-sample
         quantities. scDblFinder in particular must not see pooled samples (it would
         simulate doublets from cells that never shared a droplet). Checkpointed, so
         a crash 50 samples in does not cost the first 50.
      2. **per cohort** — the MAD thresholds, which need all of a cohort's cells at
         once. Not per sample (a small sample would get a threshold fitted to noise)
         and not pooled across cohorts (the 1.9x depth gap).

    `write_back=True` updates each checkpoint's `obs` with the cohort-derived flags
    and the final `keep` mask, so the files stage 05 reads carry the thresholds the
    stage actually shipped rather than the per-sample ones from pass 1.

    **Cells are annotated, never deleted.** Every checkpoint keeps all its barcodes
    with `keep` set. For a metric that is a fraction of zeros, "does this survive a
    different QC?" is a question that will be asked, and it can only be answered if
    the filtered cells are still there.
    """
    import anndata

    manifest = manifest if manifest is not None else io.load_manifest()
    names = list(samples) if samples is not None else manifest["sample_name"].tolist()
    checkpoint_dir = Path(checkpoint_dir or config.RESULTS_DIR / "04_qc" / "samples")

    for name in names:
        run_sample_qc(
            name,
            manifest=manifest,
            checkpoint_dir=checkpoint_dir,
            doublets=doublets,
            overwrite=overwrite,
            **flag_kwargs,
        )

    obs = collect_obs(checkpoint_dir, samples=names)
    flags, thresholds = cohort_thresholds(obs, group_key=group_key, **flag_kwargs)
    obs = obs.drop(columns=flags.columns).join(flags)
    obs["keep"] = ~(obs["outlier"] | obs["is_doublet"])

    if write_back:
        for name in names:
            path = _checkpoint_path(checkpoint_dir, name)
            adata = anndata.read_h5ad(path)
            block = obs.loc[adata.obs_names]
            for column in [*flags.columns, "keep"]:
                adata.obs[column] = block[column].to_numpy()
            adata.uns["qc_thresholds"] = thresholds
            adata.uns["qc_threshold_scope"] = group_key or "pooled"
            adata.write_h5ad(path)

    return obs, thresholds


def qc_report(obs: pd.DataFrame | AnnData, *, by: str = "sample_name") -> pd.DataFrame:
    """Per-`by` summary of what QC did and what it left.

    Takes the pooled `obs` that `run_cohort_qc` returns (an AnnData is accepted and
    its `.obs` used). One row per sample or cohort, with the pre/post counts, each
    filter's own contribution, and the post-QC depth medians.

    The per-filter columns **overlap** — a cell can be flagged by several — so they
    do not sum to `n_removed`, which is their union and is the honest total. Reported
    separately anyway, because which filter is doing the work is the thing worth
    knowing when a cohort's removal rate looks wrong.
    """
    obs = obs.obs if isinstance(obs, AnnData) else obs
    grouped = obs.groupby(by, observed=True)
    report = pd.DataFrame({
        "n_cells_pre": grouped.size(),
        "n_outlier_counts": grouped["outlier_counts"].sum(),
        "n_outlier_genes": grouped["outlier_genes"].sum(),
        "n_outlier_top20": grouped["outlier_top20"].sum(),
        "n_outlier_mt": grouped["outlier_mt"].sum(),
        "n_outlier_min_genes": grouped["outlier_min_genes"].sum(),
        "n_doublet": grouped["is_doublet"].sum(),
        "n_kept": grouped["keep"].sum(),
    })
    report["n_removed"] = report["n_cells_pre"] - report["n_kept"]
    report["pct_removed"] = 100 * report["n_removed"] / report["n_cells_pre"]

    kept = obs.loc[obs["keep"]].groupby(by, observed=True)
    report["median_genes_post"] = kept["n_genes_by_counts"].median()
    report["median_counts_post"] = kept["total_counts"].median()
    report["median_pct_mt_post"] = kept["pct_counts_mt"].median()

    for column in ("cohort", "chemistry", "sample_type"):
        if column in obs and column != by:
            first = grouped[column].agg(lambda s: s.astype(str).iat[0])
            report.insert(0, column, first)
    return report.reset_index()


def load_checkpoints(
    checkpoint_dir: Path | None = None,
    samples: Iterable[str] | None = None,
) -> AnnData:
    """Concatenate the per-sample checkpoints written by `run_sample_qc`."""
    import anndata

    checkpoint_dir = Path(checkpoint_dir or config.RESULTS_DIR / "04_qc" / "samples")
    paths = (
        [_checkpoint_path(checkpoint_dir, name) for name in samples]
        if samples is not None
        else sorted(checkpoint_dir.glob("*.h5ad"))
    )
    missing = [str(p) for p in paths if not p.exists()]
    if missing:
        raise FileNotFoundError(
            f"{len(missing)} checkpoint(s) absent, e.g. {missing[:3]}. "
            f"Run run_sample_qc() for them first."
        )
    return anndata.concat(
        [anndata.read_h5ad(p) for p in paths], join="inner", index_unique=None
    )
