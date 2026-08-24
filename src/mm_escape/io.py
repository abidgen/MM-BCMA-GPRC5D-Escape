"""
Loading GSE223060 from disk, and the sample metadata derivable from filenames.

WHY THIS MODULE EXISTS
----------------------
`scanpy.read_10x_mtx()` cannot read this deposit. The archive is old Cell Ranger
v2-style, uncompressed, one directory level deeper than the GSM folder, and the
matrix file is not called `matrix.mtx`:

    raw/samples/<GSM>_<sample>/<sample>/counts.mtx     (NOT matrix.mtx)
                                       /genes.tsv      (NOT features.tsv; ONE column)
                                       /barcodes.tsv   (uncompressed, no .gz)

`read_10x_mtx()` hardcodes the filenames it looks for, so pointing it at the right
directory fails to find anything. Everything here is built on explicit paths taken
from `raw/sample_manifest.csv` (produced by `scripts/03_build_manifest.py`, whose
column schema is the contract) instead of on any 10x-directory auto-detection.

`genes.tsv` is a single column of gene symbols — there is no `gene_id` column and no
`ENSG` string anywhere in the deposit. IDs are recovered separately and positionally
by `gene_space.attach_ensembl_ids`, which is why this module must return `var_names`
as the **deposited symbols in deposited order** and must never sort, dedupe or
otherwise touch the gene axis. `var_names_make_unique()` is banned for the same
reason — see `gene_space`.

ORIENTATION
-----------
The .mtx is genes x cells (Cell Ranger's orientation), e.g. `33538 1007 2604146`.
AnnData is cells x genes, so the matrix is transposed on read. Getting this backwards
produces an object that still looks plausible, which is why the dimensions are checked
against the two .tsv line counts rather than trusted.

PROVISIONAL PATIENT MAPPING
---------------------------
`patient_id` here is the naive rule (strip a trailing `_<digits>` only when the stem
is purely numeric), NOT the real mapping — that needs Supplementary Table S1, which is
still unresolved. Every object carries `obs["patient_id_source"] == "naive"` so a
provisional aggregate can never be mistaken for a final one, and so the switch to S1
is a one-line change with a visible marker. See the S1 policy in CLAUDE.md.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable, Iterator

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from anndata import AnnData

from . import config

__all__ = [
    "load_manifest",
    "parse_sample_name",
    "naive_patient_id",
    "classify_sample",
    "read_sample",
    "read_samples",
    "SampleLoadError",
]


class SampleLoadError(RuntimeError):
    """Raised when a sample cannot be loaded as the pipeline requires.

    Always names the offending sample and file — a sample that loads into a
    plausible-looking but wrong object is the failure mode this module guards.
    """


# `GSM6939028_MMRF_1695` -> ("GSM6939028", "MMRF_1695")
_GSM_PREFIX_RE = re.compile(r"^(GSM\d+)_(.+)$")

#: Normal bone marrow controls. Confirmed by name convention; used by stage 07 as the
#: malignant-caller negative control (polyclonal marrow must yield no clone).
_NORMAL_BM_RE = re.compile(r"^BM\d+$")

#: `ND_083017`, `ND_090617`, `ND_170531`, `ND_170607`. The suffixes are collection
#: DATES, not patient identifiers, which is how donor samples are usually labelled —
#: and CLAUDE.md's stage 07 explicitly lists `ND_*` alongside `BM*` as the normal-BM
#: controls. Treated as controls here, but flagged `sample_type_certain == False`:
#: `ND` could also read as "newly diagnosed", and the deposit does not say. S1 settles
#: it. This choice moves the naive disease-sample count 57 -> 53 and the provisional
#: patient count 47 -> 43; see `load_manifest`'s docstring.
_NORMAL_DONOR_RE = re.compile(r"^ND_\d+$")

#: The naive patient rule. Strips a trailing `_<digits>` ONLY when the stem is purely
#: numeric: `27522_1` -> `27522`, but `MMRF_1695` stays whole (its stem is `MMRF`).
_NUMERIC_SUFFIX_RE = re.compile(r"^(\d+)_(\d+)$")


# ---------------------------------------------------------------------------
# Sample metadata derivable from names alone
# ---------------------------------------------------------------------------

def parse_sample_name(sample_id: str) -> tuple[str, str]:
    """Split a manifest `sample_id` into (gsm_id, sample_name).

    `GSM6939028_MMRF_1695` -> `("GSM6939028", "MMRF_1695")`. Splitting on the first
    underscore is not enough in general, so the GSM accession is matched explicitly.
    """
    match = _GSM_PREFIX_RE.match(sample_id)
    if match is None:
        raise SampleLoadError(
            f"sample_id {sample_id!r} does not look like <GSM accession>_<sample name>. "
            f"The manifest's sample_id column is the directory name under raw/samples/."
        )
    return match.group(1), match.group(2)


def naive_patient_id(sample_name: str) -> str:
    """Provisional patient ID — NOT the real mapping.

    Strips a trailing `_<digits>` only when the stem is purely numeric, so the
    six `27522_*` samples collapse to one patient while `MMRF_1695` is left alone.

    This yields 43 provisional patients from 53 disease samples against the paper's
    reported 41 / 53, so roughly two name collapses are still being missed
    (`83942` / `MMY83942` is the obvious candidate). Do NOT hand-fix those here —
    guessing a patient mapping silently changes the denominator of the headline
    metric. Resolve against Supplementary Table S1; see CLAUDE.md's S1 policy.

    The `_N` suffixes are also of *unknown meaning* (timepoint vs. fraction vs. sort
    vs. replicate) — the bulk/sc suffix misalignment is evidence against the naive
    timepoint reading. This function only assumes they belong to one patient.
    """
    match = _NUMERIC_SUFFIX_RE.match(sample_name)
    return match.group(1) if match else sample_name


def classify_sample(sample_name: str) -> tuple[str, bool]:
    """Return (sample_type, certain) for a sample name.

    `sample_type` is `"normal_bm"` or `"myeloma"`. `certain` is False for `ND_*`,
    whose classification is a reading of the naming convention rather than a
    documented fact — see `_NORMAL_DONOR_RE`.
    """
    if _NORMAL_BM_RE.match(sample_name):
        return "normal_bm", True
    if _NORMAL_DONOR_RE.match(sample_name):
        return "normal_bm", False
    return "myeloma", True


# ---------------------------------------------------------------------------
# The manifest
# ---------------------------------------------------------------------------

def load_manifest(
    path: Path | None = None,
    *,
    drop_excluded: bool = True,
    require_complete: bool = True,
) -> pd.DataFrame:
    """Load `raw/sample_manifest.csv` and add the metadata derivable from names.

    The manifest stores **repo-root-relative** paths (it is committed and read on
    other machines), so they are resolved to absolute here and their existence is
    checked. Columns added on top of the script's schema:

        gsm_id, sample_name, patient_id, patient_id_source,
        sample_type, sample_type_certain, excluded

    `drop_excluded` removes `56203_1` (22184-gene reference, missing `TNFRSF17`
    entirely — every cell would read BCMA-negative for a purely technical reason).
    Patient 56203 is fully covered by `56203_2`, so no patient coverage is lost.
    Pass False only to inspect the exclusion, never to analyse it.

    `require_complete` hard-fails on any row not classified `triplet-ok` by the
    manifest script, rather than discovering the problem mid-load.

    NOTE on counts: 62 rows in, 61 after the exclusion, of which 8 are normal-BM
    controls (`BM2/4/5/6` plus the four `ND_*`, the latter uncertain) and 53 are
    disease samples mapping to 43 provisional patients. CLAUDE.md's inherited
    "47 patients / 57 samples" counts the four `ND_*` as disease; if `ND` does mean
    normal donor, the naive mapping is only ~2 collapses short of the paper's 41,
    not ~6.
    """
    manifest_path = Path(path) if path is not None else config.MANIFEST_CSV
    if not manifest_path.exists():
        raise SampleLoadError(
            f"No manifest at {manifest_path}. Build it with "
            f"`python scripts/03_build_manifest.py raw/samples` or "
            f"`notebooks/03_build_manifest.py`."
        )

    frame = pd.read_csv(manifest_path)

    if require_complete:
        bad = frame.loc[frame["format"] != "triplet-ok", "sample_id"].tolist()
        if bad:
            raise SampleLoadError(
                f"{len(bad)} sample(s) are not 'triplet-ok' in {manifest_path}: {bad}. "
                f"This loader reads the counts.mtx/genes.tsv/barcodes.tsv triplet only; "
                f"re-run the manifest builder or inspect those directories."
            )

    parsed = frame["sample_id"].map(parse_sample_name)
    frame["gsm_id"] = [gsm for gsm, _ in parsed]
    frame["sample_name"] = [name for _, name in parsed]
    frame["patient_id"] = frame["sample_name"].map(naive_patient_id)
    frame["patient_id_source"] = "naive"

    classified = frame["sample_name"].map(classify_sample)
    frame["sample_type"] = [t for t, _ in classified]
    frame["sample_type_certain"] = [c for _, c in classified]

    frame["excluded"] = frame["sample_name"].isin(config.EXCLUDED_SAMPLES)

    # The manifest holds repo-root-relative paths; resolve and verify them here so a
    # missing file surfaces now rather than 40 samples into a loop.
    root = config.REPO_ROOT
    for column in ("matrix_path", "barcodes_path", "genefeat_path"):
        frame[column] = [
            "" if not value else str((root / value).resolve())
            for value in frame[column].fillna("")
        ]

    if drop_excluded:
        frame = frame.loc[~frame["excluded"]].reset_index(drop=True)

    missing = [
        (row.sample_id, column)
        for row in frame.itertuples()
        for column in ("matrix_path", "barcodes_path", "genefeat_path")
        if not Path(getattr(row, column)).exists()
    ]
    if missing:
        raise SampleLoadError(
            f"{len(missing)} manifest path(s) do not exist on disk, e.g. {missing[:3]}. "
            f"Was raw/ extracted on this machine? See notebooks/01_download_data.py."
        )

    return frame


# ---------------------------------------------------------------------------
# Reading one sample
# ---------------------------------------------------------------------------

def _read_single_column(path: Path, what: str, sample_id: str) -> list[str]:
    """Read a one-column .tsv into a list of stripped strings.

    `genes.tsv` here is genuinely one column (no Ensembl ID), and `barcodes.tsv` is
    too. A file that turns out to have more columns is a different deposit than the
    one this pipeline was verified against, so it raises instead of silently taking
    the first field.
    """
    values: list[str] = []
    with open(path) as handle:
        for lineno, line in enumerate(handle, start=1):
            value = line.rstrip("\n").rstrip("\r")
            if "\t" in value:
                raise SampleLoadError(
                    f"{sample_id}: {what} file {path} has more than one column at line "
                    f"{lineno} ({value!r}). This deposit's {what} file is single-column; "
                    f"a multi-column file means a different reference format."
                )
            values.append(value)
    if not values:
        raise SampleLoadError(f"{sample_id}: {what} file {path} is empty.")
    return values


def read_sample(
    sample: str | pd.Series,
    manifest: pd.DataFrame | None = None,
    *,
    dtype: str = "float32",
) -> AnnData:
    """Read one sample's counts.mtx/genes.tsv/barcodes.tsv triplet into an AnnData.

    `sample` is either a manifest `sample_id` (`GSM6939028_MMRF_1695`), a bare
    `sample_name` (`MMRF_1695`), or a manifest row. Passing `manifest` avoids
    re-reading the CSV in a loop.

    Returns cells x genes with:
        X            raw integer counts, CSR, as `dtype` (float32 by default — the
                     values stay integral, but int64 would triple the memory and
                     scanpy normalizes to float on the first pp step anyway)
        var_names    the DEPOSITED symbols in DEPOSITED ORDER, untouched. This is a
                     precondition of gene_space.attach_ensembl_ids, which verifies
                     them position-for-position.
        obs_names    `<sample_name>_<barcode>` — barcodes repeat across samples, so
                     they must be disambiguated before any concat. The original is
                     kept in obs["barcode"].
        obs          sample_id, gsm_id, sample_name, patient_id, patient_id_source,
                     sample_type, sample_type_certain, n_genes_ref
        uns          the source paths and the reference row count

    `n_genes_ref` is carried in `.obs` (not only `.uns`) because it survives
    concatenation and stage 05 needs it as a Harmony covariate.
    """
    if manifest is None:
        manifest = load_manifest()

    if isinstance(sample, pd.Series):
        row = sample
    else:
        hits = manifest.loc[
            (manifest["sample_id"] == sample) | (manifest["sample_name"] == sample)
        ]
        if len(hits) == 0:
            if sample in config.EXCLUDED_SAMPLES:
                raise SampleLoadError(
                    f"{sample!r} is on the exclusion list and is not in the manifest. "
                    f"It was processed against a 22184-gene reference that lacks "
                    f"TNFRSF17 (BCMA) — every cell would read BCMA-negative for a "
                    f"technical reason. Patient 56203 is covered by 56203_2."
                )
            raise SampleLoadError(f"{sample!r} matches no manifest row.")
        if len(hits) > 1:
            raise SampleLoadError(
                f"{sample!r} matches {len(hits)} manifest rows "
                f"({hits['sample_id'].tolist()}); pass the full sample_id."
            )
        row = hits.iloc[0]

    sample_id = str(row["sample_id"])
    if bool(row.get("excluded", False)):
        raise SampleLoadError(
            f"{sample_id} is excluded (see config.EXCLUDED_SAMPLES) and must not be "
            f"loaded for analysis."
        )

    matrix_path = Path(row["matrix_path"])
    genes_path = Path(row["genefeat_path"])
    barcodes_path = Path(row["barcodes_path"])

    symbols = _read_single_column(genes_path, "genes", sample_id)
    barcodes = _read_single_column(barcodes_path, "barcodes", sample_id)

    matrix = scipy.io.mmread(matrix_path)
    if matrix.shape != (len(symbols), len(barcodes)):
        raise SampleLoadError(
            f"{sample_id}: {matrix_path.name} is {matrix.shape[0]}x{matrix.shape[1]} but "
            f"genes.tsv has {len(symbols)} rows and barcodes.tsv has {len(barcodes)}. "
            f"Cell Ranger writes genes x cells; a swap here would transpose the whole "
            f"object into something that still looks plausible."
        )

    # genes x cells -> cells x genes. CSR after the transpose, not before.
    counts = sp.csr_matrix(matrix.T, dtype=np.dtype(dtype))

    # The deposited column already went through R's make.unique, so duplicates mean
    # the file is not what gene_space's positional reconstruction was verified against.
    if len(set(symbols)) != len(symbols):
        duplicated = sorted({s for s in symbols if symbols.count(s) > 1})[:5]
        raise SampleLoadError(
            f"{sample_id}: genes.tsv has duplicate symbols (e.g. {duplicated}). The "
            f"deposited column is make.unique'd and must be unique. Do NOT call "
            f"var_names_make_unique() to paper over this — see gene_space."
        )
    if len(set(barcodes)) != len(barcodes):
        raise SampleLoadError(f"{sample_id}: barcodes.tsv has duplicate barcodes.")

    sample_name = str(row["sample_name"])
    obs = pd.DataFrame(
        {
            "barcode": barcodes,
            "sample_id": sample_id,
            "gsm_id": str(row["gsm_id"]),
            "sample_name": sample_name,
            "patient_id": str(row["patient_id"]),
            "patient_id_source": str(row["patient_id_source"]),
            "sample_type": str(row["sample_type"]),
            "sample_type_certain": bool(row["sample_type_certain"]),
            "n_genes_ref": len(symbols),
        },
        index=pd.Index([f"{sample_name}_{bc}" for bc in barcodes], name="cell_id"),
    )
    for column in (
        "sample_id", "gsm_id", "sample_name", "patient_id",
        "patient_id_source", "sample_type", "n_genes_ref",
    ):
        obs[column] = obs[column].astype("category")

    adata = AnnData(
        X=counts,
        obs=obs,
        var=pd.DataFrame(index=pd.Index(symbols, name="deposited_symbol")),
    )
    adata.uns["source"] = {
        "matrix": str(matrix_path),
        "genes": str(genes_path),
        "barcodes": str(barcodes_path),
        "n_genes_ref": len(symbols),
    }
    return adata


def read_samples(
    samples: Iterable[str] | None = None,
    manifest: pd.DataFrame | None = None,
    *,
    verbose: bool = True,
    **kwargs,
) -> Iterator[AnnData]:
    """Yield one AnnData per sample, in manifest order.

    A generator on purpose: the full cohort is ~181k cells across 61 samples and
    stage 04 checkpoints each sample's post-QC object individually, so nothing needs
    all 61 raw matrices resident at once.
    """
    if manifest is None:
        manifest = load_manifest()
    if samples is not None:
        wanted = set(samples)
        manifest = manifest.loc[
            manifest["sample_id"].isin(wanted) | manifest["sample_name"].isin(wanted)
        ]
        found = set(manifest["sample_id"]) | set(manifest["sample_name"])
        unknown = sorted(wanted - found)
        if unknown:
            raise SampleLoadError(f"Unknown sample(s): {unknown}")

    total = len(manifest)
    for i, (_, row) in enumerate(manifest.iterrows(), start=1):
        adata = read_sample(row, manifest=manifest, **kwargs)
        if verbose:
            print(
                f"[{i:>2}/{total}] {row['sample_id']:<26} "
                f"{adata.n_obs:>5} cells x {adata.n_vars} genes  "
                f"({row['sample_type']})"
            )
        yield adata
