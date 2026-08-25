"""
Stage 05 — gene-space harmonization, integration, clustering.

WHAT THIS STAGE PRODUCES
------------------------
One harmonized AnnData over every QC-passing cell, with a batch-corrected embedding
and Leiden clusters. Stage 06 annotates it; stage 07 subsets it.

THE BLAST RADIUS OF INTEGRATION IS DELIBERATELY CONSTRAINED
-----------------------------------------------------------
Harmony keyed on `patient_id` is correct for the immune compartment and carries a
real risk for the tumor: **the malignant clone is patient-private by definition**, so
forcing patients together can blend genuinely distinct clones into one blob and erase
the heterogeneity this project exists to measure. Therefore:

- The integrated embedding (`obsm["X_pca_harmony"]`) is for **immune-compartment
  annotation and clustering only** — stages 06 and 11.
- **All malignant subclustering is per patient and un-integrated** (stage 10), and
  must not read this embedding.
- **Per-cell antigen calls are raw counts and are therefore
  integration-independent.** This is what contains the risk, and it is the answer to
  "did Harmony distort your escape fractions?" — it cannot, because the calls never
  touch the embedding.

WHY THE GENE SPACE IS INTERSECTED ON ENSEMBL IDs
------------------------------------------------
Two Cell Ranger references are in play (33538 and 33694 genes) with different HGNC
symbol vintages. A symbol join keeps 22,164 genes and silently *mis-pairs* some of
them — `TBCE` is a different Ensembl entry in each build. The ID join keeps 32,991
and pairs them correctly. `gene_space` does the work and verifies itself; this module
only drives it and reports what it recovered.

Intersect, never union: a union would make ~11k genes structurally zero in whole
sample cohorts, which downstream is indistinguishable from a true biological zero —
exactly the quantity this project measures.

CELLS ARE FILTERED HERE, NOT AT STAGE 04
----------------------------------------
The stage-04 checkpoints keep every barcode with `obs["keep"]` set, so that "does
this result survive a different QC?" stays answerable. `load_qc_checkpoints` applies
the mask; pass `keep_only=False` to load everything for such a re-run.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from . import config, gene_space

__all__ = [
    "load_qc_checkpoints",
    "build_gene_space",
    "normalize_and_hvg",
    "run_pca_harmony",
    "cluster_and_embed",
    "batch_mixing",
    "composition_table",
]

#: Harmony's batch key, and the covariates carried alongside it.
#:
#: `patient_id` is the batch. `n_genes_ref` and `cohort` are **different axes and
#: neither substitutes for the other** — the reference-build split cuts across
#: cohorts (two WU1 samples sit on 33538, the four `ND_*` donors on 33694), and the
#: ~1.9x depth gap follows cohort rather than build. Dropping either leaves a real
#: batch structure uncorrected.
HARMONY_KEYS = ("patient_id", "n_genes_ref", "cohort")


def load_qc_checkpoints(
    checkpoint_dir: Path | None = None,
    samples: Iterable[str] | None = None,
    *,
    keep_only: bool = True,
    attach_ids: bool = True,
    verbose: bool = True,
) -> list[AnnData]:
    """Load the stage-04 checkpoints, one per sample, ready for the gene-space join.

    Returns a list rather than a concatenated object because `attach_ensembl_ids`
    must run **per sample, before any concat** — it verifies the deposited symbol
    column position-for-position against the committed map for that sample's build,
    and after a concat there is no single build to verify against.

    `keep_only=True` (the default) applies `obs["keep"]`; the checkpoints themselves
    still hold every barcode.
    """
    import anndata

    checkpoint_dir = Path(checkpoint_dir or config.RESULTS_DIR / "04_qc" / "samples")
    paths = (
        [checkpoint_dir / f"{name}.h5ad" for name in samples]
        if samples is not None
        else sorted(checkpoint_dir.glob("*.h5ad"))
    )
    missing = [str(p) for p in paths if not p.exists()]
    if missing:
        raise FileNotFoundError(
            f"{len(missing)} stage-04 checkpoint(s) absent, e.g. {missing[:3]}. "
            f"Run notebooks/04_qc.ipynb first."
        )

    blocks: list[AnnData] = []
    n_pre = n_post = 0
    for path in paths:
        adata = anndata.read_h5ad(path)
        n_pre += adata.n_obs
        if keep_only:
            if "keep" not in adata.obs:
                raise ValueError(
                    f"{path.name} has no obs['keep'] — it predates the stage-04 "
                    f"cohort pass. Re-run qc.run_cohort_qc()."
                )
            adata = adata[adata.obs["keep"].to_numpy()].copy()
        n_post += adata.n_obs
        if attach_ids:
            gene_space.attach_ensembl_ids(adata)
        blocks.append(adata)

    if verbose:
        tail = (f" -> {n_post:,} after obs['keep']" if keep_only else " (unfiltered)")
        print(f"loaded {len(blocks)} checkpoints: {n_pre:,} cells{tail}")
    return blocks


def build_gene_space(
    blocks: Sequence[AnnData],
    *,
    verbose: bool = True,
) -> AnnData:
    """Intersect on Ensembl IDs, concatenate, switch to symbols, assert the panel.

    The order is load-bearing and is the one `gene_space`'s docstrings specify:
    IDs are the index **through the merge only**, because that is where identity is
    at risk; once there is a single harmonized matrix the mis-pairing risk is gone
    and every downstream consumer (`score_genes`, dotplots, `celltypist`,
    `decoupler`, `liana`) is symbol-native.

    Reports how many genes the ID join recovered over a symbol join, because a
    regression there would otherwise be silent — and it is worth ~10,800 genes.
    """
    import anndata

    symbol_sets = [set(a.var["deposited_symbol"]) if "deposited_symbol" in a.var
                   else None for a in blocks]

    subset = gene_space.intersect_gene_space(blocks, verbose=verbose)
    adata = anndata.concat(subset, join="inner", index_unique=None, merge="same")
    if not adata.obs_names.is_unique:
        raise ValueError(
            "cell ids collided across samples — read_sample prefixes them with "
            "sample_name, so two samples must share a name."
        )

    gene_space.to_canonical_symbols(adata)
    gene_space.assert_required_genes(adata)

    if verbose and all(s is not None for s in symbol_sets):
        shared_symbols = set.intersection(*symbol_sets)  # type: ignore[arg-type]
        drifted = int(adata.var["symbol_drift"].sum())
        print(f"  a raw-symbol join would have kept {len(shared_symbols):,} genes; "
              f"the Ensembl-ID join keeps {adata.n_vars:,} "
              f"(+{adata.n_vars - len(shared_symbols):,})")
        print(f"  {drifted:,} of them carry a DIFFERENT symbol in each build and "
              f"were invisible to a symbol join")
    return adata


def normalize_and_hvg(
    adata: AnnData,
    *,
    n_top_genes: int = 2000,
    batch_key: str = "patient_id",
    target_sum: float = 1e4,
) -> AnnData:
    """Store raw counts, normalize, log1p, and select HVGs **within batch**.

    Raw counts are kept in `layers["counts"]` and are what stage 08 reads. That is
    not a convenience: antigen positivity is called on counts precisely so it is
    independent of everything this stage does.

    HVGs are selected with `batch_key` set. Without it the selection is dominated by
    genes that separate cohorts — a 1.9x depth gap and two chemistry generations
    produce plenty — and the integration would then be asked to correct a batch
    effect that the feature selection had just amplified.
    """
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    # No `adata.raw`. It exists so plotting and DE can reach genes that a subset to
    # HVGs removed — but this object is never subset: `run_pca_harmony` scales a
    # throwaway copy and writes back only the embedding, so `.X` keeps all 32,991
    # log-normalized genes. Setting `.raw` would duplicate a ~3 GB matrix to no end.

    sc.pp.highly_variable_genes(
        adata, n_top_genes=n_top_genes, batch_key=batch_key, flavor="seurat"
    )
    return adata


def run_pca_harmony(
    adata: AnnData,
    *,
    n_comps: int = 50,
    keys: Sequence[str] = HARMONY_KEYS,
    max_iter: int = 20,
    random_state: int = 0,
) -> AnnData:
    """Scale, PCA on HVGs, then Harmony on `keys`.

    Leaves both embeddings in place — `obsm["X_pca"]` uncorrected and
    `obsm["X_pca_harmony"]` corrected — so the diagnostic in `batch_mixing` can
    compare them rather than assert that integration worked.
    """
    from harmonypy import run_harmony

    missing = [key for key in keys if key not in adata.obs]
    if missing:
        raise ValueError(f"obs has no {missing}; see io.load_manifest().")

    subset = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(subset, max_value=10)
    sc.tl.pca(subset, n_comps=n_comps, svd_solver="arpack", random_state=random_state)
    adata.obsm["X_pca"] = subset.obsm["X_pca"]
    adata.uns["pca"] = subset.uns["pca"]

    meta = adata.obs[list(keys)].astype(str)
    harmony = run_harmony(
        adata.obsm["X_pca"], meta, list(keys), max_iter_harmony=max_iter,
        random_state=random_state,
    )
    # harmonypy 2.x returns `Z_corr` already as cells x dims; 1.x returned it
    # transposed. Assert rather than trusting either, because a wrong orientation
    # here produces an embedding that clusters happily and means nothing.
    corrected = np.asarray(harmony.Z_corr)
    if corrected.shape != (adata.n_obs, n_comps):
        corrected = corrected.T
    if corrected.shape != (adata.n_obs, n_comps):
        raise ValueError(
            f"Harmony returned {np.asarray(harmony.Z_corr).shape}, which is neither "
            f"({adata.n_obs}, {n_comps}) nor its transpose."
        )
    adata.obsm["X_pca_harmony"] = corrected
    adata.uns["harmony_keys"] = list(keys)
    return adata


def cluster_and_embed(
    adata: AnnData,
    *,
    use_rep: str = "X_pca_harmony",
    n_neighbors: int = 15,
    resolution: float = 1.0,
    key_added: str = "leiden",
    random_state: int = 0,
) -> AnnData:
    """Neighbor graph on the corrected embedding, then Leiden and UMAP."""
    sc.pp.neighbors(
        adata, n_neighbors=n_neighbors, use_rep=use_rep, random_state=random_state
    )
    sc.tl.leiden(
        adata, resolution=resolution, key_added=key_added,
        flavor="igraph", n_iterations=2, directed=False, random_state=random_state,
    )
    sc.tl.umap(adata, random_state=random_state)
    return adata


def batch_mixing(
    adata: AnnData,
    *,
    batch_key: str = "cohort",
    cluster_key: str = "leiden",
) -> pd.DataFrame:
    """Per-cluster batch composition and its normalized entropy.

    A cheap, honest integration diagnostic: entropy 1.0 means a cluster's batch mix
    matches the cohort's overall mix, 0.0 means it is one batch only. It answers
    "did the correction mix batches?" without claiming to be a benchmarked metric
    like kBET or iLISI.

    **A low value is not automatically a failure**, and that is the point of
    reporting it per cluster rather than as one number: a patient-private malignant
    clone SHOULD be low-entropy on `patient_id`. Read it against `cohort`, where
    low entropy has no biological excuse, and against `n_genes_ref`, where it has
    none either.
    """
    if batch_key not in adata.obs or cluster_key not in adata.obs:
        raise ValueError(f"need obs[{batch_key!r}] and obs[{cluster_key!r}].")

    counts = pd.crosstab(adata.obs[cluster_key], adata.obs[batch_key].astype(str))
    proportions = counts.div(counts.sum(axis=1), axis=0)
    n_batches = (counts.sum(axis=0) > 0).sum()

    with np.errstate(divide="ignore", invalid="ignore"):
        terms = np.where(proportions > 0, proportions * np.log(proportions), 0.0)
    entropy = -terms.sum(axis=1)
    normalized = entropy / np.log(n_batches) if n_batches > 1 else entropy * 0.0

    report = counts.copy()
    report["n_cells"] = counts.sum(axis=1)
    report["entropy"] = normalized
    report["dominant"] = proportions.idxmax(axis=1)
    report["dominant_pct"] = 100 * proportions.max(axis=1)
    return report.reset_index()


def composition_table(
    adata: AnnData,
    *,
    cluster_key: str = "leiden",
    by: str = "sample_name",
) -> pd.DataFrame:
    """Per-`by` cluster proportions — the input to stage 06's composition outputs.

    Proportions, not counts, because sample cell yields vary ~15x across this cohort
    and a raw count table would read as biology.
    """
    counts = pd.crosstab(adata.obs[by], adata.obs[cluster_key])
    return counts.div(counts.sum(axis=1), axis=0)
