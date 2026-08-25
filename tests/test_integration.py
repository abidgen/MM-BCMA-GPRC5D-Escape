"""Stage 05 — gene-space harmonization, integration, clustering.

Two tiers per `tests/conftest.py`. The gene-space contract and the batch-mixing
arithmetic need no deposit; anything reading a checkpoint sits behind `requires_data`.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from mm_escape import config, integration

from conftest import requires_data


# ---------------------------------------------------------------------------
# Configuration invariants — no data
# ---------------------------------------------------------------------------

def test_harmony_keys_carry_build_and_cohort_as_separate_axes():
    """Neither substitutes for the other, and dropping either leaves batch behind.

    The reference-build split cuts across cohorts (two WU1 samples on 33538, the
    four ND_* donors on 33694) and the ~1.9x depth gap follows cohort rather than
    build. See the constant's own note.
    """
    assert integration.HARMONY_KEYS == ("patient_id", "n_genes_ref", "cohort")


def test_run_pca_harmony_rejects_a_missing_covariate():
    adata = AnnData(
        X=np.random.default_rng(0).normal(size=(20, 5)).astype(np.float32),
        obs=pd.DataFrame({"patient_id": ["a"] * 20}, index=[f"c{i}" for i in range(20)]),
    )
    with pytest.raises(ValueError, match="n_genes_ref|cohort"):
        integration.run_pca_harmony(adata)


# ---------------------------------------------------------------------------
# batch_mixing — synthetic, no data
# ---------------------------------------------------------------------------

def _clustered(assignments):
    obs = pd.DataFrame(assignments, index=[f"c{i}" for i in range(len(assignments["leiden"]))])
    return AnnData(X=np.zeros((len(obs), 0), dtype=np.float32), obs=obs)


def test_batch_mixing_is_one_for_a_perfectly_mixed_cluster():
    adata = _clustered({"leiden": ["0"] * 4, "cohort": ["A", "B", "C", "D"]})
    row = integration.batch_mixing(adata).iloc[0]
    assert row["entropy"] == pytest.approx(1.0)
    assert row["n_cells"] == 4


def test_batch_mixing_is_zero_for_a_single_batch_cluster():
    adata = _clustered({
        "leiden": ["0"] * 4 + ["1"] * 4,
        "cohort": ["A"] * 4 + ["A", "B", "C", "D"],
    })
    report = integration.batch_mixing(adata).set_index("leiden")
    assert report.loc["0", "entropy"] == pytest.approx(0.0)
    assert report.loc["0", "dominant"] == "A"
    assert report.loc["0", "dominant_pct"] == pytest.approx(100.0)
    # Normalization is against the number of batches in the DATASET, not the
    # cluster, so a mixed cluster still reads 1.0 alongside a pure one.
    assert report.loc["1", "entropy"] == pytest.approx(1.0)


def test_batch_mixing_requires_the_columns_it_names():
    adata = _clustered({"leiden": ["0"], "cohort": ["A"]})
    with pytest.raises(ValueError, match="patient_id"):
        integration.batch_mixing(adata, batch_key="patient_id")


def test_composition_table_returns_proportions_not_counts():
    """Sample cell yields vary ~15x in this cohort; raw counts would read as biology."""
    adata = _clustered({
        "leiden": ["0", "1", "0", "0", "1", "1"],
        "sample_name": ["A", "A", "A", "B", "B", "B"],
        "cohort": ["X"] * 6,
    })
    table = integration.composition_table(adata)
    assert np.allclose(table.sum(axis=1), 1.0)
    assert table.loc["A", "0"] == pytest.approx(2 / 3)
    assert table.loc["B", "0"] == pytest.approx(1 / 3)


# ---------------------------------------------------------------------------
# Against the stage-04 checkpoints
# ---------------------------------------------------------------------------

CHECKPOINTS = config.RESULTS_DIR / "04_qc" / "samples"

requires_stage04 = pytest.mark.skipif(
    not CHECKPOINTS.is_dir() or not any(CHECKPOINTS.glob("*.h5ad")),
    reason="needs the stage-04 checkpoints (run notebooks/04_qc.ipynb)",
)


@requires_data
@requires_stage04
def test_checkpoints_still_satisfy_the_positional_gene_join():
    """`attach_ensembl_ids` is positional; a reorder anywhere upstream breaks it.

    Stage 04 writes new `var` columns, which is exactly the kind of change that
    could reorder the axis without anyone noticing until a gene silently means a
    different gene.
    """
    blocks = integration.load_qc_checkpoints(
        samples=["MMRF_1695", "27522_1"], verbose=False
    )
    for block in blocks:
        assert block.var_names.str.startswith("ENSG").all()
        assert "deposited_symbol" in block.var


@requires_data
@requires_stage04
def test_load_qc_checkpoints_filters_to_keep_by_default():
    """Cells are filtered HERE, not at stage 04 — the checkpoints hold everything."""
    filtered = integration.load_qc_checkpoints(samples=["BM4"], verbose=False)[0]
    everything = integration.load_qc_checkpoints(
        samples=["BM4"], keep_only=False, verbose=False
    )[0]
    assert filtered.n_obs < everything.n_obs
    assert filtered.n_obs == int(everything.obs["keep"].sum())
    assert filtered.obs["keep"].all()


@requires_data
@requires_stage04
def test_gene_space_recovers_the_documented_numbers():
    """32,991 genes, +10,827 over a symbol join, 11,140 drifted symbols.

    The four canonical samples span both reference builds, which is what makes this
    the same intersection the full cohort produces — the gene space is determined by
    which builds are present, not by how many samples.
    """
    blocks = integration.load_qc_checkpoints(
        samples=["MMRF_1695", "27522_1", "BM4", "56203_1"], verbose=False
    )
    symbol_join = set.intersection(*(set(b.var["deposited_symbol"]) for b in blocks))
    adata = integration.build_gene_space(blocks, verbose=False)

    assert adata.n_vars == 32_991
    assert len(symbol_join) == 22_164
    assert int(adata.var["symbol_drift"].sum()) == 11_140

    # The genes the whole project depends on, including the one symbol drift would
    # have silently dropped.
    for gene in ("TNFRSF17", "GPRC5D", "SLAMF7", "FCRL5", "SDC1", "CD38",
                 "IGKC", "NSD2"):
        assert gene in adata.var_names, gene
    # NSD2 is WHSC1 in the older build; that it resolved is the regression guard.
    nsd2 = adata.var.loc["NSD2"]
    assert nsd2["symbol_33694"] == "WHSC1"
    assert nsd2["symbol_33538"] == "NSD2"


@requires_data
@requires_stage04
def test_normalize_keeps_raw_counts_and_sets_no_redundant_raw():
    """Stage 08 reads `layers['counts']`; antigen calls must not depend on stage 05."""
    # Both builds, deliberately: intersecting a single build against itself yields
    # 33,538 genes, which is not the cohort gene space and which
    # `to_canonical_symbols` correctly refuses — see the test below.
    blocks = integration.load_qc_checkpoints(
        samples=["BM4", "27522_1"], verbose=False
    )
    adata = integration.build_gene_space(blocks, verbose=False)
    before = adata.X.sum()

    integration.normalize_and_hvg(adata, n_top_genes=500)
    assert "counts" in adata.layers
    assert adata.layers["counts"].sum() == before
    # X is log-normalized in place and holds every gene, so `.raw` would duplicate
    # a multi-GB matrix for nothing.
    assert adata.raw is None
    assert adata.n_vars == 32_991
    assert adata.var["highly_variable"].sum() == 500


@requires_data
@requires_stage04
def test_a_single_build_subset_is_refused_rather_than_silently_wrong():
    """Guard against the easy notebook mistake of exploring on one sample.

    Intersecting one 33538-build sample against itself gives 33,538 genes — 547 of
    which exist in that build only and have no place in the cohort's harmonized
    space. Producing them would be a quietly different gene set from every other
    stage's; `to_canonical_symbols` raises instead.
    """
    from mm_escape import gene_space

    blocks = integration.load_qc_checkpoints(samples=["BM4"], verbose=False)
    with pytest.raises(gene_space.GeneSpaceError, match="intersection table"):
        integration.build_gene_space(blocks, verbose=False)


@requires_data
@requires_stage04
def test_the_harmonized_object_survives_a_round_trip_through_h5ad(tmp_path):
    """Regression: the var index name must not collide with a var column.

    `to_canonical_symbols` names the index `symbol` and keeps `canonical_symbol` as
    a column, and for the 9 collision-suffixed genes the two differ. Naming both the
    same made `write_h5ad` raise — an error that appeared only on write, after the
    object had looked correct in memory through every in-memory test. Anything that
    only fails on serialization needs a test that serializes.
    """
    blocks = integration.load_qc_checkpoints(
        samples=["MMRF_1695", "27522_1"], verbose=False
    )
    adata = integration.build_gene_space(blocks, verbose=False)
    assert adata.var.index.name == "symbol"
    assert "canonical_symbol" in adata.var
    # The 9 collided symbols are exactly where index and column disagree.
    disagree = adata.var.index.to_numpy() != adata.var["canonical_symbol"].to_numpy()
    assert disagree.sum() == 18, "9 colliding symbols, 2 Ensembl entries each"
    assert all("__ENSG" in name for name in adata.var.index[disagree])

    path = tmp_path / "roundtrip.h5ad"
    adata.write_h5ad(path)

    import anndata as ad
    back = ad.read_h5ad(path)
    assert back.shape == adata.shape
    assert list(back.var_names) == list(adata.var_names)
    for gene in ("TNFRSF17", "GPRC5D", "NSD2"):
        assert gene in back.var_names


requires_stage05 = pytest.mark.skipif(
    not (config.RESULTS_DIR / "05_integration" / "integrated.h5ad").exists(),
    reason="needs the stage-05 object (run notebooks/05_integration_clustering.ipynb)",
)


@requires_data
@requires_stage05
def test_plasma_cell_depth_gap_is_a_compartment_effect_not_a_cohort_one():
    """The measurement the stage-05 interpretation rests on.

    Harmony mixes the immune compartment (cohort entropy ~0.75) and does not mix the
    plasma-cell compartment (~0.11). The explanation is that WashU was cut at 10,000
    UMIs before deposit and MMRF was not, and plasma cells — professional secretors —
    are the cells that ceiling actually removes. If that is right, the MMRF-vs-WashU
    depth gap must be far larger among plasma cells than among everything else.

    It is: ~4.5x vs ~1.8x. Asserted loosely, since the exact values would move if the
    clustering were ever recomputed, but the asymmetry itself is the claim.
    """
    import anndata

    adata = anndata.read_h5ad(
        config.RESULTS_DIR / "05_integration" / "integrated.h5ad", backed="r"
    )
    mzb1 = np.asarray(adata[:, "MZB1"].X.todense()).ravel() > 0
    clusters = adata.obs["leiden"].astype(str)
    by_cluster = pd.Series(mzb1, index=adata.obs_names).groupby(clusters).mean()
    plasma = set(by_cluster[by_cluster > 0.40].index)
    assert plasma, "no plasma-cell-like cluster found at all"

    frame = adata.obs.assign(
        is_plasma=clusters.isin(plasma).to_numpy(),
        depth=adata.obs["total_counts"].to_numpy(),
    )
    medians = frame.pivot_table(index="cohort", columns="is_plasma", values="depth",
                                aggfunc="median", observed=True)

    immune_gap = medians.loc["MMRF", False] / medians.loc["WU1", False]
    plasma_gap = medians.loc["MMRF", True] / medians.loc["WU1", True]
    assert immune_gap < 2.5, f"immune depth gap unexpectedly large: {immune_gap:.1f}x"
    assert plasma_gap > 3.0, f"plasma depth gap unexpectedly small: {plasma_gap:.1f}x"
    assert plasma_gap > 2 * immune_gap

    # And the mechanism: MMRF plasma cells sit above a ceiling WashU cells cannot cross.
    mmrf_plasma = frame[(frame["cohort"] == "MMRF") & frame["is_plasma"]]
    wu1_plasma = frame[(frame["cohort"] == "WU1") & frame["is_plasma"]]
    assert (mmrf_plasma["depth"] > 10_000).mean() > 0.30
    assert (wu1_plasma["depth"] > 10_000).mean() < 0.05


@requires_data
@requires_stage05
def test_gprc5d_is_far_lower_abundance_than_bcma():
    """Recorded because it is evidence, not because it is a problem.

    GPRC5D fails HVG selection (mean 0.061 vs TNFRSF17's 0.492). That does not affect
    the embedding's job or stage 08 — antigen calls read raw counts — but it is this
    cohort's own confirmation that GPRC5D is a low-abundance transcript, and hence
    that GPRC5D-negative calls carry more technical-zero risk than BCMA-negative ones.
    """
    import anndata

    adata = anndata.read_h5ad(
        config.RESULTS_DIR / "05_integration" / "integrated.h5ad", backed="r"
    )
    means = adata.var["means"]
    assert means["GPRC5D"] < means["TNFRSF17"] / 4
    assert not adata.var.loc["GPRC5D", "highly_variable"]
    assert adata.var.loc["TNFRSF17", "highly_variable"]
