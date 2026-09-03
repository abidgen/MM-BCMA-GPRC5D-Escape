"""
io.py -> gene_space.py, on real count matrices.

Renamed from `test_integration.py` on 2026-08-24, when `src/mm_escape/integration.py`
arrived and claimed that name under the one-test-file-per-module convention. "Integration"
here means the software sense (two modules exercised together); in `test_integration.py`
it means the biological sense (Harmony batch correction). Different concerns, and the
old name made them collide.

This is the pairing neither module can verify alone: `gene_space` asserts that
`var_names` are the deposited column position-for-position, and `io` is the only
thing that produces them. Everything here needs the extracted deposit.
"""

from __future__ import annotations

import anndata
import numpy as np
import pytest

from mm_escape import config, gene_space, io

from conftest import CANONICAL_SAMPLES, requires_data

pytestmark = requires_data


@pytest.fixture(scope="module")
def harmonized(samples):
    """The four canonical samples, merged the way stage 05 will merge all 62."""
    objects = [gene_space.attach_ensembl_ids(samples[n].copy()) for n in CANONICAL_SAMPLES]
    objects = gene_space.intersect_gene_space(objects, verbose=False)
    merged = anndata.concat(objects, join="inner")
    return gene_space.to_canonical_symbols(merged)


def test_attach_accepts_what_the_loader_produces(samples):
    for name in CANONICAL_SAMPLES:
        attached = gene_space.attach_ensembl_ids(samples[name].copy())
        assert attached.var_names.str.match(r"^ENSG\d{11}$").all()
        assert "deposited_symbol" in attached.var


def test_attach_rejects_a_reordered_gene_axis(samples):
    # The positional join is only safe because this is impossible to do by accident.
    scrambled = samples["BM4"][:, ::-1].copy()
    with pytest.raises(gene_space.GeneSpaceError, match="position-for-position"):
        gene_space.attach_ensembl_ids(scrambled)


def test_intersection_on_real_data_matches_the_committed_table(harmonized):
    assert harmonized.n_vars == 32991


def test_intersecting_on_symbols_is_refused(samples):
    # Objects straight from the loader are symbol-keyed; intersecting them is the
    # bug the whole gene_space module exists to prevent.
    with pytest.raises(gene_space.GeneSpaceError, match="not keyed on Ensembl IDs"):
        gene_space.intersect_gene_space(
            [samples["MMRF_1695"].copy(), samples["BM4"].copy()], verbose=False
        )


def test_every_required_gene_survives_the_real_merge(harmonized):
    gene_space.assert_required_genes(harmonized)


def test_drifted_symbols_are_joined_not_dropped(harmonized):
    # The 33538 and 33694 builds spell this gene differently; both rows are the same
    # Ensembl entry and must end up as one column.
    assert harmonized.var.loc["NSD2", "symbol_33538"] == "NSD2"
    assert harmonized.var.loc["NSD2", "symbol_33694"] == "WHSC1"
    assert int(harmonized.var["symbol_drift"].sum()) == 11140


def test_cells_from_all_four_samples_survive_with_unique_names(harmonized, samples):
    assert harmonized.obs["sample_name"].nunique() == 4
    assert not harmonized.obs_names.duplicated().any()
    assert harmonized.n_obs == sum(samples[n].n_obs for n in CANONICAL_SAMPLES)


def test_metadata_survives_the_merge(harmonized):
    for column in ("cohort", "chemistry", "patient_id", "sample_type", "n_genes_ref"):
        assert column in harmonized.obs
    # 56203_1 and 27522_1 are both WashU 1; MMRF_1695 is MMRF; BM4 is a donor.
    assert set(harmonized.obs["cohort"].unique()) == {"WU1", "MMRF", "Donor"}


def test_counts_are_preserved_through_the_merge(harmonized, samples):
    # Subsetting to the intersection drops genes, so totals fall -- but no cell may
    # gain counts, and the retained genes must carry exactly what they carried.
    merged_total = float(harmonized.X.sum())
    raw_total = sum(float(samples[n].X.sum()) for n in CANONICAL_SAMPLES)
    assert 0 < merged_total <= raw_total
    one = samples["MMRF_1695"]
    row = one.obs_names[0]
    assert float(harmonized[row, "TNFRSF17"].X.sum()) == float(
        one[row, "TNFRSF17"].X.sum()
    )


@pytest.mark.slow
def test_the_whole_cohort_loads(manifest):
    """62 samples, every one of them, with the totals this deposit is known to have."""
    n_cells = 0
    per_cohort: dict[str, int] = {}
    for adata in io.read_samples(manifest=manifest, verbose=False):
        assert adata.n_vars in config.BUILDS
        assert adata.n_obs > 0
        n_cells += adata.n_obs
        cohort = str(adata.obs["cohort"].iloc[0])
        per_cohort[cohort] = per_cohort.get(cohort, 0) + adata.n_obs

    assert n_cells == 204_040
    assert sorted(per_cohort) == ["Donor", "MMRF", "WU1", "WU2"]
    assert per_cohort["Donor"] == 48_150


@pytest.mark.slow
def test_sequencing_depth_tracks_cohort_not_chemistry_version(manifest):
    """The confounder, asserted so a regression in the metadata join is visible.

    MMRF cells carry ~1.9x the genes of WashU 1's. Note v2-vs-v3 alone is only 1.38x
    with overlapping distributions -- the separation follows cohort, and the project
    plan says not to quote a 2-3x chemistry effect.
    """
    medians: dict[str, list[float]] = {}
    for adata in io.read_samples(manifest=manifest, verbose=False):
        cohort = str(adata.obs["cohort"].iloc[0])
        genes_per_cell = np.diff(adata.X.indptr)
        medians.setdefault(cohort, []).append(float(np.median(genes_per_cell)))

    by_cohort = {k: float(np.median(v)) for k, v in medians.items()}
    assert by_cohort["MMRF"] > by_cohort["WU2"] > by_cohort["Donor"] > by_cohort["WU1"]
    assert 1.7 < by_cohort["MMRF"] / by_cohort["WU1"] < 2.1
