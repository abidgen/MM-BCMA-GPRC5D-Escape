"""
The loader. Name parsing, the GEO metadata join, and the truncation repair need no
count matrices; everything under `requires_data` does.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mm_escape import config, gene_space, io

from conftest import CANONICAL_SAMPLES, requires_data


# ---------------------------------------------------------------------------
# Names
# ---------------------------------------------------------------------------

def test_parse_sample_name_splits_on_the_gsm_accession():
    # Splitting on the first underscore would give ("GSM6939028", "MMRF") here.
    assert io.parse_sample_name("GSM6939028_MMRF_1695") == ("GSM6939028", "MMRF_1695")
    assert io.parse_sample_name("GSM6939047_27522_1") == ("GSM6939047", "27522_1")


def test_parse_sample_name_rejects_a_bare_sample_name():
    with pytest.raises(io.SampleLoadError, match="GSM accession"):
        io.parse_sample_name("MMRF_1695")


@pytest.mark.parametrize(
    "sample_name,expected",
    [
        ("27522_1", "27522"),      # numeric stem -> collapse
        ("27522_6", "27522"),
        ("MMRF_1695", "MMRF_1695"),  # non-numeric stem -> leave alone
        ("MMY83942", "MMY83942"),
        ("ND_083017", "ND_083017"),  # the suffix is a date, not a sample index
        ("BM4", "BM4"),
        ("83942", "83942"),
    ],
)
def test_naive_patient_id(sample_name, expected):
    assert io.naive_patient_id(sample_name) == expected


@pytest.mark.parametrize("name", ["BM2", "BM6", "ND_083017", "ND_170607"])
def test_donor_names_classify_as_normal_marrow(name):
    assert io.classify_sample(name) == ("normal_bm", True)


@pytest.mark.parametrize("name", ["MMRF_1695", "27522_1", "MMY98423"])
def test_patient_names_classify_as_myeloma(name):
    assert io.classify_sample(name) == ("myeloma", True)


# ---------------------------------------------------------------------------
# The committed GEO metadata
# ---------------------------------------------------------------------------

def test_scrna_metadata_covers_the_whole_deposit():
    meta = io.load_sample_metadata("scrna")
    assert len(meta) == 62
    assert not meta["sample_name"].duplicated().any()


def test_the_eight_donors_are_the_samples_with_no_diagnosis():
    # This is what settled ND_* -- GEO gives them source_name "Donor BMMC" and no
    # diagnosis, while the other 54 read "Multiple myeloma (MM)".
    meta = io.load_sample_metadata("scrna")
    donors = meta.loc[meta["sample_type"] == "normal_bm", "sample_name"]
    assert sorted(donors) == [
        "BM2", "BM4", "BM5", "BM6",
        "ND_083017", "ND_090617", "ND_170531", "ND_170607",
    ]
    assert (meta.loc[meta["sample_type"] == "normal_bm", "diagnosis"] == "none").all()
    assert (meta["sample_type"] == "myeloma").sum() == 54


def test_chemistry_follows_cohort_and_only_washu_1_is_v2():
    meta = io.load_sample_metadata("scrna")
    counts = meta.groupby("cohort").size().to_dict()
    assert counts == {"Donor": 8, "MMRF": 18, "WU1": 23, "WU2": 13}
    v2 = meta.loc[meta["chemistry"] == "10x 3' v2", "cohort"].unique()
    assert list(v2) == ["WU1"]
    # WashU 1 is also the only cohort that skipped dead-cell removal.
    no_removal = meta.loc[meta["dead_cell_removal"] == "no", "cohort"].unique()
    assert list(no_removal) == ["WU1"]


def test_cohort_protocol_table_agrees_with_the_parsed_metadata():
    meta = io.load_sample_metadata("scrna")
    for cohort, protocol in io.COHORT_PROTOCOL.items():
        rows = meta.loc[meta["cohort"] == cohort]
        assert (rows["chemistry"] == protocol["chemistry"]).all()
        assert (rows["dead_cell_removal"] == protocol["dead_cell_removal"]).all()


def test_bulk_metadata_splits_sorted_from_unsorted():
    # Stage 09 must not pool these: MMRF bulk is CD138+ sorted and pairs with
    # malignant pseudobulk, WashU 1 bulk is whole marrow and pairs with whole-sample
    # pseudobulk. Pooling makes 10 of 26 comparisons measure tumour burden.
    meta = io.load_sample_metadata("bulk")
    assert len(meta) == 31
    prep = meta.groupby(["cohort", "prep"]).size().to_dict()
    assert prep == {("MMRF", "CD138+ sorted"): 18, ("WU1", "unsorted BMMC"): 13}


def test_bulk_overlap_with_the_single_cell_cohort_is_26():
    # Computed, not the inherited "~28".
    sc = set(io.load_sample_metadata("scrna")["sample_name"])
    bulk = io.load_sample_metadata("bulk")["sample_name"]
    stubs = {"MMRF_1505", "MMRF_2259"}  # 114-byte failed deposits
    matched = [n for n in bulk if n in sc and n not in stubs]
    assert len(matched) == 26
    orphans = sorted(n for n in bulk if n not in sc)
    assert orphans == ["47499", "59114_2", "98433"]


def test_rebuild_metadata_rejects_a_deposit_of_the_wrong_size(tmp_path):
    soft = tmp_path / "tiny.soft"
    soft.write_text("^SAMPLE = GSM1\n!Sample_title = x\n")
    with pytest.raises(io.SampleLoadError, match="expected 62"):
        io.rebuild_sample_metadata_from_soft(soft, "scrna")


# ---------------------------------------------------------------------------
# The truncation repair (needs the committed gene map, not the deposit)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def canonical_33694():
    return list(gene_space.load_gene_map(33694)["deposited_symbol"])


@pytest.fixture(scope="module")
def truncated(canonical_33694):
    """Reproduce the deposited damage: 22184 intact rows plus a partial 22185th."""
    rows = config.TRUNCATED_GENE_FILES["56203_1"]["deposited_rows"]
    return canonical_33694[: rows - 1] + [canonical_33694[rows - 1][:-1]]


def test_repair_restores_the_full_column(canonical_33694, truncated):
    repaired = io._repair_truncated_genes("56203_1", truncated)
    assert repaired == canonical_33694
    assert len(repaired) == 33694
    # The genes that were "missing" were simply past the cut.
    assert "TNFRSF17" not in truncated and "TNFRSF17" in repaired
    for gene in ("IGLC1", "IGLC2", "IGLC3"):
        assert gene in repaired


def test_repair_keeps_every_row_the_deposit_actually_wrote(canonical_33694, truncated):
    repaired = io._repair_truncated_genes("56203_1", truncated)
    assert repaired[: len(truncated) - 1] == truncated[:-1]


def test_repair_refuses_a_file_of_unexpected_length(truncated):
    with pytest.raises(io.SampleLoadError, match="recorded truncation"):
        io._repair_truncated_genes("56203_1", truncated[:-1])


def test_repair_refuses_when_the_written_prefix_diverges(truncated):
    damaged = ["NOTAGENE"] + truncated[1:]
    with pytest.raises(io.SampleLoadError, match="diverges at row 1"):
        io._repair_truncated_genes("56203_1", damaged)


def test_repair_refuses_when_the_last_row_is_not_a_truncated_symbol(truncated):
    with pytest.raises(io.SampleLoadError, match="not a prefix of"):
        io._repair_truncated_genes("56203_1", truncated[:-1] + ["ZZZNOTREAL"])


# ---------------------------------------------------------------------------
# The manifest, against the real deposit
# ---------------------------------------------------------------------------

@requires_data
def test_manifest_covers_all_62_samples_with_nothing_excluded(manifest):
    assert len(manifest) == 62
    assert not manifest["excluded"].any()
    assert "56203_1" in set(manifest["sample_name"])


@requires_data
def test_manifest_paths_are_absolute_and_exist(manifest):
    for column in ("matrix_path", "barcodes_path", "genefeat_path"):
        for value in manifest[column]:
            path = Path(value)
            assert path.is_absolute() and path.exists()


@requires_data
def test_manifest_carries_the_geo_metadata(manifest):
    for column in ("cohort", "chemistry", "dead_cell_removal", "diagnosis"):
        assert column in manifest
    assert manifest["sample_type_certain"].all()
    assert (manifest["sample_type"] == "myeloma").sum() == 54
    assert (manifest["sample_type"] == "normal_bm").sum() == 8


@requires_data
def test_patient_ids_are_flagged_provisional_until_s1(manifest):
    assert (manifest["patient_id_source"] == "naive").all()
    myeloma = manifest.loc[manifest["sample_type"] == "myeloma"]
    # 54 samples -> 43 naive patients, against the paper's 53 / 41.
    assert myeloma["patient_id"].nunique() == 43
    assert (manifest.loc[manifest["sample_name"].str.startswith("27522"),
                         "patient_id"] == "27522").all()


@requires_data
def test_manifest_can_be_loaded_without_the_metadata_join(manifest):
    bare = io.load_manifest(with_metadata=False)
    assert "cohort" not in bare
    assert len(bare) == len(manifest)


# ---------------------------------------------------------------------------
# Reading samples
# ---------------------------------------------------------------------------

@requires_data
@pytest.mark.parametrize("name", CANONICAL_SAMPLES)
def test_sample_is_cells_by_genes_with_integral_counts(samples, name):
    adata = samples[name]
    assert adata.n_obs < adata.n_vars           # cells x genes, never the transpose
    assert adata.n_vars in config.BUILDS
    assert adata.X.dtype == np.float32
    data = adata.X.data
    assert data.min() > 0 and np.all(data == np.floor(data))


@requires_data
@pytest.mark.parametrize("name", CANONICAL_SAMPLES)
def test_var_names_are_the_deposited_column_in_deposited_order(samples, name):
    # attach_ensembl_ids joins positionally, so any reordering here is fatal.
    adata = samples[name]
    if adata.uns["genes_repaired"]:
        pytest.skip("gene axis was repaired; covered by the repair tests")
    written = [line.rstrip("\n") for line in open(adata.uns["source"]["genes"])]
    assert list(adata.var_names) == written


@requires_data
@pytest.mark.parametrize("name", CANONICAL_SAMPLES)
def test_obs_names_are_unique_and_carry_the_sample(samples, name):
    adata = samples[name]
    barcodes = [line.rstrip("\n") for line in open(adata.uns["source"]["barcodes"])]
    assert list(adata.obs["barcode"]) == barcodes
    assert list(adata.obs_names) == [f"{name}_{bc}" for bc in barcodes]
    assert not adata.obs_names.duplicated().any()


@requires_data
def test_the_transpose_is_right_way_round(samples):
    # A transposed object still looks plausible, so this checks values, not shape.
    adata = samples["MMRF_1695"]
    triplets = pd.read_csv(
        adata.uns["source"]["matrix"], sep=r"\s+", skiprows=2, header=None,
        names=["gene", "cell", "count"],
    )
    assert len(triplets) == adata.X.nnz
    assert int(triplets["count"].sum()) == int(adata.X.sum())
    for _, row in triplets.sample(12, random_state=0).iterrows():
        assert adata.X[int(row.cell) - 1, int(row.gene) - 1] == row["count"]


@requires_data
def test_56203_1_loads_repaired_with_its_antigens_restored(samples):
    adata = samples["56203_1"]
    assert adata.uns["genes_repaired"] is True
    assert adata.n_vars == 33694 and adata.n_obs == 1837
    written = [line.rstrip("\n") for line in open(adata.uns["source"]["genes"])]
    assert len(written) == 22185 and written[-1] == "KBTBD"
    assert adata.var_names[22184] == "KBTBD7"
    assert list(adata.var_names[:22184]) == written[:22184]
    assert "TNFRSF17" in adata.var_names


@requires_data
def test_intact_samples_are_not_marked_repaired(samples):
    for name in ("MMRF_1695", "27522_1", "BM4"):
        assert samples[name].uns["genes_repaired"] is False


@requires_data
def test_legacy_symbols_appear_in_the_older_build(samples):
    # 27522_1 is a 33694 sample: it spells NSD2 as WHSC1. This is why the gene-space
    # join cannot be on symbols.
    older = samples["27522_1"].var_names
    assert "WHSC1" in older and "NSD2" not in older
    assert "NSD2" in samples["MMRF_1695"].var_names


@requires_data
@pytest.mark.parametrize("name", CANONICAL_SAMPLES)
def test_cohort_metadata_reaches_every_cell(samples, name):
    obs = samples[name].obs
    for column in ("cohort", "chemistry", "dead_cell_removal", "diagnosis",
                   "patient_id", "n_genes_ref", "sample_type"):
        assert column in obs
        assert obs[column].nunique() == 1


@requires_data
def test_unknown_sample_raises(manifest):
    with pytest.raises(io.SampleLoadError, match="matches no manifest row"):
        io.read_sample("not-a-sample", manifest=manifest)


@requires_data
def test_read_samples_yields_in_manifest_order(manifest):
    names = [a.obs["sample_name"].iloc[0]
             for a in io.read_samples(["BM5", "BM4"], manifest=manifest, verbose=False)]
    assert names == ["BM4", "BM5"]


@requires_data
def test_read_samples_rejects_an_unknown_name(manifest):
    with pytest.raises(io.SampleLoadError, match="Unknown sample"):
        list(io.read_samples(["ghost"], manifest=manifest, verbose=False))
