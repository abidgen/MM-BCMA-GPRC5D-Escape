"""
Gene-space reconstruction. Runs entirely on committed resources — no count matrices.

These cover the single most consequential correctness property in the project: that
the cross-reference join is on Ensembl ID and not on symbols. A symbol join does not
merely lose genes, it can pair the WRONG gene (TBCE is a different annotation entry
in each build), and nothing downstream would look wrong.
"""

from __future__ import annotations

import pandas as pd
import pytest

from mm_escape import config, gene_space


# ---------------------------------------------------------------------------
# The two Seurat transforms
# ---------------------------------------------------------------------------

def test_make_unique_matches_r_semantics():
    # First occurrence keeps the bare name; later ones take .1, .2, ...
    assert gene_space._make_unique(["A", "A", "B", "A"]) == ["A", "A.1", "B", "A.2"]


def test_make_unique_skips_a_candidate_that_is_already_taken():
    # R's make.unique will not hand out "A.1" twice. This is the case a naive
    # counter-based implementation gets wrong, and it changes row identity.
    assert gene_space._make_unique(["A", "A.1", "A"]) == ["A", "A.1", "A.2"]


def test_seurat_sanitize_rewrites_underscores():
    # Ensembl writes RP11-442N24__B; the deposit spells it RP11-442N24--B.
    assert gene_space._seurat_sanitize("RP11-442N24__B") == "RP11-442N24--B"
    assert gene_space._seurat_sanitize("NSD2") == "NSD2"


# ---------------------------------------------------------------------------
# The committed map
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("build", sorted(config.BUILDS))
def test_gene_map_has_exactly_the_reference_row_count(build):
    frame = gene_space.load_gene_map(build)
    assert len(frame) == build
    assert list(frame.columns) == ["deposited_symbol", "ensembl_id"]
    assert frame["ensembl_id"].str.match(r"^ENSG\d{11}$").all()
    assert not frame["ensembl_id"].duplicated().any()


def test_tbce_proves_the_symbol_join_would_mis_pair():
    # The whole argument for joining on IDs, in one assertion.
    ids = {}
    for build in (33538, 33694):
        frame = gene_space.load_gene_map(build)
        ids[build] = frame.loc[frame["deposited_symbol"] == "TBCE", "ensembl_id"].iloc[0]
    assert ids[33538] != ids[33694], "TBCE must resolve to different Ensembl entries"


def test_unknown_build_raises_and_a_truncated_row_count_gets_a_hint():
    with pytest.raises(gene_space.GeneSpaceError, match="no known reference build"):
        gene_space.detect_build(["x"] * 12345)
    with pytest.raises(gene_space.GeneSpaceError, match="TRUNCATED"):
        gene_space.detect_build(["x"] * 22185)


@pytest.mark.parametrize("build", sorted(config.BUILDS))
def test_detect_build_recognises_both_references(build):
    assert gene_space.detect_build(["x"] * build) == build


# ---------------------------------------------------------------------------
# The intersection table
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def intersection():
    path = config.GENE_SPACE_DIR / "gene_space_intersection.tsv.gz"
    return pd.read_csv(path, sep="\t")


def test_intersection_recovers_the_full_gene_space(intersection):
    # 32,991 vs 22,164 on raw symbols. A regression here means the join key moved.
    assert len(intersection) == 32991
    assert not intersection["ensembl_id"].duplicated().any()


def test_drifted_symbols_are_the_reason_the_id_join_matters(intersection):
    assert int(intersection["symbol_drift"].sum()) == 11140
    nsd2 = intersection.set_index("canonical_symbol").loc["NSD2"]
    assert nsd2["symbol_33538"] == "NSD2" and nsd2["symbol_33694"] == "WHSC1"


def test_every_required_gene_survives_the_intersection(intersection):
    present = set(intersection["canonical_symbol"])
    missing = sorted(config.all_required_genes() - present)
    assert not missing, f"required genes absent from the intersection: {missing}"


def test_nine_symbols_still_collide_and_are_suffixed_not_dropped(intersection):
    collisions = intersection["canonical_symbol"].duplicated(keep=False)
    assert sorted(set(intersection.loc[collisions, "canonical_symbol"])) == [
        "COG8", "CYB561D2", "EMG1", "LINC01238", "LINC01505",
        "MATR3", "PINX1", "RGS5", "TMSB15B",
    ]


# ---------------------------------------------------------------------------
# assert_required_genes
# ---------------------------------------------------------------------------

class _FakeVar:
    """Minimal stand-in for an AnnData, so these run without building matrices."""

    def __init__(self, symbols):
        self.var = pd.DataFrame(index=pd.Index(symbols))
        self.var_names = self.var.index


def test_assert_required_genes_passes_on_a_complete_gene_set(intersection):
    gene_space.assert_required_genes(_FakeVar(list(intersection["canonical_symbol"])))


def test_missing_gene_is_named_and_the_legacy_symbol_is_hinted(intersection):
    symbols = [s for s in intersection["canonical_symbol"] if s != "NSD2"]
    with pytest.raises(gene_space.GeneSpaceError) as excinfo:
        gene_space.assert_required_genes(_FakeVar(symbols))
    message = str(excinfo.value)
    assert "NSD2" in message and "WHSC1" in message


def test_a_missing_antigen_fails_loudly(intersection):
    symbols = [s for s in intersection["canonical_symbol"] if s != "TNFRSF17"]
    with pytest.raises(gene_space.GeneSpaceError, match="TNFRSF17"):
        gene_space.assert_required_genes(_FakeVar(symbols))
