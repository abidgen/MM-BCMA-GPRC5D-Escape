"""Constants that other modules assert against, so they get asserted on themselves."""

from __future__ import annotations

from mm_escape import config


def test_required_genes_groups_are_populated():
    assert config.REQUIRED_GENES
    for group, genes in config.REQUIRED_GENES.items():
        assert genes, f"required-gene group {group!r} is empty"


def test_all_required_genes_flattens_without_loss():
    flat = config.all_required_genes()
    assert flat == frozenset().union(*config.REQUIRED_GENES.values())
    # The two antigens the whole project is about.
    assert {"TNFRSF17", "GPRC5D"} <= flat


def test_nothing_is_excluded_and_56203_1_is_repaired_instead():
    # Reversed 2026-08-24: 56203_1 was excluded on the belief that it came from a
    # 22184-gene reference lacking TNFRSF17. It is a truncated 33694 file.
    assert config.EXCLUDED_SAMPLES == frozenset()
    assert "56203_1" in config.TRUNCATED_GENE_FILES


def test_truncation_record_matches_what_was_verified():
    spec = config.TRUNCATED_GENE_FILES["56203_1"]
    assert spec["build"] == 33694
    # 22185, not 22184 — the file has no trailing newline, so `wc -l` undercounts.
    assert spec["deposited_rows"] == 22185
    assert spec["build"] in config.BUILDS


def test_legacy_symbol_map_never_conflates_nsd2_with_nsd3():
    # NSD3/WHSC1L1 is a DIFFERENT GENE from NSD2/WHSC1. Fuzzy-matching these would
    # silently swap a t(4;14) readout for an unrelated paralog.
    assert config.LEGACY_SYMBOLS["WHSC1"] == "NSD2"
    assert config.LEGACY_SYMBOLS["WHSC1L1"] == "NSD3"
    assert len(set(config.LEGACY_SYMBOLS.values())) == len(config.LEGACY_SYMBOLS)
    assert not (set(config.LEGACY_SYMBOLS) & set(config.LEGACY_SYMBOLS.values()))


def test_canonical_build_is_the_newer_symbol_vintage():
    assert config.CANONICAL_BUILD == 33538
    assert config.BUILDS[33538]["ensembl_release"] == "93"
    assert config.BUILDS[33694]["ensembl_release"] == "84"


def test_committed_resources_are_present():
    assert (config.GENE_SPACE_DIR / "gene_space_intersection.tsv.gz").exists()
    for build in config.BUILDS:
        assert (config.GENE_SPACE_DIR / f"genes_{build}_ensembl.tsv.gz").exists()
    for assay in ("scrna", "bulk"):
        assert (config.SAMPLE_METADATA_DIR / f"{assay}_samples.tsv").exists()
