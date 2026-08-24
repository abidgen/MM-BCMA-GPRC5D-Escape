"""
Project-wide constants: thresholds, gene sets, exclusions.

STATUS: partial scaffold. Only the gene-space constants needed by `gene_space.py`
are defined so far. MARKER_PANEL, STATE_PROGRAMS, ANNOTATION_DECISION, QC thresholds
and the antigen noise floor are added by the stages that derive them (04-08) — see
CLAUDE.md. Do not duplicate any of these lists into a notebook; import from here.

Every value is env-var overridable via `_env()` (same convention as the R build's
lib/00_config.R), so a notebook can override without editing source.
"""

from __future__ import annotations

import os
from pathlib import Path

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------

def _env(name: str, default: str) -> str:
    return os.environ.get(f"MM_{name}", default)


#: Repo root, resolved from this file's location (src/mm_escape/config.py).
REPO_ROOT = Path(__file__).resolve().parents[2]

RAW_DIR = Path(_env("RAW_DIR", str(REPO_ROOT / "raw")))
SAMPLES_DIR = RAW_DIR / "samples"
MANIFEST_CSV = RAW_DIR / "sample_manifest.csv"

RESULTS_DIR = Path(_env("RESULTS_DIR", str(REPO_ROOT / "results")))

#: Committed Ensembl-ID reconstruction (see gene_space.py and CLAUDE.md).
GENE_SPACE_DIR = Path(_env("GENE_SPACE_DIR", str(REPO_ROOT / "resources" / "gene_space")))

#: Per-sample GEO metadata parsed from the *_family.soft.gz files (see io.py). Holds
#: cohort / chemistry / diagnosis, none of which is derivable from filenames.
SAMPLE_METADATA_DIR = Path(
    _env("SAMPLE_METADATA_DIR", str(REPO_ROOT / "resources" / "sample_metadata"))
)

# --------------------------------------------------------------------------
# Sample exclusions
# --------------------------------------------------------------------------

#: No samples are excluded. `56203_1` was excluded until 2026-08-24 on the belief that
#: it came from a 22184-gene reference lacking TNFRSF17; the GEO metadata and the files
#: show that was a misdiagnosis — see TRUNCATED_GENE_FILES below. It is repaired and
#: retained instead.
EXCLUDED_SAMPLES: frozenset[str] = frozenset()

# --------------------------------------------------------------------------
# Damaged deposits
# --------------------------------------------------------------------------

#: `56203_1`'s `genes.tsv` write FAILED PART-WAY: the file holds 22185 rows and ends
#: `KBTBD` with no trailing newline, where the 33694 reference has `KBTBD7` at that
#: position. It is a strict prefix of the standard 33694 list, not a different
#: reference — its `counts.mtx` header reads `33694 1837 2135520`, a normal
#: 33694-build matrix. `TNFRSF17` (canonical row 25539) and `IGLC1/2/3` (rows
#: 32548-32552) were not absent from a reference; they were past the cut.
#:
#: (The "22184 genes" in earlier versions of CLAUDE.md is a `wc -l` artifact — the
#: file has no trailing newline, so `wc -l` undercounts by one. There are 22185 rows.)
#:
#: The repair is provable rather than a guess: substitute the canonical symbols for
#: the declared build, which come from the committed, position-verified gene map, and
#: assert the truncated file is a prefix of them. `io.read_sample` does this and
#: hard-fails if the prefix check does not hold.
TRUNCATED_GENE_FILES: dict[str, dict[str, int]] = {
    "56203_1": {"build": 33694, "deposited_rows": 22185},
}

# --------------------------------------------------------------------------
# Cell Ranger reference builds
# --------------------------------------------------------------------------

#: Gene-row count -> the public reference it identifies. Verified by exact positional
#: reconstruction from the Ensembl GTF (0 mismatches, all rows) — see
#: `gene_space.rebuild_gene_map_from_gtf`.
BUILDS: dict[int, dict[str, str]] = {
    33538: {
        "reference": "refdata-cellranger-GRCh38-3.0.0",
        "ensembl_release": "93",
        "gtf_url": (
            "https://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/"
            "Homo_sapiens.GRCh38.93.gtf.gz"
        ),
    },
    33694: {
        "reference": "refdata-cellranger-GRCh38-1.2.0",
        "ensembl_release": "84",
        "gtf_url": (
            "https://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/"
            "Homo_sapiens.GRCh38.84.gtf.gz"
        ),
    },
}

#: The biotype whitelist 10x's `cellranger mkgtf` applies when building the human
#: reference. Reproducing this exactly is what yields 33538 / 33694 rows.
MKGTF_BIOTYPES = frozenset({
    "protein_coding", "lincRNA", "antisense",
    "IG_LV_gene", "IG_V_gene", "IG_V_pseudogene", "IG_D_gene",
    "IG_J_gene", "IG_J_pseudogene", "IG_C_gene", "IG_C_pseudogene",
    "TR_V_gene", "TR_V_pseudogene", "TR_D_gene",
    "TR_J_gene", "TR_J_pseudogene", "TR_C_gene",
})

#: Canonical symbol vintage for the merged object. Ensembl 93 is the newer of the two
#: builds, so its symbols are the modern HGNC ones (NSD2, not WHSC1).
CANONICAL_BUILD = 33538

# --------------------------------------------------------------------------
# Required genes — assertions, not documentation
# --------------------------------------------------------------------------
#
# These caught the NSD2/WHSC1 symbol drift that manual inspection missed across two
# prior builds of this project. They stay, and they stay loud. A missing required
# gene means "check for a legacy symbol / check the join key", never "biologically
# absent". All 65 are confirmed present in the Ensembl-ID intersection (2026-08-21).

REQUIRED_GENES: dict[str, frozenset[str]] = {
    # BCMA, GPRC5D + the coverage-matrix targets (stage 08)
    "antigens": frozenset({
        "TNFRSF17", "GPRC5D", "SLAMF7", "FCRL5", "SDC1", "CD38", "ITGB7", "NCSTN",
    }),
    # kappa/lambda restriction (stage 07). Ratio-based, so all IGLC members matter.
    "light_chain": frozenset({
        "IGKC", "IGLC1", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7",
    }),
    # seven-class annotation panel (stage 06)
    "markers": frozenset({
        "MZB1", "XBP1", "IRF4",                      # PlasmaCell (+ SDC1, CD38)
        "MS4A1", "CD79A", "CD19",                    # Bcell
        "CD3D", "CD3E", "CD8A", "CD4",               # Tcell
        "NCAM1", "NKG7", "GNLY",                     # NK
        "CD14", "LYZ", "ITGAM",                      # Myeloid
        "HBB", "GYPA",                               # Erythroid
        "CD34", "KIT",                               # HSPC
    }),
    # TC-like expression subtype (stage 10). NOT a translocation call — a proxy.
    # NSD2 is WHSC1 in the 33694 build; without a correct join this class is uncallable.
    "tc": frozenset({
        "CCND1", "CCND2", "CCND3", "NSD2", "FGFR3", "MAF", "MAFB", "CKS1B",
    }),
    # pre-registered gamma-secretase hypothesis (stage 10)
    "gamma_secretase": frozenset({
        "NCSTN", "PSEN1", "APH1A", "APH1B", "PSENEN",
    }),
    # orthogonal cell-state programs (stages 06/10) — continuous scores, never labels
    "programs": frozenset({
        "MKI67", "TOP2A", "PCNA",                    # cell cycle
        "ISG15", "IFI6", "STAT1", "MX1",             # interferon
        "B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-DRA", # antigen presentation
        "ATF4", "HSPA5", "DDIT3",                    # UPR (+ XBP1)
        "MYC",                                       # MYC program (stage 10)
    }),
}

#: Legacy-symbol pairs that a raw-symbol intersection silently drops.
#:
#: DEMOTED to a regression assertion (2026-08-21): the gene-space join is on Ensembl
#: ID, which resolves all 11,140 drifted symbols. This map covers 4 of them and was
#: never the harmonization mechanism. Keep it as a canary — if a canonical symbol here
#: ever goes missing, the join key regressed to symbols somewhere.
#:
#: NSD3/WHSC1L1 is a DIFFERENT GENE from NSD2/WHSC1. Never fuzzy-match these.
LEGACY_SYMBOLS: dict[str, str] = {
    "WHSC1": "NSD2",        # t(4;14) — highest-risk MM translocation
    "FAM46C": "TENT5C",     # recurrently deleted MM tumour suppressor (1p12)
    "WHSC1L1": "NSD3",      # NOT NSD2
    "ATP5A1": "ATP5F1A",    # OXPHOS program member
}


def all_required_genes() -> frozenset[str]:
    """Flatten REQUIRED_GENES into a single set of canonical symbols."""
    return frozenset().union(*REQUIRED_GENES.values())
