"""
Project-wide constants: thresholds, gene sets, exclusions.

STATUS: partial scaffold. Only the gene-space constants needed by `gene_space.py`
are defined so far. QC thresholds
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
# Supplementary Table S1 (clinical)
# --------------------------------------------------------------------------

#: The paper's Supplementary Table S1 (Cancer Research CAN-22-1769). Landed in the
#: repo 2026-08-24 and closes the patient mapping that had blocked stage 08 — see
#: `io.rebuild_clinical_metadata_from_s1`. Sheet 1 is per-patient clinical
#: characteristics; sheet 2 is the per-sample disease stage for WashU cohort 1, which
#: is what proves the `_N` suffixes are serial timepoints.
S1_XLSX = RAW_DIR / "can-22-1769_table_s1_suppst1.xlsx"

#: `83942` (WashU cohort 1) and `MMY83942` (WashU cohort 2) are listed by S1 as two
#: patients but carry IDENTICAL age (63), gender (Male), race (White), ISS stage (3)
#: and treatment (Unknown), and the shared numeric stem is not a coincidence. Treating
#: them as one patient sampled under both cohort protocols is what makes the deposit
#: reproduce the paper's own arithmetic exactly: 54 deposited MM samples minus
#: `25183` (below) is 53, and 43 naive patients minus `25183` minus this collapse is
#: 41 — the "53 bone marrow aspirates from 41 MM patients" the series summary states.
#: Maps the S1/deposit name -> the canonical patient id.
PATIENT_ALIASES: dict[str, str] = {"MMY83942": "83942"}

#: `25183` is deposited (scRNA *and* bulk, WashU cohort 1) but appears in NO
#: supplementary table — not the clinical summary, not the disease-stage sheet. It is
#: the sample the paper's 53-vs-54 gap is made of. It is NOT dropped here: the data is
#: real and stage 07's malignant caller can use it. It carries `clinical_source ==
#: "none"` and `in_paper_cohort == False` so any per-patient aggregate can exclude it
#: deliberately rather than by accident.
SAMPLES_WITHOUT_CLINICAL: frozenset[str] = frozenset({"25183"})

#: What the deposit must reproduce once S1 is applied. Asserted by the parser so a
#: revised S1 or a changed deposit fails loudly instead of quietly moving the
#: denominator of the headline metric.
N_MYELOMA_SAMPLES_DEPOSITED = 54
N_MYELOMA_SAMPLES_IN_PAPER = 53
N_PATIENTS_IN_PAPER = 41

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


# --------------------------------------------------------------------------
# Stage 06 — annotation
# --------------------------------------------------------------------------

#: The seven project classes and the markers that define them (Method A, manual).
#:
#: Scored at CLUSTER level, never per cell: this cohort's median cell has 1,162
#: detected genes (172,940 post-QC cells; 1,521 in the plasma compartment), and
#: clustering is what absorbs that dropout. A per-cell marker call on a gene that
#: dropped out is a wrong call, not a missing one.
#:
#: These seven are the ONLY load-bearing identity labels in the project. Everything
#: downstream reads `obs["cell_type"]`, which takes exactly these values (plus
#: `Ambiguous`). Fine subtypes live in `cell_type_fine` and are never load-bearing.
#: REVISED FOR v2 (2026-08-25) — see the stage 06 v1->v2 entry in CLAUDE.md. The v1
#: panel produced NK = 33,556 against Tcell = 19,133 in bone marrow, which is not
#: credible. Two definitional defects, both fixed here from lineage biology rather
#: than by tuning against the misassigned clusters:
#:
#:   Erythroid was ("HBB", "GYPA"). Stage 04 established that haemoglobin is the
#:   dominant ambient species in this marrow (~32% of counts in the flagged decile),
#:   and HBB is in fact detected in 61-85% of cells of EVERY class. A two-gene panel
#:   half-driven by ambient made Erythroid the runner-up class in 18 of 30 clusters.
#:   It is now a six-gene erythroid program led by lineage-specific evidence; the
#:   globins stay as supporting evidence but can no longer carry an assignment alone.
#:
#:   Tcell was ("CD3D", "CD3E", "CD8A", "CD4") and NK was ("NCAM1", "NKG7", "GNLY").
#:   NKG7 and GNLY are a CYTOTOXIC GRANULE program shared by NK cells and cytotoxic
#:   T cells, so they cannot separate the two; CD4 is also on monocytes and CD8A on
#:   subsets of NK/DC. T is now the CD3/TCR complex, which is definitional, and NK
#:   adds KLRF1 (NKp80), which unlike NKG7/GNLY is genuinely NK-restricted.
MARKER_PANEL: dict[str, tuple[str, ...]] = {
    "PlasmaCell": ("SDC1", "CD38", "MZB1", "XBP1", "IRF4"),
    "Bcell":      ("MS4A1", "CD79A", "CD19"),
    "Tcell":      ("CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2"),
    "NK":         ("NCAM1", "NKG7", "GNLY", "KLRD1", "KLRF1"),
    "Myeloid":    ("CD14", "LYZ", "ITGAM"),
    "Erythroid":  ("GYPA", "AHSP", "ALAS2", "CA1", "HBA1", "HBA2"),
    "HSPC":       ("CD34", "KIT"),
}

#: The label used when no class wins cleanly. Recorded as ambiguous rather than forced
#: into a class — CLAUDE.md is explicit that forcing is the error to avoid.
AMBIGUOUS_LABEL = "Ambiguous"

#: Per-class concordance bars, **declared before looking at any result**. Pre-declaring
#: is the whole point: otherwise "choose the best method" becomes post-hoc
#: rationalisation of whichever output looks tidier.
#:
#: The tiers follow what each class actually feeds, not how interesting it is:
#:   PlasmaCell  0.95  sets the denominator of frac_double_negative (stage 07)
#:   T/NK/Myeloid 0.90 defines stage 08's ambient noise floor — a plasma cell leaking
#:                     into the "confidently antigen-negative" population inflates the
#:                     floor and biases every antigen call
#:   the rest    0.85  nothing downstream is load-bearing on these
#:
#: F1 here measures **concordance, not accuracy**. The manual labels are a third
#: opinion from the same matrix, not ground truth. The biological evidence is the
#: marker-coverage test, which can VETO a class regardless of concordance.
CONCORDANCE_THRESHOLDS: dict[str, float] = {
    "PlasmaCell": 0.95,
    "Tcell": 0.90,
    "NK": 0.90,
    "Myeloid": 0.90,
    "Bcell": 0.85,
    "Erythroid": 0.85,
    "HSPC": 0.85,
}

#: Minimum mean scaled expression of a class's own markers, within cells that a method
#: assigned to that class, for the marker-coverage test to pass. Declared in advance
#: alongside the concordance bars.
MARKER_COVERAGE_MIN = 0.30

#: Orthogonal cell-state programs — **continuous scores, never identity labels**.
#:
#: Identity and state are different axes: a cell has one `cell_type` but can carry
#: several active programs at once. A cycling plasma cell is `PlasmaCell` PLUS a high
#: cell-cycle score, not a "Cycling" cell type. If any annotation method emits a
#: proliferation label as an identity, it is remapped to PlasmaCell + score.
#:
#: Genes absent from the harmonized space are dropped with a warning rather than
#: failing the run — unlike REQUIRED_GENES, these panels are indicative, not a
#: contract. Hypoxia in particular is a cheap standard confounder whose exact panel
#: nothing downstream depends on.
STATE_PROGRAMS: dict[str, tuple[str, ...]] = {
    # A proliferative escape subclone is a different clinical risk from a quiescent
    # one — feeds stage 10.
    "cell_cycle": ("MKI67", "TOP2A", "PCNA", "CCNB1", "CDK1"),
    # Immune-pressure marker; feeds stage 11's evasion question.
    "interferon": ("ISG15", "IFI6", "STAT1", "MX1", "IFI44L", "OAS1"),
    # B2M loss is a documented immune-escape route in myeloma. CAR-T is
    # MHC-independent so this does NOT affect the escape metric — it is a *competing*
    # evasion mechanism and belongs in the stage 11/12 interpretation.
    "antigen_presentation": ("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1"),
    # Plasma cells are professional secretors; UPR tone is core plasma-cell biology.
    "upr": ("XBP1", "ATF4", "HSPA5", "DDIT3", "EDEM1", "HERPUD1"),
    # Standard confounder — cheap to score, expensive to discover late.
    "hypoxia": ("VEGFA", "SLC2A1", "PGK1", "LDHA", "NDRG1"),
    # Stage 10 asks whether the escape population is MYC-high; scored here so the
    # program panels live in one place.
    "myc": ("MYC", "NPM1", "NCL", "RPL3"),
    # Stage 10 Level-2. A standard axis of malignant plasma-cell heterogeneity and a
    # common covariate of proliferation. NOTE: broad activity/ribosome/ETC-type sets are
    # known to track library depth as a pure technical property, so stage 10 reports the
    # score-vs-depth correlation BEFORE any DN-vs-comparator difference.
    "oxphos": ("NDUFA4", "NDUFB2", "NDUFS5", "COX5A", "COX6C", "COX7C", "UQCRB",
               "UQCRQ", "ATP5F1E", "ATP5MC2", "SDHB"),
    # Stage 10 Level-2 stress tone: heat-shock plus immediate-early.
    "stress": ("HSPA1A", "HSPA1B", "HSPB1", "DNAJB1", "JUN", "FOS", "EGR1"),
    # PRE-REGISTERED before any stage-10 biology was inspected (see CLAUDE.md stage 10).
    # gamma-secretase cleaves BCMA off the cell surface, so a gamma-secretase-high escape
    # phenotype would be directly actionable. EXACTLY these five genes: no gene may be
    # added after seeing a result, and no single member may carry the claim alone.
    "gamma_secretase": ("NCSTN", "PSEN1", "APH1A", "APH1B", "PSENEN"),
}

#: Stage 10 — TC (Translocation/Cyclin D) molecular subgroup, assigned PER PATIENT from
#: pseudobulk over malignant cells, DESCRIPTIVE ONLY. These are the genes whose
#: dysregulation defines the founder event; overexpression is *consistent with* the
#: translocation and is never a breakpoint call. Every output label reads "TC-like
#: expression subtype", never "patient has t(4;14)". S1 carries no cytogenetics, so there
#: is nothing in this deposit to validate the proxy against and it stays a proxy.
#: `NSD2` is `WHSC1` in the older reference — this depends on stage 05's symbol
#: harmonization; without it the t(4;14)-like class cannot be called at all.
TC_GENES: dict[str, tuple[str, ...]] = {
    "TC_11_14_like": ("CCND1",),
    "TC_6_14_like": ("CCND3",),
    "TC_4_14_like": ("NSD2", "FGFR3"),
    "TC_14_16_like": ("MAF",),
    "TC_14_20_like": ("MAFB",),
    "TC_D2_like": ("CCND2",),
}

#: 1q21-gain readout, reported beside the TC class (also cross-checks a CNV call on that
#: arm, which this project cannot make — stage 07 CNV is NOT_EVALUABLE).
TC_1Q21_GENE: str = "CKS1B"

#: Stage 10 Level-2 primary program set, FROZEN. Level-2 evidence may come only from
#: these; no program is added after seeing a result.
LEVEL2_PROGRAMS: tuple[str, ...] = (
    "myc", "oxphos", "stress", "interferon", "upr", "antigen_presentation",
    "gamma_secretase",
)

# --------------------------------------------------------------------------
# Stage 06 — structured lineage support (Checkpoint 2, 2026-08-25, BEFORE any result)
# --------------------------------------------------------------------------
#
# WHY A FLAT PANEL WAS NOT ENOUGH
# -------------------------------
# `Myeloid` was ("CD14", "LYZ", "ITGAM"). At this cohort's depth CD14 sits at 21-24%
# in unambiguous monocyte clusters — just under MANUAL_MARKER_DETECT_MIN — and ITGAM
# is a surface protein whose transcript is poorly captured, so only LYZ cleared and
# 6,533 monocytes failed a 2-of-3 support test. That is the same defect as the v1
# two-gene Erythroid panel: too few genes, too ambient-dependent.
#
# The fix is NOT a bigger flat list. `Myeloid` is an umbrella over biologically
# distinct programs, and "half the genes from one combined list" is not how any of
# them is recognised. Support is therefore per-SUBPROGRAM, and any one coherent
# subprogram is sufficient.

#: ONTOLOGY NOTE — pDC sits under broad `Myeloid` in this project, and this is a
#: deliberate, pre-existing decision rather than an accident of marker membership:
#: `annotation.CELLTYPIST_TO_BROAD` has mapped "pDC" -> "Myeloid" since stage 06 v1,
#: and `benchmark.py` does the same. What changes here is only that the MANUAL
#: reference stops being blind to pDC. CLAUDE.md previously described pDC as a
#: population "the seven-class panel does not cover" — true of the old flat panel,
#: no longer true of this one, and updated there accordingly.
#: Broad Myeloid has THREE INDEPENDENT ROUTES. Any one coherent route suffices; there
#: is no combined list, because "half the genes from one flat panel" is not how any of
#: these programs is recognised.
#:
#: The conventional-DC route uses ANCHOR + CONTEXT, and the split is a lineage-biology
#: statement, not a tuning device:
#:
#:   ANCHORS are myeloid/DC-RESTRICTED. They are not part of lymphoid biology, so their
#:   presence is genuine lineage evidence:
#:     FCER1G  - Fc-receptor gamma chain, the signalling adaptor of myeloid Fc
#:               receptors; absent from B and T lymphocytes.
#:     TYROBP  - DAP12, the myeloid/NK signalling adaptor (TREM2/SIRP-beta partner);
#:               not B-lineage.
#:     LST1    - myeloid-restricted, MHC class III locus, monocyte/DC expression.
#:     AIF1    - Iba1, monocyte/macrophage-restricted.
#:     FCER1A  - high-affinity IgE receptor alpha, cDC2 and basophils.
#:     CD1C    - cDC2 lineage marker.
#:
#:   CONTEXT is SHARED PROFESSIONAL-APC biology. B cells, activated T cells and
#:   progenitors all present antigen, so these genes corroborate a DC call but can
#:   never establish myeloid identity:
#:     HLA-DRA/DRB1/DPA1/DPB1, CD74 - MHC class II and its invariant chain.
#:     CST3, CTSS - cystatin C and cathepsin S. CTSS in particular is used by B cells
#:               for invariant-chain processing, so it is NOT myeloid-restricted.
#:
#: The evidence for putting MHC-II on the context axis is direct: Leiden cluster 3, a
#: T-cell cluster, scored 6/7 on the old flat DC list (0.857) purely on these shared
#: genes, which cost it its manual Tcell call. B-cell clusters 10 and 21 scored 0.86 the
#: same way. A gene that scores that high in T and B cells is not identifying myeloid.
#: C2d — cDC anchors are cDC-SPECIFIC, not broad innate.
#:
#: C2b used FCER1G/TYROBP/LST1/AIF1 as DC anchors. The C2c monocyte redesign then
#: classified those same four genes as broad innate CONTEXT, which is a contradiction:
#: a gene cannot be lineage-restricted evidence in one route and generic machinery in
#: another. The cross-route audit found exactly these four conflicts and no others.
#: They are demoted to context here; the anchors become genes that mark conventional
#: DC identity itself.
#:
#:   FCER1A   - high-affinity IgE receptor alpha. Among marrow mononuclear cells its
#:              expression is essentially confined to cDC2 (and basophils), unlike the
#:              Fc-receptor ADAPTOR FCER1G which every innate cell uses.
#:   CD1C     - BDCA-1, the defining surface marker of the cDC2 subset.
#:   CLEC10A  - CD301/MGL, a cDC2-restricted C-type lectin.
#:   CD1E     - CD1 family lipid-antigen presentation, restricted to conventional DC.
#:
#: Four genes, so the 0.5 axis needs TWO independent anchors: no single accidental or
#: shared marker can carry a cDC call. The threshold is untouched; the anchor set was
#: sized so the threshold means something.
#:
#: STATED LIMITATION: this set is cDC2-weighted. cDC1 (CLEC9A/XCR1/BATF3/CADM1) is rare
#: in bone marrow, and mixing cDC1 and cDC2 anchors into one set would leave NEITHER
#: subset able to reach 0.5 on its own — a pure cDC1 cluster would score 4/9. So cDC1
#: is under-called by design rather than by accident. Under-calling is the safe
#: direction: it costs an Ambiguous, not a false Myeloid.
#:
#: FLT3 was considered and REJECTED as an anchor: it is expressed on hematopoietic
#: progenitors, which would reintroduce exactly the progenitor false-positive this
#: revision exists to remove.
MYELOID_DC_ANCHORS: tuple[str, ...] = ("FCER1A", "CD1C", "CLEC10A", "CD1E")

#: Context: shared professional-APC and broad innate machinery. None of it may
#: establish cDC identity. FCER1G/TYROBP/LST1/AIF1 join the MHC-II genes, CD74, CST3
#: and CTSS here, which makes their role identical to their role in the monocyte
#: route and resolves the cross-route contradiction.
MYELOID_DC_CONTEXT: tuple[str, ...] = (
    "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "CST3", "CTSS",
    "FCER1G", "TYROBP", "LST1", "AIF1",
)

#: pDC is its own route and must NOT reach Myeloid through the MHC-II/DC path — pDCs
#: are MHC-II-low and their identity rests on a distinct transcription-factor and
#: receptor program.
#:   TCF4 (E2-2) - the master pDC transcription factor.
#:   LILRA4 (ILT7), CLEC4C (BDCA-2), IL3RA (CD123) - pDC-restricted surface receptors.
#:   PLD4 - endolysosomal nuclease, strongly pDC-enriched.
#: IRF7/IRF8/GZMB are CONTEXTUAL: IRF7/8 are interferon-axis factors shared with other
#: myeloid cells, and GZMB is a shared cytotoxic-granule gene (the same trap as
#: NKG7/GNLY for NK), so none of them may carry a pDC call alone.
MYELOID_PDC_CORE: tuple[str, ...] = ("TCF4", "LILRA4", "CLEC4C", "IL3RA", "PLD4")

#: The monocyte route uses the same ANCHOR + CONTEXT architecture as conventional DC,
#: for the same reason. The Checkpoint-2b flat list reached 0.75 in Leiden 24 — a
#: cluster carrying genuine CD34/SPINK2/SOX4 progenitor evidence — because six of its
#: eight genes are broad innate machinery that a myeloid-primed progenitor also runs.
#: Broad innate machinery is not committed-monocyte identity.
#:
#:   ANCHORS — committed monocyte/macrophage differentiation. Each is a product of, or
#:   receptor for, the monocyte differentiation programme itself, not of general
#:   hematopoietic or innate biology:
#:     CSF1R    - the M-CSF receptor. Signalling through CSF1R is what COMMITS a
#:                progenitor to the monocyte/macrophage lineage; it is the defining
#:                receptor of that lineage and is not expressed by uncommitted HSPC.
#:     FCN1     - M-ficolin, a secreted lectin that is a hallmark product of mature
#:                classical monocytes; a differentiated secretory product, not
#:                progenitor machinery.
#:     VCAN     - versican, the classical-monocyte proteoglycan, expressed on exit
#:                from the marrow rather than during progenitor stages.
#:     MS4A7    - membrane-spanning 4A family member restricted to non-classical
#:                monocytes and macrophages.
#:     SERPINA1 - alpha-1-antitrypsin, a monocyte/macrophage secretory product.
#:     CD300E   - monocyte-restricted immunoreceptor, absent from lymphoid and
#:                progenitor compartments.
#:     FCGR3A   - CD16, the non-classical monocyte receptor. Also carried by NK cells,
#:                which is acceptable here because NK is adjudicated by the frozen
#:                T/NK logic and does not compete for the Myeloid hypothesis.
#:
#:   CONTEXT — broad innate/hematopoietic machinery. Useful corroboration, never
#:   identity, because each is shared with progenitors and/or other lineages:
#:     FCER1G, TYROBP - innate Fc-receptor and DAP12 signalling adaptors, shared with
#:                NK cells and expressed during myeloid priming before commitment.
#:     LST1     - MHC class III locus, broadly innate.
#:     AIF1     - expressed across the myelomonocytic range including immature stages.
#:     CTSS     - lysosomal protease shared with B cells (invariant-chain processing).
#:     LYZ      - the dominant ambient species in this marrow; supporting only, and
#:                never sufficient, exactly as at Checkpoint 2.
#:
#: NOTE ON DERIVATION: these placements come from monocyte differentiation biology.
#: That FCN1/FCGR3A separate Leiden 24 from Leiden 25 in this dataset is VALIDATION
#: evidence, checked after the split was fixed — not the reason for the split.
#: LILRB1 was in the C2b list and is DROPPED from the anchor set: ILT2 is carried by
#: B cells, NK cells and T-cell subsets, so it is not committed-monocyte evidence.
MYELOID_MONO_ANCHORS: tuple[str, ...] = (
    "CSF1R", "FCN1", "VCAN", "MS4A7", "SERPINA1", "CD300E", "FCGR3A",
)
MYELOID_MONO_CONTEXT: tuple[str, ...] = (
    "FCER1G", "TYROBP", "LST1", "AIF1", "CTSS", "LYZ",
)

#: Empty: every Myeloid route is now anchor+context or an independent core, so there
#: is no flat-list route left. Kept for API stability.
MYELOID_SUBPROGRAMS: dict[str, tuple[str, ...]] = {}

#: Supporting evidence only. LYZ is among the most ambient-abundant transcripts in
#: this marrow, so it may corroborate a subprogram but can never establish Myeloid
#: identity on its own — the v1 Erythroid/HBB lesson applied to the other big ambient
#: species. Likewise GZMB is contextual for pDC and is NOT generic myeloid evidence:
#: it is a shared cytotoxic-granule gene, exactly like NKG7/GNLY for NK.
MYELOID_SUPPORTING: tuple[str, ...] = ("LYZ", "CD14", "ITGAM")
PDC_CONTEXTUAL: tuple[str, ...] = ("IRF7", "IRF8", "GZMB")

#: Primitive progenitor identity, kept separate from lineage priming.
HSPC_CORE: tuple[str, ...] = ("CD34", "HLF", "SPINK2", "GATA2", "MEIS1", "SOX4")

#: Priming context. A progenitor may be lymphoid-primed (IGLL1) or myeloid-primed
#: (MPO) or neither; NEITHER is required and they are never required together.
#: Recorded for interpretation, not used for support.
HSPC_CONTEXT: dict[str, tuple[str, ...]] = {
    "lymphoid_primed": ("IGLL1",),
    "myeloid_primed": ("MPO",),
}

#: THE MATURE-PLASMA PREDICATE, declared before Checkpoint 2 was run.
#:
#: Cluster 24 exposed the problem: a lymphoid progenitor expressing MZB1/XBP1 was
#: pulled toward PlasmaCell by secretory genes alone. A progenitor can run a secretory
#: program without being a mature plasma cell, so a PlasmaCell claim must show
#: CONCORDANT evidence on TWO INDEPENDENT AXES:
#:
#:   (a) secretory/plasma program  — ALL of PLASMA_SECRETORY
#:   (b) mature plasma identity    — at least one of PLASMA_MATURE
#:
#: Deliberately excluded from axis (b): CD38, which is expressed on activated T cells,
#: NK cells, pro-B cells and progenitors and is not plasma-specific. It remains in
#: MARKER_PANEL for descriptive scoring and cannot satisfy this predicate alone.
#: MZB1/XBP1 alone — axis (a) without (b) — also cannot satisfy it.
PLASMA_SECRETORY: tuple[str, ...] = ("MZB1", "XBP1")
PLASMA_MATURE: tuple[str, ...] = ("SDC1", "TNFRSF17")

# --------------------------------------------------------------------------
# Stage 06 — manual adjudication (v3, 2026-08-25, BEFORE any v3 result)
# --------------------------------------------------------------------------
#
# WHY THE MANUAL CLASSIFIER WAS REPLACED
# --------------------------------------
# v1 and v2 assigned each cluster by argmax over `scanpy.tl.score_genes` outputs from
# different panels. Those scores are NOT comparable across panels: score_genes bins all
# genes by mean expression and subtracts a control set drawn from each gene's own bin,
# so a panel's score carries a baseline offset determined by where its genes sit inside
# their bins. Measured on this dataset:
#
#     T panel  : mean(gene - its bin's control mean) = -0.2364   (all 6 genes in bin 24)
#     NK panel : mean(gene - its bin's control mean) = -0.0328   (spread over bins 18-24)
#     systematic offset favouring NK, before any biology  =  0.2036
#
# Clusters 3 and 12 were called NK over T by margins of 0.068 and 0.196 — both SMALLER
# than that offset. The top expression bin is wide enough that its control mean (0.813)
# exceeds every T gene, penalising the entire T panel at once, while the NK panel escapes
# by being distributed across four bins.
#
# The replacement adjudicates on DETECTION FRACTIONS, which are on a common [0, 1] scale
# with no per-panel normalisation, so a cross-panel comparison is meaningful. score_genes
# remains available as a descriptive within-program quantity and no longer decides
# identity.
#
#     manual annotation = positive lineage evidence + specificity/exclusion evidence
#     NOT                = largest independently normalised module score

#: A canonical marker counts as EXPRESSED in a cluster when detected in at least this
#: fraction of the cluster's cells.
#:
#: Set from the two competing error sources, not from any cluster's value. Ambient
#: contamination yields detection roughly proportional to a transcript's share of the
#: droplet pool, which for most genes is in the low single digits; dropout at this
#: cohort's 1,162 median genes/cell puts a genuine marker in a positive population
#: somewhere in the tens of percent. 0.25 sits between the two.
MANUAL_MARKER_DETECT_MIN = 0.25

#: A class is SUPPORTED in a cluster when at least this fraction of its MARKER_PANEL
#: genes are expressed. A majority, because single markers are shared across lineages
#: (NKG7 and GNLY are a cytotoxic-granule program carried by NK *and* cytotoxic T
#: cells) — the requirement is the program, not a gene.
MANUAL_POSITIVE_MIN = 0.5

#: When two classes both survive support and exclusion, the leader must beat the
#: runner-up by this much positive-evidence fraction to be assigned; otherwise the
#: cluster is Ambiguous.
#:
#: Panels hold 3-6 genes, so one marker is worth 0.17-0.33 of the fraction. 0.15 means
#: the leader must lead by roughly a whole marker rather than a rounding difference —
#: which is exactly the distinction the v1/v2 argmax could not make, since it would
#: assign a class on a margin of 0.068.
MANUAL_DECISION_MARGIN = 0.15

# The exclusion step reuses CONTRADICTION_PAIRS / CONTRADICTION_MAX_RATE below. No new
# exclusion threshold is introduced, and none of the v2 acceptance bars change.

# --------------------------------------------------------------------------
# Stage 07 — dominant-clone membership from immunoglobulin evidence
# (predeclared 2026-08-25, BEFORE any cell was assigned)
# --------------------------------------------------------------------------
#
# WHAT THIS AXIS IS, AND WHAT IT IS NOT
# -------------------------------------
# CNV was attempted as an INDEPENDENT malignancy axis and REJECTED, before any disease
# CNV distribution was inspected, because healthy-donor plasma cells failed the
# negative control (donor false-positive 0.0-50.6% at z>3; one donor at median z +3.03).
# It is frozen `NOT_EVALUABLE` cohort-wide.
#
# Patient-specific V-gene usage passed donor and repeated-sample validation and is
# accepted ONLY as a higher-specificity refinement of the immunoglobulin clonality
# axis — NOT as an independent line of evidence. kappa/lambda class and V-gene usage
# are the same molecule and the same biological event. Anything built on them is
# SINGLE-AXIS and must be described that way.
#
# HARD INVARIANT: V ABSENCE IS NOT NEGATIVE EVIDENCE.
# This is 10x 3' data; V segments sit 5' and are captured at ~1 UMI (J segments are
# essentially absent: IGKJ 0.00%, IGHJ 0.04% of cells). A missing dominant-V call may
# be dropout, low depth, a true non-clone cell, or something else. Only POSITIVE V
# detection may establish clone support; absence may never establish incompatibility.

#: Minimum share of a patient's V-positive plasma cells carrying one V gene.
#:
#: Set inside the empirical gap from the donor negative control: healthy donors with
#: >=20 V-positive cells sit at 0.204 / 0.232 / 0.250 / 0.378, evaluable disease
#: patients begin at 0.562. 0.50 is ~32% above the donor maximum and is NOT placed at
#: the disease minimum, so it is not tuned to admit any particular patient. Higher =
#: more specific, which is the direction this stage needs.
DOMINANT_V_MIN_FRAC = 0.50

#: A patient needs this many V-positive plasma cells before a dominant V is called.
V_EVALUABLE_MIN_CELLS = 50
V_PARTIAL_MIN_CELLS = 20

#: ...and this share of its plasma cells carrying any light-chain V.
V_EVALUABLE_MIN_PCT = 0.50
V_PARTIAL_MIN_PCT = 0.20

#: The dominant V must be this many times more often detected in the patient's own
#: plasma cells than in other patients'. Median observed is 23.4x, minimum 3.1x, so
#: 3.0 excludes nothing observed — it is a floor against a V gene so common it carries
#: no membership information, not a selection criterion.
DOMINANT_V_MIN_ENRICHMENT = 3.0

#: A cell counts as expressing an ALTERNATIVE V (positive incompatibility evidence)
#: only at >=2 UMI on a non-dominant V gene. Typical V detection is 1 UMI, so 1 UMI is
#: within noise and must never establish incompatibility.
ALT_V_MIN_UMI = 2

#: Light-chain UMI needed before a cell's kappa/lambda class is called at all.
LC_CLASS_MIN_UMI = 3

#: Share of a cell's light-chain UMI that one chain must hold for a class call.
LC_CLASS_MIN_FRAC = 0.80

# --------------------------------------------------------------------------
# Stage 06 — cytotoxic-lymphocyte lineage axes (2026-08-25, predeclared)
# --------------------------------------------------------------------------
#
# OPERATIONAL CONCLUSION FROM THE LEIDEN-23 TRBC DIAGNOSTIC
# ---------------------------------------------------------
# Isolated TRBC1/TRBC2 expression is insufficient evidence of T-lineage commitment,
# because it frequently occurs WITHOUT coordinated CD3/TRAC expression in cells that
# carry strong NK-lineage evidence.
#
# That is what the data establish, and it is what justifies the split below. Among
# 5,788 Leiden-23 cells with strong NK evidence and no CD3/TRAC, 74.8% were
# TRBC-positive at a median of 2 UMI with 32.3% at >=3 UMI; among TRBC-positive
# "mixed" cells only 12.9% carried both CD3 and TRAC while 28.2% carried neither.
#
# Candidate mechanisms — germline/unrearranged TRB transcription, ambient spillover,
# residual multiplets — are HYPOTHESES. The diagnostic does not distinguish between
# them: TRBC intensity is above a one-molecule ambient profile, TRBC positivity does
# NOT rise with sample-level T-cell abundance (Spearman r = -0.163, p = 0.3), and the
# doublet-score difference (0.139 vs 0.112) is too small to explain the compartment.
# Do not state any of these mechanisms as established.
#
# WHY THIS IS NOT A THRESHOLD CHANGE: under the diagnostic sensitivity analysis ~4,080
# cells move mixed -> NK, and only ~30.2% of those movers carry a single TRBC UMI. The
# transition is therefore driven by the requirement for COORDINATION with
# lineage-defining machinery, not by thresholding away weak signal. No numeric
# parameter is altered by this revision.

#: T-lineage IDENTITY anchors — the CD3 complex and the TCR-alpha constant region.
T_IDENTITY_ANCHORS: tuple[str, ...] = ("CD3D", "CD3E", "CD3G", "TRAC")

#: T-lineage CONTEXT — measured and reported everywhere, never sufficient alone.
T_CONTEXT: tuple[str, ...] = ("TRBC1", "TRBC2")

#: NK identity axis, unchanged from the frozen Part-B cluster-23 analysis.
NK_IDENTITY: tuple[str, ...] = ("KLRD1", "KLRF1", "NCAM1", "FCGR3A", "KLRC1")

#: gamma-delta axis, unchanged. TRDC alone is NOT promoted to identity: in Leiden 23
#: TRGC detection was near-independent of TRDC (28.2% vs 21.6%) and 92.6% of TRDC+
#: cells carried strong NK evidence.
GD_IDENTITY: tuple[str, ...] = ("TRDC", "TRGC1", "TRGC2")

#: Shared cytotoxic-effector STATE. Never establishes any lineage — the same trap as
#: NKG7/GNLY that produced the original NK over-call.
CYTOTOXIC_STATE: tuple[str, ...] = ("NKG7", "GNLY", "PRF1", "GZMB", "GZMA", "CTSW")

# --------------------------------------------------------------------------
# Stage 06 — lineage exclusivity (added v2, 2026-08-25, BEFORE any v2 result)
# --------------------------------------------------------------------------

#: Ambient-robust evidence for each lineage, used ONLY to detect contradictions.
#:
#: These are deliberately NOT the same lists as MARKER_PANEL. A panel used to
#: *identify* a class can afford to include an ambient-prone gene among several; a
#: program used to *accuse* a cell of belonging to another lineage cannot, because a
#: false accusation is created out of ambient rather than hidden by dropout.
#: So the two dominant ambient species in this marrow are excluded here:
#:   - the globins (HBB/HBA1/HBA2) are absent from `erythroid`
#:   - LYZ is absent from `myeloid`
#: Immunoglobulin is the other big ambient species, which is why there is no
#: plasma/B contradiction program at all (see CONTRADICTION_PAIRS).
LINEAGE_PROGRAMS: dict[str, tuple[str, ...]] = {
    "T":         ("CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2"),
    "B":         ("MS4A1", "CD79A", "CD79B", "CD19"),
    "erythroid": ("GYPA", "AHSP", "ALAS2", "CA1"),
    "myeloid":   ("CD14", "FCN1", "MNDA", "ITGAM"),
}

#: Which lineages are INCOMPATIBLE with each project class, from developmental
#: biology — decided before v2 ran, and not derived from any observed cluster.
#:
#: What is deliberately absent is as load-bearing as what is present:
#:   - PlasmaCell is not contradicted by B, and Bcell is not contradicted by plasma:
#:     plasma cells ARE B-lineage, and the plasmablast continuum is real biology.
#:   - HSPC has NO contradictions. Progenitors legitimately co-express lineage-priming
#:     programs; flagging that would be a biology error, not a QC finding.
#: The two pairs this revision exists for are NK<-T and Erythroid<-T.
CONTRADICTION_PAIRS: dict[str, tuple[str, ...]] = {
    "PlasmaCell": ("T", "erythroid", "myeloid"),
    "Bcell":      ("T", "erythroid", "myeloid"),
    "Tcell":      ("B", "erythroid", "myeloid"),
    "NK":         ("T", "B", "erythroid", "myeloid"),
    "Myeloid":    ("T", "B", "erythroid"),
    "Erythroid":  ("T", "B", "myeloid"),
    "HSPC":       (),
}

#: A cell carries strong positive evidence for a lineage when at least this many of
#: that lineage's LINEAGE_PROGRAMS genes are DETECTED (count > 0).
#:
#: Detection, not absence, and this is the point. The instruction that NK must not
#: require literal absence of T transcripts is respected structurally: dropout can
#: only HIDE evidence, so a detection-based rule under-calls contradictions and can
#: never manufacture one from a zero. Two genes rather than one because a single
#: transcript can be ambient or mismapped; two independent genes of the same complex
#: is far harder to explain that way. Not three, because at 1,162 median genes per
#: cell a 3-gene requirement would be defeated by dropout.
CONTRADICTION_MIN_GENES = 2

#: A class fails lineage exclusivity when more than this share of the cells assigned
#: to it carry strong evidence for an incompatible lineage.
#:
#: Set from technical expectation, not from any observed rate. 25% is deliberately
#: PERMISSIVE relative to expected residual technical contamination (post-scDblFinder
#: doublets plus ambient), so the veto triggers only when contradictory evidence
#: affects a substantial fraction of the assigned class — not mild diffuse noise.
#: A class above it is not lightly contaminated; it is the wrong class.
#:
#: The class-level rate alone cannot distinguish diffuse noise from one bad cluster,
#: so `contradiction_concentration` reports the per-cluster breakdown alongside it.
#: That is a reporting aid, not part of the veto. Value fixed before v2; not revisited.
CONTRADICTION_MAX_RATE = 0.25

#: Filled in by stage 06 from `results/06_annotation/annotation_decision.md`: the
#: per-class map of which method's label was taken. Empty until that stage runs, so
#: that no downstream module branches on annotation logic — they read
#: `obs["cell_type"]` and this map explains its provenance.
ANNOTATION_DECISION: dict[str, str] = {}


def all_required_genes() -> frozenset[str]:
    """Flatten REQUIRED_GENES into a single set of canonical symbols."""
    return frozenset().union(*REQUIRED_GENES.values())
