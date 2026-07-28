# Multiple Myeloma Dual-Antigen (BCMA/GPRC5D) Escape Risk Analysis

Single-cell RNA-seq analysis of multiple myeloma bone marrow samples, asking a
question adjacent to CAR-T target selection: **how much of a patient's tumor would
already evade a combined BCMA + GPRC5D targeting strategy, before any treatment even
starts?**

## The problem

BCMA-directed CAR-T (e.g. CARVYKTI) is effective in multiple myeloma, but a real share
of patients relapse via **antigen escape** — the surviving tumor cells are the ones
that stopped expressing BCMA. The field's response has been a second target, **GPRC5D**,
used either as a fallback or paired with BCMA in a dual-target construct.

Most antigen-escape research asks a *before/after* question: did the antigen disappear
after treatment. This project asks a *baseline* question instead: how much dual-antigen
escape risk is already present in a patient's tumor pre-treatment, driven by existing
clonal heterogeneity rather than acquired resistance.

## Approach

1. Load and QC 62 bone marrow scRNA-seq samples (41 multiple myeloma patients + normal
   controls) from a public target-discovery dataset.
2. Integrate across samples (Harmony), cluster, and annotate cell types.
3. Identify malignant plasma cells via immunoglobulin light-chain restriction
   (kappa/lambda clonality) rather than clustering alone.
4. Score each malignant cell for BCMA and GPRC5D expression, classifying cells as
   dual-positive, single-positive, or double-negative.
5. Compute a **dual-antigen escape fraction** per patient — the share of the malignant
   clone that would be invisible to a combined BCMA + GPRC5D strategy.
6. Relate escape risk to the immune microenvironment via CellChat, testing whether
   high-escape-risk patients also show weaker immune cell-cell signaling.
7. Output a per-patient risk ranking usable for single- vs. dual- vs.
   sequential-target CAR-T strategy discussions.

## Data

[GSE223060](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223060) (scRNA-seq)
and [GSE223061](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223061) (matched
bulk RNA-seq) — public data from *Single-Cell Discovery and Multiomic Characterization
of Therapeutic Targets in Multiple Myeloma* (Cancer Research, 2023), covering 53 bone
marrow samples from 41 myeloma patients across the MMRF Immune Atlas Pilot study and
two WashU cohorts, plus normal bone marrow controls (62 sample entries total in the
GEO archive).

Only processed Cell Ranger output is publicly available for this series — no raw
FASTQ/SRA accession exists. `scripts/01_download_data.sh` pulls and unpacks it
directly from GEO's FTP.

## Setup

```bash
mamba env create -f env_creation/environment_working.yml
mamba activate mm-dual-antigen
Rscript env_creation/install_cran_overrides.R
```

## Pipeline

```bash
bash scripts/01_download_data.sh        # download + unpack from GEO
bash scripts/02_check_files.sh          # confirm per-sample file structure
python scripts/03_build_manifest.py raw/samples   # build sample -> file-path manifest
# 04 onward: load/QC/integrate -> annotate -> malignant calling ->
#            antigen scoring -> escape fraction -> CellChat -> final figures
```

See `CLAUDE.md` for the full step-by-step technical plan and `mm_analysis_overview.md`
for a plain-language walkthrough of the reasoning behind each stage.

## Repo structure

```
├── CLAUDE.md                              # full technical plan / project context
├── mm_analysis_overview.md                # plain-language explanation of the approach
├── mm_dual_antigen_escape_pipeline.md     # detailed pipeline documentation
├── env_creation/                          # mamba environment specs
├── scripts/                               # numbered pipeline scripts
├── raw/                                   # data (gitignored — regenerate via script 01)
└── results/                               # pipeline outputs (large files gitignored)
```

## Status

Data download, file-structure verification, and sample manifest construction are
complete (62/62 samples confirmed). QC/integration onward is in progress — see
`CLAUDE.md` for current status and next steps.
