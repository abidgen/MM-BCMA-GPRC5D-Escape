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
   controls) from a public target-discovery dataset, using MAD-based (median absolute
   deviation) outlier filtering rather than fixed thresholds.
2. Integrate across samples (Harmony), cluster, and annotate cell types.
3. Identify malignant plasma cells via immunoglobulin light-chain restriction
   (kappa/lambda clonality) rather than clustering alone.
4. Score each malignant cell for BCMA and GPRC5D expression against an empirically
   derived detection threshold (ambient-RNA correction isn't possible for this dataset —
   see below — so the positivity threshold is set above the background noise floor
   instead), classifying cells as dual-positive, single-positive, or double-negative.
5. Compute a **dual-antigen escape fraction** per patient — the share of the malignant
   clone that would be invisible to a combined BCMA + GPRC5D strategy — reported with
   bootstrap confidence intervals and a sensitivity band across detection thresholds,
   never as a bare point estimate.
6. Test **whether the escape population is a real subclone or scattered technical
   noise** — cells clustered together in transcriptional space imply a pre-existing
   resistant clone that therapy would select for; cells scattered at random imply
   dropout. Only the first predicts relapse, so the distinction is the point.
7. Validate the antigen calls against **matched bulk RNA-seq** from the same samples,
   and derive a normal-plasma-cell expression baseline from the healthy marrow
   controls — adding an on-target/off-tumor **safety** axis alongside efficacy.
8. Extend from two antigens to a **combinatorial coverage matrix** over BCMA, GPRC5D,
   SLAMF7, FCRL5, CD38 and others: which target pair or triple covers the most of
   *this* patient's clone, traded off against normal-tissue expression.
9. Relate escape risk to the immune microenvironment via LIANA+ (a native Python
   reimplementation of CellChat's algorithm), testing whether high-escape-risk patients
   also show weaker immune cell-cell signaling — with the patient (not the cell) as the
   unit of replication, and immune-cell abundance controlled as a confounder.
10. Output a per-patient risk ranking usable for single- vs. dual- vs.
    sequential-target CAR-T strategy discussions.

A planned Phase 2 independently re-runs the same pipeline on a second, external
cohort (GSE117156) once Phase 1 is complete, to test whether the core finding
replicates beyond this one dataset and sequencing technology.

## Data

[GSE223060](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223060) (scRNA-seq)
and [GSE223061](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223061) (matched
bulk RNA-seq) — public data from *Single-Cell Discovery and Multiomic Characterization
of Therapeutic Targets in Multiple Myeloma* (Cancer Research, 2023), covering 53 bone
marrow samples from 41 myeloma patients across the MMRF Immune Atlas Pilot study and
two WashU cohorts, plus normal bone marrow controls (62 sample entries total in the
GEO archive).

Only processed, *filtered* Cell Ranger output is publicly available for this series —
no raw FASTQ/SRA accession, and no unfiltered (pre-cell-calling) matrices either, which
rules out formal ambient-RNA correction (SoupX/DecontX) — see `CLAUDE.md` for how this
is handled instead. `scripts/01_download_data.sh` pulls and unpacks the data directly
from GEO's FTP.

The matched bulk RNA-seq (GSE223061, ~28 samples overlapping the single-cell cohort)
is used as an **orthogonal check on the antigen quantification** in stage 09: if bulk
shows signal where the single-cell data reads zero, that quantifies the dropout rate
directly rather than leaving it as a caveat.

### Known limitations, stated up front

- **Two error sources, opposite directions.** Ambient RNA makes true negatives look
  positive (deflating the escape fraction); dropout makes true positives look negative
  (inflating it). Dropout is the larger effect here because GPRC5D is a low-abundance
  transcript and the median cell has ~2,000 detected genes. The headline number is
  reported as a bracketed interval, and the defensible claim is the **stability of the
  patient ranking** across detection thresholds, not any single value.
- **Transcript is not surface protein.** CAR-T binds protein; this measures mRNA. BCMA
  is actively shed from the cell surface by γ-secretase, and GPRC5D transcript
  correlates imperfectly with surface density.
- **Patient-ID mapping is provisional** pending the paper's Supplementary Table S1 (a
  naive rule yields 47 patients where the paper reports 41). Every affected number is
  labelled provisional in the output until resolved.

## Setup

Three conda/mamba environments, split by actual dependency-conflict risk
(see `CLAUDE.md` for the full reasoning):

```bash
mamba env create -f envs/env-qc.yml            # data loading, QC, scDblFinder (via rpy2)
mamba env create -f envs/env-core.yml          # integration, annotation, malignant calling, scoring, robustness
mamba env create -f envs/env-communication.yml # LIANA+ (CellChat-equivalent)
# envs/env-composition.yml (scCODA) only if the compositional analysis is run —
# it pulls TensorFlow and is kept out of mm-core deliberately

# register a Jupyter kernel per env
mamba run -n mm-qc python -m ipykernel install --user --name mm-qc
mamba run -n mm-core python -m ipykernel install --user --name mm-core
mamba run -n mm-communication python -m ipykernel install --user --name mm-communication
```

## Pipeline

Numbering is continuous and 1:1 with output directories — `notebooks/NN_*.ipynb`
always writes to `results/NN_*/`, starting right after `scripts/01-03`:

```bash
bash scripts/01_download_data.sh        # download + unpack from GEO
bash scripts/02_check_files.sh          # confirm per-sample file structure
python scripts/03_build_manifest.py raw/samples   # build sample -> file-path manifest
```

| # | Notebook | Env | Output |
|---|---|---|---|
| 04 | `notebooks/04_qc.ipynb` | `mm-qc` | `results/04_qc/` |
| 05 | `notebooks/05_integration_clustering.ipynb` | `mm-core` | `results/05_integration/` |
| 06 | `notebooks/06_annotation.ipynb` | `mm-core` | `results/06_annotation/` |
| 07 | `notebooks/07_malignant_calling.ipynb` | `mm-core` | `results/07_malignant/` |
| 08 | `notebooks/08_antigen_escape_fraction.ipynb` | `mm-core` | `results/08_escape_fraction/` |
| 09 | `notebooks/09_escape_robustness.ipynb` | `mm-core` | `results/09_robustness/` |
| 10 | `notebooks/10_escape_subclone_phenotype.ipynb` | `mm-core` | `results/10_subclone/` |
| 11 | `notebooks/11_cellchat_liana.ipynb` | `mm-communication` | `results/11_communication/` |
| 12 | `notebooks/12_decision_packet.ipynb` | `mm-core` | `results/12_decision_packet/` |

Number order is execution order throughout — 04 → 05 → 06 → 07 → 08 → 09 → 10 → 11
→ 12, no exceptions.

A planned Phase 2 independently re-runs the same shape (own `phase2_NN_*` numbered
notebooks/results, never mixed with the numbers above) against a second, external
cohort (GSE117156) once Phase 1 is complete. See `CLAUDE.md`.

See `CLAUDE.md` for the full technical plan and all confirmed data ground truth,
and `mm_analysis_overview.md` for a plain-language walkthrough of the reasoning
behind each stage.

## Repo structure

```
├── CLAUDE.md                              # full technical plan / decisions log
├── RESUME_HERE.md                         # session state — read first when resuming
├── mm_analysis_overview.md                # plain-language explanation of the approach
├── mm_dual_antigen_escape_pipeline.md     # pipeline walkthrough (narrative, not code)
├── envs/                                  # three conda env specs, split by dependency risk
├── src/mm_escape/                         # all analysis logic — importable, Codex-reviewable
├── notebooks/                             # numbered 04-12, jupytext-paired with src/
├── scripts/                               # 01-03, data acquisition (reused unchanged from the R build)
├── raw/                                   # data (gitignored — regenerate via scripts/01)
└── results/                               # numbered 04-12, matching notebooks 1:1 (gitignored)
```

## Status

A prior R/Seurat build reached: data acquisition solved and verified (62/62 samples),
QC/doublet-removal run on the full cohort, integration not yet run. **Now being
rebuilt in Python from scratch** — that R build is preserved in git history under the
`r-build-snapshot` tag and is not being ported; the dataset knowledge it earned
carries forward via `CLAUDE.md`.

Current state: `raw/` intact (62 samples), `scripts/01-03` in place, no Python
implemented yet. See `RESUME_HERE.md` for exact session state and `CLAUDE.md` for
full technical context, all confirmed data-format gotchas, and the settled
architecture decisions.
