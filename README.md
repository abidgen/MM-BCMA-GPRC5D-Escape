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
2. Integrate across samples (Harmony — benchmarked against scVI, Scanorama and an
   unintegrated baseline with `scib-metrics`, and retained), cluster, and annotate cell types **three ways —
   manual marker panel, CellTypist, and SingleR — then choose per cell type against
   agreement thresholds fixed in advance**, rather than trusting one labeller.
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
6. Test **co-negativity enrichment**: are the *same* cells losing both antigens more
   often than the two antigens' individual negative rates would predict? A 6% escape
   fraction built from two independent partial failures is a different clinical
   proposition from a 6% built by a coordinated antigen-low phenotype, and only the
   second is a problem a second binder cannot solve. The null is stratified by
   sequencing depth, because shallow cells read zero for both genes and would
   otherwise manufacture exactly this signal.
7. Test **whether the escape population is structured or scattered technical noise**,
   on three escalating levels — non-random location in transcriptional space, a shared
   expression program, and finally support from CNV-defined tumor substructure. Only
   the third licenses the word *subclone*; the first two establish an escape-associated
   *state*. Transcriptional clustering has many non-genetic causes, so the claim
   escalates with the evidence rather than ahead of it.
8. Validate the antigen calls against **matched bulk RNA-seq** from the same samples —
   as a check on antigen *abundance* and on whether single-cell zeros are credible, not
   on the escape fraction itself, which is a joint single-cell quantity bulk cannot
   see — and derive a normal-plasma-cell expression baseline from the healthy marrow
   controls for **marrow on-target/off-tumor expression context**.
9. Extend from two antigens to a **combinatorial coverage matrix** over BCMA, GPRC5D,
   SLAMF7, FCRL5, CD38 and others: which target pair or triple covers the most of
   *this* patient's clone — reported as separate columns (uncovered fraction,
   incremental gain from the second target, co-loss enrichment, normal *marrow*
   expression) rather than collapsed into a utility score, since the weights such a
   score needs would encode a clinical judgement this data cannot supply.
10. As an **exploratory** extension, relate escape risk to the immune microenvironment
    via LIANA+ (a native Python reimplementation of CellChat's algorithm) — with the
    patient (not the cell) as the unit of replication and immune-cell abundance
    controlled as a confounder. Exploratory by design: n ≈ 41 patients against hundreds
    of ligand-receptor pairs does not support a confirmatory claim.
11. Synthesize everything upstream into a **six-axis per-patient evidence matrix** —
    measurement robustness, DN structure, DN phenotype, genomic evidence, immune context
    and multi-antigen coverage — kept as **separate columns, never a composite score**.
    The provisional measurement tiers (robust-high / uncertain / robust-low) from stage 08
    are one of those six axes, not folded together with the others: inspecting the joint
    distribution across all six showed the measurement-robust and structurally-supported
    patient sets to be **disjoint**, so no single categorical risk label or rank ordering
    was produced — that decision, and the finding behind it, is itself part of the result.

**Phase 1 is complete.** A planned Phase 2 independently re-runs the same pipeline on a
second, external cohort (GSE117156), to test whether the core finding replicates beyond
this one dataset and sequencing technology.

## Data

[GSE223060](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223060) (scRNA-seq)
and [GSE223061](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223061) (matched
bulk RNA-seq) — public data from *Single-Cell Discovery and Multiomic Characterization
of Therapeutic Targets in Multiple Myeloma* (Cancer Research, 2023), covering 53 bone
marrow samples from 41 myeloma patients across the MMRF Immune Atlas Pilot study and
two WashU cohorts, plus normal bone marrow controls (62 sample entries total in the
GEO archive).

Only processed, *filtered* Cell Ranger output is publicly available for this series.
Raw reads exist but are under controlled access on dbGaP (`phs000159` for the healthy
donors, `phs000748` for MMRF bulk) rather than absent, and no unfiltered
(pre-cell-calling) matrices were deposited anywhere — which rules out formal
ambient-RNA correction (SoupX/DecontX) regardless. See the main project document for how
that is handled instead. `scripts/01_download_data.sh` pulls and unpacks the data
directly from GEO's FTP.

The 62 samples span three collection cohorts on **different 10x chemistries** (WashU
cohort 1 on 3′ v2 with no dead-cell removal; MMRF and WashU cohort 2 on v3.3/v3.2 with
it), which produces a ~1.9× spread in genes detected per cell across cohorts. Since the
headline metric is a fraction of zeros, cohort is carried as a covariate rather than
ignored — see the limitations below.

The matched bulk RNA-seq (GSE223061, **26 samples** with an exact single-cell match,
computed from the GEO metadata rather than inherited) is used as an **orthogonal check
on the antigen quantification** in stage 09. The two bulk cohorts are not the same
assay — MMRF bulk is CD138+ sorted and pairs with malignant-cell pseudobulk, WashU
cohort 1 bulk is unsorted whole marrow and pairs with whole-sample pseudobulk — so they
are compared separately, since pooling them would make a third of the comparisons
measure tumour burden instead of antigen abundance. Bulk
signal where the single-cell data reads zero is quantified evidence *consistent with*
dropout, and it bounds the false-negative rate rather than leaving it as a caveat. It
is not a direct measurement of the dropout rate — the two assays differ in cellular
composition and sensitivity, so the discordance has more than one possible cause.

### Known limitations, stated up front

- **Two error sources, opposite directions.** Ambient RNA makes true negatives look
  positive (deflating the escape fraction); dropout makes true positives look negative
  (inflating it). Dropout is the larger effect here because GPRC5D is a low-abundance
  transcript and the median cell has ~2,000 detected genes. The headline number is
  reported as a bracketed interval; ranking stability across detection thresholds is
  the robustness check that earns a patient a tier, and the tiers — not an ordinal
  ranking, and not any single value — are the deliverable.
- **Transcript is not surface protein.** CAR-T binds protein; this measures mRNA. BCMA
  is actively shed from the cell surface by γ-secretase, and GPRC5D transcript
  correlates imperfectly with surface density.
- **"Subclone" is a claim this data can only sometimes support.** Resolving CNV
  substructure *within* one patient's tumor is much harder than telling tumor from
  normal, and at ~2,000 genes per cell it is often underpowered. Where it is, the
  escape population is reported as an escape-associated *state* and the CNV level as
  *not evaluable* — never as evidence that no subclone exists.
- **Bulk RNA-seq cannot validate the escape fraction.** A tumor that is half
  BCMA-only and half GPRC5D-only looks, in bulk, like both antigens are well
  expressed, while containing no dual-positive cells at all. Bulk constrains each
  gene's abundance; the joint distribution across cells has no orthogonal check here.
- **Marrow expression is not a safety profile.** The healthy-marrow controls give
  normal plasma-cell antigen levels in marrow. GPRC5D's clinically decisive off-tumor
  site is keratinized tissue, which a bone marrow dataset cannot observe at all.
- **Sequencing depth differs by cohort, and the metric is a fraction of zeros.** The
  three cohorts span two 10x chemistry generations and differ ~1.9× in genes detected
  per cell. That is bounded the same way dropout is — as a covariate in the depth
  regression and a stratum in the permutation null — not assumed away.
- **Patient-ID mapping is resolved** (Supplementary Table S1, 2026-08-24): **41
  patients over 53 in-cohort samples**, reproducing the paper. The naive rule's 43
  differed by two independent facts — `25183` is deposited but appears in no
  supplementary table, and `83942`/`MMY83942` are one patient sampled under both
  WashU protocols. S1 carries **no cytogenetics**, so the t(4;14)/1q21 annotation
  remains unavailable and the TC subgroup call stays a transcriptional proxy.
- **Batch integration is benchmarked, and its limits are measured.** Harmony was
  compared against scVI, Scanorama and an unintegrated baseline across three batch
  definitions (`scib-metrics`, immune compartment, provisional CellTypist labels). The
  incumbent configuration was retained. The instructive part: on this dataset the arms
  that score best on conventional scIB are the ones that **merge the cohort-censored
  plasma populations** — up to 20x the incumbent's plasma mixing — because batch
  metrics cannot distinguish "correctly left apart" from "failed to merge". Scoring is
  therefore on the immune compartment only, with plasma mixing reported as a
  diagnostic that is never optimized. **No integration method can restore cells that
  were never deposited**, so this changes nothing about the censoring below.
- **The deposit is pre-filtered, differently in each cohort.** WashU cohorts 1 and 2
  were cut at 10,000 UMIs before deposit; MMRF and the donors were not (MMRF was cut
  at 10% mitochondrial instead). Plasma cells are the highest-RNA-content cells in
  marrow, so that ceiling censors a band enriched **20–70× for `GPRC5D`** — meaning 36
  of the 54 myeloma samples had part of their antigen-positive population removed
  before the data was public. This inflates the escape fraction for those cohorts and
  is carried as a covariate with a truncate-all-cohorts sensitivity analysis, not
  corrected away.

## Setup

Six conda/mamba environment specifications; five environments built, split by actual
dependency-conflict risk (see the main project document for the full reasoning):

```bash
mamba env create -f envs/env-qc.yml            # data loading, QC, scDblFinder (via rpy2)
mamba env create -f envs/env-core.yml          # integration, annotation, malignant calling, scoring, robustness
mamba env create -f envs/env-annotation.yml    # CellTypist + SingleR (isolates R, like env-qc)
mamba env create -f envs/env-communication.yml # LIANA+ (CellChat-equivalent)
mamba env create -f envs/env-integration.yml   # stage 05b only: every integration method + scib-metrics
# envs/env-composition.yml (scCODA) only if the compositional analysis is run —
# it pulls TensorFlow and is kept out of mm-core deliberately

# register a Jupyter kernel per env
mamba run -n mm-qc python -m ipykernel install --user --name mm-qc
mamba run -n mm-core python -m ipykernel install --user --name mm-core
mamba run -n mm-annotation python -m ipykernel install --user --name mm-annotation
mamba run -n mm-communication python -m ipykernel install --user --name mm-communication
mamba run -n mm-integration python -m ipykernel install --user --name mm-integration

# make the package importable from notebooks (src/ layout). --no-deps is not
# optional: envs/*.yml is the dependency manifest, and letting pip resolve into
# these envs has broken them before.
mamba run -n mm-qc   pip install -e . --no-deps
mamba run -n mm-core pip install -e . --no-deps
```

### Tests

```bash
mamba run -n mm-core pytest              # 155 pass, ~27 s, with the deposit present
mamba run -n mm-core pytest -m "not slow"   # skip the two full-cohort passes
```

Two-tier on purpose: tests covering the Ensembl-ID gene-space join, the truncated-
deposit repair and the required-gene assertions need no data at all and run on a fresh
clone (100 pass, 57 skip). Tests that need the extracted deposit are gated and **skip
rather than fail**.

## Pipeline

**The whole analysis runs in notebooks, stages 01 through 12** — every stage is
openable and steppable. Numbering is continuous and 1:1 with output directories:
`notebooks/NN_*.ipynb` writes to `results/NN_*/`.

Notebooks 01-03 wrap the verified acquisition scripts rather than reimplementing them,
so the scripts stay available as a headless CLI fallback (fresh clone, remote box, CI)
and both paths produce byte-identical output:

```bash
bash scripts/01_download_data.sh        # download + unpack from GEO
bash scripts/02_check_files.sh          # confirm per-sample file structure
python scripts/03_build_manifest.py raw/samples   # build sample -> file-path manifest
```

| # | Notebook | Env | Output |
|---|---|---|---|
| 01 | `notebooks/01_download_data.ipynb` | `mm-qc` | `raw/` |
| 02 | `notebooks/02_check_files.ipynb` | `mm-qc` | `raw/` |
| 03 | `notebooks/03_build_manifest.ipynb` | `mm-qc` | `raw/sample_manifest.csv` |
| 04 | `notebooks/04_qc.ipynb` | `mm-qc` | `results/04_qc/` |
| 05 | `notebooks/05_integration_clustering.ipynb` | `mm-core` | `results/05_integration/` |
| 05b | `notebooks/05b_integration_benchmark.ipynb` | `mm-integration` | `results/05b_benchmark/` |
| 06 | `notebooks/06_annotation.ipynb` | `mm-annotation` | `results/06_annotation/` |
| 07 | `notebooks/07_malignant_plasma.ipynb` | `mm-core` | `results/07_malignant_plasma/` |
| 08 | `notebooks/08_dual_antigen_escape.ipynb` | `mm-core` | `results/08_dual_antigen_escape/` |
| 08c | `notebooks/08c_multi_antigen_coverage.ipynb` | `mm-core` | `results/08_dual_antigen_escape/multi_antigen_coverage/` |
| 09 | `notebooks/09_bulk_validation.ipynb` | `mm-core` | `results/09_bulk_validation/` |
| 09b | `notebooks/09b_risk_tiers.ipynb` | `mm-core` | `results/08_dual_antigen_escape/risk_tier_provisional/` |
| 10 | `notebooks/10_dn_coherence.ipynb` | `mm-core` | `results/10_dn_coherence/` |
| 11 | `notebooks/11_immune_context.ipynb` | `mm-communication` | `results/11_immune_context/` |
| 11b | `notebooks/11b_liana_verification.ipynb` | `mm-communication` | `results/11_immune_context/liana_verification/` |
| 12 | `notebooks/12_final_synthesis.py` | `mm-core` | `results/12_final_synthesis/` |

Number order is execution order throughout — 04 → 05 → 05b → 06 → 07 → 08 → 08c →
09 → 09b → 10 → 11 → 11b → 12, no exceptions. Lettered stages (`05b`, `08c`, `09b`,
`11b`) are side-arms feeding the numbered stage they attach to.

The **producers** of the frozen Stage 06–10 tables are committed under `production/`
(recovered verbatim from the session transcripts, 2026-08-26); the notebooks above are
the narrative reports over those tables. See `production/README.md` and
`provenance/README.md`.

A planned Phase 2 independently re-runs the same shape (own `phase2_NN_*` numbered
notebooks/results, never mixed with the numbers above) against a second, external
cohort (GSE117156) once Phase 1 is complete. See the main project document.

See the main project document for the full technical plan and all confirmed data
ground truth, and `mm_analysis_overview.md` for a plain-language walkthrough of the
reasoning behind each stage.

## Repo structure

```
├── mm_analysis_overview.md                # plain-language explanation of the approach
├── mm_dual_antigen_escape_pipeline.md     # pipeline walkthrough (narrative, not code)
├── docs/                                  # dataset ground truth, stage results, decisions archive, envs
├── envs/                                  # six env specs, five built, split by dependency risk
├── src/mm_escape/                         # reusable/testable logic — importable, Codex-reviewable
├── tests/                                 # pytest suite; most of it runs without raw/
├── resources/                             # committed gene-space map + parsed GEO metadata
├── notebooks/                             # numbered 01-12, jupytext-paired (.ipynb gitignored)
├── production/                            # recovered historical producers of the frozen stage 06-10 tables
├── provenance/                            # frozen-artifact manifest + SHA256s + as-built env exports
├── scripts/                               # 01-03 acquisition — CLI fallback, wrapped by notebooks 01-03
├── raw/                                   # data (gitignored — regenerate via scripts/01)
└── results/                               # numbered 04-12, matching notebooks 1:1 (gitignored)
```

## Status

**Phase 1 (stages 01-12, the full GSE223060 pipeline) is complete.** A prior R/Seurat
build had reached data acquisition (62/62 samples) and cohort-wide QC/doublet-removal
before integration; it was set aside and rebuilt in Python from scratch — that build is
preserved in git history under the `r-build-snapshot` tag and was not ported, though the
dataset knowledge it earned carries forward via the main project document.

The full cohort — **62 samples, 204,040 pre-QC cells** — passes QC to **172,940 cells**
and harmonizes to **32,991 genes** with every required marker present, across 30 Leiden
clusters, seven annotated cell classes, and a resolved patient mapping of **41 patients
over 53 in-cohort samples** matching the source paper exactly. From there: malignant
plasma cells identified by light-chain/V-gene clonal evidence (32 patients, 21,906 cells
in the primary denominator); BCMA/GPRC5D antigen scoring and the dual-antigen escape
fraction (median 0.335, range 0.017-0.783), bounded by a threshold-sensitivity band,
depth regressions, bootstrap intervals and a truncate-at-10k censoring check rather than
reported as a bare point estimate; matched-bulk and normal-marrow validation; a
three-level test for whether double-negative cells are structured, phenotypically
distinct, or genomically distinct (4 of 32 patients show non-random structure; a
cohort-level but weakly patient-discriminative transcriptional phenotype; genomic
subclone evidence not evaluable for any patient); an exploratory, ultimately negative
look at the immune microenvironment; a seven-target multi-antigen coverage comparison;
and a final synthesis that assembled all of it into a six-axis per-patient evidence
matrix rather than a risk ranking, because the axes turned out to disagree with each
other in a way that made a single score dishonest. **507+ tests pass** across the suite
built alongside the pipeline (`pytest -m "not slow"`, `mm-core`).

Two data-integrity findings from early in the project are baked into every stage after
them. The cross-reference gene join is on **reconstructed Ensembl IDs**, not symbols,
which recovers 32,991 genes against 22,164 and prevents silent mis-pairing (`TBCE` is a
*different* annotation entry in each build). And `56203_1`, long excluded as an
incompatible reference missing BCMA, is actually a normal sample whose gene file was
truncated mid-write; it is repaired on read behind an assertion and retained.

**Phase 2** (independent validation on GSE117156) is the only work that remains — see
the main project document's Phase 2 section for scope and the explicit no-merge
constraint. See that document for the full technical context, all confirmed data-format
gotchas, the settled architecture decisions, and the complete Stage-12 outcome.
