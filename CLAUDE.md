# CLAUDE.md — MM Dual-Antigen (BCMA/GPRC5D) Escape Risk Analysis (Python rebuild)

> **Detail lives in `docs/`.** This file carries the objective, the live conventions,
> the forward instructions for unrun stages, and the current position on every settled
> question. Four companion files carry the evidence and the history:
>
> | file | holds |
> |---|---|
> | `docs/dataset-ground-truth.md` | archive forensics, the Ensembl-ID reconstruction, GEO series metadata, bulk inventory, the S1 patient mapping |
> | `docs/stage-results.md` | run output for stages 01–05b, module validation logs, test-suite baselines |
> | `docs/decisions-archive.md` | every superseded position, why it was wrong, and the two design reviews verbatim |
> | `docs/environments.md` | the five `envs/*.yml` bodies, channel deviations, build traps, verification |
>
> If a `docs/` file and this one disagree, **this one is current**.

## Objective

For each multiple myeloma patient in GSE223060, quantify the fraction of malignant
plasma cells that would evade BOTH a BCMA-directed and a GPRC5D-directed CAR-T/TCE
(the "dual-antigen escape fraction"), then relate that metric to the immune
microenvironment via LIANA+/CellChat. Final output: a per-patient risk **tiering**
(not a false-precision rank ordering — see stage 12) usable for a single- vs. dual-
vs. sequential-target CAR-T strategy discussion.
Built as a portfolio project (Legend Biotech context — BCMA/GPRC5D are their core
CAR-T targets), no fixed deadline — code quality and defensibility of each
analytical choice matter as much as the final numbers.

Most antigen-escape literature asks a before/after question (did the antigen
disappear after treatment). This project asks a baseline question instead: how much
dual-antigen escape risk is already present in a patient's tumor before any
treatment, due to pre-existing clonal heterogeneity.
**The plan went through two design reviews (2026-08-20 scope expansion, 2026-08-21
five-overclaim correction) and both were adopted in full.** Their conclusions are
restated at the stages they bind; the reviews themselves, and every position they
replaced, are in `docs/decisions-archive.md`. Read that before re-opening a settled
question — most "obvious improvements" to this plan have already been tried and
reversed for a stated reason.

### Scientific hierarchy — what this project claims, in order

Added 2026-08-21. The stages are not equally important, and the deliverable should not
present them as if they were. Descending order of centrality:

1. **Primary question** — how heterogeneous is baseline malignant-cell coverage by
   BCMA and GPRC5D across MM patients?
2. **Individual antigen loss** — how often is each target absent on its own?
3. **Primary metric** — threshold-robust dual-negative fraction, with uncertainty.
4. **Clinical value of the second target** — incremental coverage gain,
   `P(A⁻) − P(A⁻ ∩ B⁻)`. Promoted above co-escape (2026-08-21): it is the quantity a
   single- vs. dual-target decision actually turns on, and it stays positive and
   material even when loss is correlated.
5. **Key derived metric** — is dual-negativity enriched beyond independent antigen
   loss (depth-conditioned)? I.e. how much of the pair's expected complementarity is
   eroded by correlated loss.
6. **Reliability** — depth/downsampling, expression-matched controls, bulk antigen
   abundance, normal-PC baseline.
7. **Biological significance** — are DN cells transcriptionally coherent?
8. **Stronger clonal evidence** — does the DN population align with CNV-defined
   malignant substructure, where that is evaluable?
9. **Mechanism** — which programs distinguish DN cells (γ-secretase, MYC, OXPHOS,
   stress)?
10. **Generalization** — does the pattern reproduce independently in GSE117156?
11. **Exploratory extension** — does escape heterogeneity track immune composition and
    signaling (stage 11)?

Read as a sentence: *how much dual targeting adds, where it still fails, whether those
failures occur together beyond technical expectation, and whether the residual
population has coherent biology.*

Stage 11 is last, not co-equal with the antigen analysis. Presenting it as co-equal is
what turns this into a kitchen-sink scRNA project; it is framed as exploratory unless
the signal is strong and stable — see stage 11.

**This is a from-scratch Python rebuild of a project that was substantially built in
R first.** The R version reached: data acquisition fully solved, QC/doublet-removal
run on the full 61-sample cohort, integration/clustering not yet run. None of that
R code is being ported or reused — but every piece of ground-truth knowledge about
the data (format quirks, reference mismatches, patient-mapping gaps) below was
earned the hard way against the real files and carries over directly. Do not
re-discover any of it; read this document first.

---

## Working tree

The R build was removed from the working tree on 2026-08-20 and is preserved in git
history under the **`r-build-snapshot`** tag (`git show r-build-snapshot` to inspect,
`git checkout r-build-snapshot -- <path>` to recover a file). It is not being ported
— only the dataset knowledge in this document carries forward.

Current tree; everything else in this document describes what gets added from here:
```
.
|-- CLAUDE.md
|-- RESUME_HERE.md
|-- README.md
|-- mm_analysis_overview.md
|-- mm_dual_antigen_escape_pipeline.md
|-- raw/                      # the downloaded/extracted data — do NOT redo this
`-- scripts/
    |-- 01_download_data.sh
    |-- 02_check_files.sh
    `-- 03_build_manifest.py
```

`raw/` is the one irreplaceable thing here (`GSE223060_RAW.tar` is ~970 MB and
`GSE223061_RAW.tar` ~77 MB, both gitignored). Everything else in the tree is either
source or regenerable.

---

## Data — the rules that bind code

**Full evidence: `docs/dataset-ground-truth.md`.** Everything below was verified
against the real files; do not re-guess it, and do not re-derive it from the paper —
where the deposit's own metadata disagrees with the files, the files win.

**Source.** GSE223060 (scRNA-seq) + GSE223061 (matched bulk). 62 GEO sample entries:
**54 myeloma** (41 patients, 53 in-cohort samples under the S1 mapping) + **8 donor
marrows** (`BM*`, `ND_*`). Only Cell Ranger *filtered* matrices are public; raw data
exists but sits behind dbGaP controlled access, so **SoupX/DecontX cannot be run**.

**Loader rules** (`io.py`):
- Files are Cell Ranger v2-style — `counts.mtx`, `genes.tsv`, `barcodes.tsv`, all
  uncompressed, one directory level deeper than the GSM folder.
  **`scanpy.read_10x_mtx()` will silently fail**; use `read_mtx` + manual `.obs`/`.var`.
- `genes.tsv` is a **single column of symbols** — no `gene_id` column exists, and there
  are **zero `ENSG` strings in any of the 62 files**. Do not go looking for one.
- Matrices are genes × cells and must be transposed.
- **`56203_1` is repaired, not excluded.** Its `genes.tsv` write was truncated at row
  22185; the canonical column is substituted from the committed gene map behind a
  prefix assertion that raises rather than guessing (`config.TRUNCATED_GENE_FILES`).
  `config.EXCLUDED_SAMPLES` is empty.

**Gene-space rules** (`gene_space.py`):
- Two Cell Ranger references are in play — **33538 genes (37 samples, Ensembl 93)** and
  **33694 (25 samples, Ensembl 84)** — with different HGNC symbol vintages.
- **Join on reconstructed Ensembl ID, never on symbols.** The IDs were reconstructed
  from the public GTFs and verified position-for-position (0 mismatches on both
  builds); the map is committed at `resources/gene_space/` and does **not** need
  regenerating. Symbols retain **22,164** genes; IDs retain **32,991** (+10,827),
  because 11,140 intersected IDs carry a different symbol in each build.
- **Intersect across retained samples, never union** (`anndata.concat(join="inner")`).
  A union would make ~11k genes structurally zero in whole cohorts — indistinguishable
  from the true biological zero this project measures.
- `var_names` = Ensembl ID **through the merge only**, then canonical Ensembl-93
  symbol; the 9 symbols still colliding get `SYMBOL__ENSG…`, never
  `var_names_make_unique()`. Retain `var["ensembl_id"]`, `var["symbol_33538"]`,
  `var["symbol_33694"]`.
- **The required-gene assertions stay, and stay loud.** They caught the `NSD2`/`WHSC1`
  drift that manual inspection missed twice. A missing required gene means *check for a
  legacy symbol*, never *biologically absent*. The four-gene alias map
  (`WHSC1`→`NSD2`, `FAM46C`→`TENT5C`, `WHSC1L1`→`NSD3`, `ATP5A1`→`ATP5F1A`) is a
  **regression assertion only**, not the harmonization mechanism.
- **AnnData, not MuData** — one cell-level modality. Stage-09 bulk is sample-level and
  joins on `sample_id` as a DataFrame.

**Cohort structure, and the confounder it carries:**

| cohort | n | chemistry | dead-cell removal | deposit UMI ceiling |
|---|---:|---|---|---|
| WashU 1 | 23 | 10x 3′ v2 | no | **< 10,000** |
| WashU 2 | 13 | 10x 3′ v3.2 | yes | **< 10,000** |
| MMRF | 18 | 10x 3′ v3.3 | yes | none |
| Donors | 8 | 10x 3′ v3.2 | yes | none |

- **The deposit is pre-filtered, differently per cohort, and the 10,000-UMI ceiling is
  a first-order confounder for stage 08.** Cells above it are enriched **3-21× for
  `TNFRSF17`** and **20-70× for `GPRC5D`**, so 36 of 54 myeloma samples had the
  antigen-positive tail of their own tumours removed before deposit. This inflates
  `frac_double_negative` for WU1/WU2 — **a bias in the project's own direction of
  interest**, baked in and not undoable. Carried as a covariate; stage 08 owes a
  truncate-all-at-10k sensitivity run.
- **Carry `cohort` AND `chemistry` as covariates.** The axis that separates depth is
  cohort (MMRF ≈ 1.9× the others), of which chemistry is one component alongside site
  and protocol. **`n_genes_ref` is not a proxy** — the build split cuts across cohorts.
  Do not quote a "2-3× v2-vs-v3 chemistry effect"; this cohort shows 1.38× with
  overlapping distributions.
- **The 8 donors span both references and have no clone**, so stage 07's negative
  control doubles as a free build/chemistry control.

**Patient mapping — settled by Supplementary Table S1.** 41 patients / 53 in-cohort
samples, reproducing the paper exactly. `25183` is deposited but absent from every
supplementary table — retained, flagged `in_paper_cohort == False`. `83942` (WU1) and
`MMY83942` (WU2) are **one patient**. `io.s1_patient_id` is the mapping;
`_assert_s1_reproduces_the_paper` checks all three counts so a revised S1 fails loudly
instead of quietly moving the metric's denominator.

**The `_N` suffixes are serial disease-course timepoints**, settled by S1 sheet 2
(`27522_1` Primary → `_6` Relapse-3). Not fractions, sorts or replicates.

**Matched bulk (GSE223061): 29 usable, 26 with an exact scRNA match.** Two 114-byte
gzip stubs are excluded (`GSM6939104_MMRF_1505`, `GSM6939120_MMRF_2259`); three IDs
have no scRNA counterpart. **The two bulk cohorts are not the same assay and stage 09
must not pool them** — MMRF bulk is CD138+ sorted (pairs with malignant-cell
pseudobulk), WashU 1 bulk is unsorted BMMC (pairs with **whole-sample** pseudobulk).
Pooling would make 10 of 26 comparisons measure tumour burden instead of antigen
abundance.

**Supplementary tables S3 and S5 are committed and not yet consumed** — S3 is the
paper's 38×38 target-gene correlation matrix (`GPRC5D`×`TNFRSF17` pooled r = 0.064 but
**MMRF +0.62 / WU2 +0.54 / WU1 −0.09**, sign tracking cohort depth exactly — read
stage 08's co-escape result against this); S5 is the paper's own bulk-vs-scRNA r per
gene, a direct comparator for stage 09.

**Datasets evaluated and rejected as sources:** He et al. 2022 (no data availability
statement at all — and the `GSE124310` accession once floated for it belongs to an
unrelated glioblastoma study; do not resurrect it), GSE118900 (no healthy controls,
597 cells across 15 patients).

## Phase 2 (unchanged from the R plan) — external validation on GSE117156

**Sequencing rule: do not start until Phase 1 (the full GSE223060 pipeline,
stages 01-12 below) is fully complete.** No interview deadline — a complete
single-cohort analysis beats a rushed two-cohort one.

Ledergor et al. 2018, *Nat Med* — GSE117156, 51,840 cells, 11 healthy controls + 29
MM patients spanning asymptomatic/pre-diagnosis through NDMM through post-treatment
MRD. Real accession, confirmed data-availability statement.

**Acquisition:** check `scanpy.datasets` / `pooch`-based dataset registries first
before falling back to a manual GEO pull.

**CRITICAL constraint — MARS-seq (plate-based), not 10x droplet-based, unlike
GSE223060. Never merge/integrate the two into one object** — this is a platform
difference, not a correctable batch effect. Run as a fully independent pipeline
(own numbered stage sequence, e.g. `notebooks/phase2_01_qc.ipynb` etc. — a
separate numbering track from Phase 1, prefixed `phase2_` so the two are never
ambiguous on disk): own QC (own MAD thresholds derived from its own distribution),
own clustering/annotation, own malignant calling (kappa/lambda logic carries over
unchanged), own antigen scoring and escape-fraction computation. Comparison across
cohorts is distributional/qualitative only — never a merged per-patient ranking,
no LIANA+/CellChat cross-comparison between cohorts.

---

## QC methodology

Grounded in `sc-best-practices.org`'s QC chapter, adapted for what this dataset allows.
Thresholds are **derived, run and documented** — `results/04_qc/qc_thresholds.csv`, and
`docs/stage-results.md` for what the data forced.

**1. MAD-based outlier filtering, not fixed cutoffs**, re-derived **per cohort** (MMRF
cells carry ~1.9× the genes/cell of WU1's, so a pooled MAD would flag much of WashU as
low-quality for a batch reason). 5 MADs on `log1p_total_counts` and
`log1p_n_genes_by_counts`, plus a dataset-specific `pct_counts_mt` cap that is
**one-sided**: an unusually *low* mitochondrial fraction is not itself evidence of a
low-quality cell, so only the high tail filters.
The tutorial's numeric MAD counts and 8% mito cap are defaults for healthy PBMC/BMMC —
**never copy them onto myeloma marrow.**

**2. `pct_counts_in_top_20_genes` is computed and reported but never filters.** In
myeloma marrow an Ig-dominated library is a plasma cell's *normal* state, not a defect:
the flagged decile is 21× enriched for `TNFRSF17` and its `SDC1` detection is 18.8% vs
0.0% elsewhere. A MAD cut on it would preferentially delete antigen-**positive**
malignant cells and inflate `frac_double_negative`. It is kept because it is one of the
few ambient-Ig handles available at all. Re-enabling it is one explicit argument to
`flag_outliers` (`qc.DEFAULT_FILTERS` vs `qc.ALL_FLAGS`), for a sensitivity run.

**3. QC annotates, it does not delete.** Every stage-04 checkpoint holds every barcode
with `obs["keep"]` set; filtering happens at stage 05. For a fraction-of-zeros metric,
"does this survive a different QC?" is a question that *will* be asked, and it is only
answerable if the filtered cells are still on disk.

**4. Ambient RNA cannot be formally corrected** — SoupX/DecontX need the unfiltered
matrices this deposit lacks. Mitigate at stage 08 instead: derive an empirical noise
floor for BCMA and GPRC5D from cell types with no reason to express either (T, NK,
myeloid), where any nonzero signal is ambient by definition, and set the positivity
threshold above it. Document the value and the population it came from.

**5. Dropout is the opposite-signed bias, and it is the larger one.** It must be
**bounded, not just mentioned**. `frac_double_negative` is a fraction of zeros — the
noisiest quantity scRNA-seq produces — with two error sources pointing opposite ways:

| Source | Effect on a cell's antigen call | Effect on `frac_double_negative` |
|---|---|---|
| Ambient RNA | true negative reads as positive | **deflates** |
| Dropout / limited sensitivity | true positive reads as negative | **inflates** |
| Residual normal PCs called malignant | dilutes with antigen-lower cells | **inflates** |
| Deposit UMI ceiling (WU1/WU2) | antigen-positive tail removed | **inflates** |
| mRNA measured, protein targeted | see stage 12 limitation | direction unknown |

Dropout matters more than usual here because **`GPRC5D` is a low-abundance GPCR
transcript** and this cohort is shallow: **1,162 median genes/cell** over 172,940
post-QC cells (1,521 in the plasma compartment). Stage 05 measured the gap directly —
`GPRC5D` mean **0.061** vs `TNFRSF17` **0.492**, an 8× difference, and `GPRC5D` fails
HVG selection entirely. **BCMA-negative and GPRC5D-negative calls are therefore not
equally reliable**; a materially larger share of GPRC5D zeros are technical.
Mitigations live in stage 08 (sensitivity band, depth-matched checks,
expression-matched false-negative floor) and stage 09 (bulk cross-check). **Neither
bias may be left unquantified, and the headline number is reported as a bracketed
interval, never a bare point estimate.**

Doublet detection: **`scDblFinder` (R) via an isolated `rpy2` bridge**, per the
Xi & Li 2021 benchmark — not swapped for a pure-Python alternative for language purity.

## Environments

**Full specs, build traps and verification: `docs/environments.md`.** Five built, split
by actual dependency-conflict risk (not one per stage); two exist to quarantine R.

| env | stages | why isolated |
|---|---|---|
| `mm-qc` | 04 | `scDblFinder` + `rpy2` — keeps R out of everything else |
| `mm-core` | 05, 07-10, 12, **tests** | the shared scanpy/pydeseq2/decoupler stack |
| `mm-annotation` | 06 | `SingleR` + `celldex` — the second R quarantine |
| `mm-communication` | 11 | `liana`/`omnipath`'s version-sensitive tree |
| `mm-integration` | 05b | every integration method under comparison + `scib-metrics` |

`envs/env-composition.yml` (scCODA) is **written but deliberately not built** — it pulls
TensorFlow and must not destabilize `mm-core`. Build it only if stage 06's
compositional comparison is actually run.

**Two hard rules:**
1. **Never `pip install` into these envs casually.** Installing `decoupler==1.8.0` once
   silently downgraded numpy 2.5.2 → 1.26.4 *and* numba, breaking scanpy, scipy,
   pydeseq2 and zarr at once; the repair was to delete and recreate the env. The two
   legitimate pip entries (`infercnvpy`, `decoupler==2.2.0`) live in the yml so a
   rebuild reproduces them. If a pip install is unavoidable, rebuild afterwards.
2. **`decoupler` 2.x is an API rewrite** (`dc.mt.*` / `dc.op.*`, not `dc.run_mlm`).
   Stage 10 must be written against 2.x; do not follow 1.x tutorials, including
   `sc-best-practices`'s. Downgrading is not an option — 1.x pins `numpy<2`.

Note `mm-communication` deliberately runs a different scanpy/anndata (1.12.3 / 0.12.19)
from the others (1.11.5 / 0.13.2). Stage 11 meets `mm-core` only on disk; if an `.h5ad`
compatibility problem ever appears, that version gap is the first place to look.

## Code structure — notebooks for every stage + a modular `src/` package

**`src/mm_escape/` is an importable package, named by function, NOT numbered** —
it's a library, not a pipeline sequence, and some modules (e.g. `plotting.py`) get
used out of stage order:
```
src/mm_escape/
|-- __init__.py
|-- config.py           # thresholds, gene sets, exclusions, MARKER_PANEL,
|                        # STATE_PROGRAMS, TC_GENES, ANNOTATION_DECISION,
|                        # REQUIRED_GENES — env-var overridable, same convention
|                        # as the R build's lib/00_config.R
|-- io.py                # manifest loading, the read_mtx wrapper handling this
|                        # dataset's single-column genes.tsv / counts.mtx naming
|-- gene_space.py         # cross-sample gene intersection + hard assertions
|-- qc.py                  # MAD-based outlier calling, scDblFinder rpy2 bridge
|-- integration.py          # normalize, HVG, PCA, Harmony, clustering, UMAP
|-- annotation.py            # celltypist / marker-panel scoring
|-- malignant.py               # kappa/lambda ratio restriction, infercnvpy
|                               # (required concordance layer, not optional)
|-- antigen.py                  # antigen scoring incl. ambient-noise-floor
|                                # threshold, escape-fraction computation,
|                                # threshold sensitivity band, coverage matrix,
|                                # co-negativity enrichment (2x2 + depth-cond. null)
|-- robustness.py                 # bootstrap CIs (flat + hierarchical),
|                                  # depth/downsampling checks, expression-matched
|                                  # false-negative floor + detection curve,
|                                  # depth-stratified permutation nulls
|-- bulk.py                         # GSE223061 TPM loading, pseudobulk-vs-bulk
|                                    # correlation (stage 09)
|-- subclone.py                       # DN coherence levels 1-3 (kNN/Moran's I,
|                                      # program coherence, CNV support),
|                                      # per-patient malignant subclustering,
|                                      # pseudobulk DE + decoupler (stage 10)
|-- communication.py                    # LIANA+ wrapper, per-patient LR scores,
|                                        # continuous escape-fraction modelling
`-- plotting.py                           # shared plot helpers, consistent styling
```

Five modules are new relative to the original plan (`robustness.py`, `bulk.py`,
`subclone.py`, plus expanded `antigen.py`/`malignant.py`/`communication.py`). The
naming rule is unchanged: **named by function, never numbered** — `robustness.py`
is used by stages 08, 09 and 10, which is exactly why it cannot carry a number.

**Notebooks and results ARE numbered, matching each other 1:1** — this is the
actual tracking mechanism: notebook `NN_*.ipynb` produces `results/NN_*/`. Numbers
continue straight on from the existing `scripts/01-03`, so the whole pipeline is
one unambiguous sequence from `01` to `12`. **Every stage is a notebook** — there is
no stage you cannot open and step through:

| # | Notebook | Env | Output dir |
|---|---|---|---|
| 01 | `notebooks/01_download_data.ipynb` | `mm-qc` | `raw/` |
| 02 | `notebooks/02_check_files.ipynb` | `mm-qc` | `raw/` |
| 03 | `notebooks/03_build_manifest.ipynb` | `mm-qc` | `raw/sample_manifest.csv` |
| 04 | `notebooks/04_qc.ipynb` | `mm-qc` | `results/04_qc/` |
| 05 | `notebooks/05_integration_clustering.ipynb` | `mm-core` | `results/05_integration/` |
| 05b | `notebooks/05b_integration_benchmark.ipynb` | `mm-integration` | `results/05b_benchmark/` |
| 06 | `notebooks/06_annotation.ipynb` | `mm-annotation` | `results/06_annotation/` |
| 07 | `notebooks/07_malignant_calling.ipynb` | `mm-core` | `results/07_malignant/` |
| 08 | `notebooks/08_antigen_escape_fraction.ipynb` | `mm-core` | `results/08_escape_fraction/` |
| 09 | `notebooks/09_escape_robustness.ipynb` | `mm-core` | `results/09_robustness/` |
| 10 | `notebooks/10_escape_subclone_phenotype.ipynb` | `mm-core` | `results/10_subclone/` |
| 11 | `notebooks/11_cellchat_liana.ipynb` | `mm-communication` | `results/11_communication/` |
| 12 | `notebooks/12_decision_packet.ipynb` | `mm-core` | `results/12_decision_packet/` |

**Number order IS execution order** — 04 → 05 → 06 → 07 → 08 → 09 → 10 → 11 → 12, with
no exceptions. The scope expansion renumbered the sequence rather than appending to it,
precisely so this holds. **Never add a stage whose number and run position disagree.**

Each `results/NN_*/` holds everything that stage produces (figures, CSVs, checkpointed
`.h5ad`) — no flat ad hoc filenames in one `results/` dump like the R build had.
**Phase 2 uses a `phase2_NN_*` prefix**, never plain `NN_*`.

Notebooks are **paired with `jupytext`** (percent format); the `.ipynb` is gitignored
and generated from the committed `.py` (see `.gitignore` for how to flip that).

**Division of labour.** Notebooks are the analysis — narrative, plots, intermediate
inspection, the reasoning a reader steps through. `src/mm_escape/` holds what earns
being a library: **reusable** (more than one stage), **testable** (worth asserting on
independently), or **fiddly** (the `read_mtx` loader, symbol harmonization, the
noise-floor derivation). The test is reuse and testability, **not line count** — a
single-use plotting call belongs in the notebook; a threshold calculation duplicated
across three notebooks does not. This protects two things: Codex reviews `.py` diffs
rather than notebook JSON, and logic with one home cannot drift between copies — which
is exactly what made `mm_dual_antigen_escape_pipeline.md` go stale once already.

**`scripts/01-03` are retained as a headless CLI fallback**, not deleted — already
solved and verified, useful on a fresh clone or in CI. Notebooks 01-03 **call into
them** and must never reimplement them; the scripts' output is the contract and
byte-identical parity is verified both directions.

---

## Working with Claude Code + ChatGPT Codex on this repo

- Claude Code implements against this file and the numbered stages above; keep
  `RESUME_HERE.md` current as work proceeds (session state, not decisions —
  decisions belong here).
- Route code review through Codex on the **paired `.py` files**, not raw
  `.ipynb`. `jupytext --sync` keeps the pair in lockstep; review the `.py`, then
  re-sync into the notebook.
- Every module in `src/mm_escape/` should be reviewable and lightly testable
  independent of any notebook.

- **Most tests need no data at all.** The things most worth protecting — the
  Ensembl-ID join, the `make.unique` reimplementation, the `56203_1` truncation repair,
  the required-gene assertions, the GEO metadata join — are exercised entirely by what
  is committed under `resources/`.
- **Data-backed tests are gated and skip rather than fail** (`conftest.requires_data`,
  `requires_s1`, `requires_r`). `pytest -m "not slow"` skips the two full-cohort passes.

Data-backed tests run over the four canonical samples in `conftest.CANONICAL_SAMPLES`,
chosen because they cover the deposit's failure modes — `MMRF_1695` (33538 build),
`27522_1` (33694, legacy symbols), `BM4` (donor), `56203_1` (truncated). **Add cases
there** rather than picking new samples ad hoc.

Assertions encode the cohort's known invariants (204,040 pre-QC cells, 172,940 post-QC,
32,991 intersected genes, 11,140 drifted symbols, 26 matched bulk samples, the
per-cohort depth ordering), so a regression in the loader, the gene map or the metadata
tables fails loudly instead of quietly changing a number. Current counts and gate
details: `docs/stage-results.md`.

---

## Pipeline stages (numbers match notebook filenames and result directories)

**Stages 01-05b are RUN.** The blocks below say what each produced and what it binds
downstream. **Full run output, tables and figures: `docs/stage-results.md`**, and the
`results/*.csv` files are the source of truth over any table reproduced anywhere.

**01-03 — Data acquisition** (`mm-qc`) — **RUN.** 62/62 `triplet-ok`; manifest
byte-identical via notebook and CLI. Notebooks **wrap** `scripts/01-03` (bash ones via
`subprocess`, `03` by importing `build_manifest()`) and must never reimplement them —
the scripts are the contract and are kept as a headless fallback.
- **Only 3 distinct `genes.tsv` checksums exist across 62 files.** That is what makes
  the Ensembl-ID reconstruction possible; notebook 02 asserts exactly 3, so a new
  reference build fails loudly instead of merging silently.
- **Under a Jupyter kernel `sys.argv[1]` is `-f`**, so the script's module-level
  `RAW_DIR` evaluates to `Path('-f')`. Always pass paths explicitly; never use
  `mf.RAW_DIR`. The manifest holds **repo-root-relative** paths.

**04 — QC + doublets** (`mm-qc`) — **RUN.** 204,040 → **172,940 cells (84.8%)** over
all 62 samples. Per-cohort MAD thresholds, `scDblFinder` on every sample, one
checkpoint per sample under `results/04_qc/samples/` with **every barcode retained**
and `obs["keep"]` set. Removal is stable from 4 to 6 MADs (16.5% → 14.6%), so nothing
downstream hangs on the exact count.
- **Binds downstream:** the deposit's per-cohort pre-filtering (see Data above) is
  quantified here — `results/04_qc/umi_censoring_effect.csv`. Stage 08 **owes** the
  truncate-all-cohorts-at-10k sensitivity analysis. Not corrected at stage 04, because
  matching MMRF to WashU would discard 42% of MMRF's cells to make every cohort equally
  damaged.

**05 — Gene-space intersection, integration, clustering** (`mm-core`) — **RUN.**
172,940 × **32,991 genes**, 30 Leiden clusters → `results/05_integration/integrated.h5ad`.
Gene space landed exactly as predicted: 22,164 on symbols → 32,991 on Ensembl IDs,
11,140 drifted symbols joined correctly, `NSD2` resolving against `WHSC1`. CP10K +
`log1p` with raw integers in `layers["counts"]`; 2,000 HVGs with `batch_key="patient_id"`;
50 PCs. **Peak ~20 GB RAM — the one stage that concatenates the matrices, and therefore
the project's machine-size constraint.**
- **Harmony works on the immune compartment (cohort-mixing entropy 0.751) and does NOT
  work on the plasma compartment (0.105).** The three largest plasma clusters are one
  per cohort, each spanning ~30 patients — which rules out a patient-private clone
  (that would fragment into ~41 clusters, not three cohort-shaped ones).
- **This is the stage-04 censoring, measured rather than inferred.** Median UMIs in the
  plasma compartment: MMRF 22,477 vs WU1 5,036 (**4.5×**, against 1.8× for non-plasma).
  MMRF's two largest plasma clusters have 68% and 88% of cells above 10,000 UMIs —
  cells the WashU deposits *cannot contain*. WashU's plasma clusters instead press up
  against the ceiling. **No one-to-one population correspondence remains for any method
  to recover, and no correction restores cells that were never deposited.**
- **Forbidden:** reading any cross-cohort comparison of malignant-cell state off this
  embedding. **Contained by design:** the embedding feeds only stages 06 and 11;
  per-cell antigen calls read `layers["counts"]` and are integration-independent; stage
  10's malignant subclustering is per-patient and un-integrated. This is the answer to
  "did Harmony distort your escape fractions?" — it cannot.
- **`GPRC5D` is not a highly variable gene** (mean 0.061 vs `TNFRSF17` 0.492; HVG in 6
  of ~50 patients). This affects nothing mechanically — its value is **evidential**, as
  this cohort's own first number behind the dropout-asymmetry argument. The antigen
  panel is deliberately **not** forced into the HVG set.

**05b — Integration-method benchmark** (`mm-integration`) — **RUN. The incumbent
survived; nothing about stage 05's output changed.** Seven arms (unintegrated; Harmony
/ scVI / Scanorama on `sample_name`; the stage-05 incumbent; Harmony / scVI on
`cohort`) scored with `scib-metrics` against provisional CellTypist labels. It reads
`integrated.h5ad` read-only and asserts byte-identity afterwards.
- **It takes a letter, not a number**, because number order is execution order with no
  exceptions and this is a side-comparison feeding the stage-05 *choice*.
- **The scoring is deliberately not standard scIB, and that design point matters more
  than the leaderboard.** scIB's batch metrics cannot distinguish *correctly left apart*
  from *failed to merge*, so on this dataset a naive global ranking **structurally
  rewards overcorrection**. Therefore the immune compartment is **scored** and the
  plasma compartment is only **diagnosed** — plasma mixing never contributes positively.
- **The failure mode this design was built to catch actually occurred.** A standard
  global scIB benchmark would have selected `harmony_sample` (best batch **and** bio
  scores) — which mixes plasma **13.5×** harder than the incumbent while encoding more
  depth; `harmony_cohort` reaches **20.2×**.
- **The incumbent is the worst batch corrector of the seven (0.427, below unintegrated's
  0.450)** and by a wide margin the least depth-encoding (R² **0.369** vs 0.51-0.69) and
  least plasma-merging (**0.038**). Correcting on `patient_id` + `n_genes_ref` + `cohort`
  trades cohort mixing for exactly the properties this dataset needs. **This finding
  exists only because the incumbent was entered as its own arm rather than assumed.**
- **Stated honestly: `depth_ok` did all the gating**, at a +0.05 tolerance fixed before
  the spread was known. Two things keep that from being a threshold artifact — three
  arms independently fail `overcorrection_ok` and lose even with the depth criterion
  removed, and relaxing depth enough to admit anything admits only the two weakest batch
  performers. **The tolerance is not re-tuned after the fact**; doing so is exactly the
  post-hoc rationalisation pre-declaration exists to prevent.
- **Propagate to stage 08: no integration method restores cells that were never
  deposited.** A well-mixed latent space has not undone the ascertainment bias in the
  raw counts stage 08 reads. The truncate-all-at-10k sensitivity analysis is still owed,
  and selecting a fancier method must never create the impression the censoring was
  "handled".

---

**06 — Annotation** (`notebooks/06_annotation.ipynb`, `src/mm_escape/annotation.py`;
env: **`mm-annotation`**). Follows `sc-best-practices.org`'s annotation chapter.
**Three methods are run and compared, and the choice is made per class against
thresholds declared in advance** — not "`celltypist` and/or marker scoring", which
would leave a load-bearing decision to be settled implicitly by whatever ran first.

**What is load-bearing.** This stage feeds three things and only three:

| Downstream need | Labels that matter | Cost of getting it wrong |
|---|---|---|
| Stage 07 malignant calling | PlasmaCell vs. everything else | wrong plasma-cell set → wrong denominator for `frac_double_negative` |
| Stage 08 ambient noise floor | T / NK / myeloid purity | a plasma cell leaking into the "confidently antigen-negative" population inflates the floor and biases every antigen call |
| Stage 11 confounder control | T/NK abundance, ideally with subsets | composition artifact misread as immune evasion |

Fine-grained subtypes (CD4 Tcm vs. Tem) benefit stage 11 only and **must never be the
reason a method is chosen**.

**Method A — manual, marker-based**, at cluster level (not per cell; clustering absorbs
dropout, which matters at 1,162 median genes/cell). `scanpy.tl.score_genes` over the
project panel: PlasmaCell (SDC1/CD38/MZB1/XBP1/IRF4), Bcell (MS4A1/CD79A/CD19), Tcell
(CD3D/CD3E/CD8A/CD4), NK (NCAM1/NKG7/GNLY), Myeloid (CD14/LYZ/ITGAM), Erythroid
(HBB/GYPA), HSPC (CD34/KIT). Then `dotplot(..., standard_scale="var")` as saved
assignment evidence, and `rank_genes_groups(method="wilcoxon")` +
`filter_rank_genes_groups`. **The DE step is not optional** — it is the only thing that
can surface a population the seven-class panel does not cover (pDC, erythroid
progenitors, a doublet-driven cluster). Record ambiguous clusters as ambiguous.

**Method B — `celltypist`.** Normalize to 10,000 counts/cell then `log1p` (its stated
requirement); run on **expression, never the Harmony embedding**. Enumerate models with
`models_description()` rather than assuming names; default `Immune_All_Low.pkl` +
`Immune_All_High.pkl`, plus a healthy-BM model if one ships.
`annotate(..., majority_voting=True, over_clustering=<stage-05 leiden key>)` — passing
the existing Leiden key is what makes the methods comparable (same partition, different
labelings). Retain `conf_score`.

**Method C — `SingleR`**, chosen to cover CellTypist's blind spot rather than duplicate
its strengths: `Immune_All_*` is immune-only, so its predictable failure is erythroid and
HSPC, which `celldex`'s `NovershternHematopoieticData` covers
(`HumanPrimaryCellAtlasData` as fallback; verify both exist in the installed `celldex`).
Run with `clusters=<leiden key>`; retain `pruned.labels` and the delta/score matrix as
SingleR's own low-confidence signal.

**Two caveats about plasma cells, pointing opposite ways.** The references contain
*normal* plasma cells, so malignant PCs will be labelled "Plasma cells" — correct and
sufficient, since malignant-vs-normal is stage 07's job; do not count a method down for
failing to find the tumor. But conversely, because no malignant class exists in the
reference, a heavily aneuploid clone may be labelled something else or split across
labels. **Run the plasma-cell marker-coverage test on myeloma marrows specifically, not
only the donors** — otherwise a systematic failure on exactly the cells this project
measures passes unnoticed. Stage 05's three cohort-specific plasma clusters cost
annotation nothing (all three annotate as PlasmaCell) but make this check more important.

**The comparison** (artifacts to `results/06_annotation/`):
1. Confusion matrices — manual × CellTypist, manual × SingleR, **CellTypist × SingleR** —
   at cluster and cell level, plus ARI and per-class F1/Jaccard. ARI alone is
   insufficient: a method can score well overall while failing on plasma cells, the one
   class that must not be wrong. The two automated methods agreeing with *each other* is
   the strongest evidence available, being trained on different references.
2. **The marker-coverage test — the decisive one.** Dotplot the *manual* panel grouped by
   each *automated* method's labels. If CellTypist's T cells are CD3D/CD3E-high, its
   plasma cells MZB1/SDC1-high, and so on, the automated labels already encode what the
   manual panel encodes.
3. Confidence/coverage: `conf_score` per cluster, SingleR pruned-`NA` rate and deltas,
   the fraction of cells unassigned or labelled outside the panel.

**These are concordance scores, not validation accuracy.** The manual annotation is a
third opinion from the same expression matrix, not ground truth — "F1 against manual"
measures *agreement*, and the two automated methods agreeing is agreement between
references that share marker-biology priors and may share their blind spots. Label it
all **concordance**. The load-bearing evidence of *biological* validity is the
marker-coverage test: a label set can be perfectly self-consistent and biologically
wrong. Headline figure, in this order — **label concordance → marker-expression validity
→ uncertainty/unassigned** — concordance first because it is the weakest of the three.

**The decision rule — per class, declared before looking**, so "choose the best" cannot
become post-hoc rationalization. A class goes to an automated method when its
marker-coverage test passes, its own confidence signal is not flagging the cluster, and
concordance with manual clears a pre-set F1 bar:
- **PlasmaCell: F1 ≥ 0.95** — strictest, because it sets the metric's denominator.
- **T / NK / Myeloid: F1 ≥ 0.90** — these define stage 08's noise floor.
- **Bcell / Erythroid / HSPC: F1 ≥ 0.85** — nothing downstream is load-bearing on these.

Where both qualify, take the higher concordance and record that both passed. Where
neither does, the class falls back to the manual cluster label. **A failed
marker-coverage test vetoes a class regardless of concordance** — high agreement on a
biologically unsupported label is agreement on an error. Expected but **not** assumed:
immune classes and plasma cells from CellTypist, erythroid/HSPC from SingleR or manual.
**The numbers decide.** Outcome to `results/06_annotation/annotation_decision.md`.

**Interface contract — downstream stages read `cell_type` and nothing else:**
- `obs["cell_type"]` — the canonical load-bearing coarse label, seven project classes.
- `obs["cell_type_fine"]` — CellTypist fine label where available, else `NA`. Stage 11
  only; never load-bearing.
- `obs["annotation_source"]` — per cell: `celltypist` | `singler` | `manual`. Without it
  a mixed-provenance label column is untraceable.
- `obs["annotation_conf"]` — the winning method's confidence for that cell.
- `config.ANNOTATION_DECISION` — the per-class method map, so no downstream module
  branches on annotation logic.

This decoupling is deliberate: the comparison can be redone or reversed later without
touching stages 07-12.

**Orthogonal cell-state programs — continuous, never categorical.** A cell has one
`cell_type` but can carry several active programs at once. Score with
`scanpy.tl.score_genes`, store as float `obs` columns, **never collapse into
`cell_type`**:

| Program | Why it matters here |
|---|---|
| Cell cycle (`MKI67`, `TOP2A`, `PCNA`) | a proliferative escape subclone is a different risk from a quiescent one — feeds stage 10 |
| Interferon (`ISG15`, `IFI6`, `STAT1`, `MX1`) | immune-pressure marker; feeds stage 11 |
| Antigen presentation (`B2M`, HLA I/II) | `B2M` loss is a documented myeloma immune-escape route. CAR-T is MHC-independent so it does not affect the escape metric, but it is a *competing* evasion mechanism and belongs in the stage 11/12 interpretation |
| UPR (`XBP1`, `ATF4`, `HSPA5`, `DDIT3`) | plasma cells are professional secretors; UPR tone is core PC biology |
| Hypoxia / stress | standard confounder — cheap to score, expensive to discover late |

A cycling plasma cell is `cell_type == "PlasmaCell"` **plus** a high cell-cycle score,
not a "Cycling" identity. If any method emits a proliferation label as an identity,
remap it to PlasmaCell + score.

**Per-patient composition is a first-class output**, not a by-product: malignant-PC
fraction of the marrow (tumor burden, context for stage 12) and T/NK/myeloid abundance
(stage 11's primary confounder). Any composition *comparison* between groups uses
**`scCODA`**, not a raw proportion test — proportions are compositional (they sum to 1,
so one type rising forces others down) and naive per-type tests are anticonservative.
`scCODA` needs its own env.

**07 — Malignant plasma cell identification**
(`notebooks/07_malignant_calling.ipynb`, `src/mm_escape/malignant.py`; env:
`mm-core`). Subset to plasma cell clusters. Kappa (`IGKC`) vs. lambda
(`IGLC1-7`) restriction scoring per cell; per-patient dominant restriction class
(>90% in an involved marrow) marks malignant cells. Prefer scVDJ-seq clonotype
calls over the restriction proxy if available (check GEO supplementary files —
unconfirmed).

Three requirements, all motivated by the same fact: **this stage defines the
denominator of the headline metric, so its errors propagate straight into
`frac_double_negative`.**
- **Ratio-based restriction, not presence/absence.** IG transcripts are the most
  ambient-contaminated genes in this tissue — plasma cells secrete enormous quantities
  of Ig mRNA into the droplet background — so a presence-based kappa-or-lambda call is
  much noisier than it looks. Use the per-cell kappa:lambda **ratio**, which is robust
  to a shared additive background in a way a presence call is not.
- **`infercnvpy` is REQUIRED, not optional.** Residual *normal* plasma cells are
  antigen-lower than malignant ones, so every normal PC mistakenly called malignant
  inflates the escape fraction, and an independent CNV call is the only cross-check
  available. Use minority-restriction cells as the per-patient normal reference.
  **Report the light-chain vs. CNV agreement rate as a stage output** — a low rate
  invalidates stage 08 and must halt the pipeline, not be noted and passed over.
- **Normal-BM negative control.** Run the identical logic on `BM2`, `BM4`, `BM5`, `BM6`
  and the `ND_*` samples. Normal marrow is polyclonal, so the correct output is *no
  malignant calls*; if the method calls a clone there, it is broken and nothing
  downstream is trustworthy. The cheapest strong validation available for the most
  method-dependent step in the project, on samples already downloaded and otherwise idle.

**Malignant-evidence tiers, not a binary call.** With two
independent lines of evidence there is no reason to collapse them into one bit and
throw the disagreements away. Emit `obs["malignant_confidence"]`:

| Tier | Evidence | Used where |
|---|---|---|
| `high` | plasma-cell identity + light-chain restriction + CNV abnormality | primary analysis **and** the restricted sensitivity re-run |
| `probable` | plasma-cell identity + strong restriction, CNV inconclusive/not evaluable | primary analysis only |
| `uncertain` | the two lines disagree | excluded from the metric; counted and reported |

**Stage 08's headline ranking is then re-run on `high` cells only, as a sensitivity
analysis.** If the patient ordering survives that restriction, the metric is not an
artifact of the weaker malignant calls — which is a much stronger statement than the
agreement rate alone. Same distinction as stage 10's coherence hierarchy: `uncertain`
is a *reported quantity*, never a silent drop, and `probable` from CNV being **not
evaluable** is not the same as CNV being negative.
**08 — Antigen scoring + dual-antigen escape fraction**
(`notebooks/08_antigen_escape_fraction.ipynb`, `src/mm_escape/antigen.py` +
`robustness.py`; env: `mm-core`). Per malignant cell: positivity for BCMA
(`TNFRSF17`), GPRC5D, and backups (`SLAMF7`, `FCRL5`), using the ambient-noise-floor
threshold (**not** a naive `>0` call). Classify into `dual_positive` / `BCMA_only` /
`GPRC5D_only` / `double_negative`. Compute `frac_double_negative` per patient — **the
core novel metric**. The per-patient aggregation uses the S1 mapping (41 patients / 53
in-cohort samples), which is settled.

**This stage carries the metric's defense. A single point estimate off a single
threshold is not a defensible deliverable.**
- **Threshold sensitivity band.** Compute `frac_double_negative` under at least three
  calling rules — naive `>0`, the ambient noise floor, a stricter floor. The claim is
  **not** any one number; it is the **stability of the patient ordering** across rules.
  Report pairwise Spearman ρ. An ordering that survives all three is a real result; one
  that doesn't is an artifact of a cutoff, and is reported as such.
- **Depth / dropout checks.** Regress per-patient `frac_double_negative` on the median
  UMIs-per-cell of that patient's malignant cells. A strong negative slope means the
  metric is measuring sequencing depth rather than biology — **a falsification test, to
  be run before the ranking is presented anywhere**. Also downsample all malignant cells
  to a common depth and recompute.
- **The truncate-all-cohorts-at-10,000-UMI sensitivity run. Owed, not optional.**
  Stages 04 and 05 both found the WashU deposit ceiling censoring a band enriched
  20-70× for `GPRC5D`, in 36 of 54 myeloma samples, biased toward this project's own
  hypothesis. Truncating every cohort at 10k costs nothing here and answers it
  directly. If the ordering survives, the metric is robust to the censoring; if it does
  not, WU1/WU2's escape fractions are partly an artifact of what the depositors removed
  and the framing must say so.
- **Expression-matched false-negative floor.** Select control genes matched to
  `GPRC5D`'s mean expression in malignant cells; their zero-fraction in those same cells
  is the technical false-negative rate the antigen call cannot possibly beat. This turns
  "GPRC5D is lowly expressed" from a hand-wave into a number.
- **Uncertainty on every patient.** Bootstrap CI over malignant cells per patient (or a
  Jeffreys interval). Cell counts vary ~15× across samples, so an unqualified rank
  ordering is not defensible. Declare a **minimum malignant-cell inclusion rule up
  front** and **report the excluded patients explicitly** — never silently drop them.
  - **Derive the minimum from the resolution the claim needs, not from a round number.**
    At n = 50 one cell *is* 2%, so 1% / 2% / 3% escape are not distinguishable and an
    ordering across that range is noise with a number attached. Work backwards from the
    smallest DN fraction the project intends to call meaningful: a 5% population yields
    ~2.5 expected DN cells at n=50, ~5 at n=100, ~10 at n=200. **≥50 is a floor, not the
    answer; expect the defensible threshold nearer 100-200.** Inspect the per-patient
    distribution, fix the threshold, *then* look at the ranking. Patients below it are
    reported descriptively rather than ranked.
  - **Bootstrap at the level of the question.** Cells within a patient are not
    independent draws, and **eight patients contribute more than one sample** —
    `27522` (six), plus `47491`, `56203`, `58408`, `59114`, `60359`, `81012` and
    `83942` (two each). Those eight are where the sample-level term has anything to
    estimate; a flat cell-level bootstrap treats sample-level batch variation as
    biological spread and reports CIs that are too narrow. **Do not resample patients
    for a per-patient CI** — a CI *for patient A* is conditional on patient A, so
    patient is fixed, not random.

    **This list is not the longitudinal list, and the two must not be conflated.**
    `83942` belongs here because the patient has two samples (`83942` in WashU 1 and
    `MMY83942` in WashU 2, merged by S1), **not** because it is a serial `_N`
    trajectory — which is why the longitudinal arm's list of seven correctly excludes
    it. Multiple samples is what nesting needs; serial timepoints is a different claim.

    | quantity | resampling scheme |
    |---|---|
    | **per-patient CI** on `frac_double_negative` | **sample → cell**, within that patient |
    | **cohort-level** inference (mean escape, regression coefficients, distributions) | **patient → sample → cell** |

    Report the flat and sample-aware per-patient intervals side by side so the narrowing
    is visible. **For the many single-sample patients this reduces to a cell bootstrap,
    which cannot see sample-level variation at all** — their CIs are optimistic in a way
    multi-sample patients' are not, and that asymmetry is stated with the ranking rather
    than buried.
- **Multi-antigen combinatorial coverage matrix.** `SLAMF7`/`FCRL5` are promoted from
  "backups" to a deliverable. For every pair and triple over {`TNFRSF17`, `GPRC5D`,
  `SLAMF7`, `FCRL5`, `CD38`, `SDC1`, `ITGB7`}, compute the uncovered fraction of each
  patient's clone. This answers what a target-strategy audience actually asks — *is
  BCMA+GPRC5D the best pair for this patient, or would BCMA+FCRL5 cover more?*

  **Report it as separate columns; do NOT collapse it into a utility score.** A weighted
  `coverage − λ · exposure` needs a principled λ, and there isn't one — the weights would
  encode a clinical judgement the data cannot supply while hiding the inputs a reader
  could otherwise disagree with. Per pair/triple, per patient:

  | column | source |
  |---|---|
  | uncovered fraction | this stage |
  | incremental gain vs. the best single target | this stage, `P(A⁻) − P(A⁻ ∩ B⁻)` |
  | co-loss enrichment | this stage, depth-conditioned |
  | **normal *marrow* expression** | stage 09 |

  The last column is normal **marrow** expression specifically — not "normal tissue",
  which this dataset cannot observe. Coverage is read against it rather than maximized
  blindly: a target covering 100% of the tumor that also hits normal marrow plasma cells
  is not a better target. Extra-marrow liabilities (GPRC5D in keratinized tissue) stay a
  cited external caveat, never a measured column.
- **The bias table** (QC methodology, above) is authored as a figure here and referenced
  from stage 12.

**BCMA/GPRC5D co-negativity enrichment — the key derived metric.**
`frac_double_negative` alone cannot distinguish two clinically different tumors. Per
patient, build the 2×2 over malignant cells:

|  | GPRC5D⁺ | GPRC5D⁻ |
|---|---|---|
| **BCMA⁺** | dual-positive | BCMA-only |
| **BCMA⁻** | GPRC5D-only | **double-negative** |

and compare observed DN against the independence expectation
`E[DN] = P(BCMA⁻) × P(GPRC5D⁻)`. Report the **co-escape enrichment ratio**
`observed / expected` with Fisher's exact and a permutation CI. This separates three
facts the single metric fuses: how often each antigen is individually absent, how many
cells are DN, and whether the *same* cells are disproportionately losing both. A patient
at 6% DN ≈ 0.3 × 0.2 has two independent partial failures; a patient at 6% DN against a
1.5% independence expectation has a coordinated antigen-low phenotype, and is the one
stage 10 investigates mechanistically. Read whatever this finds against Table S3's
per-cohort sign flip (MMRF +0.62 / WU2 +0.54 / WU1 −0.09).

**What co-escape enrichment does NOT mean.** It does **not** mean dual targeting
"doesn't help", and it does not determine whether a second binder is worth adding. The
arithmetic shows why. Adding GPRC5D to BCMA moves the uncovered fraction from
`P(BCMA⁻)` to `P(BCMA⁻ ∩ GPRC5D⁻)`. At 30% BCMA⁻ / 20% GPRC5D⁻ under independence that
is 30% → 6%. With co-loss enrichment pushing DN to 15%, it is 30% → 15% — less than
independence promised, but still halving the escape population. Enrichment measures
**how much of the two targets' expected complementarity is eroded by correlated loss**,
not whether the second target is worth adding. Use that framing everywhere.

**Incremental coverage gain — reported alongside.** Co-escape enrichment is a statement
about *biology* (is loss correlated); the clinical question is about *value* (what does
the second target buy). Different quantities, both cheap off the same 2×2:

    gain from adding GPRC5D to BCMA  =  P(BCMA⁻)   − P(BCMA⁻ ∩ GPRC5D⁻)
    gain from adding BCMA to GPRC5D  =  P(GPRC5D⁻) − P(BCMA⁻ ∩ GPRC5D⁻)

Report both per patient with CIs, as separate columns. A patient can carry high
enrichment *and* a large incremental gain — not in tension, and collapsing them into one
number would hide exactly that case. This is the quantity a single- vs. dual- vs.
sequential-target discussion actually turns on, and it generalizes to the coverage matrix.

**The null must be depth-conditioned, or this test measures library size.** Dropout is a
per-*cell* property: a shallow cell is more likely to read zero for *both* genes, so
depth heterogeneity alone produces positive BCMA⁻/GPRC5D⁻ association. A permutation that
shuffles labels freely within a patient destroys the depth↔label coupling and will report
co-escape enrichment on data with no biological co-occurrence at all — an artifact
pointing in exactly the direction the project wants to find, which is the worst kind.
- Stratify cells by depth (or `n_genes_by_counts`) within patient and **permute labels
  within stratum**; equivalently compute `E[DN]` from a per-cell independence model where
  `P_i(BCMA⁻)` and `P_i(GPRC5D⁻)` are functions of cell *i*'s own depth, summed over cells.
- **Report the unconditioned ratio next to the conditioned one.** The gap between them
  *is* the depth artifact, quantified — a more convincing exhibit than the conditioned
  number alone.

**The detection curve, and what it cannot deliver.**
`Σ_i P_i(BCMA⁻) · P_i(GPRC5D⁻)` is **circular as a correction** and must never be used as
one: multiplying the marginals assumes exactly the independence the co-escape test exists
to interrogate, so a tumor with genuinely correlated loss would be "corrected" toward the
null it violates.
- **That computation is reported once, as the depth-adjusted DN expectation under
  conditional independence** — a *technical baseline* the observed value is compared
  against, never a corrected truth. The "dropout-adjusted DN" and the "expected DN under
  depth-conditioned independence" are the same number, not two deliverables.
- **No dropout-corrected DN point estimate is produced, and none is claimed.** Dropout is
  *bounded* here — threshold band, false-negative floor, depth regression, downsampling —
  not corrected. The observed DN stays the point estimate, reported as an interval.

Still build the detection curve: fit detection probability against cell depth and gene
mean on the expression-matched control genes already selected for the floor, giving each
observed zero an approximate `P(false zero)`. It is what makes the depth-conditioned null
quantitative rather than rank-based. It does not license a corrected DN.

**A genuinely dropout-corrected DN would need a joint model, and is deferred.** The
defensible version is a latent-class model over the four true states
(B⁺G⁺ / B⁺G⁻ / B⁻G⁺ / B⁻G⁻) with per-cell detection probabilities from the curve, fit by
EM over the observed 2×2 — estimating the true joint *without* assuming independence and
yielding co-escape enrichment as a by-product rather than an input. Real statistical work,
**not** on the critical path. Filed so it is not reinvented casually, and so "bounded, not
corrected" reads as a deliberate choice rather than an oversight.

**Imputation/denoising (MAGIC, scVI, ALRA, …) is forbidden for positivity calls**, and
not as a stylistic preference: imputation manufactures low-level expression by borrowing
from neighbors, and the entire question is whether a transcript is genuinely absent.
Smoothing over the zeros erases the measurement. The detection curve models the
uncertainty instead of filling it in.

**09 — Escape robustness** (`notebooks/09_escape_robustness.ipynb`,
`src/mm_escape/bulk.py` + `robustness.py`; env: `mm-core`). Everything here exists to
answer "how do you know your escape fractions are real?"
- **Matched bulk RNA-seq — orthogonal validation of antigen *abundance*, not of the
  DN fraction.** For the **26** samples with matched bulk, **split by cohort**: MMRF bulk is
  CD138+ sorted and pairs with malignant-cell pseudobulk, WashU 1 bulk is unsorted
  BMMC and pairs with **whole-sample** pseudobulk. Pooling them would make 10 of the 26
  comparisons measure tumour burden instead of antigen abundance. Then
  correlate malignant-cell pseudobulk `TNFRSF17`/`GPRC5D` against bulk TPM (Spearman,
  concordance, residuals, and the named discordant cases). The load-bearing question
  is whether scRNA zero-rates run systematically high where bulk says the transcript
  is plainly present — which is direct, quantified evidence of dropout and feeds back
  into stage 08's false-negative floor.
  **What bulk cannot do is validate `frac_double_negative`.** Bulk destroys the joint
  single-cell distribution: a tumor that is 50% BCMA⁺GPRC5D⁻ plus 50% BCMA⁻GPRC5D⁺
  shows healthy bulk expression of *both* genes while containing zero dual-positive
  cells — and the converse misreads are equally available. Bulk constrains marginal
  abundance per gene; the *joint* distribution over cells is visible only in
  single-cell data and has no orthogonal check in this project. Phrase every output
  here as **"orthogonal validation of antigen abundance and the plausibility of
  scRNA-derived antigen-negative calls"**, never as validation of the escape fraction.
  Handle the two empty 114-byte stubs and the three ID mismatches documented in the
  Data section.
- **Normal plasma-cell antigen baseline — marrow expression context, not a safety
  axis.** Do *normal* plasma cells (from `BM*`/`ND_*`
  marrow) express BCMA and GPRC5D, and what do other marrow lineages show? This is
  real and worth having: BCMA carries broad normal-PC and B-lineage expression, and the
  malignant-vs-normal-PC contrast is what makes a coverage number interpretable rather
  than absolute. It feeds the coverage matrix in stage 08.
  **It is not a safety axis and must not be called one.** GPRC5D's decisive off-tumor
  liability is keratinized tissue — the nail, skin and taste toxicity seen with
  talquetamab — which a marrow dataset cannot observe at all, and expression is not
  toxicity. Keep three things separate in the writeup:
  **(a) tumor coverage**, **(b) normal *marrow* expression** (measured here),
  **(c) known extra-marrow liabilities** (external evidence, cited, not measured).
  A genuine target-ranking utility score of the form
  `coverage − λ · normal-tissue exposure` needs a normal-tissue atlas (GTEx/HPA or a
  normal scRNA atlas) and is filed as a future extension, not claimed from this data.
- **The label-permutation null lives at stage 08**, as the co-negativity test — it
  tests independence between the two negativities, not absence of signal, so it belongs
  with the co-escape enrichment and its depth-stratified null.

**10 — Escape subclone + phenotype** (`notebooks/10_escape_subclone_phenotype.ipynb`,
`src/mm_escape/subclone.py`; env: `mm-core`). **The project's actual scientific payoff**
rather than another robustness check.
- **Is the double-negative population structured, or scattered noise?** "3% of this
  patient's cells are double-negative" and "this patient has a pre-existing 3%
  resistant subclone" are different claims, and **only the second one predicts
  selection under therapy** — which is the entire clinical premise of the project.

  **Transcriptional clustering alone does NOT establish clonality.** A
  transcriptionally coherent group can arise from cell cycle, stress, interferon tone,
  metabolic state, sequencing depth, or sample-prep batch as easily as from a genetic
  subclone; and conversely, cells of a genuine genetic clone need not form a tidy
  transcriptional island. So the question is **DN coherence**, evaluated at three
  escalating levels with the claim escalating with it:

  | Level | Question | Method | What it licenses saying |
  |---|---|---|---|
  | **1 — enrichment** | Are DN cells non-randomly located in malignant transcriptional space? | kNN-neighborhood enrichment, Moran's I on the DN label, Leiden-subcluster Fisher enrichment, **within-patient label permutation (depth-stratified, per stage 08)** | "DN cells are non-randomly distributed" |
  | **2 — transcriptional coherence** | Do DN cells share a reproducible program? | the stage-06/10 program scores — MYC, OXPHOS, stress, IFN, UPR, antigen presentation, γ-secretase | "an escape-associated **state**" |
  | **3 — genomic coherence** | Do DN cells preferentially occupy a CNV-defined malignant subclone? | `infercnvpy` substructure from stage 07 | "an escape-associated **subclone**" |

  **Only level 3 licenses the word "subclone."** Levels 1+2 without level 3 are
  reported as an escape-associated *state* — still a real and interesting finding,
  and still the thing that separates a structured 3% from a scattered 3%, but not a
  claim about pre-existing genetic clones under selection.

  **A negative at level 3 is not evidence of absence.** Resolving CNV substructure
  *within* one patient's clone is a much harder problem than separating tumor from
  normal, and at this cohort's depth (1,521 median genes/cell in the plasma
  compartment) it will often be underpowered. Report level 3 as **supported / not
  evaluable**, with the per-patient CNV resolution stated — never as "no CNV
  subclone". Treating an underpowered null as a negative result would systematically
  understate exactly the risk this project exists to measure.

  Emit the per-patient level attained alongside the escape fraction. Runs on the
  per-patient un-integrated embedding from stage 05, never the Harmony one — and note
  that the level-1 depth-stratified permutation is the same guard as stage 08's:
  shallow cells are both more likely to be DN *and* to cluster together in
  low-dimensional space, so an unconditioned enrichment test sees depth structure and
  calls it biology.
- **Phenotype of the escape cells.** Pseudobulk differential expression
  (`pydeseq2`/`decoupler`) of double-negative vs. dual-positive malignant cells,
  aggregated with **patient as the unit of replication** — `sc-best-practices` is
  explicit that per-cell DE tests treat cells as independent replicates and badly
  inflate FDR. Then pathway/TF activity via `decoupler` (Hallmark, PROGENy,
  CollecTRI).
- **Pre-registered hypothesis: the γ-secretase axis** (`NCSTN`, `PSEN1`, `APH1A`,
  `APH1B`, `PSENEN`). γ-secretase cleaves BCMA off the cell surface, and
  γ-secretase-inhibitor + BCMA CAR-T combinations are in active clinical development
  precisely to counter it, so a γ-secretase-high escape phenotype would be directly
  actionable rather than descriptive. **Registered before looking**, so it stays a
  hypothesis test and not a post-hoc story.
- **Malignant-cell program scoring.** The stage-06 orthogonal programs (cell cycle,
  IFN, antigen presentation, UPR, hypoxia) plus two that are myeloma-specific and only
  meaningful once malignant cells are isolated: the **MYC program** (`MYC` + targets —
  MYC activation is a recognized myeloma progression event, so "is the escape population
  MYC-high?" is substantive rather than generic) and **OXPHOS** (a standard axis of
  malignant PC heterogeneity and a common covariate of proliferation). All stay
  **continuous scores, never categorical labels** — the stage-06 identity/state
  separation holds here too. A cycling MYC-high escape cell is one cell carrying three
  scores, not a new cell type.
- **TC (Translocation/Cyclin D) molecular subgroup, per patient.** Assign from
  per-patient pseudobulk over malignant cells using the genes whose dysregulation
  defines the founder event: `CCND1` (t(11;14)), `CCND3` (t(6;14)), `NSD2`/`FGFR3`
  (**`NSD2` is `WHSC1` in the older reference — depends on stage 05's symbol
  harmonization; without it this class cannot be called at all**)
  (t(4;14)), `MAF` (t(14;16)), `MAFB` (t(14;20)), `CCND2`, plus **`CKS1B` as the
  1q21-gain readout** (which also cross-checks `infercnvpy`'s CNV call on that arm from
  stage 07). Two reasons it earns its place: it is **cheap** (~8 bimodal genes off
  pseudobulk, versus reconstructing the bulk-array UAMS 7-group signatures), and it asks
  a **target-strategy question** — does dual-antigen escape risk concentrate in a
  molecular subgroup? If t(4;14) or t(11;14) patients carry systematically higher
  `frac_double_negative`, that speaks directly to who needs a different construct.
  A hypothesis the data can test, not a known result.

  **It is a transcriptional proxy for the translocation, not a detection of it.**
  `NSD2`/`FGFR3` overexpression is *consistent with* t(4;14); it is not a breakpoint
  call, and expression can be driven by other things. Every output label reads
  **"TC-like expression subtype"** or **"transcriptionally inferred TC class"** — never
  "patient has t(4;14)". This costs nothing and removes the single easiest claim in the
  project to attack. **S1 carries no cytogenetics column**, so there is nothing in this
  deposit to validate the proxy against and it stays a proxy — do not write it up as if
  a confirmation is pending.

  **Assigned per patient from pseudobulk, never per cell.** These signatures come
  from bulk arrays of purified plasma cells; per-cell assignment at 1,521 median
  genes/cell would be over-claiming. The founder translocation is clonal, so per-patient
  uniformity is the expectation — a patient splitting across two TC classes is more
  likely a doublet or patient-mapping artifact than biology, and is flagged as a QC
  signal rather than reported as heterogeneity.

  **Use it descriptively, not as a statistical stratifier.** At n ≈ 41 patients
  across ~5 TC classes, an association test is underpowered; the class belongs in
  the ranked table so a reader can see the high-escape group's composition. Run a
  formal test only where a group is large enough, and say plainly when it is not.
  **UAMS 7-group is deliberately not adopted** — it splits an already-small cohort
  into seven bins to support a test the cohort cannot power, and its signatures need
  supplement sourcing, while TC delivers most of the interpretive value for far less.

**11 — Cell-cell communication** (`notebooks/11_cellchat_liana.ipynb`,
`src/mm_escape/communication.py`; env: `mm-communication`). LIANA+ using the
CellChat-algorithm method specifically (or the full consensus rank-aggregate).
Tests whether high-escape-risk patients show weaker NK/T-cell engagement toward
malignant plasma cells.

**This stage is explicitly EXPLORATORY and is presented as such (decided
2026-08-21).** It is ninth in the scientific hierarchy near the top of this document,
not co-equal with the antigen analysis. The reason is power, not interest: n ≈ 41
patients before the stage-08 malignant-cell inclusion rule and fewer after, against a
ligand-receptor search space of hundreds of interactions with a confounder (T/NK
abundance) that is itself correlated with the predictor. That combination does not
support a confirmatory claim. Label its outputs exploratory, report effect sizes with
CIs rather than a filtered list of significant hits, and promote a finding to
foreground only if it is strong, stable across the LIANA+ methods, and survives the
abundance covariate. Left un-demoted, this is the stage that turns a focused antigen
project into a kitchen-sink scRNA project.

**Statistical design corrected from the original plan.** The original "split into
high/low tertiles and compare" approach pools cells across patients into two
groups, which is **pseudoreplication** — it treats thousands of cells from one
patient as thousands of independent observations, inflating significance
arbitrarily. Instead:
- Compute **per-patient** LR interaction scores; the patient is the unit of
  replication, n ≈ 41, not n ≈ 172,940.
- Model `frac_double_negative` as a **continuous predictor** across patients rather
  than binarizing into tertiles. This avoids an arbitrary cutoff and has strictly
  more power than discarding the middle tertile.
- **Control the obvious confounder.** High-escape patients may simply have fewer
  T/NK cells, in which case "weaker signaling" is a composition artifact, not
  immune evasion. Include cell-type abundance (from stage 06) as a covariate, and/or
  downsample to equal per-cell-type counts per patient before running LIANA+. This
  confounder is fatal to the stage's claim if unaddressed, so it is not optional.

**12 — Decision packet** (`notebooks/12_decision_packet.ipynb`; env: `mm-core`).
The final stage; consumes the output of everything upstream. Assembles:
- Escape-fraction table, annotated from S1 with ISS stage, treatment, TTPD, and
  disease stage **for WashU cohort 1 only** — MMRF and WashU 2 have no per-sample stage
  and it is *not* imputed. **S1 carries no cytogenetics**, so no karyotype column exists
  and stage 10's TC class enters as an explicitly transcriptional proxy.
- **Caterpillar plot with confidence intervals**, never a ranked bar chart — a bar
  chart implies a precision this metric does not have.
- **Risk tiers, not a rank ordering.** Printing "#1 patient 123, #2 patient 456" is
  false precision when the CIs overlap heavily, and the ordinal positions are not
  stable quantities. Assign tiers instead, and only where the evidence justifies one:
  - **Robust high escape** — DN estimate elevated, stable across the threshold band,
    sufficient malignant cells, low depth dependence.
  - **Uncertain** — wide CI, threshold-sensitive, or below the malignant-cell
    inclusion rule.
  - **Robust low escape** — consistently low across every sensitivity assumption.

  Then carry **co-escape enrichment** (stage 08) and **DN coherence level** (stage 10)
  as independent columns rather than folding them into one score, because they answer
  different questions and a patient can be high on one and low on another:

  | Patient | DN fraction | uncertainty | co-escape | DN coherence | interpretation |
  |---|---|---|---|---|---|
  | A | 8.2% | narrow | enriched | level 3 (CNV-supported) | strong baseline escape signal |
  | B | 7.5% | wide | neutral | level 1 only | uncertain |
  | C | 2.1% | narrow | enriched | level 2 (state) | small but structured population |

  Patient C may be the more interesting case despite the smaller raw fraction, and a
  bare ranking would bury that below patient B. **The ranking-stability check (stage
  08's Spearman ρ across thresholds) stays** — it is a robustness diagnostic and is
  what earns a patient the "robust" label; it is not itself the deliverable.
- The multi-antigen coverage matrix (stage 08), the DN-coherence level (stage 10),
  the co-escape enrichment ratio (stage 08), the bulk-validation correlation
  (stage 09), and the bias-direction table.
- UMAP of malignant cells by `coverage_class` faceted by patient; LIANA+
  differential interaction plot.
- **The mRNA-vs-protein limitation, stated explicitly and mechanistically.** CAR-T
  binds surface protein; this analysis measures transcript. BCMA is actively shed
  from the cell surface by γ-secretase, and GPRC5D transcript correlates
  imperfectly with surface density. This is the first question a target-strategy
  audience will ask, so it is answered in the deliverable rather than waited for.
  Calibrate against published CITE-seq/flow mRNA-protein correlations for these two
  antigens if obtainable.
- Short written interpretation with **decision rules stated in advance** (which
  patients are poor candidates for a BCMA/GPRC5D dual-target construct alone, and
  on what threshold), not fitted after seeing the ranking.

---

## Status / immediate next step

**Stages 01-05b are run and verified. Stage 06 (annotation) is the next artifact.**
`raw/` is intact at 62 samples, the five `envs/*.yml` are built with kernels registered,
and `src/mm_escape/` holds `config.py`, `gene_space.py`, `io.py`, `qc.py`,
`integration.py` and `benchmark.py`. Run output and module validation logs are in
`docs/stage-results.md`; `RESUME_HERE.md` carries exact session state.

**Supplementary Table S1 is in the repo and parsed.** The patient mapping is resolved
(41 patients / 53 in-cohort samples), the `_N` suffixes are confirmed as serial
timepoints, and clinical covariates reach every cell. **It carries no cytogenetics** —
so stage 10's TC class has nothing in this deposit to validate against and stays a
proxy, and stage 12's karyotype annotation is not available. The full account of what
S1 did and did not close is in `docs/decisions-archive.md`.

First actions, in order:
1. Write `src/mm_escape/annotation.py` and `notebooks/06_annotation.ipynb` in
   **`mm-annotation`** (needs `pip install -e . --no-deps` in that env first — it was
   only done for `mm-qc` and `mm-core`).
2. Run all three methods and compare per class against the F1 thresholds **declared in
   advance** (PlasmaCell 0.95 / T-NK-myeloid 0.90 / rest 0.85). Do not settle this
   implicitly by whichever runs first. Write
   `results/06_annotation/annotation_decision.md`.
3. The marker-coverage test is the load-bearing evidence, not the concordance numbers —
   a label set can be perfectly self-consistent and biologically wrong, and it must be
   run **on myeloma marrows specifically**, not only the donors.
4. Emit the interface contract exactly: `cell_type`, `cell_type_fine`,
   `annotation_source`, `annotation_conf`, plus `config.ANNOTATION_DECISION`. Stages
   07-12 read `cell_type` and nothing else.
5. Score the orthogonal state programs as **continuous** `obs` floats — never let them
   leak into `cell_type`.

### S1-gated additions, now unblocked

- **Within-patient longitudinal escape trajectory.** `27522_1`…`_6` is six timepoints
  from one patient, plus `47491_1/2`, `56203_1/2`, `58408_1/2`, `59114_1/4`,
  `60359_1/2`, `81012_1/2`. Confirmed serial by S1 sheet 2, so this is a real
  longitudinal arm at zero extra data cost — does escape fraction rise over time within
  a patient? Note `37692_2` and `57075_3` are lone later timepoints whose earlier draws
  were not deposited.
- **Escape fraction vs. clinical covariates.** NDMM vs. RRMM, ISS stage, treatment,
  time-to-progression. Descriptive at n ≈ 41, but "is baseline escape risk higher in
  relapsed/refractory disease?" is a real testable hypothesis on this cohort as-is.
  **Cytogenetic covariates (1q21 gain, t(4;14)) are NOT available** — S1 has no such
  column.

### Execution ordering

Chosen so the project has several presentable stopping points rather than being
all-or-nothing:

1. ~~Re-confirm `scripts/01-03` → scaffold `src/mm_escape/` + envs.~~ **Done.**
2. Stages 04-08 core path. **04, 05, 05b done.**
   **First presentable state: escape fractions with co-escape enrichment.**
3. Stage 08's defense layer (sensitivity band, dropout checks, CIs) + stage 09.
   **Second presentable state: a ranking that survives hostile questioning.**
4. Stage 10 (subclone test + phenotype). **Third: the actual scientific finding.**
5. Stage 11 (communication), then stage 12 (decision packet) last.
6. Phase 2 (GSE117156) — strictly last.

---

## Things to not re-litigate

A scan list, not an argument. Each line is settled; **the reasoning is at the stage it
binds, and every position these replaced is in `docs/decisions-archive.md`.** Entries
that merely restate the Data or Environments rules above are not repeated here.

**Infrastructure**
- No Biowulf/HPC, no SINCLAIR — ~970 MB of data, no alignment step exists.
- Five envs split by dependency-conflict risk, not one per stage; R stays quarantined in
  `mm-qc` and `mm-annotation`; `scCODA` gets its own env if used.
- Never `pip install` casually into a built env — rebuild from the yml.
- The whole analysis runs in notebooks 01-12; `scripts/01-03` are a wrapped CLI fallback,
  never reimplemented.
- Notebooks carry the analysis; `src/` carries what is reusable, testable or fiddly.
  Paired via jupytext, review on the `.py`.
- Notebooks and `results/` numbered 04-12, matching 1:1; `src/` modules named by
  function, never numbered. **Number order is execution order, no exceptions.**

**Data and gene space** — full rules under "Data" above.
- Custom `read_mtx` loader, not `scanpy.read_10x_mtx()`.
- `56203_1` is repaired and retained; zero-filling is wrong.
- Ensembl-ID join, verified and committed — **do not re-open, do not regenerate, and do
  not search the raw files for an `ENSG` column.** Intersect, never union. The alias map
  is a regression assertion only; the "drop ~52 ambiguous symbols" interim is superseded
  and **must not be implemented**.
- The patient mapping is S1's (41/53), not the naive rule's. The `_N` suffixes are serial
  timepoints, so the longitudinal arm is real.
- Ambient RNA correction is not attemptable — no unfiltered matrices are public.
- The deposit is pre-filtered per cohort; the WashU 10,000-UMI ceiling is a real
  confounder biased toward this project's hypothesis, carried as a covariate with a
  sensitivity run owed at stage 08.

**QC and integration**
- MAD thresholds re-derived per this cohort **and per cohort within it**.
- `pct_counts_in_top_20_genes` is reported but never filters — cutting on it deletes
  antigen-**positive** plasma cells.
- QC annotates, it does not delete. `scDblFinder` via an isolated `rpy2` bridge.
- Harmony on `patient_id` + `n_genes_ref` + `cohort` survived a seven-arm benchmark and
  is not re-opened. **On this dataset a standard global scIB ranking selects the wrong
  method** — the immune compartment is scored, the plasma compartment only diagnosed, and
  plasma mixing never contributes positively.
- `R²(depth ~ latent)` is the depth statistic, for rotation-invariance. BBKNN is excluded
  (graph, no embedding); scANVI/scGen are deferred, not rejected (they need stage-06
  labels, which come from the embedding under selection).
- Malignant subclustering is per-patient and un-integrated. **No integration method
  restores cells that were never deposited.**

**Annotation and cell state**
- Annotation is decided empirically, per class, against F1 thresholds declared in advance
  (0.95 / 0.90 / 0.85). Those numbers are **concordance, not accuracy**; the
  marker-coverage test is the biological evidence and can veto a class.
- `obs["cell_type"]` is the **only** load-bearing annotation output — which is what lets
  the decision be revisited without touching stages 07-12.
- Identity and state are separate axes: programs are continuous scores, never labels, and
  never leak into `cell_type`.
- No custom `celltypist` model for malignant states — it would force continuous tumor
  substructure into discrete bins and hide the intermediates.

**The metric and its defense**
- Light-chain restriction by **ratio**, not presence/absence — IG genes are the most
  ambient-contaminated in this tissue.
- `infercnvpy` is required, not optional, with the agreement rate reported as a stage
  output; a low rate halts the pipeline.
- Normal-BM samples are controls, not filler — but that baseline is **marrow expression
  context, not a safety axis**.
- Dropout is **bounded, not corrected**, and is the larger of the two biases. The headline
  metric is a bracketed interval, never a bare point estimate. `Σ P(A⁻)·P(B⁻)` is the
  independence baseline, **not** a correction; the latent-class/EM model is the honest
  correction and is deliberately deferred.
- Binary calls stay primary; **imputation is forbidden for positivity calls** — whether a
  transcript is genuinely absent is the entire question.
- Co-negativity enrichment is a first-class result and its null is **depth-stratified**;
  an unconditioned null would report enrichment from library size alone, biased toward
  this project's own hypothesis.
- Co-escape enrichment measures **eroded complementarity, not futility** — never write
  that enrichment means dual targeting "doesn't help"; report incremental gain next to it.
- Bootstrap at the level of the question: sample → cell within patient for a per-patient
  CI; patient → sample → cell for cohort-level quantities. Single-sample patients get
  optimistic CIs and this is stated.
- **Risk tiers are the deliverable; ranking stability is the diagnostic.**
- **No composite risk score, anywhere** — the inputs stay separate columns so a reader can
  disagree with them.
- Matched bulk validates **antigen abundance**, never the dual-negative fraction, which is
  a joint quantity bulk destroys by construction. The two bulk cohorts are never pooled.
- **"Subclone" requires CNV support**; level 3 reports *supported / not evaluable*, never
  "no CNV subclone".
- Pseudobulk DE with patient as replicate, never per-cell DE tests.
- TC yes, UAMS 7-group no. TC is per-patient, descriptive, and labelled a **TC-like
  expression subtype** — never a translocation call.
- The γ-secretase hypothesis is pre-registered.
- Stage 11 is exploratory and ranks ninth; patient is the unit of replication, escape
  fraction is a continuous predictor, T/NK abundance is controlled.

**Phase 2**
- GSE117156 is confirmed; GSE118900 was evaluated and rejected; He et al. 2022 has no
  data availability statement and is unusable.
- **GSE117156 must never be merged with GSE223060** (MARS-seq vs 10x) — separate
  `phase2_`-prefixed pipeline, distributional comparison only.

## Open questions to resolve during implementation

Closed questions (patient mapping, `_N` suffix meaning, MAD thresholds) have moved to
`docs/decisions-archive.md`. What remains genuinely open:

- **How to treat the WashU 10,000-UMI censoring in the headline metric.** Stage 08 owes
  the truncate-all-cohorts-at-10k sensitivity analysis. If the patient ordering survives
  it, the metric is robust to the censoring; if it does not, WU1/WU2's escape fractions
  are partly an artifact of what the depositors removed and the framing must say so.
  **Not optional — the bias points toward the hypothesis.**
- **The minimum malignant-cell inclusion threshold** for stage 08 — ≥50 is a floor, not
  the answer. Re-derive from the smallest DN fraction the project intends to call
  meaningful once the per-patient distribution is known; expect 100-200. Fixed **before**
  the ranking is looked at.
- **How many patients survive that threshold.** Post-QC cell counts vary ~15× (min 480,
  median 2,555, max 7,937 per sample); if many patients fall below the minimum the
  ranking's usable n may be well under 41, and the framing must adjust honestly rather
  than quietly.
- **The ambient-noise-floor antigen cutoff** — derived at stage 08 from T/NK/myeloid
  cells, still to be computed and documented with the population it came from.
- **Whether `infercnvpy` resolves sub-clonal CNV within a single patient's clone at this
  depth (1,521 median genes/cell in the plasma compartment)** — gates level 3 of the
  stage-10 coherence hierarchy, and therefore whether the project can use the word
  "subclone" at all. Determine per patient and report resolution; assume neither outcome.
- **How much of any observed co-escape enrichment survives depth stratification** — the
  gap between the unconditioned and conditioned ratios is itself a reportable number, and
  if the conditioned enrichment collapses to ~1 across the cohort that is a real
  (negative) finding about dual-antigen escape, not a failed analysis.
- **Whether the stage-08 ranking survives restriction to `high`-confidence malignant
  cells** (stage 07 tiers). If it does, the metric is robust to the weakest link in its
  own denominator; if not, the CNV-inconclusive cells are driving it and the framing must
  say so.
- **Whether any samples have paired scVDJ-seq** for a stronger malignant call than the
  kappa/lambda proxy — check GEO supplementary files.
- **Whether a published CITE-seq/flow calibration exists** for BCMA and GPRC5D
  mRNA-vs-surface-protein correlation in myeloma — determines whether stage 12's protein
  limitation can be quantified or only stated.
- **Per-sample disease stage for MMRF and WashU 2, and cytogenetics for anyone.** Neither
  is in S1; neither is imputed. Cytogenetic risk annotation would need a source outside
  this deposit.
