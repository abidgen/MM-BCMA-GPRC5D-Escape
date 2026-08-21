# CLAUDE.md — MM Dual-Antigen (BCMA/GPRC5D) Escape Risk Analysis (Python rebuild)

## Objective

For each multiple myeloma patient in GSE223060, quantify the fraction of malignant
plasma cells that would evade BOTH a BCMA-directed and a GPRC5D-directed CAR-T/TCE
(the "dual-antigen escape fraction"), then relate that metric to the immune
microenvironment via LIANA+/CellChat. Final output: a per-patient risk ranking
usable for a single- vs. dual- vs. sequential-target CAR-T strategy discussion.
Built as a portfolio project (Legend Biotech context — BCMA/GPRC5D are their core
CAR-T targets), no fixed deadline — code quality and defensibility of each
analytical choice matter as much as the final numbers.

Most antigen-escape literature asks a before/after question (did the antigen
disappear after treatment). This project asks a baseline question instead: how much
dual-antigen escape risk is already present in a patient's tumor before any
treatment, due to pre-existing clonal heterogeneity.

**Scope expansion, 2026-08-20.** A design review of this plan added a robustness
layer and two new pipeline stages, all adopted. The motivating observation:
`frac_double_negative` is a *fraction of zeros*, the noisiest quantity scRNA-seq
produces, and the original plan bounded only one of its two error directions. The
additions fall into four groups, each documented in place below rather than
duplicated here:
1. **Defend the metric** — threshold sensitivity band, dropout/depth checks,
   bootstrap CIs, matched bulk RNA-seq validation, normal-BM controls (stages 08, 09).
2. **New science** — is the double-negative population a real subclone or scattered
   noise, what phenotype do escape cells have, and a multi-antigen coverage matrix
   (stages 08, 10).
3. **Rigor upgrades to planned stages** — `infercnvpy` required, ratio-based
   light-chain calls, per-patient un-integrated subclustering, patient-level
   statistics in stage 11, `scCODA` for composition (stages 05-07, 11).
4. **Framing** — the mRNA-vs-protein limitation and a two-sided bias table (stage 12).

The sequence was **renumbered, not appended to**, so that stage numbers still run in
true execution order: the two new stages took 09 and 10, pushing cell-cell
communication to 11 and the decision packet to 12.

The single highest-value addition is the subclone-vs-noise test in stage 10: "3% of
cells are double-negative" and "there is a pre-existing 3% resistant subclone" are
different claims, and only the second one predicts selection under therapy — which
is the project's entire clinical premise.

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

## Data — confirmed ground truth (do not re-guess this, it's been verified)

**Source:** GSE223060 (scRNA-seq) / GSE223061 (matched bulk RNA-seq — **no longer
"optional"; it is now the stage 09 orthogonal validation set**), NCBI GEO. Paper:
"Single-Cell Discovery and Multiomic Characterization of Therapeutic Targets in
Multiple Myeloma" (Cancer Research, 2023). 53 BM samples from 41 patients per the
paper's main analysis (MMRF Immune Atlas Pilot + 2 WashU cohorts); the GEO archive
itself bundles 62 sample entries (includes normal/control BM beyond the 53 disease
samples).

**No raw FASTQ/SRA exists for this series.** Confirmed directly against NCBI. Only
processed Cell Ranger *filtered* output is public — not raw, not unfiltered
(matters for ambient-RNA correction, see below). Do not attempt to find or process
raw reads for this dataset; this is a closed question.

**Archive structure (confirmed by direct extraction):**
```
GSE223060_RAW.tar                              # outer archive, ~970 MB
  |-- <GSM_ID>_<sample_name>.tar.gz             # one per-sample archive (62 total)
        |-- (extract to) raw/samples/<GSM_ID>_<sample_name>/
              |-- <sample_name>/                # ONE EXTRA NESTING LEVEL
                    |-- barcodes.tsv             # uncompressed, no .gz
                    |-- genes.tsv                # single column, gene symbols only
                    |-- counts.mtx               # NOT "matrix.mtx"
```

**Critical implications for the Python loader:**
- Filenames are old Cell Ranger v2-style (`counts.mtx`, `genes.tsv`), not the
  `matrix.mtx`/`features.tsv` that `scanpy.read_10x_mtx()` expects — that function
  will silently fail to find these files. Load explicitly with
  `scanpy.read_mtx()` on `counts.mtx`, then attach `.obs`/`.var` from
  `barcodes.tsv`/`genes.tsv` manually, rather than relying on any 10x-directory
  auto-detection helper.
- `genes.tsv` is a **single column of gene symbols**, not the usual 2-3 column
  Ensembl+symbol file — don't write code that assumes a `gene_id` column exists.
- Data sits one directory level deeper than the GSM-named folder.

Sample naming observed: mix of `MMRF_####` (MMRF Immune Atlas Pilot cohort), bare
numeric IDs with `_N` suffixes (e.g. `27522_1` through `27522_6` — likely multiple
timepoints/fractions per patient), `MMY#####` IDs, `ND_######`, and `BM#` (likely
normal bone marrow controls). **Disease stage/cytogenetic metadata per sample is
NOT in the filenames** — pull it from the paper's Supplementary Table S1 before
building the analysis manifest.

**CRITICAL — mixed Cell Ranger references across samples (confirmed):** The 62
samples were processed against three different references, distinguishable by row
count in `genes.tsv`:
- 33538 genes — 37 samples
- 33694 genes — 24 samples
- 22184 genes — 1 sample only: `56203_1`

`56203_1` is missing `TNFRSF17` (BCMA) entirely, plus `IGLC1`/`IGLC2`/`IGLC3`. If
merged naively, every cell in that sample would read as BCMA-negative for a purely
technical reason. **Decision: exclude `56203_1`.** Patient 56203 is still fully
covered by `56203_2` (complete 33694-gene reference) — zero patient coverage lost.

**The 33538- and 33694-gene reference sets share only 22164 genes.** A union merge
across retained samples would make ~11k genes structurally zero in whole sample
cohorts — indistinguishable downstream from a true biological zero, which is
exactly the quantity this project measures. **Decision: intersect gene sets across
retained samples, never union**, before any concatenation (`anndata.concat` with
`join="inner"`, not `join="outer"`). Assert that all markers needed downstream
survive the intersection — BCMA (`TNFRSF17`), GPRC5D, backup antigens (`SLAMF7`,
`FCRL5`), the full annotation marker panel, and the light-chain genes (`IGKC`,
`IGLC1-7`) — hard-fail with the specific missing gene name(s) if any assertion
fails, rather than silently proceeding on a partial marker set.

**Verified 2026-08-20 — every target gene survives the intersection.** Spot-checked
directly against `genes.tsv` in a 33538-gene sample (`MMRF_1695`), a 33694-gene
sample (`27522_1`), and a normal-BM control (`BM4`): `TNFRSF17`, `GPRC5D`,
`SLAMF7`, `FCRL5`, `SDC1`, `CD38`, `ITGB7`, `NCSTN`, `IGKC` are present in **both**
reference builds. The 33538 build additionally carries `GPRC5D-AS1` (an antisense
transcript, absent from 33694) — harmless, dropped by the intersection, and must
**not** be substituted for `GPRC5D` anywhere. The hard assertions in
`gene_space.py` are still required as a regression guard, but are expected to pass.

**Matched bulk RNA-seq (GSE223061) — already downloaded, previously unused.**
`raw/unpacked_bulk/` holds **30 usable bulk samples**: 18 MMRF samples as
`<GSM>_<sample>_tpm.tsv.gz` (gene × TPM tables) and 12 WashU samples as
`<GSM>_<sample>.tar.gz`. Overlap with the scRNA cohort is ~28 samples — enough for
a real orthogonal check on the antigen quantification (stage 09). Gotchas confirmed
by direct inspection:
- **Two files are empty 114-byte gzip stubs and must be excluded**:
  `GSM6939104_MMRF_1505_tpm.tsv.gz`, `GSM6939120_MMRF_2259_tpm.tsv.gz`. Do not
  treat a 114-byte read as "zero expression" — it is a failed deposit.
- **Three bulk/sc sample-ID mismatches to reconcile against Supplementary Table
  S1**, alongside the patient-mapping fix below — do not guess these pairings:
  bulk `47499` vs. sc `47491_1`/`47491_2`; bulk `98433` vs. sc `MMY98423`; bulk
  `59114_2` vs. sc `59114_1`/`59114_4` (the suffix numbering does not align across
  assays, which is itself evidence the suffixes are not simple timepoint indices).

**Patient mapping is a known unresolved gap, not yet fixed in the R build either.**
A naive rule (strip a trailing `_<digits>` only when the stem is purely numeric,
e.g. `27522_1` -> `27522`) yields 47 disease patients from 57 disease samples, vs.
the paper's reported 41 patients / 53 samples — roughly six sample-name collapses
are being missed (`83942`/`MMY83942` is a likely pair). **This must be resolved
against Supplementary Table S1 before any per-patient aggregation step** — not
optional, and not yet done in either language's build.

**Ambient RNA correction (SoupX/DecontX) is not possible for this dataset.** Both
require the *unfiltered* Cell Ranger matrix (including empty droplets) to estimate
the background contamination profile; GEO only hosts the *filtered* per-sample
matrices. This is a real, stated limitation, not an oversight — see "QC
methodology" below for the mitigation actually used instead.

**Datasets explicitly evaluated and rejected as data sources:**
- **He et al. 2022** (`ctm2.757`, paired NDMM/RRMM scVDJ-seq) — confirmed by
  directly reading the paper (PMC8926895): no data availability statement, no
  accession of any kind. Not usable as a data source under any framing. A
  fabricated-looking accession (`GSE124310`) was floated for this paper at one
  point and independently verified to actually belong to an unrelated
  glioblastoma study — do not resurrect that number.
- **GSE118900** (MGUS/SMM/NDMM/RRMM staging cohort) — real and open, but rejected
  as the Phase 2 validation dataset: zero healthy controls, only 597 cells across
  15 patients (as low as 7-24 per patient), too statistically thin for a per-cell
  classification metric. `GSE117156` (below) is superior on every axis that
  matters here.

---

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

## QC methodology (upgraded from the R build's provisional fixed thresholds)

Three changes, grounded in `sc-best-practices.org`'s QC chapter
(https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html),
adapted for what this specific dataset actually allows:

**1. MAD-based (median absolute deviation) outlier filtering, not fixed cutoffs.**
Flag a cell as an outlier if it's more than 5 MADs from the median on
`log1p_total_counts`, `log1p_n_genes_by_counts`, or `pct_counts_in_top_20_genes`
(fraction of a cell's counts from its 20 highest-expressed genes — not tracked in
the R build, add it), plus a tighter, dataset-specific check on `pct_counts_mt`.
**The tutorial's exact numeric MAD counts and 8% mitochondrial cap are defaults
for their demo dataset (healthy PBMC/BMMC), not something to copy verbatim onto
bone marrow myeloma samples** — recompute and document per this cohort.

**2. Ambient RNA cannot be formally corrected (see Data section) — mitigate at the
antigen-scoring stage instead of leaving a naive `>0` positivity call in place.**
Derive an empirical background/noise floor for BCMA and GPRC5D from cell types
with no biological reason to express either gene (T cells, NK cells, myeloid) —
any nonzero signal there is ambient contamination by definition. Set the
malignant-cell antigen-positivity threshold above that floor. Document the
derived threshold value and the population it was derived from. Uncorrected
ambient contamination would tend to make truly double-negative cells register as
false-positive, **underestimating** escape fraction — a stated, real bias
direction.

**3. Dropout is the opposite-signed bias, and it is the larger one. It must be
bounded, not just mentioned.** The earlier version of this document accounted only
for ambient RNA. That was half the picture. `frac_double_negative` is a *fraction
of zeros*, which is the noisiest quantity scRNA-seq produces, and it has two error
sources pointing in opposite directions:

| Source | Effect on a cell's antigen call | Effect on `frac_double_negative` |
|---|---|---|
| Ambient RNA | true negative reads as positive | **deflates** |
| Dropout / limited sensitivity | true positive reads as negative | **inflates** |
| Residual normal PCs called malignant | dilutes with antigen-lower cells | **inflates** |
| mRNA measured, protein targeted | see stage 12 limitation | direction unknown |

Dropout matters more than usual here specifically because **`GPRC5D` is a
low-abundance GPCR transcript** and this cohort's median cell has only ~2,044
detected genes (R-build QC, 61 samples, 181,336 cells). A large share of
"GPRC5D-negative" calls will be technical zeros. The mitigations live in stage 08
(threshold sensitivity band, depth-matched checks, expression-matched
false-negative floor) and stage 09 (bulk RNA-seq cross-check); the point of this
entry is that **neither bias may be left unquantified, and the headline number is
reported as a bracketed interval, not a bare point estimate.**

Doublet detection: `scDblFinder` (R) remains the benchmarked best choice per
`sc-best-practices` (highest accuracy in the Xi & Li 2021 benchmark) — used from
Python via a deliberately isolated `rpy2` bridge (see Environment below), not
switched to a pure-Python alternative purely for language purity.

---

## Environment / tooling — modular envs per pipeline stage

Four environments, split by actual dependency-conflict risk (not one per file). Two of
them exist to quarantine R — `env-qc` for `scDblFinder`, `env-annotation` for
`SingleR` — so the pure-Python stack never carries an R dependency it doesn't need:

**`envs/env-qc.yml`** — stage 04 only (loading, QC, doublet detection). Isolates
the R/rpy2 bridge so no other environment needs to carry R at all.
```yaml
name: mm-qc
channels: [conda-forge, bioconda]
dependencies:
  - python=3.12
  - scanpy=1.11
  - anndata
  - pandas
  - numpy
  - scipy
  - seaborn
  - matplotlib
  - jupyterlab
  - ipykernel
  - r-base=4.3.3
  - rpy2=3.5.11
  - bioconductor-scdblfinder=1.16.0
  - bioconductor-singlecellexperiment
  - anndata2ri
```
(Adapted directly from `sc-best-practices`'s own published environment for this
exact chapter.)

**`envs/env-core.yml`** — stages 05-10 and 12 (integration, annotation, malignant
calling, antigen scoring/escape fraction, robustness, subclone/phenotype, decision
packet — everything except stage 04's QC and stage 11's LIANA+). Shared
dependencies, no conflicts between them.
```yaml
name: mm-core
channels: [conda-forge, bioconda]
dependencies:
  - python=3.12
  - scanpy=1.11
  - anndata
  - harmonypy
  - leidenalg
  - python-igraph
  - celltypist
  - infercnvpy
  - pydeseq2          # stage 10 pseudobulk DE (patient as unit of replication)
  - decoupler         # stage 10 pathway/TF activity (Hallmark/PROGENy/CollecTRI)
  - pandas
  - numpy
  - scipy
  - scikit-learn
  - statsmodels       # stages 08/09 regression, bootstrap CIs, confounder models
  - seaborn
  - matplotlib
  - jupyterlab
  - ipykernel
```

**`scCODA` is deliberately NOT in `env-core`.** It pulls TensorFlow, which is a
heavyweight and version-brittle dependency with no relationship to the rest of the
core stack — exactly the conflict-risk criterion these env splits exist to respect.
If the compositional analysis (stage 06) is actually run, give it
`envs/env-composition.yml` with `sccoda` + `anndata` + `ipykernel` only. Do not
destabilize `mm-core` to save an environment.

**`envs/env-communication.yml`** — stage 11 only (LIANA+). Isolated because
`liana`/`omnipath`'s dependency tree is version-sensitive and unrelated to the
core scanpy stack.
```yaml
name: mm-communication
channels: [conda-forge, bioconda]
dependencies:
  - python=3.12
  - liana
  - omnipath
  - anndata
  - scanpy
  - jupyterlab
  - ipykernel
```

**`envs/env-annotation.yml`** — stage 06 only. Isolated for the same reason `env-qc`
is: `SingleR` is R, and R stays quarantined in the environments that actually need it
rather than being pulled into `mm-core`.
```yaml
name: mm-annotation
channels: [conda-forge, bioconda]
dependencies:
  - python=3.12
  - scanpy=1.11
  - anndata
  - celltypist
  - r-base=4.3.3
  - rpy2=3.5.11
  - bioconductor-singler
  - celldex
  - anndata2ri
  - scikit-learn        # ARI / F1 for the annotation comparison
  - pandas
  - seaborn
  - matplotlib
  - jupyterlab
  - ipykernel
```
`celltypist` also remains in `env-core` — pure Python, no conflict, and convenient if
labels ever need re-deriving outside stage 06.

Register a distinct Jupyter kernel per env: `python -m ipykernel install --user
--name mm-qc` (and `mm-core`, `mm-annotation`, `mm-communication`).

**If scVI-based integration is ever considered as an alternative to Harmony**, it
gets its own separate env (`env-scvi.yml`) — not created yet, only if actually needed.

---

## Code structure — modular `src/` package + numbered thin notebooks

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
|                                # threshold sensitivity band, coverage matrix
|-- robustness.py                 # bootstrap CIs, depth/downsampling checks,
|                                  # expression-matched false-negative floor,
|                                  # label-permutation null
|-- bulk.py                         # GSE223061 TPM loading, pseudobulk-vs-bulk
|                                    # correlation (stage 09)
|-- subclone.py                       # DN neighborhood enrichment / Moran's I,
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
one unambiguous sequence from `01` to `12`, split across two directories only
because 01-03 are plain file-handling scripts with no need for a notebook:

| # | Notebook | Env | Output dir |
|---|---|---|---|
| 01-03 | (`scripts/`, not notebooks — data acquisition) | — | `raw/` |
| 04 | `notebooks/04_qc.ipynb` | `mm-qc` | `results/04_qc/` |
| 05 | `notebooks/05_integration_clustering.ipynb` | `mm-core` | `results/05_integration/` |
| 06 | `notebooks/06_annotation.ipynb` | `mm-annotation` | `results/06_annotation/` |
| 07 | `notebooks/07_malignant_calling.ipynb` | `mm-core` | `results/07_malignant/` |
| 08 | `notebooks/08_antigen_escape_fraction.ipynb` | `mm-core` | `results/08_escape_fraction/` |
| 09 | `notebooks/09_escape_robustness.ipynb` | `mm-core` | `results/09_robustness/` |
| 10 | `notebooks/10_escape_subclone_phenotype.ipynb` | `mm-core` | `results/10_subclone/` |
| 11 | `notebooks/11_cellchat_liana.ipynb` | `mm-communication` | `results/11_communication/` |
| 12 | `notebooks/12_decision_packet.ipynb` | `mm-core` | `results/12_decision_packet/` |

**Number order IS execution order** — 04 → 05 → 06 → 07 → 08 → 09 → 10 → 11 → 12,
with no exceptions. When the 2026-08-20 scope expansion added two stages, the
sequence was renumbered rather than appended to, precisely so this invariant holds:
escape robustness took 09 and subclone/phenotype took 10 (both run right after the
escape fraction they interrogate), pushing cell-cell communication to 11 and the
decision packet to 12, which consumes everything upstream and therefore stays last.
Do not reintroduce a stage whose number and run position disagree.

Each `results/NN_*/` subdirectory holds everything that stage produces (figures,
CSVs, checkpointed AnnData `.h5ad` files) — no more flat, ad hoc filenames in a
single `results/` dump like the R build had. **Phase 2's notebooks/results use a
`phase2_NN_*` prefix**, never plain `NN_*`, so the two phases are never
ambiguous on disk or in `results/`.

Notebooks are **paired with `jupytext`** (percent format) for git-diffability and
Codex review — commit the paired `.py` alongside (or instead of, with `.ipynb` in
`.gitignore`) the notebook, and strip heavy binary outputs (`nbstripout` or a
pre-commit hook) before committing either way. **All real logic lives in
`src/mm_escape/`; notebooks import and call, they don't contain the logic
themselves** — this is what makes Codex review practical (it reviews `.py` diffs,
not notebook JSON) and avoids the drift that made `mm_dual_antigen_escape_pipeline.md`
go stale once already in the R build.

`scripts/01_download_data.sh`, `02_check_files.sh`, `03_build_manifest.py` are
reused as-is from the R build — pure bash/Python, zero R or scanpy dependency,
already solved and verified. Do not rewrite; do not notebook-ify them.

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

---

## Pipeline stages (numbers match notebook filenames and result directories)

**01-03 — Data acquisition.** Unchanged, reused as-is from the R build. Status:
complete, confirmed working, all 62 samples verified.

**04 — QC + doublets** (`notebooks/04_qc.ipynb`, `src/mm_escape/qc.py`; env:
`mm-qc`). Custom loader (handles `counts.mtx`/`genes.tsv` naming). QC metrics via
`scanpy.pp.calculate_qc_metrics` including `pct_counts_in_top_20_genes`. MAD-based
outlier filtering, re-derived per this cohort. `scDblFinder` via the `rpy2`
bridge. `56203_1` excluded here. Checkpoint each sample's post-QC AnnData
individually (mirrors the R build's resumable-per-sample design).

**05 — Gene-space intersection, integration, clustering**
(`notebooks/05_integration_clustering.ipynb`, `src/mm_escape/gene_space.py` +
`integration.py`; env: `mm-core`). Intersect gene sets across retained samples
(hard-fail with specifics if required genes don't survive).
`anndata.concat(join="inner")`. Normalize, HVG, PCA, `harmonypy` keyed on
`patient_id` **with `n_genes_ref` as an additional covariate**, Leiden clustering,
UMAP. Diagnostic UMAP colored by reference version (`n_genes_ref`) to confirm the
intersection actually neutralized processing batch.

**Constrain integration's blast radius (new).** Harmony keyed on `patient_id` is
correct for the immune compartment but carries a real risk for the tumor: the
malignant clone is *patient-private by definition*, so forcing patients together
can blend genuinely distinct clones into one blob and erase the heterogeneity this
project exists to measure. Therefore:
- Integrated embedding is used for **immune-compartment annotation and clustering
  only** (stages 06, 11).
- **All malignant subclustering is done per patient, un-integrated** (stage 10
  depends on this and must not read the Harmony embedding).
- State explicitly in the notebook that **per-cell antigen calls are raw counts and
  are therefore integration-independent** — this is what contains the risk, and it
  is the answer to "did Harmony distort your escape fractions?" (it cannot; the
  calls never touch the embedding).

**06 — Annotation** (`notebooks/06_annotation.ipynb`, `src/mm_escape/annotation.py`;
env: **`mm-annotation`**). **Three methods are run and compared, and the choice is
made per class against thresholds declared in advance** — the earlier "`celltypist`
and/or marker-panel scoring" left a load-bearing methodological decision unmade in
the middle of the pipeline, where it would have been settled implicitly by whatever
ran first. Follows `sc-best-practices.org`'s annotation chapter
(https://www.sc-best-practices.org/cellular_structure/annotation.html).

**What is load-bearing here.** The comparison is weighted by what this stage actually
feeds, which is three things and only three:

| Downstream need | Labels that matter | Cost of getting it wrong |
|---|---|---|
| Stage 07 malignant calling | PlasmaCell vs. everything else | Wrong plasma-cell set → wrong denominator for `frac_double_negative` |
| Stage 08 ambient noise floor | T / NK / myeloid purity | A plasma cell leaking into the "confidently antigen-negative" population inflates the floor and biases every antigen call |
| Stage 11 confounder control | T/NK abundance, ideally with subsets | Composition artifact misread as immune evasion |

Fine-grained subtypes (CD4 Tcm vs. Tem, etc.) are a bonus only stage 11 benefits from
and **must never be the reason a method is chosen**.

**Method A — manual, marker-based**, at cluster level (not per cell; clustering
absorbs dropout, which matters at ~2,044 median genes/cell). `scanpy.tl.score_genes`
over the project panel: PlasmaCell (SDC1/CD38/MZB1/XBP1/IRF4), Bcell
(MS4A1/CD79A/CD19), Tcell (CD3D/CD3E/CD8A/CD4), NK (NCAM1/NKG7/GNLY), Myeloid
(CD14/LYZ/ITGAM), Erythroid (HBB/GYPA), HSPC (CD34/KIT). Then
`scanpy.pl.dotplot(..., standard_scale="var")` as saved assignment evidence, and
`scanpy.tl.rank_genes_groups(method="wilcoxon")` + `filter_rank_genes_groups`.
**The DE step is not optional** — it is the only thing that can surface a population
the seven-class panel does not cover (pDC, erythroid progenitors, a doublet-driven
cluster). Record ambiguous clusters as ambiguous rather than forcing them into a class.

**Method B — automated #1, `celltypist`.** Input normalized to 10,000 counts/cell then
`log1p` (CellTypist's stated requirement); run on **expression**, never the Harmony
embedding. Enumerate models with `celltypist.models.models_description()` rather than
assuming names; default pair `Immune_All_Low.pkl` + `Immune_All_High.pkl`, and
evaluate a healthy-bone-marrow model too if the installed version ships one.
`celltypist.annotate(..., majority_voting=True, over_clustering=<stage-05 leiden key>)`
— passing the existing Leiden key is what makes the methods directly comparable (same
partition, different labelings). Retain `conf_score`.

**Method C — automated #2, `SingleR`** with a hematopoietic reference. Chosen to cover
CellTypist's expected blind spot rather than duplicate its strengths: `Immune_All_*`
is immune-only, so its predictable failure is erythroid and HSPC — non-immune marrow
populations the reference cannot represent. `celldex`'s `NovershternHematopoieticData`
is built from sorted hematopoietic populations including erythroid and progenitor
compartments (`HumanPrimaryCellAtlasData` as broader fallback; verify both exist in the
installed `celldex` rather than assuming). Run with `clusters=<leiden key>` for
comparability, keeping per-cell scores for pruning diagnostics; retain `pruned.labels`
and the delta/score matrix as SingleR's own low-confidence signal.

**Two caveats about plasma cells, pointing opposite ways.** The automated references
contain *normal* plasma cells, not malignant ones, so malignant PCs will be labelled
"Plasma cells" — correct and sufficient, since separating malignant from normal is
stage 07's job. Do not expect an auto method to find the tumor or count that against
it. But conversely: because no malignant class exists in the reference, a heavily
aneuploid clone with an unusual transcriptome may be labelled something else, or split
across labels. **Check the plasma-cell marker-coverage test on myeloma marrows
specifically, not only on the normal-BM controls** — otherwise a systematic failure on
exactly the cells this project measures passes unnoticed.

**The comparison** (all artifacts to `results/06_annotation/`):
1. Confusion matrices — manual × CellTypist, manual × SingleR, **CellTypist × SingleR**
   — at cluster and cell level, plus adjusted Rand index and per-class F1/Jaccard.
   ARI alone is insufficient: a method can score well overall while failing on plasma
   cells specifically, the one class that must not be wrong. The two automated methods
   agreeing with *each other* is the strongest evidence available, since they are
   independently trained on different references.
2. **The marker-coverage test — the decisive one.** Dotplot the *manual* panel grouped
   by each *automated* method's labels. If CellTypist's T cells are CD3D/CD3E-high, its
   plasma cells MZB1/SDC1-high, its NK cells NKG7/GNLY-high, and so on, the automated
   labels already encode what the manual panel encodes and manual annotation adds only
   labor for those classes.
3. Confidence/coverage report: `conf_score` per cluster, SingleR pruned-`NA` rate and
   deltas, and the fraction of cells left unassigned or labelled outside the panel.

**The decision rule — per class, declared before looking.** Pre-declaring is the point;
otherwise "choose the best" becomes post-hoc rationalization of whichever result looks
tidier. A class goes to an automated method when its marker-coverage test passes, its
own confidence signal is not flagging the cluster, and agreement with manual clears a
pre-set bar:
- **PlasmaCell: F1 ≥ 0.95** — strictest, because it sets the metric's denominator.
- **T / NK / Myeloid: F1 ≥ 0.90** — these define stage 08's noise floor.
- **Bcell / Erythroid / HSPC: F1 ≥ 0.85** — nothing downstream is load-bearing on these.

Where both automated methods qualify, take the higher agreement and record that both
passed. Where neither qualifies, that class falls back to the manual cluster label.
Expected (to be confirmed, not assumed): immune classes and plasma cells from
CellTypist, erythroid/HSPC from SingleR or manual. **The numbers decide, not the
expectation.** Outcome written to `results/06_annotation/annotation_decision.md` with
the per-class table and the numbers behind it, so the choice is auditable.

**Interface contract — downstream stages read `cell_type` and nothing else:**
- `obs["cell_type"]` — canonical load-bearing coarse label, the seven project classes.
- `obs["cell_type_fine"]` — CellTypist fine label where available, else `NA`. Stage 11
  only; never load-bearing.
- `obs["annotation_source"]` — per cell: `celltypist` | `singler` | `manual`. Required
  under the hybrid; without it a mixed-provenance label column is untraceable.
- `obs["annotation_conf"]` — the winning method's confidence for that cell.
- `config.ANNOTATION_DECISION` — the per-class method map, so no downstream module
  branches on annotation logic.

This decoupling is deliberate: the comparison can be redone or reversed later without
touching stages 07-12.

**Orthogonal cell-state programs — continuous, never categorical.** Identity and state
are different axes: a cell has one `cell_type` but can carry several active programs at
once. Score these with `scanpy.tl.score_genes` and store as float `obs` columns —
**never collapse them into `cell_type`**:

| Program | Why it matters here |
|---|---|
| Cell cycle (`MKI67`, `TOP2A`, `PCNA`) | A proliferative escape subclone is a different risk from a quiescent one — feeds stage 10 |
| Interferon response (`ISG15`, `IFI6`, `STAT1`, `MX1`) | Immune-pressure marker; feeds stage 11's evasion question |
| Antigen presentation (`B2M`, HLA class I/II) | `B2M` loss is a documented immune-escape route in myeloma. CAR-T is MHC-independent so it does not affect the escape metric, but it is a *competing* evasion mechanism and belongs in the stage 11/12 interpretation |
| Unfolded-protein response (`XBP1`, `ATF4`, `HSPA5`, `DDIT3`) | Plasma cells are professional secretors; UPR tone is core plasma-cell biology |
| Hypoxia / stress | Standard confounder — cheap to score, expensive to discover late |

A cycling plasma cell is `cell_type == "PlasmaCell"` **plus** a high cell-cycle score,
not a separate "Cycling" identity. If any method emits a proliferation label as an
identity, remap it to PlasmaCell + score.

**Per-patient composition is a first-class output of this stage (new)**, not a
by-product: malignant-PC fraction of the marrow (tumor burden), and T/NK/myeloid
abundance. Both are needed downstream — tumor burden is context for stage 12, and
T/NK abundance is the primary confounder for stage 11's central claim. Any
composition *comparison* between groups uses **`scCODA`**, not a raw proportion
test: cell-type proportions are compositional data (they sum to 1, so one cell type
rising forces others down) and naive per-type tests on them are anticonservative.
Note `scCODA` needs its own env — see the Environment section.

**07 — Malignant plasma cell identification**
(`notebooks/07_malignant_calling.ipynb`, `src/mm_escape/malignant.py`; env:
`mm-core`). Subset to plasma cell clusters. Kappa (`IGKC`) vs. lambda
(`IGLC1-7`) restriction scoring per cell; per-patient dominant restriction class
(>90% in an involved marrow) marks malignant cells. Prefer scVDJ-seq clonotype
calls over the restriction proxy if available (check GEO supplementary files —
unconfirmed).

Three upgrades over the original plan, all motivated by the same fact: **this stage
defines the denominator of the headline metric, so its errors propagate directly
into `frac_double_negative`.**
- **Ratio-based restriction, not presence/absence.** IG transcripts are the single
  most ambient-contaminated genes in this tissue (plasma cells secrete enormous
  quantities of Ig mRNA into the droplet background — this document's own ambient
  discussion says so). A presence-based kappa-or-lambda call is therefore much
  noisier than it looks. Use the per-cell kappa:lambda **ratio**, which is robust to
  a shared additive background in a way a presence call is not.
- **`infercnvpy` is promoted from optional to REQUIRED.** Residual *normal* plasma
  cells are antigen-lower than malignant ones, so every normal PC mistakenly called
  malignant inflates the escape fraction. An independent CNV-based call is the only
  cross-check available. Use minority-restriction cells as the per-patient normal
  reference. **Report the light-chain vs. CNV agreement rate as a stage output** —
  a low agreement rate invalidates stage 08 and must halt the pipeline, not be
  noted and passed over.
- **Normal-BM negative control (new).** Run the identical calling logic on `BM2`,
  `BM4`, `BM5`, `BM6` and the `ND_*` samples. Normal marrow is polyclonal, so the
  correct output is *no malignant calls*. If the method calls a malignant clone in
  normal marrow, the method is broken and nothing downstream is trustworthy. This
  is the cheapest strong validation available for the most method-dependent step in
  the project, and it uses samples already downloaded and otherwise idle.

**08 — Antigen scoring + dual-antigen escape fraction**
(`notebooks/08_antigen_escape_fraction.ipynb`, `src/mm_escape/antigen.py` +
`robustness.py`; env: `mm-core`). Per malignant cell: positivity for BCMA
(`TNFRSF17`), GPRC5D, and backups (`SLAMF7`, `FCRL5`), using the ambient-noise-floor
threshold (not a naive `>0` call). Classify into `dual_positive`/`BCMA_only`/
`GPRC5D_only`/`double_negative`. Compute `frac_double_negative` per patient — **the
core novel metric**. **Blocked on the patient-mapping fix from Supplementary Table
S1** before this aggregation runs for real (see the S1 policy under "Status").

This stage now also carries the metric's defense. A single point estimate off a
single threshold is not a defensible deliverable:
- **Threshold sensitivity band.** Compute `frac_double_negative` under at least
  three calling rules — naive `>0`, the ambient noise floor, and a stricter floor.
  The claim being made is **not** any one number; it is the **stability of the
  patient ranking** across rules. Report pairwise Spearman ρ between the rankings.
  A ranking that survives all three thresholds is a real result; one that doesn't
  is an artifact of a cutoff choice, and must be reported as such.
- **Depth / dropout checks.** Regress per-patient `frac_double_negative` on the
  median UMIs-per-cell of that patient's malignant cells. A strong negative slope
  means the metric is measuring sequencing depth rather than biology — a
  falsification test, and it must be run before the ranking is presented anywhere.
  Additionally downsample all malignant cells to a common depth and recompute.
- **Expression-matched false-negative floor.** Select control genes matched to
  `GPRC5D`'s mean expression in malignant cells; their zero-fraction in those same
  cells is the technical false-negative rate the antigen call cannot possibly beat.
  This converts "GPRC5D is lowly expressed" from a hand-wave into a number.
- **Uncertainty on every patient.** Bootstrap CI over malignant cells per patient
  (or a Jeffreys binomial interval). Cell counts vary ~15× across samples in this
  cohort, so an unqualified rank ordering is not defensible. Declare a **minimum
  malignant-cell inclusion rule up front** (start at ≥50 cells) and **report the
  excluded patients explicitly** — never silently drop them.
- **Multi-antigen combinatorial coverage matrix (new).** `SLAMF7`/`FCRL5` are
  promoted from "backups" to a deliverable. For every pair and triple over
  {`TNFRSF17`, `GPRC5D`, `SLAMF7`, `FCRL5`, `CD38`, `SDC1`, `ITGB7`}, compute the
  uncovered fraction of each patient's clone. This answers the question a
  target-strategy audience actually asks — *is BCMA+GPRC5D the best pair for this
  patient, or would BCMA+FCRL5 cover more?* — which the two-antigen metric alone
  cannot. Coverage is traded off against normal-cell expression from stage 09, not
  maximized blindly: a target that covers 100% of the tumor and also hits normal
  tissue is not a better target.
- **The bias table** (in the QC methodology section above) is authored as a figure
  here and referenced from stage 12.

**09 — Escape robustness** (`notebooks/09_escape_robustness.ipynb`,
`src/mm_escape/bulk.py` + `robustness.py`; env: `mm-core`). New stage. Everything
here exists to answer "how do you know your escape fractions are real?"
- **Matched bulk RNA-seq validation (GSE223061).** For the ~28 samples with matched
  bulk, correlate malignant-cell pseudobulk `TNFRSF17`/`GPRC5D` against bulk TPM.
  Agreement means the antigen quantification is technically credible. Disagreement
  in the specific direction of *bulk-positive where single-cell reads zero* is
  direct, quantified evidence of dropout, and feeds back into stage 08's
  false-negative floor. This turns an already-downloaded, entirely unused dataset
  into the project's strongest technical-credibility exhibit at near-zero cost.
  Handle the two empty 114-byte stubs and the three ID mismatches documented in the
  Data section.
- **Normal plasma-cell antigen baseline.** Do *normal* plasma cells (from `BM*`/
  `ND_*` marrow) express BCMA and GPRC5D? This converts the project from
  efficacy-only to **efficacy plus on-target/off-tumor safety** — the axis that
  actually separates the two antigens in the clinic: BCMA has broad normal-PC and
  B-lineage expression, while GPRC5D is more tumor-restricted in marrow but carries
  the keratinized-tissue expression behind talquetamab's nail and skin toxicity.
  Feeds the coverage matrix's risk trade-off in stage 08.
- **Label-permutation null.** Permute antigen labels within patient and recompute
  the escape-fraction distribution, establishing what the metric looks like under
  no signal. Observed values must be well-separated from this null.

**10 — Escape subclone + phenotype** (`notebooks/10_escape_subclone_phenotype.ipynb`,
`src/mm_escape/subclone.py`; env: `mm-core`). New stage, and **the project's actual
scientific payoff** rather than another robustness check.
- **Is the double-negative population a subclone, or scattered noise?** "3% of this
  patient's cells are double-negative" and "this patient has a pre-existing 3%
  resistant subclone" are different claims, and **only the second one predicts
  selection under therapy** — which is the entire clinical premise of the project.
  They are distinguishable: per patient, test whether DN cells are clustered in
  transcriptional space (kNN-neighborhood enrichment or Moran's I on the DN label,
  and/or per-patient malignant subclustering with Fisher enrichment of DN per
  subcluster). **Random scatter is the signature of dropout; spatial clustering is
  the signature of a real subclone.** Emit a per-patient **clonality-of-escape**
  score reported alongside the escape fraction — a patient with 3% DN cells
  concentrated in one subclone is a materially different risk from a patient with
  3% scattered at random. Runs on the per-patient un-integrated embedding from
  stage 05, never the Harmony one.
- **Phenotype of the escape cells.** Pseudobulk differential expression
  (`pydeseq2`/`decoupler`) of double-negative vs. dual-positive malignant cells,
  aggregated with **patient as the unit of replication** — `sc-best-practices` is
  explicit that per-cell DE tests treat cells as independent replicates and badly
  inflate FDR. Then pathway/TF activity via `decoupler` (Hallmark, PROGENy,
  CollecTRI).
- **Pre-registered hypothesis: the γ-secretase axis** (`NCSTN`, `PSEN1`, `APH1A`,
  `APH1B`, `PSENEN`). γ-secretase cleaves BCMA off the cell surface, and
  γ-secretase-inhibitor + BCMA CAR-T combinations are in active clinical
  development precisely to counter it. A γ-secretase-high escape phenotype would be
  a directly actionable, mechanistically grounded finding rather than a descriptive
  one. Registered here, before looking, so it stays a hypothesis test and not a
  post-hoc story.
- **Malignant-cell program scoring.** Score the malignant compartment on the
  orthogonal programs defined at stage 06 (cell cycle, IFN, antigen presentation,
  UPR, hypoxia) plus two that are specifically myeloma-relevant and only meaningful
  once malignant cells are isolated:
  - **MYC program** (`MYC` + MYC targets). MYC rearrangement/activation is a
    recognized progression event in myeloma, which makes "is the escape subclone
    MYC-high?" a substantive question rather than a generic one.
  - **Oxidative/metabolic (OXPHOS)**. Standard axis of malignant plasma-cell
    heterogeneity and a common covariate of proliferation.
  All remain **continuous scores on malignant cells, never categorical labels** —
  the identity/state separation established at stage 06 holds here too. A cycling
  MYC-high escape cell is one cell carrying three scores, not a new cell type.
- **TC (Translocation/Cyclin D) molecular subgroup, per patient.** Assign from
  per-patient pseudobulk over malignant cells using the genes whose dysregulation
  defines the founder event: `CCND1` (t(11;14)), `CCND3` (t(6;14)), `NSD2`/`FGFR3`
  (t(4;14)), `MAF` (t(14;16)), `MAFB` (t(14;20)), `CCND2`, plus **`CKS1B` as the
  1q21-gain readout** (which also cross-checks `infercnvpy`'s CNV call on that arm
  from stage 07). Three reasons this earns its place:
  1. **Cheap** — ~8 bimodal genes off pseudobulk, versus reconstructing the
     bulk-array-derived UAMS 7-group signatures (Zhan et al., *Blood* 2006).
  2. **S1-independent** — gives patient stratification without Supplementary
     Table S1, which is still unresolved; and when S1 lands the two cross-validate
     (TC-inferred t(4;14) should match S1's reported t(4;14)).
  3. **Asks a target-strategy question** — does dual-antigen escape risk concentrate
     in a molecular subgroup? If t(4;14) or t(11;14) patients carry systematically
     higher `frac_double_negative`, that speaks directly to who needs a different
     construct. **A hypothesis the data can test, not a known result.**

  **Assigned per patient from pseudobulk, never per cell.** These signatures come
  from bulk arrays of purified plasma cells; per-cell assignment on ~2,044-gene cells
  would be over-claiming. The founder translocation is clonal, so per-patient
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

**Statistical design corrected from the original plan.** The original "split into
high/low tertiles and compare" approach pools cells across patients into two
groups, which is **pseudoreplication** — it treats thousands of cells from one
patient as thousands of independent observations, inflating significance
arbitrarily. Instead:
- Compute **per-patient** LR interaction scores; the patient is the unit of
  replication, n ≈ 41, not n ≈ 181,336.
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
- Ranked escape-fraction table, annotated with disease stage/cytogenetics from
  Supplementary Table S1 where available.
- **Caterpillar plot with confidence intervals**, replacing the original ranked bar
  chart — a bar chart implies a precision this metric does not have.
- The multi-antigen coverage matrix (stage 08), the clonality-of-escape score
  (stage 10), the bulk-validation correlation (stage 09), and the bias-direction
  table.
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

**No Python has been implemented yet.** The working tree is clean (R build removed,
preserved under `r-build-snapshot`), `raw/` is intact at 62 samples, and
`scripts/01-03` are in place. See `RESUME_HERE.md` for exact session state as work
proceeds.

First actions, in order:
1. Re-run `scripts/01-03` (unchanged) and confirm `raw/sample_manifest.csv` still
   comes out clean (62 samples, no INCOMPLETE entries) — should be a no-op
   confirmation, not new debugging.
2. Scaffold `src/mm_escape/` and the four `envs/*.yml` files.
3. Build `env-qc`, register its kernel, start `notebooks/04_qc.ipynb` against
   `src/mm_escape/io.py` + `qc.py` — validate the loader against 2-3 real sample
   directories before scaling to all 61.

### Supplementary Table S1 policy (decided 2026-08-20)

S1 is still not in the repo and still blocks the patient mapping. **Do not stall the
pipeline on it.** The policy is: build everything S1-independent first, running on
the provisional naive mapping, and **label every S1-dependent number as provisional
in the output itself** (filename suffix and in-figure annotation, not just a comment)
so a provisional number can never be mistaken for a final one. Attempt S1 retrieval
when stage 08's aggregation is reached. S1 gates only:
- the 47-vs-41 patient mapping (stage 08 aggregation),
- the within-patient longitudinal arm (below),
- disease stage / cytogenetic annotation (stage 12),
- the three bulk/sc ID mismatches (stage 09).

Everything else — QC, integration, annotation, malignant calling, per-cell antigen
scoring, the robustness suite, the subclone test — runs to completion without it.

### S1-gated additions (flagged provisional until S1 lands)

- **Within-patient longitudinal escape trajectory.** `27522_1` through `27522_6` is
  six samples from one patient, plus `47491_1/2`, `58408_1/2`, `59114_1/4`,
  `60359_1/2`, `81012_1/2`. If those suffixes are timepoints, this is a longitudinal
  arm at zero extra data cost — does escape fraction rise over time within a patient?
  **Do not assume the suffixes are timepoints**: they may be fractions, sorts, or
  replicates, and the bulk/sc suffix misalignment noted in the Data section is
  evidence against a naive timepoint reading. S1 settles it.
- **Escape fraction vs. clinical/genomic covariates.** NDMM vs. RRMM, ISS stage,
  1q21 gain, t(4;14). Even descriptive at n ≈ 41, "is baseline escape risk higher in
  relapsed/refractory disease?" is a real, testable hypothesis on this cohort as-is.

### Execution ordering

Chosen so the project has several presentable stopping points rather than being
all-or-nothing:

1. Re-confirm `scripts/01-03` → scaffold `src/mm_escape/` + envs.
2. Stages 04-08 core path on the provisional mapping.
   **First presentable state: a working escape-fraction ranking.**
3. Stage 08's defense layer (sensitivity band, dropout checks, CIs) + stage 09.
   **Second presentable state: a ranking that survives hostile questioning.**
4. Stage 10 (subclone test + phenotype).
   **Third: the actual scientific finding.**
5. Stage 11 (communication), then stage 12 (decision packet) last.
6. S1 retrieval → un-flag provisional numbers, add longitudinal + cytogenetics.
7. Phase 2 (GSE117156) — unchanged, and still strictly last.

---

## Things to not re-litigate (already decided, with reasoning)

- **No Biowulf/HPC needed.** Data is small (~970MB), no raw alignment step
  exists to run.
- **No SINCLAIR.** Was for raw FASTQ processing; moot, no raw reads exist.
- **Custom `read_mtx`-based loader, not `scanpy.read_10x_mtx()`.** The latter
  hardcodes filenames this archive doesn't use.
- **`56203_1` excluded, gene sets intersected (not unioned).** Full reasoning in
  Data section — don't revisit by trying to "recover" it via zero-filling.
- **Ambient RNA (SoupX/DecontX) is not attemptable** — no unfiltered matrices
  exist for this dataset. Mitigated via an empirical antigen-positivity noise
  floor instead, not left uncorrected silently.
- **MAD-based QC thresholds, re-derived per this cohort** — not a straight copy
  of `sc-best-practices`'s tutorial numbers.
- **`scDblFinder` (R) via an isolated `rpy2` bridge in `env-qc`**, not a
  pure-Python doublet-detection swap.
- **Annotation is decided empirically, per class, at stage 06** — manual +
  `celltypist` + `SingleR`, compared, with F1 thresholds declared before looking
  (PlasmaCell 0.95 / T-NK-myeloid 0.90 / rest 0.85). Not a preference, not an
  `and/or`. The decision and its numbers live in
  `results/06_annotation/annotation_decision.md`.
- **`obs["cell_type"]` (seven coarse classes) is the ONLY load-bearing annotation
  output.** Everything downstream reads it and nothing else; `cell_type_fine` is for
  stage 11's convenience and is never load-bearing. This is what lets the annotation
  decision be revisited without touching stages 07-12.
- **Identity and state are separate axes.** Cell-cycle / IFN / antigen-presentation /
  UPR / MYC / OXPHOS are **continuous scores**, never categorical labels, and never
  leak into `cell_type`. A cycling plasma cell is PlasmaCell + a high score, not a
  "Cycling" cell type.
- **No custom `celltypist` model for malignant states.** A regularized linear
  classifier forces plastic, continuous tumor substructure into discrete bins and
  hides the intermediates. Malignant identity stays at stage 07 (light-chain +
  `infercnvpy`); malignant substructure stays at stage 10 (per-patient un-integrated
  subclustering + pseudobulk DE + `decoupler`). Both are score-and-cluster, not
  classification.
- **TC molecular subgroup yes, UAMS 7-group no.** TC is ~8 bimodal genes off
  per-patient pseudobulk and is S1-independent; UAMS-7 needs bulk-array signature
  sourcing and splits n≈41 into unpowerable bins. TC is assigned **per patient**, and
  used **descriptively** — it is not a statistical stratifier at this cohort size.
- **R stays isolated in its own environments** (`env-qc` for `scDblFinder`,
  `env-annotation` for `SingleR`) — never merged into `mm-core`.
- **Malignant calling via light-chain restriction, not clustering alone** — and by
  **ratio**, not presence/absence, because IG genes are the most ambient-contaminated
  in this tissue.
- **`infercnvpy` is required, not optional**, with the agreement rate reported. This
  stage sets the metric's denominator; its errors propagate straight into the
  headline number.
- **Normal-BM samples are controls, not filler.** They validate the malignant caller
  (polyclonal marrow must yield no clone) and provide the normal-PC antigen baseline
  for the safety axis. They are not to be dropped as "not myeloma."
- **Dropout is bounded, not just mentioned.** It is the opposite-signed and larger
  counterpart to ambient RNA, and `GPRC5D` is a low-abundance transcript. The
  headline metric is reported as a bracketed interval with a threshold sensitivity
  band, never as a bare point estimate.
- **The patient ranking's stability across thresholds is the claim** — not any
  single threshold's value.
- **Matched bulk RNA-seq (GSE223061) is used, not shelved.** It was already
  downloaded and previously unused; it is the only orthogonal check available on the
  antigen quantification.
- **Stage 11 uses patient as the unit of replication and escape fraction as a
  continuous predictor** — the original tertile split was pseudoreplication. T/NK
  abundance is controlled as a confounder.
- **Pseudobulk DE with patient as replicate, never per-cell DE tests** (stage 10) —
  per-cell tests inflate FDR by treating cells as independent replicates.
- **Malignant subclustering is per-patient and un-integrated.** Harmony is for the
  immune compartment; the malignant clone is patient-private and must not be blended
  across patients. Per-cell antigen calls are raw counts and are unaffected by
  integration either way.
- **`scCODA` gets its own env if used** — it pulls TensorFlow and must not
  destabilize `mm-core`.
- **The γ-secretase hypothesis is pre-registered** (stage 10), so it stays a
  hypothesis test rather than a post-hoc narrative.
- **Number order is execution order, with no exceptions** (04 → 12). The
  2026-08-20 scope expansion renumbered the sequence rather than appending new
  stages out of position, specifically to preserve this. Never add a stage whose
  number and run position disagree.
- **GSE117156 is the confirmed Phase 2 dataset; GSE118900 was evaluated and
  rejected** (no healthy controls, too few cells per patient). He et al. 2022
  has no data availability statement at all and is unusable as a data source.
- **GSE117156 must never be merged with GSE223060** (platform mismatch) — runs
  as a fully separate, `phase2_`-prefixed pipeline; comparison is
  distributional, not a merged ranking.
- **Logic lives in `src/mm_escape/`, notebooks are thin orchestration**, paired
  via `jupytext`.
- **Four environments split by actual dependency-conflict risk**
  (`mm-qc`/`mm-core`/`mm-annotation`/`mm-communication`), not one per pipeline stage.
  Two of them exist purely to quarantine R.
- **Notebooks and `results/` subdirectories are numbered `04`-`12`, matching
  1:1, continuing straight on from `scripts/01-03`.** `src/mm_escape/` modules
  are named by function, not numbered — different rule for a library vs. a
  pipeline sequence, on purpose.

## Open questions to resolve during implementation

- **Patient mapping is unresolved** — 47 vs. 41 patients, needs Supplementary
  Table S1. Blocks stage 08's per-patient aggregation specifically. See the S1
  policy above: proceed provisionally, label provisional output as such.
- Whether any samples have paired scVDJ-seq for a stronger malignant-cell call
  than the kappa/lambda proxy — check GEO supplementary files.
- Disease stage (NDMM/RRMM/normal) and cytogenetic risk annotation per sample —
  also needs Supplementary Table S1.
- Exact MAD thresholds and the derived ambient-noise-floor antigen cutoff — both
  need to be computed against the real cohort and documented with their values.
- **What the `_N` sample suffixes actually mean** (timepoint vs. fraction vs. sort
  vs. replicate). Gates the longitudinal arm. The bulk/sc suffix misalignment
  (bulk `59114_2` vs. sc `59114_1`/`59114_4`) is evidence against the naive
  timepoint reading — needs S1, do not assume.
- **The minimum malignant-cell inclusion threshold** for stage 08 — starts at ≥50
  cells, but should be re-derived once the per-patient malignant-cell distribution
  is known, and fixed before the ranking is looked at.
- **How many patients survive that threshold.** Cohort cell counts vary ~15×
  (min 480, median 2,555, max 7,937 cells/sample post-QC); if a large share of
  patients fall below the malignant-cell minimum, the ranking's usable n may be
  well under 41 and the framing must adjust honestly rather than quietly.
- **Whether a published CITE-seq/flow calibration exists** for BCMA and GPRC5D
  mRNA-vs-surface-protein correlation in myeloma — determines whether stage 12's
  protein limitation can be quantified or only stated.
