# CLAUDE.md — MM Dual-Antigen (BCMA/GPRC5D) Escape Risk Analysis (Python rebuild)

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

**Second design review, 2026-08-21 — five overclaims corrected, one analysis added.**
A second review read the plan as a whole and found no method errors, but flagged that
several *claims* outrun what the data can support. All five corrections are adopted,
plus two caveats the review itself did not raise. Each is documented in place below:

1. **Clustering ≠ clonality (stage 10).** Transcriptional coherence of the DN
   population has many possible causes — cell cycle, stress, IFN, depth, batch — not
   only a genetic subclone; and a genuine genetic subclone need not form a tidy
   transcriptional island. `clonality-of-escape` is replaced by a three-level
   **DN-coherence hierarchy** in which the word *subclone* requires CNV support.
2. **Bulk RNA-seq validates antigen *abundance*, not the DN fraction (stage 09).**
   Bulk destroys the joint single-cell distribution: a tumor that is 50%
   BCMA⁺GPRC5D⁻ plus 50% BCMA⁻GPRC5D⁺ looks, in bulk, like both antigens are well
   expressed, while containing *zero* dual-positive cells. Reframed accordingly.
3. **A probabilistic dropout-adjusted DN estimate joins the binary one (stage 08).**
   The binary call stays primary and is what stage 10 consumes; imputation/denoising
   is still forbidden for positivity calls.
4. **The "label-permutation null" was mis-specified (moves stage 09 → 08).** Permuting
   antigen labels within patient does not test "no signal" — it tests *independence*
   between BCMA-negativity and GPRC5D-negativity. That is a better question than the
   one it was written to answer, and it becomes the added analysis below.
5. **Normal-BM controls give marrow expression context, not a safety axis (stage 09).**
   GPRC5D's clinically decisive off-tumor site is keratinized tissue, which a bone
   marrow dataset cannot see at all.

**The added analysis — BCMA/GPRC5D co-negativity enrichment (stage 08).** Per patient,
build the 2×2 of BCMA± × GPRC5D± over malignant cells and test whether double-negativity
exceeds the product of the two marginal negative rates. Two patients with the same 6% DN
fraction are different clinical propositions if one is 6% ≈ 0.3 × 0.2 (independent
antigen heterogeneity, two separate partial failures) and the other is 6% at 4× the
independence expectation (the *same* cells suppressing both targets). Only the second is
a coordinated antigen-low phenotype, which is precisely what a dual-target construct
cannot escape by adding a second binder — and it is what stage 10 then goes looking for
a mechanism behind. This is the sharpest question the dataset can answer and it costs
one contingency table per patient.

**Two caveats the review did not raise, both load-bearing:**
- **The independence null must be depth-conditioned, or "co-escape enrichment" is a
  library-size artifact.** Dropout is a per-*cell* property: a shallow cell is more
  likely to read zero for *both* genes, manufacturing exactly the positive association
  the test is looking for. A naive within-patient permutation destroys the depth↔label
  coupling and will therefore report enrichment on pure noise, biased in the project's
  own direction of interest. The null is permuted **within depth strata** (or computed
  against a per-cell independence model where each cell's P(BCMA⁻) and P(GPRC5D⁻) are
  functions of its own depth). This is not a refinement; the unconditioned test is
  simply wrong here.
- **`infercnvpy` may not resolve *sub*clonal CNV, and a null there is not evidence of
  absence.** Separating tumor from normal by CNV is far easier than resolving structure
  *within* a single patient's clone at ~2,044 genes/cell. Level 3 of the coherence
  hierarchy therefore reports **supported / not evaluable**, never "no CNV subclone",
  and per-patient CNV resolution is stated rather than assumed.

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

**CRITICAL — the two references also use different HGNC symbol vintages, so a naive
symbol intersection silently drops genes present in BOTH builds.** Confirmed
2026-08-20 by direct comparison, and caught by the stage-03 notebook's required-gene
assertions rather than by inspection — which is the argument for keeping those
assertions:

| 33538 (newer symbols) | 33694 (older symbols) | consequence if unharmonized |
|---|---|---|
| `NSD2`    | `WHSC1`   | **t(4;14) becomes uncallable** — the highest-risk MM translocation, and the `MS` class of the TC subgrouping at stage 10 |
| `TENT5C`  | `FAM46C`  | loses a recurrently-deleted MM tumour suppressor (1p12) |
| `NSD3`    | `WHSC1L1` | a **different gene** from NSD2 — do not conflate the two |
| `ATP5F1A` | `ATP5A1`  | OXPHOS program member (stage 10) |

**Decision (superseded 2026-08-21 — see "Ensembl-ID reconstruction" below): join on
Ensembl ID, not on symbols at all.** Symbol canonicalization recovered these four
(22,164 → 22,168 genes); the ID join recovers 32,991 and makes the alias map a
regression assertion rather than the mechanism. The table above still matters — it is
what the assertions check and why they exist.
The four above are the ones this project's gene lists depend on; a fuller HGNC-based
reconciliation is worth doing in `gene_space.py` since the same drift certainly affects
genes outside the required set. **Never assume a missing required gene is biologically
absent — check for a legacy symbol first.**

Note also that the ~11.4k/11.5k symbols unique to each build are dominated by
**annotation-version noise** — versioned clone identifiers such as `AC000032.1` vs
`AC000032.2` — not genuinely absent genes. So 22,164 *understates* the recoverable gene
space. That is tolerable for lncRNA/clone entries (HVG selection never reaches them) but
was not tolerable for `NSD2`, which is exactly why the assertions exist.

### Ensembl-ID reconstruction — SOLVED and verified 2026-08-21

The deposit has no ID column: `genes.tsv` is symbol-only and contains **zero `ENSG`
strings across all 62 samples**. Do not go looking for one. But the IDs are fully
**reconstructible**, and this was verified end-to-end — it is the method, not a plan.

There are only **three distinct `genes.tsv` files in the cohort** (checksum-verified,
byte-identical within group), each a positional dump of a public reference:

| rows | samples | reference | source GTF |
|---|---|---|---|
| 33538 | 37 | `refdata-cellranger-GRCh38-3.0.0` | Ensembl 93 |
| 33694 | 24 | `refdata-cellranger-GRCh38-1.2.0` | Ensembl 84 |
| 22184 | 1 | — | `56203_1`, excluded on other grounds |

**Reconstruction recipe (reproduces both files exactly):**
1. Fetch the Ensembl GTF (release 93 / 84) — ~41-44 MB each, not the 11 GB 10x tarball.
2. Take `feature == "gene"` rows **in GTF order**, keep first occurrence per `gene_id`,
   filtered to 10x's `mkgtf` biotype list (`protein_coding`, `lincRNA`, `antisense`,
   the 8 `IG_*` and 6 `TR_*` classes). This yields exactly 33538 / 33694 rows.
3. Apply the depositor's two transforms, in order: **`gsub("_", "-")`** then
   **R `make.unique`**. Both are Seurat artifacts — the deposited files were written
   through Seurat, which is what produced the `TBCE`/`TBCE.1` suffixes *and* the
   `RP11-442N24--B` (Ensembl writes `__`) spellings.
4. **Assert the result equals the deposited column position-for-position.** This is
   what makes the reconstruction self-certifying: a wrong biotype filter or a wrong
   release changes the row count or the order, and the assertion fails loudly.

**Result: 0 mismatches / 33538 rows and 0 mismatches / 33694 rows.** The mapping is
committed at `resources/gene_space/` (~1 MB gzipped, three TSVs: per-build
`row_index → deposited_symbol → ensembl_id`, plus the intersection table). It does not
need regenerating; the GTFs do not need re-downloading.

**The payoff is large and was not anticipated by either review:**

| join key | genes retained |
|---|---|
| raw symbols | 22,164 |
| symbols + 4-gene alias map | 22,168 |
| **Ensembl IDs** | **32,991** |

**+10,827 genes — a ~49% larger gene space**, because 11,140 intersected IDs carry a
*different symbol* in each build and were silently invisible to a symbol join. The
earlier note that the unique-symbol remainder was "dominated by annotation-version
noise" was too optimistic: most of it is real, recoverable genes.

**It also proves the mis-pairing risk was real, not theoretical:**

```
33538 build:  TBCE -> ENSG00000285053
33694 build:  TBCE -> ENSG00000116957     <- a DIFFERENT annotation entry
```

A symbol join silently merges those two rows. Nine canonical symbols remain duplicated
*within* the ID intersection (`COG8`, `CYB561D2`, `EMG1`, `LINC01238`, `LINC01505`,
`MATR3`, `PINX1`, `RGS5`, `TMSB15B`) — disambiguate as `SYMBOL__ENSGxxxxxxxxxxx`, never
by dropping or by `var_names_make_unique()`, which would re-introduce exactly the
positional mangling this section undoes.

**Implemented and tested 2026-08-21** in `src/mm_escape/gene_space.py` (with a partial
`config.py` carrying only the gene-space constants; the rest is filled in by the stages
that derive it). Public API: `attach_ensembl_ids` (per sample, pre-concat, verifies
position-for-position) -> `intersect_gene_space` (IDs only, refuses symbol-keyed input)
-> `to_canonical_symbols` (once, post-merge) -> `assert_required_genes`.
`rebuild_gene_map_from_gtf` regenerates the committed map and re-verifies it. Verified
against the real files: 33538+33694 -> 32,991 genes, 11,140 drifted symbols joined
correctly, all 64 required genes present, and four failure paths confirmed to raise —
reordered `var_names`, symbol-keyed intersection, a missing required gene (with the
legacy-symbol hint), and `var_names_make_unique()`-style mangling.

**Decision — `gene_space.py` joins on Ensembl ID; the alias map is demoted to an
assertion.** All 21 spot-checked required genes survive the ID intersection. The
four-gene alias dictionary (`WHSC1`→`NSD2` etc.) is kept **only** as a regression
assertion, not as the harmonization mechanism — it addressed 4 of 11,140 drifted
symbols. The earlier "drop the ~52 `make.unique`-ambiguous symbols" interim is
**superseded and must not be implemented**: those genes are now resolved correctly
rather than discarded.

**Index convention.** Use `ensembl_id` as `var_names` **through the merge only** — that
is where identity is load-bearing. Once the object is a single harmonized matrix the
mis-pairing risk is gone, so switch `var_names` to the canonical Ensembl-93 symbol
(with the 9 collisions suffixed as above) and retain `var["ensembl_id"]`,
`var["symbol_33538"]`, `var["symbol_33694"]`. Reason: every downstream consumer is
symbol-native — `score_genes` marker panels, dotplots, `celltypist`, `decoupler`,
`liana` — and keeping IDs as the index buys no further correctness while making every
figure unreadable and every gene-set call a translation step.

**AnnData, not MuData.** MuData exists for multiple modalities over a shared `obs`
axis (RNA + ADT + ATAC). This project has one cell-level modality. The stage-09 bulk
RNA-seq is *sample*-level, so it shares no `obs` axis with the single-cell object and
belongs as a plain DataFrame joined on `sample_id`, not as a MuData modality. Adding
MuData here buys nothing and costs compatibility across the scanpy stack.

**The 33538- and 33694-gene reference sets share only 22164 genes** *on raw symbols*
(22,168 after harmonization). A union merge
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
`raw/unpacked_bulk/` holds **29 usable bulk samples** (inventory corrected
2026-08-21 by direct count in `notebooks/01_download_data.py`, which now asserts
these numbers): **18** MMRF samples as `<GSM>_<sample>_tpm.tsv.gz` (gene × TPM
tables, GSM6939103-120) of which 2 are empty stubs, plus **13** WashU samples as
`<GSM>_<sample>.tar.gz` (GSM6939090-102, all 4.5-5.4 MB). The earlier "18 + 12 = 30
usable" was wrong twice over — the WashU count was 12 and the total did not subtract
the stubs. Correct arithmetic: (18 - 2) + 13 = 29.
Overlap with the scRNA cohort is ~28 samples — **this figure is inherited and not
yet verified**; it depends on the three ID mismatches below and on the S1 patient
mapping, so recompute it at stage 09 rather than quoting it. Still ample for a real
orthogonal check on the antigen quantification. Gotchas confirmed by direct
inspection:
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

**The 47/57 above counts the four `ND_*` samples as disease, and most of the gap
disappears if they are not (found 2026-08-24 while building `io.py`).** This document
already treats `ND_*` as normal-BM controls in two places — the sample-naming note
above ("`ND_######`, and `BM#` (likely normal bone marrow controls)") and stage 07's
negative control, which runs the malignant caller on "`BM2`, `BM4`, `BM5`, `BM6` and
the `ND_*` samples". Counting them consistently as controls:

| | samples | naive patients |
|---|---|---|
| `ND_*` counted as disease (the inherited figure) | 57 | 47 |
| **`ND_*` counted as controls** (+ `56203_1` excluded) | **53** | **43** |

53 is the paper's disease-sample count **exactly**, and 43 is two collapses from its
41 rather than six — with `83942`/`MMY83942` the one obvious remaining pair. The four
names are `ND_083017`, `ND_090617`, `ND_170531`, `ND_170607`: the suffixes are
collection **dates** (MMDDYY / YYMMDD), not patient identifiers, which is how donor
samples are usually labelled and is weak independent support.

**This is a hypothesis, not a resolution — `ND` could also read "newly diagnosed",
and the deposit does not say.** S1 settles it. `io.py` therefore classifies `ND_*` as
`sample_type == "normal_bm"` with `sample_type_certain == False`, so the uncertainty
is carried in the data rather than in a comment. Note the 53 also depends on dropping
`56203_1`, which the paper had no reason to exclude for *our* BCMA-reference reason —
so the exact sample-count match is suggestive, not proof.

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

**Written and built 2026-08-21** as `envs/env-{qc,core,annotation,communication}.yml`,
plus `envs/env-composition.yml` which is **written but deliberately not built** — per
the scCODA note above, it is created on demand only if stage 06's compositional
comparison is actually run. Three deviations from the specs above, all forced by what
the channels actually carry (verified against `conda-forge`/`bioconda`, not assumed):

| spec said | reality | resolution |
|---|---|---|
| `infercnvpy` (conda) | **not packaged on any conda channel** | installed via `pip:` inside `env-core.yml` — still REQUIRED, not optional |
| `celldex` | no such conda package | the R package is `bioconductor-celldex` |
| `decoupler` (conda) | bioconda is stuck at **1.5.0 (2023)**, pins `numpy<2`, and fails to import under numba 0.67 | `pip: decoupler==2.2.0` in `env-core.yml` |
| `liana` (conda) | bioconda recipe **under-declares**: `import liana` fails outright | pin `pydantic` + `mudata<0.4` in `env-communication.yml` |
| (unstated) | notebooks are jupytext-paired in every env | `jupytext` added to all four |

Every pinned version in the specs does exist and was confirmed before building:
`scanpy=1.11`, `r-base=4.3.3`, `rpy2=3.5.11`, `bioconductor-scdblfinder=1.16.0`.

**`decoupler` 2.x is an API rewrite** (`dc.mt.*` / `dc.op.*`, not `dc.run_mlm`). Stage 10
must be written against 2.x — do **not** follow 1.x tutorials, including
`sc-best-practices`'s, without checking the call names. Downgrading is not an option:
1.x pins `numpy<2`, which collides with `pydeseq2`/`scipy`/`zarr` in the same env.

**Never `pip install` into these envs casually.** Installing `decoupler==1.8.0` during
setup silently downgraded `numpy` 2.5.2 → 1.26.4 *and* `numba`, breaking `scanpy`,
`scipy`, `pydeseq2` and `zarr` at once; the repair was to delete and recreate the env
from the yml. The two legitimate pip entries (`infercnvpy`, `decoupler`) live in the yml
so a rebuild reproduces them. If a pip install is unavoidable, rebuild the env
afterwards rather than patching it.

**Verified 2026-08-21 — all four envs import their key packages, both R bridges work:**

| env | verified |
|---|---|
| `mm-qc` | scanpy 1.11.5, rpy2 3.5.11, anndata2ri 2.0.1; **R 4.3.3 + scDblFinder 1.16.0 + SingleCellExperiment 1.24.0**, `scDblFinder()` callable |
| `mm-core` | scanpy 1.11.5, numpy 2.5.2, harmonypy 2.0.0, celltypist 1.7.1, infercnvpy 0.6.1, pydeseq2 0.5.4, decoupler 2.2.0; `cnv.tl.infercnv` + `run_harmony` callable |
| `mm-annotation` | celltypist 1.7.1, rpy2 3.5.11; **R 4.3.3 + SingleR 2.4.0 + celldex 1.12.0**, both `NovershternHematopoieticData` and `HumanPrimaryCellAtlasData` present |
| `mm-communication` | liana 1.8.1 with **both `liana.mt.cellchat` and `liana.mt.rank_aggregate`**, omnipath 1.0.12 |

Kernels registered for all four (`mm-qc`, `mm-core`, `mm-annotation`, `mm-communication`),
each kernelspec confirmed to point at its own interpreter.

**Note `mm-communication` deliberately runs a different scanpy/anndata** (1.12.3 /
0.12.19) from the other three (1.11.5 / 0.13.2) — the spec leaves scanpy unpinned there,
and pinning it would over-constrain liana's already version-sensitive tree. Stage 11
reads `.h5ad` written by `mm-core`, so the two only meet on disk; if a forward/backward
`.h5ad` compatibility problem ever appears, that version gap is the first place to look.

**If scVI-based integration is ever considered as an alternative to Harmony**, it
gets its own separate env (`env-scvi.yml`) — not created yet, only if actually needed.

---

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
Codex review. The `.ipynb` is gitignored and generated from the committed `.py` — see
`.gitignore`, which documents how to flip that if you'd rather commit notebooks with
`nbstripout`.

**Division of labour between notebook and package.** Notebooks are the analysis: they
carry the narrative, the plots, the intermediate inspection, and the reasoning a reader
steps through. `src/mm_escape/` holds the parts that earn being a library —
**reusable** (used by more than one stage), **testable** (worth asserting on
independently), or **fiddly** (the `read_mtx` loader, symbol harmonization, the
noise-floor derivation). The test is reuse and testability, not line count: a
single-use plotting call belongs in the notebook, and duplicating a threshold
calculation across three notebooks does not.

Two things this protects, and the reason the split exists at all: Codex reviews `.py`
diffs rather than notebook JSON, and logic with one home cannot drift between copies —
which is exactly what made `mm_dual_antigen_escape_pipeline.md` go stale once already
during the R build.

**The whole analysis runs in notebooks, stages 01 through 12.** Decided 2026-08-20:
every stage is openable and steppable, with nothing hidden behind a bare CLI script.

`scripts/01_download_data.sh`, `02_check_files.sh`, `03_build_manifest.py` are
**retained as a CLI fallback**, not deleted — they are already solved and verified, they
work headlessly (useful for a fresh clone, a remote box, or CI), and the stage 01-03
notebooks call into them rather than reimplementing them. Do not rewrite the scripts;
the CLI remains a valid path and its output is the contract.

**Notebooks 01-03 wrap the scripts; they must not duplicate their logic.**
`notebooks/03_build_manifest.py` is the reference pattern: it imports `build_manifest()`
from `scripts/03_build_manifest.py` by file path and writes the manifest on the script's
exact column schema, so the two produce byte-identical output (verified both
directions). It adds what a notebook is *for* — the Cell Ranger reference split, a
symbol-harmonized intersection preview, and the required-gene assertions that caught the
`NSD2` drift. Two implementation traps it documents, both hit for real:
- Under a Jupyter kernel `sys.argv[1]` is `-f`, so the script's module-level `RAW_DIR`
  evaluates to `Path('-f')`. **Always pass paths explicitly**; never use `mf.RAW_DIR`.
- The manifest must hold **repo-root-relative** paths (it is committed and read on other
  machines), so the notebook normalizes the absolute paths it must pass in.

Its `n_genes_ref` diagnostic is intentionally **not** persisted: adding a column would
break schema parity with the script.

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

**01-03 — Data acquisition** (`notebooks/01_download_data.ipynb`,
`02_check_files.ipynb`, `03_build_manifest.ipynb`; env: `mm-qc`). Notebooks wrapping the
verified `scripts/01-03`, which remain as a CLI fallback. Status: **all three written,
executed and passing (01/02 added 2026-08-21)** — 62/62 `triplet-ok`, manifest
byte-identical via both paths, all assertions green.

Scripts 01 and 02 are **bash**, so their notebooks `subprocess` out to them rather than
importing a function as notebook 03 does — that is what "wrap, don't duplicate" means
for a shell script. Notebook 02 parses the script's classification table into a
DataFrame for assertions, but never re-derives the classification itself.

Each adds what a notebook is for. **01**: a pre-flight stating whether the run will
download or skip (so a 1 GB transfer is never started by accident), post-conditions on
the 62 sample dirs and their `.extraction_complete` markers, and the bulk inventory
assertion that caught the 12-vs-13 WashU error above. **02**: assertions on the
`triplet-ok` classification, the reference split, a manifest-vs-disk cross-check, and
the **checksum finding** — 62 `genes.tsv` files, only **3 distinct checksums**, which is
the fact that makes the Ensembl-ID reconstruction in `gene_space.py` possible at all. It
asserts exactly 3, so a new reference build fails loudly instead of being merged
silently. **03** reports the reference split and runs the required-gene assertions
(which is how the `NSD2`/`WHSC1` symbol drift was found).

Note notebook 03's kernelspec said `mm-dual-antigen` (an R-build leftover) and was
repointed to `mm-qc` on 2026-08-21, matching the env this stage actually declares.

**04 — QC + doublets** (`notebooks/04_qc.ipynb`, `src/mm_escape/qc.py`; env:
`mm-qc`). Custom loader (handles `counts.mtx`/`genes.tsv` naming). QC metrics via
`scanpy.pp.calculate_qc_metrics` including `pct_counts_in_top_20_genes`. MAD-based
outlier filtering, re-derived per this cohort. `scDblFinder` via the `rpy2`
bridge. `56203_1` excluded here. Checkpoint each sample's post-QC AnnData
individually (mirrors the R build's resumable-per-sample design).

**05 — Gene-space intersection, integration, clustering**
(`notebooks/05_integration_clustering.ipynb`, `src/mm_escape/gene_space.py` +
`integration.py`; env: `mm-core`). **Canonicalize gene symbols first, then** intersect
gene sets across retained samples (hard-fail with specifics if required genes don't
survive). The harmonization step is not optional — see the Data section's symbol-drift
table; without it `NSD2` is dropped and stage 10 loses t(4;14) entirely. Report how many
genes the harmonization recovers, so a regression here is visible rather than silent.
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

**These are concordance scores, not validation accuracy (clarified 2026-08-21).** The
manual annotation is not a ground-truth labelled dataset — it is a third opinion derived
from the same expression matrix. So "F1 against manual" measures *agreement*, not
correctness, and CellTypist agreeing with SingleR is agreement between two references
that share canonical marker biology and may share its blind spots. Report and label all
of it as **concordance**. The load-bearing evidence of *biological* validity is the
marker-coverage test (#2 above), not the concordance numbers — a label set can be
perfectly self-consistent and biologically wrong. The stage's headline figure is
therefore three panels in this order: **label concordance → marker-expression validity →
uncertainty / unassigned rate**, with concordance first because it is the weakest of
the three, not the strongest.

**The decision rule — per class, declared before looking.** Pre-declaring is the point;
otherwise "choose the best" becomes post-hoc rationalization of whichever result looks
tidier. A class goes to an automated method when its marker-coverage test passes, its
own confidence signal is not flagging the cluster, and concordance with manual clears a
pre-set bar (F1 used as the concordance statistic, per above):
- **PlasmaCell: F1 ≥ 0.95** — strictest, because it sets the metric's denominator.
- **T / NK / Myeloid: F1 ≥ 0.90** — these define stage 08's noise floor.
- **Bcell / Erythroid / HSPC: F1 ≥ 0.85** — nothing downstream is load-bearing on these.

Where both automated methods qualify, take the higher concordance and record that both
passed. Where neither qualifies, that class falls back to the manual cluster label.
**A failed marker-coverage test vetoes a class regardless of concordance** — high
agreement on a biologically unsupported label is agreement on an error.
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

**Malignant-evidence tiers, not a binary call (added 2026-08-21).** With two
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
agreement rate alone. Note this is the same distinction the 2026-08-21 review makes at
stage 10: `uncertain` here is a *reported quantity*, never a silent drop, and
`probable` from CNV being **not evaluable** is not the same as CNV being negative.

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
  malignant-cell inclusion rule up front** and **report the excluded patients
  explicitly** — never silently drop them.
  - **Derive the minimum from the resolution the claim needs, not from a round
    number (revised 2026-08-21).** At n = 50 malignant cells one cell *is* 2%, so
    1% / 2% / 3% escape are not distinguishable and a rank ordering across that
    range is noise with a number attached. Work backwards from the smallest DN
    fraction the project intends to call meaningful: a 5% population yields ~2.5
    expected DN cells at n=50, ~5 at n=100, ~10 at n=200. **≥50 is a floor, not the
    answer; expect the defensible threshold to land nearer 100-200.** The rule stays
    as written — inspect the per-patient malignant-cell distribution, fix the
    threshold, *then* look at the ranking — and patients below it are reported
    descriptively rather than ranked.
  - **Bootstrap hierarchically — but at the right level (corrected 2026-08-21).**
    Cells within a patient are not independent draws, and several patients contribute
    multiple samples (`27522_1`…`_6`, `47491_1/2`, and the rest); a flat cell-level
    bootstrap treats sample-level batch variation as biological spread and reports CIs
    that are too narrow. An earlier draft said to resample **patient → sample → cell**
    for the per-patient CI. **That is the wrong level**: a CI *for patient A* is
    conditional on patient A, so patient is fixed, not random, and resampling patients
    answers a different question. The correct split:

    | quantity | resampling scheme |
    |---|---|
    | **per-patient CI** on `frac_double_negative` | **sample → cell**, within that patient |
    | **cohort-level** inference (mean escape, regression coefficients, distributions) | **patient → sample → cell** |

    Report the flat and sample-aware per-patient intervals side by side so the
    narrowing is visible. **For the many patients with a single sample this reduces to
    a cell bootstrap, which cannot see sample-level technical or biological variation
    at all** — their CIs are therefore optimistic in a way multi-sample patients' are
    not, and that asymmetry is stated with the ranking rather than buried. Correct
    nesting is S1-gated (it needs the real sample→patient map), so this is provisional
    until S1 lands, like every other per-patient aggregate here.
- **Multi-antigen combinatorial coverage matrix (new).** `SLAMF7`/`FCRL5` are
  promoted from "backups" to a deliverable. For every pair and triple over
  {`TNFRSF17`, `GPRC5D`, `SLAMF7`, `FCRL5`, `CD38`, `SDC1`, `ITGB7`}, compute the
  uncovered fraction of each patient's clone. This answers the question a
  target-strategy audience actually asks — *is BCMA+GPRC5D the best pair for this
  patient, or would BCMA+FCRL5 cover more?* — which the two-antigen metric alone
  cannot.

  **Report it as separate columns; do NOT collapse it into a utility score
  (2026-08-21).** The same discipline that keeps DN fraction, co-escape and coherence
  apart at stage 12 applies here. A weighted `coverage − λ · exposure` needs a
  principled λ, and there isn't one — the weights would encode a clinical judgement
  the data cannot supply, while hiding the inputs that a reader could otherwise
  disagree with. Per pair/triple, per patient:

  | column | source |
  |---|---|
  | uncovered fraction | this stage |
  | incremental gain vs. the best single target | this stage, `P(A⁻) − P(A⁻ ∩ B⁻)` |
  | co-loss enrichment | this stage, depth-conditioned |
  | **normal *marrow* expression** | stage 09 |

  The last column is normal **marrow** expression specifically — not "normal tissue",
  which this dataset cannot observe. Coverage is read against it rather than maximized
  blindly: a target covering 100% of the tumor that also hits normal marrow plasma
  cells is not a better target. Extra-marrow liabilities (GPRC5D in keratinized tissue)
  stay a cited external caveat, never a measured column.
- **The bias table** (in the QC methodology section above) is authored as a figure
  here and referenced from stage 12.

**BCMA/GPRC5D co-negativity enrichment — the key derived metric (added 2026-08-21).**
`frac_double_negative` alone cannot distinguish two clinically different tumors. Per
patient, build the 2×2 contingency table over malignant cells:

|  | GPRC5D⁺ | GPRC5D⁻ |
|---|---|---|
| **BCMA⁺** | dual-positive | BCMA-only |
| **BCMA⁻** | GPRC5D-only | **double-negative** |

and compare the observed DN fraction against the independence expectation
`E[DN] = P(BCMA⁻) × P(GPRC5D⁻)`. Report the **co-escape enrichment ratio**
`observed(DN) / expected(DN)` with a Fisher's exact test and a permutation CI. This
separates three facts the single metric fuses into one: how often each antigen is
individually absent, how many cells are DN, and whether the *same* cells are
disproportionately losing both. A patient at 6% DN ≈ 0.3 × 0.2 has two independent
partial failures; a patient at 6% DN against a 1.5% independence expectation has a
coordinated antigen-low phenotype, and is the one stage 10 then investigates
mechanistically.

**What co-escape enrichment does NOT mean (corrected 2026-08-21).** An earlier draft
said an enriched patient "is the one dual targeting doesn't help", and that co-escape
"determines whether a second binder helps at all". **Both overstate it, and the
arithmetic shows why.** Adding GPRC5D to BCMA moves the uncovered fraction from
`P(BCMA⁻)` to `P(BCMA⁻ ∩ GPRC5D⁻)`. At 30% BCMA⁻ / 20% GPRC5D⁻ under independence that
is 30% → 6%. With co-loss enrichment pushing DN to 15%, it is 30% → 15% — less than
independence promised, but still halving the escape population. Enrichment measures
**how much of the two targets' expected complementarity is eroded by correlated loss**,
not whether the second target is worth adding. Use that framing everywhere: it is both
more precise and more useful to a target-strategy audience than the binary claim it
replaces.

**Incremental coverage gain — reported alongside (added 2026-08-21).** Co-escape
enrichment is a statement about *biology* (is loss correlated); the clinical question is
a statement about *value* (what does the second target buy). Different quantities, both
cheap off the same 2×2:

    gain from adding GPRC5D to BCMA  =  P(BCMA⁻)   − P(BCMA⁻ ∩ GPRC5D⁻)
    gain from adding BCMA to GPRC5D  =  P(GPRC5D⁻) − P(BCMA⁻ ∩ GPRC5D⁻)

Report both per patient with CIs, as separate columns. A patient can carry high
enrichment *and* a large incremental gain — those are not in tension, and collapsing
them into one number would hide exactly that case. This is the quantity a single- vs.
dual- vs. sequential-target discussion actually turns on, and it generalizes directly
to the coverage matrix below.

**The null must be depth-conditioned, or this test measures library size.** Dropout is
a per-*cell* property: a shallow cell is more likely to read zero for *both* genes, so
depth heterogeneity alone produces positive BCMA⁻/GPRC5D⁻ association. A permutation
that shuffles antigen labels freely within a patient destroys the depth↔label coupling
and will report co-escape enrichment on data with no biological co-occurrence at all —
an artifact pointing in exactly the direction the project wants to find, which is the
worst kind. Therefore:
- Stratify cells by sequencing depth (or `n_genes_by_counts`) within patient and
  **permute labels within stratum**; equivalently, compute `E[DN]` from a per-cell
  independence model where `P_i(BCMA⁻)` and `P_i(GPRC5D⁻)` are functions of cell *i*'s
  own depth, and sum over cells.
- **Report the unconditioned ratio next to the conditioned one.** The gap between them
  *is* the depth artifact, quantified — which is a more convincing exhibit than the
  conditioned number on its own.
- This supersedes the "label-permutation null" originally filed under stage 09: that
  test, as written, was already this test — it just wasn't labelled as one.

**The detection curve, and what it can and cannot deliver (corrected 2026-08-21).**
An earlier draft of this stage proposed a "dropout-adjusted expected DN" computed as
`Σ_i P_i(BCMA⁻) · P_i(GPRC5D⁻)`. **That formula is circular and is not used as a
corrected estimate**: multiplying the two marginals assumes exactly the independence
that the co-escape test above exists to interrogate, so a tumor with genuinely
correlated antigen loss would be "corrected" toward the very null it violates.

Two consequences, and the first is a simplification worth having:

- **The "dropout-adjusted DN" and the "expected DN under depth-conditioned
  independence" are the same number.** They were specified as two separate
  deliverables; they are one computation serving one purpose. It is reported once, as
  the **depth-adjusted DN expectation under conditional independence** — a *technical
  baseline* the observed value is compared against, never a corrected truth. Merging
  them also removes a deliverable that would have invited exactly the misreading above.
- **No dropout-corrected DN point estimate is produced, and none is claimed.** Dropout
  is *bounded* here — by the threshold sensitivity band, the expression-matched
  false-negative floor, the depth regression and the downsampling check — not corrected.
  The observed DN stays the point estimate, reported as an interval. Saying so plainly
  is stronger than shipping a number whose correction rests on an assumption the
  project is simultaneously testing.

Still build the detection curve: fit detection probability against cell depth and gene
mean on the expression-matched control genes already selected for the false-negative
floor, giving each observed zero an approximate `P(false zero)`. It is what makes the
depth-conditioned null above quantitative rather than rank-based, and it is what turns
"GPRC5D is lowly expressed" into a number. It just does not license a corrected DN.

**A genuinely dropout-corrected DN would need a joint model, and is deferred.** The
defensible version is a latent-class model over the four true states
(B⁺G⁺ / B⁺G⁻ / B⁻G⁺ / B⁻G⁻) with per-cell detection probabilities from the curve above,
fit by EM over the observed 2×2 — which estimates the true joint *without* assuming
independence, and yields the co-escape enrichment as a by-product rather than an input.
That is a real piece of statistical work and it is **not** on the critical path: it is
filed here so it is not reinvented casually, and so that the current position ("bounded,
not corrected") is understood as a deliberate choice rather than an oversight.

**Imputation/denoising (MAGIC, scVI, ALRA, …) is forbidden for positivity calls**, and
this is not a stylistic preference: imputation manufactures low-level expression by
borrowing from neighbors, and the entire scientific question here is whether a
transcript is genuinely absent. Smoothing over the zeros erases the measurement. The
detection-curve approach models the uncertainty instead of filling it in.

**09 — Escape robustness** (`notebooks/09_escape_robustness.ipynb`,
`src/mm_escape/bulk.py` + `robustness.py`; env: `mm-core`). New stage. Everything
here exists to answer "how do you know your escape fractions are real?"
- **Matched bulk RNA-seq — orthogonal validation of antigen *abundance*, not of the
  DN fraction (scope corrected 2026-08-21).** For the ~28 samples with matched bulk,
  correlate malignant-cell pseudobulk `TNFRSF17`/`GPRC5D` against bulk TPM (Spearman,
  concordance, residuals, and the named discordant cases). The load-bearing question
  is whether scRNA zero-rates run systematically high where bulk says the transcript
  is plainly present — which is direct, quantified evidence of dropout and feeds back
  into stage 08's false-negative floor.
  **What bulk cannot do is validate `frac_double_negative`**, and the earlier wording
  claimed it could. Bulk destroys the joint single-cell distribution: a tumor that is
  50% BCMA⁺GPRC5D⁻ plus 50% BCMA⁻GPRC5D⁺ shows healthy bulk expression of *both*
  genes while containing zero dual-positive cells — and the converse misreads are
  equally available. Bulk constrains marginal abundance per gene; the *joint*
  distribution over cells is visible only in single-cell data and has no orthogonal
  check in this project. Phrase every output here as **"orthogonal validation of
  antigen abundance and the plausibility of scRNA-derived antigen-negative calls"**,
  never as validation of the escape fraction.
  Handle the two empty 114-byte stubs and the three ID mismatches documented in the
  Data section.
- **Normal plasma-cell antigen baseline — marrow on-target/off-tumor expression
  context (scope corrected 2026-08-21).** Do *normal* plasma cells (from `BM*`/`ND_*`
  marrow) express BCMA and GPRC5D, and what do other marrow lineages show? This is
  real and worth having: BCMA carries broad normal-PC and B-lineage expression, and
  the malignant-vs-normal-PC contrast is what makes a coverage number interpretable
  rather than absolute. It feeds the coverage matrix's risk trade-off in stage 08.
  **It is not a safety axis, and must not be called one.** GPRC5D's clinically
  decisive off-tumor liability is keratinized tissue — the nail, skin and taste
  toxicity seen with talquetamab — and a bone marrow dataset cannot observe that at
  all. Expression is also not toxicity. Keep three things separate in the writeup:
  **(a) tumor coverage**, **(b) normal *marrow* expression** (measured here),
  **(c) known extra-marrow liabilities** (external evidence, cited, not measured).
  A genuine target-ranking utility score of the form
  `coverage − λ · normal-tissue exposure` needs a normal-tissue atlas (GTEx/HPA or a
  normal scRNA atlas) and is filed as a future extension, not claimed from this data.
- **Label-permutation null → moved to stage 08 as the co-negativity test.** Permuting
  antigen labels within patient while preserving marginals does not establish "what
  the metric looks like under no signal" — it holds each antigen's negative rate fixed
  and tests only whether the two negativities *co-occur* beyond independence. That is
  a better question than the one it was written for, so it moves to stage 08 as the
  co-escape enrichment test, with the depth-stratified null documented there.

**10 — Escape subclone + phenotype** (`notebooks/10_escape_subclone_phenotype.ipynb`,
`src/mm_escape/subclone.py`; env: `mm-core`). New stage, and **the project's actual
scientific payoff** rather than another robustness check.
- **Is the double-negative population structured, or scattered noise?** "3% of this
  patient's cells are double-negative" and "this patient has a pre-existing 3%
  resistant subclone" are different claims, and **only the second one predicts
  selection under therapy** — which is the entire clinical premise of the project.

  **Transcriptional clustering alone does NOT establish clonality (corrected
  2026-08-21).** The earlier wording — "random scatter is the signature of dropout;
  spatial clustering is the signature of a real subclone" — was too strong in both
  directions. A transcriptionally coherent group can arise from cell cycle, stress,
  interferon tone, metabolic state, sequencing depth, or sample-prep batch as easily
  as from a genetic subclone; and conversely, cells of a genuine genetic clone need
  not form a tidy transcriptional island. The name `clonality-of-escape` prejudged
  exactly the question the analysis is supposed to ask, so it is retired in favour of
  **DN coherence**, evaluated at three escalating levels with the claim escalating
  with it:

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
  normal, and at this cohort's ~2,044 median genes/cell it will often be
  underpowered. Report level 3 as **supported / not evaluable**, with the per-patient
  CNV resolution stated — never as "no CNV subclone". Treating an underpowered null
  as a negative result would systematically understate exactly the risk this project
  exists to measure.

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
  (**`NSD2` is `WHSC1` in the older reference — depends on stage 05's symbol
  harmonization; without it this class cannot be called at all**)
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

  **It is a transcriptional proxy for the translocation, not a detection of it
  (terminology fixed 2026-08-21).** `NSD2`/`FGFR3` overexpression is *consistent with*
  t(4;14); it is not a breakpoint call, and expression can be driven by other things.
  Every output label reads **"TC-like expression subtype"** or **"transcriptionally
  inferred TC class"** — never "patient has t(4;14)". This costs nothing and removes
  the single easiest claim in the project to attack. When S1 lands with real
  cytogenetics, the proxy becomes testable rather than assumed, and the reportable
  result is *"the expression-based TC proxy agreed with clinical cytogenetics in X/Y
  evaluable patients"* — which is a better finding than the proxy alone ever was.

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
- Escape-fraction table, annotated with disease stage/cytogenetics from
  Supplementary Table S1 where available.
- **Caterpillar plot with confidence intervals**, replacing the original ranked bar
  chart — a bar chart implies a precision this metric does not have.
- **Risk tiers, not a rank ordering (decided 2026-08-21).** The caterpillar plot fixed
  the chart but not the deliverable: printing "#1 patient 123, #2 patient 456" is
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

**Stages 01-03 are run and verified. Stage 04's loader exists; its QC does not yet.**
The working tree is clean (R build removed, preserved under `r-build-snapshot`), `raw/`
is intact at 62 samples, `scripts/01-03` are confirmed (62/62 `triplet-ok`), notebooks
01-03 are written and executed, the four `envs/*.yml` are built with kernels
registered, and `src/mm_escape/` holds `config.py`, `gene_space.py` and `io.py`. See
`RESUME_HERE.md` for exact session state as work proceeds.

**`io.py` is written and validated against the real files (2026-08-24).** All 61
retained samples load in ~4 s; the three-sample failure-mode set (`MMRF_1695` = 33538
build, `27522_1` = 33694 build, `BM4` = normal-BM control) round-trips through
`gene_space.py` to 32,991 genes with all 65 required genes present. Verified: the
genes x cells -> cells x genes transpose (against random raw `.mtx` triplets plus exact
total-count and nnz equality), deposited gene order preserved untouched for
`attach_ensembl_ids`'s positional join, `<sample_name>_<barcode>` obs_names unique
across the concat, and five failure paths raising. Cohort **pre-QC** totals: 202,203
cells over 61 samples (154,053 disease / 48,150 normal-BM), 509 / 2,767 / 9,328
min/median/max cells per sample (18.3x spread).

First actions, in order:
1. Write `src/mm_escape/qc.py` (MAD outlier calling + the `scDblFinder` rpy2 bridge)
   and `notebooks/04_qc.ipynb` against it, in `mm-qc`. The loader it builds on is
   done — do not re-litigate `io.py`.
2. Re-derive the MAD thresholds and the `pct_counts_mt` cap against THIS cohort's
   distributions; do not copy `sc-best-practices`'s PBMC/BMMC demo numbers.
3. Checkpoint each sample's post-QC AnnData individually (resumable per sample), then
   stage 05 does the gene-space intersection on those checkpoints.

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
   **First presentable state: escape fractions with co-escape enrichment.**
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
- **The gene-space join is on Ensembl ID, reconstructed and verified (2026-08-21).**
  The deposit is symbol-only, but the three distinct `genes.tsv` files are positional
  dumps of public references and were reproduced exactly (0 mismatches / 33538 and
  33694 rows) from the Ensembl 93 / 84 GTFs plus Seurat's `gsub("_","-")` +
  `make.unique`. Mapping committed at `resources/gene_space/`. Recovers **32,991 genes
  vs. 22,164 on symbols (+10,827)**. Do not re-open, do not regenerate, and do not
  search for an `ENSG` column in the raw files — there isn't one.
- **The four-gene alias dictionary is now only a regression assertion.** It addressed
  4 of the 11,140 symbols that drift between builds; it was never a harmonization
  method and must not be treated as one.
- **The "drop ~52 `make.unique`-ambiguous symbols" interim is superseded** — those
  genes resolve correctly under the ID join. The `.N` suffixes encode row order, which
  is why a symbol join could pair the *wrong* gene (`TBCE` is a different Ensembl entry
  in each build); IDs remove the ambiguity rather than working around it.
- **`var_names` = Ensembl ID through the merge, canonical symbol afterwards.** Identity
  matters at the join; readability matters everywhere downstream, and the whole scanpy
  stack is symbol-native. The 9 symbols still colliding after the ID intersection get
  `SYMBOL__ENSG...`, never `var_names_make_unique()`.
- **AnnData, not MuData.** One cell-level modality; the stage-09 bulk is sample-level
  and joins on `sample_id` as a DataFrame. MuData would add no capability here.
- **Gene symbols are canonicalized before the gene-space intersection.** The two
  Cell Ranger references use different HGNC vintages; intersecting raw symbols drops
  `NSD2`/`WHSC1`, `TENT5C`/`FAM46C`, `NSD3`/`WHSC1L1`, `ATP5F1A`/`ATP5A1`. A missing
  required gene means "check for a legacy symbol", not "biologically absent".
- **The required-gene assertions stay, and stay loud.** They caught the `NSD2` drift
  that manual inspection had missed across two prior builds of this project.
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
  It is labelled a **TC-like expression subtype**, never a translocation call — it is
  a transcriptional proxy, and S1's cytogenetics is what would test it.
- **R stays isolated in its own environments** (`env-qc` for `scDblFinder`,
  `env-annotation` for `SingleR`) — never merged into `mm-core`.
- **Malignant calling via light-chain restriction, not clustering alone** — and by
  **ratio**, not presence/absence, because IG genes are the most ambient-contaminated
  in this tissue.
- **`infercnvpy` is required, not optional**, with the agreement rate reported. This
  stage sets the metric's denominator; its errors propagate straight into the
  headline number.
- **Normal-BM samples are controls, not filler.** They validate the malignant caller
  (polyclonal marrow must yield no clone) and provide the normal-PC antigen baseline.
  They are not to be dropped as "not myeloma." But that baseline is **marrow
  expression context, not a safety axis** — GPRC5D's decisive off-tumor site is
  keratinized tissue, which this dataset cannot see. Tumor coverage, normal *marrow*
  expression, and known extra-marrow liabilities stay three separate things.
- **Dropout is bounded, not just mentioned.** It is the opposite-signed and larger
  counterpart to ambient RNA, and `GPRC5D` is a low-abundance transcript. The
  headline metric is reported as a bracketed interval with a threshold sensitivity
  band, never as a bare point estimate.
- **Ranking stability across thresholds is the robustness *diagnostic*; risk tiers
  are the deliverable.** No single threshold's value is the claim, and neither is an
  ordinal rank when the CIs overlap — stage 12 emits robust-high / uncertain /
  robust-low, with co-escape enrichment and DN coherence as separate columns.
- **Matched bulk RNA-seq (GSE223061) is used, not shelved.** It was already
  downloaded and previously unused; it is the only orthogonal check available on the
  antigen quantification. It validates **antigen abundance and the plausibility of
  antigen-negative calls** — never the dual-negative fraction itself, which is a
  joint single-cell quantity that bulk destroys by construction.
- **No composite risk score, anywhere.** DN fraction, incremental coverage gain,
  co-escape enrichment, DN coherence, malignant confidence, threshold sensitivity and
  bulk concordance stay separate columns — in the patient table and in the coverage
  matrix alike. A weighted utility would need principled weights that do not exist, and
  would hide the inputs a reader could otherwise disagree with.
- **Co-escape enrichment measures eroded complementarity, not futility.** Adding a
  second target moves the uncovered fraction from `P(A⁻)` to `P(A⁻ ∩ B⁻)`, which is a
  real gain even under strong correlated loss. Never write that enrichment means dual
  targeting "doesn't help"; report the incremental gain next to it.
- **No dropout-corrected DN point estimate is claimed.** `Σ P(A⁻)·P(B⁻)` is the
  independence baseline the co-escape test compares against, *not* a correction —
  using it as one would assume away the dependence being measured. Dropout is bounded
  (sensitivity band, false-negative floor, depth regression, downsampling), not
  corrected. A joint latent-class/EM model over the four true states is the honest way
  to correct it and is deliberately deferred, not forgotten.
- **Bootstrap at the level of the question**: sample → cell within patient for a
  per-patient CI (patient is conditioned on, not random); patient → sample → cell for
  cohort-level quantities. Single-sample patients get optimistic CIs and this is stated.
- **Co-negativity enrichment is a first-class result, and its null is
  depth-stratified.** Whether the *same* cells lose both antigens is a different and
  sharper question than how many are double-negative. An unconditioned permutation
  null would report enrichment from library-size variation alone, biased toward the
  project's own hypothesis — so the null is permuted within depth strata and the
  unconditioned value is reported next to it as the size of the artifact.
- **"Subclone" requires CNV support (stage 10).** Transcriptional clustering of DN
  cells establishes an escape-associated *state*, not a clone — coherence has many
  non-genetic causes. `clonality-of-escape` is retired for a three-level coherence
  hierarchy, and level 3 reports *supported / not evaluable*, never "no CNV
  subclone", because within-clone CNV resolution is often underpowered here.
- **Binary antigen calls stay primary; a dropout-adjusted estimate runs alongside;
  imputation is forbidden.** Denoising manufactures low-level expression, and whether
  a transcript is genuinely absent is the entire question. Model the uncertainty with
  a detection curve, never fill it in.
- **Annotation numbers are concordance, not validation accuracy** (stage 06). Manual
  labels are a third opinion, not ground truth, and the two automated methods share
  marker-biology priors. The marker-coverage test is the biological evidence and can
  veto a class regardless of concordance.
- **Stage 11 is exploratory and ranks ninth**, not co-equal with the antigen
  analysis — n ≈ 41 against hundreds of LR pairs with a correlated confounder does
  not support a confirmatory claim.
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
- **The whole analysis runs in notebooks, stages 01-12** — every stage is openable and
  steppable. `scripts/01-03` are kept as a headless CLI fallback and are *wrapped* by
  notebooks 01-03, never reimplemented; their output is the contract and byte-identical
  parity is verified.
- **Notebooks carry the analysis; `src/mm_escape/` carries what is reusable, testable,
  or fiddly.** Paired via `jupytext`, `.ipynb` gitignored. This is a deliberate
  relaxation of the earlier "notebooks are thin orchestration" rule — the goal was never
  thinness, it was avoiding duplicated logic and keeping review on `.py` diffs.
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
- **The minimum malignant-cell inclusion threshold** for stage 08 — ≥50 cells is a
  floor, not the answer (at n=50, one cell is 2%, so 1%/2%/3% escape are not
  separable). Re-derive it from the smallest DN fraction the project intends to call
  meaningful once the per-patient malignant-cell distribution is known; expect
  100-200. Fixed before the ranking is looked at.
- **How many patients survive that threshold.** Cohort cell counts vary ~15×
  (min 480, median 2,555, max 7,937 cells/sample post-QC); if a large share of
  patients fall below the malignant-cell minimum, the ranking's usable n may be
  well under 41 and the framing must adjust honestly rather than quietly.
- **Whether a published CITE-seq/flow calibration exists** for BCMA and GPRC5D
  mRNA-vs-surface-protein correlation in myeloma — determines whether stage 12's
  protein limitation can be quantified or only stated.
- **Whether `infercnvpy` resolves sub-clonal CNV structure within a single patient's
  clone at this depth** — gates level 3 of the stage-10 coherence hierarchy, and
  therefore whether the project can use the word "subclone" at all. Determine per
  patient and report resolution; do not assume either outcome.
- **How much of any observed co-escape enrichment survives depth stratification** —
  the gap between the unconditioned and depth-conditioned ratios is itself a
  reportable number, and if the conditioned enrichment collapses to ~1 across the
  cohort, that is a real (negative) finding about dual-antigen escape, not a failed
  analysis.
- **Whether the stage-08 ranking survives restriction to `high`-confidence malignant
  cells** (stage 07 tiers). If it does, the metric is robust to the weakest link in
  its own denominator; if it does not, the CNV-inconclusive cells are driving it and
  the framing must say so.
