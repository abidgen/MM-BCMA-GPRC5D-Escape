# Multiple Myeloma Dual-Antigen (BCMA/GPRC5D) Coverage & Escape Risk

> **Naming note (updated for Stage 12's completion).** Several stage headings in this
> document once used planning-era notebook names. They now match the files on disk:
> `07_malignant_plasma`, `08_dual_antigen_escape`, `09_bulk_validation`,
> `10_dn_coherence`, `11_immune_context`, plus the lettered arms `05b`, `08c`,
> `09b`, `11b`, and finally `12_final_synthesis` — renamed from the planning-era
> "decision packet" for the same reason Stage 11 was renamed from `11_cellchat_liana`:
> the executed stage turned out narrower and more disciplined than the original name
> implied. **Phase 1 (stages 01-12) is complete as of 2026-08-27.**
>
> The **producers** of the frozen Stage 06-10 tables are committed under
> `production/`; the main project document remains current where this one disagrees.
## Pipeline Walkthrough (Python rebuild)

**Objective:** For each multiple myeloma patient in GSE223060, quantify the fraction
of malignant plasma cells that would evade BOTH a BCMA-directed and a
GPRC5D-directed CAR-T/TCE, then relate that "dual-antigen escape fraction" to the
immune microenvironment via LIANA+ (CellChat's algorithm, natively reimplemented in
Python), and lay out the resulting evidence for a single- vs. dual- vs.
sequential-target CAR-T strategy discussion. **The output is not a per-patient
risk-tier ranking** — that was the original plan, and Stage 12 explicitly declined it
once the frozen evidence was assembled, because the patients who looked strongest on
measurement robustness and the patients who looked strongest on structural evidence
turned out to be disjoint sets. The actual output is a **six-axis per-patient evidence
matrix** (measurement, structure, phenotype, genomic, immune, coverage) kept as separate
columns, plus a claim ladder stating exactly what is and is not supported. See Stage 12
below for the full reasoning.

**This document is a narrative walkthrough of the pipeline's logic and reasoning,
not a specification.** It answers *"what does the pipeline do, and why is it built
that way"*. It deliberately does not restate thresholds, gene lists, resampling
schemes or pre-declared decision rules — those live in the main project document, which
is authoritative, and duplicating them here is exactly what made this file go stale once
already during the R build. For exact implementation read `src/mm_escape/`; for
superseded positions and how the plan evolved read `docs/decisions-archive.md`.

**Where the numbers come from.** Run output for stages 01-05b is described in
`docs/stage-results.md` and stored in the `results/*.csv` files it points at. Frozen
scientific outputs are local or externally archived artifacts authenticated by the
committed pre-Stage-12 provenance manifest; an ignored/local CSV alone is not durable
provenance and is authoritative only when its hash matches that manifest. Committed
production code and provenance metadata document how those artifacts were generated and
verified. Dataset facts (file-naming quirks, the two Cell Ranger references, the patient
mapping) are in
`docs/dataset-ground-truth.md`. This is a from-scratch Python rebuild: no R code carries
over, all data-format knowledge does, and the R build is preserved in git under the
`r-build-snapshot` tag.

**Numbering:** stage numbers match notebook filenames and their `results/NN_*/` output
directory one-to-one. **Number order is execution order** — 04 → 05 → 06 → 07 → 08 →
09 → 10 → 11 → 12, with no exceptions, and every stage from 01 to 12 is a notebook you
can open and step through.

**The through-line.** The headline metric is a *fraction of zeros* — the noisiest
thing scRNA-seq measures — and it has two error sources pointing in opposite
directions. Most of the pipeline past Stage 08 exists to bound both of them rather
than to add analyses.

---

## Stages 01-03 — Data acquisition (`notebooks/01_download_data.ipynb`, `02_check_files.ipynb`, `03_build_manifest.ipynb`; env: `mm-qc`)

Downloads and unpacks `GSE223060_RAW.tar` / `GSE223061_RAW.tar` from GEO's FTP, verifies
per-sample file structure, and builds `raw/sample_manifest.csv` mapping each sample ID to
its exact `barcodes.tsv`/`genes.tsv`/`counts.mtx` paths.

These three run as notebooks like every other stage, but they **wrap** the original
bash/Python scripts rather than reimplementing them. The scripts were already solved and
verified against the real archive, they run headlessly (which a notebook does not — useful
for a fresh clone or a remote box), and having one implementation means the two paths
cannot disagree. Byte-identical output is verified in both directions, so `scripts/01-03`
stay in the repo as a CLI fallback.

**Status: run and confirmed.** All 62 samples `triplet-ok`, manifest reproduced exactly.

Notebook 03 does more than the script it wraps, which is the point of having it: it
reports the Cell Ranger reference split, previews the gene intersection, and runs the
required-gene assertions. Those assertions are what surfaced the `NSD2` symbol drift
described under Stage 05 — a defect two earlier builds of this project walked straight
past. Notebook 02 contributes the observation the whole gene-space solution rests on:
across 62 samples there are only **three distinct `genes.tsv` checksums**, which is why
the deposited symbol columns could be reconstructed from public references at all.

---

## Stage 04 — QC + doublet detection (`notebooks/04_qc.ipynb`, `src/mm_escape/qc.py`; env: `mm-qc`)

Loads each sample via a custom loader (`src/mm_escape/io.py`) rather than
`scanpy.read_10x_mtx()`, since this archive's filenames (`counts.mtx`, single-column
`genes.tsv`) don't match what that convenience function expects.

QC uses MAD-based outlier filtering rather than fixed thresholds, following
`sc-best-practices.org`'s QC chapter — but with the tutorial's numeric cutoffs re-derived
against this cohort, since they were set for healthy PBMC/BMMC and this is myeloma
marrow. **Run 2026-08-24: 204,040 → 172,940 cells (84.8% kept).** Two things the data
changed:

- **`pct_counts_in_top_20_genes` is computed and reported but does not filter.** A 5-MAD
  band on it flags 17% of MMRF and 15% of WashU 1 against 3% of WashU 2, and the flagged
  cells are two quite different populations: `IGKC` at ~25% of counts (plasma cells) and
  haemoglobin at ~32% (erythroid debris). The plasma-cell half is the project's subject —
  `TNFRSF17` is detected in **21.8%** of that decile against **0.8%** elsewhere. An
  Ig-dominated library is a plasma cell's *normal* state, because a plasma cell is a
  professional secretor; filtering on it would preferentially delete
  antigen-**positive** malignant cells and inflate the escape fraction. It is kept as an
  ambient-Ig handle for Stage 08, which needs one, but is not allowed to remove cells.
- **Thresholds are derived per cohort, not pooled.** The cohorts ran different 10x
  chemistries and differ ~1.9× in genes detected per cell; a single pooled MAD across
  that spread would flag much of WashU cohort 1 as low-quality for what is a batch
  difference, not a cell-quality one.

**And the deposit turns out to be pre-filtered, differently in each cohort:** WashU 1
and 2 at <10,000 UMIs and <20% mt, MMRF uncensored on UMIs but cut at <10% mt, donors
uncensored on both. That WashU UMI ceiling censors a band enriched **20-70× for
`GPRC5D`** and is a live confounder for the headline metric — it recurs at Stage 05 and
is owed a sensitivity analysis at Stage 08.

Doublet detection uses `scDblFinder` (R), called from Python via an isolated `rpy2`
bridge contained entirely within the `mm-qc` environment — kept as R because it's the
benchmarked best-performing method (Xi & Li 2021), not swapped to a pure-Python
alternative purely for language purity.

`56203_1` is **repaired on load, not excluded**. It was long believed to come from an
incompatible 22184-gene reference missing BCMA; it is in fact a normal 33694-build
sample whose `genes.tsv` write stopped mid-symbol at row 22185, putting `TNFRSF17` and
`IGLC1/2/3` past the cut rather than absent from a reference.

**Ambient RNA correction (SoupX/DecontX) is not run — not an oversight, a hard
constraint.** Both need the unfiltered Cell Ranger matrix, including empty droplets, to
estimate the background contamination profile, and only the filtered per-sample matrices
were ever deposited. (Raw reads do exist under controlled access on dbGaP, but no
unfiltered matrices exist anywhere, so the constraint holds either way.) The mitigation
lives at Stage 08 instead — an empirically derived noise floor rather than formal
correction. This matters more for plasma-cell data than most: plasma cells secrete
enormous quantities of immunoglobulin transcript, making ambient contamination of
`IGKC`/`IGLC` — and potentially of the antigen genes themselves — a real and
non-hypothetical concern, which is also why Stage 07's light-chain call is ratio-based
rather than presence-based.

**Ambient RNA is only half the problem, though.** It biases in one direction: true
negatives read as faintly positive, which *deflates* the escape fraction. Dropout biases
the other way — a cell that genuinely expresses the antigen reads as zero because the
transcript was never captured, which *inflates* it. For this project dropout is the
larger of the two, because `GPRC5D` is a low-abundance GPCR transcript and this cohort
is shallow: **1,162 median genes per cell** across the 172,940 that survive QC. Neither
bias may be left unquantified; Stages 08 and 09 do the bounding.

Each sample's post-QC AnnData is checkpointed individually, and **QC annotates rather
than deletes** — every barcode survives with `obs["keep"]` set, so "would this hold
under a different QC?" stays an answerable question rather than a lost one. (The load
itself is no longer the slow part: the Python loader reads all 62 samples, 204,040
pre-QC cells, in about two seconds.)

---

## Stage 05 — Gene-space intersection, integration, clustering (`notebooks/05_integration_clustering.ipynb`; env: `mm-core`)

Before any concatenation: intersect gene sets across all retained samples — the 62
samples split across 33538- and 33694-gene Cell Ranger references, and a union merge
would make ~11k genes structurally zero in whole cohorts, indistinguishable downstream
from the true biological zero this project measures. Assert all required marker and
antigen genes survive, hard-failing with the specific missing gene names otherwise.

**The reason that assertion exists is worth the space.** The two references were built
against different HGNC symbol vintages, so some genes exist in **both** builds under
different names, and a naive symbol intersection throws them away silently. This was
found the hard way when Stage 03's assertions failed on `NSD2` — which the older build
calls `WHSC1`. The others are `TENT5C`/`FAM46C` (a recurrently deleted myeloma tumour
suppressor), `NSD3`/`WHSC1L1` (a genuinely different gene from NSD2 — do not conflate
them), and `ATP5F1A`/`ATP5A1`.

The `NSD2` case is the one that matters most. NSD2 is how t(4;14) is read
transcriptionally, and t(4;14) is the highest-risk myeloma translocation. Left
unharmonized, Stage 10's molecular-subgroup call would quietly lose its highest-risk
class — not error out, just never report it.

The wider lesson those assertions bought: **a required gene coming up missing means
"check for a legacy symbol", not "biologically absent".** And the initial guess that the
~11k symbols unique to each build were mostly annotation-version noise (`AC000032.1`
versus `AC000032.2`) was far too optimistic. Joining on reconstructed Ensembl IDs rather
than symbols recovers **32,991 genes against 22,164**, because 11,140 intersected IDs
simply carry a different symbol in each build. The symbol join was also *unsafe*, not
merely lossy: `TBCE` is a different Ensembl entry in each build, so a symbol join
silently merges two unrelated rows.

Then normalize, select highly variable genes, PCA, and `harmonypy` integration keyed on
`patient_id` with `n_genes_ref` and `cohort` as additional covariates; Leiden
clustering; UMAP. **The two covariates are not interchangeable, which is easy to get
wrong.** Reference build and collection cohort cut across each other — two WashU cohort
1 samples sit on the newer 33538 build while the four `ND_*` donors sit on the older
33694 one. `n_genes_ref` captures which reference processed a sample; `cohort` captures
which chemistry and protocol generated it, and that is where the ~1.9× depth difference
lives. Neither substitutes for the other.

**The result splits by compartment, and the split is informative.** Harmony mixes the
immune compartment well (cohort-mixing entropy 0.751) and the plasma compartment barely
at all (0.105) — with the three largest plasma clusters coming out one per cohort, each
spanning ~30 patients, which rules out the benign explanation: a patient-private clone
would fragment into ~41 clusters, not three cohort-shaped ones.

**This is the Stage 04 censoring surfacing a second time.** Median UMIs in the plasma
compartment run MMRF 22,477 against WashU 1's 5,036 — cells the WashU deposits *cannot
contain*, because WashU was cut at 10,000 before deposit. So this is not Harmony failing:
what separates the compartments is a non-recoverable sampling asymmetry, and no
correction restores cells that were never deposited. The compartment specificity follows
directly — T/NK/myeloid/B cells sit well below 10,000 UMIs everywhere, so the ceiling
never touched them. (Numbers: `docs/stage-results.md`.)

**One incidental finding worth recording: `GPRC5D` is not a highly variable gene** —
mean 0.061 against `TNFRSF17`'s 0.492, an 8× gap, and HVG in only 6 of ~50 patients.
Mechanically this changes nothing, since the embedding does not need it and Stage 08
reads raw counts. Its value is *evidential*: the claim that dropout matters more for
GPRC5D than for BCMA had been an inference from the literature, and this is the first
number from this cohort supporting it.

**Integration is deliberately confined to the immune compartment.** Harmonizing on
`patient_id` is right for T/NK/myeloid cells, which should look alike across patients.
It is actively risky for the tumor: the malignant clone is patient-private by
definition, so forcing patients together can blend distinct clones and erase the
heterogeneity the project exists to measure. So the integrated embedding is used for
immune annotation and clustering only; **all malignant subclustering happens per
patient, un-integrated** (Stage 10 relies on this). The reassuring part, worth stating
out loud because it is the obvious objection: **per-cell antigen calls are raw counts
and never touch the embedding**, so no integration choice can distort the escape
fractions themselves.

### Stage 05b — the integration method is benchmarked, not assumed

Harmony was chosen above because it was the obvious default.
`notebooks/05b_integration_benchmark.ipynb` (env `mm-integration`) tests that, running
seven arms — unintegrated, Harmony / scVI / Scanorama on `sample_name`, the stage-05
configuration, and Harmony / scVI on `cohort` — scored with `scib-metrics` against
provisional CellTypist labels. It reads stage 05's output **read-only** and asserts the
file is byte-identical afterwards.

**Scoring is on the immune compartment only, and that is the design's whole point.**
scIB's batch metrics cannot distinguish *"correctly left apart"* from *"failed to
merge"*. Because WashU was cut at 10,000 UMIs before deposit, an aggressive method can
fuse the plasma populations and score better while manufacturing correspondence that is
not recoverable. So plasma mixing is **reported as a diagnostic and never optimized**.

**Result: no arm qualified; the stage-05 configuration stays.** And the reason is the
one the design predicted:

The arms that win on conventional scIB are exactly the arms that merge the censored
plasma populations. `harmony_sample` posts the best batch *and* bio scores while mixing
plasma **13.5×** harder than the incumbent; `harmony_cohort` reaches **20.2×**. The
incumbent, meanwhile, is the worst batch corrector of the seven and by a wide margin the
least depth-encoding (R² 0.369 against 0.51–0.69) and the least plasma-merging. **A
standard global scIB ranking would have selected `harmony_sample`.** That the incumbent
was entered as its own arm rather than assumed is what made this visible at all. Full
table: `docs/stage-results.md`.

Two caveats recorded rather than buried: the depth criterion did all the gating (though
three arms fail the overcorrection criterion independently of it), and **no integration
method can restore cells that were never deposited** — Stage 08 still owes its
truncate-all-cohorts-at-10,000 sensitivity analysis whatever wins here.

---

## Stage 06 — Cell type annotation (`notebooks/06_annotation.ipynb`; env: `mm-annotation`)

**Three annotation methods are run, compared, and chosen between per cell-type class.**
An earlier plan said "`celltypist` and/or marker-panel scoring" — an `and/or` sitting in
the middle of a pipeline is an unmade decision, and left alone it gets settled implicitly
by whichever method happens to run first. This stage makes it explicitly and on evidence.

**What the comparison is actually judging** is not "which annotation is better in
general". This stage feeds exactly three things — the plasma-cell boundary for Stage 07,
which sets the metric's denominator; T/NK/myeloid *purity* for Stage 08, because one
plasma cell leaking into the "confidently antigen-negative" population inflates the
ambient noise floor and biases every antigen call downstream; and T/NK abundance for
Stage 11. Fine-grained subtypes are a bonus and must never be why a method wins.

**Method A, manual:** marker-panel scoring at cluster level, not per cell — clustering
absorbs dropout, which matters at 1,162 median genes per cell. Plus a Wilcoxon DE pass
per cluster, **not optional**, because it is the only step that can reveal a population
the seven-class panel doesn't cover at all. **Method B, CellTypist**, run over the
*existing* Leiden partition so the methods are directly comparable — same clusters,
different labels. **Method C, SingleR** against a sorted-hematopoietic reference, chosen
to cover CellTypist's predictable blind spot rather than duplicate its strengths: the
`Immune_All_*` models are immune-only, so erythroid and HSPC — real marrow populations
simply absent from an immune reference — are where automated annotation is most likely to
fail here. That costs an R bridge and its own environment, the same isolation `env-qc`
already applies to scDblFinder.

**Two caveats about plasma cells, pointing opposite ways.** The automated references
contain *normal* plasma cells only, so malignant PCs get labelled "plasma cell" — which
is fine, because telling malignant from normal is Stage 07's job, not this one's. Don't
expect an automated method to find the tumour, and don't hold it against it. But the same
fact cuts the other way: with no malignant class in the reference, a heavily aneuploid
clone with an odd transcriptome could be labelled something else entirely, or split
across labels. So the plasma-cell check has to be run **on myeloma marrows
specifically**, not just on the healthy controls — otherwise a systematic failure on
precisely the cells this project measures would sail through unnoticed. Stage 05's three
cohort-specific plasma clusters make that check more important, not less.

**The decisive test** is the one a reader can check at a glance: take the *manual*
marker panel and plot it grouped by each *automated* method's labels. If CellTypist's
T cells are CD3D-high and its plasma cells MZB1/SDC1-high straight down the panel, the
automated labels already encode everything the manual panel encodes, and manual
annotation is adding labour rather than information.

**The agreement numbers are concordance, not accuracy** — a distinction that changes how
they should be read. The manual annotation is a third opinion derived from the same
expression matrix, not ground truth, so "F1 against manual" measures agreement. And the
two automated methods agreeing with each other, while the strongest evidence available
because they were trained independently on different references, is still agreement
between two references that share canonical marker biology and may share its blind
spots. That is why the marker-coverage test above is the load-bearing evidence and can
veto a class regardless of concordance: a label set can be perfectly self-consistent and
biologically wrong.

**The thresholds are declared before looking**, which is the whole point; otherwise
"pick the best" quietly becomes "justify whichever looked tidier". The choice is made
**per class**, not once for the whole stage — the methods are expected to fail on
*different* classes, so a single verdict would throw away good labels to punish an
unrelated weakness. The per-class bars are in the main project document; the outcome and
its numbers go to `results/06_annotation/annotation_decision.md`. Everything downstream
then reads one column, `cell_type`, never needing to know which method produced it —
which is what makes the comparison reversible later.

**Identity and state are kept as separate axes.** A cell has one identity but can be
running several programs at once, so cell cycle, interferon response, antigen
presentation (`B2M` and HLA — `B2M` loss being a real, *competing* immune-escape route
in myeloma, though one CAR-T is indifferent to since it is MHC-independent), UPR and
stress are scored as **continuous values**, never folded into the identity label. A
cycling plasma cell is a plasma cell with a high cell-cycle score, not a "Cycling" cell
type. Collapsing those into categories would throw away exactly the intermediate cells
Stage 10 later needs.

Per-patient composition — malignant-PC fraction of the marrow (tumor burden) and
T/NK/myeloid abundance — is a **first-class output**, not a by-product: tumor burden is
context for the final packet, and T/NK abundance is the main confounder for Stage 11's
central claim. Any group comparison uses `scCODA` rather than a per-cell-type proportion
test, because proportions are compositional (they sum to one, so one population rising
mechanically pushes the others down) and naive tests on them are anticonservative.

---

## Stage 07 — Malignant plasma cell identification (`notebooks/07_malignant_plasma.ipynb`; env: `mm-core`)

Subset to plasma cell clusters. Clustering alone can't separate malignant from residual
normal plasma cells — clonality can. Score kappa (`IGKC`) vs. lambda (`IGLC1-7`)
restriction per cell; the per-patient dominant restriction class (>90% in an involved
marrow) marks malignant cells, minority-restriction cells are residual normal plasma
cells. Prefer actual scVDJ-seq/BCR clonotype calls over the restriction proxy if any
sample turns out to have them.

Everything else about this stage is driven by one fact: **it defines the denominator of
the headline metric**, so its mistakes don't stay local — they land directly in
`frac_double_negative`.

**The restriction call is ratio-based, not presence/absence.** Immunoglobulin
transcripts are the most ambient-contaminated genes in this entire tissue — plasma cells
pour Ig mRNA into the droplet background, which is the same reason Stage 04 flags
ambient RNA as a real concern here. A presence-based "does this cell have `IGKC`?" call
is therefore far noisier than it appears, because a chunk of that signal isn't the
cell's. A kappa:lambda **ratio** is robust to a shared additive background in a way a
presence call simply isn't.

**`infercnvpy` is required rather than optional**, using minority-restriction cells as
the per-patient normal reference, with the light-chain/CNV **agreement rate reported as
a stage output**. The reason it can't stay optional: residual normal plasma cells
express *less* BCMA and GPRC5D than malignant ones, so every normal plasma cell wrongly
called malignant pushes the escape fraction up. A second, independent call is the only
way to catch that. A poor agreement rate invalidates Stage 08 and stops the pipeline —
it isn't something to note and move past.

**And the two lines of evidence are not collapsed into one bit.** With two independent
calls there is no reason to throw the disagreements away, so cells carry a confidence
tier instead of a boolean, and Stage 08's headline result is re-run on the
highest-confidence cells alone as a sensitivity check — a much stronger statement than
the agreement rate by itself. The tier definitions are in the main project document; the
one that matters conceptually is that CNV being *not evaluable* is not the same thing as
CNV being negative.

**The normal bone marrow samples become a negative control.** `BM2/4/5/6` and the `ND_*`
samples get run through the identical calling logic. (That these eight are genuinely
donor marrow is confirmed, not inferred from their names: the GEO metadata gives all
eight `source_name = "Donor BMMC"` and no `diagnosis` characteristic, while the other 54
read "Multiple myeloma (MM)". They also span both reference builds, so the control
doubles as a build check in a population carrying no clone to confound it.) Normal
marrow is polyclonal, so the correct answer is *no malignant cells found*. If the method
discovers a clone in healthy marrow, the method is broken and everything downstream is
worthless. This is the cheapest strong validation available for the project's most
method-dependent step, and it uses samples that were already downloaded and doing
nothing.

---

## Stage 08 — Antigen scoring + dual-antigen escape fraction (`notebooks/08_dual_antigen_escape.ipynb`; env: `mm-core`)

Per malignant cell, positivity for BCMA (`TNFRSF17`), GPRC5D, and backup candidates
(`SLAMF7`, `FCRL5`) — using the empirical ambient-noise-floor threshold derived from
confidently antigen-negative cell types (T/NK/myeloid), not a naive `>0` call. Each cell
classifies into `dual_positive`/`BCMA_only`/`GPRC5D_only`/`double_negative`.

**The core novel metric**: per patient, `frac_double_negative` = the fraction of
malignant cells negative for both antigens at once, computed at baseline before any
treatment — reframing the antigen-escape question from the literature's usual
before/after-relapse framing into a pre-treatment risk score.

Per-patient aggregation runs on the S1 mapping: **41 patients over 53 in-cohort
samples**, reproducing the paper exactly. This was a genuine blocker until S1 arrived —
a naive rule gave 43 patients from 54 samples, and the risk was one patient's cells
splitting across duplicate entries, each yielding its own wrong partial escape
fraction.

### Defending the metric

A fraction of zeros is the most fragile thing you can compute from scRNA-seq, and the
whole of the rest of this stage exists because of it.

The two errors point **opposite ways**. Ambient RNA makes a truly negative cell look
faintly positive, *deflating* the escape fraction. Dropout does the reverse — a cell that
genuinely expresses the antigen reads as zero because the transcript was never captured,
*inflating* it. Dropout is the bigger problem here for a specific reason: **`GPRC5D` is a
low-abundance GPCR transcript** and this cohort's median cell carries only ~1,162
detected genes, so a large share of "GPRC5D-negative" calls are going to be technical
rather than biological. BCMA-negative and GPRC5D-negative calls are therefore not equally
trustworthy, and the writeup has to say so.

So the stage reports a **range, not a number**. Four checks do the bounding — a
threshold sensitivity band, a depth-regression falsification test, an
expression-matched false-negative floor, and per-patient bootstrap intervals; the
specifications are in the main project document. Two of them have reasoning that isn't
obvious from the specification.

**Why the minimum cell count is not a round number.** Cell counts vary about fifteen-fold
across this cohort, so some patients simply cannot support the claim. At n = 50 malignant
cells a *single cell is 2%* — so 1%, 2% and 3% escape are not distinguishable, and
ordering patients across that range is noise with a decimal point on it. A 5% population
yields ~2.5 expected DN cells at n=50, ~5 at n=100, ~10 at n=200, so 50 is a floor rather
than an answer. The procedure matters as much as the number: inspect the distribution,
fix the threshold, *then* look at the results, and name the excluded patients rather than
quietly dropping them.

**Why the bootstrap level is a real question.** A flat cell-level bootstrap treats
sample-level batch variation as biological spread and reports intervals that are too
narrow. But the fix is not "resample at every level": a confidence interval *for patient
A* conditions on patient A, so patient is fixed rather than random there, and resampling
patients would answer a different question. One consequence is worth stating with the
results rather than burying — for the many single-sample patients this collapses to a
plain cell bootstrap, blind to sample-level variation, so their intervals are optimistic
relative to the eight patients who contribute more than one sample.

### Co-negativity enrichment — the key derived metric

The escape fraction alone can't distinguish two clinically different tumors. Per patient,
build the 2×2 of BCMA± × GPRC5D± over malignant cells and compare observed
double-negatives against the independence expectation `E[DN] = P(BCMA⁻) × P(GPRC5D⁻)`,
reporting the **co-escape enrichment ratio** `observed / expected`.

This splits three facts the single metric fuses: how often each antigen is individually
absent, how many cells are DN, and whether the *same* cells are disproportionately losing
both. A patient at 6% DN against a 6% independence expectation has two unrelated partial
failures; one at 6% DN against a 1.5% expectation has a coordinated antigen-low
phenotype, and that second patient is the one Stage 10 goes looking for a mechanism
behind.

**What it does not mean.** Enrichment does *not* say a second binder fails to help. The
arithmetic: adding GPRC5D to BCMA moves the uncovered fraction from `P(BCMA⁻)` to
`P(BCMA⁻ ∩ GPRC5D⁻)` — at 30% BCMA⁻ / 20% GPRC5D⁻ under independence, 30% → 6%; with
enrichment pushing DN to 15%, 30% → 15%. Less than independence promised, but still
halving the escape population. Enrichment measures **how much of the pair's expected
complementarity is eroded by correlated loss**, not whether the second target is worth
adding.

**So report incremental coverage gain beside it**, which is the quantity a target
decision actually turns on and is free off the same 2×2:

    gain from adding GPRC5D to BCMA  =  P(BCMA⁻)   − P(BCMA⁻ ∩ GPRC5D⁻)
    gain from adding BCMA to GPRC5D  =  P(GPRC5D⁻) − P(BCMA⁻ ∩ GPRC5D⁻)

Enrichment is a claim about biology (is loss correlated); incremental gain is a claim
about value (what the second target buys). A patient can score high on both — they are
not in tension — which is exactly why they stay separate columns.

**The null must be depth-conditioned, or this measures library size.** Dropout is a
per-*cell* property: a shallow cell reads zero for *both* genes, so depth heterogeneity
alone produces positive BCMA⁻/GPRC5D⁻ association. A permutation that shuffles labels
freely within a patient destroys the depth↔label coupling and will report co-escape
enrichment on data with no biological co-occurrence whatsoever — an artifact pointing in
exactly the direction the project hopes to find, which is the worst kind. So the null is
permuted **within depth strata**, and the unconditioned ratio is reported alongside: the
gap between the two *is* the artifact, quantified.

Read whatever this finds against Supplementary Table S3, where the sample-level
`GPRC5D`×`TNFRSF17` correlation is +0.62 in MMRF, +0.54 in WashU 2 and −0.09 in
WashU 1 — a sign flip tracking the cohort depth ordering exactly.

### The detection curve — and what it cannot deliver

It is tempting to compute a "dropout-adjusted expected DN" as
`Σ_i P_i(BCMA⁻) · P_i(GPRC5D⁻)` and call it corrected. **That is circular.** Multiplying
the marginals assumes exactly the independence the co-escape test exists to interrogate,
so a tumor with genuinely correlated loss would be "corrected" toward the null it
violates. The computation is still worth having — as the *technical baseline* the
observed value is compared against, never as a corrected truth.

Hence the stage's most important negative statement: **no dropout-corrected DN point
estimate is produced.** Dropout is *bounded* by the four checks above, not corrected, and
observed DN remains the point estimate reported as an interval. That is a stronger
position than shipping a number whose correction rests on an assumption the project is
simultaneously testing.

The detection curve itself is still built — detection probability against cell depth and
gene mean, fit on the expression-matched controls — because it makes the
depth-conditioned null quantitative rather than rank-based. A genuinely corrected DN
would need a latent-class model over the four true states, fit by EM over the observed
2×2, which estimates the joint *without* assuming independence. Real work, not on the
critical path, filed so it is neither reinvented casually nor mistaken for an
oversight.

**Imputation and denoising (MAGIC, scVI, ALRA) are forbidden for positivity calls.** They
manufacture low-level expression by borrowing from neighboring cells, and whether a
transcript is genuinely absent is the entire question being asked. Model the uncertainty;
don't fill it in.

### From two antigens to a coverage matrix

`SLAMF7` and `FCRL5` get promoted from "backups" to a real deliverable. For every pair
and triple across a seven-antigen panel, compute how much of each patient's clone would
be left uncovered. This answers the question a target-strategy audience actually asks and
the two-antigen metric structurally cannot: *is BCMA+GPRC5D the best pair for this
patient, or would BCMA+FCRL5 cover more of their tumor?*

Coverage is read against normal *marrow* expression from Stage 09 rather than maximized
on its own — a target covering the whole tumor that also hits normal marrow plasma cells
is not a better target. The columns stay separate rather than collapsing into a utility
score: a weighted `coverage − λ · exposure` needs a principled λ that does not exist, and
inventing one would encode a clinical judgement the data cannot supply while hiding the
inputs a reader could disagree with.

### The confounder to watch here

The deposit UMI ceiling from Stage 04 lands squarely on this metric: it censors a band
enriched 20–70× for `GPRC5D` in 36 of 54 myeloma samples, inflating
`frac_double_negative` for WashU — in the project's own direction of interest. So this
stage owes a **truncate-all-cohorts-at-10,000 sensitivity analysis** alongside the
threshold band and the depth regression. If the patient ordering survives that
truncation the metric is robust to the censoring; if it does not, the framing has to say
so plainly.

---

## Stage 09 — Escape robustness (`notebooks/09_bulk_validation.ipynb`; env: `mm-core`)

Everything in this stage exists to answer one question: *how do you know your escape
fractions are real?*

**Matched bulk RNA-seq.** GSE223061 — the matched bulk data — was downloaded at the
start of the project and then never used by the plan. It holds 29 usable bulk samples,
of which **26 have an exact single-cell match**. For those, compare pseudobulk
`TNFRSF17` and `GPRC5D` against bulk TPM, **cohort by cohort**: MMRF bulk is CD138+
sorted and pairs with *malignant-cell* pseudobulk, while WashU cohort 1 bulk is unsorted
whole marrow and pairs with *whole-sample* pseudobulk. Pooling the two would be a real
error rather than a nicety — correlating malignant-PC pseudobulk against unsorted bulk
measures tumour burden, since dilution by non-plasma cells scales with how much tumour is
present, and that would corrupt 10 of the 26 comparisons in a direction correlated with
the metric itself.

Agreement means the antigen quantification is technically credible: an orthogonal assay,
generated independently, agrees with the per-cell calls. Disagreement in one specific
direction — **bulk-positive where single-cell reads zero** — is direct quantified
evidence of dropout, and feeds straight back into Stage 08's false-negative floor. Either
outcome is informative, which is what makes this worth doing.

**But bulk validates antigen *abundance*, not the escape fraction**, and the distinction
is not pedantic. Bulk averages every cell in the sample together, destroying the joint
per-cell distribution the metric is built on. A tumor that is 50% BCMA⁺GPRC5D⁻ and 50%
BCMA⁻GPRC5D⁺ shows healthy bulk expression of *both* genes while containing zero
dual-positive cells — and the converse misreads are equally available. So every output
here is phrased as orthogonal validation of antigen abundance and of whether the
single-cell antigen-negative calls are plausible, never as validation of
`frac_double_negative`, which has no orthogonal check available anywhere in this project.

**Normal plasma-cell antigen baseline — marrow expression context, not a safety axis.**
Do *normal* plasma cells, from the healthy `BM*`/`ND_*` marrow, express BCMA and GPRC5D?
This is genuinely worth having: BCMA is broadly expressed on normal plasma cells and
B-lineage cells, and the malignant-versus-normal-PC contrast is what makes a coverage
number interpretable rather than absolute.

What it is *not* is a safety axis. GPRC5D's clinically decisive off-tumor site is
keratinized tissue — the nail, skin and taste toxicity seen with talquetamab — which a
marrow dataset cannot observe at all, and expression is not toxicity in any case. Three
things stay separate in the writeup: **tumor coverage**, **normal marrow expression**
(measured here), and **known extra-marrow liabilities** (cited, not measured). A real
`coverage − λ · exposure` utility score would need a normal-tissue atlas such as GTEx or
HPA, and is a future extension rather than a claim from this data.

**The label-permutation null lives at Stage 08, not here.** Shuffling antigen labels
within a patient while preserving each antigen's marginal negative rate cannot produce
"the metric with no signal in it" — the marginals are exactly what it holds fixed. What
it actually tests is whether BCMA-negativity and GPRC5D-negativity *co-occur* beyond
independence, which is a sharper question than the one it was written for, and it belongs
with the co-negativity enrichment test.

---

## Stage 10 — Escape subclone + phenotype (`notebooks/10_dn_coherence.ipynb`; env: `mm-core`)

The stage carrying the project's actual scientific payoff rather than another robustness
check.

**Is the double-negative population structured, or just scattered noise?** This is the
most important question the project can ask of its own metric. "3% of this patient's
malignant cells are double-negative" and "this patient carries a pre-existing 3%
resistant subclone" sound like the same statement and are not — only the second predicts
that therapy will *select* for those cells, which is the entire clinical premise of
measuring baseline escape in the first place.

**But transcriptional clustering does not establish clonality**, and it is worth being
precise about why, because the intuitive version of this test is wrong in both
directions. A transcriptionally coherent group of cells can come from cell cycle, stress,
interferon tone, metabolic state, sequencing depth or sample-prep batch just as easily as
from a genetic subclone. And conversely, a genuine genetic clone need not form a tidy
transcriptional island at all. So the question is **DN coherence**, tested at three
escalating levels with the claim escalating alongside the evidence:

| Level | Question | Licenses saying |
|---|---|---|
| **1 — enrichment** | Are DN cells non-randomly located in malignant transcriptional space? | "non-randomly distributed" |
| **2 — transcriptional coherence** | Do DN cells share a reproducible program? | "an escape-associated **state**" |
| **3 — genomic coherence** | Do DN cells preferentially occupy a CNV-defined subclone? | "an escape-associated **subclone**" |

**Only level 3 licenses the word "subclone."** Levels 1+2 alone are an escape-associated
*state* — still the thing that separates a structured 3% from a scattered 3%, and still
worth reporting, but not a claim about pre-existing genetic clones under selection.

**And a negative at level 3 is not evidence of absence.** Resolving CNV substructure
*within* a single patient's clone is far harder than separating tumor from normal, and at
this depth — ~1,521 median genes per cell in the plasma compartment — it will often be
underpowered. Level 3 reports **supported / not evaluable**, with per-patient CNV
resolution stated, never "no CNV subclone". Treating an underpowered null as a negative
would systematically understate exactly the risk this project exists to measure.

The per-patient level attained is reported next to the escape fraction. This runs on the
per-patient un-integrated embedding, never the Harmony one, for the reason given at Stage
05. And the level-1 permutation carries the same depth-stratification requirement as
Stage 08's co-negativity test: shallow cells are both likelier to read double-negative
*and* to sit together in low-dimensional space, so an unconditioned enrichment test will
see depth structure and report it as biology.

**What else is different about the escape cells?** Stage-10 differential expression used
depth-matched patient pseudobulks. For each gene, paired patient-level DN-versus-comparator
log-fold changes were tested with a two-sided Wilcoxon signed-rank test, followed by
Benjamini–Hochberg correction. **Patient is the biological unit.** `decoupler` was used
separately for pathway and TF activity. Per-cell DE
tests treat cells as independent replicates and badly inflate the false discovery rate.
Program scores cover MYC and OXPHOS alongside the Stage 06 programs — MYC because its activation is
a recognized myeloma progression event, which makes "is the escape population MYC-high?"
a substantive question rather than a generic one.

**One hypothesis is pre-registered: the γ-secretase axis** (`NCSTN`, `PSEN1`, `APH1A`,
`APH1B`, `PSENEN`). γ-secretase physically cleaves BCMA off the cell surface, which is
why γ-secretase-inhibitor plus BCMA CAR-T combinations are in active clinical
development. If escape cells turn out to be γ-secretase-high, that's a directly
actionable, mechanistically grounded finding rather than a descriptive one. It's written
down here, before looking, specifically so it stays a hypothesis test instead of becoming
a post-hoc story.

**The TC (Translocation/Cyclin D) subgroup is assigned per patient** from pseudobulk,
asking whether dual-antigen escape risk concentrates in a molecular subgroup. It is a
**transcriptional proxy, never a translocation call** — `NSD2`/`FGFR3` overexpression is
*consistent with* t(4;14), not a breakpoint detection — and since S1 carries no
cytogenetics column, nothing in this deposit can confirm it. Used descriptively; at ~41
patients across ~5 classes an association test would be underpowered.

---

## Stage 11 — Exploratory immune context (`notebooks/11_immune_context.ipynb`; env: `mm-communication`)

Runs LIANA+ using the CellChat-algorithm method specifically (LIANA+ natively
reimplements it), or the full consensus rank-aggregate across
CellPhoneDB/CellChat/Connectome/NATMI/SingleCellSignalR — arguably stronger evidence than
any single method, since it is cross-validated across independent approaches. The question
is whether high-escape-risk patients also show weaker NK/T-cell engagement signaling
toward malignant plasma cells.

**This stage is explicitly exploratory**, and ranks last in the project's scientific
hierarchy — for reasons of power, not interest: ~41 patients against hundreds of
ligand-receptor pairs, with a confounder correlated with the predictor. Left un-demoted,
this is the stage that turns a focused antigen project into a kitchen-sink scRNA project.

The design also has to avoid **pseudoreplication**. Splitting patients into
high/low tertiles and pooling their cells treats thousands of cells from one person as
thousands of independent observations, inflating significance essentially arbitrarily —
the real sample size is the number of patients (~41), not the number of cells
(~173,000). So interaction scores are computed **per patient**, `frac_double_negative`
enters as a **continuous predictor** rather than being binned (no arbitrary cutoff, and
strictly more power than discarding the middle third), and **the confounder is controlled
explicitly**: high-escape patients might simply have fewer T and NK cells to begin with,
in which case "weaker immune signaling" is a composition artifact and the finding
evaporates. Cell-type abundance from Stage 06 goes in as a covariate, and/or cells are
downsampled to equal per-type counts per patient. That one is fatal to the stage's claim
if ignored, so it isn't treated as optional polish.

---

## Stage 12 — Final synthesis (`notebooks/12_final_synthesis.py`; env: `mm-core`; output `results/12_final_synthesis/`)

**RUN, COMPLETE (2026-08-27T07:13:56Z, commit `5bbecbc`).** The final stage, assembling
everything upstream of it — 29 frozen artifacts, hash-verified before use — into a
synthesis a target-strategy discussion can actually use. It is deliberately
**synthesis only**: no new statistical test, no fitted model, no new threshold, no score.
Every number it reports already existed in a frozen artifact.

**No risk tiers, and no ranking — a decision made after inspecting the evidence, not
before.** The plan going in was "robust high escape / uncertain / robust low escape."
Stage 12 built the full cross-tabulation of the frozen tier (08/09b), Level-1 structure
and Level-2 phenotype axes first, and found only 5 of the 18 possible combinations
occupied — with the four measurement-`robust-high` patients and the four
Level-1-`SUPPORTED` patients forming **disjoint sets**. Any single categorical label over
that structure would either relabel three existing columns with risk-sounding language
that the genomic axis (all 32 patients `NOT_EVALUABLE`) cannot support, or merge cells
arbitrarily, including a one-patient singleton. **The decision: no categorical patient
label, no rank ordering, no composite score of any kind.**

**What replaced it: a 32-row evidence matrix with six axes kept as separate columns** —
measurement robustness (08/09b), DN structure (10 Level 1), DN phenotype (10 Level 2),
genomic evidence (10 Level 3), immune context (11/11b), and multi-antigen coverage (08c)
— plus deterministically generated interpretation text (never freeform prose) and a
12-claim ladder stating exactly which claims the frozen evidence does and does not
support. Keeping the axes separate is the point, and it is now a demonstrated one rather
than a design preference: a patient can be strong on measurement and unremarkable on
structure, or the reverse, and folding them into one number would have erased exactly
that disagreement — which turned out to be the project's central Stage-12 finding.

**Clinical annotation comes from S1**, which supplies ISS stage, treatment and
time-to-progression per patient, plus disease stage for WashU cohort 1 only — MMRF and
WashU 2 have no per-sample stage and it is not imputed. **S1 carries no cytogenetics**,
so there is no karyotype column and Stage 10's TC class enters explicitly as a
transcriptional proxy, unchanged through to the final synthesis.

**The bias-direction table** — ambient, dropout, malignant-call error, deposit UMI
censoring, mRNA-vs-protein, each with its sign on the metric — is included as one of the
29 consumed inputs, feeding the uncertainty register rather than left implicit.

**The mRNA-versus-protein limitation is stated explicitly and mechanistically, not
quantified.** CAR-T binds surface protein; this analysis measures transcript. BCMA is
actively shed from the surface by γ-secretase, and GPRC5D transcript correlates
imperfectly with surface density. No published CITE-seq/flow calibration for these two
antigens was incorporated, so the limitation sits in the uncertainty register as a stated
gap — this is the first question a target-strategy audience will ask, and it is answered
directly rather than left for them to catch, but it is not resolved.

**The headline result, verbatim from `results/12_final_synthesis/stage12_summary.md`:**
observed transcript-level BCMA/GPRC5D double-negativity is common at baseline (median
0.335 across 32 patients), but its magnitude is dominated by measurement limitations
rather than demonstrable biology. Depth conditioning removes most apparent co-loss
enrichment and most apparent DN structure; the pre-registered γ-secretase hypothesis is
not supported; genomic subclone evidence is not evaluable for any patient; and neither
immune composition nor ligand-receptor communication survives correction. No patient
shows convergent measurement-robust, structurally-supported and phenotype-supported
evidence, and the evidence axes are substantially discordant — which is reported as the
result, not smoothed over into a single number.

---

## Phase 2 — External validation on GSE117156 (not started; sequenced strictly after Phase 1 completes)

Independently re-runs the same pipeline shape (QC → malignant calling → antigen scoring →
escape fraction) against GSE117156 (Ledergor et al. 2018, *Nat Med* — 51,840 cells, 11
healthy controls + 29 MM patients spanning asymptomatic disease through post-treatment
MRD) as a second, independent cohort, testing whether the core finding replicates beyond
this one dataset and technology.

Full reasoning, acquisition method, and the explicit no-merge constraint (MARS-seq vs.
10x — a platform difference, not a correctable batch effect) are in the main project
document's Phase 2 section, not duplicated here.
