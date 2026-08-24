# Multiple Myeloma Dual-Antigen (BCMA/GPRC5D) Coverage & Escape Risk
## Pipeline Walkthrough (Python rebuild)

**Objective:** For each multiple myeloma patient in GSE223060, quantify the fraction
of malignant plasma cells that would evade BOTH a BCMA-directed and a
GPRC5D-directed CAR-T/TCE, then relate that "dual-antigen escape fraction" to the
immune microenvironment via LIANA+ (CellChat's algorithm, natively reimplemented in
Python). Output is a per-patient **risk-tier table** — robust-high / uncertain /
robust-low, with co-escape enrichment, incremental coverage gain and DN coherence as
separate columns — usable for a single- vs. dual- vs. sequential-target CAR-T strategy
discussion. Deliberately not a rank ordering: the confidence intervals overlap too
heavily for ordinal positions to carry meaning (see Stage 12).

**This document is a narrative walkthrough of the pipeline's logic and reasoning,
not a copy of the code.** For exact implementation, read `src/mm_escape/` directly
— this doc intentionally does not duplicate code inline, since exactly that
duplication caused this file to go stale once already during the R build (an
earlier version described a monolithic script that had already been split and
deleted). For current execution state, read `RESUME_HERE.md`. For the settled
technical decisions and their reasoning, read `CLAUDE.md`. This file answers "what
does the pipeline do and why," not "what state is it in right now."

**This is a from-scratch Python rebuild.** An earlier R/Seurat version of this
project reached: data acquisition fully solved, QC/doublet-removal run on the full
cohort, integration not yet run. No R code carries over; all data-format knowledge
does — see `CLAUDE.md`'s Data section for the full detail (file-naming quirks,
reference mismatch, patient-mapping gap). The R build has been removed from the
working tree (preserved in git under the `r-build-snapshot` tag); only `raw/` and
`scripts/01-03` carry over.

**Numbering note:** stage numbers below match notebook filenames and their
`results/NN_*/` output directory one-to-one (`notebooks/04_qc.ipynb` ->
`results/04_qc/`, etc.). This is deliberate — it's the actual mechanism for tracking
pipeline order end to end. **Every stage from 01 to 12 is a notebook**; there is no
stage you cannot open and step through.
**Number order is execution order** — 04 → 05 → 06 → 07 → 08 → 09 → 10 → 11 → 12,
with no exceptions. Nothing runs out of sequence.

**Scope expansion (2026-08-20).** A design review added a robustness layer and two
new stages, and the sequence was renumbered so the new stages sit where they
actually run: escape robustness at 09 and the subclone/phenotype analysis at 10,
which pushed cell-cell communication to 11 and the decision packet to 12 (the
packet consumes everything upstream, so it stays last). The reasoning behind each
addition lives in `CLAUDE.md`; what follows is the narrative version. The
through-line: the headline metric is a *fraction of zeros*, the noisiest thing
scRNA-seq measures, and the original plan bounded only one of its two error
directions.

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

**Status: run and confirmed 2026-08-20.** All 62 samples `triplet-ok`, nothing
re-downloaded, manifest reproduced exactly.

Notebook 03 does more than the script it wraps, which is the point of having it: it
reports the Cell Ranger reference split, previews the symbol-harmonized gene intersection,
and runs the required-gene assertions. Those assertions are what surfaced the `NSD2`
symbol drift described under Stage 05 — a defect two earlier builds of this project had
walked straight past.

---

## Stage 04 — QC + doublet detection (`notebooks/04_qc.ipynb`, `src/mm_escape/qc.py`; env: `mm-qc`)

Loads each sample via a custom loader (`src/mm_escape/io.py`) rather than
`scanpy.read_10x_mtx()`, since this archive's filenames (`counts.mtx`, single-column
`genes.tsv`) don't match what that convenience function expects.

Computes QC metrics via `scanpy.pp.calculate_qc_metrics`, including
`pct_counts_in_top_20_genes` — a covariate the R build didn't track. Applies
MAD-based (median absolute deviation) outlier filtering rather than the R build's
fixed thresholds (`nFeature_RNA` 200-6000, `percent.mt` < 15), following
`sc-best-practices.org`'s QC chapter method: flag outliers beyond 5 MADs on
`log1p_total_counts`, `log1p_n_genes_by_counts`, and `pct_counts_in_top_20_genes`,
with a separate, tighter check on `pct_counts_mt`. The tutorial's exact numeric
cutoffs are for a different (healthy PBMC/BMMC) dataset and are not copied
verbatim — the actual thresholds get re-derived and documented against this
cohort's own distributions.

**Run 2026-08-24: 204,040 → 172,940 cells (84.8% kept),** thresholds in
`results/04_qc/qc_thresholds.csv`, one checkpoint per sample keeping **every** barcode
with `obs["keep"]` set. Two things the data changed, both departures from the
paragraph above:

- **`pct_counts_in_top_20_genes` is computed and reported but does not filter.**
  A 5-MAD band on it flags 17% of MMRF and 15% of WashU 1 against 3% of WashU 2, and
  the flagged cells are two populations: `IGKC` at ~25% of counts (plasma cells) and
  haemoglobin at ~32% (erythroid debris). The plasma-cell half is the project's
  subject — `TNFRSF17` is detected in **21.8%** of that decile against **0.8%**
  elsewhere. An Ig-dominated library is a plasma cell's *normal* state; filtering on
  it would delete antigen-**positive** malignant cells and inflate the escape
  fraction. The metric is kept as an ambient-Ig handle for Stage 08 (which needs one,
  since SoupX cannot run on this deposit) but is not allowed to remove cells.
- **`pct_counts_mt` is one-sided.** A cell with unusually *few* mitochondrial reads is
  not low quality, so only the upper bound applies.

**And the deposit turns out to be pre-filtered, differently in each cohort:** WashU 1
and 2 at <10,000 UMIs and <20% mt, MMRF uncensored on UMIs but cut at <10% mt, donors
uncensored on both. See Stage 08's note — the WashU UMI ceiling censors a band
enriched 20–70× for `GPRC5D` and is a live confounder for the headline metric.

**Thresholds are derived per cohort, not pooled.** The three collection cohorts ran
different 10x chemistries (WashU 1 on 3′ v2 with no dead-cell removal; MMRF on v3.3
and WashU 2 / donors on v3.2, both with it) and differ ~1.9× in genes detected per
cell — MMRF 1916 vs. WashU 1 1023, sample-level medians. A single pooled MAD across
that spread would flag much of WashU cohort 1 as low-quality for what is a batch
difference, not a cell-quality one.

Doublet detection uses `scDblFinder` (R), called from Python via an isolated
`rpy2` bridge contained entirely within the `mm-qc` environment — kept as R
because it's the benchmarked best-performing method per `sc-best-practices`'s
cited benchmark (Xi & Li 2021), not swapped to a pure-Python alternative purely
for language purity.

`56203_1` is **repaired on load, not excluded** (reversed 2026-08-24). It was long
believed to come from an incompatible 22184-gene reference missing BCMA; it is in
fact a normal 33694-build sample whose `genes.tsv` write stopped mid-symbol at row
22185, putting `TNFRSF17` and `IGLC1/2/3` past the cut rather than absent from a
reference. `io.read_sample` substitutes the canonical column behind a prefix
assertion that raises rather than guessing. See `CLAUDE.md`'s Data section.

**Ambient RNA correction (SoupX/DecontX) is not run — not an oversight, a hard
constraint.** Both require the unfiltered Cell Ranger matrix (including empty
droplets) to estimate the background contamination profile, and only the filtered
per-sample matrices were ever deposited. (Raw reads do exist, under controlled access
on dbGaP — `phs000159` / `phs000748` — but no unfiltered matrices exist anywhere, so
the constraint holds either way.) The mitigation lives at Stage 08
(antigen scoring) instead — an empirically-derived noise floor rather than
formal ambient correction. Worth understanding for plasma-cell data specifically:
plasma cells secrete enormous quantities of immunoglobulin transcript, making
ambient contamination of `IGKC`/`IGLC` and potentially the antigen genes
themselves a real, non-hypothetical concern here — which is also why Stage 07's
light-chain call is ratio-based rather than presence-based.

**Ambient RNA is only half the problem, though.** It biases in one direction (true
negatives read as positive, *deflating* the escape fraction). Dropout biases the
other way — true positives read as zero, *inflating* it — and for this project
dropout is the larger of the two, because `GPRC5D` is a low-abundance GPCR transcript
and this cohort's median cell carries only ~2,044 detected genes. Neither bias may be
left unquantified; Stage 08 and Stage 09 handle the bounding.

Each sample's post-QC AnnData is checkpointed individually — mirrors the R
build's resumable-per-sample-checkpoint design, which is worth keeping regardless
of language — QC and doublet detection across 62 samples benefit from resumability
the same way in either stack. (The load itself is no longer the slow part: the Python
loader reads all 62 samples, 204,040 pre-QC cells, in about two seconds.)

---

## Stage 05 — Gene-space intersection, integration, clustering (`notebooks/05_integration_clustering.ipynb`; env: `mm-core`)

Before any concatenation: intersect gene sets across all retained samples (the
62 samples split across 33538- and 33694-gene Cell Ranger references — see
`CLAUDE.md`'s Data section for why this must be an intersection, never a union).
Assert all required marker/antigen genes survive; hard-fail naming the specific
missing gene(s) otherwise. `anndata.concat(..., join="inner")` performs the merge.

**But intersect *canonicalized* symbols, not raw ones.** The two references were built
against different HGNC symbol vintages, so some genes exist in **both** builds under
different names and a naive set intersection throws them away. Found the hard way when
stage 03's assertions failed on `NSD2` — which the older build calls `WHSC1`. The
others are `TENT5C`/`FAM46C` (a recurrently-deleted myeloma tumour suppressor),
`NSD3`/`WHSC1L1` (a genuinely different gene from NSD2 — don't conflate them), and
`ATP5F1A`/`ATP5A1`.

That `NSD2` case is the one that matters most: NSD2 is how t(4;14) is detected, and
t(4;14) is the highest-risk myeloma translocation. Left unharmonized, Stage 10's
molecular-subgroup call would quietly lose its highest-risk class — not error out, just
never report it. Harmonizing first recovers all four genes (22,164 → 22,168).

The wider lesson the assertions bought: **a required gene coming up missing means
"check for a legacy symbol", not "biologically absent".** The guess that the ~11k
symbols unique to each build were mostly annotation-version noise (`AC000032.1` versus
`AC000032.2`) turned out to be far too optimistic — joining on reconstructed Ensembl
IDs instead of symbols recovers **32,991 genes against 22,164**, because 11,140
intersected IDs simply carry a different symbol in each build. The symbol join was
also unsafe, not merely lossy: `TBCE` is a *different* Ensembl entry in each build, so
a symbol join silently merges two unrelated rows.

Normalize, select highly variable genes, PCA, then `harmonypy` integration keyed
on `patient_id` (with **both `n_genes_ref` and `cohort`** as additional covariates) to
correct for patient-of-origin batch effects. Leiden clustering, UMAP. Includes
diagnostic UMAPs colored by reference version and by cohort — if clusters visibly
track either instead of biology, the gene-space intersection failed to neutralize
processing batch and every downstream antigen call is suspect.

**The two covariates are not interchangeable, which is easy to get wrong.** The
reference build and the collection cohort cut across each other: two WashU cohort 1
samples sit on the newer 33538 build while the four `ND_*` donors sit on the older
33694 one. `n_genes_ref` captures which reference processed a sample; `cohort`
captures which chemistry and protocol generated it, and that is where the ~1.9× depth
difference lives. Neither substitutes for the other.

**Integration is deliberately confined to the immune compartment.** Harmonizing on
`patient_id` is right for T/NK/myeloid cells, which should look alike across
patients. It is actively risky for the tumor: the malignant clone is patient-private
by definition, so forcing patients together can blend distinct clones and erase the
heterogeneity the project exists to measure. So the integrated embedding is used for
immune annotation and clustering only; **all malignant subclustering happens per
patient, un-integrated** (stage 10 relies on this). The reassuring part, worth
stating out loud because it's the obvious objection: **per-cell antigen calls are
raw counts and never touch the embedding**, so no integration choice can distort the
escape fractions themselves.

---

## Stage 06 — Cell type annotation (`notebooks/06_annotation.ipynb`; env: `mm-annotation`)

**Three annotation methods are run, compared, and chosen between per cell-type
class.** The earlier plan said "`celltypist` and/or marker-panel scoring" — an `and/or`
sitting in the middle of the pipeline is an unmade decision, and left alone it gets
settled implicitly by whichever method happens to run first. This stage makes it
explicitly and on evidence, following `sc-best-practices.org`'s annotation chapter.
Note the env change: stage 06 runs under `mm-annotation`, not `mm-core`.

**What the comparison is actually judging.** Not "which annotation is better in
general" — this stage feeds exactly three things downstream, and the acceptance test is
weighted accordingly. Stage 07 needs the plasma-cell boundary right, because it sets
the denominator of the headline metric. Stage 08 needs T/NK/myeloid *purity*, because
those cells define the ambient noise floor — one plasma cell leaking into that
"confidently antigen-negative" population inflates the floor and biases every antigen
call downstream. Stage 11 needs T/NK abundance. Fine-grained subtypes are a bonus and
must never be why a method wins.

**Method A, manual.** Marker-panel scoring at cluster level (not per cell — clustering
absorbs dropout, which matters at ~2,044 median genes/cell) against the same panel used
in the R build: PlasmaCell (SDC1/CD38/MZB1/XBP1/IRF4), Bcell (MS4A1/CD79A/CD19), Tcell
(CD3D/CD3E/CD8A/CD4), NK (NCAM1/NKG7/GNLY), Myeloid (CD14/LYZ/ITGAM), Erythroid
(HBB/GYPA), HSPC (CD34/KIT). Plus a Wilcoxon DE pass per cluster — **not optional**,
because it is the only step that can reveal a population the seven-class panel doesn't
cover at all.

**Method B, CellTypist.** The `Immune_All_Low`/`Immune_All_High` models with majority
voting, run over the *existing* Leiden partition so the methods are directly
comparable — same clusters, different labels.

**Method C, SingleR** against a sorted-hematopoietic reference. Chosen deliberately to
cover CellTypist's predictable blind spot rather than duplicate its strengths: the
`Immune_All_*` models are immune-only, so erythroid and HSPC — real bone-marrow
populations that simply aren't in an immune reference — are where automated annotation
is most likely to fail here. This costs an R bridge and therefore its own environment,
the same isolation `env-qc` already applies to scDblFinder.

**Two caveats about plasma cells, pointing opposite ways.** The automated references
contain *normal* plasma cells only, so malignant PCs get labelled "plasma cell" — which
is fine, because telling malignant from normal is Stage 07's job, not this one's. Don't
expect an automated method to find the tumour, and don't hold it against it. But the
same fact cuts the other way: with no malignant class in the reference, a heavily
aneuploid clone with an odd transcriptome could be labelled something else entirely, or
split across labels. So the plasma-cell check has to be run **on myeloma marrows
specifically**, not just on the healthy controls — otherwise a systematic failure on
precisely the cells this project measures would sail through unnoticed.

**The decisive test** is the one a reader can check at a glance: take the *manual*
marker panel and plot it grouped by each *automated* method's labels. If CellTypist's
T cells are CD3D-high and its plasma cells MZB1/SDC1-high straight down the panel, then
the automated labels already encode everything the manual panel encodes, and manual
annotation is adding labour rather than information. Alongside that: confusion matrices
between all three methods, and per-class agreement scores. **The two automated methods
agreeing with each other is the strongest evidence available**, since they were trained
independently on different references — much harder to dismiss than either one matching
the manual panel.

**The thresholds are declared before looking**, which is the whole point; otherwise
"pick the best" quietly becomes "justify whichever looked tidier." Plasma cells need
F1 ≥ 0.95 (strictest, because they set the metric's denominator), T/NK/myeloid ≥ 0.90
(the noise floor depends on them), everything else ≥ 0.85. The choice is made **per
class**, not once for the whole stage — the methods are expected to fail on *different*
classes, so a single verdict would throw away good labels to punish an unrelated
weakness. Where neither automated method clears the bar for a class, that class falls
back to manual. The result, with its numbers, is written to
`results/06_annotation/annotation_decision.md`.

Everything downstream then reads one column, `cell_type`, and never needs to know which
method produced it — which is what makes the whole comparison reversible later.

**Identity and state are kept as separate axes.** A cell has one identity but can be
running several programs at once, so cell cycle, interferon response, antigen
presentation (`B2M` and HLA — `B2M` loss being a real, *competing* immune-escape route
in myeloma), UPR and stress are scored as **continuous values**, never folded into the
identity label. A cycling plasma cell is a plasma cell with a high cell-cycle score,
not a "Cycling" cell type. Collapsing those into categories would throw away exactly
the intermediate cells that Stage 10 later needs.

Per-patient composition — malignant-PC fraction of the marrow (tumor burden) and
T/NK/myeloid abundance — is a **first-class output of this stage**, not a
by-product. Tumor burden is context for the final packet, and T/NK abundance turns
out to be the main confounder for Stage 11's central claim, so it has to be measured
here. Any group comparison of composition uses `scCODA` rather than a per-cell-type
proportion test: proportions are compositional (they sum to one, so one population
rising mechanically pushes the others down), and naive tests on them are
anticonservative.

---

## Stage 07 — Malignant plasma cell identification (`notebooks/07_malignant_calling.ipynb`; env: `mm-core`)

Subset to plasma cell clusters. Clustering alone can't separate malignant from
residual normal plasma cells — clonality can. Score kappa (`IGKC`) vs. lambda
(`IGLC1-7`) restriction per cell; the per-patient dominant restriction class
(>90% in an involved marrow) marks malignant cells, minority-restriction cells
are residual normal plasma cells. Prefer actual scVDJ-seq/BCR clonotype calls over
the restriction proxy if any sample turns out to have them (check GEO supplementary
files — status unconfirmed either way).

This stage got three upgrades in the scope expansion, all driven by one fact: **it
defines the denominator of the headline metric**, so its mistakes don't stay local —
they land directly in `frac_double_negative`.

**The restriction call is ratio-based, not presence/absence.** Immunoglobulin
transcripts are the most ambient-contaminated genes in this entire tissue — plasma
cells pour Ig mRNA into the droplet background, which is the same reason Stage 04
flags ambient RNA as a real concern here. A presence-based "does this cell have
`IGKC`?" call is therefore far noisier than it appears, because a chunk of that
signal isn't the cell's. A kappa:lambda **ratio** is robust to a shared additive
background in a way a presence call simply isn't.

**`infercnvpy` is now required rather than optional**, using minority-restriction
cells as the per-patient normal reference, with the light-chain/CNV **agreement rate
reported as a stage output**. The reason it can't stay optional: residual normal
plasma cells express *less* BCMA and GPRC5D than malignant ones, so every normal
plasma cell wrongly called malignant pushes the escape fraction up. A second,
independent call is the only way to catch that. A poor agreement rate invalidates
Stage 08 and stops the pipeline — it isn't something to note and move past.

**The normal bone marrow samples become a negative control.** `BM2/4/5/6` and the
`ND_*` samples get run through the identical calling logic. (That these eight are
genuinely donor marrow is confirmed, not inferred from their names: the GEO metadata
gives all eight `source_name = "Donor BMMC"` and no `diagnosis` characteristic, while
the other 54 samples read "Multiple myeloma (MM)". They also happen to span both
reference builds, so the control doubles as a build check in a population that carries
no clone to confound it.) Normal marrow is
polyclonal, so the correct answer is *no malignant cells found*. If the method
discovers a clone in healthy marrow, the method is broken and everything downstream
is worthless. This is the cheapest strong validation available for the project's
most method-dependent step, and it uses samples that were already downloaded and
doing nothing.

---

## Stage 08 — Antigen scoring + dual-antigen escape fraction (`notebooks/08_antigen_escape_fraction.ipynb`; env: `mm-core`)

Per malignant cell, positivity for BCMA (`TNFRSF17`), GPRC5D, and backup
candidates (`SLAMF7`, `FCRL5`) — using the empirical ambient-noise-floor threshold
derived from confidently antigen-negative cell types (T/NK/myeloid), not a naive
`>0` call (see Stage 04's ambient-RNA note for why this matters more than usual for
this tissue type). Classifies each cell into `dual_positive`/`BCMA_only`/
`GPRC5D_only`/`double_negative`.

**The core novel metric**: per patient, `frac_double_negative` = the fraction of
malignant cells negative for both antigens at once, computed at baseline before
any treatment — reframing the antigen-escape question from the literature's usual
before/after-relapse framing into a pre-treatment risk score.

### Defending the metric

A fraction of zeros is the most fragile thing you can compute from scRNA-seq, and
the original plan bounded only one of the two ways it can go wrong. Both directions
now get handled here.

The two errors point **opposite ways**. Ambient RNA makes a truly negative cell look
faintly positive, which *deflates* the escape fraction — that one was already
documented. Dropout does the reverse: a cell that genuinely expresses the antigen
reads as zero because the transcript was never captured, which *inflates* it. Dropout
is the bigger problem here for a specific reason — **`GPRC5D` is a low-abundance GPCR
transcript**, and this cohort's median cell has only ~2,044 detected genes. A large
share of "GPRC5D-negative" calls are going to be technical, not biological.

So the stage reports a range, not a number:

- **A threshold sensitivity band.** Compute the metric under at least three calling
  rules — naive `>0`, the ambient noise floor, and a stricter floor. The claim being
  made isn't any single value; it's whether the **patient ranking holds** across all
  three (reported as pairwise Spearman ρ). A ranking that survives every threshold is
  a finding. One that doesn't is an artifact of a cutoff, and gets said so.
- **A falsification test for dropout.** Regress each patient's escape fraction
  against the median UMIs-per-cell of their malignant cells. If that slope is
  strongly negative, the metric is measuring sequencing depth rather than biology —
  and that has to be checked *before* the ranking is shown to anyone, not after.
  Malignant cells are also downsampled to common depth and the metric recomputed.
- **An expression-matched false-negative floor.** Pick control genes whose expression
  matches `GPRC5D`'s in malignant cells; their zero-fraction in those same cells is
  the technical false-negative rate the antigen call can't beat. This turns "GPRC5D
  is lowly expressed" from a hand-wave into an actual number.
- **Confidence intervals on every patient.** Bootstrap over each patient's malignant
  cells. Sample cell counts in this cohort vary about fifteen-fold, so a bare rank
  ordering claims a precision it doesn't have. A **minimum malignant-cell rule** is
  declared up front and excluded patients are **named in the output**, never quietly
  dropped.
  - **≥50 cells is a floor, not the answer (revised 2026-08-21).** At n=50 a single
    cell *is* 2%, so 1% / 2% / 3% escape are not separable and ordering patients across
    that range is noise with a decimal point on it. Derive the threshold from the
    smallest DN fraction the project intends to call meaningful — a 5% population gives
    ~2.5 expected DN cells at n=50, ~5 at n=100, ~10 at n=200 — so expect it to land
    nearer **100-200**. The procedure is unchanged: inspect the distribution, fix the
    threshold, *then* look at the results.
  - **Bootstrap at the right level (corrected 2026-08-21).** Cells within a patient
    aren't independent draws, and several patients contribute multiple samples
    (`27522_1`…`_6` and others); a flat cell-level bootstrap treats sample-level batch
    variation as biological spread and reports intervals that are too narrow. An earlier
    draft said **patient → sample → cell** for the per-patient CI — wrong level, because
    a CI *for patient A* conditions on patient A, so patient is fixed rather than
    random. Use **sample → cell within patient** for per-patient CIs, and
    **patient → sample → cell** only for cohort-level quantities (mean escape,
    regression coefficients, distributions). For single-sample patients this collapses
    to a cell bootstrap that cannot see sample-level variation at all, so their
    intervals are optimistic relative to multi-sample patients — stated with the
    results, not buried. The nesting was S1-gated; S1 landed 2026-08-24, so the real
    sample→patient map is now available and this is no longer provisional. Eight
    patients contribute more than one sample (`27522` six, and `47491`, `56203`,
    `58408`, `59114`, `60359`, `81012`, `83942` two each), which is where the
    sample-level term actually has something to estimate.

### Co-negativity enrichment — the key derived metric (added 2026-08-21)

The escape fraction alone can't distinguish two clinically different tumors. Per
patient, build the 2×2 of BCMA± × GPRC5D± over malignant cells and compare observed
double-negatives against the independence expectation
`E[DN] = P(BCMA⁻) × P(GPRC5D⁻)`, reporting the **co-escape enrichment ratio**
`observed / expected` with a Fisher's exact test and a permutation CI.

This splits three facts the single metric fuses together: how often each antigen is
individually absent, how many cells are DN, and whether the *same* cells are
disproportionately losing both.

**What it does not mean (corrected 2026-08-21).** An earlier draft said an enriched
patient is one "a second binder does not solve". That overstates it. Adding GPRC5D to
BCMA moves the uncovered fraction from `P(BCMA⁻)` to `P(BCMA⁻ ∩ GPRC5D⁻)`: at 30%
BCMA⁻ / 20% GPRC5D⁻ under independence, 30% → 6%; with enrichment pushing DN to 15%,
30% → 15%. Less than independence promised, but still halving the escape population.
Enrichment measures **how much of the pair's expected complementarity is eroded by
correlated loss** — not whether the second target is worth adding.

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
enrichment on data with no biological co-occurrence whatsoever — an artifact pointing
in exactly the direction the project hopes to find, which is the worst kind. So permute
**within depth strata** (or compute `E[DN]` from a per-cell independence model where
each cell's `P(BCMA⁻)` and `P(GPRC5D⁻)` are functions of its own depth), and **report
the unconditioned ratio alongside** — the gap between the two *is* the artifact,
quantified.

### The detection curve — and what it cannot deliver (corrected 2026-08-21)

An earlier draft proposed a "dropout-adjusted expected DN" as
`Σ_i P_i(BCMA⁻) · P_i(GPRC5D⁻)`. **That formula is circular** — multiplying the two
marginals assumes exactly the independence the co-escape test exists to interrogate, so
a tumor with genuinely correlated loss would be "corrected" toward the null it violates.

Two consequences. First, a simplification: the "dropout-adjusted DN" and the "expected
DN under depth-conditioned independence" are **the same number**, specified twice. It is
reported once, as the **depth-adjusted DN expectation under conditional independence** —
a technical baseline to compare against, never a corrected truth.

Second, and stated plainly: **no dropout-corrected DN point estimate is produced.**
Dropout is *bounded* here — threshold sensitivity band, expression-matched
false-negative floor, depth regression, downsampling — not corrected. Observed DN
remains the point estimate, reported as an interval. That is a stronger position than
shipping a number whose correction rests on an assumption the project is simultaneously
testing.

The detection curve itself is still built (detection probability vs. cell depth and gene
mean, fit on the expression-matched controls): it is what makes the depth-conditioned
null quantitative rather than rank-based, and what turns "GPRC5D is lowly expressed"
into a number. A genuinely dropout-corrected DN would need a latent-class model over the
four true states with per-cell detection probabilities, fit by EM over the observed 2×2 —
which estimates the joint without assuming independence. Real work, not on the critical
path, filed so it is neither reinvented casually nor mistaken for an oversight.

**Imputation and denoising (MAGIC, scVI, ALRA) are forbidden for positivity calls.**
They manufacture low-level expression by borrowing from neighboring cells, and whether
a transcript is genuinely absent is the entire question being asked. Model the
uncertainty; don't fill it in.

### From two antigens to a coverage matrix

`SLAMF7` and `FCRL5` get promoted from "backups" to a real deliverable. For every
pair and triple across {`TNFRSF17`, `GPRC5D`, `SLAMF7`, `FCRL5`, `CD38`, `SDC1`,
`ITGB7`}, compute how much of each patient's clone would be left uncovered. This
answers the question a target-strategy audience actually asks and the two-antigen
metric structurally cannot: *is BCMA+GPRC5D the best pair for this patient, or would
BCMA+FCRL5 cover more of their tumor?* Coverage gets traded off against normal-cell
expression from Stage 09 rather than maximized on its own — a target that covers the
whole tumor and also hits healthy tissue is not a better target.

**Blocking prerequisite — DISCHARGED 2026-08-24.** `patient_id` was provisional (a
naive rule gave 43 patients from 54 myeloma samples vs. the paper's 41/53), and the
risk was that one patient's cells split across duplicate entries, each producing its
own wrong, partial escape fraction. Supplementary Table S1 resolved it: **41 patients
over 53 in-cohort samples**, exactly. Two independent corrections — `25183` is
deposited but listed in no supplementary table (`in_paper_cohort == False`, retained
but excludable on purpose), and `83942`/`MMY83942` are one patient sampled under both
WashU protocols. `io.s1_patient_id` is the mapping and the parser asserts all three
counts, so a revised S1 fails loudly rather than quietly moving this denominator.

S1 also confirmed the `_N` suffixes are **serial disease-course timepoints**
(`27522_1` Primary → `_6` Relapse-3), which unblocks the longitudinal arm. It carries
**no cytogenetics**, so the t(4;14)/1q21 annotation this stage also wanted is still
unavailable from the deposit.

**A different confounder now takes its place as the thing to watch here.** Stage 04
found the deposit is pre-filtered differently per cohort: WashU 1 and 2 were cut at
10,000 UMIs before deposit, MMRF and the donors were not. That ceiling censors a band
enriched **20–70× for `GPRC5D`**, in 36 of the 54 myeloma samples, and it inflates
`frac_double_negative` for those cohorts — in the project's own direction of interest.
It is carried as a covariate, and this stage owes a **truncate-all-cohorts-at-10,000
sensitivity analysis** alongside the threshold band and depth regression.

---

## Stage 09 — Escape robustness (`notebooks/09_escape_robustness.ipynb`; env: `mm-core`)

New stage. Everything in it exists to answer one question: *how do you know your
escape fractions are real?*

**Matched bulk RNA-seq validation.** GSE223061 — the matched bulk data — was
downloaded at the very start of the project and then never used by the plan.
`raw/unpacked_bulk/` holds 29 usable bulk samples, of which **26 have an exact
single-cell match** — computed from the GEO metadata, not inherited. For those,
compare pseudobulk `TNFRSF17` and `GPRC5D` against bulk TPM, **cohort by cohort**:
MMRF bulk is CD138+ sorted and pairs with *malignant-cell* pseudobulk, while WashU
cohort 1 bulk is unsorted whole marrow and pairs with *whole-sample* pseudobulk.
Pooling the two would be a real error rather than a nicety — correlating malignant-PC
pseudobulk against unsorted bulk measures tumour burden, since the dilution by
non-plasma cells scales with how much tumour is present, and that would corrupt 10 of
the 26 comparisons in a direction correlated with the metric itself. Agreement means the antigen quantification is technically
credible — an orthogonal assay, generated independently, agrees with the per-cell
calls. Disagreement in one specific direction, **bulk-positive where single-cell
reads zero**, is direct quantified evidence of dropout, and feeds straight back into
Stage 08's false-negative floor. Either outcome is informative, which is what makes
this worth doing. Two of the bulk files are empty 114-byte stubs and three sample
IDs (`47499`, `59114_2`, `98433`) have no single-cell counterpart — both documented
in `CLAUDE.md`, both to be handled rather than silently absorbed.

**Scope correction, 2026-08-21 — bulk validates antigen abundance, not the escape
fraction.** Bulk averages every cell in the sample together, which destroys the joint
per-cell distribution the metric is built on. A tumor that is 50% BCMA⁺GPRC5D⁻ and 50%
BCMA⁻GPRC5D⁺ shows healthy bulk expression of both genes while containing zero
dual-positive cells. So the correct framing of every output here is **orthogonal
validation of antigen abundance and of whether the single-cell antigen-negative calls
are plausible** — never validation of `frac_double_negative`, which has no orthogonal
check available in this project.

**Normal plasma-cell antigen baseline — marrow expression context, not a safety axis
(scope corrected 2026-08-21).** Do *normal* plasma cells, from the healthy `BM*`/`ND_*`
marrow, express BCMA and GPRC5D? This is genuinely worth having: BCMA is broadly
expressed on normal plasma cells and B-lineage cells, and the malignant-versus-normal-PC
contrast is what makes a coverage number interpretable rather than absolute. It feeds
the coverage-matrix trade-off in Stage 08, because coverage alone is the wrong thing to
maximize. What it is *not* is a safety axis. GPRC5D's clinically decisive off-tumor site
is keratinized tissue — the nail, skin and taste toxicity seen with talquetamab — and a
bone marrow dataset cannot observe that at all; expression is also not toxicity. Keep
three things separate in the writeup: **tumor coverage**, **normal marrow expression**
(measured here), and **known extra-marrow liabilities** (cited from external evidence,
not measured). A real utility score of the form `coverage − λ · normal-tissue exposure`
needs a normal-tissue atlas (GTEx/HPA) and is a future extension, not a claim from this
data.

**Label-permutation null — moved to Stage 08 (corrected 2026-08-21).** Shuffling antigen
labels within a patient while preserving each antigen's marginal negative rate does not
produce "the metric with no signal in it" — the marginals are exactly what it holds
fixed. What it actually tests is whether BCMA-negativity and GPRC5D-negativity
*co-occur* beyond independence, which is a sharper and more interesting question than
the one it was written for. It therefore moves to Stage 08 as the co-negativity
enrichment test, with the depth-stratified null documented there.

---

## Stage 10 — Escape subclone + phenotype (`notebooks/10_escape_subclone_phenotype.ipynb`; env: `mm-core`)

New stage, and the one carrying the project's actual scientific payoff rather than
another robustness check.

**Is the double-negative population structured, or just scattered noise?** This is
the most important question the project can ask of its own metric. "3% of this
patient's malignant cells are double-negative" and "this patient carries a
pre-existing 3% resistant subclone" sound like the same statement and are not — only
the second one predicts that therapy will *select* for those cells, which is the
entire clinical premise of measuring baseline escape in the first place.

**Corrected 2026-08-21 — transcriptional clustering does not establish clonality.** The
earlier formulation ("cells scattered at random is the signature of dropout; cells
clustered together is the signature of a real subclone") was too strong in both
directions. A transcriptionally coherent group can come from cell cycle, stress,
interferon tone, metabolic state, sequencing depth or sample-prep batch just as easily
as from a genetic subclone — and a genuine genetic clone need not form a tidy
transcriptional island. The name `clonality-of-escape` prejudged precisely the question
the analysis exists to ask, so it is retired in favour of **DN coherence**, tested at
three escalating levels with the claim escalating alongside the evidence:

| Level | Question | Method | Licenses saying |
|---|---|---|---|
| **1 — enrichment** | Are DN cells non-randomly located in malignant transcriptional space? | kNN-neighborhood enrichment, Moran's I on the DN label, per-patient subclustering + Fisher, **depth-stratified within-patient permutation** | "non-randomly distributed" |
| **2 — transcriptional coherence** | Do DN cells share a reproducible program? | the program scores below — MYC, OXPHOS, stress, IFN, UPR, antigen presentation, γ-secretase | "an escape-associated **state**" |
| **3 — genomic coherence** | Do DN cells preferentially occupy a CNV-defined subclone? | `infercnvpy` substructure from Stage 07 | "an escape-associated **subclone**" |

**Only level 3 licenses the word "subclone."** Levels 1+2 alone are an escape-associated
*state* — still the thing that separates a structured 3% from a scattered 3%, and still
worth reporting, but not a claim about pre-existing genetic clones under selection.

**And a negative at level 3 is not evidence of absence.** Resolving CNV substructure
*within* a single patient's clone is far harder than separating tumor from normal, and
at ~2,044 median genes/cell it will often be underpowered. Level 3 reports **supported /
not evaluable**, with per-patient CNV resolution stated — never "no CNV subclone".
Treating an underpowered null as a negative would systematically understate exactly the
risk this project exists to measure.

The per-patient level attained is reported next to the escape fraction. This runs
on the per-patient un-integrated embedding from Stage 05, never the Harmony one — for
the reason given there. Note that the level-1 permutation carries the same
depth-stratification requirement as Stage 08's co-negativity test: shallow cells are
both likelier to read double-negative *and* to sit together in low-dimensional space,
so an unconditioned enrichment test will see depth structure and report it as biology.

**What else is different about the escape cells?** Pseudobulk differential expression
(`pydeseq2`/`decoupler`) of double-negative vs. dual-positive malignant cells, with
**patient as the unit of replication** — `sc-best-practices` is blunt that per-cell DE
tests treat cells as independent replicates and badly inflate the false discovery
rate. Then pathway and TF activity via `decoupler` (Hallmark, PROGENy, CollecTRI).

**One hypothesis is pre-registered: the γ-secretase axis** (`NCSTN`, `PSEN1`,
`APH1A`, `APH1B`, `PSENEN`). γ-secretase physically cleaves BCMA off the cell
surface, which is why γ-secretase-inhibitor plus BCMA CAR-T combinations are in
active clinical development. If escape cells turn out to be γ-secretase-high, that's
a directly actionable, mechanistically grounded finding rather than a descriptive
one. It's written down here, before looking, specifically so it stays a hypothesis
test instead of becoming a post-hoc story.

---

## Stage 11 — Cell-cell communication (`notebooks/11_cellchat_liana.ipynb`; env: `mm-communication`)

Runs LIANA+ using the CellChat-algorithm method specifically (LIANA+ natively
reimplements it), or optionally the full consensus rank-aggregate across
CellPhoneDB/CellChat/Connectome/NATMI/SingleCellSignalR — arguably stronger evidence
than any single method alone, since it's cross-validated across independent
approaches. The question is whether high-escape-risk patients also show weaker
NK/T-cell engagement signaling toward malignant plasma cells.

**The statistical design changed in the scope expansion.** The original plan split
patients into high/low tertiles and compared the two pools. That's
**pseudoreplication**: pooling cells across patients treats thousands of cells from
one person as thousands of independent observations, which inflates significance
essentially arbitrarily. The real sample size is the number of patients (~41), not
the number of cells (~181,000). So instead:

- Interaction scores are computed **per patient**, making the patient the unit of
  replication.
- `frac_double_negative` enters as a **continuous predictor** rather than being
  binned into tertiles — no arbitrary cutoff, and strictly more power than throwing
  away the middle third of the cohort.
- **The confounder gets controlled explicitly.** High-escape patients might simply
  have fewer T and NK cells to begin with, in which case "weaker immune signaling"
  is a composition artifact and the finding evaporates. Cell-type abundance from
  Stage 06 goes in as a covariate, and/or cells are downsampled to equal per-type
  counts per patient. This one is fatal to the stage's claim if ignored, so it isn't
  treated as optional polish.

---

## Stage 12 — Final decision packet (`notebooks/12_decision_packet.ipynb`; env: `mm-core`)

The final stage, assembling everything upstream of it. Assembles: the ranked escape-fraction table (annotated with disease stage and
cytogenetic risk from Supplementary Table S1, once resolved), a UMAP of malignant
cells colored by `coverage_class` faceted by patient, the LIANA+ differential
interaction plot, and a short written interpretation of which patients are poor
candidates for a BCMA/GPRC5D dual-target construct alone.

Changes from the scope expansion:

- **A caterpillar plot with confidence intervals replaces the ranked bar chart.** A
  bar chart asserts a precision this metric doesn't have; the CIs from Stage 08 are
  part of the result, not a footnote to it.
- **Risk tiers replace the rank ordering (2026-08-21).** The caterpillar plot fixed
  the chart but not the deliverable — "#1 patient 123, #2 patient 456" is false
  precision when the intervals overlap. Patients get **robust high escape** /
  **uncertain** / **robust low escape**, and the ranking-stability check from Stage 08
  becomes what *earns* a tier rather than being the output itself.
- **The multi-antigen coverage matrix, the co-escape enrichment ratio, the DN-coherence
  level, and the bulk-validation correlation** all join the packet as separate columns
  — respectively: is there a better target pair for this patient, are the same cells
  losing both antigens, is that population structured, and does an orthogonal assay
  agree with the single-cell antigen levels. They are kept as separate columns rather
  than folded into one score, because a patient can be low on the escape fraction and
  high on co-escape enrichment, and that patient is more interesting than a bare
  ranking would ever show.
- **The bias-direction table** (ambient, dropout, malignant-call error,
  mRNA-vs-protein, each with its sign on the metric) is included rather than left
  implicit.
- **The mRNA-versus-protein limitation is stated explicitly and mechanistically.**
  CAR-T binds surface protein; this analysis measures transcript. BCMA is actively
  shed from the surface by γ-secretase, and GPRC5D transcript correlates imperfectly
  with surface density. This is the first question a target-strategy audience will
  ask, so it gets answered in the deliverable rather than waited for. If a published
  CITE-seq or flow calibration exists for these two antigens, quantify against it;
  otherwise state it plainly.
- **Decision rules are declared in advance** — which escape-fraction threshold makes
  a patient a poor dual-target candidate is fixed before the tiers are looked at,
  not fitted to them afterwards.

---

## Phase 2 — External validation on GSE117156 (not started; sequenced strictly after Phase 1 completes)

Independently re-runs the same pipeline shape (QC -> malignant calling -> antigen
scoring -> escape fraction) against GSE117156 (Ledergor et al. 2018, *Nat Med* —
51,840 cells, 11 healthy controls + 29 MM patients spanning asymptomatic disease
through post-treatment MRD) as a second, independent cohort, testing whether the
core finding replicates beyond this one dataset and technology.

Full reasoning, acquisition method, and the explicit no-merge constraint (MARS-seq
vs. 10x — a platform difference, not a correctable batch effect) are documented in
`CLAUDE.md`'s Phase 2 section, not duplicated here.
