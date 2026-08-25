# Decisions archive — superseded positions and how the plan evolved

Split out of `CLAUDE.md` on 2026-08-24. `CLAUDE.md` states the **current** position
only, one line per decision. This file keeps what each position replaced and why, so a
reversal is auditable and a settled question is not accidentally re-opened from an
older draft.

**Nothing in this file is an instruction.** Every position recorded here as "the old
reading" is wrong or obsolete. If this file and `CLAUDE.md` disagree, `CLAUDE.md` is
current.

---

## Superseded positions, at a glance

| # | The old position | Why it was wrong | Current position | Corrected |
|---|---|---|---|---|
| 1 | `56203_1` is a third, 22184-gene reference build missing `TNFRSF17`; exclude it | Misdiagnosis. The matrix header declares a normal 33694 build; `genes.tsv` is a truncated **write** (22185 rows, no trailing newline — `wc -l` undercounted by one), and the written rows are a strict prefix of the canonical list. `TNFRSF17` and `IGLC1/2/3` were past the cut, never absent from a reference | Repaired and retained. `io.read_sample` substitutes the canonical column behind a prefix assertion that raises rather than guessing. `config.EXCLUDED_SAMPLES` is empty; the mechanism is `config.TRUNCATED_GENE_FILES`. Recovers 1,837 cells and a second timepoint for patient 56203 | 2026-08-24 |
| 2 | Harmonize the gene space by canonicalizing HGNC symbols, with a four-gene alias map (`WHSC1`→`NSD2` etc.) | It addressed 4 of the **11,140** symbols that drift between the two builds, and a symbol join can silently pair the *wrong* gene — `TBCE` is `ENSG00000285053` in one build and `ENSG00000116957` in the other | Join on reconstructed Ensembl ID. 32,991 genes vs 22,164 on raw symbols (+10,827). The alias map is demoted to a **regression assertion**, never the mechanism | 2026-08-21 |
| 3 | Drop the ~52 `make.unique`-ambiguous symbols as an interim | Those genes resolve correctly under the ID join; dropping them discards real data | **Must not be implemented.** The 9 symbols still colliding *within* the ID intersection get `SYMBOL__ENSG…` suffixes, never `var_names_make_unique()` | 2026-08-21 |
| 4 | The cohort is 47 patients / 57 disease samples | It counted the four `ND_*` samples as disease. `!Sample_source_name_ch1` reads `Donor BMMC` for all four `ND_*` and all four `BM*`; the other 54 read `diagnosis: Multiple myeloma (MM)` | 43 naive patients over 54 myeloma samples; **41 patients / 53 in-cohort samples** under the S1 mapping, reproducing the paper exactly | 2026-08-24 |
| 5 | The depositors' stated Seurat filter "was not applied to what is deposited" | Reasoned from a **cohort-wide average** UMI count that pooled MMRF with WashU and hid per-cohort structure. Read per cohort the boundaries are unmistakable — a max of exactly 9,999 is a cutoff, not a distribution | The deposit **is** pre-filtered, differently per cohort. WU1/WU2 were cut at 10,000 UMIs; MMRF and Donor were not. A first-order confounder for stage 08, carried as a covariate with a truncate-all-at-10k sensitivity analysis owed | 2026-08-24 |
| 6 | No raw data exists for this series | Too strong. The series design sends raw data to dbGaP for patient privacy — `ND_*` at `phs000159`, MMRF bulk at `phs000748` | Raw data exists **under controlled access**. The practical conclusion is unchanged (no unfiltered matrices without a DAC application, so SoupX/DecontX stays unavailable), but "doesn't exist" and "exists behind a data-access committee" are different claims | 2026-08-24 |
| 7 | ~28 bulk samples match scRNA; 18 + 12 = 30 usable bulk | Wrong twice over — the WashU count was 13, not 12, and the total did not subtract the two empty 114-byte stubs | **29 usable bulk** ((18−2)+13); **26** with an exact scRNA match, computed from the GEO titles | 2026-08-21 / 2026-08-24 |
| 8 | The `_N` suffixes are probably fractions, sorts or replicates; bulk/scRNA suffix misalignment argues against timepoints | The misalignment argument rested on `59114` alone, and reflects incomplete bulk coverage of a shared index. `37692_2` and `57075_3` are **lone** samples with non-`_1` suffixes — a fraction label starts at 1 for a single sample, a serial event index does not | **Serial disease-course timepoints**, settled outright by S1 sheet 2 (`27522_1` Primary → `_6` Relapse-3). The longitudinal arm is real, not speculative | 2026-08-24 |
| 9 | `clonality-of-escape`: random scatter is the signature of dropout, spatial clustering is the signature of a real subclone | Too strong in both directions. Transcriptional coherence has many non-genetic causes (cell cycle, stress, IFN, depth, batch), and a genuine genetic clone need not form a tidy transcriptional island. The name prejudged the question | Retired for the three-level **DN-coherence hierarchy**. Only level 3 (CNV support) licenses the word *subclone*; levels 1+2 report an escape-associated **state**. Level 3 reports *supported / not evaluable*, never "no CNV subclone" | 2026-08-21 |
| 10 | An enriched co-escape patient "is the one dual targeting doesn't help"; co-escape "determines whether a second binder helps at all" | The arithmetic contradicts it. Adding GPRC5D moves the uncovered fraction from `P(BCMA⁻)` to `P(BCMA⁻ ∩ GPRC5D⁻)`: 30% → 6% under independence, 30% → 15% under strong correlated loss — still halving the escape population | Enrichment measures **how much of the pair's expected complementarity is eroded by correlated loss**, not whether the second target is worth adding. Incremental coverage gain is reported next to it as a separate column | 2026-08-21 |
| 11 | Report a dropout-adjusted DN estimate as `Σ_i P_i(BCMA⁻) · P_i(GPRC5D⁻)` | **Circular.** Multiplying the two marginals assumes exactly the independence the co-escape test exists to interrogate, so a tumour with genuinely correlated loss would be "corrected" toward the null it violates | No dropout-corrected DN point estimate is produced or claimed. The same computation is reported once, as the **depth-adjusted DN expectation under conditional independence** — a technical baseline, never a corrected truth. Dropout is *bounded*, not corrected. A joint latent-class/EM model is the honest correction and is deliberately deferred | 2026-08-21 |
| 12 | Bulk RNA-seq validates the dual-negative fraction | Bulk destroys the joint single-cell distribution. A tumour that is 50% BCMA⁺GPRC5D⁻ plus 50% BCMA⁻GPRC5D⁺ shows healthy bulk expression of *both* genes while containing zero dual-positive cells | Bulk is **orthogonal validation of antigen abundance and the plausibility of antigen-negative calls**, never of the escape fraction. The joint distribution over cells has no orthogonal check in this project | 2026-08-21 |
| 13 | Normal-BM controls give a safety axis for GPRC5D | GPRC5D's clinically decisive off-tumour site is keratinized tissue (the nail/skin/taste toxicity seen with talquetamab), which a bone marrow dataset cannot observe at all. Expression is also not toxicity | They give **normal *marrow* expression context**. Tumour coverage, normal marrow expression, and known extra-marrow liabilities stay three separate things | 2026-08-21 |
| 14 | The "label-permutation null" tests what the metric looks like under no signal (stage 09) | Permuting antigen labels within patient holds each marginal negative rate fixed — it tests **independence** between BCMA-negativity and GPRC5D-negativity, not absence of signal | A better question than the one it was written for. Moved to stage 08 as the **co-escape enrichment test**, with a depth-stratified null | 2026-08-21 |
| 15 | Harmony fails on plasma cells because "the populations genuinely differ" between cohorts | Implied biological divergence between cohorts and claimed more than the data supports | The separation is a **non-recoverable sampling/censoring asymmetry**: WashU's observed plasma-cell distribution is missing its high-RNA portion, so no one-to-one population correspondence remains for any method to recover | 2026-08-24 |
| 16 | There is a 2-3× v2-vs-v3 chemistry effect on depth | This cohort does not show one. v2 vs all-v3 is 1023 vs 1408 genes/cell = **1.38×**, and the sample distributions overlap (v2 max 1602 > v3 min 793) | The axis that separates is **cohort** (MMRF ≈ 1.9× the others), of which chemistry is one component alongside site and protocol. Carry `cohort` **and** `chemistry` as covariates; `n_genes_ref` is not a proxy for either | 2026-08-24 |
| 17 | The cohort's median cell has ~2,044 detected genes | An **R-build** figure (61 samples, 181,336 cells, fixed thresholds) carried into the Python rebuild without re-deriving it, then propagated to five separate arguments | Recomputed from `results/05_integration/integrated.h5ad`: **1,162** genes/cell over 172,940 post-QC cells; **1,521** in the plasma-like compartment. The cohort is about half as deep as the project had been arguing from, which *strengthens* every argument that depended on it. See `docs/stage-results.md` | 2026-08-24 |
| 18 | Split patients into high/low escape tertiles and compare (stage 11) | **Pseudoreplication** — pools cells across patients into two groups, treating thousands of cells from one patient as independent observations | Patient is the unit of replication (n ≈ 41), `frac_double_negative` is a **continuous** predictor, and T/NK abundance is controlled as a confounder | 2026-08-20 |
| 19 | Resample patient → sample → cell for the per-patient CI | Wrong level. A CI *for patient A* is conditional on patient A, so patient is fixed, not random | **sample → cell** within patient for a per-patient CI; **patient → sample → cell** for cohort-level quantities | 2026-08-21 |
| 20 | Rank patients 1..N by escape fraction | False precision when the CIs overlap heavily; ordinal positions are not stable quantities | **Risk tiers** (robust-high / uncertain / robust-low), with co-escape enrichment and DN coherence as independent columns. Ranking stability across thresholds stays as the robustness *diagnostic* that earns a "robust" label — it is not the deliverable | 2026-08-21 |
| 21 | `envs/env-scvi.yml` will be created if scVI is ever considered | Superseded by the stage-05b design, which compares **every** integration method on one footing | `envs/env-integration.yml`, named for what it holds rather than one member: harmonypy, scvi-tools, scanorama, bbknn, scib-metrics, celltypist | 2026-08-24 |
| 22 | Notebooks are thin orchestration over `src/` | The goal was never thinness — it was avoiding duplicated logic and keeping review on `.py` diffs | Notebooks carry the analysis; `src/mm_escape/` carries what is **reusable, testable, or fiddly**. The test is reuse and testability, not line count | 2026-08-20 |

---

## The design reviews, verbatim

These describe how the plan changed. The resulting rules are restated at each stage in
`CLAUDE.md`; they are kept here because the *reasoning* is what stops a settled
question being re-opened. Two figures below were corrected after the fact and are
annotated inline.

### First design review — scope expansion, 2026-08-20

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

### Second design review — 2026-08-21

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
  *within* a single patient's clone at this depth (~1,521 median genes/cell in the
  plasma compartment; the review said ~2,044, an R-build figure — see row 17 above). Level 3 of the coherence
  hierarchy therefore reports **supported / not evaluable**, never "no CNV subclone",
  and per-patient CNV resolution is stated rather than assumed.



---

## The Supplementary Table S1 policy, discharged 2026-08-24

S1 is in the repo and parsed. The policy below is **dead** and is kept only because the
"label provisional output as provisional" discipline still applies to anything computed
before it landed.

**What S1 closed:** the patient mapping (41/53), the `_N` suffix meaning (serial
timepoints), disease stage for WashU cohort 1, and ISS / treatment / TTPD per patient.

**What it did NOT close:** per-sample disease stage for MMRF and WashU 2 (absent, not
imputed), **cytogenetics** — S1 carries no t(4;14), 1q21 or other karyotype column, so
stage 10 has nothing in this deposit to validate its TC proxy against and it stays a
proxy — and the three bulk/sc ID mismatches, of which `47499`/`98433` are now explained
as bulk-only patients that S1 does list.

### The original policy, verbatim

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
