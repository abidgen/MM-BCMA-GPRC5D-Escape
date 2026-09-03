# Decisions archive — superseded positions and how the plan evolved

Split out of the main project document on 2026-08-24. That document states the
**current** position only, one line per decision. This file keeps what each position
replaced and why, so a reversal is auditable and a settled question is not accidentally
re-opened from an older draft.

**Nothing in this file is an instruction.** Every position recorded here as "the old
reading" is wrong or obsolete. If this file and the main project document disagree, that
document is current.

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

## Stage 06 v1 — executed as declared, not accepted (2026-08-25)

Preserved at `results/06_annotation_v1/`. Recorded here because the correction it
prompted is exactly the kind that must not look, later, like threshold-shopping.

**What v1 did right.** It executed the pre-declared rules without retuning anything.
Four classes cleared their bars on CellTypist (PlasmaCell F1 0.984, Bcell 1.000,
Myeloid 0.982, HSPC 1.000). Three did not, and took the pre-declared fallback to manual.

**What it got wrong.** The fallback produced assignments that are not credible for bone
marrow:

    NK        33,556      Tcell   19,133      Erythroid  35,855
    (CellTypist, for comparison:  Tcell 60,896   NK 11,424   Erythroid 16,224)

Erythroid was the runner-up class in **18 of 30 clusters**, and in the final manual
labelling 47% of "Erythroid" cells and 63% of "NK" cells were CD3D-positive.

**Diagnosis — two definitional defects, no threshold defect.**

1. *Ambient haemoglobin.* The manual Erythroid panel was `("HBB", "GYPA")`. Stage 04 had
   already established haemoglobin as the dominant ambient species in this marrow
   (~32% of counts in the flagged top-20 decile), and direct measurement confirmed HBB
   detected in **61-85% of cells of every class** while GYPA sat at 1-3% outside
   erythroid cells. A two-gene panel half-driven by ambient will attach itself to
   everything, which is precisely the 18-of-30 pattern.
2. *No exclusivity check.* `NKG7`/`GNLY` are a cytotoxic-granule program shared by NK
   cells **and** cytotoxic T cells. The validation framework asked only a
   precision-like question — do cells labelled NK express NK markers? — which those
   cells pass honestly, since they really are cytotoxic. Marker coverage for manual NK
   was **1.00**. Nothing in the framework could ask whether they *also* carried strong
   T-lineage evidence, so a class swallowing another lineage was structurally invisible.

**A third, structural flaw the first two exposed.** Manual was simultaneously the
concordance reference and the fallback. When manual is wrong for a class, both automated
methods "fail" by disagreeing with a bad reference, and the class then falls back to that
same bad reference. The fallback must clear the same vetoes as any other method.

**What was NOT changed.** `CONCORDANCE_THRESHOLDS` (0.95 / 0.90 / 0.85) and
`MARKER_COVERAGE_MIN = 0.30` are unchanged, and remain unchanged in v2. The correction is
to the **reference specification and the validation framework**, not to the acceptance
bars. Revising a definition that measurement showed to be ambient-driven is not the same
act as lowering a bar until the desired labels appear, and the distinction is the whole
reason this entry exists.

**Justification wording corrected mid-run (2026-08-25), value unchanged.** The original
rationale for `CONTRADICTION_MAX_RATE` claimed 25% sat "roughly an order of magnitude
above what contamination can explain". That is quantitatively loose — a few percent of
doublets plus a few percent of ambient is not an order of magnitude below 25%. The
threshold **value was not changed**; only the prose defending it, which now reads that
25% is deliberately *permissive* relative to expected residual technical contamination,
so the veto triggers only when contradictory evidence affects a substantial fraction of
the assigned class. Recorded here because correcting a justification after a run has
started is exactly the kind of edit that must be visible rather than quiet.

**What v2 changes** (all fixed in `config.py` and the main project document before v2
was run):
revised `MARKER_PANEL` for Erythroid, Tcell and NK from lineage biology; new
`LINEAGE_PROGRAMS`, `CONTRADICTION_PAIRS`, `CONTRADICTION_MIN_GENES = 2` and
`CONTRADICTION_MAX_RATE = 0.25`; and the fallback subjected to both vetoes.

**One v1 result that must not be misread.** SingleR/Novershtern's PlasmaCell F1 of 0.000
is **not evidence against the plasma-cell annotation**. `NovershternHematopoieticData`
contains no plasma-cell label at any level — its B-lineage stops at "Mature B cells class
switched" — so the class is *not evaluable* against that reference rather than failed.

---

## Stage 06 v2 — the exclusivity veto worked; the manual classifier was still broken (2026-08-25)

Preserved at `results/06_annotation_v2/`. v2 fixed what v1 got wrong and exposed a third,
separate defect underneath it.

**What v2 fixed.** Erythroid collapsed from 35,855 to 16,224 cells once the panel stopped
being half-driven by ambient haemoglobin, and manual/CellTypist agreement went F1
0.623 -> 1.000. Manual Tcell coverage went 0.71 -> 0.987. PlasmaCell, Bcell, Myeloid and
HSPC re-derived byte-identically.

**What v2 caught.** The new lineage-exclusivity veto fired on NK: every method exceeded
`CONTRADICTION_MAX_RATE`, the fallback included, so NK was reported **unresolved** and
33,556 cells became `Ambiguous` rather than being shipped as NK. The six other classes
sat at 0.06-0.13 contradiction, which is the technical floor — NK at 0.56-0.83 was 4-7x
that, so the signal was real and not ambient TCR.

**The checkpoint that changed the diagnosis.** A hypothesis that all the NK contradiction
was localised to Leiden 23 (and that the veto was therefore over-firing at class level)
was **refuted**. Manual-NK spanned three clusters:

    leiden 3   11,819 cells   T-contradiction 0.958   CellTypist: Tcell   SingleR: Tcell
    leiden 12  10,313 cells   T-contradiction 0.985   CellTypist: Tcell   SingleR: Tcell
    leiden 23  11,424 cells   T-contradiction 0.528   CellTypist: NK      SingleR: NK

Cluster 23 carried only **21.9%** of the contradictions. Clusters 3 and 12 — 78.1% —
are T cells by every automated method, with CD3D at 86.7%/90.6% and KLRF1 at 8.1%/13.2%.
The veto was right; the manual labels feeding it were wrong.

**The third defect: cross-panel `score_genes` argmax is not calibrated.** v1 and v2 both
assigned clusters by argmax over module scores from different panels. `score_genes`
subtracts, per gene, the mean of a control set drawn from that gene's own expression bin,
so each panel carries a baseline offset set by where its genes sit inside their bins.
Measured here:

    T panel  : mean(gene - its bin's control mean) = -0.2364   (all 6 genes in bin 24)
    NK panel : mean(gene - its bin's control mean) = -0.0328   (spread across bins 18-24)
    systematic offset favouring NK                 =  0.2036

Clusters 3 and 12 were called NK over T by **0.068 and 0.196 — both smaller than that
offset**, which is therefore sufficient on its own to produce and reverse the calls.
The mechanism is not that NK's controls are lower in absolute terms; it is that the top
expression bin is wide, its control mean (0.813) exceeds *every* T gene, and the entire
T panel is penalised at once while the NK panel escapes by spanning four bins. Cluster
23, a genuine NK call, had a margin of 2.03 — an order of magnitude clear of the offset.

**v3 replaces the manual classifier, and nothing else.** Module-score magnitude no longer
decides identity. Clusters are adjudicated on **detection fractions**, which share a
common [0, 1] scale with no per-panel normalisation, under the rule
*positive lineage evidence + specificity/exclusion evidence*:
`MANUAL_MARKER_DETECT_MIN = 0.25`, `MANUAL_POSITIVE_MIN = 0.5`,
`MANUAL_DECISION_MARGIN = 0.15`, exclusion reusing the existing
`CONTRADICTION_PAIRS`/`CONTRADICTION_MAX_RATE`. `score_genes` survives as a descriptive
within-program quantity only.

**Unchanged across all three versions:** the v2 marker panels, the F1 bars
(0.95/0.90/0.85), `MARKER_COVERAGE_MIN = 0.30`, `CONTRADICTION_MIN_GENES = 2` and
`CONTRADICTION_MAX_RATE = 0.25`. Each revision has been to a *definition* or a *method*
that measurement showed to be broken; none has been to an acceptance bar.

**Still open, deliberately not addressed in v3.** Leiden 23 remains a biological
resolution question — `KLRD1` 89.6% / `KLRF1` 75.0% / `FCGR3A` 65.8% is NK-like, but
`TRDC` 66.4% and `TRBC2` 46.6% indicate a γδ-T or NKT-like component, and its 0.528
contradiction rate would fail the veto on its own. And whether clean NK populations
outside cluster 23 should survive a class-level veto triggered by one cluster is a
separate question, held for review.

---

## Stage 06 — TRBC1/2 demoted from T identity to T context (2026-08-25)

**Isolated `TRBC1`/`TRBC2` expression is insufficient evidence of T-lineage
commitment, because it frequently occurs without coordinated `CD3`/`TRAC` expression in
cells with strong NK-lineage evidence.** Validated on Leiden 23: among 5,788 cells with
strong NK evidence and no CD3/TRAC, 74.8% were TRBC-positive; among TRBC-positive mixed
cells only 12.9% carried both CD3 and TRAC. Candidate mechanisms — germline/unrearranged
TRB transcription, ambient spillover, residual multiplets — are **hypotheses**; the
diagnostic did not distinguish them (TRBC intensity exceeds a one-molecule ambient
profile, TRBC positivity does *not* rise with sample T-cell abundance, r = -0.163,
p = 0.3, and the doublet-score gap is too small to explain the compartment). Do not
state any mechanism as established.

**This is not a threshold change.** Of the 4,080 cells that move mixed -> NK under the
revision, only 30.2% carry a single TRBC UMI and 100% remain TRBC-positive at a median
of 2 UMI. The transition is driven by the requirement for coordination with
lineage-defining machinery, not by thresholding away weak signal. `config.T_IDENTITY_ANCHORS`
(CD3D/E/G, TRAC) and `config.T_CONTEXT` (TRBC1/2) are predeclared; `MARKER_PANEL`,
`LINEAGE_PROGRAMS` and every numeric parameter are unchanged, and the revision is so far
**validated on cluster 23 only** — global C2d is untouched.

Cluster-23 validation: NK 4,871 -> 8,951 (78.4%, 49 patients, top-patient 8.6%);
NKT_like_mixed 5,460 -> T_NK_mixed 1,380 (now 100% any-CD3, 55.6% CD3+TRAC);
unresolved 524 -> 827; T_ab 416 -> 170; T_gd 153 -> 96. No reverse transitions.
gamma-delta conclusion unchanged: no robust cohort-wide population.
Artifacts: `results/06_annotation/cluster23_local/trbc_context_revision/`.

---

## The design reviews, verbatim

These describe how the plan changed. The resulting rules are restated at each stage in
the main project document; they are kept here because the *reasoning* is what stops a
settled question being re-opened. Two figures below were corrected after the fact and
are annotated inline.

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

---

## Risk tiers are PROVISIONAL (recorded 2026-08-25)

**The Stage-08/09b risk tiers are PROVISIONAL, and the word "final" has been retired from
them (2026-08-25).** They live at `results/08_dual_antigen_escape/risk_tier_provisional/`
and are **measurement-robust provisional tiers under the frozen Stage-08/09b rule** — a
`robust-high` label means the observed DN fraction survived the denominator, depth,
repeated-sample, null-scheme and threshold sensitivity analyses, and **nothing more**.
> **Stage 10 sits between provisional measurement-robust tiering and any final biological
> risk classification.** No provisional tier may be cited as a final biological escape
> classification until stage 10's coherence states are frozen.
Provisional outcome: 4 `robust-high` (`MMRF_1267`, `MMY18273`, `MMY74196`, `MMY98423`),
28 `uncertain`, **0 `robust-low`**, all 32 `THRESHOLD_ROBUST` across
`TAU_HIGH` ∈ {0.20, 0.25, 0.33}. The relabelling changed **no** patient membership, no
threshold, no measurement and not the algorithm — provenance and terminology only.

---

## Stage 10 frozen; the two axes stay separate (2026-08-25)

**Stage 10 is RUN, ACCEPTED and FROZEN exactly as generated (2026-08-25).**
4 `DN_COHERENCE_SUPPORTED` (`MMRF_1720`, `MMRF_2038`, `MMY34339`, `MMY80649`),
23 `DN_COHERENCE_NOT_SUPPORTED`, 5 `DN_COHERENCE_NOT_EVALUABLE`, and
**`CNV_SUBCLONE_NOT_EVALUABLE` for all 32** — the failed stage-07 CNV method was not
reused, and `CNV_SUBCLONE_NOT_SUPPORTED` is never emitted because an underpowered null is
not a negative. Coherence licenses **"escape-associated transcriptional state"**, never
"subclone". Full results: `results/10_dn_coherence/`.

### The two frozen axes stay separate — no composite score, ever

> **Elevated DN burden and coherent DN-associated transcriptional organization are
> separable properties in this cohort.**

| axis | states |
|---|---|
| **measurement** (stage 08/09b) | `measurement_robust_high` · `measurement_uncertain` · `measurement_robust_low` |
| **biological coherence** (stage 10) | `DN_COHERENCE_SUPPORTED` · `DN_COHERENCE_NOT_SUPPORTED` · `DN_COHERENCE_NOT_EVALUABLE` |
| **genetic** | `CNV_SUBCLONE_NOT_EVALUABLE` (all patients) |

**Do not invent a post-hoc rule collapsing these into one scalar escape-risk score.** The
observed cross-tab **is itself a result and must not be "fixed"**: 4 measurement-robust-high
patients, **all** `DN_COHERENCE_NOT_SUPPORTED`; 4 `DN_COHERENCE_SUPPORTED` patients, **all**
measurement-uncertain; **zero patients with convergent high-measurement + coherence-supported
evidence.** The final synthesis displays measurement robustness, DN coherence and
exploratory immune context as **three distinct evidence layers**.

### `MMRF_1640` — the project's primary methodological example

Moran's I **0.47**, unconditioned permutation p **0.001**, depth-stratified p **0.499**.
Use it to illustrate the project-level lesson, not as a patient anecdote:

> **A visually and statistically strong DN-coherence signal can be generated by sequencing
> depth alone, and conditioning on depth can eliminate the apparent biology.**

Its force is that **the same technical confound reproduces across two different
statistical frameworks** — stage 08's count-based co-negativity enrichment (unconditioned
median 1.052 → depth-conditioned 1.009; `MMRF_1505` 4.61 → 1.42) and stage 10's
neighbourhood coherence. Same artifact, different mathematics, same direction — toward the
project's own hypothesis. **No stage-10 result is altered because of this example.**

---

## Stage 10 COMPLETE and FROZEN — three-level hierarchy (2026-08-25)

**Stage 10 is COMPLETE and FROZEN (2026-08-25), through the full three-level evidence
hierarchy.** Results: `results/10_dn_coherence/`.

### The three levels are independent axes and are never combined

| level | states | licenses |
|---|---|---|
| **1 — DN structure** | `DN_STRUCTURE_SUPPORTED` **4** · `NOT_SUPPORTED` **23** · `NOT_EVALUABLE` **5** | *non-random DN organization* — **nothing more** |
| **2 — DN phenotype** | `DN_STATE_SUPPORTED` **26** · `NOT_SUPPORTED` **1** (`MMY67868`) · `NOT_EVALUABLE` **5** | *escape-associated transcriptional state* |
| **3 — genomic** | `CNV_SUBCLONE_NOT_EVALUABLE` **32 (all)** | nothing |

**Level 1 alone is NOT an escape-associated transcriptional state.**
**`CNV_SUBCLONE_NOT_SUPPORTED` is never emitted or implied, and no patient may be called a
genetic escape subclone.**

> **The predeclared per-patient Level-2 criterion proved weakly discriminative, with 26/27
> evaluable patients satisfying it. `DN_STATE_SUPPORTED` therefore indicates compatibility
> with the cohort-level DN-associated program and should not be interpreted as strong
> patient-specific evidence of a distinct escape state.** The rule was **not** retuned
> after seeing that result. **The scientifically useful Level-2 result is the cohort-level
> transcriptional phenotype.**

### The frozen Level-2 biological conclusion

> DN cells exhibit a reproducible **cohort-level** transcriptional shift characterised by
> reduced plasma-cell secretory/differentiation, antigen-presentation, OXPHOS and
> interferon programs. **This supports a biological state associated with observed DN
> status, but does not establish antigen-specific escape.** Patient-level Level-2
> classification is weakly discriminative, and independent genomic subclone evidence is not
> evaluable.

**Keep "biology associated with DN status" and "antigen-specific escape mechanism"
distinct. The latter is not claimed.**

Reproducible DN-**lower** programs (BH < 0.10 under **both** denominators, consistent
sign): **antigen presentation · OXPHOS · interferon**. Not reproducible: **MYC · stress ·
UPR · γ-secretase**. Do not search additional programs.

**Pseudobulk DE — 190 genes significant under both denominators**, paired over patients on
depth-matched cells. Many DN-lower genes are ER/secretory-machinery and mature plasma-cell
identity genes (`SPCS1`, `SPCS2`, `SEC61B`, `UBE2J1`, `TMBIM6`, `MZB1`, `B2M`).
**The DN phenotype is compatible with a less secretory / less differentiated plasma-cell
state.** Never write that DN cells "specifically evolved an antigen-escape program".
**The present data cannot separate a broad plasma-cell-state shift from an antigen-specific
escape mechanism**, and BCMA/GPRC5D are themselves secretory-pathway-dependent surface
proteins, which is exactly why the two are confounded here.

### The pre-registered γ-secretase result is a clean NEGATIVE

Frozen five-gene hypothesis, **no gene added after seeing results**: `NCSTN`, `PSEN1`,
`APH1A`, `APH1B`, `PSENEN`. **Not supported** under the frozen both-denominators rule
(depth-matched BH 0.387 primary / 0.070 sensitivity). **The direction is opposite to the
pre-registered γ-secretase-high prediction.** This stays visible in the stage-10 summary.

### The depth lesson — now the THIRD independent instance

Depth conditioning altered the biological interpretation again at Level 2: matching moved
the DN/comparator depth ratio from **~0.470 to ~0.992**; **γ-secretase looked strongly
negative before matching and not after**; **MYC changed sign**; OXPHOS was the strongest
depth-tracking program tested (ρ 0.274).

> **Three independent places where naive/unconditioned analysis was materially misleading:
> (1) stage 08 co-negativity enrichment, (2) stage 10 Level-1 DN structure, (3) stage 10
> Level-2 transcriptional phenotype.**

`MMRF_1640` remains the project-level example — **Moran's I 0.47, unconditioned p 0.001,
depth-stratified p 0.499**. This is part of the frozen methodological result, **not** a
prompt for another correction or rerun.

**decoupler** (Hallmark/PROGENy/CollecTRI, full tested space retained): E2F targets / G2M /
MAPK up, p53 / OXPHOS / JAK-STAT down. **Limitation, retained: decoupler ULM significance
is based on the fitted contrast vector and is not an independent patient-replicated
hypothesis test.** These p-values may not create patient states or upgrade Level-2 evidence.

**TC-like subtype** stays per-patient, a transcriptional proxy, **descriptive only**, never
a translocation diagnosis, and never an input to structure/state calls. No association
testing.

### Cross-axis — four separate axes, no combined score, ever

Measurement **4 / 28 / 0** · Level-1 **4 / 23 / 5** · Level-2 **26 / 1 / 5** ·
Level-3 **32 not evaluable**.

> **No patient is simultaneously measurement-robust-high, Level-1 structure-supported and
> Level-2 state-supported.**
>
> **The apparent measurement-high × Level-2 overlap of four patients is not strongly
> discriminative, because Level-2 support occurs in 26 of 27 evaluable patients.**

**Do not invent a combination rule or a scalar risk score.**
