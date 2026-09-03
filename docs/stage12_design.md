# Stage 12 — final synthesis / decision packet: DESIGN

**Status: DESIGN ONLY. Stage 12 has not been executed and no notebook exists.**
**Design date:** 2026-08-27
**Design pinned to:** commit `5bbecbc`, tag lineage `pre-stage12-audit` → `pre-stage12-remediated`
**Authority rule:** where historical planning text disagrees with the recovered production
code under `production/`, **the recovered implementation is authoritative.**

Nothing in this pass modified a frozen result, threshold, tier, state or artifact. Every
number below was re-derived from the frozen artifacts while writing this document, and all
393 manifest hashes verify.

---

## Executive objective

> **Stage 12 synthesizes the already-frozen evidence axes into a scientifically defensible
> patient-level and cohort-level decision packet, while preserving uncertainty,
> non-evaluability, and disagreement between evidence types.**

Stage 12 is **not a discovery stage.** It runs no new statistical test, fits no model,
creates no threshold, and produces no score. Its entire job is to lay six frozen evidence
axes side by side, per patient and per cohort, in a form that makes their **disagreement**
legible — and to state precisely which claims that evidence does and does not license.

### The Stage-12 scientific question

> Given the frozen measurement, structural, phenotypic, genomic, immune and multi-antigen
> evidence for each patient, **what conclusions can be supported about baseline
> dual-antigen vulnerability, and how certain are those conclusions?**

Three questions Stage 12 explicitly does **not** ask, because each exceeds the evidence:

| forbidden framing | why it exceeds the evidence |
|---|---|
| "Who is high risk?" | risk requires biological escape evidence; CNV is `NOT_EVALUABLE` for all 32 patients and Level-2 is near-vacuous |
| "Which patient should get dual CAR-T?" | requires efficacy, protein-level target density, toxicity and manufacturability — none measured here |
| "What is the best target pair?" | the observed pair ordering is a **detection-rate artifact**; no target in the panel is depth-robust |

### The shape of the answer, stated up front

The honest synthesis Stage 12 will deliver is:

> Observed transcript-level BCMA/GPRC5D double-negativity is common at baseline (median
> 0.335 over 32 patients), but its magnitude is dominated by measurement limitations rather
> than by demonstrable biology. Depth conditioning removes most apparent co-loss enrichment
> and most apparent DN structure; the pre-registered γ-secretase hypothesis is not
> supported; genomic subclone evidence is not evaluable; and neither immune composition nor
> ligand–receptor communication survives correction. **No patient shows convergent
> measurement-robust, structurally-supported and phenotype-supported evidence, and the
> evidence axes are substantially discordant.**

---

## Non-negotiable frozen constraints

These are inherited, not re-decided. Every one is enforced by an existing or planned test.

| # | constraint |
|---|---|
| 1 | **Patient is the biological unit.** One Stage-12 row per patient. No cell-level inference. Repeated samples are within-patient evidence, never independent patients. |
| 2 | **Raw antigen measurement stays distinct from biological interpretation.** The observed DN fraction is a measurement, never an escape probability. |
| 3 | **Stage-09b tiers are measurement tiers only** — provisional by design. `robust-high` means the estimate survived the sensitivity analyses, and nothing more. |
| 4 | **Level 1 = non-random DN structure only.** Not a state, not a subclone, not antigen-specific. |
| 5 | **Level 2 = DN-associated phenotype / compatibility only**, and it is **weakly discriminative at patient level** (26 of 27 evaluable). |
| 6 | **Level 3 requires genomic evidence and is `NOT_EVALUABLE` for all 32 patients.** `CNV_SUBCLONE_NOT_SUPPORTED` is never emitted or implied. |
| 7 | **Stages 11 and 11b are exploratory, negative/confounded, and non-tier-changing.** They may not create or alter any classification. |
| 8 | **Bulk validates marginal antigen abundance only** — never the joint DN state, which bulk destroys by construction. |
| 9 | **08c coverage is descriptive transcript coverage, not therapeutic utility.** No combination is optimal, recommended or best. |
| 10 | **Normal marrow is expression context, never safety.** |
| 11 | **No hidden retuning.** No frozen threshold may be adjusted, and none may be re-derived from Stage-12 outputs. |
| 12 | **No future-derived thresholds.** Stage 12 may not invent a cutoff by looking at the synthesis it is producing. |
| 13 | **`NOT_EVALUABLE` is never converted to negative** on any axis. |
| 14 | **No composite scalar classifier**, no weighted score, no ranking, no learned model. |

Two further constraints are added by this design, from the recovered code:

| # | constraint |
|---|---|
| 15 | **Stage-10 DE must be described as a paired patient-level Wilcoxon signed-rank test on pseudobulk log-fold changes with BH correction.** The recovered producer never imports pydeseq2; the `~ patient + group` string in `pseudobulk_de_evaluability.csv` describes intent, not a fitted model. |
| 16 | **Level-2 programme results may only be quoted under the frozen both-denominators rule** (BH < 0.10 under *both* primary and sensitivity, with consistent sign). Single-denominator BH values must never be quoted as the result — γ-secretase, MYC and UPR all sit at BH ≈ 0.056–0.069 under the sensitivity denominator alone. |

---

## Inputs and provenance

### Provenance protocol (mandatory, before any input is read)

1. Load `provenance/frozen_artifacts_pre_stage12.tsv`.
2. For every Stage-12 input, recompute SHA256 and compare to the manifest.
3. **Abort on the first mismatch.** Do not proceed, do not "note and continue", and never
   regenerate the manifest to agree with a file.
4. Record commit, tag, environment name, and the resolved hash of every input into
   `stage12_input_manifest.tsv`.
5. Open every input **read-only**. Stage 12 writes only to `results/12_final_synthesis/`.

### Verified input set — 29 artifacts, all present in the committed manifest

Confirmed during this design pass: **all 29 are manifested**, and all currently verify.

| # | artifact | stage | role in synthesis |
|---|---|---|---|
| 1 | `08_.../patient_antigen_states_primary.csv` | 08 | DN fraction, marginal detection, enrichment, bootstrap, perm p (primary) |
| 2 | `08_.../patient_antigen_states_sensitivity.csv` | 08 | same under the sensitivity denominator |
| 3 | `08_.../patient_bootstrap_intervals.csv` | 08 | CI bounds, `low_n`, `single_sample` |
| 4 | `08_.../patient_conegativity_enrichment.csv` | 08 | conditioned/unconditioned enrichment |
| 5 | `08_.../patient_evidence_states.csv` | 08 | per-patient uncertainty flags |
| 6 | `08_.../repeated_sample_antigen_consistency.csv` | 08 | per-sample DN within repeated patients |
| 7 | `08_.../truncate10k_sensitivity.csv` | 08 | censoring sensitivity |
| 8 | `08_.../primary_vs_sensitivity_denominator_comparison.csv` | 08 | denominator sensitivity |
| 9 | `08_.../noise_floor_technical_zero.csv` | 08 | dropout burden (cohort × stratum) |
| 10 | `08_.../noise_floor_ambient.csv` | 08 | ambient floor, opposite-signed bias |
| 11 | `08_.../depth_stratified_null.csv` | 08 | null-scheme sensitivity |
| 12 | `08_.../risk_tier_provisional/risk_tiers_provisional.csv` | 09b | provisional measurement tier + flags |
| 13 | `08_.../risk_tier_design/patient_evidence_matrix.csv` | 09b | tier evidence inputs |
| 14 | `10_.../stage10_evidence_levels.csv` | 10 | Level-1/2/3 states, licensed language |
| 15 | `10_.../stage10_evaluability.csv` | 10 | richest per-patient L1/L2 source: n_dn, Moran's I, perm p, depth ratios, DE counts |
| 16 | `10_.../dn_local_structure_by_patient.csv` | 10 | L1 statistics per denominator |
| 17 | `10_.../dn_program_scores_by_patient.csv` | 10 | per-patient programme deltas, raw and matched |
| 18 | `10_.../level2_program_cohort_tests.csv` | 10 | cohort-level programme results (the interpretable Level-2 result) |
| 19 | `10_.../pseudobulk_de_results.csv` | 10 | 190-gene both-denominator DE set |
| 20 | `10_.../gamma_secretase_hypothesis.csv` | 10 | the pre-registered negative |
| 21 | `10_.../repeated_sample_dn_coherence.csv` | 10 | L1 repeated-sample consistency |
| 22 | `10_.../tc_like_subtype.csv` | 10 | descriptive TC-like proxy |
| 23 | `09_.../bulk_vs_sc_by_cohort.csv` | 09 | marginal validation, cohort-split |
| 24 | `09_.../normal_marrow_antigen_context.csv` | 09 | marrow expression context |
| 25 | `11_.../immune_vs_dn_measurement.csv` | 11 | composition associations (0/28 at BH<0.10) |
| 26 | `11_.../liana_verification/liana_vs_dn_associations.csv` | 11b | LIANA associations + circularity flag |
| 27 | `08_.../multi_antigen_coverage/stage12_multi_antigen_interface.csv` | 08c | the per-patient coverage interface |
| 28 | `08_.../multi_antigen_coverage/target_measurement_qc.csv` | 08c | target eligibility, technical-zero, depth sensitivity |
| 29 | `04_qc/umi_censoring_effect.csv` | 04 | the bias-direction table |

**Denominator policy for Stage 12:** the **primary** denominator is the reporting default;
the **sensitivity** denominator is carried in adjacent columns for every quantity where it
exists. **They are never averaged, and neither is selected for looking better.**

---

## Evidence axes

Six axes, kept structurally separate. Each answers a different question and has a different
failure mode.

### A. Measurement robustness

**Question:** *Is the observed DN measurement technically credible?*
**Does not answer:** whether a biological escape clone exists.

Inputs: observed DN (both denominators), marginal BCMA/GPRC5D detection, Stage-09b tier,
bootstrap CI, denominator sensitivity, truncate-10k sensitivity, technical-zero/dropout
burden, repeated-sample instability, null-scheme sensitivity, cell-count adequacy.

Frozen anchors: median DN **0.335** (0.017–0.783); tiers **4 / 28 / 0**; GPRC5D detection
spans **9.4×** across cohorts; pooled technical-zero floor BCMA 0.276 vs GPRC5D **0.620**;
truncate-10k ΔDN WU1/WU2 **0.000** vs MMRF **+0.059**; primary→sensitivity median
**+0.032** with 12/32 moving > 5 points.

### B. DN structure (Level 1)

**Question:** *Are DN cells non-randomly organized within the patient beyond the implemented
depth null?*
**Does not establish:** a phenotype or a clone.

Inputs: Level-1 state, Moran's I, kNN DN fraction, best-cluster enrichment,
depth-stratified permutation p (and the unconditioned p, for contrast), repeated-sample
consistency, evaluability.

Frozen anchors: **4 supported / 23 not supported / 5 not evaluable.**
Worked example to carry forward: `MMRF_1640` — Moran's I **0.470**, unconditioned
p **0.001**, depth-stratified p **0.499**.

### C. DN-associated phenotype (Level 2)

**Question:** *Are DN cells compatible with the cohort-level DN-associated transcriptional
phenotype?*

Inputs: Level-2 state, per-patient programme directions (matched), cohort-level programme
tests under the both-denominators rule, pseudobulk DE summary, secretory/differentiation
interpretation.

Frozen anchors: **26 supported / 1 not supported / 5 not evaluable.**
Reproducible DN-**lower** programmes: **antigen presentation** (BH 0.0011 / 0.0001),
**OXPHOS** (0.0044 / 0.0563), **interferon** (0.0899 / 0.0111). **190 DE genes** significant
under both denominators.

> **Stage 12 must state that the per-patient Level-2 state is weakly discriminative because
> 26 of 27 evaluable patients are supported, and must give the cohort-level phenotype more
> interpretive weight than the near-universal patient label.**

### D. Genomic evidence (Level 3)

**Question:** *Is there genomic evidence for a distinct DN subclone?*
**Answer:** *Not evaluable* — all 32 patients.
**Never downgrade to "no subclone."**

### E. Immune context

**Question:** *Is there independent evidence that immune context explains or tracks DN
burden?*
**Frozen answer:** *no robust independent evidence.*

Frozen anchors: **0 of 28** composition tests at BH < 0.10 (smallest 0.49); receptor-only
falsification shows **11 of 15 receptors** move down with DN burden; LIANA tested 87
interactions with **1** at BH < 0.10 — `Myeloid TNFSF13B → TNFRSF17`, **antigen-circular**;
decomposition of all 12 raw hits gives **3 confounded / 9 not reproduced / 0 LIANA-only**.

**This axis stays cohort-level and exploratory. It generates no patient labels.** Its
per-patient columns exist only to record that no patient-level immune support exists.

### F. Multi-antigen coverage

**Question:** *How does observed transcript-level uncovered fraction change under alternative
target combinations?*
**Does not answer:** therapeutic superiority, clinical recommendation, safety, efficacy.

Frozen anchors: 5 eligible / 2 not evaluable (`GPRC5D` on dropout at technical-zero
**0.620** vs a 0.50 gate fixed beforehand; `SDC1` on circularity); **no target is
depth-robust** (ρ 0.32–0.48, 3–16× stratum spread); anchor uncovered median **0.335**;
gain from adding GPRC5D **0.011**; alternative pairs lower the uncovered fraction in
**32/32** patients, median advantage 0.098 — **a detection-rate artifact.**

**Observation to carry into the write-up:** every patient's lowest-uncovered eligible pair
is `TNFRSF17 + X` (SLAMF7 ×11, FCRL5 ×11, ITGB7 ×5, CD38 ×5). The anchor's BCMA half is
never displaced; only the GPRC5D half is. That is consistent with the coverage difference
being about **GPRC5D's detection rate specifically**, not about BCMA being a poor target.

---

## Patient-level synthesis schema

**Output:** `results/12_final_synthesis/stage12_patient_evidence_matrix.csv`
**One row per patient, 32 rows.** Every field is either copied verbatim from a frozen
artifact or derived by **predeclared deterministic logic** stated in this document.

### Identity

| column | source | notes |
|---|---|---|
| `patient_id` | 08 primary | string dtype, always — int coercion produced an empty join once before |
| `cohort` | 08 primary | MMRF / WU1 / WU2 |
| `n_samples` | 08 primary | samples contributing primary-denominator cells |
| `repeated_sample_flag` | derived | `n_samples > 1` |
| `in_paper_cohort` | manifest/io | `25183` is `False` |

### Measurement (axis A)

`n_primary_cells`, `n_sensitivity_cells`, `observed_dn_primary`, `observed_dn_sensitivity`,
`bcma_detection`, `gprc5d_detection`, `provisional_measurement_tier`,
`bootstrap_ci_lower`, `bootstrap_ci_upper`, `enrichment_depth_conditioned`,
`enrichment_unconditioned`, `perm_p_depth_stratified`, plus the frozen flags
`denominator_unstable`, `depth_sensitive`, `low_n`, `repeated_sample_unstable`,
`null_scheme_sensitive`, `dropout_compatible`, `intermediate_dn`, and
`n_uncertainty_flags` (a **count**, never a score).

> `enrichment_unconditioned` is carried **only** as the paired contrast that quantifies the
> depth artifact. It is never interpreted on its own.

### Structure (axis B)

`level1_state`, `level1_evaluable`, `level1_not_evaluable_reason`, `morans_i_primary`,
`morans_i_sensitivity`, `knn_dn_frac_primary`, `best_cluster_enr_primary`,
`depth_stratified_p_primary`, `depth_stratified_p_sensitivity`,
`unconditioned_p_primary` (contrast only), `n_depth_bins`, `repeated_structure_consistent`.

### Phenotype (axis C)

`level2_state`, `level2_evaluable`, `phenotype_compatibility_note`,
`antigen_presentation_direction`, `oxphos_direction`, `interferon_direction`,
`secretory_direction`, `n_de_padj05_primary`, `depth_ratio_matched`.

> **Direction fields are signs of per-patient matched deltas, not tests.** They must be
> labelled as such and must not be treated as independent hypothesis tests — the frozen
> design tests programmes **at cohort level**, and the per-patient values feed only the
> weakly-discriminative Level-2 rule. `phenotype_compatibility_note` carries the 26/27
> disclosure on every row where `level2_state == DN_STATE_SUPPORTED`.

### Genomic (axis D)

`level3_state` (constant `CNV_SUBCLONE_NOT_EVALUABLE`), `cnv_evaluable` (constant `False`),
`cnv_not_evaluable_reason` (constant: donor negative control failed; method rejected before
disease interpretation).

### Immune (axis E)

`immune_context_summary`, `composition_supported` (expected `False` for all),
`communication_supported` (expected `False`), `liana_independent_support` (expected
`False`), `immune_confound_flag`, `liana_evaluable` (`False` for `25183`).

> These columns exist to record a **cohort-level negative**, not to differentiate patients.
> The header block of the CSV and the axis documentation must say so.

### Coverage (axis F)

`anchor_pair_uncovered`, `lowest_observed_uncovered_eligible_pair`,
`lowest_observed_uncovered_eligible_pair_value`,
`lowest_observed_uncovered_eligible_triple`,
`lowest_observed_uncovered_eligible_triple_value`, `uncovered_bcma_alone`,
`uncovered_gprc5d_alone`, `gain_from_adding_gprc5d`, `anchor_vs_alternative_note`,
`coverage_depth_sensitive` (constant `True` — no target is depth-robust),
`coverage_repeat_sensitive`, `coverage_qc_note`, `eligible_targets`,
`not_evaluable_targets`.

> **Column naming is load-bearing.** The 08c source columns are named
> `greatest_coverage_pair_descriptive` / `..._uncovered`. Stage 12 renames them to
> `lowest_observed_uncovered_eligible_pair` to remove any residual suggestion of "greatest"
> being a merit ranking. **Never `best_pair`, `optimal_pair`, or `recommended_pair`.**

### Synthesis text — generated by predeclared logic, never freeform prose

Five short text fields, each produced by a deterministic rule stated here and implemented as
a pure function with unit tests. **No model-generated sentences may enter the CSV.**

| field | rule sketch |
|---|---|
| `measurement_interpretation` | template selected by `(tier, n_uncertainty_flags, low_n)`; e.g. tier `robust-high` + 0 flags → *"observed DN estimate survived all frozen sensitivity analyses"* |
| `biological_interpretation` | template selected by `(level1_state, level2_state)`; Level-3 is always appended as *"genomic subclone evidence not evaluable"* |
| `main_uncertainty` | the highest-priority flag under a **fixed priority order** declared below — a selection, not a ranking of patients |
| `allowed_claim` | derived from `subclone.licensed_language(level1, level2, level3)`, the existing frozen function |
| `prohibited_claim` | constant per (level1, level2) combination, always including *"no genomic subclone claim"* and *"no clinical recommendation"* |

**Fixed `main_uncertainty` priority order** (declared now, before execution, and never
reordered after seeing output): `not_evaluable` → `low_n` → `dropout_compatible` →
`denominator` → `repeated_sample` → `depth` → `null_scheme` → `intermediate_dn` → `none`.

This is a **display priority for a text field**, not an ordering of patients and not a
weighting of evidence.

---

## Cohort-level synthesis schema

**Output:** `results/12_final_synthesis/stage12_cohort_summary.csv` — seven rows, one per
evidence domain, with a fixed vocabulary.

| column | content |
|---|---|
| `domain` | measurement · co-loss · structure · phenotype · genomic · immune · coverage |
| `question` | the axis question, verbatim from this design |
| `result` | the frozen quantitative result |
| `strength_of_evidence` | one of the five fixed values below |
| `major_control` | the control that produced the result |
| `major_limitation` | the binding limitation |
| `allowed_interpretation` | permitted wording |
| `prohibited_interpretation` | forbidden wording |

### Pre-populated content (all values verified against frozen artifacts)

| domain | result | strength | major control | major limitation | allowed | prohibited |
|---|---|---|---|---|---|---|
| **measurement** | observed DN common; median **0.335** (0.017–0.783), 32 patients, 21,906 primary cells | `STRONG` | denominator, depth, truncate-10k, repeated-sample and threshold sensitivity | GPRC5D technical-zero **0.620**; WashU 10k censoring | "observed transcript-level DN burden" | "true target-negative clone prevalence" |
| **co-loss** | unconditioned median **1.052** (max 4.606) → conditioned **1.009** (max 1.750); **4/32** significant, all MMRF | `NOT_SUPPORTED` (as biology) | cohort-specific depth-stratified permutation null | significance tracks depth/power, not necessarily biology | "most apparent co-negativity enrichment is attributable to depth" | "coordinated biological antigen co-loss" |
| **structure** | **4 / 23 / 5** (supported / not / not evaluable) | `SUPPORTED_WITH_CAVEATS` | depth-stratified permutation + adaptive bins | conservative by design; 5 not evaluable; `MMRF_1640` shows the failure mode | "non-random DN organization in a minority of patients" | "subclone", "pre-existing resistant clone" |
| **phenotype** | cohort-level DN-lower antigen presentation / OXPHOS / interferon; 190 DE genes under both denominators; per-patient **26/27** | `SUPPORTED_WITH_CAVEATS` | depth-matched cells before pseudobulk; both-denominators rule | confounded with a broad secretory/differentiation shift; BCMA and GPRC5D are themselves secretory-pathway-dependent | "cohort-level DN-associated transcriptional state" | "antigen-specific escape mechanism"; "individualized phenotype risk" |
| **genomic** | **32 / 32 `CNV_SUBCLONE_NOT_EVALUABLE`** | `NOT_EVALUABLE` | donor negative control, run before disease interpretation | donor false-positive span 0.0–50.6% at z>3 | "genomic subclone evidence is not evaluable" | "absence of a subclone"; any use of the word *subclone* |
| **immune** | **0/28** composition at BH<0.10; LR receiver-state confounded; LIANA **1/87** and antigen-circular | `NOT_SUPPORTED` | receptor-only falsification; antigen-circularity test | n≈32 with a confounder correlated with the predictor | "no robust independent immune association" | "immune-evasion mechanism"; "PDCD1→CD274 axis" |
| **coverage** | alternatives lower observed uncovered fraction in **32/32**, median 0.098; **no target depth-robust**; 2/7 not evaluable | `EXPLORATORY` | target-specific technical-zero floors; predeclared eligibility gate | detection differs 1.8–2.8×; GPRC5D dropout; SDC1/TNFRSF17 selection dependence | "greatest observed transcript-level malignant-cell coverage among evaluated combinations" | "optimal / recommended / best pair"; "GPRC5D is redundant"; any safety claim |

### Evidence-strength vocabulary — exact semantics, fixed before execution

| value | means |
|---|---|
| `STRONG` | Directly measured, reproduced under **every** frozen sensitivity analysis that applies to it, and not dependent on a contested assumption. |
| `SUPPORTED_WITH_CAVEATS` | Reproduced under the frozen rule (e.g. both denominators), but conditional on a stated technical assumption, or holding in only a minority of patients. |
| `EXPLORATORY` | Descriptive; predeclared but not inferentially tested, or tested without surviving correction; may not license a mechanism or a decision. |
| `NOT_SUPPORTED` | Tested under a predeclared control and failed it. This is a **negative result about the hypothesis**, not about evaluability. |
| `NOT_EVALUABLE` | The measurement or assay could not be assessed at all. **Carries no evidential weight in either direction and may never be read as negative.** |

**No patient-level numeric confidence value may be invented on top of this vocabulary.**

---

## Concordance / discordance analysis

This is the analytical heart of Stage 12: making the axes' disagreement visible.

### The joint contingency, computed from frozen artifacts during this design pass

**Tier × Level-1**

| | L1 NOT_EVALUABLE | L1 NOT_SUPPORTED | L1 SUPPORTED |
|---|---:|---:|---:|
| robust-high | 0 | **4** | **0** |
| uncertain | 5 | 19 | 4 |

**Tier × Level-2**

| | L2 NOT_EVALUABLE | L2 NOT_SUPPORTED | L2 SUPPORTED |
|---|---:|---:|---:|
| robust-high | 0 | 0 | 4 |
| uncertain | 5 | 1 | 22 |

**Level-1 × Level-2**

| | L2 NOT_EVALUABLE | L2 NOT_SUPPORTED | L2 SUPPORTED |
|---|---:|---:|---:|
| L1 NOT_EVALUABLE | **5** | 0 | 0 |
| L1 NOT_SUPPORTED | 0 | 1 | 22 |
| L1 SUPPORTED | 0 | 0 | 4 |

### Four structural facts, verified

1. **Only 5 of 18 possible (tier × L1 × L2) cells are occupied:**

   | tier | L1 | L2 | n | patients |
   |---|---|---|---:|---|
   | uncertain | NOT_SUPPORTED | SUPPORTED | **18** | 25183, 27522, 37692, 47491, 56203, 57075, 58408, 77570, 81012, 83942, MMRF_1413, MMRF_1505, MMRF_1640, MMY21940, MMY22933, MMY34600, MMY40511, MMY47218 |
   | uncertain | NOT_EVALUABLE | NOT_EVALUABLE | **5** | 59114, MMRF_1325, MMRF_1424, MMRF_1641, MMRF_1777 |
   | robust-high | NOT_SUPPORTED | SUPPORTED | **4** | MMRF_1267, MMY18273, MMY74196, MMY98423 |
   | uncertain | SUPPORTED | SUPPORTED | **4** | MMRF_1720, MMRF_2038, MMY34339, MMY80649 |
   | uncertain | NOT_SUPPORTED | NOT_SUPPORTED | **1** | MMY67868 |

2. **Level-1 and Level-2 evaluability are perfectly coupled.** The same 5 patients are
   `NOT_EVALUABLE` on both; **zero** patients are non-evaluable on exactly one. Evaluability
   is therefore a single joint property, not two independent ones — a fact Stage 12 should
   state rather than implying two separate evaluability checks.

3. **Level-2 is nearly a function of evaluability.** Among the 27 evaluable patients, 26 are
   `SUPPORTED`. The axis carries almost no cross-patient information.

4. **Claim K holds: no patient is simultaneously measurement-robust-high, Level-1 supported
   and Level-2 supported.** The four `robust-high` patients are all `L1 NOT_SUPPORTED`; the
   four `L1 SUPPORTED` patients are all measurement-`uncertain`. **The two sets are
   disjoint.**

### Supporting overlap facts

- **Uncertainty-flag load:** 4 patients carry 0 flags (exactly the `robust-high` four), 14
  carry 1, 8 carry 2, 4 carry 3, 2 carry 4. Most frequent flags: `dropout_compatible` (15),
  `intermediate_dn` (14), `denominator` (9).
- **Repeated-sample status:** 26 `not_assessable_for_agreement`, 4 `agree`, 2 `disagree`.
  All 4 `robust-high` patients are `not_assessable_for_agreement` (single-sample) — so their
  intervals are the optimistic kind, and Stage 12 must say so beside the tier.
- **Repeated patients in the primary denominator: 7**, not 8 — `60359` has **zero
  primary-denominator cells** and does not appear in the Stage-08 primary table at all.

### Output

`stage12_concordance_matrix.csv` (a tidy long-form rendering of the above) plus Figure 2.
The matrix is **categorical only**. No axis may be encoded as an ordinal integer, because an
integer encoding invites summation — which is the composite score by another route.

---

## Repeated-patient handling

**Output:** `results/12_final_synthesis/stage12_repeated_patient_summary.csv`

Seven patients contribute more than one sample to the primary denominator. Per patient:
`n_samples`, per-sample `n_cells`, observed DN **range** (never a mean),
`measurement_tier_stability`, `level1_repeated_consistency`, coverage range, and an explicit
`small_denominator_flag` for any sample below the frozen 20-cell floor.

Frozen values to reproduce:

| patient | n samples | primary cells | observed DN (primary) | note |
|---|---:|---:|---:|---|
| `27522` | 6 | 1,250 | 0.297 | per-sample DN 0.275 / 0.000 / 0.500 / 0.509 / 0.571 / 0.511 on n = 1,134 / **1** / **10** / 53 / **7** / 45 |
| `47491` | 2 | 486 | 0.202 | |
| `56203` | 2 | 360 | 0.103 | second timepoint exists only because of the `56203_1` truncation repair |
| `58408` | 2 | 168 | 0.357 | |
| `59114` | 2 | 59 | 0.339 | also `NOT_EVALUABLE` on both Stage-10 levels |
| `81012` | 2 | 373 | 0.386 | |
| `83942` | 2 | 4,705 | 0.053 | one patient, two **protocols** (`83942` WU1 + `MMY83942` WU2), not a serial timepoint |

**Design rules:**
- **Never average across samples within a patient.** Report the range and the per-sample n.
- **Show tiny denominators explicitly.** `27522`'s apparent 0.275 → 0.571 "trajectory" rests
  on samples of 1, 10, 7 and 45 cells; the design requires those n values to be printed
  adjacent to every DN value so the variation is not read as disease evolution.
- **Keep the two lists distinct.** Multi-sample patients (bootstrap-relevant) vs serial `_N`
  timepoint patients (longitudinal-relevant). `83942` is in the first and not the second.
- `60359` is reported as **no primary-denominator cells**, not as a missing row.

---

## Claim ladder

**Output:** `results/12_final_synthesis/stage12_claim_ladder.csv` — a central deliverable.

All statuses below were **verified against frozen artifacts during this design pass.**

| # | Claim | Evidence required | Current status | Can Stage 12 say it? |
|---|---|---|---|---|
| 1 | Observed DN burden exists at baseline | per-patient DN over a defended denominator | median **0.335**, 32 patients, 21,906 cells | **YES** — as *observed transcript-level* burden |
| 2 | The DN measurement is robust | survival of all frozen sensitivity analyses | **4 / 32** patients (0 uncertainty flags); 0 `robust-low` | **YES, for 4 patients only**, and only as a *measurement* statement |
| 3 | DN co-loss is enriched beyond depth | conditioned enrichment > 1 with a depth-stratified null | median **1.009**; 4/32 significant, all in the deepest cohort | **NO** as a cohort claim; **only as "4 patients retain significance, in the cohort where the test has power"** |
| 4 | DN cells are structurally organized | Level-1 support under the depth-stratified null, both denominators | **4 / 23 / 5** | **YES for those 4 patients**, as *non-random organization* only |
| 5 | A DN-associated phenotype exists | cohort-level programmes reproducible under both denominators | AP / OXPHOS / interferon DN-lower; 190 DE genes | **YES at cohort level**; **NO** as individualized evidence (26/27) |
| 6 | DN reflects an antigen-specific escape state | separation of the phenotype from a general secretory shift | not separable — BCMA and GPRC5D are secretory-pathway-dependent | **NO** |
| 7 | A genomic DN subclone exists | CNV evidence | **32/32 `NOT_EVALUABLE`** | **NO** — and its negation is equally prohibited |
| 8 | An immune-evasion mechanism operates | independent immune association surviving correction | 0/28 composition; LR receiver-confounded; LIANA 1/87 and circular | **NO** |
| 9 | An alternative target combination is superior | comparable detection + efficacy/safety evidence | alternatives lower uncovered fraction 32/32, but detection differs 1.8–2.8× and no target is depth-robust | **NO** — only "greatest observed transcript-level coverage among evaluated combinations" |
| 10 | A clinical recommendation follows | protein-level target density, efficacy, toxicity, manufacturability | none measured | **NO** |
| 11 | γ-secretase mediates DN escape | pre-registered programme reproducible under both denominators | BH **0.3865 / 0.0694** — fails the both-denominators rule; **direction negative**, opposite to prediction | **NO** — report as a clean pre-registered negative |
| 12 | Convergent multi-axis evidence exists in any patient | measurement-high ∧ L1 ∧ L2 | **n = 0** | **NO** |

### Statements A–K, verification status

Every statement in the brief was checked against the frozen data:

| | statement | verified |
|---|---|---|
| A | DN cells common at baseline | ✔ median 0.335 |
| B | magnitude highly sensitive to depth, especially GPRC5D | ✔ 9.4× cohort spread; technical-zero 0.620 |
| C | co-loss enrichment largely collapses after depth conditioning | ✔ 1.052 → 1.009; max 4.606 → 1.750 |
| D | only a minority show non-random organization | ✔ 4/32 |
| E | cohort-level less-secretory / lower antigen-presentation state | ✔ AP, OXPHOS, IFN all DN-lower under the both-denominators rule |
| F | phenotype is not evidence of antigen-specific escape | ✔ (interpretive; grounded in the secretory-dependence confound) |
| G | γ-secretase not supported | ✔ **with a required nuance — see below** |
| H | genomic subclone evidence not evaluable | ✔ 32/32 |
| I | immune analyses give no robust independent support | ✔ 0/28; 1/87 circular |
| J | alternatives reduce uncovered fraction but ranking is precluded | ✔ 32/32; detection differs 1.8–2.8× |
| K | no patient has convergent measurement-high + L1 + L2 | ✔ n = 0 |

> **Required precision on G.** γ-secretase is **not** "p > 0.10 everywhere". It is BH
> **0.3865** under the primary denominator and **0.0694** under the sensitivity denominator,
> so it **fails the frozen both-denominators rule**. MYC (0.3865 / 0.0563) and UPR
> (0.3865 / 0.0563) sit in exactly the same position. Stage 12 must quote the rule, not a
> single denominator — quoting the sensitivity value alone would misrepresent three
> programmes as near-significant findings.

---

## Uncertainty hierarchy

Three **separate** hierarchies. They are reported in three distinct blocks and **never
combined into one uncertainty score.**

### 1. Measurement uncertainty (per patient, from frozen flags)

dropout / technical-zero burden · depth sensitivity · WashU 10k censoring · denominator
instability · low n · repeated-sample instability · null-scheme sensitivity.

### 2. Biological uncertainty (per patient + cohort)

absence of Level-1 structure in 23 of 32 · non-specificity of the Level-2 phenotype
(inseparable from a secretory/differentiation shift) · Level-3 non-evaluability in all 32 ·
the near-vacuity of the per-patient Level-2 label.

### 3. External-validity uncertainty (cohort-level only)

**transcript vs protein** — CAR-T binds surface protein; BCMA is actively shed by
γ-secretase and GPRC5D transcript correlates imperfectly with surface density ·
**normal marrow vs whole body** — keratinized-tissue liability is structurally unobservable
here · **untreated baseline vs post-therapy selection** — this cohort is baseline, so no
selection has occurred and no post-treatment inference is available ·
**cohort/protocol confounding** — depth, site, chemistry and censoring all move together
and cannot be separated in this design.

---

## Figures

Five designed; **if effort must be cut, produce 2, 4 and 5**, per the brief.

### Figure 1 — Evidence architecture (schematic)

Six axes drawn as **parallel lanes, not a ladder**, with an explicit annotation that they are
not a linear escalation and are never summed. Genomic is drawn greyed-out with
`NOT_EVALUABLE`. Immune and coverage sit as contextual side panels.

### Figure 2 — Patient evidence matrix *(priority)*

Categorical heatmap: rows = 32 patients, columns = axes + uncertainty flags.
**Sort rule, declared now: provisional measurement tier, then patient_id ascending.**
No other sort. No clustering, no seriation, no optimisation of visual pattern.
Colour encodes **category identity only** — supported / not supported / not evaluable / flag
present — using a palette with no light-to-dark ordinal ramp, so no axis reads as ordinal.

### Figure 3 — Measurement vs Level-1 structure

x = observed DN fraction (primary), y = Moran's I (primary), point shape = Level-1 state,
annotation = provisional tier. **Purpose is to show discordance**, and the caption states
that the four `robust-high` patients are all Level-1 `NOT_SUPPORTED`. `MMRF_1640` is
labelled with both its p-values (0.001 unconditioned, 0.499 depth-stratified).

### Figure 4 — Stage-10 cohort phenotype *(priority)*

Compact dot/forest plot of the **seven predeclared programmes**, matched effects, both
denominators, with BH values. **No programme may be added.** Raw (unmatched) effects are
shown in a paler adjacent series — this is the depth lesson made visual, and it is where
γ-secretase's collapse is visible.

### Figure 5 — Multi-antigen coverage *(priority)*

Descriptive uncovered-fraction summary: anchor and eligible alternatives. Mandatory
annotations: **GPRC5D `COVERAGE_NOT_EVALUABLE` (technical-zero 0.620)**, **SDC1 excluded on
circularity**, **TNFRSF17 selection-dependence caveat**, **no target is depth-robust**.
Ordering of combinations on the axis is **alphabetical**, never by uncovered fraction, so the
figure cannot be read as a ranking. Caption must carry the detection-artifact explanation.

---

## Output files

All under a new namespace: **`results/12_final_synthesis/`**. No upstream file is touched.

| file | content |
|---|---|
| `stage12_patient_evidence_matrix.csv` | 32 rows; the authoritative patient table |
| `stage12_cohort_summary.csv` | 7 rows; one per evidence domain |
| `stage12_claim_ladder.csv` | 12 rows; claim / evidence required / status / permitted |
| `stage12_concordance_matrix.csv` | tidy long-form axis cross-tab |
| `stage12_repeated_patient_summary.csv` | 7 repeated patients + the `60359` note |
| `stage12_uncertainty_register.csv` | the three uncertainty hierarchies, kept separate |
| `stage12_input_manifest.tsv` | input · stage · SHA256 · producer · role · verified |
| `stage12_design_snapshot.md` | this design, frozen at execution time, with commit/tag/env |
| `figures/` | Figures 1–5 |
| `stage12_summary.md` | the narrative write-up |

### `stage12_summary.md` — narrative structure

1. What we measured
2. How robust the measurement is
3. Whether DN cells are structured
4. Whether they share a phenotype
5. Whether a genomic subclone is demonstrated
6. Whether immune context adds support
7. What alternative transcript-level target coverage shows
8. What remains unresolved
9. What claims are justified
10. What claims are explicitly rejected

Sections 9 and 10 are **generated from `stage12_claim_ladder.csv`**, not written freehand,
so the narrative cannot drift from the table.

---

## Tests

Written **before** execution, in `tests/test_stage12_synthesis.py`.

### A. Freeze / integrity
- every Stage-12 input hash matches `provenance/frozen_artifacts_pre_stage12.tsv`;
- Stage 12 writes **only** into `results/12_final_synthesis/`;
- no Stage 04–11b or 08c artifact changes across a Stage-12 run (hash before/after).

### B. No hidden classifier (AST/source-level, docstrings stripped)
Ban in Stage-12 source: weighted sums over axis columns, `.sum(axis=1)` over evidence
columns, `np.dot`/`@` against a weight vector, `rank`/`argsort`/`sort_values` on any
evidence or score column, `PCA`/`TSNE`/`UMAP`/`KMeans`, any `sklearn` estimator, any
`score`/`utility`/`weight`/`points`/`composite` identifier bound to a patient-level value.
Mirror the existing `test_g_no_weighted_aggregate_exists_anywhere_in_the_module` pattern
from `tests/test_coverage.py`.

### C. Tier invariance
`stage12_patient_evidence_matrix.csv` reproduces **4 robust-high / 28 uncertain /
0 robust-low**, with the same four patient IDs.

### D. Level invariance
Level-1 **4 / 23 / 5**; Level-2 **26 / 1 / 5**; Level-3 **32 `CNV_SUBCLONE_NOT_EVALUABLE`**;
the specific L1-supported set `{MMRF_1720, MMRF_2038, MMY34339, MMY80649}` and the single
L2-not-supported patient `MMY67868` are unchanged.

### E. Biological-unit invariance
Exactly one row per patient (32 rows, 32 unique IDs); no output has a cell-level index; the
repeated-patient table never collapses samples into an average.

### F. Coverage invariance
No output column or text field contains `optimal`, `best`, `recommended`, `superior`, or
`preferred` applied to a combination; GPRC5D and SDC1 QC caveats present on every coverage
row; `coverage_depth_sensitive` is `True` for all patients.

### G. Claim guards (behavioural where possible, literal where wording *is* the invariant)
Fail if any generated text contains: `confirmed subclone`, `genomic subclone`,
`immune escape mechanism`, `immune evasion` (unqualified), `optimal target pair`,
`GPRC5D redundant`, `clinically safe`, `high-risk patient`, `risk classifier`,
`IMMUNE_EVASION_CONFIRMED`, `CNV_SUBCLONE_NOT_SUPPORTED`.
Behavioural companions: assert `licensed_language()` output for every patient never exceeds
what `(level1, level2, level3)` permits; assert no row asserts genomic evidence while
`cnv_evaluable == False`.

### H. Evaluability guards
`NOT_EVALUABLE` never maps to `False`/`0`/`not supported` in any derived column; the five
jointly-non-evaluable patients keep `NOT_EVALUABLE` on both L1 and L2.

### I. Text-generation determinism
Every synthesis text field is produced by a pure function of frozen inputs; running twice
yields byte-identical output; the functions are unit-tested against fixed inputs.

### J. Denominator separation
Primary and sensitivity columns both present for every quantity that has them; no column is
a mean of the two; no code path selects between them by comparing outcomes.

---

## Prohibited analyses and claims

**Prohibited analyses.** New p-values · new regressions · new clustering · new dimensionality
reduction · new pathway or programme tests · new patient subgroup tests · new thresholds ·
re-running any frozen stage · re-deriving depth strata · imputation of any kind · averaging
across denominators · averaging repeated samples · any ordinal encoding of an evidence axis ·
any weighted or learned aggregation.

**Prohibited claims.** The twelve carried forward from the frozen record — no confirmed
genomic subclone (and no assertion of its absence) · no antigen-specific escape mechanism ·
no immune-evasion mechanism · no clinical target-combination recommendation · no safety
claim from marrow · no claim that GPRC5D is redundant · no composite risk classifier · no
equating transcript coverage with therapeutic superiority · no cross-cohort DN comparison ·
no dropout-corrected DN · no treatment of `DN_STATE_SUPPORTED` as strong individual evidence ·
no description of Stage-10 DE as pydeseq2.

---

## Decision on patient synthesis categories

> ### **Recommendation: NO. Stage 12 produces the evidence matrix only, with no categorical patient synthesis label.**

This decision was made **after** inspecting the frozen cross-tabulation, as required, and is
recorded here before execution.

### The evidence

Only **5 of 18** possible (tier × L1 × L2) cells are occupied, and the occupancy is extreme:
**18 / 5 / 4 / 4 / 1**.

Any category scheme over this structure would be one of two things:

**(a) A lossless relabelling.** The five cells are already fully described by three frozen
columns. Inventing names like `MEASUREMENT_HIGH_BIOLOGY_UNRESOLVED` for the `robust-high`
cell adds no information, and adds real risk: a compound name reads as a **rank** (it sounds
like a tier), invites sorting, and hides the fact that its members are `L1 NOT_SUPPORTED` —
which is the most important thing about them.

**(b) An arbitrary aggregation.** Any scheme with fewer than five categories must merge
cells, and every available merge is unprincipled — most obviously merging `MMY67868`
(the single `L2 NOT_SUPPORTED` patient) into either neighbour, which would either erase the
one patient who fails the Level-2 rule or create a **singleton category**, the textbook sign
of an over-fitted taxonomy.

### Three further reasons

1. **Level-2 carries almost no information.** 26 of 27 evaluable patients are supported, so
   any category incorporating L2 is effectively a category about *evaluability*. The
   contingency confirms this: L1 and L2 non-evaluability are **perfectly coupled** (5
   patients, both axes, zero patients non-evaluable on exactly one).
2. **The axes are disjoint where it matters.** `robust-high` and `L1 SUPPORTED` share **no**
   patients. A single label per patient must therefore either privilege one axis over the
   other or fuse two incompatible meanings — and the disjointness *is the finding*. A
   category would conceal exactly what Stage 12 exists to reveal.
3. **Categories drift toward risk language.** Every plausible name for the `robust-high`
   cell implies elevated risk. With CNV `NOT_EVALUABLE` for all 32 patients, no patient-level
   label can be given a biological meaning, and a label with no biological meaning that
   *sounds* biological is worse than no label.

### What replaces categories

An optional **`evidence_profile`** column containing the **verbatim concatenation** of the
three frozen states, e.g.
`tier=uncertain|L1=NOT_SUPPORTED|L2=SUPPORTED`.

This is explicitly **not** a category: it is lossless, mechanically derived, creates no new
vocabulary, cannot be ordered, and takes exactly the five values already enumerated above.
Its only purpose is grouping rows in a table or figure. If even this proves to invite
misreading in review, drop it — the three source columns are sufficient.

**The five profiles, with membership, are already documented in the concordance section
above and require no new naming.**

---

## Exact execution order for Stage 12

1. **Verify provenance.** Hash all 29 inputs against the committed manifest. **Abort on any
   mismatch.** Write `stage12_input_manifest.tsv`.
2. **Snapshot the design.** Copy this document to
   `results/12_final_synthesis/stage12_design_snapshot.md`, stamped with commit, tag,
   environment name and UTC time. Nothing after this point may change the design.
3. **Load frozen inputs read-only**, with `dtype={"patient": str}` everywhere.
4. **Assert the frozen invariants immediately** — tier 4/28/0, L1 4/23/5, L2 26/1/5, L3 32
   NE, 32 unique patients. **Abort on any deviation** before doing any work.
5. **Assemble the patient evidence matrix** (axes A–F, verbatim copies plus predeclared
   derivations). No transformation beyond renaming and sign extraction.
6. **Generate synthesis text** via the predeclared pure functions.
7. **Build the concordance matrix** from the assembled table.
8. **Build the repeated-patient summary**, ranges only, with per-sample n printed.
9. **Build the cohort summary** from the pre-populated content in this design.
10. **Build the claim ladder** from the pre-populated content in this design.
11. **Build the uncertainty register**, three separate blocks.
12. **Render figures** under the declared sort/order rules.
13. **Generate `stage12_summary.md`**, with sections 9–10 derived from the claim ladder.
14. **Re-verify** that no upstream artifact changed (re-hash all 393).
15. **Run the full test suite**, including the new Stage-12 tests.
16. **Report**; do not stage or commit without review.

---

## Stop conditions

Execution halts immediately, with nothing written, if any of these occur:

1. any upstream hash fails verification;
2. tier or state counts differ from 4/28/0 · 4/23/5 · 26/1/5 · 32 NE;
3. patient IDs change, or the row count is not 32;
4. a synthesis category would require a post-hoc threshold;
5. any output would imply genomic evidence while CNV is `NOT_EVALUABLE`;
6. any algorithm ranks target combinations;
7. any upstream artifact is modified;
8. a new statistical test appears necessary;
9. a text field would contain a prohibited claim;
10. `NOT_EVALUABLE` would be coerced to a negative or numeric value.

**In every case: stop and report. Do not adapt the design to proceed.**

---

## New statistics: none required

Stage 12 needs **no new statistical test**. Every quantity it reports already exists in a
frozen artifact. Specifically:

- concordance/discordance is **counting**, not testing — and deliberately so. A formal test
  of association between axes (e.g. Fisher on tier × L1) is **not proposed**: with 4 and 4
  positives among 32 it would be badly underpowered, it would invite a p-value where the
  contingency table is the honest presentation, and it would be a **new test in a synthesis
  stage**. The disjointness is reported as an observed fact, not a significant one.
- the repeated-patient summary reports **ranges**, not tests.
- the claim ladder cites **existing** BH values.

If a new test ever seems indispensable, the design requires documenting: the exact question,
why frozen results cannot answer it, the predeclared null, the replication unit, the
correction method, and the post-hoc overfitting risk — **and stopping for review before
implementing it.**

---

## Two reconciled narrative points

Incorporated into the design as required.

### 08c chronology

**Logical pipeline position** and **historical execution chronology** are different things
and must be stated separately:

> 08c is a **supplemental arm attached to Stage 08**, and it consumes frozen Stage-08
> artifacts, which is why it carries a letter rather than a number. **Historically it was
> executed later, during the pre-Stage-12 completion pass (2026-08-26), after Stages 09,
> 09b, 10 and 11.**

**Its stage number is not evidence that it ran before Stage 09.** What guarantees
correctness is not chronology but isolation: 08c hashes every frozen artifact it consumes
into `frozen_upstream_digests.json`, and `tests/test_coverage.py` re-checks those digests,
so a later edit to a frozen stage fails a test. Stage 12 must describe 08c this way.

### Plasma mixing

The permitted interpretation is the cautious one:

> **Plasma mixing is strongly confounded by cohort-specific depth, protocol, and
> non-recoverable WashU censoring; therefore cross-cohort plasma-state comparisons from the
> integrated embedding are not trusted.**

Stage 12 must **not** state that poor plasma Harmony mixing is simply expected because
malignant plasma cells are patient-private. That framing was explicitly retired: the three
largest plasma clusters are **one per cohort, each spanning ~30 patients**, which is
inconsistent with patient-private structure — a patient-private clone would fragment into
~41 clusters. The evidence points to a sampling/censoring asymmetry (MMRF plasma median
22,477 UMIs vs WU1 5,036; MMRF's two largest plasma clusters have 68% and 88% of cells above
a ceiling the WashU deposits cannot contain), which no method can recover.

---

## Final readiness checklist

| # | item | status |
|---|---|---|
| 1 | All 29 Stage-12 inputs identified and present in the committed manifest | ✔ verified |
| 2 | All 393 frozen artifacts currently verify against the manifest | ✔ verified |
| 3 | Joint contingency inspected **before** any category decision | ✔ 5 of 18 cells occupied |
| 4 | Category decision made and justified | ✔ **NO categories** |
| 5 | Claim-ladder statements A–K verified against frozen artifacts | ✔ all 11, with the G nuance recorded |
| 6 | Evidence-strength vocabulary defined with exact semantics | ✔ 5 values |
| 7 | Figure sort/order rules declared before execution | ✔ tier → patient_id; alphabetical combinations |
| 8 | Text-generation logic predeclared and deterministic | ✔ template rules + fixed priority order |
| 9 | Output namespace defined and isolated | ✔ `results/12_final_synthesis/` |
| 10 | Test plan written before execution | ✔ suites A–J |
| 11 | Stop conditions enumerated | ✔ 10 conditions |
| 12 | New statistics required | ✔ **none** |
| 13 | Recovered implementation treated as authoritative (Wilcoxon, not pydeseq2) | ✔ constraint 15 |
| 14 | 08c chronology reconciled | ✔ |
| 15 | Plasma-mixing interpretation reconciled | ✔ |
| 16 | No composite score anywhere in the design | ✔ and test-enforced |
| 17 | Stage-12 notebook created | ✘ **deliberately not yet** |

**Design status: COMPLETE. Stage 12 is ready to be implemented against this design.
It has not been executed.**
