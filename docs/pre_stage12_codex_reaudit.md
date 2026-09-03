# Focused pre-Stage-12 Codex re-audit

> **Editorial note, added after the fact:** this re-audit cleared Stage 12 to begin
> ("safe to begin, with non-blocking caveats"). Stage 12 subsequently ran to completion
> on 2026-08-27 (`results/12_final_synthesis/`, design in `docs/stage12_design.md`).
> This document is preserved verbatim as the historical checkpoint record and has not
> been edited to reflect that outcome.

Audit date: 2026-08-26  
Repository: `/media/wrath/CART_mm_dual_antigen`  
Audited commit: `ea26a949c40bc070c1cf4983ea936252fa36d145`  
Audited tag: `pre-stage12-remediated`  

This was a read-only verification pass except for creation of this report. Stage 12 was
not started; no scientific code, threshold, tier, patient state, frozen artifact, or
environment was changed or regenerated.

## Executive summary

The remediation resolves the original Stage-12 blockers sufficiently to begin final
synthesis. The formerly missing Stage 07-10 producers are now present, committed,
syntactically valid, and reviewable. Direct code tracing supports the documented stage
order and shows patient-level inference, antigen exclusion, separate primary/sensitivity
denominators, and report-only use of bulk in Stage 09b. The recovered code is best called
`RECOVERED_WITH_PROVENANCE_CAVEAT`: transcript heredocs, execution timestamps, parsing,
and a byte-identical surviving Stage-07 control provide strong support, but no historical
hash exists for every original scratchpad producer.

All 393 manifested frozen artifacts exist and pass both SHA256 and recorded-size checks.
The standalone SHA256 file and TSV contain the same 393 paths and hashes. The committed
manifest therefore authenticates the current local frozen state. It does not prove
historical immutability, and the only artifact copies remain local/ignored pending an
external immutable deposit. That residual H1 limitation is explicit and is acceptable
for starting Stage 12, but it remains a real archival caveat.

H2 is resolved. Runtime validation is bound to `config.LEVEL2_PROGRAMS`, invalid names
raise on all evaluability branches, and behavioral tests invoke the function with invalid
names. Frozen calls remain 26 supported of 27 evaluable patients (plus five not evaluable),
and the documentation correctly calls this a weak compatibility label rather than a
patient-level risk classifier.

The full suite passed: **570 passed, 2 skipped, 0 deselected, 6 warnings**. No new CRITICAL
or HIGH defect was found. Three non-blocking MEDIUM precision issues remain: the main
pipeline still contains one unqualified “CSVs are the source of truth” sentence; the
Stage-10 narrative/environment record says pydeseq2 although the recovered frozen producer
uses paired patient-level Wilcoxon tests; and 197 manifest rows retain the placeholder
`pending-commit(this-repair-pass)` rather than the now-known recovery commit.

## Stage-12 readiness

`YES, WITH NON-BLOCKING CAVEATS`

The remaining H1 external-storage limitation does not prevent current local artifact
authentication: the manifest is committed at the audited tag, all local files match it,
and no evidence of mutation was found. Stage 12 must consume the hashed artifacts and
must preserve the documented weak Level-2 semantics. The three MEDIUM metadata/narrative
issues below should be corrected before public final reporting, but none changes an input
number or patient-level result.

## Original finding closure

| Finding | Previous status | Current status | Evidence | Blocking? |
|---|---|---|---|---|
| C1 — missing Stage 07-10 production paths | CRITICAL / open | **RESOLVED** | Forty-one recovered files under `production/`; all are tracked, all Python payloads parse, all shell payloads pass `bash -n`, and the artifact manifest maps frozen outputs to producers. Direct traces are summarized below. | No |
| H1 — frozen-artifact provenance not independently anchored | HIGH / open | **PARTIALLY_RESOLVED** | The committed 393-row TSV and SHA256 manifests exactly authenticate all current local artifacts and record sizes, stages, producers, dates, environments, and reproducibility classes. Immutable external storage is still absent, and 197 producer-commit cells retain a placeholder. | No for Stage 12; still an archival caveat |
| H2 — Level-2 validation / semantics | HIGH / open | **RESOLVED** | `src/mm_escape/subclone.py:validate_program_names()` binds to `config.LEVEL2_PROGRAMS`; `level2_state()` validates before every branch. `tests/test_production_paths.py:test_a1`-`test_a6` call the behavior. Frozen states remain unchanged. | No |
| Pseudobulk behavioral-test weakness | MEDIUM / open | **RESOLVED WITH TEST-SCOPE CAVEAT** | `tests/test_pseudobulk_production.py` uses synthetic counts and the same committed binning/matching primitives to demonstrate one row per patient/group, repeated-sample pooling, pre-sum matching, antigen removal before normalization, and deterministic seeds. It reproduces the recovered loop rather than importing the top-level producer; several linkage guards are source inspection. | No |
| Tautological Level-2 test | MEDIUM / open | **RESOLVED** | Tests now pass invented and plausible-but-unfrozen names to the runtime function and require `UnknownProgramError`, including the not-evaluable branch. | No |
| Documentation/source-of-truth contradictions | MEDIUM / open | **MOSTLY RESOLVED** | The main project document and `docs/stage-results.md` now explain committed code + manifest + ignored/archived artifacts and explicit checksum verification; Jupytext direction, title, notebook names, and environment count were corrected. One stale sentence remains in the main pipeline narrative. | No |

## Production-driver recovery assessment

All files below are committed at `ea26a94`. Independent syntax checks used
`PYTHONPYCACHEPREFIX` under `/tmp` for Python and `bash -n` for shell, without writing to
the repository. The drivers are historical producers and deliberately overwrite frozen
namespaces if run; `production/README.md` correctly says not to execute them merely to
check a number. Several chains depend on ephemeral `/tmp/*.pkl`, so reproducibility
requires running each stage chain continuously from its first step.

### Stage 07

- **Producer present/outputs:** `production/stage07/s07a`-`s07j` cover reference
  availability and the CNV gate, light-chain dominance, antigen-excluded CNV gene space,
  donor CNV calibration, V/J availability and donor/repeated-sample checks, clone
  membership, primary/sensitivity denominators, and antigen-perturbation invariance. Their
  declared and actual writes resolve to manifested `results/07_malignant_plasma/**`
  artifacts.
- **Inputs/order:** the code reads accepted Stage-06 annotation/integration inputs and no
  future stage. `s07h_clone_membership.py` groups by `patient_id`; repeated samples are
  audited separately rather than counted as independent patients.
- **Scientific invariants:** clone membership uses light-chain/V evidence, not antigen
  expression. `s07c` excludes `TNFRSF17` and `GPRC5D` from CNV features. `s07j` checks
  antigen perturbation invariance. CNV remains `CNV_SUBCLONE_NOT_EVALUABLE`; no downstream
  normal/negative interpretation was found. Primary and sensitivity denominators are
  emitted separately.
- **Conclusion:** producer is present, committed, valid, correctly ordered, and
  patient-disciplined.

### Stage 08

- **Producer present/outputs:** `production/stage08/s08a`-`s08e` produce patient antigen
  states, cohort/global depth-stratified nulls, bootstrap intervals, technical-zero and
  ambient controls, truncate-10k sensitivity, repeated-sample tables, denominator
  comparison, evidence states, and the Stage-08 checkpoint.
- **Inputs/order:** the chain reads accepted Stage-05/06 cell data and Stage-07 clone
  membership/evaluability, with no Stage-09/10/11 input. Writes remain in the Stage-08
  namespace.
- **Scientific invariants:** raw antigen counts are used; patient is the inferential unit;
  repeated samples contribute within patient; fixed seeds are present; primary and
  sensitivity outputs remain distinct. Shared `mm_escape.antigen` and
  `mm_escape.subclone` depth machinery is called where documented.
- **Conclusion:** producer is present, committed, valid, correctly ordered, and
  patient-disciplined.

### Stage 09

- **Producer present/outputs:** `production/stage09/s09a`-`s09d` cover bulk TPM loading,
  exact scRNA/bulk pairing, marginal TNFRSF17/GPRC5D correlations, truncate-10k checks,
  repeated samples, and normal-marrow context. Writes match manifested Stage-09 paths.
- **Inputs/order:** Stage 09 uses Stage-08 state carried through the historical chain plus
  raw deposited bulk files and the frozen GTF; it does not read Stage 09b/10/11.
- **Scientific invariants:** MMRF runs are averaged rather than summed; one assessable
  matched timepoint per patient enters correlations; bulk tests marginal abundance only;
  normal marrow is donor-level context. No joint-DN validation occurs in bulk.
- **Conclusion:** producer is present, committed, valid, correctly ordered, and
  patient-disciplined.

### Stage 09b

- **Producer present/outputs:** `production/stage09b/s09b1`-`s09b4` cover the Stage-07/08
  evidence matrix, frozen provisional-tier construction, threshold sensitivity, cohort
  diagnostic, and historical final-to-provisional relabel. The resulting manifested
  membership remains **4 robust-high, 28 uncertain, 0 robust-low**.
- **Inputs/order:** tier evidence is constructed before Stage 10/11 and no Stage-10/11
  value is read. The driver does read Stage-09 tables at load time, but, as verified below,
  they are not visible to the decision function and are attached only after tiers exist.
- **Scientific invariants:** `src/mm_escape/risk_tiers.py` accepts no bulk, cohort,
  coherence, or evidence-level argument. Threshold variants remain sensitivity outputs,
  not a best-result selection.
- **Conclusion:** producer is present, committed, valid, correctly ordered, and
  patient-disciplined.

### Stage 10

- **Producer present/outputs:** `production/stage10/s10a`-`s10g` cover antigen-excluded
  Level-1 embeddings, adaptive depth bins and stratified permutations, depth matching,
  patient pseudobulks, paired DE, frozen programs and gamma-secretase, decoupler full
  space, TC-like subtype, and evidence-level construction. Writes map to manifested
  `results/10_dn_coherence/**` outputs.
- **Inputs/order:** the scientific drivers read Stage-08 and earlier products, plus
  same-stage intermediates, never Stage 11. Antigen columns are removed before
  normalization/features. The recovered DE loop iterates eligible patients, pools all
  their samples, matches cells within adaptive depth bins, sums one DN and one comparator
  vector per patient, and tests patient-paired log-fold changes.
- **Scientific invariants:** Level 1, Level 2, and Level 3 terminology is separated; no
  patient is called a genomic subclone. Frozen programs come only from
  `config.LEVEL2_PROGRAMS`. Gamma-secretase is not retuned. The cohort-level phenotype is
  distinct from the weak 26/27 per-patient compatibility rule.
- **Operational caveat:** `s10g` and the Stage-09b relabel script are historical
  freeze/migration steps and can edit narrative/source paths as well as results. They are
  reviewable historical evidence, not safe rerun entry points; the README warning is
  therefore load-bearing.
- **Conclusion:** producer is present, committed, valid, correctly ordered, and
  patient-disciplined. The pydeseq2-versus-Wilcoxon narrative mismatch is documented under
  new issues.

## Recovered-code provenance confidence

**Classification: `RECOVERED_WITH_PROVENANCE_CAVEAT`.**

Evidence supporting faithful recovery:

- every file identifies the transcript UUID, original scratchpad path/form, and execution
  UTC time;
- transcript Bash heredocs preserve the full producer payloads;
- all extracted payloads parse;
- recorded execution windows agree with artifact mtimes;
- the surviving historical
  `results/07_malignant_plasma/v_clone_membership/antigen_circularity_invariance.py` is
  byte-identical (4,517 bytes) to the recovered copy; and
- the recovered paths generate the expected schemas/names and reproduce the documented
  dependency structure.

This is strong evidence of exact historical recovery, but it is not a cryptographic proof
for every producer. The original scratchpad files were deleted and no original hash was
taken for each one. `EXACT_HISTORICAL_RECOVERY_WITH_STRONG_SUPPORT` would overstate what
can now be proven; `RECONSTRUCTED` would understate the transcript evidence.

## Frozen-artifact verification

- Expected manifest rows: **393**; unique paths: **393**.
- SHA256 entries: **393**; the SHA file and TSV contain identical path/hash mappings.
- Files present: **393/393**.
- Recorded size matches: **393/393**.
- `sha256sum -c provenance/frozen_artifacts_pre_stage12.sha256`: **393/393 OK**, exit 0.
- Coverage by stage: 04 (70), 05 (20), 05b (10), 06 (74), 07 (44), 08 (16), 08c
  (19), 09 (12), 09b (12), 10 (78), 11 (26), 11b (12).
- Reproducibility classes: 191 `REPRODUCIBLE_FROM_COMMITTED_CODE`, 180
  `REPRODUCIBLE_WITH_EXTERNAL_INPUT`, 22 `ARCHIVED_ONLY`, 0
  `NOT_FULLY_REPRODUCIBLE`.

No manifested frozen scientific artifact changed during this re-audit, and no mismatch
suggests that remediation rewrote one. Because no independently committed pre-remediation
hash existed, historical immutability before this checkpoint cannot be proven
retrospectively; that is the residual H1 caveat, not evidence that a frozen number changed.

H1 status is **PARTIALLY_RESOLVED**. The committed manifest authenticates current local
files and distinguishes accepted Stage-06 outputs plus Stage-08c/11b artifacts. External
immutable storage is still absent. This is sufficient for Stage 12 provided all inputs are
verified against the committed manifest before synthesis and the limitation remains
disclosed.

## Environment provenance

`provenance/environments/` contains nonempty `conda list --explicit`, no-build environment
YAML, and `pip freeze` exports for all five built environments. R package/version tables
are present for `mm-qc` and `mm-annotation`. `ENVIRONMENT_SUMMARY.md` records capture date,
OS/kernel, GPU model/driver, stage-to-environment mapping, Python, scanpy/anndata, NumPy,
pandas, SciPy, pydeseq2, decoupler, LIANA, scDblFinder, SingleCellExperiment, R/rpy2,
scVI, torch, and the relevant cross-environment version differences. CUDA is diagnostically
identified for the GPU-dependent 05b arm; Stages 07-10 are recorded as CPU work.

Assessment: **exact enough for reasonable reconstruction and strong enough for diagnosis**.
The explicit conda URLs/builds are the strongest freeze record; no container is required
for this checkpoint. External resources (raw GEO/dbGaP availability, celldex/SingleR
downloads, LIANA resource) remain correctly represented as external reproducibility
caveats.

One metadata error remains: the environment summary labels pydeseq2 0.5.4 as “used by
stage 10 pseudobulk DE,” while the recovered producer imports SciPy Wilcoxon and never
imports pydeseq2. The package was installed, but this frozen DE result did not use it.

## Stage-09b ordering check

Direct trace of `production/stage09b/s09b2_provisional_tiers.py` confirms the invariant:

1. `EV` is the Stage-07/08 patient evidence matrix.
2. The loop builds an `ev` dictionary solely from cell count, DN fraction, truncation,
   enrichment, confidence interval, null-scheme, and repeated-sample fields.
3. `RT.final_tier(ev)` is called before `T = pd.DataFrame(rows)`.
4. Only after `T` exists does the explicit `Stage-09 context: REPORT-ONLY` block append
   `s09_*` bulk columns.

`RT.final_tier()` cannot see bulk, cohort, Stage 10, Stage 11, or coherence fields. No
Stage-10/11 path occurs in Stage 09b. The hoisted Stage-09 reads are therefore not future
leakage or tier evidence. Frozen membership is unchanged at 4/28/0.

## Test-suite result

Command:

```text
NUMBA_CACHE_DIR=/tmp/mm_numba_cache_reaudit conda run -n mm-core pytest -q -rs
```

Result: **570 passed, 2 skipped, 0 deselected, 6 warnings** in 41.25 seconds on the
confirmatory run (the first equivalent run completed in 51.08 seconds).

Skips:

- one repaired-gene-axis branch in `tests/test_io.py`, explicitly covered by repair tests;
- one scDblFinder bridge test requiring the `mm-qc` R environment rather than `mm-core`.

Warnings were four upstream Scanpy/Python deprecations, one intentional empty-marker
runtime warning, and one upstream pandas/Scanpy deprecation. The only workaround was a
writable Numba cache under `/tmp`; this is already documented and is not a scientific
failure.

Production-path test quality:

- **Substantive behavioral:** invalid Level-2 vocabulary, frozen tier membership,
  pseudobulk replication unit, repeated-sample pooling, depth matching effect and order,
  antigen perturbation before normalization, and deterministic matching.
- **Schema/provenance:** manifest schema, stage coverage, path existence, producer
  presence, and full slow artifact hash verification.
- **Useful literal/AST guards:** stage namespace, future-stage path references, fixed seed
  literals, denominator separation, CNV label vocabulary, producer headers, and actual
  driver-loop linkage. These are valuable freeze guards but are not dynamic execution of
  the historical top-level scripts.

The pseudobulk synthetic helper recreates the recovered loop using the exact committed
`adaptive_depth_bins` and `depth_matched_indices` functions. It does not import or execute
the top-level `s10e` producer, so the “driver iterates patients,” call ordering, seed, and
feature-drop linkage are source-string checks. This limitation is disclosed in the test
docstring and does not make the behavioral assertions tautological.

The production-path namespace test is also primarily literal/regex based; it does not
sandbox-run dangerous historical drivers. That is appropriate for frozen isolation here,
but it cannot prove every computed/dynamic path. Direct code tracing found no future-stage
scientific input or unintended prior-stage result write.

## New issues introduced by remediation

### N1 — producer commit placeholder remains in 197 manifest rows

- **File/section:** `provenance/frozen_artifacts_pre_stage12.tsv`,
  `producer_code_commit`; `provenance/README.md`, manifest schema.
- **Evidence:** 197 rows say `pending-commit(this-repair-pass)` after the repair was
  committed as `ea26a949...`. Other values are 74 rows at `625b7be`, 70 at `d6a85c6`, 30
  at `65487cb`, and 22 archived documents with no code commit.
- **Why it matters:** the audited tag couples the manifest and producer paths, so current
  authentication remains effective, but the per-row field does not fulfill its stated
  purpose of naming the exact commit that carries each recovered producer.
- **Severity/type:** **MEDIUM — provenance/documentation issue**.
- **Recommended action:** in a provenance-only follow-up, replace the placeholder with the
  exact recovery commit or explicitly define it as “same commit as this manifest” and pin
  the tag/commit in `provenance/README.md`. Do not regenerate artifact hashes.

### N2 — Stage-10 DE implementation is Wilcoxon, not pydeseq2

- **Files/sections:** `production/stage10/s10e_pseudobulk_de_decoupler_tc.py`, pseudobulk
  DE block; `production/stage10/s10d_hypotheses_and_evidence_levels.py`, evaluability
  `design` string; `mm_dual_antigen_escape_pipeline.md`, Stage 10;
  `provenance/environments/ENVIRONMENT_SUMMARY.md`, analysis-critical packages;
  `docs/environments.md`, `mm-core` dependencies.
- **Evidence:** the recovered producer imports `scipy.stats.wilcoxon`, forms paired
  patient CPM log-fold changes, and applies per-gene Wilcoxon plus BH correction. It never
  imports or calls pydeseq2. The evaluability CSV nevertheless labels the design
  `~ patient + group`, and narrative/environment text calls the analysis pydeseq2.
- **Why it matters:** the actual method is still a valid patient-paired, depth-matched
  analysis and the frozen output is authenticated; this is not pseudoreplication and does
  not show that any result is numerically wrong. It does matter for exact methods reporting
  and clean-room reproduction.
- **Severity/type:** **MEDIUM — documentation/scientific-method provenance issue**.
- **Recommended action:** before public synthesis text is finalized, describe the frozen
  implementation exactly as paired Wilcoxon over patient pseudobulk log-fold changes and
  remove claims that pydeseq2 produced these frozen values. Do not rerun or substitute a
  method during Stage 12.

No newly introduced CRITICAL or HIGH scientific/provenance issue was found.

## Stage-by-stage reproducibility update

| Stage | Rating | Basis |
|---|---|---|
| 07 | `REPRODUCIBLE_WITH_DOCUMENTED_EXTERNAL_ARTIFACT_CAVEAT` | Complete committed chain and manifested Stage-06 inputs/Stage-07 outputs; transcript-recovery and local-only artifact caveats remain. |
| 08 | `REPRODUCIBLE_WITH_DOCUMENTED_EXTERNAL_ARTIFACT_CAVEAT` | Complete committed chain for all requested core/sensitivity/control outputs; manifested inputs and outputs; local-only freeze. |
| 09 | `REPRODUCIBLE_WITH_DOCUMENTED_EXTERNAL_ARTIFACT_CAVEAT` | Complete committed marginal-bulk/context chain; rerun also needs deposited raw bulk/GTF resources and continuous temporary-state chain. |
| 09b | `REPRODUCIBLE_WITH_DOCUMENTED_EXTERNAL_ARTIFACT_CAVEAT` | Complete committed tier/evidence chain; correct report-only bulk ordering; historical relabel step and local-only artifacts are documented. |
| 10 | `REPRODUCIBLE_WITH_DOCUMENTED_EXTERNAL_ARTIFACT_CAVEAT` | Complete committed Level-1/2/DE/program/decoupler/TC/evidence chain; patient and antigen invariants hold; method-label and transcript-recovery caveats remain. |
| 08c | `REPRODUCIBLE_WITH_DOCUMENTED_EXTERNAL_ARTIFACT_CAVEAT` | Assessment materially improves because 19 artifacts, including the Stage-12 interface and QC, are now authenticated by the committed manifest; producer remains the committed notebook. |
| 11 | `REPRODUCIBLE_WITH_DOCUMENTED_EXTERNAL_ARTIFACT_CAVEAT` | Twenty-six custom-analysis artifacts are authenticated; producer was already committed. External/local artifact archival caveat remains. |
| 11b | `REPRODUCIBLE_WITH_DOCUMENTED_EXTERNAL_ARTIFACT_CAVEAT` | Twelve LIANA verification artifacts, environment/API versions, and producer are authenticated; bundled/external resource and local archive caveats remain. |

These ratings are deliberately more conservative than the manifest's per-artifact
`REPRODUCIBLE_FROM_COMMITTED_CODE` class: they include the repository-level fact that the
large frozen inputs/outputs are local and not yet in immutable external storage.

## Remaining pre-Stage-12 requirements

No blocking scientific recomputation or retuning is required. Before Stage 12 reads any
input:

1. verify the required artifact subset against the committed SHA256 manifest;
2. pin the synthesis work to `ea26a949...` / `pre-stage12-remediated` plus this re-audit;
3. preserve the current 4/28/0 tier membership and 26/27 weak Level-2 compatibility
   semantics; and
4. do not execute recovered historical drivers into frozen namespaces.

Non-blocking corrections recommended before public final reporting:

- update the 197 placeholder producer-commit entries without changing artifact hashes;
- correct pydeseq2/`~ patient + group` wording to the implemented paired Wilcoxon method;
- correct `mm_dual_antigen_escape_pipeline.md:35`, which still says ignored CSVs alone
  are “the source of truth.” The accurate rule already present in the main project
  document and `docs/stage-results.md` is that frozen local/archived artifacts are
  authoritative only after matching the committed manifest; and
- provision immutable external storage when available and record its versioned pointer
  plus the manifest's own digest.

## Final verdict

**Stage 12 is safe to begin with non-blocking caveats.** C1 is resolved. H2 is resolved.
H1 is partially resolved but sufficient for this checkpoint because all current local
frozen artifacts are authenticated by a committed manifest at the audited tag. The
recovered Stage 07-10 paths are complete and reviewable enough for synthesis, with strong
but non-cryptographic historical-recovery support. Patient-level replication,
antigen-circularity controls, tier ordering, and frozen-result integrity remain intact.

Finding counts:

- **Unresolved original findings:** CRITICAL 0; HIGH 1 (H1 partial, non-blocking for Stage
  12); MEDIUM 1 (one residual source-of-truth sentence); LOW 0 in this focused scope.
- **Newly introduced findings:** CRITICAL 0; HIGH 0; MEDIUM 2; LOW 0.
- **Frozen numerical result appears invalid:** no.
- **Any frozen number changed:** no evidence of change; 393/393 current frozen artifacts
  match the committed freeze exactly.
- **New CRITICAL/HIGH issue:** none.
- **Stage-12 readiness:** `YES, WITH NON-BLOCKING CAVEATS`.
- **Git status before report creation:** clean, branch `main` ahead of `origin/main` by one
  commit.
- **Expected git status after this report:** only untracked
  `docs/pre_stage12_codex_reaudit.md`; nothing staged or committed.
