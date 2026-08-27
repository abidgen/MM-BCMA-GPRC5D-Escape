# Pre-Stage-12 Codex audit — remediation record

Remediation date: 2026-08-26
Audit remediated: `docs/pre_stage12_codex_audit.md` (1 CRITICAL, 2 HIGH, 8 MEDIUM, 6 LOW)
Scope: **reproducibility and provenance only.**

> **No scientific result changed.** No threshold was retuned, no patient moved tier, no
> Level-1/Level-2/CNV state changed, no Stage-08c target eligibility changed, no
> Stage-11/LIANA interpretation changed, and Stage 12 was not started. Nothing under
> `results/` was written during this pass — verified by mtime and by all 393 manifest
> hashes matching (`sha256sum -c provenance/frozen_artifacts_pre_stage12.sha256`).

## Closure table

| Finding | Status | Action taken | Evidence | Remaining caveat |
|---|---|---|---|---|
| **C1 — Stage 07–10 frozen outputs have no committed producer** (CRITICAL) | `RESOLVED` | The **exact historical** producers were recovered verbatim from the Claude Code session transcripts and committed under `production/` (39 drivers across stages 06–10, each with a provenance header naming session, UTC execution time and outputs). Recovery was preference-order 1, not reconstruction: **no file is a `RECONSTRUCTED_PRODUCTION_DRIVER`.** Nothing was re-executed. | All 293 frozen Stage 06–11 artifacts resolve to a producer (`production/RECOVERY_INVENTORY.tsv`, `provenance/frozen_artifacts_pre_stage12.tsv`). Three independent validations: (a) `results/07_malignant_plasma/v_clone_membership/antigen_circularity_invariance.py` survived on disk and the transcript-recovered copy is **byte-identical**, 4,517 bytes both; (b) every driver's execution timestamp matches the mtime of the artifacts it writes within its runtime (e.g. `s08a` ran 16:06:57Z, `patient_antigen_states_primary.csv` mtime 16:29:37Z); (c) all 188 extracted Python payloads parse, none truncated. `tests/test_production_paths.py` (35 tests). | The session scratchpad files themselves were deleted; the transcripts are the surviving record and **no hash of the originals was ever taken**. The three checks above are strong evidence of faithful recovery, not a cryptographic guarantee. Also: the recovered chain passes state between steps through `/tmp/*.pkl` files that no longer exist, so a rerun must go start-to-finish rather than resume midway. |
| **H1 — frozen-artifact digests are not independently anchored** (HIGH) | `PARTIALLY_RESOLVED` | Committed `provenance/frozen_artifacts_pre_stage12.tsv` (393 rows, stages 04–11b) and `provenance/frozen_artifacts_pre_stage12.sha256`, **outside** ignored `results/`. Each row carries stage, path, size, SHA256, artifact role, committed producer path, producer commit, freeze date, environment and reproducibility class. Committed as-built environment exports for all five conda environments (`conda list --explicit`, `conda env export --no-builds`, `pip freeze`, R package tables) plus `ENVIRONMENT_SUMMARY.md`. `provenance/README.md` documents verification and the archival policy. | `sha256sum -c` passes on all 393 files; `test_f6_manifest_hashes_match_the_local_frozen_artifacts` (slow marker) reproduces it in-suite; `test_f1`–`test_f5`, `test_f7` pin schema, resolution, stage coverage, producer existence and that `results/` stays ignored. | **Immutable external storage is not yet provisioned.** The manifest anchors the freeze *in Git*, which is what was missing, but the artifacts themselves remain local. `provenance/README.md` states the four requirements for the external deposit and asks that the bucket/DOI and the manifest's own hash be recorded there once it exists. Until then the honest position is: artifacts local, present state checkable. |
| **H2 — Level-2 classification is near-vacuous; the guard test is tautological** (HIGH) | `RESOLVED` | Two separate things, both fixed without touching the frozen rule. **(1) Code:** `subclone.level2_state()` now calls a new `validate_program_names()` which raises `UnknownProgramError` for any name outside `config.LEVEL2_PROGRAMS`; validation runs on every branch, including the non-evaluable one. **(2) Semantics:** the weak-discriminator disclosure is now carried in the function docstring and in `CLAUDE.md` as binding on Stage 12 — *"the per-patient Level-2 state is a weak compatibility label, not a discriminative risk classifier; the cohort-level phenotype is the scientifically informative result."* **The 26/27 rule was not retuned.** | `tests/test_production_paths.py::test_a1`–`test_a7` — `test_a1`/`test_a2`/`test_a4` **actually call** `level2_state` with an invented name and assert the raise, replacing the tautological set-membership assertion. `test_a7` records why the frozen numbers were already safe: every program loop in `production/stage10/s10c_*.py` and `s10d_*.py` iterates `config.LEVEL2_PROGRAMS`, so no post-hoc program had a route into `repro` or `hits`. | **The library function was not the production call site.** The recovered `s10d` driver inlines the Level-2 rule rather than calling `subclone.level2_state()`. The new validation therefore hardens the library for future callers; it did not retroactively guard the frozen run. What guarded the frozen run is the structural bound in `test_a7`. This is stated rather than glossed. |
| **MEDIUM #3 — pseudobulk test asserts array shape, not the replication unit** | `RESOLVED` | Added `tests/test_pseudobulk_production.py` (17 tests): an end-to-end synthetic cohort, deliberately depth-confounded (DN cells shallower), run through the recovered driver's loop shape using the committed primitives. Tests behaviour: exactly one pseudobulk per patient per group; duplicating a patient's cells adds no replicate; repeated samples collapse into their patient and a multi-sample patient's pseudobulk draws from both samples; depth matching happens on cells **before** summation and demonstrably shrinks the depth gap; an unmatchable patient contributes nothing rather than a biased row; antigen columns are dropped **before** normalisation (shown by perturbing an antigen by 5,000 counts and asserting CPM is unchanged when dropped first and visibly changed when dropped after). | `tests/test_pseudobulk_production.py` A1–A3, B1–B3, C1–C4, D1–D4, E1–E3, all passing. | The synthetic cohort exercises the design, not the real data; it cannot and does not re-derive any frozen DE number. |
| **MEDIUM #4 — `test_g2_level2_state_accepts_only_predeclared_program_names` is tautological** | `RESOLVED` | Superseded by behavioural tests that invoke the function (see H2). The original test is left in place as a cheap freeze guard on the constant. | `test_a1`, `test_a2`, `test_a4` in `tests/test_production_paths.py`. | None. |
| **MEDIUM — stale notebook/result names in README and pipeline doc** | `RESOLVED` | `README.md`: the stage table now names `07_malignant_plasma`, `08_dual_antigen_escape`, `09_bulk_validation`, `10_dn_coherence`, `11_immune_context` and the lettered arms `05b`, `08c`, `09b`, `11b`, with correct result directories. `mm_dual_antigen_escape_pipeline.md`: 10 stale notebook/result references corrected and a dated naming banner added at the top. `CLAUDE.md`: the two remaining stage-09/10 prose references corrected. | `grep` for the five planning-era names returns only the explanatory banner. | `mm_dual_antigen_escape_pipeline.md` remains a planning-era document; its banner says `CLAUDE.md` is current where they disagree. Only naming was corrected — no scientific prose was rewritten. |
| **MEDIUM — contradictory `jupytext --sync` instruction** | `RESOLVED` | `CLAUDE.md:441` no longer says sync keeps the pair in lockstep; it now points at the notebook-pairing rule and mandates explicit `jupytext --to notebook` in the `.py → .ipynb` direction. | `grep -n jupytext CLAUDE.md` — the forbidding rule and the review workflow now agree. | None. |
| **MEDIUM — source-of-truth contradiction (`results/*.csv` called the source of truth while ignored)** | `RESOLVED` | Replaced in both places (`CLAUDE.md`, `docs/stage-results.md`) with: *frozen scientific outputs are local/archived artifacts authenticated by the committed freeze manifest; committed production code and provenance metadata define how those artifacts were produced and verified.* Both now state that an ignored local file is not by itself durable provenance and outranks prose only once its hash matches. | `grep "source of truth"` in `CLAUDE.md` and `docs/stage-results.md`. | Full closure depends on H1's external deposit. |
| **LOW — `docs/stage-results.md` title says "stages 01 through 05b"** | `RESOLVED` | Retitled "run output for stages 01 through 11b (incl. 05b, 08c, 09b)". | `head -1 docs/stage-results.md`. | None. |
| **LOW — environment-count phrasing** | `RESOLVED` | `docs/environments.md` now states **"Five environments are built; six YAML specifications exist"**, naming `envs/env-composition.yml` as written-but-unbuilt, and gained an "As-built exports" section pointing at `provenance/environments/`. | `head docs/environments.md`, tail of the same file. | None. |
| **LOW — writable Numba cache undocumented** | `RESOLVED` | Recorded in `provenance/environments/ENVIRONMENT_SUMMARY.md` and in the new `docs/environments.md` section, with the exact `NUMBA_CACHE_DIR=…` invocation. | Both files. | None. |
| **MEDIUM — no committed manifest declares which annotation tree is accepted** | `RESOLVED` | The freeze manifest hashes `results/06_annotation/**` (74 rows) as the accepted tree and names its producer. The nine rejected `06_annotation_*` trees are deliberately **not** manifested, which is what distinguishes them from the accepted endpoint on a clean clone. | `provenance/frozen_artifacts_pre_stage12.tsv`, stage `06`. | The rejected trees remain on disk, ignored and unmanifested. Quarantining them under a versioned archive manifest stays an after-Stage-12 improvement. |
| **MEDIUM — duplicated depth-binning wrappers (08 local vs `antigen.py` vs Stage-10 adaptive)** | `OPEN` (deferred by design) | Not touched. The audit itself recommends consolidating only **after** the freeze, with a regression-equivalence proof; doing it now would edit frozen production paths. | — | Carried forward as an after-Stage-12 item. The Stage-10 adaptive amendment is documented and intentional, not drift. |
| **MEDIUM — ambiguous `total_counts` naming (full-reference vs intersected space)** | `OPEN` (deferred by design) | Not touched; renaming fields would change code that produced frozen results. Stage 11 already repairs the distinction explicitly (Amendment 3). | — | After-Stage-12 rename to `total_counts_full_ref` / `total_counts_intersection`. |
| **MEDIUM — `communication.ols_association` uses nominal df after a pseudoinverse** | `OPEN` (deferred by design) | Not touched. The audit found no rank defect in current results and recommends adding rank/condition diagnostics in future verification **without** retrospectively changing frozen values. | — | After-Stage-12. Stage 11 is exploratory and non-tier-changing, so this cannot move a classification. |
| **MEDIUM — `PATIENT_ALIASES` rests on concordant S1 characteristics, not a direct identifier** | `OPEN` (documentation, Stage 12) | Not touched; it is a limitation to carry, not a defect to fix. | — | Stage 12 must carry the collapse rationale and the split-ID sensitivity in its limitations. |
| **LOW — implicit index alignment / `.astype(str)` on categoricals** | `OPEN` (deferred by design) | Not touched in the recovered drivers — editing them would destroy their value as historical evidence. | — | Applies to any **new** code (Stage 12 included): assert one-to-one merges and non-missing values before conversion. |
| **LOW — stale/rejected annotation and checkpoint artifacts coexist locally** | `PARTIALLY_RESOLVED` | The manifest now distinguishes accepted from rejected by inclusion (see above). | — | Physical quarantine deferred to after Stage 12. |
| **LOW — literal source-substring tests test naming more than invariants** | `PARTIALLY_RESOLVED` | The new suites are behavioural where it matters (H2 rejection, pseudobulk replication unit, antigen-drop ordering, depth-matching order, tier isolation by signature and by construction order). One pre-existing string test (`test_module_never_calls_transcriptional_coherence_a_subclone`) was narrowed to inspect **string constants** rather than `repr`, because a class `repr` carries the module path `mm_escape.subclone` and made the test pass or fail on where code lives rather than what it claims. | `tests/test_production_paths.py`, `tests/test_pseudobulk_production.py`, `tests/test_dn_coherence.py:219`. | Literal freeze guards are retained alongside, as the audit recommends. |
| **LOW — `.ipynb`-vs-`.py` execution path references** | `OPEN` | Not addressed beyond the jupytext fix. | — | Minor; after Stage 12. |

## One finding raised and cleared during remediation

A new test initially failed because the recovered `production/stage09b/s09b2_provisional_tiers.py`
**does** read two Stage-09 bulk tables — which appeared to contradict "no bulk may enter the
tier driver". This was investigated before anything was changed.

**It is not a defect.** Tier assignment happens in the `for _,r in EV.iterrows()` loop via
`RT.final_tier(ev)`, where `ev` is built only from Stage-07/08 evidence-matrix columns.
Bulk is attached **after** `T=pd.DataFrame(rows)` under the driver's own comment
`# ---- Stage-09 context: REPORT-ONLY, appended AFTER tiers exist ----`, producing four
descriptive columns (`s09_bulk_available`, `s09_bulk_TNFRSF17_tpm`, `s09_bulk_GPRC5D_tpm`,
`s09_marginal_context`). The decision module `src/mm_escape/risk_tiers.py` contains no
bulk, cohort, coherence or p-value term at all.

The test was rewritten to pin the **real** invariant rather than ban a filename:
`test_c3c_stage09_bulk_enters_09b_only_as_report_only_context_after_tiers_exist` asserts the
ordering and that the `ev` dict carries no bulk term; `test_c3a`/`test_c3a2` assert the
decision module and its signatures are bulk-, cohort- and coherence-free;
`test_c3d_frozen_tier_membership_is_unchanged` pins **4 robust-high / 28 uncertain /
0 robust-low**. This matches CLAUDE.md exactly: Stage 09 is interpretation context, never
an additional scoring axis.

## Verification performed

| check | result |
|---|---|
| full suite, `pytest -q -m "not slow"` (`mm-core`) | **566 passed, 1 skipped, 5 deselected** |
| slow manifest-hash test | **1 passed** |
| independent `sha256sum -c` over the manifest | **393/393 OK, 0 mismatches** |
| frozen artifacts modified during this pass | **none** (`find results -newermt '2026-08-26 19:50'` empty) |
| required Stage-12 upstream artifacts represented in the manifest | **yes**, all layers (see `provenance/README.md`) |
| recovered drivers use no future-stage input | enforced by `test_c2`, `test_c3b`, `test_c4` |
| recovered Python drivers parse | 39/39 (`test_b2`) |

## Status against the audit's five required fixes

1. **Resolve C1** — done, by verbatim recovery rather than reimplementation.
2. **Resolve H1** — manifest and environment provenance committed; external immutable
   storage still to be provisioned.
3. **Bind Stage-12 semantics to H2** — the weak-discriminator wording is now binding in
   `CLAUDE.md` and in the function's own docstring.
4. **Re-run tests; add end-to-end synthetic production tests** — done; 52 new tests.
5. **Reconcile source-of-truth and stale paths** — done.

**Recommendation: READY FOR RE-AUDIT**, with the two disclosed residual caveats (external
immutable storage not yet provisioned; the recovered drivers are transcript-derived and
validated by a byte-identical control plus timing rather than by a hash of the originals).

## 2026-08-26 non-blocking re-audit closure addendum

This addendum closes only the three non-blocking findings in
`docs/pre_stage12_codex_reaudit.md`. No scientific output, threshold, tier, patient state,
frozen artifact, or analysis result was changed, and Stage 12 was not started.

| Re-audit finding | Closure | Correction |
|---|---|---|
| **N1 — recovered-producer commit placeholders** | `RESOLVED` | Replaced all 197 `pending-commit(this-repair-pass)` values in `provenance/frozen_artifacts_pre_stage12.tsv` with `ea26a949c40bc070c1cf4983ea936252fa36d145`, the commit containing those recovered producers. `provenance/README.md` now distinguishes that repository commit from historical execution-time provenance and retains the missing-original-hash caveat. Artifact paths, sizes, hashes, dates, environments and reproducibility classes were not changed. |
| **N2 — frozen Stage-10 DE mislabeled as pydeseq2 / `~ patient + group`** | `RESOLVED IN CURRENT DOCUMENTATION` | Current methods and environment documentation now state that frozen Stage-10 DE used depth-matched patient pseudobulks, paired patient-level DN-versus-comparator log-fold changes, two-sided Wilcoxon signed-rank tests, and BH correction. Patient remains the biological unit. pydeseq2 remains recorded as installed but unused for this result. The inaccurate historical label remains untouched in two frozen artifacts and their verbatim recovered producer; `provenance/README.md` records that fact explicitly. |
| **Residual “CSVs are the source of truth” sentence** | `RESOLVED` | `mm_dual_antigen_escape_pipeline.md` now states that ignored/local outputs are authoritative only after matching the committed provenance manifest and are not durable provenance by themselves. |
