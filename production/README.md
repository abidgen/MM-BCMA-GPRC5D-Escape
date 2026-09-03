# `production/` — the recovered producers of the frozen Stage 06–10 results

**Status: recovery pass, 2026-08-26. Nothing here was executed during recovery, and no
frozen artifact was regenerated, recomputed or altered.**

## Why this directory exists

The pre-Stage-12 Codex audit (`docs/pre_stage12_codex_audit.md`) raised one CRITICAL
finding, **C1**: substantial Stage 07–10 frozen outputs had no producer in the repository.
The authoritative notebooks (`notebooks/07_*.py` … `notebooks/10_*.py`) are *reports* —
they read the frozen tables and narrate them. The code that actually computed those
tables was written into per-session scratchpad files, executed, and then lost when the
scratchpads were cleaned. A clean checkout could not recreate or code-review any of it.

The audit was explicit that this is an archival failure, not evidence of a wrong number:
it found **no direct evidence that any frozen numerical result is invalid**.

## What was recovered, and how

The producers were recovered **verbatim** from the Claude Code session transcripts under
`~/.claude/projects/-media-wrath-CART-mm-dual-antigen/*.jsonl`. Those transcripts record
the full text of every `Bash` tool call, including the `cat > … <<'PYEOF'` heredocs that
created each production script and the inline `python - <<EOF` blocks that ran directly.

This is **recovery preference order 1 — exact historical production code**. It is not a
reimplementation from memory, and none of these files is a `RECONSTRUCTED_PRODUCTION_DRIVER`.

Three independent checks support that claim:

1. **A surviving control.** `results/07_malignant_plasma/v_clone_membership/antigen_circularity_invariance.py`
   was written into `results/` at the time and still exists on disk. The transcript-recovered
   copy of it is **byte-identical** (4,517 bytes, both). The same extraction produced every
   other file here.
2. **Timing.** Every recovered script's execution timestamp matches the mtime of the
   artifacts it writes, to within the script's runtime (UTC in the transcript, EDT on disk).
   Example: `s08a` ran 2026-08-25T16:06:57Z; `patient_antigen_states_primary.csv` has
   mtime 2026-08-25 12:29:37 EDT = 16:29:37Z, 22 minutes later.
3. **Syntax.** All 188 extracted Python payloads parse; none is truncated.

Every one of the 293 frozen Stage 06–11 artifacts resolves to a producer. The mapping is
`RECOVERY_INVENTORY.tsv` here, and the per-artifact record is
`../provenance/frozen_artifacts_pre_stage12.tsv`.

## Layout

| directory | stage | steps |
|---|---|---|
| `stage06/` | 06 annotation (supporting) | `s06a`–`s06i` — the accepted driver is `notebooks/06_annotation.py`; these are the cluster-23, TRBC and benchmark arms |
| `stage07/` | 07 malignant plasma | `s07a`–`s07j` — CNV gate, light-chain dominance, CNV gene set, donor CNV calibration, IG V/J audit, clone membership, denominators, antigen invariance |
| `stage08/` | 08 dual-antigen escape | `s08a`–`s08e` — patient antigen states, null/truncate/denominator/ambient, per-cell states + technical-zero floor, evidence states, checkpoint |
| `stage09/` | 09 bulk validation | `s09a`–`s09d` — bulk abundance, pairing, marginal correlations, repeated samples + normal marrow |
| `stage09b/` | 09b provisional tiers | `s09b1`–`s09b4` — evidence matrix, tiers, cohort diagnostic, the `final → provisional` relabel |
| `stage10/` | 10 DN coherence | `s10a`–`s10g` — Level-1 structure + per-patient DE, state calls, program scores, hypotheses + evidence levels, pseudobulk DE + decoupler + TC, decoupler full space, Level-2 interpretation + freeze |

> **Naming trap.** `stage08/s08c_*.py` is **Stage-08 step c**. It is unrelated to
> **notebook `08c`**, the supplemental multi-antigen coverage arm, which is committed at
> `notebooks/08c_multi_antigen_coverage.py` and is not part of this recovery.

Stages **11**, **11b**, **08c** and **12** are not recovered here: their producers were
committed as executable notebooks from the start
(`notebooks/11_immune_context.py`, `notebooks/11b_liana_verification.py`,
`notebooks/08c_multi_antigen_coverage.py`, `notebooks/12_final_synthesis.py`). The freeze
manifest records the artifacts for 11/11b/08c and points at those notebooks; Stage 12
(run 2026-08-27, after this recovery pass) consumed the manifest as an input rather than
contributing rows to it — see `provenance/README.md` for what that means for Stage 12's
own outputs.

## Reading order

Number order is execution order, exactly as elsewhere in this project:

```
06 → 07 → 08 → 09 → 09b → 10          (then 11, 11b, 08c from notebooks/)
```

Within a stage, the letter suffix is the run order. Several steps hand state to the next
through a pickle in `/tmp` (`/tmp/s08.pkl`, `/tmp/s10c.pkl`, …); those files are gone, so
the chain must be run start-to-finish rather than resumed midway.

## Rules for this directory

1. **Do not run these to "check" a number.** Every one of them writes into `results/`, so
   running any of them overwrites a frozen artifact and destroys the freeze record. The
   frozen state is authenticated by `../provenance/frozen_artifacts_pre_stage12.sha256`;
   verify against that instead.
2. **Do not edit them to improve them.** They are historical evidence. Their value is that
   they are what actually ran. A better implementation belongs in `src/mm_escape/`, after
   Stage 12, with a documented equivalence check.
3. **Do not treat them as the narrative.** The notebooks remain the analysis a reader steps
   through. These are the derivation a reviewer audits.
4. If a genuine reproducibility defect is ever found, the fix is a **new** numbered arm with
   its own namespace and its own predeclaration — never an in-place edit here.

## Known limitation, stated plainly

The scratchpad copies themselves are gone; the transcripts are the surviving record. The
byte-identical control above, the timestamp correspondence and the clean parse are strong
evidence that the recovery is faithful, but they are evidence rather than a cryptographic
guarantee — no hash of the original scratchpad files was ever taken. This is recorded as
a residual caveat in `../docs/pre_stage12_audit_remediation.md` rather than papered over.
