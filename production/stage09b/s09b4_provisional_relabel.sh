# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 09b
# Step         : s09b4_provisional_relabel.sh
# What it does : terminological relabel final -> provisional (rename + README only; no number changed)
# Writes       : results/08_dual_antigen_escape/risk_tier_provisional/ (directory and file renames, README_PROVISIONAL.md)
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : inline shell/python block
#   Original path  : (inline heredoc executed directly from the shell)
#   Executed (UTC) : 2026-08-25T22:25:05Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
cd /media/wrath/CART_mm_dual_antigen
D=results/08_dual_antigen_escape
git mv 2>/dev/null; mv $D/risk_tier_final $D/risk_tier_provisional
P=$D/risk_tier_provisional
mv $P/risk_tiers_final.csv $P/risk_tiers_provisional.csv
mv $P/risk_tier_final_summary.md $P/risk_tier_provisional_summary.md
HDR='> **PROVISIONAL — Stage-08/09b measurement-robust classification only. Stage 10
> biological coherence assessment is required before final biological interpretation.**
>
> These tiers say a DN measurement survived the technical and threshold sensitivity
> analyses. They do **not** say the double-negative population is biologically coherent.
> Stage 10 sits between provisional measurement-robust tiering and any final biological
> risk classification.

'
for f in $P/*.md; do printf '%s' "$HDR" | cat - "$f" > /tmp/h && mv /tmp/h "$f"; done
# CSVs carry the marker in the filename; a README covers the directory
cat > $P/README_PROVISIONAL.md << 'MDEOF'
# PROVISIONAL — read before using anything in this directory

> **PROVISIONAL — Stage-08/09b measurement-robust classification only. Stage 10
> biological coherence assessment is required before final biological interpretation.**

Every file here describes a **measurement-robust provisional tier under the frozen
Stage-08/09b rule**. A `robust-high` label means the observed double-negative fraction
survived the denominator, depth, repeated-sample, null-scheme and threshold sensitivity
analyses. **It does not mean the DN population is a biologically coherent escape
state**, and it never meant a genetic subclone.

**Stage 10 sits between provisional measurement-robust tiering and any final biological
risk classification.** Until Stage 10's coherence states are frozen, no file here may be
cited as a final biological escape classification.

Patient membership, thresholds, the tier algorithm and every Stage-08/09 measurement are
**unchanged** by this relabelling — it is provenance and terminology only.

| file | holds |
|---|---|
| `risk_tiers_provisional.csv` | per-patient provisional tier, all evidence columns, machine-readable uncertainty reasons |
| `risk_tier_provisional_summary.md` | the written result and its limitations |
| `risk_tier_policy.md` | the frozen rule, written before application |
| `risk_tiers_tau{020,025,033}.csv` | the rule applied at each threshold |
| `risk_tier_threshold_comparison.csv` | tier under each threshold + threshold-robustness state |
| `risk_tier_by_cohort_diagnostic.csv` | post-assignment cohort diagnostic (diagnostic only) |
| `risk_tier_evidence_long.csv` | long-form uncertainty reasons |
MDEOF
sed -i 's|risk_tier_final|risk_tier_provisional|g; s|risk_tiers_final\.csv|risk_tiers_provisional.csv|g; s|risk_tier_final_summary\.md|risk_tier_provisional_summary.md|g' notebooks/09b_risk_tiers.py
sed -i 's|^# # 09b — Final risk tiers|# # 09b — Provisional risk tiers (measurement-robust only)|' notebooks/09b_risk_tiers.py
grep -n "risk_tier_provisional\|Provisional risk tiers" notebooks/09b_risk_tiers.py | head -5
ls $P
