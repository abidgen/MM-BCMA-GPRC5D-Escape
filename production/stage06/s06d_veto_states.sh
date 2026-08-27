# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 06
# Step         : s06d_veto_states.sh
# What it does : per-class veto states
# Writes       : results/06_annotation/veto_states.csv
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : inline shell/python block
#   Original path  : (inline heredoc executed directly from the shell)
#   Executed (UTC) : 2026-08-25T05:58:58Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
cd /media/wrath/CART_mm_dual_antigen
sed -i 's|V_PREV = REPO / "results" / "06_annotation_v2" / "annotated.h5ad"|V_PREV = REPO / "results" / "06_annotation_v3" / "annotated.h5ad"|; s|print(f"reusing v2 automated predictions|print(f"reusing prior automated predictions|' notebooks/06_annotation.py
python3 - <<'PY'
p='notebooks/06_annotation.py'
s=open(p).read()
s=s.replace('''print(decision[["chosen_method", "reason", "f1_threshold", "f1_chosen",
                "coverage_chosen", "contradiction_chosen",
                "vetoed_by_coverage", "vetoed_by_contradiction"]].to_string())''',
'''print(decision[["chosen_method", "reason", "f1_chosen", "coverage_chosen",
                "contradiction_chosen", "vetoed_by_coverage",
                "vetoed_by_contradiction", "not_evaluable"]].to_string())

# Explicit NOT_EVALUABLE audit across all seven classes — no combination may be
# silently read as a pass or a veto.
state_cols = [c for c in decision.columns if c.startswith(("cov_state_", "contra_state_"))]
states = decision[state_cols]
states.to_csv(OUT / "veto_states.csv")
print("\\n--- NOT_EVALUABLE combinations (all seven classes) ---")
ne = [(cls, c) for cls in decision.index for c in state_cols
      if decision.loc[cls, c] == ann.NOT_EVALUABLE]
for cls, c in ne:
    print(f"  {cls:11s} {c}")
print(f"  total: {len(ne)}")''',1)
open(p,'w').write(s)
PY
/home/abid/miniforge3/envs/mm-annotation/bin/python -c "import ast;ast.parse(open('notebooks/06_annotation.py').read());print('syntax OK')"
/home/abid/miniforge3/envs/mm-annotation/bin/python -u notebooks/06_annotation.py > /tmp/claude-1000/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646/scratchpad/s06/run06_v4.log 2>&1
