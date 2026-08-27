# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 10
# Step         : s10g_level2_interpretation_and_freeze.sh
# What it does : adds level2_interpretation column and writes the Stage-10 freeze README
# Writes       : results/10_dn_coherence/stage10_evidence_levels.csv (level2_interpretation column); README_STAGE10_FROZEN.md
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : inline shell/python block
#   Original path  : (inline heredoc executed directly from the shell)
#   Executed (UTC) : 2026-08-26T00:42:30Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
cd /media/wrath/CART_mm_dual_antigen
# Level-2 caveat carried on the artifact itself
/home/abid/miniforge3/envs/mm-core/bin/python - << 'PY'
import pandas as pd
p='results/10_dn_coherence/stage10_evidence_levels.csv'
d=pd.read_csv(p,dtype={'patient':str})
d['level2_interpretation']=d.level2_state.map({
 'DN_STATE_SUPPORTED':'compatible with the cohort-level DN-associated program; NOT strong '
   'patient-specific evidence of a distinct escape state (26/27 evaluable patients qualify)',
 'DN_STATE_NOT_SUPPORTED':'evaluable, no reproducible program in the same direction',
 'DN_STATE_NOT_EVALUABLE':'insufficient DN or comparator cells; NOT negative evidence'})
d.to_csv(p,index=False); print('evidence-levels columns:',list(d.columns))
PY
cat > results/10_dn_coherence/README_STAGE10_FROZEN.md << 'MDEOF'
# Stage 10 — COMPLETE and FROZEN (2026-08-25)

Three **independent** evidence axes. They are never combined into a scalar score.

| level | file(s) | states | licenses |
|---|---|---|---|
| **1 structure** | `dn_coherence_final_states.csv`, `dn_local_structure_by_patient.csv` | 4 / 23 / 5 | *non-random DN organization* |
| **2 state** | `dn_program_scores_by_patient.csv`, `level2_program_cohort_tests.csv`, `pseudobulk_de_results.csv` | 26 / 1 / 5 | *escape-associated transcriptional state* |
| **3 genomic** | — | `CNV_SUBCLONE_NOT_EVALUABLE` ×32 | nothing |

`stage10_evidence_levels.csv` is the joined view and carries `level2_interpretation`.

## Compatibility alias — documented, not silently changed

`dn_coherence_final_states.csv` is the **original Level-1 raw table** and keeps its
original schema and column name `dn_coherence_state`. Its values map to the corrected
vocabulary as:

| original value | Level-1 vocabulary |
|---|---|
| `DN_COHERENCE_SUPPORTED` | `DN_STRUCTURE_SUPPORTED` |
| `DN_COHERENCE_NOT_SUPPORTED` | `DN_STRUCTURE_NOT_SUPPORTED` |
| `DN_COHERENCE_NOT_EVALUABLE` | `DN_STRUCTURE_NOT_EVALUABLE` |

The original file was **not** rewritten, so provenance for the Level-1 run is intact.
`subclone.level1_structure_state()` performs the mapping in code.

## Read before citing any of it

- **Level 1 alone is not an escape-associated state.**
- **`DN_STATE_SUPPORTED` is weakly discriminative** — 26 of 27 evaluable patients qualify.
  It means *compatible with the cohort-level program*, not patient-specific evidence.
- **No patient is a genetic subclone**; `CNV_SUBCLONE_NOT_SUPPORTED` is never emitted.
- **No patient is simultaneously measurement-robust-high, Level-1 supported and Level-2
  supported.**
- The **pre-registered γ-secretase hypothesis is a clean negative**, in the direction
  opposite to prediction.
- The DN phenotype is **compatible with a less secretory / less differentiated plasma-cell
  state**; this data cannot separate that from an antigen-specific escape mechanism.
MDEOF
echo "--- jupytext note into CLAUDE.md:"
python3 - << 'PY'
from pathlib import Path
p=Path('CLAUDE.md'); t=p.read_text()
old="Notebooks are **paired with `jupytext`** (percent format); the `.ipynb` is gitignored\nand generated from the committed `.py` (see `.gitignore` for how to flip that)."
new=("Notebooks are **paired with `jupytext`** (percent format); the `.ipynb` is gitignored\n"
 "and generated from the committed `.py` (see `.gitignore` for how to flip that).\n"
 "**Do NOT use `jupytext --sync` in this repo.** It has repeatedly treated the `.ipynb` as\n"
 "the authoritative representation and silently dropped cells appended to the `.py` — it did\n"
 "so three times in one session (notebooks 08, 09b, 10), each time leaving a `.py` that\n"
 "compiled and a `.ipynb` missing the new work. **The `.py` is authoritative; always convert\n"
 "explicitly in that direction:** `jupytext --to notebook notebooks/NN_x.py -o notebooks/NN_x.ipynb`,\n"
 "then execute and verify the expected counts appear in the output. Never adjust scientific\n"
 "output to make a notebook reproduce.")
assert old in t; p.write_text(t.replace(old,new,1)); print('ok',len(t.replace(old,new,1)))
PY
