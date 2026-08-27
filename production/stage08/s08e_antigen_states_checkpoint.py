# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 08
# Step         : s08e_antigen_states_checkpoint.py
# What it does : antigen-state AnnData checkpoint
# Writes       : results/08_dual_antigen_escape/antigen_states.h5ad
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : heredoc-written script, run with `conda run -n mm-core python <path>`
#   Original path  : $S/s08e.py  (session scratchpad, since deleted)
#   Executed (UTC) : 2026-08-25T16:09:51Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
import numpy as np, pandas as pd, scanpy as sc, pickle
from pathlib import Path
OUT=Path('results/08_dual_antigen_escape')
S=pickle.load(open('/tmp/s08.pkl','rb')); prim,sens=S['prim'],S['sens']
cell=pd.read_csv(OUT/'cell_antigen_states.csv.gz',index_col=0)
g=sc.read_h5ad('results/05_integration/integrated.h5ad')
sub=g[(prim|sens)].copy(); del g
for c in ['clone_state','in_primary','in_sensitivity','depth_stratum_cohort',
          'observed_state','observed_state_umi2','bcma','gprc5d','depth_ex_antigen','total_umi']:
    sub.obs[c]=cell.reindex(sub.obs_names)[c].values
for k in list(sub.obsm): del sub.obsm[k]
sub.write(OUT/'antigen_states.h5ad')
print('checkpoint', sub.shape, 'primary', int(sub.obs.in_primary.sum()))
