# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 09b
# Step         : s09b3_cohort_diagnostic.py
# What it does : tier-by-cohort diagnostic
# Writes       : results/08_dual_antigen_escape/risk_tier_provisional/risk_tier_by_cohort_diagnostic.csv
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/9ccb4002-678a-444f-b2e7-87985ef626d2.jsonl
#   Original form  : heredoc-written script, run with `conda run -n mm-core python <path>`
#   Original path  : $SP/diag.py  (session scratchpad, since deleted)
#   Executed (UTC) : 2026-08-26T22:08:20Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
import sys; sys.path.insert(0,'src')
import numpy as np, pandas as pd
from mm_escape import communication as CM
L=pd.read_csv("results/06_annotation/per_cell_labels.csv.gz",dtype={"patient_id":str})
R=pd.read_csv("results/06_annotation/cluster23_local/trbc_context_revision/revised_lineage_calls.csv.gz",usecols=["cell_id","call"])
C23={"NK":"NK_core","T_NK_mixed":"cytotoxic_mixed","unresolved":"cytotoxic_mixed","T_ab":"Tcell","T_gd":"Tcell"}
L=L.merge(R,on="cell_id",how="left"); L["lineage"]=np.where(L.call.notna(),L.call.map(C23),L.cell_type)
A=pd.read_csv("results/08_dual_antigen_escape/cell_antigen_states.csv.gz",dtype={"patient_id":str})
prim=set(A.loc[A.in_primary,"cell_id"])
P=pd.read_csv("results/08_dual_antigen_escape/patient_antigen_states_primary.csv",dtype={"patient":str})
L["role"]=np.where(L.cell_id.isin(prim),"clone_primary",
          np.where(L.lineage=="PlasmaCell","other_plasma",
          np.where(L.lineage.isin(CM.LR_SENDERS),L.lineage,"drop")))
c=L[(L.role!="drop")&L.patient_id.isin(set(P.patient))]
tab=c.groupby(["patient_id","role"]).size().unstack(fill_value=0)
ge20=(tab>=20).sum(axis=1)
print("n groups with >=20 cells, per patient:")
print(ge20.value_counts().sort_index().to_string())
bad=ge20[ge20<2]
print("\nPATIENTS WITH <2 GROUPS >=20 CELLS:")
print(tab.loc[bad.index].to_string() if len(bad) else "  none")
print("\nrows where some group is 0 (empty 'rest' risk):")
print(tab[(tab==0).any(axis=1)].to_string())
