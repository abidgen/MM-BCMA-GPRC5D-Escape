# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 08
# Step         : s08c_cell_states_repeated_technical_zero.py
# What it does : per-cell antigen states, repeated-sample consistency, expression-matched technical-zero floor
# Writes       : results/08_dual_antigen_escape/cell_antigen_states.csv.gz; repeated_sample_antigen_consistency.csv; noise_floor_technical_zero.csv
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : heredoc-written script, run with `conda run -n mm-core python <path>`
#   Original path  : $S/s08c.py  (session scratchpad, since deleted)
#   Executed (UTC) : 2026-08-25T16:08:46Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
import numpy as np, pandas as pd, pickle, scanpy as sc, scipy.sparse as sp
from pathlib import Path
from mm_escape import antigen as A
OUT=Path('results/08_dual_antigen_escape'); pd.set_option('display.width',250)
S=pickle.load(open('/tmp/s08.pkl','rb')); B=pickle.load(open('/tmp/s08b.pkl','rb'))
prim,sens,strata_c=S['prim'],S['sens'],S['strata_c']
bc,gc,tot,D,elig,obs=S['bc'],S['gc'],S['tot'],S['D'],S['elig'],S['obs']
npm=B['npm']; SEED=20260825

# ---- cell-level artifact ----
cell=obs.loc[prim|sens,['sample_name','patient_id','cohort','clone_state','bcma','gprc5d',
                        'total_umi','depth_ex_antigen']].copy()
cell['in_primary']=prim[prim|sens]; cell['in_sensitivity']=True
cell['depth_stratum_cohort']=strata_c[prim|sens]
cell['observed_state']=A.observed_states(cell.bcma.values>0,cell.gprc5d.values>0)
cell['observed_state_umi2']=A.observed_states(cell.bcma.values>=2,cell.gprc5d.values>=2)
cell.to_csv(OUT/'cell_antigen_states.csv.gz',compression='gzip')
print(f"cell-level rows {len(cell)}")
print("\nobserved states, PRIMARY denominator (>0 rule):")
print((cell[cell.in_primary].observed_state.value_counts(normalize=True)*100).round(2).to_string())
print("declared sensitivity rule (>=2 UMI), primary denominator:")
print((cell[cell.in_primary].observed_state_umi2.value_counts(normalize=True)*100).round(2).to_string())

# ---- 6. repeated-sample consistency ----
REP=['27522','47491','56203','58408','59114','60359','81012','83942']
r=[]
for p in REP:
    m=prim&(obs.patient_id==p).values
    if m.sum()==0:
        r.append({'patient':p,'sample':'--','n_cells':0,'note':'no primary denominator cells'}); continue
    for s in sorted(obs.sample_name[m].unique()):
        q=m&(obs.sample_name==s).values
        if q.sum()==0: continue
        r.append({'patient':p,'sample':s,'n_cells':int(q.sum()),
            'median_depth':float(np.median(D[q])),
            'bcma_detect':round(float((bc[q]>0).mean()),4),
            'gprc5d_detect':round(float((gc[q]>0).mean()),4),
            'obs_dn':round(float(((bc[q]==0)&(gc[q]==0)).mean()),4),'note':''})
RS=pd.DataFrame(r); RS.to_csv(OUT/'repeated_sample_antigen_consistency.csv',index=False)
print("\n=== 10. REPEATED-SAMPLE CONSISTENCY (primary) ==="); print(RS.to_string(index=False))
w=RS[RS.n_cells>=20].groupby('patient').agg(n_samples=('sample','size'),
    dn_min=('obs_dn','min'),dn_max=('obs_dn','max'),depth_min=('median_depth','min'),
    depth_max=('median_depth','max'))
w['dn_spread']=w.dn_max-w.dn_min; w['depth_ratio']=w.depth_max/w.depth_min
print("\n  within-patient spread across samples:"); print(w.round(3).to_string())

# ---- 8. expression-matched technical-zero floor ----
g=sc.read_h5ad('results/05_integration/integrated.h5ad')
C=sp.csc_matrix(g.layers['counts']); V=np.array(g.var_names); del g
mu_den=np.asarray(C[prim].mean(axis=0)).ravel()
mu_np=np.asarray(C[npm].mean(axis=0)).ravel()
iB,iG=int(np.where(V=='TNFRSF17')[0][0]),int(np.where(V=='GPRC5D')[0][0])
bad=np.array([str(x).startswith(('IGK','IGL','IGH')) for x in V])
bad[[iB,iG]]=True
rows=[]
for nm,i in [('TNFRSF17',iB),('GPRC5D',iG)]:
    target=mu_den[i]
    cand=np.flatnonzero((~bad)&(mu_np>0))
    sel=cand[np.argsort(np.abs(mu_np[cand]-target))[:100]]
    X=C[:,sel]
    for ch in ['MMRF','WU1','WU2']:
        e=S['EDGES'][ch]; sa=A.assign_strata(D,e)
        for k in range(len(e)-1):
            m=npm&(obs.cohort==ch).values&(sa==k)
            if m.sum()<50: continue
            z=float((np.asarray(X[m].todense())==0).mean())
            rows.append({'antigen':nm,'target_mean_in_denominator':round(float(target),4),
                'cohort':ch,'stratum':k,'n_reference_cells':int(m.sum()),
                'median_depth':float(np.median(D[m])),
                'matched_control_genes':len(sel),
                'technical_zero_fraction':round(z,4)})
TZ=pd.DataFrame(rows); TZ.to_csv(OUT/'noise_floor_technical_zero.csv',index=False)
print("\n=== TECHNICAL-ZERO FLOOR (expression-matched control genes, non-plasma cells) ===")
print(TZ.to_string(index=False))
print("\nSTEP3 DONE",flush=True)
