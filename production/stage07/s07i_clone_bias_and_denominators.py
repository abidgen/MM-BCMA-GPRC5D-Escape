# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 07
# Step         : s07i_clone_bias_and_denominators.py
# What it does : V-detection antigen bias, repeated-sample clone consistency, primary/sensitivity denominator sizes
# Writes       : results/07_malignant_plasma/v_clone_membership/v_detection_antigen_bias_*.csv; repeated_sample_clone_consistency.csv; primary_vs_sensitivity_denominator_sizes.csv
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : heredoc-written script, run with `conda run -n mm-core python <path>`
#   Original path  : /tmp/claude-1000/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646/scratchpad/s06/clone2.py  (session scratchpad, since deleted)
#   Executed (UTC) : 2026-08-25T15:52:11Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
"""Post-freeze diagnostics. Rules are already fixed; nothing here may change them."""
import numpy as np, pandas as pd, scanpy as sc
from pathlib import Path
from mm_escape import config, malignant as M
OUT=Path('results/07_malignant_plasma/v_clone_membership'); pd.set_option('display.width',250)

o = pd.read_csv(OUT/'clone_membership_per_cell.csv.gz', index_col=0, dtype={'patient_id':str,'sample_name':str})
P = pd.read_csv(OUT/'patient_v_evaluability.csv', dtype={'patient':str})
g = sc.read_h5ad('results/05_integration/integrated.h5ad')
lab = pd.read_csv('results/06_annotation/per_cell_labels.csv.gz', index_col=0)
g.obs['cell_type']=lab.reindex(g.obs_names)['cell_type'].astype(str).values
pl = g[(g.obs['cell_type']=='PlasmaCell').values].copy(); del g
X = pl[:, ['TNFRSF17','GPRC5D']].layers['counts']
X = np.asarray(X.todense()) if hasattr(X,'todense') else np.asarray(X)
o = o.reindex(pl.obs_names)
o['BCMA_det']=X[:,0]>0; o['GPRC5D_det']=X[:,1]>0

cl=P[(P.clonality_state=='CLONAL_STRONG')&(P.v_state=='V_EVALUABLE')].patient.tolist()
q=o[o.patient_id.isin(cl)&o.clone_state.isin([M.CLONE_SUPPORTED,M.CLONE_COMPATIBLE_V_UNOBSERVED])]

print("=== 9. ANTIGEN-DETECTION BIAS (DIAGNOSTIC ONLY — cannot change any rule) ===")
pool=q.groupby('clone_state')[['BCMA_det','GPRC5D_det']].mean()
pool['n']=q.groupby('clone_state').size(); pool.round(4).to_csv(OUT/'v_detection_antigen_bias_pooled.csv')
print(pool.round(3).to_string())
coh=q.groupby(['cohort','clone_state'],observed=True)[['BCMA_det','GPRC5D_det']].mean()
coh['n']=q.groupby(['cohort','clone_state'],observed=True).size()
coh.round(4).to_csv(OUT/'v_detection_antigen_bias_by_cohort.csv')
print("\nby cohort:"); print(coh.round(3).to_string())
rows=[]
for p,gg in q.groupby('patient_id',observed=True):
    a=gg[gg.clone_state==M.CLONE_SUPPORTED]; b=gg[gg.clone_state==M.CLONE_COMPATIBLE_V_UNOBSERVED]
    if len(a)<10 or len(b)<10: continue
    rows.append({'patient':p,'cohort':gg.cohort.iloc[0],'n_sup':len(a),'n_unobs':len(b),
      'BCMA_sup':a.BCMA_det.mean(),'BCMA_unobs':b.BCMA_det.mean(),
      'GPRC5D_sup':a.GPRC5D_det.mean(),'GPRC5D_unobs':b.GPRC5D_det.mean()})
A=pd.DataFrame(rows); A.round(4).to_csv(OUT/'v_detection_antigen_bias_by_patient.csv',index=False)
print(f"\npatient-level (n={len(A)}): BCMA higher in SUPPORTED for {int((A.BCMA_sup>A.BCMA_unobs).sum())}/{len(A)}"
      f" | GPRC5D higher for {int((A.GPRC5D_sup>A.GPRC5D_unobs).sum())}/{len(A)}")
print(f"  median BCMA  {A.BCMA_sup.median():.3f} vs {A.BCMA_unobs.median():.3f}")
print(f"  median GPRC5D {A.GPRC5D_sup.median():.3f} vs {A.GPRC5D_unobs.median():.3f}")

# repeated samples
REP=['27522','47491','56203','58408','59114','60359','81012','83942']
rows=[]
for p in REP:
    gg=o[o.patient_id==p]
    if not len(gg): continue
    r=P[P.patient==p].iloc[0]
    for s,hh in gg.groupby('sample_name',observed=True):
        vc=hh.clone_state.value_counts()
        rows.append({'patient':p,'sample':s,'n_plasma':len(hh),
          'patient_clonality':r.clonality_state,'patient_v_state':r.v_state,
          'dominant_class':r.dominant_class,'dominant_V':r.dominant_V,
          'pct_lc_dominant':round(float((hh.lc_class==r.dominant_class).mean()),3),
          'pct_dominant_V_det':round(float(hh.dominant_v_detected.mean()),3),
          'frac_SUPPORTED':round(float(vc.get(M.CLONE_SUPPORTED,0)/len(hh)),3),
          'frac_V_UNOBSERVED':round(float(vc.get(M.CLONE_COMPATIBLE_V_UNOBSERVED,0)/len(hh)),3)})
R=pd.DataFrame(rows); R.to_csv(OUT/'repeated_sample_clone_consistency.csv',index=False)
print("\n=== 6. REPEATED-SAMPLE CONSISTENCY ==="); print(R.to_string(index=False))

# denominators
d=[]
for p,gg in o[o.patient_id.isin(cl)].groupby('patient_id',observed=True):
    ns=int((gg.clone_state==M.CLONE_SUPPORTED).sum()); nu=int((gg.clone_state==M.CLONE_COMPATIBLE_V_UNOBSERVED).sum())
    d.append({'patient':p,'cohort':gg.cohort.iloc[0],'n_plasma':len(gg),
              'primary_CLONE_SUPPORTED':ns,'sensitivity_broader':ns+nu,
              'ratio':round((ns+nu)/max(ns,1),3)})
D=pd.DataFrame(d).sort_values('primary_CLONE_SUPPORTED',ascending=False)
D.to_csv(OUT/'primary_vs_sensitivity_denominator_sizes.csv',index=False)
print(f"\n=== 11. DENOMINATORS ({len(D)} patients) ===")
print(f"  primary CLONE_SUPPORTED total     : {D.primary_CLONE_SUPPORTED.sum():,}")
print(f"  sensitivity broader total         : {D.sensitivity_broader.sum():,}")
print(f"  patients >=100 primary            : {int((D.primary_CLONE_SUPPORTED>=100).sum())}/{len(D)}")
print(f"  patients >=100 sensitivity        : {int((D.sensitivity_broader>=100).sum())}/{len(D)}")
print(f"  median expansion factor           : {D.ratio.median():.2f}x  (range {D.ratio.min():.2f}-{D.ratio.max():.2f})")
print(D.head(8).to_string(index=False))
print("\nDONE",flush=True)
