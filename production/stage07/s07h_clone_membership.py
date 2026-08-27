# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 07
# Step         : s07h_clone_membership.py
# What it does : dominant-V rule, per-cell clone-state construction, V-detection depth bias
# Writes       : results/07_malignant_plasma/v_clone_membership/patient_v_evaluability.csv; patient_dominant_v_rule.csv; clone_membership_per_cell.csv.gz; clone_membership_summary_by_patient.csv; v_detection_depth_bias_*.csv
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : heredoc-written script, run with `conda run -n mm-core python <path>`
#   Original path  : /tmp/claude-1000/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646/scratchpad/s06/clone.py  (session scratchpad, since deleted)
#   Executed (UTC) : 2026-08-25T15:50:49Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
import numpy as np, pandas as pd, scanpy as sc
from pathlib import Path
from mm_escape import config, malignant as M
OUT=Path('results/07_malignant_plasma/v_clone_membership'); pd.set_option('display.width',250)

g = sc.read_h5ad('results/05_integration/integrated.h5ad')
lab = pd.read_csv('results/06_annotation/per_cell_labels.csv.gz', index_col=0)
g.obs['cell_type']=lab.reindex(g.obs_names)['cell_type'].astype(str).values
pl = g[(g.obs['cell_type']=='PlasmaCell').values].copy(); del g
V=list(pl.var_names); C=pl.layers['counts']
assert 'TNFRSF17' in V and 'GPRC5D' in V   # present, but must not touch the classifier
KV=[x for x in V if x.startswith('IGKV')]; LV=[x for x in V if x.startswith('IGLV')]
KC=['IGKC']; LC=[x for x in V if x in ('IGLC1','IGLC2','IGLC3','IGLC7')]
def sub(gs):
    idx=[V.index(x) for x in gs]; X=C[:,idx]
    return np.asarray(X.todense()) if hasattr(X,'todense') else np.asarray(X)
XKV,XLV,XKC,XLC = sub(KV),sub(LV),sub(KC),sub(LC)

o = pl.obs[['sample_name','patient_id','cohort','sample_type','total_counts',
            'n_genes_by_counts','pct_counts_mt']].copy()
o['patient_id']=o['patient_id'].astype(str); o['sample_name']=o['sample_name'].astype(str)
o['kappa_umi']=XKC.sum(axis=1); o['lambda_umi']=XLC.sum(axis=1)
o['lc_class']=M.light_chain_class(o.kappa_umi, o.lambda_umi)
o['ig_umi']=o.kappa_umi+o.lambda_umi+XKV.sum(axis=1)+XLV.sum(axis=1)
o['has_LCV']=((XKV>0).any(axis=1))|((XLV>0).any(axis=1))

# ---------- patient level ----------
rows=[]
for p, gg in o.groupby('patient_id', observed=True):
    m=o.index.get_indexer(gg.index)
    called=gg[gg.lc_class.isin(['kappa','lambda'])]
    if len(called)>=20:
        fk=(called.lc_class=='kappa').mean(); D=max(fk,1-fk)
        dom='kappa' if fk>=0.5 else 'lambda'
        cl = 'CLONAL_STRONG' if D>=0.95 else 'CLONAL_WEAK' if D>=0.85 else 'NO_RESTRICTION'
    else:
        D=np.nan; dom=None; cl='NOT_EVALUABLE'
    X,gs = (XKV,KV) if dom=='kappa' else (XLV,LV) if dom=='lambda' else (XKV,KV)
    Dt=X[m]>0; npos=int(Dt.any(axis=1).sum())
    if npos:
        cnt=Dt.sum(axis=0); i=int(np.argmax(cnt)); gene=gs[i]; frac=float(cnt[i]/npos)
        gi=V.index(gene); col=sub([gene]).ravel()
        own=(o.patient_id==p).values
        fo=float((col[own]>0).mean()); fx=float((col[~own]>0).mean())
        enr = fo/fx if fx>0 else np.inf
    else:
        gene=None; frac=np.nan; enr=np.nan
    vs = M.patient_v_evaluability(npos, float(gg.has_LCV.mean()), frac, enr)
    rows.append({'patient':p,'sample_type':gg.sample_type.iloc[0],'cohort':gg.cohort.iloc[0],
        'n_plasma':len(gg),'n_lc_called':len(called),'D':D,'clonality_state':cl,
        'dominant_class':dom,'pct_LCV':round(float(gg.has_LCV.mean()),4),
        'n_V_positive':npos,'dominant_V':gene,'top_V_frac':frac,
        'enrichment':round(enr,2) if np.isfinite(enr) else np.inf,'v_state':vs})
P=pd.DataFrame(rows)
P.to_csv(OUT/'patient_v_evaluability.csv',index=False)
P[['patient','clonality_state','dominant_class','dominant_V','top_V_frac','enrichment','v_state']]\
 .to_csv(OUT/'patient_dominant_v_rule.csv',index=False)
print("=== PATIENT STATES ==="); print(pd.crosstab(P.clonality_state,P.v_state).to_string())
print(f"\n myeloma CLONAL_STRONG & V_EVALUABLE: "
      f"{int(((P.sample_type=='myeloma')&(P.clonality_state=='CLONAL_STRONG')&(P.v_state=='V_EVALUABLE')).sum())}"
      f" of {int((P.sample_type=='myeloma').sum())}")

# ---------- cell level ----------
pmap=P.set_index('patient')
dom_v_det=np.zeros(len(o),bool); alt_v_det=np.zeros(len(o),bool)
for p, gg in o.groupby('patient_id', observed=True):
    m=o.index.get_indexer(gg.index); r=pmap.loc[p]
    if not isinstance(r.dominant_V,str): continue
    X,gs=(XKV,KV) if r.dominant_class=='kappa' else (XLV,LV)
    j=gs.index(r.dominant_V)
    dom_v_det[m]=X[m][:,j]>0
    other=np.delete(X[m],j,axis=1)
    alt_v_det[m]=(other>=config.ALT_V_MIN_UMI).any(axis=1)
o['dominant_v_detected']=dom_v_det; o['alt_v_detected']=alt_v_det
o['clone_state']=[M.clone_membership(pmap.loc[p,'clonality_state'],pmap.loc[p,'v_state'],
                    lc, pmap.loc[p,'dominant_class'] or '', dv, av)
                  for p,lc,dv,av in zip(o.patient_id,o.lc_class,o.dominant_v_detected,o.alt_v_detected)]
o.to_csv(OUT/'clone_membership_per_cell.csv.gz',compression='gzip')
print("\n=== CELL STATES ==="); print(o.clone_state.value_counts().to_string())
print("\n  myeloma only:"); print(o[o.sample_type=='myeloma'].clone_state.value_counts().to_string())

by=o.groupby(['patient_id'],observed=True).clone_state.value_counts().unstack(fill_value=0)
by=by.join(P.set_index('patient')[['cohort','sample_type','clonality_state','v_state','dominant_V','n_plasma']])
by.to_csv(OUT/'clone_membership_summary_by_patient.csv')

# ---------- 7/8. informative missingness ----------
cl_pat=P[(P.clonality_state=='CLONAL_STRONG')&(P.v_state=='V_EVALUABLE')].patient.tolist()
q=o[o.patient_id.isin(cl_pat)&o.clone_state.isin([M.CLONE_SUPPORTED,M.CLONE_COMPATIBLE_V_UNOBSERVED])].copy()
METR=['total_counts','n_genes_by_counts','ig_umi','pct_counts_mt']
pooled=q.groupby('clone_state')[METR].median()
pooled.loc['ratio_supported_over_unobserved']=(pooled.loc[M.CLONE_SUPPORTED]/pooled.loc[M.CLONE_COMPATIBLE_V_UNOBSERVED])
pooled.round(3).to_csv(OUT/'v_detection_depth_bias_pooled.csv')
print("\n=== 7. DEPTH BIAS, POOLED (median) ==="); print(pooled.round(2).to_string())

coh=q.groupby(['cohort','clone_state'],observed=True)[METR].median()
coh['n']=q.groupby(['cohort','clone_state'],observed=True).size()
coh.round(2).to_csv(OUT/'v_detection_depth_bias_by_cohort.csv')
print("\n=== 8. DEPTH BIAS BY COHORT (median) ==="); print(coh.round(1).to_string())

pp=[]
for p,gg in q.groupby('patient_id',observed=True):
    a=gg[gg.clone_state==M.CLONE_SUPPORTED]; b=gg[gg.clone_state==M.CLONE_COMPATIBLE_V_UNOBSERVED]
    if len(a)<10 or len(b)<10: continue
    pp.append({'patient':p,'cohort':gg.cohort.iloc[0],'n_supported':len(a),'n_unobserved':len(b),
      'med_umi_sup':a.total_counts.median(),'med_umi_unobs':b.total_counts.median(),
      'umi_ratio':a.total_counts.median()/max(b.total_counts.median(),1),
      'med_genes_sup':a.n_genes_by_counts.median(),'med_genes_unobs':b.n_genes_by_counts.median(),
      'genes_ratio':a.n_genes_by_counts.median()/max(b.n_genes_by_counts.median(),1)})
PP=pd.DataFrame(pp); PP.round(3).to_csv(OUT/'v_detection_depth_bias_by_patient.csv',index=False)
print(f"\n=== 9. PATIENT-LEVEL DIRECTIONALITY (n={len(PP)} patients with >=10 cells in both) ===")
print(f"  UMI ratio  (supported/unobserved): median {PP.umi_ratio.median():.2f}  "
      f"patients >1: {int((PP.umi_ratio>1).sum())}/{len(PP)}  range {PP.umi_ratio.min():.2f}-{PP.umi_ratio.max():.2f}")
print(f"  genes ratio: median {PP.genes_ratio.median():.2f}  "
      f"patients >1: {int((PP.genes_ratio>1).sum())}/{len(PP)}  range {PP.genes_ratio.min():.2f}-{PP.genes_ratio.max():.2f}")
print("\nDONE",flush=True)
