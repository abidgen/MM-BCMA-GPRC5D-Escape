# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 07
# Step         : s07j_antigen_circularity_invariance.py
# What it does : six-perturbation antigen-invariance proof for the clone call
# Writes       : results/07_malignant_plasma/v_clone_membership/antigen_circularity_invariance.{py,md}
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : heredoc-written script, run with `conda run -n mm-core python <path>`
#   Original path  : $S/invar.py  (session scratchpad, since deleted)
#   Executed (UTC) : 2026-08-25T15:57:13Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
import numpy as np, pandas as pd, scanpy as sc, hashlib, scipy.sparse as sp
from mm_escape import config, malignant as M

g = sc.read_h5ad('results/05_integration/integrated.h5ad')
lab = pd.read_csv('results/06_annotation/per_cell_labels.csv.gz', index_col=0)
g.obs['cell_type']=lab.reindex(g.obs_names)['cell_type'].astype(str).values
pl = g[(g.obs['cell_type']=='PlasmaCell').values].copy(); del g
V=list(pl.var_names); C0=sp.csc_matrix(pl.layers['counts'])
AG=[V.index('TNFRSF17'), V.index('GPRC5D')]
base_obs = pl.obs[['sample_name','patient_id','cohort','sample_type','total_counts']].copy()
base_obs['patient_id']=base_obs.patient_id.astype(str)
KV=[x for x in V if x.startswith('IGKV')]; LV=[x for x in V if x.startswith('IGLV')]
LC=[x for x in V if x in ('IGLC1','IGLC2','IGLC3','IGLC7')]

def run(C, tot):
    def sub(gs):
        X=C[:,[V.index(x) for x in gs]]
        return np.asarray(X.todense()) if hasattr(X,'todense') else np.asarray(X)
    XKV,XLV,XKC,XLC = sub(KV),sub(LV),sub(['IGKC']),sub(LC)
    o=base_obs.copy(); o['total_counts']=tot
    o['lc_class']=M.light_chain_class(XKC.sum(axis=1), XLC.sum(axis=1))
    o['has_LCV']=((XKV>0).any(axis=1))|((XLV>0).any(axis=1))
    rows=[]
    for p,gg in o.groupby('patient_id',observed=True):
        m=o.index.get_indexer(gg.index); called=gg[gg.lc_class.isin(['kappa','lambda'])]
        if len(called)>=20:
            fk=(called.lc_class=='kappa').mean(); D=max(fk,1-fk); dom='kappa' if fk>=0.5 else 'lambda'
            cl='CLONAL_STRONG' if D>=0.95 else 'CLONAL_WEAK' if D>=0.85 else 'NO_RESTRICTION'
        else: dom=None; cl='NOT_EVALUABLE'
        X,gs=(XKV,KV) if dom=='kappa' else (XLV,LV) if dom=='lambda' else (XKV,KV)
        Dt=X[m]>0; npos=int(Dt.any(axis=1).sum())
        if npos:
            cnt=Dt.sum(axis=0); i=int(np.argmax(cnt)); gene=gs[i]; frac=float(cnt[i]/npos)
            col=sub([gene]).ravel(); own=(o.patient_id==p).values
            fx=float((col[~own]>0).mean()); enr=float((col[own]>0).mean())/fx if fx>0 else np.inf
        else: gene=None; frac=np.nan; enr=np.nan
        rows.append({'patient':p,'clonality_state':cl,'dominant_class':dom,'dominant_V':gene,
            'top_V_frac':frac,'enrichment':enr,
            'v_state':M.patient_v_evaluability(npos,float(gg.has_LCV.mean()),frac,enr)})
    P=pd.DataFrame(rows).sort_values('patient').reset_index(drop=True); pmap=P.set_index('patient')
    dv=np.zeros(len(o),bool); av=np.zeros(len(o),bool)
    for p,gg in o.groupby('patient_id',observed=True):
        m=o.index.get_indexer(gg.index); r=pmap.loc[p]
        if not isinstance(r.dominant_V,str) or r.dominant_class not in ('kappa','lambda'): continue
        X,gs=(XKV,KV) if r.dominant_V in KV else (XLV,LV); j=gs.index(r.dominant_V)
        dv[m]=X[m][:,j]>0; av[m]=(np.delete(X[m],j,axis=1)>=config.ALT_V_MIN_UMI).any(axis=1)
    st=np.array([M.clone_membership(pmap.loc[p,'clonality_state'],pmap.loc[p,'v_state'],lc,
                 pmap.loc[p,'dominant_class'] or '',a,b)
                 for p,lc,a,b in zip(o.patient_id,o.lc_class,dv,av)])
    h=lambda s: hashlib.sha256(s.encode()).hexdigest()[:16]
    return {'patient_states':h(P.to_csv(index=False)),
            'v_evaluability':h(','.join(P.v_state)),
            'dominant_V':h(','.join(map(str,P.dominant_V))),
            'cell_states':h(','.join(st)), 'n_supported':int((st==M.CLONE_SUPPORTED).sum())}

tot0=base_obs.total_counts.values.astype(float)
ag_sum=np.asarray(C0[:,AG].sum(axis=1)).ravel()
def perturb(zero_idx=(), extreme_idx=()):
    C=C0.copy().tolil(); d=np.zeros(C0.shape[0])
    for i in zero_idx:
        d-=np.asarray(C0[:,i].todense()).ravel(); C[:,i]=0
    for i in extreme_idx:
        old=np.asarray(C0[:,i].todense()).ravel(); C[:,i]=10000; d+=10000-old
    return sp.csc_matrix(C), tot0+d

base=run(C0,tot0); print('BASELINE',base)
cases={
 'TNFRSF17->0':        perturb(zero_idx=[AG[0]]),
 'GPRC5D->0':          perturb(zero_idx=[AG[1]]),
 'TNFRSF17->10000':    perturb(extreme_idx=[AG[0]]),
 'GPRC5D->10000':      perturb(extreme_idx=[AG[1]]),
 'BOTH->0 (joint)':    perturb(zero_idx=AG),
 'BOTH->10000 (joint)':perturb(extreme_idx=AG),
}
ok=True
for name,(C,t) in cases.items():
    r=run(C,t); same=all(r[k]==base[k] for k in base)
    ok &= same
    print(f"{'PASS' if same else 'FAIL'}  {name:22s} "
          f"dtot_median={np.median(t-tot0):+.0f}  " +
          ' '.join(f"{k}={'=' if r[k]==base[k] else 'DIFF'}" for k in base))
print("\nANTIGEN-CIRCULARITY INVARIANCE:", "ALL PASS" if ok else "FAILED")
