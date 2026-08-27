# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 09b
# Step         : s09b1_patient_evidence_matrix.py
# What it does : per-patient evidence matrix from Stage-07/08 evidence only
# Writes       : results/08_dual_antigen_escape/risk_tier_design/patient_evidence_matrix.csv
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : heredoc-written script, run with `conda run -n mm-core python <path>`
#   Original path  : $S/mat.py  (session scratchpad, since deleted)
#   Executed (UTC) : 2026-08-25T16:25:19Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
import numpy as np, pandas as pd, pickle
from pathlib import Path
from mm_escape import antigen as A
OUT=Path('results/08_dual_antigen_escape'); RT=OUT/'risk_tier_design'
S=pickle.load(open('/tmp/s08.pkl','rb')); SEED=20260825
prim,sens,strata_c,strata_g=S['prim'],S['sens'],S['strata_c'],S['strata_g']
bc,gc,tot,D,elig,obs,EDGES=S['bc'],S['gc'],S['tot'],S['D'],S['elig'],S['obs'],S['EDGES']
pd.set_option('display.width',300)

# truncated counts / strata (same frozen procedure)
bt=A.downsample_gene_counts(tot,bc,10000,SEED); gt=A.downsample_gene_counts(tot,gc,10000,SEED)
Dt=A.depth_ex_antigen(np.minimum(tot,10000),bt,gt)
st_t=np.zeros(len(obs),np.int64)
for ch in sorted(obs.cohort[prim].unique()):
    m=(obs.cohort==ch).values
    st_t[m]=A.assign_strata(Dt[m],A.quantile_edges(Dt[prim&m],5))

def core(mask,strata,b_,g_,boot=False):
    r={}
    for p in sorted(elig):
        m=mask&(obs.patient_id==p).values
        if m.sum()==0: continue
        b0=b_[m]==0; g0=g_[m]==0; st=A.merge_sparse_strata(strata[m],20)
        o=float((b0&g0).mean()); e=A.stratified_expected_dn(b0,g0,st)
        d={'n':int(m.sum()),'n_samples':int(len(np.unique(obs.sample_name.values[m]))),
           'bcma_neg':float(b0.mean()),'gprc5d_neg':float(g0.mean()),'obs_dn':o,
           'exp_dn':e,'enr':o/e if e>0 else np.nan,'excess':o-e}
        if boot:
            sid=obs.sample_name.values[m]; idx=np.arange(m.sum())
            db=A.hierarchical_bootstrap(sid,idx,2000,SEED,lambda i:float((b0[i]&g0[i]).mean()))
            def rr(i):
                ee=A.stratified_expected_dn(b0[i],g0[i],st[i])
                return float((b0[i]&g0[i]).mean()/ee) if ee>0 else np.nan
            eb=A.hierarchical_bootstrap(sid,idx,2000,SEED,rr)
            d['dn_lo'],d['dn_hi']=np.percentile(db,[2.5,97.5])
            d['enr_lo'],d['enr_hi']=np.nanpercentile(eb,[2.5,97.5])
            _,pv=A.permutation_null_dn(b0,g0,st,2000,SEED); d['perm_p']=pv
        r[p]=d
    return r

P=core(prim,strata_c,bc,gc,True);  Sn=core(sens,strata_c,bc,gc,True)
Pg=core(prim,strata_g,bc,gc);      Sg=core(sens,strata_g,bc,gc)
Pt=core(prim,st_t,bt,gt);          St=core(sens,st_t,bt,gt)

# per-sample estimates (>=20 cells), both denominators
def per_sample(mask,b_,g_):
    out={}
    for p in sorted(elig):
        m=mask&(obs.patient_id==p).values
        est=[]
        for s in sorted(obs.sample_name[m].unique()):
            q=m&(obs.sample_name==s).values
            if q.sum()>=20:
                est.append((s,int(q.sum()),float(((b_[q]==0)&(g_[q]==0)).mean())))
        out[p]=est
    return out
PS=per_sample(prim,bc,gc); SS=per_sample(sens,bc,gc)

rows=[]
for p in sorted(elig):
    a,b=P[p],Sn.get(p,{})
    ps=PS[p]; ss=SS[p]
    vals=[v for _,_,v in ps]
    r={'patient':p,'cohort':obs.cohort[(obs.patient_id==p).values].iloc[0]}
    for tag,src in [('primary',a),('sensitivity',b)]:
        for k in ['n','n_samples','bcma_neg','gprc5d_neg','obs_dn','exp_dn','enr','excess',
                  'dn_lo','dn_hi','enr_lo','enr_hi','perm_p']:
            r[f'{k}_{tag}']=src.get(k,np.nan)
    r['obs_dn_trunc10k_primary']=Pt[p]['obs_dn']; r['obs_dn_trunc10k_sensitivity']=St[p]['obs_dn']
    r['delta_trunc10k_primary']=Pt[p]['obs_dn']-a['obs_dn']
    r['delta_trunc10k_sensitivity']=St[p]['obs_dn']-b.get('obs_dn',np.nan)
    r['delta_primary_vs_sensitivity']=b.get('obs_dn',np.nan)-a['obs_dn']
    r['enr_globalbins_primary']=Pg[p]['enr']; r['enr_globalbins_sensitivity']=Sg[p]['enr']
    r['delta_nullscheme_primary']=Pg[p]['enr']-a['enr']
    r['n_samples_assessable_primary']=len(ps)
    r['per_sample_dn_primary']='|'.join(f'{s}:{n}:{v:.3f}' for s,n,v in ps)
    r['repeated_dn_range_primary']=(max(vals)-min(vals)) if len(vals)>1 else np.nan
    r['repeated_dn_min_primary']=min(vals) if vals else np.nan
    r['repeated_dn_max_primary']=max(vals) if vals else np.nan
    sv=[v for _,_,v in ss]
    r['repeated_dn_range_sensitivity']=(max(sv)-min(sv)) if len(sv)>1 else np.nan
    r['ci_width_primary']=a.get('dn_hi',np.nan)-a.get('dn_lo',np.nan)
    r['ci_width_sensitivity']=b.get('dn_hi',np.nan)-b.get('dn_lo',np.nan)
    r['flag_low_n']=bool(a['n']<100)
    r['flag_single_sample']=bool(a['n_samples']==1)
    r['flag_denominator_unstable']=bool(abs(r['delta_primary_vs_sensitivity'])>0.05)
    r['flag_depth_sensitive']=bool(abs(r['delta_trunc10k_primary'])>0.05)
    r['flag_null_scheme_sensitive']=bool((Pg[p]['enr']-1)*(a['enr']-1)<0)
    rows.append(r)
M=pd.DataFrame(rows)
M.round(5).to_csv(RT/'patient_evidence_matrix.csv',index=False)
print(f"evidence matrix: {M.shape[0]} patients x {M.shape[1]} columns")
show=['patient','cohort','n_primary','n_samples_primary','obs_dn_primary','dn_lo_primary','dn_hi_primary',
 'enr_primary','enr_lo_primary','enr_hi_primary','perm_p_primary','obs_dn_sensitivity',
 'obs_dn_trunc10k_primary','repeated_dn_range_primary','ci_width_primary']
print(M[show].round(3).to_string(index=False))
print("\nCI width (primary): median %.3f  min %.3f  max %.3f"%(M.ci_width_primary.median(),M.ci_width_primary.min(),M.ci_width_primary.max()))
print("enr_lo>1 (primary):",int((M.enr_lo_primary>1).sum()),"| enr_lo>1 (sensitivity):",int((M.enr_lo_sensitivity>1).sum()))
print("patients with >1 assessable sample:",int((M.n_samples_assessable_primary>1).sum()))
print("\nflag counts:"); print(M[[c for c in M if c.startswith('flag_')]].sum().to_string())
