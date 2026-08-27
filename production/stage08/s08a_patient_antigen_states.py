# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 08
# Step         : s08a_patient_antigen_states.py
# What it does : raw-count antigen calls, cohort depth strata, per-patient DN + depth-conditioned null + hierarchical bootstrap
# Writes       : results/08_dual_antigen_escape/depth_strata_definition.csv; patient_antigen_states_primary.csv; patient_antigen_states_sensitivity.csv; patient_conegativity_enrichment.csv; patient_bootstrap_intervals.csv
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : heredoc-written script, run with `conda run -n mm-core python <path>`
#   Original path  : $S/s08.py  (session scratchpad, since deleted)
#   Executed (UTC) : 2026-08-25T16:06:57Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
from pathlib import Path
from mm_escape import antigen as A, malignant as M
OUT=Path('results/08_dual_antigen_escape'); OUT.mkdir(exist_ok=True)
pd.set_option('display.width',260)
SEED=20260825

# ---------- load ----------
g=sc.read_h5ad('results/05_integration/integrated.h5ad')
lab=pd.read_csv('results/06_annotation/per_cell_labels.csv.gz',index_col=0)
g.obs['cell_type']=lab.reindex(g.obs_names)['cell_type'].astype(str).values
V=list(g.var_names); C=sp.csc_matrix(g.layers['counts'])
iB,iG=V.index('TNFRSF17'),V.index('GPRC5D')
bc=np.asarray(C[:,iB].todense()).ravel().astype(np.int64)
gc=np.asarray(C[:,iG].todense()).ravel().astype(np.int64)
tot=np.asarray(C.sum(axis=1)).ravel().astype(np.int64)
obs=g.obs[['sample_name','patient_id','cohort','sample_type','cell_type']].copy()
obs['patient_id']=obs.patient_id.astype(str); obs['sample_name']=obs.sample_name.astype(str)
obs['cohort']=obs.cohort.astype(str); obs['bcma']=bc; obs['gprc5d']=gc; obs['total_umi']=tot
obs['depth_ex_antigen']=A.depth_ex_antigen(tot,bc,gc)

# ---------- 1. denominators, frozen from stage 07 ----------
cm=pd.read_csv('results/07_malignant_plasma/v_clone_membership/clone_membership_per_cell.csv.gz',index_col=0)
obs['clone_state']=cm.reindex(obs.index)['clone_state'].astype(str).values
P7=pd.read_csv('results/07_malignant_plasma/v_clone_membership/patient_v_evaluability.csv',dtype={'patient':str})
elig=set(P7[(P7.clonality_state=='CLONAL_STRONG')&(P7.v_state=='V_EVALUABLE')].patient)
inel=obs.patient_id.isin(elig).values
prim = inel & (obs.clone_state==M.CLONE_SUPPORTED).values
sens = inel & obs.clone_state.isin([M.CLONE_SUPPORTED,M.CLONE_COMPATIBLE_V_UNOBSERVED]).values
assert not (obs.clone_state[sens].isin([M.CLONE_INCOMPATIBLE,M.CLONE_UNCERTAIN])).any()
print(f"eligible patients {len(elig)} | primary {prim.sum()} | sensitivity {sens.sum()}")

# ---------- 4. cohort-specific strata, edges from PRIMARY denominator only ----------
def build_edges(mask, depth, by_cohort=True, nb=5):
    ed={}
    if by_cohort:
        for ch in sorted(obs.cohort[mask].unique()):
            m=mask&(obs.cohort==ch).values
            ed[ch]=A.quantile_edges(depth[m],nb)
    else:
        ed['__global__']=A.quantile_edges(depth[mask],nb)
    return ed
def apply_edges(ed, depth, by_cohort=True):
    s=np.zeros(len(obs),dtype=np.int64)
    if by_cohort:
        for ch,e in ed.items():
            m=(obs.cohort==ch).values; s[m]=A.assign_strata(depth[m],e)
    else:
        s=A.assign_strata(depth,ed['__global__'])
    return s
D=obs.depth_ex_antigen.values
EDGES=build_edges(prim,D,True); GEDGES=build_edges(prim,D,False)
strata_c=apply_edges(EDGES,D,True); strata_g=apply_edges(GEDGES,D,False)
rows=[]
for ch,e in EDGES.items():
    for i in range(len(e)-1):
        m=prim&(obs.cohort==ch).values&(strata_c==i)
        rows.append({'scheme':'cohort_specific','cohort':ch,'stratum':i,
                     'lo':e[i],'hi':e[i+1],'n_primary_cells':int(m.sum())})
for i in range(len(GEDGES['__global__'])-1):
    e=GEDGES['__global__']
    rows.append({'scheme':'global_secondary','cohort':'ALL','stratum':i,
                 'lo':e[i],'hi':e[i+1],'n_primary_cells':int((prim&(strata_g==i)).sum())})
SD=pd.DataFrame(rows); SD.to_csv(OUT/'depth_strata_definition.csv',index=False)
print("\n=== 3. DEPTH STRATA (frozen, cohort-specific primary) ==="); print(SD.to_string(index=False))

# ---------- per-patient engine ----------
def patient_table(mask, strata, bcma, gprc, tag, boot=True, nperm=2000):
    res=[]
    for p in sorted(elig):
        m=mask&(obs.patient_id==p).values
        n=int(m.sum())
        if n==0:
            res.append({'patient':p,'n_cells':0}); continue
        b0=bcma[m]==0; g0=gprc[m]==0
        st=A.merge_sparse_strata(strata[m],20)
        obs_dn=float((b0&g0).mean())
        exp_s=A.stratified_expected_dn(b0,g0,st)
        exp_u=A.unconditioned_expected_dn(b0,g0)
        sids=obs.sample_name.values[m]
        r={'patient':p,'cohort':obs.cohort.values[m][0],'n_cells':n,
           'n_samples':int(len(np.unique(sids))),'n_strata':int(len(np.unique(st))),
           'bcma_detect':float((bcma[m]>0).mean()),'gprc5d_detect':float((gprc[m]>0).mean()),
           'bcma_neg':float(b0.mean()),'gprc5d_neg':float(g0.mean()),
           'obs_double_positive':float(((bcma[m]>0)&(gprc[m]>0)).mean()),
           'obs_BCMA_only':float(((bcma[m]>0)&g0).mean()),
           'obs_GPRC5D_only':float((b0&(gprc[m]>0)).mean()),
           'observed_double_negative_fraction':obs_dn,
           'expected_dn_stratified':exp_s,'expected_dn_unconditioned':exp_u,
           'enrichment_stratified':obs_dn/exp_s if exp_s>0 else np.nan,
           'enrichment_unconditioned':obs_dn/exp_u if exp_u>0 else np.nan,
           'excess_dn':obs_dn-exp_s,
           'median_depth_ex_antigen':float(np.median(D[m])),
           'low_n':bool(n<100),'single_sample':bool(len(np.unique(sids))==1)}
        if boot and n>=20:
            idx=np.arange(n)
            dnb=A.hierarchical_bootstrap(sids,idx,2000,SEED,lambda i:float((b0[i]&g0[i]).mean()))
            def enr(i):
                e=A.stratified_expected_dn(b0[i],g0[i],st[i])
                return float((b0[i]&g0[i]).mean()/e) if e>0 else np.nan
            eb=A.hierarchical_bootstrap(sids,idx,2000,SEED,enr)
            r['dn_ci_lo'],r['dn_ci_hi']=np.percentile(dnb,[2.5,97.5])
            r['enr_ci_lo'],r['enr_ci_hi']=np.nanpercentile(eb,[2.5,97.5])
            r['ci_width']=r['dn_ci_hi']-r['dn_ci_lo']
        if nperm and n>=20:
            null,pv=A.permutation_null_dn(b0,g0,st,nperm,SEED)
            r['perm_p']=pv; r['null_dn_lo'],r['null_dn_hi']=np.percentile(null,[2.5,97.5])
        res.append(r)
    return pd.DataFrame(res)

PRI=patient_table(prim,strata_c,bc,gc,'primary')
SEN=patient_table(sens,strata_c,bc,gc,'sensitivity')
PRI.to_csv(OUT/'patient_antigen_states_primary.csv',index=False)
SEN.to_csv(OUT/'patient_antigen_states_sensitivity.csv',index=False)
cols=['patient','cohort','n_cells','n_samples','bcma_detect','gprc5d_detect',
      'observed_double_negative_fraction','dn_ci_lo','dn_ci_hi','expected_dn_stratified',
      'enrichment_stratified','enr_ci_lo','enr_ci_hi','perm_p','low_n']
print("\n=== PRIMARY (CLONE_SUPPORTED) ==="); print(PRI[cols].round(4).to_string(index=False))
print("\n=== SENSITIVITY (broader clone-compatible) ==="); print(SEN[cols].round(4).to_string(index=False))
PRI.assign(scheme='cohort_specific').to_csv(OUT/'patient_conegativity_enrichment.csv',index=False)
PRI[['patient','cohort','n_cells','n_samples','observed_double_negative_fraction',
     'dn_ci_lo','dn_ci_hi','enrichment_stratified','enr_ci_lo','enr_ci_hi','low_n','single_sample']]\
  .to_csv(OUT/'patient_bootstrap_intervals.csv',index=False)
np.save('/tmp/s08_state.npy',np.array([0]))
import pickle; pickle.dump(dict(prim=prim,sens=sens,strata_c=strata_c,strata_g=strata_g,
    bc=bc,gc=gc,tot=tot,D=D,elig=elig,obs=obs,EDGES=EDGES),open('/tmp/s08.pkl','wb'))
print("\nSTEP1 DONE",flush=True)
