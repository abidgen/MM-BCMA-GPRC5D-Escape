# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 08
# Step         : s08b_null_scheme_truncate10k_denominator_ambient.py
# What it does : cohort vs global null scheme, truncate-all-at-10k, primary-vs-sensitivity, ambient floor
# Writes       : results/08_dual_antigen_escape/depth_stratified_null.csv; truncate10k_sensitivity.csv; primary_vs_sensitivity_denominator_comparison.csv; noise_floor_ambient.csv
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : heredoc-written script, run with `conda run -n mm-core python <path>`
#   Original path  : $S/s08b.py  (session scratchpad, since deleted)
#   Executed (UTC) : 2026-08-25T16:08:02Z
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
OUT=Path('results/08_dual_antigen_escape'); pd.set_option('display.width',260)
S=pickle.load(open('/tmp/s08.pkl','rb')); SEED=20260825
prim,sens,strata_c,strata_g=S['prim'],S['sens'],S['strata_c'],S['strata_g']
bc,gc,tot,D,elig,obs=S['bc'],S['gc'],S['tot'],S['D'],S['elig'],S['obs']

def summarize(mask,strata,bcma,gprc):
    r=[]
    for p in sorted(elig):
        m=mask&(obs.patient_id==p).values
        if m.sum()==0: continue
        b0=bcma[m]==0; g0=gprc[m]==0; st=A.merge_sparse_strata(strata[m],20)
        o=float((b0&g0).mean()); e=A.stratified_expected_dn(b0,g0,st)
        u=A.unconditioned_expected_dn(b0,g0)
        r.append({'patient':p,'cohort':obs.cohort.values[m][0],'n':int(m.sum()),
            'bcma_detect':float((bcma[m]>0).mean()),'gprc5d_detect':float((gprc[m]>0).mean()),
            'obs_dn':o,'exp_strat':e,'exp_uncond':u,
            'enr_strat':o/e if e>0 else np.nan,'enr_uncond':o/u if u>0 else np.nan})
    return pd.DataFrame(r)

# ---- 4/6. cohort-specific vs global bins ----
CS=summarize(prim,strata_c,bc,gc); GL=summarize(prim,strata_g,bc,gc)
cmp=CS[['patient','cohort','n','obs_dn','exp_strat','enr_strat']].merge(
    GL[['patient','exp_strat','enr_strat']],on='patient',suffixes=('_cohortbins','_globalbins'))
cmp['enr_delta']=cmp.enr_strat_globalbins-cmp.enr_strat_cohortbins
cmp=cmp.merge(CS[['patient','enr_uncond']],on='patient')
cmp.round(4).to_csv(OUT/'depth_stratified_null.csv',index=False)
print("=== NULL SCHEME COMPARISON (primary denominator) ===")
print(cmp.round(3).to_string(index=False))
print(f"\n  |enrichment delta| cohort-bins vs global-bins: median {cmp.enr_delta.abs().median():.4f} "
      f"max {cmp.enr_delta.abs().max():.4f}")
print(f"  stratified enrichment: median {cmp.enr_strat_cohortbins.median():.3f}  "
      f"range {cmp.enr_strat_cohortbins.min():.3f}-{cmp.enr_strat_cohortbins.max():.3f}")
print(f"  UNCONDITIONED enrichment: median {cmp.enr_uncond.median():.3f}  "
      f"range {cmp.enr_uncond.min():.3f}-{cmp.enr_uncond.max():.3f}")
print(f"  --> depth artifact = median {(cmp.enr_uncond-cmp.enr_strat_cohortbins).median():.4f}")

# ---- 10. truncate all cohorts at 10k ----
bt=A.downsample_gene_counts(tot,bc,10000,SEED); gt=A.downsample_gene_counts(tot,gc,10000,SEED)
tt=np.minimum(tot,10000); Dt=A.depth_ex_antigen(tt,bt,gt)
st_t=np.zeros(len(obs),dtype=np.int64)
for ch in sorted(obs.cohort[prim].unique()):
    m=(obs.cohort==ch).values
    e=A.quantile_edges(Dt[prim&m],5); st_t[m]=A.assign_strata(Dt[m],e)
TR=summarize(prim,st_t,bt,gt)
tc=CS[['patient','cohort','n','bcma_detect','gprc5d_detect','obs_dn','enr_strat']].merge(
    TR[['patient','bcma_detect','gprc5d_detect','obs_dn','enr_strat']],on='patient',
    suffixes=('_orig','_trunc10k'))
tc['dn_delta']=tc.obs_dn_trunc10k-tc.obs_dn_orig
tc['gprc5d_delta']=tc.gprc5d_detect_trunc10k-tc.gprc5d_detect_orig
tc.round(4).to_csv(OUT/'truncate10k_sensitivity.csv',index=False)
print("\n=== 11. TRUNCATE-ALL-AT-10k (primary) ==="); print(tc.round(3).to_string(index=False))
print("\n  by cohort, mean delta:"); print(tc.groupby('cohort')[['dn_delta','gprc5d_delta']].agg(['mean','max']).round(4).to_string())
print(f"  n cells above 10k UMI: {int((tot>10000).sum())} of {len(tot)} "
      f"({int((tot[prim]>10000).sum())} of {int(prim.sum())} in primary denominator)")
print("  spearman rank stability of obs_dn orig vs trunc:",
      round(tc.obs_dn_orig.corr(tc.obs_dn_trunc10k,method='spearman'),4))

# ---- 11/12. primary vs sensitivity ----
SS=summarize(sens,strata_c,bc,gc)
ps=CS[['patient','cohort','n','bcma_detect','gprc5d_detect','obs_dn','enr_strat']].merge(
    SS[['patient','n','bcma_detect','gprc5d_detect','obs_dn','enr_strat']],on='patient',
    suffixes=('_primary','_sensitivity'))
ps['dn_delta']=ps.obs_dn_sensitivity-ps.obs_dn_primary
ps['expansion']=ps.n_sensitivity/ps.n_primary
ps.round(4).to_csv(OUT/'primary_vs_sensitivity_denominator_comparison.csv',index=False)
print("\n=== 12. PRIMARY vs SENSITIVITY ==="); print(ps.round(3).to_string(index=False))
print(f"\n  DN delta: median {ps.dn_delta.median():+.4f}  range {ps.dn_delta.min():+.4f} to {ps.dn_delta.max():+.4f}")
print(f"  patients with |delta| > 0.05: {int((ps.dn_delta.abs()>0.05).sum())} of {len(ps)}")
print(f"  spearman(primary, sensitivity) = {ps.obs_dn_primary.corr(ps.obs_dn_sensitivity,method='spearman'):.4f}")
print("  by cohort mean DN delta:"); print(ps.groupby('cohort').dn_delta.agg(['mean','max']).round(4).to_string())

# ---- 4/8. noise floor from NON-PLASMA ----
NP=['Tcell','Myeloid','Bcell','HSPC','NK']
npm=obs.cell_type.isin(NP).values & obs.cohort.isin(['MMRF','WU1','WU2']).values
print(f"\n=== 8. TECHNICAL-ZERO / AMBIENT REFERENCE ===")
print(f"  populations used: {NP}; 'Ambiguous' (Leiden 23) excluded by name")
print(f"  reference cells: {int(npm.sum())}  (excluded Ambiguous: {int((obs.cell_type=='Ambiguous').sum())})")
rows=[]
for ch in ['MMRF','WU1','WU2']:
    e=S['EDGES'][ch]
    for i in range(len(e)-1):
        m=npm&(obs.cohort==ch).values&(A.assign_strata(D,e)==i)
        if m.sum()<50: continue
        rows.append({'cohort':ch,'stratum':i,'n_reference_cells':int(m.sum()),
            'median_depth':float(np.median(D[m])),
            'ambient_bcma_detect':float((bc[m]>0).mean()),
            'ambient_gprc5d_detect':float((gc[m]>0).mean())})
NF=pd.DataFrame(rows)
print(NF.round(4).to_string(index=False))
NF.to_csv(OUT/'noise_floor_ambient.csv',index=False)
pickle.dump(dict(CS=CS,cmp=cmp,tc=tc,ps=ps,NF=NF,npm=npm),open('/tmp/s08b.pkl','wb'))
print("\nSTEP2 DONE",flush=True)
