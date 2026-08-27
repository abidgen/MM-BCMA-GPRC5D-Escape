# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 09b
# Step         : s09b2_provisional_tiers.py
# What it does : tier assignment + tau threshold sensitivity
# Writes       : results/08_dual_antigen_escape/risk_tier_provisional/risk_tiers_provisional.csv (written as risk_tiers_final.csv, later renamed); risk_tiers_tau0{20,25,33}.csv; risk_tier_threshold_comparison.csv; risk_tier_evidence_long.csv
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : heredoc-written script, run with `conda run -n mm-core python <path>`
#   Original path  : $S/tier.py  (session scratchpad, since deleted)
#   Executed (UTC) : 2026-08-25T21:03:45Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
import numpy as np, pandas as pd
from pathlib import Path
from mm_escape import risk_tiers as RT
RT_DIR=Path('results/08_dual_antigen_escape/risk_tier_final'); RT_DIR.mkdir(exist_ok=True)
EV=pd.read_csv('results/08_dual_antigen_escape/risk_tier_design/patient_evidence_matrix.csv')
MAP=pd.read_csv('results/09_bulk_validation/sc_bulk_sample_mapping.csv')
SC9=pd.read_csv('results/09_bulk_validation/sc_marginal_antigen_primary.csv')
pd.set_option('display.width',300)

def parse_ps(s):
    if not isinstance(s,str) or not s: return []
    out=[]
    for tok in s.split('|'):
        a,b,c=tok.split(':'); out.append((a,int(b),float(c)))
    return out

rows=[]
for _,r in EV.iterrows():
    ev=dict(n_primary=int(r.n_primary), dn_primary=r.obs_dn_primary,
            dn_sensitivity=r.obs_dn_sensitivity, dn_trunc10k=r.obs_dn_trunc10k_primary,
            enr_ci_lo=r.enr_lo_primary, enr_cohortbins=r.enr_primary,
            enr_globalbins=r.enr_globalbins_primary,
            per_sample=parse_ps(r.per_sample_dn_primary))
    final,robust,tiers,reasons=RT.final_tier(ev)
    agree={t:RT.sample_agreement(ev['per_sample'],t) for t in RT.TAU_HIGH_SET}
    d={'patient':r.patient,'cohort':r.cohort,'final_tier':final,
       'threshold_robustness':robust,
       'tier_tau020':tiers[0.20],'tier_tau025':tiers[0.25],'tier_tau033':tiers[0.33],
       'n_primary':int(r.n_primary),'n_samples_assessable':int(r.n_samples_assessable_primary),
       'repeated_sample_status':agree[0.25],
       'dn_primary':r.obs_dn_primary,'dn_sensitivity':r.obs_dn_sensitivity,
       'dn_trunc10k':r.obs_dn_trunc10k_primary,'dn_ci_lo':r.dn_lo_primary,'dn_ci_hi':r.dn_hi_primary,
       'enr_cohortbins':r.enr_primary,'enr_ci_lo':r.enr_lo_primary,'enr_ci_hi':r.enr_hi_primary,
       'excess_dn':r.excess_primary,'enr_globalbins':r.enr_globalbins_primary,
       'perm_p_reported_not_used':r.perm_p_primary,
       **reasons}
    d['uncertainty_reasons']=';'.join([k.replace('uncertain_','') for k,v in reasons.items() if v]) or 'none'
    rows.append(d)
T=pd.DataFrame(rows)

# ---- Stage-09 context: REPORT-ONLY, appended AFTER tiers exist ----
bulkpat=set(SC9.sc_patient.astype(str))
b=SC9.set_index(SC9.sc_patient.astype(str))
T['s09_bulk_available']=T.patient.astype(str).isin(bulkpat)
T['s09_bulk_TNFRSF17_tpm']=T.patient.astype(str).map(b.bulk_TNFRSF17_tpm)
T['s09_bulk_GPRC5D_tpm']=T.patient.astype(str).map(b.bulk_GPRC5D_tpm)
T['s09_marginal_context']=np.where(~T.s09_bulk_available,'no matched bulk',
    np.where(T.cohort=='MMRF','sorted-bulk cohort (marginal rho=0.93 both antigens)',
             'unsorted-bulk cohort (n=5, underpowered; WashU GPRC5D dropout context)'))

for tau,col in [(0.20,'tier_tau020'),(0.25,'tier_tau025'),(0.33,'tier_tau033')]:
    T[['patient','cohort',col,'n_primary','dn_primary','dn_sensitivity','dn_trunc10k',
       'enr_ci_lo','uncertainty_reasons']].rename(columns={col:'tier'}).round(4)\
      .to_csv(RT_DIR/f'risk_tiers_tau{str(tau).replace("0.","0")}.csv',index=False)
T[['patient','cohort','tier_tau020','tier_tau025','tier_tau033','threshold_robustness','final_tier']]\
  .to_csv(RT_DIR/'risk_tier_threshold_comparison.csv',index=False)
T.round(5).to_csv(RT_DIR/'risk_tiers_final.csv',index=False)
long=T.melt(id_vars=['patient','cohort','final_tier'],
    value_vars=list(RT.UNCERTAINTY_REASONS),var_name='reason',value_name='fired')
long[long.fired].to_csv(RT_DIR/'risk_tier_evidence_long.csv',index=False)

print("=== TIER AT EACH TAU_HIGH ===")
for c in ['tier_tau020','tier_tau025','tier_tau033']:
    print(f"  {c}: {dict(T[c].value_counts())}")
print(f"\n  threshold robustness: {dict(T.threshold_robustness.value_counts())}")
print(f"  FINAL: {dict(T.final_tier.value_counts())}")
print("\n=== FINAL TIERS ===")
print(T[['patient','cohort','n_primary','dn_primary','dn_sensitivity','dn_trunc10k','enr_ci_lo',
         'tier_tau020','tier_tau025','tier_tau033','threshold_robustness','final_tier',
         'uncertainty_reasons']].round(3).to_string(index=False))
print("\n=== UNCERTAINTY REASON FREQUENCIES ===")
print(T[list(RT.UNCERTAINTY_REASONS)].sum().sort_values(ascending=False).to_string())
T.to_pickle('/tmp/tiers.pkl')
