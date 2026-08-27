# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 10
# Step         : s10d_hypotheses_and_evidence_levels.py
# What it does : pre-registered gamma-secretase test, MYC/OXPHOS, stress/IFN/UPR/AP, Level-2 state calls, evidence levels
# Writes       : results/10_dn_coherence/gamma_secretase_hypothesis.csv; myc_oxphos_programs.csv; stress_ifn_upr_antigen_presentation.csv; stage10_evidence_levels.csv; pseudobulk_de_evaluability.csv
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : heredoc-written script, run with `conda run -n mm-core python <path>`
#   Original path  : $S/s10d.py  (session scratchpad, since deleted)
#   Executed (UTC) : 2026-08-25T22:57:42Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
import numpy as np, pandas as pd, pickle, warnings
from pathlib import Path
from scipy.stats import wilcoxon
from mm_escape import subclone as SC, config, communication as CM
warnings.filterwarnings('ignore')
OUT=Path('results/10_dn_coherence'); pd.set_option('display.width',300); SEED=20260825
PP,R,repro,DS,PROGS,scores,obs,is_dn,depth,C,genes=pickle.load(open('/tmp/s10c.pkl','rb'))

# ---- gamma-secretase, pre-registered: reported in full regardless of outcome ----
gs=R[(R.program=='gamma_secretase')]
gs.round(5).to_csv(OUT/'gamma_secretase_hypothesis.csv',index=False)
print("=== PRE-REGISTERED GAMMA-SECRETASE HYPOTHESIS ===")
print(f"  genes (frozen, unchanged): {config.STATE_PROGRAMS['gamma_secretase']}")
print(gs.round(4).to_string(index=False))
per=PP.pivot(index='patient',columns='denominator',values='gamma_secretase_delta_matched')
print(f"  patients with positive matched effect: primary {int((per['primary']>0).sum())}/{per['primary'].notna().sum()},"
      f" sensitivity {int((per['sensitivity']>0).sum())}/{per['sensitivity'].notna().sum()}")
print(f"  VERDICT: gamma_secretase in cohort-reproducible set? {'gamma_secretase' in repro}")

PP[[c for c in PP.columns if c.startswith(('patient','denominator','myc','oxphos'))]+['n_matched_per_group']]\
  .round(5).to_csv(OUT/'myc_oxphos_programs.csv',index=False)
PP[[c for c in PP.columns if c.startswith(('patient','denominator','stress','interferon','upr','antigen_presentation'))]]\
  .round(5).to_csv(OUT/'stress_ifn_upr_antigen_presentation.csv',index=False)

# ---- per-patient LEVEL 2 ----
EV=pd.read_csv(OUT/'stage10_evaluability.csv'); F1=pd.read_csv(OUT/'dn_coherence_final_states.csv')
rs=pd.read_csv(OUT/'repeated_sample_dn_coherence.csv')[['patient','repeated_sample_status']].drop_duplicates()
rs['patient']=rs.patient.astype(str)
piv={k:PP.pivot(index='patient',columns='denominator',values=f'{k}_delta_matched') for k in config.LEVEL2_PROGRAMS}
rows=[]
for _,f in F1.iterrows():
    p=str(f.patient)
    ev=bool(f.evaluable_primary) and bool(f.evaluable_sensitivity)
    l1=('DN_STRUCTURE_SUPPORTED' if f.dn_coherence_state=='DN_COHERENCE_SUPPORTED' else
        'DN_STRUCTURE_NOT_SUPPORTED' if f.dn_coherence_state=='DN_COHERENCE_NOT_SUPPORTED' else
        'DN_STRUCTURE_NOT_EVALUABLE')
    if not ev or not repro:
        l2='DN_STATE_NOT_EVALUABLE' if not ev else 'DN_STATE_NOT_SUPPORTED'; hits=[]
    else:
        hits=[k for k in sorted(repro)
              if p in piv[k].index and np.isfinite(piv[k].loc[p,'primary']) and np.isfinite(piv[k].loc[p,'sensitivity'])
              and np.sign(piv[k].loc[p,'primary'])==np.sign(piv[k].loc[p,'sensitivity'])
              and np.sign(piv[k].loc[p,'primary'])==np.sign(
                  R[(R.program==k)&(R.effect=='matched')&(R.denominator=='primary')].median_delta.iloc[0])]
        l2='DN_STATE_SUPPORTED' if hits else 'DN_STATE_NOT_SUPPORTED'
    rows.append({'patient':p,'cohort':f.cohort,'level1_structure':l1,'level2_state':l2,
        'level3_cnv':SC.CNV_NOT_EVALUABLE,'reproducible_programs_matched':';'.join(hits) or 'none',
        'n_dn_primary':f.n_dn_primary,'repeated_sample_status':f.repeated_sample_status,
        'licensed_language':('escape-associated transcriptional state' if l2=='DN_STATE_SUPPORTED'
            else 'non-random DN organization' if l1=='DN_STRUCTURE_SUPPORTED' else 'none')})
L=pd.DataFrame(rows); L.to_csv(OUT/'stage10_evidence_levels.csv',index=False)
print("\n=== STAGE-10 EVIDENCE LEVELS ===")
print(pd.crosstab(L.level1_structure,L.level2_state).to_string())
print(L[['patient','cohort','level1_structure','level2_state','reproducible_programs_matched','licensed_language']].to_string(index=False))

# ---- pseudobulk DE evaluability audit ----
elig=PP[PP.n_matched_per_group>=SC.MIN_GROUP_CELLS]
aud=[]
for tag in ['primary','sensitivity']:
    s=elig[elig.denominator==tag]
    aud.append({'denominator':tag,'n_patients_with_both_groups_ge20_matched':len(s),
        'min_required':10,'estimable':bool(len(s)>=10),
        'design':'~ patient + group (paired, patient as blocking factor)',
        'aggregation':'sum raw counts over depth-matched cells',
        'note':'per-cell DE is NOT used as a substitute if this is not estimable'})
AU=pd.DataFrame(aud); AU.to_csv(OUT/'pseudobulk_de_evaluability.csv',index=False)
print("\n=== PSEUDOBULK DE EVALUABILITY ===" ); print(AU.to_string(index=False))
pickle.dump((L,elig,repro),open('/tmp/s10d.pkl','wb'))
