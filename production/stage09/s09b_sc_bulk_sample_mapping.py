# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 09
# Step         : s09b_sc_bulk_sample_mapping.py
# What it does : matched bulk/scRNA pairing, cohort-split assay rule
# Writes       : results/09_bulk_validation/sc_bulk_sample_mapping.csv
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : heredoc-written script, run with `conda run -n mm-core python <path>`
#   Original path  : $S/s09b.py  (session scratchpad, since deleted)
#   Executed (UTC) : 2026-08-25T16:37:23Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
import numpy as np, pandas as pd, pickle
from pathlib import Path
OUT=Path('results/09_bulk_validation')
St=pickle.load(open('/tmp/s08.pkl','rb')); B=pickle.load(open('/tmp/s09_bulk.pkl','rb'))
obs=St['obs']; elig=St['elig']; prim,sens=St['prim'],St['sens']
pd.set_option('display.width',250)

sc_samples=sorted(obs.sample_name.unique())
sc_by_sample={s:obs.patient_id[(obs.sample_name==s).values].iloc[0] for s in sc_samples}
print(f"scRNA samples {len(sc_samples)}; eligible (frozen-denominator) patients {len(elig)}")

rows=[]
for _,r in B.iterrows():
    bs=r.bulk_sample
    if not r.usable:
        rows.append({'sc_patient':None,'sc_sample':None,'bulk_sample':bs,'timepoint':None,
            'match_status':'NOT_EVALUABLE','notes':'empty bulk file (114-byte stub)'}); continue
    if bs in sc_by_sample:
        p=sc_by_sample[bs]
        tp=bs.split('_')[-1] if ('_' in bs and bs.split('_')[-1].isdigit()) else 'single'
        ok=p in elig
        rows.append({'sc_patient':p,'sc_sample':bs,'bulk_sample':bs,'timepoint':tp,
            'match_status':'MATCHED_EXACT' if ok else 'NOT_EVALUABLE',
            'notes':'' if ok else 'sc patient has no frozen Stage-07 denominator '
                                  '(not CLONAL_STRONG & V_EVALUABLE)'})
    else:
        near=[s for s in sc_samples if s.replace('MMY','')==bs or bs.replace('MMY','')==s]
        cand=[s for s in sc_samples if abs(len(s)-len(bs))<=3 and s[:3]==bs[:3]]
        rows.append({'sc_patient':None,'sc_sample':None,'bulk_sample':bs,'timepoint':None,
            'match_status':'NOT_EVALUABLE',
            'notes':f'no scRNA sample with this identifier; nearest by prefix: '
                    f'{cand[:4] if cand else "none"} — NOT matched (identifiers differ)'})
M=pd.DataFrame(rows)
M=M.merge(B[['bulk_sample','bulk_cohort','specimen','TNFRSF17_tpm','GPRC5D_tpm','n_runs']],
          on='bulk_sample',how='left')
M.to_csv(OUT/'sc_bulk_sample_mapping.csv',index=False)
print("\n=== MAPPING TABLE ===")
print(M[['sc_patient','sc_sample','bulk_sample','bulk_cohort','timepoint','match_status','notes']]
      .to_string(index=False))
print("\n", M.match_status.value_counts().to_string())
mt=M[M.match_status=='MATCHED_EXACT']
print(f"\nusable matched sc-bulk pairs: {len(mt)} across {mt.sc_patient.nunique()} patients")
print("  by cohort:"); print(mt.groupby('bulk_cohort').agg(pairs=('bulk_sample','size'),
      patients=('sc_patient','nunique')).to_string())
print("  one-to-many (patients with >1 matched bulk sample):",
      sorted(mt.sc_patient.value_counts()[lambda x:x>1].index.tolist()))
pickle.dump(M,open('/tmp/s09_map.pkl','wb'))
