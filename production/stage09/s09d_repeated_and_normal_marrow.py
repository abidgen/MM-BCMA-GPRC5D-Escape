# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 09
# Step         : s09d_repeated_and_normal_marrow.py
# What it does : repeated-sample bulk consistency and normal-marrow expression context
# Writes       : results/09_bulk_validation/repeated_sample_bulk_consistency.csv; normal_marrow_antigen_context.csv
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : heredoc-written script, run with `conda run -n mm-core python <path>`
#   Original path  : $S/s09d.py  (session scratchpad, since deleted)
#   Executed (UTC) : 2026-08-25T16:38:54Z
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
St=pickle.load(open('/tmp/s08.pkl','rb')); D,R,mt=pickle.load(open('/tmp/s09_res.pkl','rb'))
B=pickle.load(open('/tmp/s09_bulk.pkl','rb')); M=pickle.load(open('/tmp/s09_map.pkl','rb'))
obs=St['obs']; prim=St['prim']; bc,gc,tot=St['bc'],St['gc'],St['tot']
pd.set_option('display.width',250)

# ---- 8. WashU GPRC5D interpretation (context only; changes nothing) ----
w=D[D.cohort=='WU1'].copy()
w['bulk_GPRC5D_rank']=w.bulk_GPRC5D_tpm.rank(ascending=False)
print("=== 11. WASHU GPRC5D: bulk abundance vs single-cell detection ===")
print(w[['sc_patient','sc_sample','n_cells_primary','bulk_GPRC5D_tpm','sc_GPRC5D_detect_primary',
         'sc_GPRC5D_cpm_primary','bulk_TNFRSF17_tpm','sc_TNFRSF17_detect_primary']]
      .round(4).to_string(index=False))
print("\n  MMRF for contrast:")
print(D[D.cohort=='MMRF'][['sc_patient','bulk_GPRC5D_tpm','sc_GPRC5D_detect_primary',
      'sc_GPRC5D_cpm_primary']].round(4).to_string(index=False))

# ---- 13. repeated samples / timepoints ----
rs=[]
for _,r in mt.iterrows():
    m=prim&(obs.sample_name==r.sc_sample).values
    rs.append({'sc_patient':r.sc_patient,'sc_sample':r.sc_sample,'bulk_sample':r.bulk_sample,
        'timepoint':r.timepoint,'n_primary_cells':int(m.sum()),
        'assessable_ge20':bool(m.sum()>=20),
        'used_in_correlation':bool(r.sc_sample in set(D.sc_sample)),
        'bulk_TNFRSF17_tpm':r.TNFRSF17_tpm,'bulk_GPRC5D_tpm':r.GPRC5D_tpm,
        'sc_TNFRSF17_cpm':1e6*bc[m].sum()/tot[m].sum() if m.sum() and tot[m].sum()>0 else np.nan,
        'sc_GPRC5D_cpm':1e6*gc[m].sum()/tot[m].sum() if m.sum() and tot[m].sum()>0 else np.nan})
RS=pd.DataFrame(rs).sort_values(['sc_patient','timepoint'])
RS.round(3).to_csv(OUT/'repeated_sample_bulk_consistency.csv',index=False)
print("\n=== 13. REPEATED SAMPLES / TIMEPOINTS (never pooled, never counted as patients) ===")
print(RS.round(2).to_string(index=False))

# ---- 11. normal marrow, DONOR as the unit ----
don=(obs.sample_type=='donor').values if 'donor' in set(obs.sample_type) else \
    (~obs.sample_type.isin(['myeloma'])).values
pl=(obs.cell_type=='PlasmaCell').values
rows=[]
for s in sorted(obs.sample_name[don].unique()):
    for pop,mask in [('PlasmaCell',pl),('Tcell',(obs.cell_type=='Tcell').values),
                     ('Myeloid',(obs.cell_type=='Myeloid').values)]:
        m=don&mask&(obs.sample_name==s).values
        if m.sum()<10: continue
        rows.append({'donor':s,'population':pop,'n_cells':int(m.sum()),
            'median_total_umi':float(np.median(tot[m])),
            'TNFRSF17_detect':round(float((bc[m]>0).mean()),4),
            'GPRC5D_detect':round(float((gc[m]>0).mean()),4),
            'TNFRSF17_cpm':round(1e6*bc[m].sum()/tot[m].sum(),2),
            'GPRC5D_cpm':round(1e6*gc[m].sum()/tot[m].sum(),2)})
N=pd.DataFrame(rows); N.to_csv(OUT/'normal_marrow_antigen_context.csv',index=False)
print("\n=== 12. NORMAL MARROW CONTEXT (donor = biological unit) ===")
print(N.to_string(index=False))
npc=N[N.population=='PlasmaCell']
print(f"\n  normal plasma cells: {len(npc)} donors, {int(npc.n_cells.sum())} cells total")
print(f"  TNFRSF17 detection across donors: median {npc.TNFRSF17_detect.median():.3f} "
      f"range {npc.TNFRSF17_detect.min():.3f}-{npc.TNFRSF17_detect.max():.3f}")
print(f"  GPRC5D   detection across donors: median {npc.GPRC5D_detect.median():.3f} "
      f"range {npc.GPRC5D_detect.min():.3f}-{npc.GPRC5D_detect.max():.3f}")
print(f"  depth dependence (Spearman detection ~ median UMI across donors): ", end='')
from scipy.stats import spearmanr
print(f"BCMA {spearmanr(npc.median_total_umi,npc.TNFRSF17_detect)[0]:.3f}, "
      f"GPRC5D {spearmanr(npc.median_total_umi,npc.GPRC5D_detect)[0]:.3f}")
print("\n  malignant (primary denominator) for contrast: "
      f"BCMA detect median {D.sc_TNFRSF17_detect_primary.median():.3f}, "
      f"GPRC5D {D.sc_GPRC5D_detect_primary.median():.3f}")
