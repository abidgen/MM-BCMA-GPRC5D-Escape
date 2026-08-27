# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 07
# Step         : s07f_ig_vj_donor_and_repeated.py
# What it does : donor V/J summary, repeated-sample consistency, V evaluability
# Writes       : results/07_malignant_plasma/ig_clone_feasibility/normal_donor_vj_summary.csv; repeated_sample_vj_consistency.csv; vj_evaluability.csv
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : heredoc-written script, run with `conda run -n mm-core python <path>`
#   Original path  : /tmp/claude-1000/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646/scratchpad/s06/vj_audit2.py  (session scratchpad, since deleted)
#   Executed (UTC) : 2026-08-25T15:34:35Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
"""V/J feasibility part 2: donor negative control, repeated samples, evaluability."""
import numpy as np, pandas as pd, scanpy as sc, re
from pathlib import Path
OUT = Path('results/07_malignant_plasma/ig_clone_feasibility'); pd.set_option('display.width',250)

g = sc.read_h5ad('results/05_integration/integrated.h5ad')
lab = pd.read_csv('results/06_annotation/per_cell_labels.csv.gz', index_col=0)
g.obs['cell_type'] = lab.reindex(g.obs_names)['cell_type'].astype(str).values
pl = g[(g.obs['cell_type']=='PlasmaCell').values].copy(); del g
V = list(pl.var_names); C = pl.layers['counts']
KV=[x for x in V if x.startswith('IGKV')]; LV=[x for x in V if x.startswith('IGLV')]
def sub(gs):
    idx=[V.index(x) for x in gs]; X=C[:,idx]
    return np.asarray(X.todense()) if hasattr(X,'todense') else np.asarray(X)
XKV, XLV = sub(KV), sub(LV)
o = pd.read_csv(OUT/'ig_vj_detection_per_cell.csv.gz', index_col=0).reindex(pl.obs_names)

def dominant_V(mask, cls):
    X, gs = (XKV, KV) if cls=='kappa' else (XLV, LV)
    D = X[mask]>0
    npos = int(D.any(axis=1).sum())
    if npos == 0: return None, np.nan, 0
    cnt = D.sum(axis=0); i = int(np.argmax(cnt))
    return gs[i], float(cnt[i]/npos), npos

# ---- 7. DONOR NEGATIVE CONTROL: does the procedure invent clones in normal marrow? ----
print("=== 7. NORMAL-DONOR NEGATIVE CONTROL (per donor, not pooled) ===")
rows=[]
for s, gg in o[o.sample_type=='normal_bm'].groupby('sample_name', observed=True):
    m = o.index.get_indexer(gg.index)
    cls = gg.lc_class.value_counts().reindex(['kappa','lambda']).fillna(0).idxmax()
    gene, frac, npos = dominant_V(m, cls)
    rows.append({'donor':s,'n_plasma':len(gg),'dominant_class':cls,
                 'pct_with_LCV':round(100*gg.has_LCV.mean(),2),
                 'top_V_gene':gene,'top_V_frac_of_Vpos':round(frac,3) if frac==frac else np.nan,
                 'n_V_positive':npos})
D=pd.DataFrame(rows).sort_values('n_plasma',ascending=False)
D.to_csv(OUT/'normal_donor_vj_summary.csv',index=False)
print(D.to_string(index=False))
q=D[D.n_V_positive>=20]
print(f"\n  donors with >=20 V-positive cells: {len(q)}")
print(f"  their top_V_frac range: {q.top_V_frac_of_Vpos.min():.3f} - {q.top_V_frac_of_Vpos.max():.3f}"
      f"  median {q.top_V_frac_of_Vpos.median():.3f}")
P = pd.read_csv(OUT/'patient_vj_summary.csv')
dis = P[(P.sample_type=='myeloma') & (P.n_V_positive>=20)]
print(f"  myeloma patients top_V_frac: min {dis.top_V_frac_of_Vpos.min():.3f} "
      f"q25 {dis.top_V_frac_of_Vpos.quantile(.25):.3f} median {dis.top_V_frac_of_Vpos.median():.3f}")
print(f"  >>> SEPARATION: donor max {q.top_V_frac_of_Vpos.max():.3f} vs myeloma q25 {dis.top_V_frac_of_Vpos.quantile(.25):.3f}")

# ---- 6. REPEATED-SAMPLE CONSISTENCY ----
REPEAT=['27522','47491','56203','58408','59114','60359','81012','83942']
print("\n=== 6. REPEATED-SAMPLE CONSISTENCY (does the same V gene recur per sample?) ===")
rows=[]
for p in REPEAT:
    gg = o[o.patient_id.astype(str)==p]
    if not len(gg): continue
    cls = gg.lc_class.value_counts().reindex(['kappa','lambda']).fillna(0).idxmax()
    pg,_,_ = dominant_V(o.index.get_indexer(gg.index), cls)
    for s, hh in gg.groupby('sample_name', observed=True):
        m=o.index.get_indexer(hh.index); sg, sf, sn = dominant_V(m, cls)
        rows.append({'patient':p,'sample':s,'n_plasma':len(hh),'dominant_class':cls,
                     'patient_top_V':pg,'sample_top_V':sg,
                     'sample_top_V_frac':round(sf,3) if sf==sf else np.nan,
                     'n_V_pos':sn,'matches_patient':sg==pg})
R=pd.DataFrame(rows); R.to_csv(OUT/'repeated_sample_vj_consistency.csv',index=False)
print(R.to_string(index=False))
ev=R[R.n_V_pos>=20]
print(f"\n  samples with >=20 V-positive cells: {len(ev)} | matching the patient's top V gene: "
      f"{int(ev.matches_patient.sum())} ({100*ev.matches_patient.mean():.1f}%)")

# ---- 5/8. CELL-LEVEL COORDINATION + EVALUABILITY (distributions only) ----
print("\n=== 8. EVALUABILITY (distributions; gate predeclared below, not tuned) ===")
rows=[]
for p, gg in o.groupby('patient_id', observed=True):
    cls = gg.lc_class.value_counts().reindex(['kappa','lambda']).fillna(0).idxmax()
    m=o.index.get_indexer(gg.index); gene, frac, npos = dominant_V(m, cls)
    rows.append({'patient':p,'sample_type':gg.sample_type.iloc[0],'n_plasma':len(gg),
      'pct_LCV':round(100*gg.has_LCV.mean(),2),'n_V_positive':npos,
      'top_V_frac':round(frac,3) if frac==frac else np.nan})
E=pd.DataFrame(rows)
E['state']=np.where((E.n_V_positive>=50)&(E.pct_LCV>=50),'V_EVALUABLE',
             np.where((E.n_V_positive>=20)&(E.pct_LCV>=20),'PARTIALLY_EVALUABLE','NOT_EVALUABLE'))
E.to_csv(OUT/'vj_evaluability.csv',index=False)
print(E.groupby(['sample_type','state']).size().to_string())
print(f"\n  myeloma patients V_EVALUABLE: {int(((E.sample_type=='myeloma')&(E.state=='V_EVALUABLE')).sum())} "
      f"of {int((E.sample_type=='myeloma').sum())}")
print("\nDONE", flush=True)
