# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 07
# Step         : s07e_ig_vj_audit.py
# What it does : IG V/J capture audit; J segments shown uncapturable at 3-prime
# Writes       : results/07_malignant_plasma/ig_clone_feasibility/ig_gene_availability.csv; ig_vj_detection_per_cell.csv.gz; patient_vj_summary.csv; cross_patient_specificity.csv
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : heredoc-written script, run with `conda run -n mm-core python <path>`
#   Original path  : /tmp/claude-1000/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646/scratchpad/s06/vj_audit.py  (session scratchpad, since deleted)
#   Executed (UTC) : 2026-08-25T15:33:06Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
"""IG V/J clone-membership FEASIBILITY audit. No cell is assigned malignant."""
import numpy as np, pandas as pd, scanpy as sc, re
from pathlib import Path
OUT = Path('results/07_malignant_plasma/ig_clone_feasibility'); pd.set_option('display.width',240)

g = sc.read_h5ad('results/05_integration/integrated.h5ad')
lab = pd.read_csv('results/06_annotation/per_cell_labels.csv.gz', index_col=0)
g.obs['cell_type'] = lab.reindex(g.obs_names)['cell_type'].astype(str).values
pl = g[(g.obs['cell_type']=='PlasmaCell').values].copy(); del g
assert 'counts' in pl.layers
V = list(pl.var_names)
assert not ({'TNFRSF17','GPRC5D'} & set(V[:0])), "n/a"

GRP = {
 'IGKV': [x for x in V if x.startswith('IGKV')], 'IGKJ': [x for x in V if x.startswith('IGKJ')],
 'IGKC': [x for x in V if x=='IGKC'],
 'IGLV': [x for x in V if x.startswith('IGLV')], 'IGLJ': [x for x in V if x.startswith('IGLJ')],
 'IGLC': [x for x in V if re.match(r'^IGLC\d',x)],
 'IGHV': [x for x in V if x.startswith('IGHV')], 'IGHJ': [x for x in V if x.startswith('IGHJ')],
 'IGHD_seg': [x for x in V if re.match(r'^IGHD\d',x)],
 'IGHC': [x for x in V if x in ('IGHG1','IGHG2','IGHG3','IGHG4','IGHA1','IGHA2','IGHM','IGHD','IGHE')],
}
C = pl.layers['counts']
def sub(gs):
    if not gs: return np.zeros((pl.n_obs,0))
    idx=[V.index(x) for x in gs]; X=C[:,idx]
    return np.asarray(X.todense()) if hasattr(X,'todense') else np.asarray(X)

rows=[]
for k,gs in GRP.items():
    X = sub(gs)
    det = (X>0)
    rows.append({'group':k,'n_genes_in_var':len(gs),
      'n_genes_ever_detected':int((det.sum(axis=0)>0).sum()),
      'pct_cells_any_detected':round(100*det.any(axis=1).mean(),2),
      'median_genes_detected_per_cell':float(np.median(det.sum(axis=1))),
      'median_umi_per_cell':float(np.median(X.sum(axis=1))),
      'pct_cells_ge2_genes':round(100*(det.sum(axis=1)>=2).mean(),2)})
av=pd.DataFrame(rows); av.to_csv(OUT/'ig_gene_availability.csv',index=False)
print("=== 1. IG GENE AVAILABILITY vs ACTUAL DETECTION (35,474 plasma cells) ===")
print(av.to_string(index=False), flush=True)

# per-cell V/J evidence
o = pl.obs[['sample_name','patient_id','cohort','sample_type','total_counts','n_genes_by_counts']].copy()
for k in ['IGKV','IGKJ','IGLV','IGLJ','IGHV','IGHJ']:
    X=sub(GRP[k]); o[f'n_{k}']= (X>0).sum(axis=1); o[f'umi_{k}']=X.sum(axis=1)
kc, lc = sub(GRP['IGKC']).sum(axis=1), sub(GRP['IGLC']).sum(axis=1)
o['kappa_umi'], o['lambda_umi'] = kc, lc
tot = kc+lc
o['lc_class'] = np.where(tot<3,'insufficient',
                  np.where(kc/np.maximum(tot,1)>=0.8,'kappa',
                    np.where(lc/np.maximum(tot,1)>=0.8,'lambda','ambiguous')))
o['has_LCV'] = (o.n_IGKV+o.n_IGLV)>0
o['has_LCJ'] = (o.n_IGKJ+o.n_IGLJ)>0
o['has_LCVJ']= o.has_LCV & o.has_LCJ
o['has_IGHV']= o.n_IGHV>0
o.to_csv(OUT/'ig_vj_detection_per_cell.csv.gz', compression='gzip')
print(f"\n=== 2. PER-CELL V/J DETECTION ===")
for c,lbl in [('has_LCV','any light-chain V'),('has_LCJ','any light-chain J'),
              ('has_LCVJ','light V AND J'),('has_IGHV','any IGHV')]:
    print(f"  {lbl:24s} {100*o[c].mean():5.2f}%  ({int(o[c].sum()):,} cells)")

# patient-level structure
def top_frac(X, gs, mask):
    if not gs or mask.sum()==0: return (np.nan, np.nan, 0)
    D=(X[mask]>0); n=D.any(axis=1).sum()
    if n==0: return (np.nan, np.nan, 0)
    counts=D.sum(axis=0); i=int(np.argmax(counts))
    return (gs[i], counts[i]/max(n,1), int(n))
KV, LV = sub(GRP['IGKV']), sub(GRP['IGLV'])
pats=[]
for p, idx in o.groupby('patient_id', observed=True).groups.items():
    m = o.index.get_indexer(idx)
    sub_o = o.loc[idx]
    dom = sub_o.lc_class.value_counts().reindex(['kappa','lambda']).fillna(0).idxmax()
    Xv, gsv = (KV, GRP['IGKV']) if dom=='kappa' else (LV, GRP['IGLV'])
    gene, frac, nv = top_frac(Xv, gsv, m)
    pats.append({'patient':p,'sample_type':sub_o.sample_type.iloc[0],'n_plasma':len(idx),
      'dominant_class':dom,'n_with_LCV':int(sub_o.has_LCV.sum()),
      'pct_with_LCV':round(100*sub_o.has_LCV.mean(),2),
      'pct_with_LCVJ':round(100*sub_o.has_LCVJ.mean(),2),
      'pct_with_IGHV':round(100*sub_o.has_IGHV.mean(),2),
      'top_V_gene':gene,'top_V_frac_of_Vpos':round(frac,3) if frac==frac else np.nan,
      'n_V_positive':nv})
P=pd.DataFrame(pats).sort_values('n_plasma',ascending=False)
P.to_csv(OUT/'patient_vj_summary.csv',index=False)
print(f"\n=== 3. PATIENT-LEVEL V/J STRUCTURE (top 12 by plasma count) ===")
print(P.head(12).to_string(index=False), flush=True)
print(f"\n  patients with >=50% of plasma cells LCV-positive: {int((P.pct_with_LCV>=50).sum())} / {len(P)}")
print(f"  patients with >=20% LCV-positive: {int((P.pct_with_LCV>=20).sum())} / {len(P)}")
print(f"  median pct_with_LCV across patients: {P.pct_with_LCV.median():.2f}%")
print(f"  median pct_with_LCVJ: {P.pct_with_LCVJ.median():.2f}%  | median pct IGHV: {P.pct_with_IGHV.median():.2f}%")

# cross-patient specificity of top V genes
print(f"\n=== 4. CROSS-PATIENT SPECIFICITY of each patient's top V gene ===")
spec=[]
for _,r in P.iterrows():
    gene=r.top_V_gene
    if not isinstance(gene,str): continue
    gi=V.index(gene); col=np.asarray(C[:,[gi]].todense()).ravel() if hasattr(C,'todense') else np.asarray(C[:,gi]).ravel()
    own = o.patient_id.astype(str)==str(r.patient)
    spec.append({'patient':r.patient,'top_V_gene':gene,
      'frac_own_patient_cells': round(float((col[own.values]>0).mean()),4),
      'frac_other_patients':    round(float((col[~own.values]>0).mean()),4),
      'n_other_patients_with_it': int(o.loc[(col>0)&(~own.values),'patient_id'].nunique()),
      'frac_donor_plasma': round(float((col[(o.sample_type=='normal_bm').values]>0).mean()),4)})
S=pd.DataFrame(spec); S.to_csv(OUT/'cross_patient_specificity.csv',index=False)
print(S.head(12).to_string(index=False))
print(f"\n  median own-patient detection of top V gene : {S.frac_own_patient_cells.median():.3f}")
print(f"  median detection in OTHER patients         : {S.frac_other_patients.median():.3f}")
print(f"  median n other patients carrying it        : {S.n_other_patients_with_it.median():.0f} of 49")
print("\nDONE", flush=True)
