# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 07
# Step         : s07d_donor_cnv_calibration.py
# What it does : infercnvpy donor negative-control calibration -> CNV rejected as NOT_EVALUABLE
# Writes       : results/07_malignant_plasma/cnv_negative_result/donor_cnv_calibration.csv; donor_cnv_cells_<sample>.csv.gz
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : heredoc-written script, run with `conda run -n mm-core python <path>`
#   Original path  : /tmp/claude-1000/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646/scratchpad/s06/donor_cnv.py  (session scratchpad, since deleted)
#   Executed (UTC) : 2026-08-25T11:53:30Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
"""Donor-ONLY CNV calibration. No disease CNV is computed or inspected."""
import numpy as np, pandas as pd, scanpy as sc, infercnvpy as cnv
from pathlib import Path
OUT = Path('results/07_malignant_plasma'); pd.set_option('display.width',240)

genes = pd.read_csv(OUT/'cnv_input_gene_set.tsv.gz', sep='\t')
keep = set(genes.symbol)
assert not ({'TNFRSF17','GPRC5D'} & keep), "ANTIGEN LEAK"
assert not any(g.startswith(('IGHV','IGKV','IGLV','IGHJ','IGKJ','IGLJ','IGHG','IGKC','IGLC')) for g in keep), "IG LEAK"

g = sc.read_h5ad('results/05_integration/integrated.h5ad')
lab = pd.read_csv('results/06_annotation/per_cell_labels.csv.gz', index_col=0)
g.obs['cell_type'] = lab.reindex(g.obs_names)['cell_type'].astype(str).values
donors = g.obs.loc[g.obs['sample_type'].astype(str)=='normal_bm','sample_name'].unique()
print("donor samples:", list(donors), flush=True)

pos = genes.set_index('symbol')
rows = []
for s in sorted(donors):
    m = (g.obs['sample_name'].astype(str)==s) & g.obs['cell_type'].isin(['PlasmaCell','Tcell','Myeloid'])
    a = g[m.values].copy()
    a = a[:, [x for x in a.var_names if x in keep]].copy()
    a.var['chromosome'] = ['chr'+str(pos.loc[x,'chromosome']) for x in a.var_names]
    a.var['start'] = [int(pos.loc[x,'start']) for x in a.var_names]
    a.var['end']   = [int(pos.loc[x,'end']) for x in a.var_names]
    a = a[:, a.var.sort_values(['chromosome','start']).index].copy()
    npl = int((a.obs.cell_type=='PlasmaCell').sum()); nref = int(a.obs.cell_type.isin(['Tcell','Myeloid']).sum())
    if npl == 0 or nref < 200:
        print(f"  {s}: SKIP (plasma {npl}, ref {nref})", flush=True); continue
    cnv.tl.infercnv(a, reference_key='cell_type', reference_cat=['Tcell','Myeloid'],
                    window_size=100, step=10)
    X = a.obsm['X_cnv']
    X = np.asarray(X.todense()) if hasattr(X,'todense') else np.asarray(X)
    burden = np.abs(X).mean(axis=1)                      # per-cell CNV burden
    a.obs['cnv_burden'] = burden
    ref = burden[a.obs.cell_type.isin(['Tcell','Myeloid']).values]
    pl  = burden[(a.obs.cell_type=='PlasmaCell').values]
    center = float(np.median(ref)); mad = float(np.median(np.abs(ref-center)))
    scale = 1.4826*mad if mad > 0 else np.nan
    z_pl = (pl-center)/scale; z_ref = (ref-center)/scale
    rows.append({'donor':s,'n_plasma':npl,'n_ref':nref,'n_T':int((a.obs.cell_type=='Tcell').sum()),
        'n_Mye':int((a.obs.cell_type=='Myeloid').sum()),
        'ref_center':center,'ref_mad_scale':scale,'n_windows':X.shape[1],
        'plasma_z_median':float(np.median(z_pl)),'plasma_z_p90':float(np.percentile(z_pl,90)),
        'plasma_z_p95':float(np.percentile(z_pl,95)),'plasma_z_max':float(z_pl.max()),
        'ref_z_p95':float(np.percentile(z_ref,95)),'ref_z_p99':float(np.percentile(z_ref,99)),
        'frac_plasma_z_gt3':float((z_pl>3).mean()),'frac_plasma_z_gt5':float((z_pl>5).mean()),
        'frac_ref_z_gt3':float((z_ref>3).mean())})
    pd.DataFrame({'cell':a.obs_names,'cell_type':a.obs.cell_type.values,'burden':burden,
                  'z':(burden-center)/scale,'donor':s}).to_csv(
        OUT/f'donor_cnv_cells_{s}.csv.gz', index=False, compression='gzip')
    print(f"  {s}: plasma {npl} ref {nref} windows {X.shape[1]} | "
          f"plasma z med {np.median(z_pl):+.2f} p95 {np.percentile(z_pl,95):+.2f} | "
          f"ref z p95 {np.percentile(z_ref,95):+.2f} | frac plasma z>3 {100*(z_pl>3).mean():.1f}%", flush=True)

d = pd.DataFrame(rows); d.to_csv(OUT/'donor_cnv_calibration.csv', index=False)
print("\n=== DONOR CNV CALIBRATION ==="); print(d.round(3).to_string(index=False), flush=True)
print("\nreference lineage skew vs plasma z_p95:")
d['skew'] = d[['n_T','n_Mye']].max(axis=1)/d.n_ref
print(d[['donor','n_T','n_Mye','skew','plasma_z_median','plasma_z_p95','frac_plasma_z_gt3']].round(3).to_string(index=False))
print("\nleave-one-donor-out max plasma_z_p95:")
for s in d.donor:
    print(f"  drop {s:12s} -> max {d[d.donor!=s].plasma_z_p95.max():+.3f}")
print("\nDONE", flush=True)
