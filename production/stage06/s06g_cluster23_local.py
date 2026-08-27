# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 06
# Step         : s06g_cluster23_local.py
# What it does : cluster-23 local subclustering
# Writes       : results/06_annotation/cluster23_local/*
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : heredoc-written script, run with `conda run -n mm-core python <path>`
#   Original path  : /tmp/claude-1000/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646/scratchpad/s06/c23_local.py  (session scratchpad, since deleted)
#   Executed (UTC) : 2026-08-25T10:57:54Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
"""Cluster 23 — local, un-integrated analysis. Diagnostic only; changes no global label."""
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from pathlib import Path
from sklearn.metrics import adjusted_rand_score as ari, normalized_mutual_info_score as nmi

OUT = Path('results/06_annotation/cluster23_local'); OUT.mkdir(parents=True, exist_ok=True)
sc.settings.verbosity = 1

# --- 1. subset from RAW COUNTS; no global embedding, graph or PCA is carried over ---
g = sc.read_h5ad('results/05_integration/integrated.h5ad')
m = (g.obs['leiden'].astype(str) == '23').to_numpy()
a = ad.AnnData(X=g.layers['counts'][m].copy(), obs=g.obs[m].copy(),
               var=pd.DataFrame(index=g.var_names))
del g
a.layers['counts'] = a.X.copy()
print(f"cluster 23 cells: {a.n_obs}, genes: {a.n_vars}", flush=True)

# --- 2. fresh local workflow ---
sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
sc.pp.highly_variable_genes(a, n_top_genes=2000, batch_key=None)
print("local HVGs:", int(a.var.highly_variable.sum()), flush=True)
sub = a[:, a.var.highly_variable].copy()
sc.pp.scale(sub, max_value=10); sc.tl.pca(sub, n_comps=30, svd_solver='arpack')
a.obsm['X_pca_local'] = sub.obsm['X_pca']; del sub
sc.pp.neighbors(a, n_neighbors=15, use_rep='X_pca_local')
sc.tl.umap(a)

# --- 3. resolution sweep, no cherry-picking ---
RES = [0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.5,2.0]
for r in RES:
    sc.tl.leiden(a, resolution=r, key_added=f'local_r{r}', flavor='igraph',
                 n_iterations=2, directed=False)
rows = []
for i, r in enumerate(RES):
    lab = a.obs[f'local_r{r}']
    row = {'resolution': r, 'n_clusters': lab.nunique(),
           'min_cluster_size': int(lab.value_counts().min()),
           'max_cluster_size': int(lab.value_counts().max())}
    if i:
        prev = a.obs[f'local_r{RES[i-1]}']
        row['ARI_vs_prev'] = round(ari(prev, lab), 3)
        row['NMI_vs_prev'] = round(nmi(prev, lab), 3)
    rows.append(row)
sweep = pd.DataFrame(rows); sweep.to_csv(OUT/'local_resolution_sweep.csv', index=False)
print(sweep.to_string(index=False), flush=True)

# --- 4. lineage axes, declared before interpretation ---
AX = {
 'T_ab':      ['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2'],
 'T_gd':      ['TRDC','TRGC1','TRGC2'],
 'NK':        ['KLRD1','KLRF1','NCAM1','FCGR3A','KLRC1'],
 'cytotoxic': ['NKG7','GNLY','PRF1','GZMB','GZMA','CTSW'],   # STATE, not lineage
}
MIN_GENES = 2   # same positive-evidence rule as the global framework
ev = {}
for name, gs in AX.items():
    gs = [x for x in gs if x in a.var_names]
    X = a[:, gs].layers['counts']
    X = np.asarray(X.todense()) if hasattr(X, 'todense') else np.asarray(X)
    ev[name] = (X > 0).sum(axis=1) >= MIN_GENES
    a.obs[f'ev_{name}'] = ev[name]
ev = pd.DataFrame(ev, index=a.obs_names)

def call(r):
    t, gd, nk = r['T_ab'], r['T_gd'], r['NK']
    if t and not nk:  return 'T_gd' if gd else 'T_ab'
    if nk and not t:  return 'NK'
    if t and nk:      return 'NKT_like_mixed'
    return 'unresolved'
a.obs['lineage_call'] = pd.Categorical([call(r) for _, r in ev.iterrows()])
print("\n=== cell-level lineage evidence ===", flush=True)
print(a.obs['lineage_call'].value_counts().to_string(), flush=True)
print("\ncross-tab T x NK:", flush=True)
print(pd.crosstab(ev['T_ab'], ev['NK']).to_string(), flush=True)
ev.assign(lineage_call=a.obs['lineage_call'].values).to_csv(OUT/'local_lineage_evidence.csv.gz',
                                                            compression='gzip')

# --- 5. per-local-cluster tables at a stable resolution chosen AFTER the sweep ---
RES_PICK = 0.3
key = f'local_r{RES_PICK}'
a.obs['local_cluster'] = a.obs[key]
det_genes = [x for gs in AX.values() for x in gs if x in a.var_names]
X = a[:, det_genes].layers['counts']
X = np.asarray(X.todense()) if hasattr(X, 'todense') else np.asarray(X)
det = pd.DataFrame(X > 0, columns=det_genes).groupby(a.obs['local_cluster'].values).mean()
det.round(4).to_csv(OUT/'local_marker_detection.csv')

tech = a.obs.groupby('local_cluster', observed=True).agg(
    n_cells=('local_cluster','size'),
    n_patients=('patient_id','nunique'), n_samples=('sample_name','nunique'),
    median_counts=('total_counts','median'), median_genes=('n_genes_by_counts','median'),
    median_pct_mt=('pct_counts_mt','median'))
tech['dominant_cohort'] = a.obs.groupby('local_cluster', observed=True)['cohort'].agg(lambda s: s.mode().iat[0])
tech['top_patient_frac'] = a.obs.groupby('local_cluster', observed=True)['patient_id'].agg(
    lambda s: round(s.value_counts(normalize=True).iat[0], 3))
for lc in ('T_ab','T_gd','NK','NKT_like_mixed','unresolved'):
    tech[f'frac_{lc}'] = a.obs.groupby('local_cluster', observed=True)['lineage_call'].agg(
        lambda s, lc=lc: round((s == lc).mean(), 3))
tech.to_csv(OUT/'local_technical_covariates.csv')
print(f"\n=== local clusters at resolution {RES_PICK} ===", flush=True)
print(tech.to_string(), flush=True)

pd.crosstab(a.obs['local_cluster'], a.obs['patient_id']).to_csv(OUT/'local_cluster_composition_by_patient.csv')
pd.crosstab(a.obs['local_cluster'], a.obs['sample_name']).to_csv(OUT/'local_cluster_composition_by_sample.csv')
pd.crosstab(a.obs['local_cluster'], a.obs['cohort']).to_csv(OUT/'local_cluster_composition_by_cohort.csv')
a.obs[['sample_name','patient_id','cohort','local_cluster','lineage_call',
       'total_counts','n_genes_by_counts']].to_csv(OUT/'local_cluster_assignments.csv.gz',
                                                   compression='gzip')

sc.tl.rank_genes_groups(a, 'local_cluster', method='wilcoxon')
pd.DataFrame({c: [x for x in a.uns['rank_genes_groups']['names'][c][:15]]
              for c in a.obs['local_cluster'].cat.categories}).T.to_csv(OUT/'local_de_markers.csv')

# patient-level representation per proposed lineage
pl = pd.crosstab(a.obs['lineage_call'], a.obs['patient_id'])
pd.DataFrame({'n_cells': a.obs['lineage_call'].value_counts(),
              'n_patients': (pl > 0).sum(axis=1),
              'n_samples': pd.crosstab(a.obs['lineage_call'], a.obs['sample_name']).gt(0).sum(axis=1),
              'top_patient_frac': pl.div(pl.sum(axis=1), axis=0).max(axis=1).round(3)}
             ).to_csv(OUT/'lineage_patient_representation.csv')
print("\n=== patient representation per lineage ===", flush=True)
print(pd.read_csv(OUT/'lineage_patient_representation.csv', index_col=0).to_string(), flush=True)

a.write_h5ad(OUT/'cluster23_local.h5ad', compression='gzip')
print("\nDONE", flush=True)
