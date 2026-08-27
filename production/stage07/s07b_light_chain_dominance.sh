# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 07
# Step         : s07b_light_chain_dominance.sh
# What it does : ratio-based kappa/lambda restriction; donor calibration
# Writes       : results/07_malignant_plasma/plasma_light_chain_percell.csv.gz; patient_light_chain_D.csv; donor_light_chain_calibration.csv
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : inline shell/python block
#   Original path  : (inline heredoc executed directly from the shell)
#   Executed (UTC) : 2026-08-25T11:30:58Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
cd /media/wrath/CART_mm_dual_antigen
/home/abid/miniforge3/envs/mm-core/bin/python << 'EOF' 2>&1 | grep -v Warning
import scanpy as sc, numpy as np, pandas as pd
pd.set_option('display.width',220)
g = sc.read_h5ad('results/05_integration/integrated.h5ad')
lab = pd.read_csv('results/06_annotation/per_cell_labels.csv.gz', index_col=0)
g.obs['cell_type'] = lab.reindex(g.obs_names)['cell_type'].values
pl = g[(g.obs['cell_type'].astype(str)=='PlasmaCell')].copy()
del g
V = list(pl.var_names)
KAPPA = ['IGKC'] + [x for x in V if x.startswith(('IGKV','IGKJ'))]
LAMBDA = [x for x in V if x in ('IGLC1','IGLC2','IGLC3','IGLC7')] + [x for x in V if x.startswith(('IGLV','IGLJ'))]
assert not ({'TNFRSF17','GPRC5D'} & set(KAPPA+LAMBDA)), "ANTIGEN LEAK"
C = pl.layers['counts']
def tot(gs):
    idx=[V.index(x) for x in gs]; X=C[:,idx]
    return np.asarray(X.todense()).sum(axis=1) if hasattr(X,'todense') else np.asarray(X).sum(axis=1)
k,l = tot(KAPPA), tot(LAMBDA); lc = k+l
EV=3; ev = lc>=EV
fk = np.divide(k, np.maximum(lc,1))
call = np.where(~ev,'insufficient', np.where(fk>=0.8,'kappa', np.where(fk<=0.2,'lambda','ambiguous')))
o = pl.obs[['sample_name','patient_id','cohort','sample_type']].copy()
o['lc_call']=call; o['lc_umi']=lc
print(f"RAW counts. plasma {len(o)} | evaluable(>={EV} UMI) {int(ev.sum())} ({100*ev.mean():.1f}%) "
      f"| called {int((o.lc_call.isin(['kappa','lambda'])).sum())}")
o.to_csv('results/07_malignant_plasma/plasma_light_chain_percell.csv.gz', compression='gzip')

d = o[o.sample_type.astype(str)=='normal_bm']
rows=[]
for s,gg in d.groupby('sample_name', observed=True):
    c = gg[gg.lc_call.isin(['kappa','lambda'])]
    D = max((c.lc_call=='kappa').mean(), 1-(c.lc_call=='kappa').mean()) if len(c) else np.nan
    rows.append({'donor':s,'n_plasma':len(gg),'n_evaluable':int((gg.lc_umi>=EV).sum()),
                 'n_called':len(c),'D':round(D,3) if len(c) else np.nan,
                 'median_lc_umi':float(gg.lc_umi.median())})
dd=pd.DataFrame(rows).sort_values('n_plasma',ascending=False)
dd.to_csv('results/07_malignant_plasma/donor_light_chain_calibration.csv',index=False)
print("\n=== DONOR-LEVEL D (raw counts; donors NOT pooled) ===")
print(dd.to_string(index=False))
q=dd[dd.n_called>=20]
print(f"\n  qualifying donors (>=20 called): {len(q)}  D range {q.D.min():.3f}-{q.D.max():.3f}  median {q.D.median():.3f}")
print("  leave-one-donor-out max D:", {s: float(q[q.donor!=s].D.max()) for s in q.donor})

print("\n=== patient-level D across all 50 patients (for context only, thresholds NOT set from this) ===")
pr=[]
for p,gg in o.groupby('patient_id', observed=True):
    c=gg[gg.lc_call.isin(['kappa','lambda'])]
    D=max((c.lc_call=='kappa').mean(),1-(c.lc_call=='kappa').mean()) if len(c) else np.nan
    pr.append({'patient':p,'n_plasma':len(gg),'n_called':len(c),'D':D,
               'sample_type':gg.sample_type.iloc[0]})
pp=pd.DataFrame(pr)
pp.to_csv('results/07_malignant_plasma/patient_light_chain_D.csv',index=False)
ok=pp[pp.n_called>=20]
print(f"  patients with >=20 called cells: {len(ok)} of {len(pp)}")
print(f"  D distribution: min {ok.D.min():.3f} q25 {ok.D.quantile(.25):.3f} median {ok.D.median():.3f} "
      f"q75 {ok.D.quantile(.75):.3f} max {ok.D.max():.3f}")
print(f"  myeloma patients D>=0.90: {int((ok[ok.sample_type!='normal_bm'].D>=0.90).sum())} of {len(ok[ok.sample_type!='normal_bm'])}")
print(f"  donors D>=0.75: {int((ok[ok.sample_type=='normal_bm'].D>=0.75).sum())} of {len(ok[ok.sample_type=='normal_bm'])}")
EOF
