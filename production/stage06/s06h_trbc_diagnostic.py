# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 06
# Step         : s06h_trbc_diagnostic.py
# What it does : TRBC1/2 context diagnostic
# Writes       : results/06_annotation/cluster23_local/trbc_diagnostic/*
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : heredoc-written script, run with `conda run -n mm-core python <path>`
#   Original path  : /tmp/claude-1000/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646/scratchpad/s06/trbc_diag.py  (session scratchpad, since deleted)
#   Executed (UTC) : 2026-08-25T11:06:36Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
"""TRBC1/TRBC2 diagnostic in Leiden 23. Interpretive only; changes no label."""
import numpy as np, pandas as pd, scanpy as sc
from pathlib import Path
from scipy.stats import spearmanr, mannwhitneyu

OUT = Path('results/06_annotation/cluster23_local/trbc_diagnostic')
a = sc.read_h5ad('results/06_annotation/cluster23_local/cluster23_local.h5ad')
print("cells:", a.n_obs, "| obs has doublet info:",
      [c for c in a.obs.columns if 'doublet' in c.lower()], flush=True)

T   = ['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2']
GD  = ['TRDC','TRGC1','TRGC2']
NK  = ['KLRD1','KLRF1','NCAM1','FCGR3A','KLRC1']
CYT = ['NKG7','GNLY','PRF1','GZMB','GZMA','CTSW']
ALL = T + GD + NK + CYT
X = a[:, ALL].layers['counts']
X = np.asarray(X.todense()) if hasattr(X,'todense') else np.asarray(X)
c = pd.DataFrame(X, columns=ALL, index=a.obs_names)          # RAW UMI
o = a.obs
call = o['lineage_call'].astype(str)

# ---- 3. TRBC intensity by category -------------------------------------------------
rows=[]
for g in ['TRBC1','TRBC2']:
    for k in ['NK','NKT_like_mixed','T_ab','T_gd','unresolved']:
        v = c.loc[(call==k).values, g]; pos = v[v>0]
        rows.append({'gene':g,'category':k,'n_cells':len(v),
            'frac_detected':round((v>0).mean(),4),
            'median_umi_pos': float(pos.median()) if len(pos) else np.nan,
            'mean_umi_pos': round(float(pos.mean()),3) if len(pos) else np.nan,
            **{f'p{q}_umi_pos': float(np.percentile(pos,q)) if len(pos) else np.nan
               for q in (25,75,90,95)},
            'frac_pos_exactly_1': round(float((pos==1).mean()),4) if len(pos) else np.nan,
            'frac_all_ge2': round(float((v>=2).mean()),4),
            'frac_all_ge3': round(float((v>=3).mean()),4),
            'frac_all_ge5': round(float((v>=5).mean()),4)})
summ = pd.DataFrame(rows); summ.to_csv(OUT/'trbc_raw_count_summary.csv', index=False)
print("\n=== TRBC intensity by category ===", flush=True)
print(summ[['gene','category','n_cells','frac_detected','median_umi_pos','frac_pos_exactly_1',
            'frac_all_ge2','frac_all_ge3']].to_string(index=False), flush=True)

# ---- 4. coordinated evidence among TRBC+ cells --------------------------------------
trbc_pos = (c[['TRBC1','TRBC2']] > 0).any(axis=1)
anycd3 = (c[['CD3D','CD3E','CD3G']] > 0).any(axis=1)
cd3_2  = (c[['CD3D','CD3E','CD3G']] > 0).sum(axis=1) >= 2
trac   = c['TRAC'] > 0
gdprog = (c[GD] > 0).sum(axis=1) >= 2
rows=[]
for k in ['NK','NKT_like_mixed','T_ab','T_gd','unresolved']:
    m = (call==k).values & trbc_pos.values
    n = int(m.sum())
    if not n: continue
    rows.append({'category':k,'n_TRBC_pos':n,
      'TRBC_only_no_CD3_no_TRAC': round(float((~anycd3[m] & ~trac[m]).mean()),4),
      'TRBC_plus_anyCD3': round(float(anycd3[m].mean()),4),
      'TRBC_plus_2CD3': round(float(cd3_2[m].mean()),4),
      'TRBC_plus_TRAC': round(float(trac[m].mean()),4),
      'TRBC_plus_CD3_and_TRAC': round(float((anycd3[m] & trac[m]).mean()),4),
      'TRBC_plus_gd_program': round(float(gdprog[m].mean()),4)})
co = pd.DataFrame(rows); co.to_csv(OUT/'trbc_coexpression_patterns.csv', index=False)
print("\n=== coordinated T evidence among TRBC+ cells ===", flush=True)
print(co.to_string(index=False), flush=True)

# ---- 5. the decisive population: strong NK, no CD3, no TRAC -------------------------
nk_strong = (c[NK] > 0).sum(axis=1) >= 2
clean_nk = nk_strong & ~anycd3 & ~trac
sel = clean_nk & trbc_pos
tot = c.loc[sel.values, ['TRBC1','TRBC2']].max(axis=1)
other = c.loc[(trbc_pos & anycd3 & trac).values, ['TRBC1','TRBC2']].max(axis=1)
u = mannwhitneyu(tot, other, alternative='less') if len(tot) and len(other) else None
print(f"\n=== strong-NK / CD3- / TRAC- cells: n={int(clean_nk.sum())}, "
      f"TRBC+ among them={int(sel.sum())} ({100*sel.sum()/max(clean_nk.sum(),1):.1f}%) ===", flush=True)
print(f"  their max TRBC UMI: median {tot.median():.0f}, mean {tot.mean():.2f}, "
      f"frac==1 {100*(tot==1).mean():.1f}%, frac>=3 {100*(tot>=3).mean():.1f}%", flush=True)
print(f"  coordinated (CD3+ and TRAC+) TRBC+ cells: median {other.median():.0f}, "
      f"mean {other.mean():.2f}, frac==1 {100*(other==1).mean():.1f}%, frac>=3 {100*(other>=3).mean():.1f}%", flush=True)
if u: print(f"  Mann-Whitney (clean-NK < coordinated): p={u.pvalue:.3g}", flush=True)

pd.DataFrame({'group':['cleanNK_TRBCpos','coordinated_TRBCpos'],
  'n':[len(tot),len(other)],'median_umi':[tot.median(),other.median()],
  'mean_umi':[round(tot.mean(),3),round(other.mean(),3)],
  'frac_eq1':[round((tot==1).mean(),4),round((other==1).mean(),4)],
  'frac_ge3':[round((tot>=3).mean(),4),round((other>=3).mean(),4)]}
  ).to_csv(OUT/'trbc_by_lineage_call.csv', index=False)

# ---- 6. ambient pattern: depth + sample T-cell abundance ----------------------------
trbc_tot = c[['TRBC1','TRBC2']].sum(axis=1)
r_d, p_d = spearmanr(o['total_counts'], trbc_tot)
r_g, p_g = spearmanr(o['n_genes_by_counts'], trbc_tot)
print(f"\n=== depth association ===\n  Spearman TRBC UMI vs total_counts  r={r_d:.3f} p={p_d:.2g}"
      f"\n  Spearman TRBC UMI vs n_genes       r={r_g:.3f} p={p_g:.2g}", flush=True)

glob = pd.read_csv('results/06_annotation_c2d_accepted/per_cell_labels.csv.gz', index_col=0)
tfrac = (glob.groupby('sample_name')['cell_type']
         .agg(lambda s: (s=='Tcell').mean()).rename('sample_T_fraction'))
loc = pd.DataFrame({'sample_name':o['sample_name'].astype(str).values,
                    'patient_id':o['patient_id'].astype(str).values,
                    'cohort':o['cohort'].astype(str).values,
                    'clean_nk':clean_nk.values,'trbc_pos':trbc_pos.values})
bys = (loc[loc.clean_nk].groupby('sample_name')
       .agg(n_cleanNK=('trbc_pos','size'), frac_TRBCpos=('trbc_pos','mean')))
bys = bys.join(tfrac).dropna()
bys = bys[bys.n_cleanNK >= 20]
r_t, p_t = spearmanr(bys.sample_T_fraction, bys.frac_TRBCpos)
bys.round(4).to_csv(OUT/'trbc_vs_tcell_abundance.csv')
bys.round(4).to_csv(OUT/'trbc_by_sample.csv')
print(f"\n=== ambient check: sample T-cell fraction vs TRBC+ rate in clean-NK cells ===", flush=True)
print(f"  n samples (>=20 clean-NK cells) = {len(bys)}; Spearman r={r_t:.3f} p={p_t:.2g}", flush=True)

byp = (loc[loc.clean_nk].groupby('patient_id')
       .agg(n_cleanNK=('trbc_pos','size'), frac_TRBCpos=('trbc_pos','mean')))
byp.round(4).to_csv(OUT/'trbc_by_patient.csv')

# ---- 7. doublet / depth ------------------------------------------------------------
dcols = [x for x in ['doublet_score','is_doublet','doublet_class'] if x in o.columns]
if dcols:
    dd = o.groupby(call, observed=True).agg(
        n=('total_counts','size'), median_counts=('total_counts','median'),
        median_genes=('n_genes_by_counts','median'),
        median_doublet_score=('doublet_score','median'))
    dd.to_csv(OUT/'trbc_doublet_depth_diagnostic.csv')
    print("\n=== doublet / depth by category ===", flush=True); print(dd.to_string(), flush=True)
else:
    print("\nNO doublet metric retained in the stage-05 object", flush=True)

# ---- 9. gamma-delta ----------------------------------------------------------------
trdc_hi = c['TRDC'] > 0
gd = pd.DataFrame({
 'group':['TRDC_positive','TRDC_neg'],
 'n':[int(trdc_hi.sum()), int((~trdc_hi).sum())],
 'median_TRDC_umi_pos':[float(c.loc[trdc_hi,'TRDC'].median()), 0.0],
 'frac_TRDC_eq1':[round(float((c.loc[trdc_hi,'TRDC']==1).mean()),4), np.nan],
 'frac_TRGC1_pos':[round(float((c.loc[trdc_hi,'TRGC1']>0).mean()),4),
                   round(float((c.loc[~trdc_hi,'TRGC1']>0).mean()),4)],
 'frac_TRGC2_pos':[round(float((c.loc[trdc_hi,'TRGC2']>0).mean()),4),
                   round(float((c.loc[~trdc_hi,'TRGC2']>0).mean()),4)],
 'frac_anyCD3':[round(float(anycd3[trdc_hi.values].mean()),4),
                round(float(anycd3[(~trdc_hi).values].mean()),4)],
 'frac_TRAC':[round(float(trac[trdc_hi.values].mean()),4),
              round(float(trac[(~trdc_hi).values].mean()),4)],
 'frac_NKstrong':[round(float(nk_strong[trdc_hi.values].mean()),4),
                  round(float(nk_strong[(~trdc_hi).values].mean()),4)]})
gd.to_csv(OUT/'gamma_delta_diagnostic.csv', index=False)
print("\n=== gamma-delta diagnostic ===", flush=True); print(gd.to_string(index=False), flush=True)

# ---- 10. sensitivity: T evidence = CD3/TRAC only (REPORT ONLY) ----------------------
t_strict = (c[['CD3D','CD3E','CD3G','TRAC']] > 0).sum(axis=1) >= 2
def scall(i):
    t, nk_, g_ = t_strict.iloc[i], nk_strong.iloc[i], gdprog.iloc[i]
    if t and not nk_: return 'T_gd' if g_ else 'T_ab'
    if nk_ and not t: return 'NK'
    if t and nk_:     return 'NKT_like_mixed'
    return 'unresolved'
sens = pd.Series([scall(i) for i in range(len(c))], index=c.index, name='sensitivity_call')
sc_tab = pd.crosstab(call.values, sens.values)
pd.DataFrame({'frozen':call.value_counts(),'sensitivity':sens.value_counts()}).fillna(0).astype(int)\
  .to_csv(OUT/'trbc_sensitivity_calls.csv')
sc_tab.to_csv(OUT/'trbc_sensitivity_crosstab.csv')
print("\n=== sensitivity (CD3/TRAC only) — REPORT ONLY ===", flush=True)
print(pd.DataFrame({'frozen':call.value_counts(),'sensitivity':sens.value_counts()}).fillna(0).astype(int).to_string(), flush=True)
print("\ncrosstab frozen x sensitivity:", flush=True); print(sc_tab.to_string(), flush=True)
moved = (call.values=='NKT_like_mixed') & (sens.values=='NK')
print(f"\nmixed->NK movers: {int(moved.sum())}; their max TRBC UMI median "
      f"{c.loc[moved,['TRBC1','TRBC2']].max(axis=1).median():.0f}, "
      f"frac==1 {100*(c.loc[moved,['TRBC1','TRBC2']].max(axis=1)==1).mean():.1f}%", flush=True)

# ---- 11. patient robustness of the clean NK core -----------------------------------
for name, mask in [('frozen_NK', (call=='NK').values), ('sensitivity_NK', (sens=='NK').values)]:
    s = o[mask]
    print(f"\n{name}: n={mask.sum()}, patients={s['patient_id'].nunique()}, "
          f"samples={s['sample_name'].nunique()}, "
          f"top_patient_frac={s['patient_id'].value_counts(normalize=True).iat[0]:.3f}", flush=True)
    print("  cohort:", s['cohort'].value_counts().to_dict(), flush=True)
print("\nDONE", flush=True)
