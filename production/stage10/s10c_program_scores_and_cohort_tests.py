# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 10
# Step         : s10c_program_scores_and_cohort_tests.py
# What it does : frozen seven program scores, depth association, cohort-level program tests
# Writes       : results/10_dn_coherence/program_score_vs_depth.csv; dn_program_scores_by_patient.csv; level2_program_cohort_tests.csv
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : heredoc-written script, run with `conda run -n mm-core python <path>`
#   Original path  : $S/s10c.py  (session scratchpad, since deleted)
#   Executed (UTC) : 2026-08-25T22:56:53Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp, warnings, pickle
from pathlib import Path
from scipy.stats import spearmanr, wilcoxon
from mm_escape import subclone as SC, config, communication as CM
warnings.filterwarnings('ignore'); sc.settings.verbosity=0
OUT=Path('results/10_dn_coherence'); pd.set_option('display.width',300); SEED=20260825

A=sc.read_h5ad('results/08_dual_antigen_escape/antigen_states.h5ad')
A.obs['patient_id']=A.obs['patient_id'].astype(str); A.obs['sample_name']=A.obs['sample_name'].astype(str)
genes=list(A.var_names); C=sp.csr_matrix(A.layers['counts'])
C,genes=SC.drop_antigen_features(C,genes)          # antigens gone before anything else
assert not (set(genes)&set(SC.ANTIGEN_FEATURES))
obs=A.obs.copy(); del A
is_dn=(obs.observed_state=='double_negative').values
depth=obs.depth_ex_antigen.values.astype(float)

ad=sc.AnnData(C.copy()); ad.var_names=genes; ad.obs=obs.reset_index(drop=True)
sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad)
PROGS={}
for name in config.LEVEL2_PROGRAMS:
    gs=[g for g in config.STATE_PROGRAMS[name] if g in set(genes)]
    miss=[g for g in config.STATE_PROGRAMS[name] if g not in set(genes)]
    sc.tl.score_genes(ad,gs,score_name=name,random_state=SEED)
    PROGS[name]=(gs,miss)
    print(f"{name:22s} {len(gs)}/{len(config.STATE_PROGRAMS[name])} genes"+(f"  MISSING {miss}" if miss else ""))
scores=ad.obs[list(config.LEVEL2_PROGRAMS)].values

# ---- MANDATORY depth screen, BEFORE any DN-vs-comparator difference ----
ds=[]
for i,name in enumerate(config.LEVEL2_PROGRAMS):
    rho,p=spearmanr(scores[:,i],depth)
    ds.append({'program':name,'n_cells':len(depth),'spearman_rho_vs_depth':round(float(rho),4),
               'p':float(p),'strong_depth_tracking':bool(abs(rho)>=0.3)})
DS=pd.DataFrame(ds); DS.to_csv(OUT/'program_score_vs_depth.csv',index=False)
print("\n=== DEPTH SCREEN (raw program score vs depth_ex_antigen, all clone cells pooled) ===")
print(DS.to_string(index=False))
print("  MYC and OXPHOS are reported here FIRST by requirement; any strong tracker is read\n"
      "  only through the depth-matched lens below.")

# ---- per-patient DN vs comparator, depth-matched and unmatched, both denominators ----
EV=pd.read_csv('results/10_dn_coherence/stage10_evaluability.csv')
ok=EV[EV.evaluable].set_index(['patient','denominator'])
rows=[]
for (pat,tag),_ in ok.iterrows():
    dm=obs.in_primary.values.astype(bool) if tag=='primary' else obs.in_sensitivity.values.astype(bool)
    m=dm&(obs.patient_id==pat).values; idx=np.flatnonzero(m)
    dn=is_dn[idx]; bins=SC.adaptive_depth_bins(depth[idx])
    di,pi=SC.depth_matched_indices(dn,bins,SEED)
    r={'patient':pat,'denominator':tag,'n_dn':int(dn.sum()),'n_pos':int((~dn).sum()),
       'n_matched_per_group':int(di.size),
       'depth_ratio_raw':float(np.median(depth[idx][dn])/max(np.median(depth[idx][~dn]),1e-9)),
       'depth_ratio_matched':float(np.median(depth[idx][di])/max(np.median(depth[idx][pi]),1e-9))}
    for i,name in enumerate(config.LEVEL2_PROGRAMS):
        s=scores[idx,i]
        r[f'{name}_delta_raw']=float(s[dn].mean()-s[~dn].mean())
        r[f'{name}_delta_matched']=float(s[di].mean()-s[pi].mean()) if di.size else np.nan
    rows.append(r)
PP=pd.DataFrame(rows); PP.round(5).to_csv(OUT/'dn_program_scores_by_patient.csv',index=False)
print(f"\npatient x denominator rows: {len(PP)}; median matched cells/group {PP.n_matched_per_group.median():.0f}")
print(f"depth ratio DN/pos: raw median {PP.depth_ratio_raw.median():.3f} -> matched median {PP.depth_ratio_matched.median():.3f}")

# ---- cohort-level: Wilcoxon signed-rank over PATIENTS, BH across the 7 programs ----
res=[]
for tag in ['primary','sensitivity']:
    sub=PP[PP.denominator==tag]
    for name in config.LEVEL2_PROGRAMS:
        for kind in ['matched','raw']:
            v=sub[f'{name}_delta_{kind}'].dropna().values
            if len(v)<6: res.append({'denominator':tag,'program':name,'effect':kind,'n_patients':len(v),
                'median_delta':np.nan,'W':np.nan,'p':np.nan}); continue
            W,p=wilcoxon(v)
            res.append({'denominator':tag,'program':name,'effect':kind,'n_patients':len(v),
                'median_delta':float(np.median(v)),'frac_positive':float((v>0).mean()),
                'W':float(W),'p':float(p)})
R=pd.DataFrame(res)
for (tag,kind),g in R.groupby(['denominator','effect']):
    R.loc[g.index,'p_BH']=CM.benjamini_hochberg(g.p.values)
R.round(5).to_csv(OUT/'level2_program_cohort_tests.csv',index=False)
print("\n=== LEVEL-2 COHORT TESTS (Wilcoxon signed-rank over patients; BH across 7 programs) ===")
print(R[R.effect=='matched'].round(4).to_string(index=False))
print("\n  (unmatched, for the depth comparison only):")
print(R[R.effect=='raw'][['denominator','program','median_delta','p','p_BH']].round(4).to_string(index=False))

repro=set()
for name in config.LEVEL2_PROGRAMS:
    a=R[(R.program==name)&(R.effect=='matched')&(R.denominator=='primary')].iloc[0]
    b=R[(R.program==name)&(R.effect=='matched')&(R.denominator=='sensitivity')].iloc[0]
    if (a.p_BH<0.10 and b.p_BH<0.10 and np.sign(a.median_delta)==np.sign(b.median_delta)):
        repro.add(name)
print(f"\n  COHORT-REPRODUCIBLE programs (BH<0.10 both denominators, consistent sign): "
      f"{sorted(repro) if repro else 'NONE'}")
pickle.dump((PP,R,repro,DS,PROGS,scores,obs,is_dn,depth,C,genes),open('/tmp/s10c.pkl','wb'))
