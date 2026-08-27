# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 10
# Step         : s10e_pseudobulk_de_decoupler_tc.py
# What it does : cohort pseudobulk DE (patient as replicate), decoupler, TC-like subtype
# Writes       : results/10_dn_coherence/level2_program_cross_correlation.csv; pseudobulk_de_results.csv; decoupler_pathway_results.csv; tc_like_subtype.csv
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : heredoc-written script, run with `conda run -n mm-core python <path>`
#   Original path  : $S/s10e.py  (session scratchpad, since deleted)
#   Executed (UTC) : 2026-08-25T22:58:37Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
import numpy as np, pandas as pd, pickle, warnings, scipy.sparse as sp
from pathlib import Path
from scipy.stats import spearmanr
from mm_escape import subclone as SC, config, communication as CM
warnings.filterwarnings('ignore')
OUT=Path('results/10_dn_coherence'); pd.set_option('display.width',300); SEED=20260825
PP,R,repro,DS,PROGS,scores,obs,is_dn,depth,C,genes=pickle.load(open('/tmp/s10c.pkl','rb'))
L,elig,_=pickle.load(open('/tmp/s10d.pkl','rb'))

# ---- ALTERNATIVE EXPLANATION: is "three programs all down" one global factor? ----
pr=[p for p in sorted(repro)]
sub=PP[PP.denominator=='primary'].set_index('patient')
Mx=sub[[f'{k}_delta_matched' for k in pr]]
cm=Mx.corr(method='spearman')
print("=== ALTERNATIVE EXPLANATION CHECK: cross-program correlation of per-patient effects ===")
print(cm.round(3).to_string())
print(f"  patients with ALL {len(pr)} reproducible programs shifted DOWN: "
      f"{int((Mx<0).all(axis=1).sum())}/{len(Mx)}")
n_neg=(Mx<0).sum(axis=1)
print("  distribution of how many of the three are down per patient:", dict(n_neg.value_counts().sort_index()))
cm.round(4).to_csv(OUT/'level2_program_cross_correlation.csv')

# ---- per-patient rule vacuity diagnostic (reported, rule NOT retuned) ----
tot=int((L.level2_state!='DN_STATE_NOT_EVALUABLE').sum())
sup=int((L.level2_state=='DN_STATE_SUPPORTED').sum())
print(f"\n=== PER-PATIENT LEVEL-2 RULE DIAGNOSTIC ===")
print(f"  {sup} of {tot} evaluable patients satisfy the predeclared rule "
      f"({sup/tot:.0%}) -> the rule barely discriminates between patients.")
print("  The rule is NOT retuned after the fact; the cohort-level test is the real result.")

# ---- pseudobulk DE (paired, patient as blocking factor) ----
from scipy.stats import wilcoxon
de_rows=[]
for tag in ['primary','sensitivity']:
    dm=obs.in_primary.values.astype(bool) if tag=='primary' else obs.in_sensitivity.values.astype(bool)
    pats=elig[elig.denominator==tag].patient.tolist()
    DNm=[];POSm=[]
    for p in pats:
        m=dm&(obs.patient_id==p).values; idx=np.flatnonzero(m)
        dn=is_dn[idx]; bins=SC.adaptive_depth_bins(depth[idx])
        di,pi=SC.depth_matched_indices(dn,bins,SEED)
        DNm.append(np.asarray(C[idx[di]].sum(axis=0)).ravel())
        POSm.append(np.asarray(C[idx[pi]].sum(axis=0)).ravel())
    DNm=np.vstack(DNm); POSm=np.vstack(POSm)
    cpm_dn=1e6*DNm/DNm.sum(1,keepdims=True); cpm_pos=1e6*POSm/POSm.sum(1,keepdims=True)
    keep=((DNm>0).sum(0)>=len(pats)//2)&((POSm>0).sum(0)>=len(pats)//2)
    lfc=np.log2((cpm_dn[:,keep]+1)/(cpm_pos[:,keep]+1))
    W=np.array([wilcoxon(lfc[:,j])[1] if np.any(lfc[:,j]!=0) else 1.0 for j in range(lfc.shape[1])])
    gn=np.array(genes)[keep]
    d=pd.DataFrame({'gene':gn,'median_log2FC':np.median(lfc,0),'frac_up':(lfc>0).mean(0),
                    'p':W,'n_patients':len(pats),'denominator':tag})
    d['p_BH']=CM.benjamini_hochberg(d.p.values)
    de_rows.append(d)
DE=pd.concat(de_rows); DE.sort_values('p').round(5).to_csv(OUT/'pseudobulk_de_results.csv',index=False)
print(f"\n=== PSEUDOBULK DE (paired over patients, depth-matched cells) ===")
for tag in ['primary','sensitivity']:
    s=DE[DE.denominator==tag]
    print(f"  {tag}: {len(s)} genes tested, BH<0.05: {int((s.p_BH<0.05).sum())}, BH<0.10: {int((s.p_BH<0.10).sum())}")
both=set(DE[(DE.denominator=='primary')&(DE.p_BH<0.05)].gene)&set(DE[(DE.denominator=='sensitivity')&(DE.p_BH<0.05)].gene)
print(f"  significant under BOTH denominators: {len(both)}")
top=DE[(DE.denominator=='primary')&(DE.gene.isin(both))].nsmallest(15,'p')
print(top[['gene','median_log2FC','frac_up','p','p_BH']].round(4).to_string(index=False))

# ---- decoupler 2.x ----
try:
    import decoupler as dc
    print(f"\ndecoupler {dc.__version__}; API check dc.mt:", hasattr(dc,'mt'))
    try:
        net=dc.op.hallmark(organism='human')
        print("  hallmark resource fetched:",net.shape)
        dcs='evaluable'
    except Exception as e:
        dcs=f'NOT_EVALUABLE — resource unavailable offline: {type(e).__name__}'
        print("  ",dcs)
except Exception as e:
    dcs=f'NOT_EVALUABLE — {type(e).__name__}: {e}'
    print("\ndecoupler:",dcs)
pd.DataFrame([{'analysis':'decoupler Hallmark/PROGENy/CollecTRI','status':dcs,
   'note':'decoupler 2.x API (dc.mt.*/dc.op.*) only; weak pathway evidence is never converted into Level-2 support'}])\
  .to_csv(OUT/'decoupler_pathway_results.csv',index=False)

# ---- TC-like subtype, per patient, descriptive only ----
gi={g:i for i,g in enumerate(genes)}
rows=[]
for p in sorted(obs.patient_id.unique()):
    m=obs.in_primary.values.astype(bool)&(obs.patient_id==p).values
    if m.sum()<SC.MIN_PATIENT_CELLS: continue
    tot_=C[m].sum(); d={'patient':p,'n_cells':int(m.sum())}
    for cls,gl in config.TC_GENES.items():
        vals=[1e6*C[m][:,gi[g]].sum()/tot_ for g in gl if g in gi]
        d[cls]=float(np.max(vals)) if vals else np.nan
    d['CKS1B_cpm']=float(1e6*C[m][:,gi[config.TC_1Q21_GENE]].sum()/tot_) if config.TC_1Q21_GENE in gi else np.nan
    cl={k:d[k] for k in config.TC_GENES if np.isfinite(d.get(k,np.nan))}
    d['TC_like_subtype']=max(cl,key=cl.get) if cl else 'NOT_EVALUABLE'
    d['top_cpm']=cl.get(d['TC_like_subtype'],np.nan)
    rows.append(d)
TC=pd.DataFrame(rows); TC.round(3).to_csv(OUT/'tc_like_subtype.csv',index=False)
print("\n=== TC-LIKE EXPRESSION SUBTYPE (per patient, DESCRIPTIVE ONLY, never a translocation call) ===")
print(TC.TC_like_subtype.value_counts().to_string())
print(TC[['patient','TC_like_subtype','top_cpm','CKS1B_cpm']].round(1).to_string(index=False))
pickle.dump((DE,TC,cm),open('/tmp/s10e.pkl','wb'))
