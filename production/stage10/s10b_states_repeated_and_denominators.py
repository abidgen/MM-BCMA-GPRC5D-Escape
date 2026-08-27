# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 10
# Step         : s10b_states_repeated_and_denominators.py
# What it does : Level-1 state calls, repeated-sample structure, primary-vs-sensitivity coherence
# Writes       : results/10_dn_coherence/repeated_sample_dn_coherence.csv; dn_coherence_final_states.csv; primary_vs_sensitivity_coherence.csv; dn_program_summary_by_patient.csv; dn_de_by_patient.csv
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : heredoc-written script, run with `conda run -n mm-core python <path>`
#   Original path  : $S/s10b.py  (session scratchpad, since deleted)
#   Executed (UTC) : 2026-08-25T22:30:56Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp, warnings, pickle
from pathlib import Path
from mm_escape import subclone as SC
warnings.filterwarnings('ignore'); sc.settings.verbosity=0
OUT=Path('results/10_dn_coherence'); pd.set_option('display.width',280)
SEED=20260825; NPERM=1000
E,de_store=pickle.load(open('/tmp/s10.pkl','rb'))
A=sc.read_h5ad('results/08_dual_antigen_escape/antigen_states.h5ad')
A.obs['patient_id']=A.obs['patient_id'].astype(str); A.obs['sample_name']=A.obs['sample_name'].astype(str)
genes=list(A.var_names); C=sp.csr_matrix(A.layers['counts'])
keep=[i for i,g in enumerate(genes) if g not in set(SC.ANTIGEN_FEATURES)]
C=C[:,keep]; genes=[genes[i] for i in keep]; obs=A.obs.copy(); del A
is_dn_all=(obs.observed_state=='double_negative').values; depth_all=obs.depth_ex_antigen.values.astype(float)

def local_p(idx,is_dn,depth):
    ad=sc.AnnData(C[idx].copy()); ad.var_names=genes
    sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad)
    sc.pp.highly_variable_genes(ad,n_top_genes=min(2000,ad.n_vars-1))
    ad=ad[:,ad.var.highly_variable].copy()
    npc=int(min(30,ad.n_obs-1,ad.n_vars-1)); sc.pp.scale(ad,max_value=10)
    sc.tl.pca(ad,n_comps=npc,svd_solver='arpack',random_state=SEED)
    k=int(min(15,ad.n_obs-1)); sc.pp.neighbors(ad,n_neighbors=k,use_rep='X_pca',random_state=SEED)
    g=ad.obsp['connectivities']
    nn=np.asarray([np.argsort(-np.asarray(g[i].todense()).ravel())[:k] for i in range(ad.n_obs)])
    stat=lambda l: SC.knn_dn_fraction(l,nn)
    _,p=SC.depth_stratified_permutation(is_dn,SC.adaptive_depth_bins(depth),stat,NPERM,SEED)
    return float(p), float(stat(is_dn))

# ---- 13. repeated samples, each analysed separately first ----
rows=[]
for pat in sorted(obs.patient_id.unique()):
    pm=(obs.patient_id==pat).values&obs.in_primary.values.astype(bool)
    samples=sorted(obs.sample_name[pm].unique())
    if len(samples)<2: continue
    for s in samples:
        m=pm&(obs.sample_name==s).values; idx=np.flatnonzero(m)
        dn=is_dn_all[idx]
        r={'patient':pat,'sample':s,'n_cells':len(idx),'n_dn':int(dn.sum()),
           'n_antigen_pos':int((~dn).sum())}
        if r['n_dn']>=SC.MIN_GROUP_CELLS and r['n_antigen_pos']>=SC.MIN_GROUP_CELLS:
            p,f=local_p(idx,dn,depth_all[idx]); r['perm_p_depth_stratified']=p; r['knn_dn_frac']=f
            r['sample_evaluable']=True
        else:
            r['perm_p_depth_stratified']=np.nan; r['knn_dn_frac']=np.nan; r['sample_evaluable']=False
        rows.append(r)
RS=pd.DataFrame(rows)
def rstat(g):
    ev=g[g.sample_evaluable]
    if len(ev)<2: return 'not_evaluable'
    sig=(ev.perm_p_depth_stratified<0.05)
    return 'coherent_across_samples' if sig.all() else ('discordant' if sig.any() else 'no_support_any_sample')
RSS=RS.groupby('patient').apply(rstat).rename('repeated_sample_status').reset_index()
RS=RS.merge(RSS,on='patient'); RS.round(4).to_csv(OUT/'repeated_sample_dn_coherence.csv',index=False)
print("=== 7. REPEATED-SAMPLE DN COHERENCE (each sample analysed separately) ===")
print(RS.round(4).to_string(index=False))

# ---- final states ----
P=E[E.denominator=='primary'].set_index('patient'); Sx=E[E.denominator=='sensitivity'].set_index('patient')
rmap=dict(zip(RSS.patient,RSS.repeated_sample_status))
out=[]
for pat in sorted(E.patient.unique()):
    pr=P.loc[pat].to_dict() if pat in P.index else {'evaluable':False}
    se=Sx.loc[pat].to_dict() if pat in Sx.index else {'evaluable':False}
    rs=rmap.get(pat)
    state=SC.coherence_state(
        {'evaluable':bool(pr.get('evaluable')),'perm_p':pr.get('perm_p_depth_stratified',np.nan)},
        {'evaluable':bool(se.get('evaluable')),'perm_p':se.get('perm_p_depth_stratified',np.nan)},
        repeated_sample_status='discordant' if rs=='discordant' else None)
    out.append({'patient':pat,'cohort':obs.cohort[(obs.patient_id==pat).values].iloc[0],
        'dn_coherence_state':state,'cnv_subclone_state':SC.CNV_NOT_EVALUABLE,
        'evaluable_primary':bool(pr.get('evaluable')),'evaluable_sensitivity':bool(se.get('evaluable')),
        'n_dn_primary':pr.get('n_dn'),'n_antigen_pos_primary':pr.get('n_antigen_pos'),
        'perm_p_primary':pr.get('perm_p_depth_stratified'),'perm_p_sensitivity':se.get('perm_p_depth_stratified'),
        'perm_p_unconditioned_primary':pr.get('perm_p_unconditioned'),
        'knn_dn_frac_primary':pr.get('knn_dn_frac'),'dn_rate_primary':pr.get('dn_rate'),
        'morans_i_primary':pr.get('morans_i'),'depth_ratio_dn_over_pos':pr.get('ratio_depth_ex_antigen'),
        'repeated_sample_status':rs or 'single_sample','not_evaluable_reason':pr.get('reason','')})
F=pd.DataFrame(out); F.to_csv(OUT/'dn_coherence_final_states.csv',index=False)
F[['patient','cohort','evaluable_primary','evaluable_sensitivity','perm_p_primary',
   'perm_p_sensitivity','dn_coherence_state']].to_csv(OUT/'primary_vs_sensitivity_coherence.csv',index=False)
print("\n=== STATES ==="); print(F.dn_coherence_state.value_counts().to_string())
print(F[['patient','cohort','n_dn_primary','perm_p_primary','perm_p_sensitivity',
         'perm_p_unconditioned_primary','depth_ratio_dn_over_pos','repeated_sample_status',
         'dn_coherence_state']].round(4).to_string(index=False))

# ---- 6. cross-patient program consistency (non-antigen) ----
sup=set(F[F.dn_coherence_state==SC.SUPPORTED].patient)
lfc={k[0]:v for k,v in de_store.items() if k[1]=='primary' and k[0] in sup}
if len(lfc)>=3:
    M=pd.DataFrame(lfc).dropna(thresh=max(3,int(0.8*len(lfc))))
    cons=(np.sign(M).sum(axis=1).abs()/M.notna().sum(axis=1))
    top=M.assign(mean_lfc=M.mean(1),frac_consistent=cons,n_pat=M.notna().sum(1))
    top=top[(top.n_pat>=3)&(top.frac_consistent==1.0)].sort_values('mean_lfc',key=abs,ascending=False)
    print(f"\n=== 6. CROSS-PATIENT DE PROGRAM (supported patients, n={len(lfc)}) ===")
    print(f"  genes tested {len(M)}; unanimous-direction genes {len(top)}")
    print(top[['mean_lfc','n_pat']].head(25).round(3).to_string())
    top.to_csv(OUT/'dn_program_summary_by_patient.csv')
    pd.concat([v.rename(k) for k,v in lfc.items()],axis=1).to_csv(OUT/'dn_de_by_patient.csv')
else:
    pd.DataFrame().to_csv(OUT/'dn_program_summary_by_patient.csv')
    pd.concat([v.rename(k) for k,v in de_store.items() if k[1]=='primary'],axis=1).to_csv(OUT/'dn_de_by_patient.csv')
    print("\ncross-patient program: <3 supported patients with DE")
pickle.dump((F,RS),open('/tmp/s10b.pkl','wb'))
