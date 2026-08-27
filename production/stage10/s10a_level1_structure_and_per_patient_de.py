# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 10
# Step         : s10a_level1_structure_and_per_patient_de.py
# What it does : patient-local unintegrated embedding, Moran I, depth-stratified permutation, per-patient depth-matched DE
# Writes       : results/10_dn_coherence/stage10_evaluability.csv; dn_local_structure_by_patient.csv; dn_vs_nondn_depth_diagnostic.csv; de_per_patient/_de_<patient>_<denominator>.csv
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : heredoc-written script, run with `conda run -n mm-core python <path>`
#   Original path  : $S/s10.py  (session scratchpad, since deleted)
#   Executed (UTC) : 2026-08-25T22:29:27Z
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
warnings.filterwarnings('ignore')
sc.settings.verbosity=0
OUT=Path('results/10_dn_coherence'); pd.set_option('display.width',280)
SEED=20260825; NPERM=1000

A=sc.read_h5ad('results/08_dual_antigen_escape/antigen_states.h5ad')
A.obs['patient_id']=A.obs['patient_id'].astype(str); A.obs['sample_name']=A.obs['sample_name'].astype(str)
genes=list(A.var_names); C=sp.csr_matrix(A.layers['counts'])
keep=[i for i,g in enumerate(genes) if g not in set(SC.ANTIGEN_FEATURES)]
assert len(keep)==len(genes)-2
C=C[:,keep]; genes=[genes[i] for i in keep]          # antigens gone BEFORE anything else
obs=A.obs.copy(); del A
is_dn_all=(obs.observed_state=='double_negative').values
depth_all=obs.depth_ex_antigen.values.astype(float)
patients=sorted(obs.patient_id.unique())
print(f"cells {C.shape[0]}, features {C.shape[1]} (antigens dropped), patients {len(patients)}")

def local(idx, is_dn, depth):
    """Fresh un-integrated local embedding. No Harmony, no global PCA, no stage-05 Leiden."""
    ad=sc.AnnData(C[idx].copy())
    ad.var_names=genes
    sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad)
    nhvg=min(2000, ad.n_vars-1)
    sc.pp.highly_variable_genes(ad,n_top_genes=nhvg); ad=ad[:,ad.var.highly_variable].copy()
    npc=int(min(30, ad.n_obs-1, ad.n_vars-1))
    sc.pp.scale(ad,max_value=10); sc.tl.pca(ad,n_comps=npc,svd_solver='arpack',random_state=SEED)
    k=int(min(15, ad.n_obs-1))
    sc.pp.neighbors(ad,n_neighbors=k,use_rep='X_pca',random_state=SEED)
    sc.tl.leiden(ad,resolution=1.0,key_added='local',flavor='igraph',n_iterations=2,
                 directed=False,random_state=SEED)
    g=ad.obsp['connectivities']
    nn=np.asarray([np.argsort(-np.asarray(g[i].todense()).ravel())[:k] for i in range(ad.n_obs)])
    stat=lambda l: SC.knn_dn_fraction(l,nn)
    bins=SC.adaptive_depth_bins(depth)
    pure=all(len(np.unique(is_dn[bins==b]))==1 for b in np.unique(bins))
    _,p_strat=SC.depth_stratified_permutation(is_dn,bins,stat,NPERM,SEED)
    _,p_unc=SC.depth_stratified_permutation(is_dn,np.zeros(len(is_dn),int),stat,NPERM,SEED)
    return dict(knn_dn_frac=stat(is_dn), dn_rate=float(is_dn.mean()),
        morans_i=SC.morans_i(is_dn,g), best_cluster_enr=SC.best_cluster_enrichment(is_dn,ad.obs['local'].values),
        n_local_clusters=int(ad.obs['local'].nunique()), n_depth_bins=int(len(np.unique(bins))),
        strata_pure=bool(pure), perm_p_depth_stratified=p_strat, perm_p_unconditioned=p_unc,
        n_pcs=npc, n_hvg=nhvg), ad, nn, bins

def de(idx, is_dn, bins):
    d,p=SC.depth_matched_indices(is_dn,bins,SEED)
    if d.size<SC.MIN_GROUP_CELLS or p.size<SC.MIN_GROUP_CELLS: return None,None
    sel=np.concatenate([d,p])
    ad=sc.AnnData(C[idx[sel]].copy()); ad.var_names=genes
    ad.obs['grp']=pd.Categorical(['DN']*d.size+['AgPos']*p.size,categories=['AgPos','DN'])
    sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad)
    sc.pp.filter_genes(ad,min_cells=3)
    sc.tl.rank_genes_groups(ad,'grp',groups=['DN'],reference='AgPos',method='wilcoxon')
    r=sc.get.rank_genes_groups_df(ad,'DN')
    return r, int(d.size)

rows=[]; de_store={}; struct=[]
for pat in patients:
    pm=(obs.patient_id==pat).values
    for tag,dm in [('primary',obs.in_primary.values.astype(bool)),
                   ('sensitivity',obs.in_sensitivity.values.astype(bool))]:
        m=pm&dm; idx=np.flatnonzero(m)
        dn=is_dn_all[idx]; dep=depth_all[idx]
        n,ndn,npos=len(idx),int(dn.sum()),int((~dn).sum())
        base={'patient':pat,'denominator':tag,'n_cells':n,'n_dn':ndn,'n_antigen_pos':npos}
        evaluable=(ndn>=SC.MIN_GROUP_CELLS and npos>=SC.MIN_GROUP_CELLS
                   and int((pm&obs.in_primary.values.astype(bool)).sum())>=SC.MIN_PATIENT_CELLS)
        base['evaluable']=bool(evaluable)
        if not evaluable:
            base['reason']=('dn<20' if ndn<SC.MIN_GROUP_CELLS else
                            'antigen_pos<20' if npos<SC.MIN_GROUP_CELLS else 'n_primary<100')
            rows.append(base); continue
        base['reason']=''
        st,ad,nn,bins=local(idx,dn,dep)
        base.update(st)
        # depth diagnostic
        for nm,v in [('depth_ex_antigen',dep),('n_genes',np.asarray((C[idx]>0).sum(1)).ravel()),
                     ('total_umi',np.asarray(C[idx].sum(1)).ravel())]:
            base[f'med_{nm}_dn']=float(np.median(v[dn])); base[f'med_{nm}_pos']=float(np.median(v[~dn]))
            base[f'ratio_{nm}']=float(np.median(v[dn])/max(np.median(v[~dn]),1e-9))
        base['n_samples']=int(obs.sample_name.values[idx].__len__() and len(np.unique(obs.sample_name.values[idx])))
        r,nmatch=de(idx,dn,bins)
        if r is not None:
            base['n_de_matched_cells_per_group']=nmatch
            sig=r[(r.pvals_adj<0.05)]
            base['n_de_padj05']=int(len(sig))
            base['n_de_up']=int((sig.logfoldchanges>0).sum()); base['n_de_down']=int((sig.logfoldchanges<0).sum())
            de_store[(pat,tag)]=r.set_index('names')['logfoldchanges']
            r.head(400).assign(patient=pat,denominator=tag).to_csv(
                OUT/f'_de_{pat}_{tag}.csv',index=False)
        rows.append(base)
        struct.append({'patient':pat,'denominator':tag,**st})
    print('.',end='',flush=True)
E=pd.DataFrame(rows); E.to_csv(OUT/'stage10_evaluability.csv',index=False)
pd.DataFrame(struct).to_csv(OUT/'dn_local_structure_by_patient.csv',index=False)
E[[c for c in E.columns if c.startswith(('patient','denominator','med_','ratio_','n_dn','n_antigen_pos','n_samples'))]]\
 .to_csv(OUT/'dn_vs_nondn_depth_diagnostic.csv',index=False)
pickle.dump((E,de_store),open('/tmp/s10.pkl','wb'))
print(f"\nevaluable rows {int(E.evaluable.sum())} of {len(E)}")
print(E[E.evaluable][['patient','denominator','n_cells','n_dn','n_antigen_pos','knn_dn_frac','dn_rate',
    'morans_i','perm_p_depth_stratified','perm_p_unconditioned','n_depth_bins','strata_pure',
    'ratio_depth_ex_antigen','n_de_padj05']].round(4).to_string(index=False))
