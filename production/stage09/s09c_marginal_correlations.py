# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 09
# Step         : s09c_marginal_correlations.py
# What it does : marginal TNFRSF17/GPRC5D pseudobulk-vs-bulk correlations, both denominators, truncate-10k
# Writes       : results/09_bulk_validation/sc_marginal_antigen_{primary,sensitivity}.csv; bulk_vs_sc_tnfrsf17.csv; bulk_vs_sc_gprc5d.csv; bulk_vs_sc_by_cohort.csv; bulk_vs_sc_truncate10k_sensitivity.csv
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : heredoc-written script, run with `conda run -n mm-core python <path>`
#   Original path  : $S/s09c.py  (session scratchpad, since deleted)
#   Executed (UTC) : 2026-08-25T16:38:07Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
import numpy as np, pandas as pd, pickle
from pathlib import Path
from scipy.stats import spearmanr
from mm_escape import antigen as A
OUT=Path('results/09_bulk_validation')
St=pickle.load(open('/tmp/s08.pkl','rb')); M=pickle.load(open('/tmp/s09_map.pkl','rb'))
obs=St['obs']; prim,sens=St['prim'],St['sens']; bc,gc,tot=St['bc'],St['gc'],St['tot']
pd.set_option('display.width',260); SEED=20260825
bt=A.downsample_gene_counts(tot,bc,10000,SEED); gt=A.downsample_gene_counts(tot,gc,10000,SEED)
tt=np.minimum(tot,10000)

mt=M[M.match_status=='MATCHED_EXACT'].copy()
# FROZEN rule, declared before any correlation: one pair per patient = the EARLIEST
# matched timepoint that has >=20 denominator cells (the pre-existing per-sample floor).
def n_prim(s): return int((prim&(obs.sample_name==s).values).sum())
mt['n_primary_cells']=mt.sc_sample.map(n_prim)
mt['assessable']=mt.n_primary_cells>=20
sel=(mt[mt.assessable].sort_values(['sc_patient','timepoint'])
       .groupby('sc_patient',as_index=False).first())
print(f"matched pairs {len(mt)}; assessable (>=20 primary cells) {int(mt.assessable.sum())}; "
      f"patient-level observations {len(sel)}")
print(mt[['sc_patient','sc_sample','bulk_sample','n_primary_cells','assessable']].to_string(index=False))

def cpm(mask, sample, gene_counts, totals):
    m=mask&(obs.sample_name==sample).values
    if m.sum()==0: return np.nan, 0
    T=totals[m].sum()
    return (1e6*gene_counts[m].sum()/T if T>0 else np.nan), int(m.sum())

rows=[]
for _,r in sel.iterrows():
    d={'sc_patient':r.sc_patient,'sc_sample':r.sc_sample,'bulk_sample':r.bulk_sample,
       'cohort':r.bulk_cohort,'specimen':r.specimen,'timepoint':r.timepoint,
       'bulk_TNFRSF17_tpm':r.TNFRSF17_tpm,'bulk_GPRC5D_tpm':r.GPRC5D_tpm}
    for tag,mask in [('primary',prim),('sensitivity',sens)]:
        for gname,gv,gt_ in [('TNFRSF17',bc,bt),('GPRC5D',gc,gt)]:
            v,n=cpm(mask,r.sc_sample,gv,tot)
            d[f'sc_{gname}_cpm_{tag}']=v
            d[f'sc_{gname}_detect_{tag}']=float((gv[mask&(obs.sample_name==r.sc_sample).values]>0).mean())
            d[f'sc_{gname}_cpm_{tag}_trunc10k']=cpm(mask,r.sc_sample,gt_,tt)[0]
        d[f'n_cells_{tag}']=n
    # whole-sample pseudobulk: the composition-appropriate pairing for unsorted WU1 bulk
    ws=(obs.sample_name==r.sc_sample).values
    d['sc_TNFRSF17_cpm_wholesample']=1e6*bc[ws].sum()/tot[ws].sum()
    d['sc_GPRC5D_cpm_wholesample']=1e6*gc[ws].sum()/tot[ws].sum()
    rows.append(d)
D=pd.DataFrame(rows)
for tag in ['primary','sensitivity']:
    cols=['sc_patient','sc_sample','bulk_sample','cohort',f'n_cells_{tag}',
          'bulk_TNFRSF17_tpm',f'sc_TNFRSF17_cpm_{tag}',f'sc_TNFRSF17_detect_{tag}',
          'bulk_GPRC5D_tpm',f'sc_GPRC5D_cpm_{tag}',f'sc_GPRC5D_detect_{tag}']
    D[cols].round(4).to_csv(OUT/f'sc_marginal_antigen_{tag}.csv',index=False)
print("\n=== MATCHED PATIENT-LEVEL DATA (primary denominator) ===")
print(D[['sc_patient','cohort','n_cells_primary','bulk_TNFRSF17_tpm','sc_TNFRSF17_cpm_primary',
         'bulk_GPRC5D_tpm','sc_GPRC5D_cpm_primary']].round(2).to_string(index=False))

def sp(x,y,label,n_min=5):
    ok=np.isfinite(x)&np.isfinite(y)
    if ok.sum()<n_min: return {'comparison':label,'n':int(ok.sum()),'rho':np.nan,'p':np.nan,
                               'status':'NOT_EVALUABLE (n<5)'}
    rho,p=spearmanr(x[ok],y[ok])
    return {'comparison':label,'n':int(ok.sum()),'rho':round(float(rho),4),
            'p':round(float(p),4),'status':'evaluable'}

res=[]
for gname in ['TNFRSF17','GPRC5D']:
    for tag in ['primary','sensitivity']:
        for sub,lab in [(D,'pooled'),(D[D.cohort=='MMRF'],'MMRF'),(D[D.cohort=='WU1'],'WU1')]:
            res.append({**sp(sub[f'bulk_{gname}_tpm'].values,sub[f'sc_{gname}_cpm_{tag}'].values,
                            f'{gname} | {tag} | {lab}'),'antigen':gname,'denominator':tag,'stratum':lab})
    for sub,lab in [(D,'pooled'),(D[D.cohort=='WU1'],'WU1')]:
        res.append({**sp(sub[f'bulk_{gname}_tpm'].values,sub[f'sc_{gname}_cpm_wholesample'].values,
                        f'{gname} | whole-sample pseudobulk | {lab}'),
                    'antigen':gname,'denominator':'whole_sample','stratum':lab})
R=pd.DataFrame(res)
R[R.antigen=='TNFRSF17'].to_csv(OUT/'bulk_vs_sc_tnfrsf17.csv',index=False)
R[R.antigen=='GPRC5D'].to_csv(OUT/'bulk_vs_sc_gprc5d.csv',index=False)
R.to_csv(OUT/'bulk_vs_sc_by_cohort.csv',index=False)
print("\n=== BULK vs SINGLE-CELL, Spearman (patient-level) ===")
print(R[['comparison','n','rho','p','status']].to_string(index=False))

tr=[]
for gname in ['TNFRSF17','GPRC5D']:
    for tag in ['primary','sensitivity']:
        a=sp(D[f'bulk_{gname}_tpm'].values,D[f'sc_{gname}_cpm_{tag}'].values,'orig')
        b=sp(D[f'bulk_{gname}_tpm'].values,D[f'sc_{gname}_cpm_{tag}_trunc10k'].values,'trunc')
        tr.append({'antigen':gname,'denominator':tag,'n':a['n'],'rho_original':a['rho'],
                   'rho_truncate10k':b['rho'],'delta':round((b['rho'] or 0)-(a['rho'] or 0),4)})
T=pd.DataFrame(tr); T.to_csv(OUT/'bulk_vs_sc_truncate10k_sensitivity.csv',index=False)
print("\n=== 10. TRUNCATE-10k SENSITIVITY (frozen Stage-08 procedure, unchanged) ===")
print(T.to_string(index=False))
pickle.dump((D,R,mt),open('/tmp/s09_res.pkl','wb'))
