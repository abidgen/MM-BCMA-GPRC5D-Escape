# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 06
# Step         : s06i_trbc_context_revision.sh
# What it does : accepted TRBC-context revision of cluster 23
# Writes       : results/06_annotation/cluster23_local/trbc_context_revision/*
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : inline shell/python block
#   Original path  : (inline heredoc executed directly from the shell)
#   Executed (UTC) : 2026-08-25T11:15:17Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
cd /media/wrath/CART_mm_dual_antigen
mkdir -p results/06_annotation/cluster23_local/trbc_context_revision
/home/abid/miniforge3/envs/mm-annotation/bin/python << 'EOF' 2>&1 | grep -v Warning | tail -60
import numpy as np, pandas as pd, scanpy as sc
from pathlib import Path
from mm_escape import annotation as ann, config
OUT=Path('results/06_annotation/cluster23_local/trbc_context_revision')
pd.set_option('display.width',250)

a = sc.read_h5ad('results/06_annotation/cluster23_local/cluster23_local.h5ad')
old = a.obs['lineage_call'].astype(str)          # frozen Part-B call
r = ann.cytotoxic_lineage_calls(a)               # revised: TRBC = context
new = r['call']
o = a.obs

r.join(o[['sample_name','patient_id','cohort','local_cluster','total_counts',
          'n_genes_by_counts','doublet_score']]).assign(frozen_call=old.values)\
 .to_csv(OUT/'revised_lineage_calls.csv.gz', compression='gzip')

trans = pd.crosstab(old.values, new.values)
trans.to_csv(OUT/'old_vs_new_transition_table.csv')
print("=== OLD -> NEW transition ==="); print(trans.to_string())
print("\n=== counts ===")
print(pd.DataFrame({'frozen_PartB': old.value_counts(), 'revised': new.value_counts()}
                   ).fillna(0).astype(int).to_string())

# per-category evidence summary
T=list(config.T_IDENTITY_ANCHORS); CTX=list(config.T_CONTEXT)
X=a[:,T+CTX].layers['counts']; X=np.asarray(X.todense()) if hasattr(X,'todense') else np.asarray(X)
c=pd.DataFrame(X,columns=T+CTX,index=a.obs_names)
rows=[]
for k in new.unique():
    m=(new==k).values; s=o[m]
    rows.append({'call':k,'n_cells':int(m.sum()),'pct':round(100*m.mean(),2),
      'n_patients':s['patient_id'].nunique(),'n_samples':s['sample_name'].nunique(),
      'top_patient_frac':round(s['patient_id'].value_counts(normalize=True).iat[0],3),
      'frac_anyCD3':round(float((c.loc[m,['CD3D','CD3E','CD3G']]>0).any(axis=1).mean()),3),
      'frac_2CD3':round(float((c.loc[m,['CD3D','CD3E','CD3G']]>0).sum(axis=1).ge(2).mean()),3),
      'frac_TRAC':round(float((c.loc[m,'TRAC']>0).mean()),3),
      'frac_CD3_and_TRAC':round(float(((c.loc[m,['CD3D','CD3E','CD3G']]>0).any(axis=1)&(c.loc[m,'TRAC']>0)).mean()),3),
      'frac_TRBC1':round(float((c.loc[m,'TRBC1']>0).mean()),3),
      'frac_TRBC2':round(float((c.loc[m,'TRBC2']>0).mean()),3),
      'median_TRBC_umi_pos':float(c.loc[m,CTX].max(axis=1).replace(0,np.nan).median()),
      'mean_NK_anchors':round(float(r.loc[m,'n_NK_identity'].mean()),2),
      'mean_cytotoxic':round(float(r.loc[m,'n_cytotoxic_state'].mean()),2),
      'median_counts':float(s['total_counts'].median()),
      'median_genes':float(s['n_genes_by_counts'].median()),
      'median_doublet':round(float(s['doublet_score'].median()),4),
      'cohorts':dict(s['cohort'].value_counts())})
summ=pd.DataFrame(rows).sort_values('n_cells',ascending=False)
summ.to_csv(OUT/'revised_lineage_evidence_summary.csv',index=False)
print("\n=== revised categories ===")
print(summ[['call','n_cells','pct','n_patients','n_samples','top_patient_frac',
            'frac_anyCD3','frac_CD3_and_TRAC','frac_TRBC1','mean_NK_anchors','median_doublet']].to_string(index=False))

# the movers
mv=((old.values=='NKT_like_mixed')&(new.values=='NK'))
print(f"\n=== movers mixed -> NK: {int(mv.sum())} ===")
print(f"  any CD3 {100*(c.loc[mv,['CD3D','CD3E','CD3G']]>0).any(axis=1).mean():.1f}% | "
      f"TRAC {100*(c.loc[mv,'TRAC']>0).mean():.1f}% | "
      f"TRBC+ {100*(c.loc[mv,CTX]>0).any(axis=1).mean():.1f}% | "
      f"median TRBC UMI {c.loc[mv,CTX].max(axis=1).median():.0f} | "
      f"frac==1 UMI {100*(c.loc[mv,CTX].max(axis=1)==1).mean():.1f}%")
print(f"  mean NK anchors {r.loc[mv,'n_NK_identity'].mean():.2f} | "
      f"patients {o[mv]['patient_id'].nunique()} | samples {o[mv]['sample_name'].nunique()} | "
      f"top-patient {o[mv]['patient_id'].value_counts(normalize=True).iat[0]:.3f}")

pd.crosstab(new.values,o['patient_id']).to_csv(OUT/'revised_patient_representation.csv')
pd.crosstab(new.values,o['sample_name']).to_csv(OUT/'revised_sample_representation.csv')
ctx=pd.DataFrame({'call':new.values,'TRBC_max':c[CTX].max(axis=1).values,
                  'TRBC_pos':(c[CTX]>0).any(axis=1).values}).groupby('call').agg(
   n=('TRBC_pos','size'),frac_TRBC_pos=('TRBC_pos','mean'),
   median_TRBC_umi=('TRBC_max','median'),mean_TRBC_umi=('TRBC_max','mean'))
ctx.round(3).to_csv(OUT/'revised_trbc_context_summary.csv')
print("\n=== TRBC context preserved per revised call ==="); print(ctx.round(3).to_string())
EOF
