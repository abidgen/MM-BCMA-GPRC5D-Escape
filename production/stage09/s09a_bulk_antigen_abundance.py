# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 09
# Step         : s09a_bulk_antigen_abundance.py
# What it does : GSE223061 TPM load, stub exclusion, MMRF multi-run averaging
# Writes       : results/09_bulk_validation/bulk_antigen_abundance.csv
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : heredoc-written script, run with `conda run -n mm-core python <path>`
#   Original path  : $S/s09a.py  (session scratchpad, since deleted)
#   Executed (UTC) : 2026-08-25T16:36:54Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
import numpy as np, pandas as pd, tarfile, gzip, re, pickle
from pathlib import Path
OUT=Path('results/09_bulk_validation'); RAW=Path('raw/unpacked_bulk')
pd.set_option('display.width',250)
AG={'TNFRSF17':'ENSG00000048462','GPRC5D':'ENSG00000111291'}

t2g={}
with gzip.open('raw/gtf/Homo_sapiens.GRCh38.93.gtf.gz','rt') as fh:
    for line in fh:
        if line[0]=='#': continue
        f=line.split('\t',9)
        if f[2]!='transcript': continue
        gn=re.search(r'gene_name "([^"]+)"',f[8])
        if gn and gn.group(1) in AG:
            t2g[re.search(r'transcript_id "([^"]+)"',f[8]).group(1).split('.')[0]]=gn.group(1)

rows=[]
for f in sorted(RAW.glob('GSM*')):
    gsm, rest = f.name.split('_',1)
    if f.name.endswith('.tar.gz'):
        sid=rest[:-7]
        with tarfile.open(f,'r:gz') as tf:
            m=[x for x in tf.getnames() if x.endswith('abundance.tsv')]
            d=pd.read_csv(tf.extractfile(m[0]),sep='\t')
        d['gene']=d.target_id.str.split('.').str[0].map(t2g)
        agg=d.dropna(subset=['gene']).groupby('gene').tpm.sum()
        rows.append({'gsm':gsm,'bulk_sample':sid,'bulk_cohort':'WU1','specimen':'unsorted BMMC',
            'assay':'kallisto transcript TPM -> gene sum','n_features':len(d),'n_runs':1,
            'usable':True,'TNFRSF17_tpm':float(agg.get('TNFRSF17',np.nan)),
            'GPRC5D_tpm':float(agg.get('GPRC5D',np.nan)),'total_tpm':float(d.tpm.sum()),'notes':''})
    else:
        sid=rest[:-len('_tpm.tsv.gz')]
        try: d=pd.read_csv(f,sep='\t')
        except Exception: d=pd.DataFrame()
        if len(d)==0:
            rows.append({'gsm':gsm,'bulk_sample':sid,'bulk_cohort':'MMRF','specimen':'CD138+ sorted',
                'assay':'gene TPM','n_features':0,'n_runs':0,'usable':False,'TNFRSF17_tpm':np.nan,
                'GPRC5D_tpm':np.nan,'total_tpm':np.nan,'notes':'114-byte header-only stub'}); continue
        runs=sorted(d.Sample_ID.unique())
        # TPM is already per-run normalized -> average across runs, never sum
        per=d[d.Gene_Name.isin(AG)].groupby(['Sample_ID','Gene_Name']).TPM.sum().unstack()
        rows.append({'gsm':gsm,'bulk_sample':sid,'bulk_cohort':'MMRF','specimen':'CD138+ sorted',
            'assay':'gene TPM','n_features':int(len(d)/len(runs)),'n_runs':len(runs),'usable':True,
            'TNFRSF17_tpm':float(per['TNFRSF17'].mean()),'GPRC5D_tpm':float(per['GPRC5D'].mean()),
            'total_tpm':float(d.TPM.sum()/len(runs)),
            'notes':f'{len(runs)} runs averaged: {",".join(runs)}' if len(runs)>1 else ''})
B=pd.DataFrame(rows)
B.to_csv(OUT/'bulk_antigen_abundance.csv',index=False)
print(f"bulk GSM files {len(B)}; usable {int(B.usable.sum())}; "
      f"stubs {sorted(B[~B.usable].bulk_sample)}")
print(B[['bulk_sample','bulk_cohort','n_features','n_runs','TNFRSF17_tpm','GPRC5D_tpm','notes']]
      .round(2).to_string(index=False))
pickle.dump(B,open('/tmp/s09_bulk.pkl','wb'))
