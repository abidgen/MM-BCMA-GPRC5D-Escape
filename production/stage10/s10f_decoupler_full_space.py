# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 10
# Step         : s10f_decoupler_full_space.py
# What it does : decoupler re-run retaining the full tested space
# Writes       : results/10_dn_coherence/decoupler_pathway_results.csv
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : heredoc-written script, run with `conda run -n mm-core python <path>`
#   Original path  : $S/s10f.py  (session scratchpad, since deleted)
#   Executed (UTC) : 2026-08-25T22:59:14Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
import numpy as np, pandas as pd, pickle, warnings
from pathlib import Path
import decoupler as dc
from mm_escape import communication as CM
warnings.filterwarnings('ignore')
OUT=Path('results/10_dn_coherence'); pd.set_option('display.width',260)
DE,TC,cm=pickle.load(open('/tmp/s10e.pkl','rb'))
print("=== cross-program correlation of per-patient matched effects (primary) ===")
print(cm.round(3).to_string())

rows=[]
for tag in ['primary','sensitivity']:
    d=DE[DE.denominator==tag].copy()
    stat=d.set_index('gene')[['median_log2FC']].T
    stat.index=['DN_vs_AgPos']
    for nm,fn in [('hallmark',dc.op.hallmark),('progeny',dc.op.progeny),('collectri',dc.op.collectri)]:
        try:
            net=fn(organism='human')
            src,tgt=('source','target')
            res=dc.mt.ulm(data=stat,net=net,tmin=5,verbose=False)
            sc_,pv=res[0] if isinstance(res,tuple) else (stat.mul(0),None)
            act=stat.copy()
        except Exception as e:
            rows.append({'denominator':tag,'resource':nm,'status':f'FAILED {type(e).__name__}: {e}'}); continue
        try:
            est=dc.pp.get_obsm if False else None
        except Exception: pass
        rows.append({'denominator':tag,'resource':nm,'status':'ran','n_sources':int(net['source'].nunique())})
print()
# proper 2.x call returns an AnnData-like; do it explicitly
out=[]
for tag in ['primary','sensitivity']:
    d=DE[DE.denominator==tag]
    mat=pd.DataFrame([d.set_index('gene').median_log2FC],index=['DN_vs_AgPos'])
    for nm,fn in [('hallmark',dc.op.hallmark),('progeny',dc.op.progeny),('collectri',dc.op.collectri)]:
        net=fn(organism='human')
        try:
            r=dc.mt.ulm(data=mat,net=net,tmin=5)
            score=r[0] if isinstance(r,tuple) else mat
            pval=r[1] if isinstance(r,tuple) else None
        except Exception as e:
            print(f"{tag} {nm}: FAILED {type(e).__name__}: {e}"); continue
        s=score.iloc[0]; p=pval.iloc[0] if pval is not None else pd.Series(np.nan,index=s.index)
        t=pd.DataFrame({'denominator':tag,'resource':nm,'source':s.index,'score':s.values,'p':p.values})
        t['p_BH']=CM.benjamini_hochberg(t.p.values)
        out.append(t)
        print(f"{tag:12s} {nm:10s} sources tested {len(t):5d}  BH<0.05: {int((t.p_BH<0.05).sum())}")
P=pd.concat(out) if out else pd.DataFrame()
P.round(5).to_csv(OUT/'decoupler_pathway_results.csv',index=False)
if len(P):
    sig=P[(P.p_BH<0.05)]
    both=set(sig[(sig.denominator=='primary')].source)&set(sig[sig.denominator=='sensitivity'].source)
    print(f"\n  significant under BOTH denominators: {len(both)}")
    show=P[(P.denominator=='primary')&(P.source.isin(both))].nsmallest(15,'p')
    print(show[['resource','source','score','p','p_BH']].round(4).to_string(index=False))
