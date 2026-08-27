# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 08
# Step         : s08d_patient_evidence_states.py
# What it does : per-patient uncertainty/evidence flags
# Writes       : results/08_dual_antigen_escape/patient_evidence_states.csv
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : heredoc-written script, run with `conda run -n mm-core python <path>`
#   Original path  : $S/s08d.py  (session scratchpad, since deleted)
#   Executed (UTC) : 2026-08-25T16:09:34Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
import numpy as np, pandas as pd
from pathlib import Path
OUT=Path('results/08_dual_antigen_escape'); pd.set_option('display.width',250)
P=pd.read_csv(OUT/'patient_antigen_states_primary.csv')
N=pd.read_csv(OUT/'depth_stratified_null.csv')
T=pd.read_csv(OUT/'truncate10k_sensitivity.csv')
V=pd.read_csv(OUT/'primary_vs_sensitivity_denominator_comparison.csv')
d=P.merge(T[['patient','dn_delta','gprc5d_delta']],on='patient').merge(
    V[['patient','dn_delta','expansion']],on='patient',suffixes=('_trunc','_denom')).merge(
    N[['patient','enr_strat_globalbins','enr_uncond']],on='patient')
med=d.observed_double_negative_fraction.median()
def cat(r):
    fl=[]
    if r.n_cells<100: fl.append('low_n')
    if abs(r.dn_delta_denom)>0.05: fl.append('denominator_unstable')
    if abs(r.dn_delta_trunc)>0.05: fl.append('depth_sensitive')
    if (r.enr_strat_globalbins-1)*(r.enrichment_stratified-1)<0: fl.append('null_scheme_sensitive')
    enr = r.perm_p<0.05 and r.enrichment_stratified>1
    base = ('high_observed_DN' if r.observed_double_negative_fraction>=med else 'low_observed_DN')
    base += '/enriched' if enr else '/not_enriched'
    return base, ';'.join(fl) if fl else 'none'
d[['evidence_state','uncertainty_flags']]=d.apply(lambda r: pd.Series(cat(r)),axis=1)
d['stable_across_all']=(d.uncertainty_flags=='none')
keep=['patient','cohort','n_cells','n_samples','observed_double_negative_fraction','dn_ci_lo','dn_ci_hi',
 'bcma_detect','gprc5d_detect','expected_dn_stratified','enrichment_stratified','enr_ci_lo','enr_ci_hi',
 'perm_p','enr_uncond','dn_delta_trunc','dn_delta_denom','evidence_state','uncertainty_flags','stable_across_all']
d[keep].round(4).to_csv(OUT/'patient_evidence_states.csv',index=False)
print("=== 15. PATIENT EVIDENCE STATES (descriptive; NO risk tiers assigned) ===")
print(d[['patient','cohort','n_cells','observed_double_negative_fraction','enrichment_stratified',
         'perm_p','evidence_state','uncertainty_flags']].round(4).to_string(index=False))
print("\n  counts:"); print(d.evidence_state.value_counts().to_string())
print(f"\n  stable across every analysis (no flags): {int(d.stable_across_all.sum())} of {len(d)}")
print("  flag frequencies:")
for f in ['low_n','denominator_unstable','depth_sensitive','null_scheme_sensitive']:
    print(f"    {f}: {int(d.uncertainty_flags.str.contains(f).sum())}")
print("\n=== COHORT SUMMARY (primary) ===")
cs=d.groupby('cohort').agg(n_patients=('patient','size'),cells=('n_cells','sum'),
    bcma=('bcma_detect','median'),gprc5d=('gprc5d_detect','median'),
    dn=('observed_double_negative_fraction','median'),enr=('enrichment_stratified','median'),
    enr_uncond=('enr_uncond','median'),med_depth=('median_depth_ex_antigen','median'))
print(cs.round(4).to_string())
print(f"\ncohort-wide median DN {med:.4f}; range {d.observed_double_negative_fraction.min():.4f}-{d.observed_double_negative_fraction.max():.4f}")
print(f"patients with stratified enrichment significant (p<0.05 & >1): "
      f"{int(((d.perm_p<0.05)&(d.enrichment_stratified>1)).sum())} of {len(d)} -> "
      f"{sorted(d[(d.perm_p<0.05)&(d.enrichment_stratified>1)].patient)}")
print(f"patients with UNCONDITIONED enrichment > 1.2: {int((d.enr_uncond>1.2).sum())} "
      f"-> {sorted(d[d.enr_uncond>1.2].patient)}")
print(f"\nmarginal negativity, pooled over primary cells: see per-patient; median BCMA- "
      f"{(1-d.bcma_detect).median():.3f}, median GPRC5D- {(1-d.gprc5d_detect).median():.3f}")
