# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 07
# Step         : s07g_candidate_clone_signatures.sh
# What it does : candidate per-patient clone signatures
# Writes       : results/07_malignant_plasma/ig_clone_feasibility/candidate_clone_signatures.csv
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : inline shell/python block
#   Original path  : (inline heredoc executed directly from the shell)
#   Executed (UTC) : 2026-08-25T15:35:26Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
cd /media/wrath/CART_mm_dual_antigen
/home/abid/miniforge3/envs/mm-core/bin/python << 'EOF' 2>&1 | grep -v Warning | tail -12
import pandas as pd, numpy as np
from pathlib import Path
OUT=Path('results/07_malignant_plasma/ig_clone_feasibility')
rd=lambda f: pd.read_csv(OUT/f, dtype={'patient':str})
P,S,E,R = rd('patient_vj_summary.csv'), rd('cross_patient_specificity.csv'), rd('vj_evaluability.csv'), rd('repeated_sample_vj_consistency.csv')
c=P.merge(S[['patient','top_V_gene','frac_own_patient_cells','frac_other_patients',
             'n_other_patients_with_it','frac_donor_plasma']],on=['patient','top_V_gene'],how='left')
c=c.merge(E[['patient','state']],on='patient',how='left')
rep=R[R.n_V_pos>=20].groupby('patient')['matches_patient'].agg(['size','sum'])
c=c.merge(rep.rename(columns={'size':'n_evaluable_samples','sum':'n_samples_matching'}),
          left_on='patient',right_index=True,how='left')
c['enrichment_vs_others']=(c.frac_own_patient_cells/c.frac_other_patients.replace(0,np.nan)).round(1)
c.sort_values('n_plasma',ascending=False).to_csv(OUT/'candidate_clone_signatures.csv',index=False)
d=c[(c.sample_type=='myeloma')&(c.state=='V_EVALUABLE')]
print(f"V_EVALUABLE myeloma patients: {len(d)}")
print(f"  top_V_frac_of_Vpos: min {d.top_V_frac_of_Vpos.min():.3f} q25 {d.top_V_frac_of_Vpos.quantile(.25):.3f} median {d.top_V_frac_of_Vpos.median():.3f}")
print(f"  enrichment vs other patients: median {d.enrichment_vs_others.median():.1f}x  min {d.enrichment_vs_others.min():.1f}x")
print(f"  distinct top V genes across {len(d)} patients: {d.top_V_gene.nunique()} (most shared: {d.top_V_gene.value_counts().iloc[0]})")
print(f"  top V gene >5% detected in donor plasma: {int((d.frac_donor_plasma>0.05).sum())} patients")
print(f"  median pct of plasma cells that are LCV-positive: {d.pct_with_LCV.median():.1f}%")
EOF
ls results/07_malignant_plasma/ig_clone_feasibility/
