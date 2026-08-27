# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 07
# Step         : s07a_cnv_reference_availability_and_gate.sh
# What it does : CNV reference availability per sample and the per-patient CNV gate
# Writes       : results/07_malignant_plasma/cnv_gate_by_patient.csv; cnv_reference_availability_by_sample.csv; cnv_negative_result/donor_cnv_calibration_feasibility.csv
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : inline shell/python block
#   Original path  : (inline heredoc executed directly from the shell)
#   Executed (UTC) : 2026-08-25T11:29:44Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
cd /media/wrath/CART_mm_dual_antigen
/home/abid/miniforge3/envs/mm-core/bin/python << 'EOF' 2>&1 | grep -v Warning
import anndata as ad, pandas as pd, numpy as np
pd.set_option('display.width',220)
a = ad.read_h5ad('results/06_annotation/annotated.h5ad', backed='r')
o = a.obs[['sample_name','patient_id','cohort','sample_type','cell_type',
           'total_counts','n_genes_by_counts']].copy()
o['cell_type']=o['cell_type'].astype(str)

ref = o[o.cell_type.isin(['Tcell','Myeloid'])]
pl  = o[o.cell_type=='PlasmaCell']
t = pd.DataFrame({
 'n_plasma': pl.groupby('sample_name').size(),
 'n_T':      ref[ref.cell_type=='Tcell'].groupby('sample_name').size(),
 'n_Myeloid':ref[ref.cell_type=='Myeloid'].groupby('sample_name').size()}).fillna(0).astype(int)
t['n_ref'] = t.n_T + t.n_Myeloid
t = t[t.n_plasma > 0]
t['ref_ge200'] = t.n_ref >= 200
t['lineage_skew'] = (t[['n_T','n_Myeloid']].max(axis=1) / t.n_ref.replace(0,np.nan)).round(3)
t['grossly_dominated'] = t.lineage_skew > 0.90
meta = o.groupby('sample_name')[['cohort','sample_type']].first()
t = t.join(meta)
t.to_csv('results/07_malignant_plasma/cnv_reference_availability_by_sample.csv')

print("=== reference-cell availability, per sample with plasma cells ===")
print(f"samples: {len(t)}")
print(f"  ref >=200 PASS : {int(t.ref_ge200.sum())}  FAIL: {int((~t.ref_ge200).sum())}")
print(f"  grossly lineage-dominated (>90% one lineage): {int(t.grossly_dominated.sum())}")
print("\n  FAILING samples:")
f = t[~t.ref_ge200]
print(f[['n_plasma','n_T','n_Myeloid','n_ref','cohort','sample_type']].to_string() if len(f) else "    none")
print("\n  dominated samples:")
d = t[t.grossly_dominated]
print(d[['n_plasma','n_T','n_Myeloid','n_ref','lineage_skew','cohort']].to_string() if len(d) else "    none")

print("\n=== per-PATIENT CNV gates (G2 >=50 plasma, G3 >=200 same-sample ref, G4 median genes >=500) ===")
pat = pd.DataFrame({'n_plasma': pl.groupby('patient_id').size(),
                    'median_genes': pl.groupby('patient_id')['n_genes_by_counts'].median(),
                    'n_samples': pl.groupby('patient_id')['sample_name'].nunique(),
                    'cohort': pl.groupby('patient_id')['cohort'].first(),
                    'sample_type': pl.groupby('patient_id')['sample_type'].first()})
smap = pl.groupby('patient_id')['sample_name'].unique()
pat['n_ref_total'] = [int(t.reindex(s)['n_ref'].fillna(0).sum()) for s in smap]
pat['all_samples_ref_ge200'] = [bool(t.reindex(s)['ref_ge200'].fillna(False).all()) for s in smap]
pat['G2_plasma_ge50'] = pat.n_plasma >= 50
pat['G4_genes_ge500'] = pat.median_genes >= 500
pat['gates_pass'] = pat.G2_plasma_ge50 & pat.all_samples_ref_ge200 & pat.G4_genes_ge500
pat.sort_values('n_plasma', ascending=False).to_csv('results/07_malignant_plasma/cnv_gate_by_patient.csv')
print(f"  patients: {len(pat)}   gates PASS: {int(pat.gates_pass.sum())}   FAIL: {int((~pat.gates_pass).sum())}")
print(f"    G2 fail (<50 plasma): {int((~pat.G2_plasma_ge50).sum())}")
print(f"    G3 fail (a sample <200 ref): {int((~pat.all_samples_ref_ge200).sum())}")
print(f"    G4 fail (median genes <500): {int((~pat.G4_genes_ge500).sum())}")
print("\n  by cohort:")
print(pat.groupby('cohort')['gates_pass'].agg(['size','sum']).to_string())

print("\n=== NORMAL-DONOR calibration feasibility (8 donors) ===")
don = pl[pl.sample_type=='normal_bm']
dd = pd.DataFrame({'n_plasma': don.groupby('sample_name').size(),
                   'median_genes': don.groupby('sample_name')['n_genes_by_counts'].median(),
                   'median_counts': don.groupby('sample_name')['total_counts'].median()})
dd['n_ref'] = t.reindex(dd.index)['n_ref']
dd['G2_ge50'] = dd.n_plasma>=50; dd['G3_ge200ref'] = dd.n_ref>=200; dd['G4_ge500'] = dd.median_genes>=500
dd['gates_pass'] = dd.G2_ge50 & dd.G3_ge200ref & dd.G4_ge500
dd.to_csv('results/07_malignant_plasma/donor_cnv_calibration_feasibility.csv')
print(dd.to_string())
print(f"\n  donors passing all gates: {int(dd.gates_pass.sum())} of {len(dd)}")
EOF
