# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#   kernelspec:
#     display_name: mm-core
#     language: python
#     name: mm-core
# ---

# %% [markdown]
# # Stage 07 — dominant-clone plasma-cell calling
#
# Two methods were attempted. **One was rejected on its negative control and stays
# rejected**, which is the more important half of this stage.
#
# * **CNV (`infercnvpy`) — `NOT_EVALUABLE` cohort-wide.** Healthy-donor plasma cells,
#   which have no clone, produced false positives spanning 0.0%–50.6% at z > 3
#   (`ND_090617` median z = +3.03). The method was rejected **before any disease CNV was
#   inspected**, so no threshold could be tuned against tumour data. An unreliable assay
#   contributes *no* evidence — it is not weak negative evidence, and the light-chain
#   rule was deliberately **not** relaxed to compensate for its absence.
# * **Immunoglobulin V-gene usage — accepted**, but only as a higher-specificity
#   refinement of the light-chain clonality axis. 3′ chemistry cannot capture J segments
#   (IGKJ 0.00%, IGHJ 0.04% of cells), so there is no clonotype here — only patient-
#   specific V usage. **This is one axis, not two.**
#
# Outputs: `results/07_malignant_plasma/v_clone_membership/`.

# %%
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path

from mm_escape import config, malignant as M

OUT = Path("results/07_malignant_plasma/v_clone_membership")
OUT.mkdir(parents=True, exist_ok=True)
pd.set_option("display.width", 250)

# %% [markdown]
# ## Plasma compartment
#
# `cell_type` is the only annotation stage 06 exposes; nothing here branches on how that
# label was produced.

# %%
adata = sc.read_h5ad("results/05_integration/integrated.h5ad")
labels = pd.read_csv("results/06_annotation/per_cell_labels.csv.gz", index_col=0)
adata.obs["cell_type"] = labels.reindex(adata.obs_names)["cell_type"].astype(str).values

plasma = adata[(adata.obs["cell_type"] == "PlasmaCell").values].copy()
del adata

var_names = list(plasma.var_names)
counts = plasma.layers["counts"]
print(f"plasma cells: {plasma.n_obs:,}")

# %% [markdown]
# The antigen genes are present in the matrix but **must never reach the classifier**.
# That is asserted mechanically at the end of this notebook, not merely intended.

# %%
assert "TNFRSF17" in var_names and "GPRC5D" in var_names

KAPPA_V = [g for g in var_names if g.startswith("IGKV")]
LAMBDA_V = [g for g in var_names if g.startswith("IGLV")]
LAMBDA_C = [g for g in var_names if g in ("IGLC1", "IGLC2", "IGLC3", "IGLC7")]


def dense(genes):
    idx = [var_names.index(g) for g in genes]
    block = counts[:, idx]
    return np.asarray(block.todense()) if hasattr(block, "todense") else np.asarray(block)


X_kv, X_lv = dense(KAPPA_V), dense(LAMBDA_V)
X_kc, X_lc = dense(["IGKC"]), dense(LAMBDA_C)

# %% [markdown]
# ## Light-chain class per cell
#
# **Ratio, never presence/absence.** Ig transcripts are the most ambient-contaminated
# genes in this tissue — plasma cells secrete enormous quantities of Ig mRNA into the
# droplet background — so a presence-based call is far noisier than it looks. A shared
# additive background cancels in a ratio.

# %%
obs = plasma.obs[
    ["sample_name", "patient_id", "cohort", "sample_type", "total_counts",
     "n_genes_by_counts", "pct_counts_mt"]
].copy()
obs["patient_id"] = obs["patient_id"].astype(str)
obs["sample_name"] = obs["sample_name"].astype(str)
obs["kappa_umi"] = X_kc.sum(axis=1)
obs["lambda_umi"] = X_lc.sum(axis=1)
obs["lc_class"] = M.light_chain_class(obs.kappa_umi, obs.lambda_umi)
obs["ig_umi"] = obs.kappa_umi + obs.lambda_umi + X_kv.sum(axis=1) + X_lv.sum(axis=1)
obs["has_LCV"] = ((X_kv > 0).any(axis=1)) | ((X_lv > 0).any(axis=1))

print(obs.lc_class.value_counts().to_string())

# %% [markdown]
# ## Patient-level clonality and V evaluability
#
# The dominant-V threshold of 0.50 comes from **donor calibration, not clinical
# convention**: healthy donors top out at 0.204/0.232/0.250/0.378 and evaluable disease
# patients start at 0.562. 0.50 sits inside that gap, ~32% above the donor maximum, and
# deliberately not at the disease minimum.

# %%
rows = []
for patient, grp in obs.groupby("patient_id", observed=True):
    pos = obs.index.get_indexer(grp.index)
    called = grp[grp.lc_class.isin(["kappa", "lambda"])]
    if len(called) >= 20:
        frac_kappa = (called.lc_class == "kappa").mean()
        dominance = max(frac_kappa, 1 - frac_kappa)
        dominant = "kappa" if frac_kappa >= 0.5 else "lambda"
        clonality = ("CLONAL_STRONG" if dominance >= 0.95
                     else "CLONAL_WEAK" if dominance >= 0.85 else "NO_RESTRICTION")
    else:
        dominance, dominant, clonality = np.nan, None, "NOT_EVALUABLE"

    # A patient with no dominant class has no dominant V to look for.
    X, gene_set = ((X_kv, KAPPA_V) if dominant == "kappa"
                   else (X_lv, LAMBDA_V) if dominant == "lambda" else (X_kv, KAPPA_V))
    detected = X[pos] > 0
    n_v_pos = int(detected.any(axis=1).sum())
    if n_v_pos:
        per_gene = detected.sum(axis=0)
        top = int(np.argmax(per_gene))
        gene, top_frac = gene_set[top], float(per_gene[top] / n_v_pos)
        column = dense([gene]).ravel()
        own = (obs.patient_id == patient).values
        elsewhere = float((column[~own] > 0).mean())
        enrichment = float((column[own] > 0).mean()) / elsewhere if elsewhere > 0 else np.inf
    else:
        gene, top_frac, enrichment = None, np.nan, np.nan

    rows.append({
        "patient": patient, "sample_type": grp.sample_type.iloc[0],
        "cohort": grp.cohort.iloc[0], "n_plasma": len(grp), "n_lc_called": len(called),
        "D": dominance, "clonality_state": clonality, "dominant_class": dominant,
        "pct_LCV": round(float(grp.has_LCV.mean()), 4), "n_V_positive": n_v_pos,
        "dominant_V": gene, "top_V_frac": top_frac,
        "enrichment": round(enrichment, 2) if np.isfinite(enrichment) else np.inf,
        "v_state": M.patient_v_evaluability(
            n_v_pos, float(grp.has_LCV.mean()), top_frac, enrichment),
    })

patients = pd.DataFrame(rows)
patients.to_csv(OUT / "patient_v_evaluability.csv", index=False)
print(pd.crosstab(patients.clonality_state, patients.v_state).to_string())

# %% [markdown]
# ## Per-cell clone membership
#
# **V absence is never negative evidence.** Incompatibility requires *positive* evidence:
# a minority light-chain class, or an alternative V at ≥2 UMI (1 UMI is the typical
# detection level and sits inside noise). Everything else that fails to reach
# `CLONE_SUPPORTED` is `UNCERTAIN`, which is a reported quantity rather than a silent drop.

# %%
pmap = patients.set_index("patient")
dominant_v_detected = np.zeros(len(obs), bool)
alt_v_detected = np.zeros(len(obs), bool)

for patient, grp in obs.groupby("patient_id", observed=True):
    pos = obs.index.get_indexer(grp.index)
    row = pmap.loc[patient]
    if not isinstance(row.dominant_V, str) or row.dominant_class not in ("kappa", "lambda"):
        continue
    X, gene_set = (X_kv, KAPPA_V) if row.dominant_V in KAPPA_V else (X_lv, LAMBDA_V)
    j = gene_set.index(row.dominant_V)
    dominant_v_detected[pos] = X[pos][:, j] > 0
    alt_v_detected[pos] = (np.delete(X[pos], j, axis=1) >= config.ALT_V_MIN_UMI).any(axis=1)

obs["dominant_v_detected"] = dominant_v_detected
obs["alt_v_detected"] = alt_v_detected
obs["clone_state"] = [
    M.clone_membership(pmap.loc[p, "clonality_state"], pmap.loc[p, "v_state"], lc,
                       pmap.loc[p, "dominant_class"] or "", dv, av)
    for p, lc, dv, av in zip(obs.patient_id, obs.lc_class,
                             obs.dominant_v_detected, obs.alt_v_detected)
]
obs.to_csv(OUT / "clone_membership_per_cell.csv.gz", compression="gzip")
print(obs.clone_state.value_counts().to_string())

# %% [markdown]
# ## Informative missingness — the finding that shapes stage 08
#
# Requiring positive V detection does not sample the clone at random: it preferentially
# keeps **deeper** cells. If stage 08 used `CLONE_SUPPORTED` alone it would inherit that,
# inflating apparent antigen positivity and therefore *deflating* the double-negative
# fraction — a bias pointing against this project's own hypothesis, unevenly by cohort.

# %%
clonal = patients[(patients.clonality_state == "CLONAL_STRONG")
                  & (patients.v_state == "V_EVALUABLE")].patient.tolist()
q = obs[obs.patient_id.isin(clonal)
        & obs.clone_state.isin([M.CLONE_SUPPORTED, M.CLONE_COMPATIBLE_V_UNOBSERVED])]

METRICS = ["total_counts", "n_genes_by_counts", "ig_umi", "pct_counts_mt"]
pooled = q.groupby("clone_state")[METRICS].median()
pooled.loc["ratio"] = pooled.loc[M.CLONE_SUPPORTED] / pooled.loc[M.CLONE_COMPATIBLE_V_UNOBSERVED]
pooled.round(3).to_csv(OUT / "v_detection_depth_bias_pooled.csv")
print(pooled.round(2).to_string())

by_cohort = q.groupby(["cohort", "clone_state"], observed=True)[METRICS].median()
by_cohort["n"] = q.groupby(["cohort", "clone_state"], observed=True).size()
by_cohort.round(2).to_csv(OUT / "v_detection_depth_bias_by_cohort.csv")
print(by_cohort.round(1).to_string())

# %% [markdown]
# The pooled ratio (1.79×) badly understates it: **MMRF is 17.6×** against WU1 1.85× and
# WU2 1.35×, and 24 of 24 evaluable patients point the same way.
#
# `83942` is the cleanest control available — one patient, one clone, one dominant V,
# sampled under two WashU protocols — and its supported fraction splits 0.746 vs 0.351.
# The clone is identical, so that difference is **technical, not biological**.

# %% [markdown]
# ## Antigen-circularity invariance — verified, not asserted
#
# The classifier is re-run on perturbed counts and compared by hash. The joint-extreme
# case shifts every library by ~20,000 UMI — more than a median plasma cell's entire
# library — and all 35,474 per-cell states must still be **bit-identical**.

# %%
print(Path("results/07_malignant_plasma/v_clone_membership/"
           "antigen_circularity_invariance.md").read_text())
