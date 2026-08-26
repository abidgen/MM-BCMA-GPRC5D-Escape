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
# # Stage 08 — BCMA/GPRC5D observed dual-negativity
#
# Every parameter below is frozen in
# `results/08_dual_antigen_escape/stage08_predeclaration.md`, written before any
# double-negative fraction existed.
#
# The distinction this notebook exists to hold: **an observed double-negative cell is a
# measurement.** Escape risk requires showing the joint zero state survives depth,
# denominator definition, repeated-sample structure and expected marginal dropout. Most
# of it does not.

# %%
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from pathlib import Path

from mm_escape import antigen as A, malignant as M

OUT = Path("results/08_dual_antigen_escape")
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20260825
pd.set_option("display.width", 260)

# %%
adata = sc.read_h5ad("results/05_integration/integrated.h5ad")
labels = pd.read_csv("results/06_annotation/per_cell_labels.csv.gz", index_col=0)
adata.obs["cell_type"] = labels.reindex(adata.obs_names)["cell_type"].astype(str).values

var_names = list(adata.var_names)
counts = sp.csc_matrix(adata.layers["counts"])
i_bcma, i_gprc = var_names.index("TNFRSF17"), var_names.index("GPRC5D")

bcma = np.asarray(counts[:, i_bcma].todense()).ravel().astype(np.int64)
gprc5d = np.asarray(counts[:, i_gprc].todense()).ravel().astype(np.int64)
total_umi = np.asarray(counts.sum(axis=1)).ravel().astype(np.int64)

obs = adata.obs[["sample_name", "patient_id", "cohort", "sample_type", "cell_type"]].copy()
for col in ("patient_id", "sample_name", "cohort"):
    obs[col] = obs[col].astype(str)
obs["bcma"], obs["gprc5d"], obs["total_umi"] = bcma, gprc5d, total_umi

# %% [markdown]
# ## Depth metric, built to be antigen-independent
#
# Bare `total_counts` would let a perturbation of either antigen shift a cell between
# depth strata — reintroducing exactly the circularity stage 07 was cleared of. Removing
# the two genes' own UMIs closes that, and
# `tests/test_antigen.py::test_depth_metric_invariant_to_antigen_counts` enforces it.

# %%
depth = A.depth_ex_antigen(total_umi, bcma, gprc5d)
obs["depth_ex_antigen"] = depth

# %% [markdown]
# ## Denominators — inherited frozen from stage 07, never redefined here
#
# Clone membership, V thresholds and CNV are not touched. `CLONE_INCOMPATIBLE` and
# `CLONE_UNCERTAIN` contribute no cells to either denominator.

# %%
membership = pd.read_csv(
    "results/07_malignant_plasma/v_clone_membership/clone_membership_per_cell.csv.gz",
    index_col=0)
obs["clone_state"] = membership.reindex(obs.index)["clone_state"].astype(str).values

p7 = pd.read_csv("results/07_malignant_plasma/v_clone_membership/patient_v_evaluability.csv",
                 dtype={"patient": str})
eligible = set(p7[(p7.clonality_state == "CLONAL_STRONG")
                  & (p7.v_state == "V_EVALUABLE")].patient)

in_eligible = obs.patient_id.isin(eligible).values
primary = in_eligible & (obs.clone_state == M.CLONE_SUPPORTED).values
sensitivity = in_eligible & obs.clone_state.isin(
    [M.CLONE_SUPPORTED, M.CLONE_COMPATIBLE_V_UNOBSERVED]).values

assert not obs.clone_state[sensitivity].isin(
    [M.CLONE_INCOMPATIBLE, M.CLONE_UNCERTAIN]).any()
print(f"patients {len(eligible)} | primary {primary.sum():,} | "
      f"sensitivity {sensitivity.sum():,}")

# %% [markdown]
# ## Depth strata — cohort-specific, frozen before any antigen result
#
# Stage 07 showed informative missingness is strongly cohort-dependent (MMRF 17.6× vs
# WU1 1.85× vs WU2 1.35×). A globally pooled bin would place deep WashU cells and shallow
# MMRF cells in one stratum despite entirely different capture regimes, leaving the null
# inadequately conditioned. Edges come from the **primary** denominator only and are
# reused unchanged for the sensitivity denominator, so the two are measured with one ruler.

# %%
def cohort_edges(mask, values, n_bins=5):
    return {c: A.quantile_edges(values[mask & (obs.cohort == c).values], n_bins)
            for c in sorted(obs.cohort[mask].unique())}


def apply_cohort_edges(edges, values):
    out = np.zeros(len(obs), dtype=np.int64)
    for cohort, e in edges.items():
        m = (obs.cohort == cohort).values
        out[m] = A.assign_strata(values[m], e)
    return out


EDGES = cohort_edges(primary, depth)
strata = apply_cohort_edges(EDGES, depth)

# The global scheme is a secondary diagnostic only, never the primary null.
global_edges = A.quantile_edges(depth[primary], 5)
strata_global = A.assign_strata(depth, global_edges)

# %% [markdown]
# ## Observed expression states
#
# Raw counts, `> 0`. No normalization, imputation, smoothing or embedding may touch the
# call — the whole question is whether a transcript is genuinely absent, and imputation
# manufactures expression by borrowing from neighbours, which erases the measurement.

# %%
def patient_metrics(mask, strata_vec, bcma_v, gprc_v, bootstrap=True, n_perm=2000):
    out = []
    for patient in sorted(eligible):
        m = mask & (obs.patient_id == patient).values
        n = int(m.sum())
        if n == 0:
            continue
        b0, g0 = bcma_v[m] == 0, gprc_v[m] == 0
        st = A.merge_sparse_strata(strata_vec[m], 20)
        observed = float((b0 & g0).mean())
        expected = A.stratified_expected_dn(b0, g0, st)
        uncond = A.unconditioned_expected_dn(b0, g0)
        samples = obs.sample_name.values[m]

        rec = {
            "patient": patient, "cohort": obs.cohort.values[m][0], "n_cells": n,
            "n_samples": int(len(np.unique(samples))), "n_strata": int(len(np.unique(st))),
            "bcma_detect": float((bcma_v[m] > 0).mean()),
            "gprc5d_detect": float((gprc_v[m] > 0).mean()),
            "bcma_neg": float(b0.mean()), "gprc5d_neg": float(g0.mean()),
            "obs_double_positive": float(((bcma_v[m] > 0) & (gprc_v[m] > 0)).mean()),
            "obs_BCMA_only": float(((bcma_v[m] > 0) & g0).mean()),
            "obs_GPRC5D_only": float((b0 & (gprc_v[m] > 0)).mean()),
            "observed_double_negative_fraction": observed,
            "expected_dn_stratified": expected,
            "expected_dn_unconditioned": uncond,
            "enrichment_stratified": observed / expected if expected > 0 else np.nan,
            "enrichment_unconditioned": observed / uncond if uncond > 0 else np.nan,
            "excess_dn": observed - expected,
            "median_depth_ex_antigen": float(np.median(depth[m])),
            "low_n": bool(n < 100),
            "single_sample": bool(len(np.unique(samples)) == 1),
        }

        if bootstrap and n >= 20:
            idx = np.arange(n)
            dn_boot = A.hierarchical_bootstrap(
                samples, idx, 2000, SEED, lambda i: float((b0[i] & g0[i]).mean()))

            def ratio(i):
                e = A.stratified_expected_dn(b0[i], g0[i], st[i])
                return float((b0[i] & g0[i]).mean() / e) if e > 0 else np.nan

            enr_boot = A.hierarchical_bootstrap(samples, idx, 2000, SEED, ratio)
            rec["dn_ci_lo"], rec["dn_ci_hi"] = np.percentile(dn_boot, [2.5, 97.5])
            rec["enr_ci_lo"], rec["enr_ci_hi"] = np.nanpercentile(enr_boot, [2.5, 97.5])
            rec["ci_width"] = rec["dn_ci_hi"] - rec["dn_ci_lo"]

        if n_perm and n >= 20:
            null, pval = A.permutation_null_dn(b0, g0, st, n_perm, SEED)
            rec["perm_p"] = pval
            rec["null_dn_lo"], rec["null_dn_hi"] = np.percentile(null, [2.5, 97.5])

        out.append(rec)
    return pd.DataFrame(out)


PRIMARY = patient_metrics(primary, strata, bcma, gprc5d)
SENSITIVITY = patient_metrics(sensitivity, strata, bcma, gprc5d)
PRIMARY.to_csv(OUT / "patient_antigen_states_primary.csv", index=False)
SENSITIVITY.to_csv(OUT / "patient_antigen_states_sensitivity.csv", index=False)

# %% [markdown]
# ## The central comparison
#
# Dropout is a per-*cell* property, so a shallow cell is more likely to read zero for
# *both* genes. Depth heterogeneity alone therefore manufactures positive BCMA⁻/GPRC5D⁻
# association — an artifact pointing in exactly the direction this project hopes to find,
# which is the worst kind. The gap between these two columns *is* that artifact.

# %%
print(PRIMARY[["patient", "cohort", "n_cells",
               "observed_double_negative_fraction",
               "enrichment_unconditioned", "enrichment_stratified",
               "perm_p"]].round(3).to_string(index=False))

# %% [markdown]
# Unconditioned median 1.052 (max **4.606**) collapses to a stratified median of
# **1.009** (max 1.750). Four of 32 patients keep significant enrichment, all of them
# MMRF — the only cohort deep enough for GPRC5D to be detectable at all, so that is where
# the test has power, not necessarily where the biology is. **Reported as the largely
# negative result it is.**
#
# The `Σ_s n_s · p_s(B) · p_s(G)` baseline is *not* a dropout correction: multiplying the
# marginals assumes exactly the independence the co-negativity test exists to interrogate.
# No dropout-corrected DN point estimate is produced, and none is claimed.

# %% [markdown]
# ## Truncate all cohorts at 10,000 UMI — owed since stage 04, not optional
#
# Cells are **downsampled, never discarded**, by exact hypergeometric subsample, so the
# comparison stays over the same cells.

# %%
bcma_t = A.downsample_gene_counts(total_umi, bcma, 10000, SEED)
gprc_t = A.downsample_gene_counts(total_umi, gprc5d, 10000, SEED)
depth_t = A.depth_ex_antigen(np.minimum(total_umi, 10000), bcma_t, gprc_t)
strata_t = apply_cohort_edges(cohort_edges(primary, depth_t), depth_t)

TRUNC = patient_metrics(primary, strata_t, bcma_t, gprc_t, bootstrap=False, n_perm=0)
print(f"cells above 10k UMI: {(total_umi > 10000).sum():,} overall; "
      f"{(total_umi[primary] > 10000).sum():,} of {primary.sum():,} in the primary denominator")

# %% [markdown]
# WashU is unchanged — it was already censored there, which is itself confirmation of the
# deposit ceiling. MMRF moves: GPRC5D detection falls by a mean of 0.186 and DN rises by a
# mean of 0.059 (max +0.182), rank stability Spearman 0.921. So MMRF's advantage in
# GPRC5D detection is substantially a **depth** advantage, and the cohort gap is partly
# what the depositors removed.

# %% [markdown]
# ## Technical-zero floor, from non-plasma cells only
#
# Malignant plasma cells may never define their own null. `Ambiguous` (Leiden 23) is
# excluded **by name** — its lineage is unresolved, and letting it in silently would
# contaminate the reference.

# %%
REFERENCE = ["Tcell", "Myeloid", "Bcell", "HSPC", "NK"]
non_plasma = (obs.cell_type.isin(REFERENCE).values
              & obs.cohort.isin(["MMRF", "WU1", "WU2"]).values)
print(f"reference cells {non_plasma.sum():,}; "
      f"Ambiguous excluded {(obs.cell_type == 'Ambiguous').sum():,}")

# %% [markdown]
# A gene at GPRC5D's abundance reads zero **79%** of the time in the shallowest WashU
# cells and still 37% in the deepest WashU stratum, against a WU2 observed
# GPRC5D-negative rate of 0.961. **The great majority of GPRC5D zeros in WashU are
# consistent with technical dropout.** BCMA's floor falls to 8–10% in deep strata, so the
# two negative calls are not equally reliable and that asymmetry is bounded, not corrected.
#
# Full results, evidence states and limitations: `stage08_escape_summary.md`.

# %%
print((OUT / "stage08_escape_summary.md").read_text())
