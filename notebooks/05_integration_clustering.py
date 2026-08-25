# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.5
#   kernelspec:
#     display_name: mm-core
#     language: python
#     name: mm-core
# ---

# %% [markdown]
# # 05 — Gene-space harmonization, integration, clustering
#
# **Env: `mm-core`.** Reads the stage-04 checkpoints, writes `results/05_integration/`.
#
# ## What this stage is for, and what it must not do
#
# It produces one harmonized object over every QC-passing cell, with a batch-corrected
# embedding and Leiden clusters. Stage 06 annotates it; stage 07 subsets it.
#
# The delicate part is not the clustering — it is **constraining what integration is
# allowed to touch**. Harmony keyed on `patient_id` is right for the immune
# compartment and carries a real risk for the tumour: the malignant clone is
# *patient-private by definition*, so forcing patients together can blend genuinely
# distinct clones into one blob and erase the heterogeneity this project exists to
# measure. So, stated here and enforced downstream:
#
# - The integrated embedding is for **immune-compartment annotation and clustering
#   only** (stages 06 and 11).
# - **All malignant subclustering is per patient and un-integrated** (stage 10). It
#   must not read `obsm["X_pca_harmony"]`.
# - **Per-cell antigen calls are raw counts** (`layers["counts"]`) and are therefore
#   integration-independent. This is what contains the risk, and it is the answer to
#   *"did Harmony distort your escape fractions?"* — it cannot, because the calls
#   never touch the embedding.
#
# ## And the gene space, which is where the real correctness risk lives
#
# Two Cell Ranger references are in play with different HGNC symbol vintages. A
# symbol join keeps 22,164 genes and silently **mis-pairs** some of them (`TBCE` is a
# different Ensembl entry in each build). The ID join keeps 32,991 and pairs them
# correctly. `gene_space.py` does the work and verifies itself position-for-position;
# this notebook drives it and reports what it recovered.

# %%
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

REPO = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
sys.path.insert(0, str(REPO / "src"))

from mm_escape import config, integration  # noqa: E402

warnings.filterwarnings("ignore", category=FutureWarning)
sns.set_theme(style="whitegrid", context="notebook")
sc.settings.verbosity = 1

OUT = config.RESULTS_DIR / "05_integration"
FIGURES = OUT / "figures"
CHECKPOINTS = config.RESULTS_DIR / "04_qc" / "samples"
for directory in (OUT, FIGURES):
    directory.mkdir(parents=True, exist_ok=True)
INTEGRATED = OUT / "integrated.h5ad"
print("writing to", OUT)

# %% [markdown]
# ## Load the stage-04 checkpoints and filter
#
# **Cells are filtered here, not at stage 04.** The checkpoints deliberately hold
# every barcode with `obs["keep"]` set, so "does this result survive a different QC?"
# stays answerable for a metric that is a fraction of zeros. `keep_only=False` loads
# everything for exactly such a re-run.
#
# One object per sample, not a concatenation: `attach_ensembl_ids` must run **per
# sample, before any concat**, because it verifies the deposited symbol column
# position-for-position against the committed map for *that sample's build*. After a
# concat there is no single build left to verify against.

# %%
if INTEGRATED.exists():
    print(f"resuming from {INTEGRATED}")
    adata = anndata.read_h5ad(INTEGRATED)
else:
    blocks = integration.load_qc_checkpoints(CHECKPOINTS)
    adata = integration.build_gene_space(blocks)
    del blocks
    integration.normalize_and_hvg(adata)
    integration.run_pca_harmony(adata)
    integration.cluster_and_embed(adata)
    adata.write_h5ad(INTEGRATED)

adata

# %% [markdown]
# ## What the Ensembl-ID join bought
#
# Reported rather than assumed, because a regression here would otherwise be silent —
# and it is worth ~10,800 genes.

# %%
recovery = pd.DataFrame([
    {"join key": "raw deposited symbols", "genes retained": 22_164},
    {"join key": "symbols + the 4-gene alias map", "genes retained": 22_168},
    {"join key": "Ensembl IDs (used)", "genes retained": adata.n_vars},
])
recovery["vs symbol join"] = recovery["genes retained"] - 22_164
display(recovery)

drifted = int(adata.var["symbol_drift"].sum())
print(f"{drifted:,} of the {adata.n_vars:,} retained genes carry a DIFFERENT symbol "
      f"in each build,\nand were therefore invisible to a symbol join.")

# %% [markdown]
# The four that this project's gene lists actually depend on — the ones the required-
# gene assertions exist to catch. Without harmonization `NSD2` is dropped and stage 10
# loses t(4;14), the highest-risk MM translocation, entirely.

# %%
watch = ["NSD2", "TENT5C", "NSD3", "ATP5F1A"]
display(adata.var.loc[[g for g in watch if g in adata.var_names],
                      ["ensembl_id", "symbol_33538", "symbol_33694", "symbol_drift"]])

# The panel every downstream stage depends on. gene_space.assert_required_genes()
# already hard-failed if any of these were missing; this is the visible receipt.
panel = ["TNFRSF17", "GPRC5D", "SLAMF7", "FCRL5", "SDC1", "CD38", "ITGB7", "NCSTN",
         "IGKC", "MZB1", "XBP1", "IRF4", "MS4A1", "CD3D", "NKG7", "LYZ", "HBB", "CD34"]
missing = [g for g in panel if g not in adata.var_names]
print(f"required-gene spot check: {len(panel) - len(missing)}/{len(panel)} present"
      + (f" — MISSING {missing}" if missing else ""))

# %% [markdown]
# ## What went in
#
# The 8 donors and `25183` are **kept in the embedding**, deliberately. They have no
# S1 patient, so Harmony treats each as its own batch — but they are stage 07's
# negative control (polyclonal marrow must yield no clone) and stage 09's normal
# plasma-cell baseline, and both need to sit in a shared space with the myeloma cells
# to be comparable at all. They are excluded at the *patient-level aggregation* in
# stage 08 instead, via `in_paper_cohort` / `sample_type`, which is where the
# distinction actually matters.

# %%
composition = pd.DataFrame({
    "n_cells": adata.obs.groupby("cohort", observed=True).size(),
    "n_samples": adata.obs.groupby("cohort", observed=True)["sample_name"].nunique(),
    "n_patients": adata.obs.groupby("cohort", observed=True)["patient_id"].nunique(),
    "median_genes": adata.obs.groupby("cohort", observed=True)["n_genes_by_counts"].median(),
})
composition.loc["TOTAL"] = [composition["n_cells"].sum(), composition["n_samples"].sum(),
                            adata.obs["patient_id"].nunique(),
                            adata.obs["n_genes_by_counts"].median()]
display(composition.round(1))
print("Harmony keys:", adata.uns["harmony_keys"])
print("HVGs:", int(adata.var['highly_variable'].sum()),
      "| selected within patient, so the choice is not driven by cohort depth")

# %% [markdown]
# ## Did integration actually mix the batches?
#
# The diagnostic UMAPs. `n_genes_ref` is the reference-build split and `cohort` is the
# collection/chemistry axis — **different axes, which is why both are Harmony
# covariates**: two WU1 samples sit on the 33538 build and the four `ND_*` donors on
# 33694, so neither substitutes for the other.

# %%
fig, axes = plt.subplots(2, 3, figsize=(18, 10))
for ax, colour in zip(axes.ravel(),
                      ["leiden", "cohort", "n_genes_ref",
                       "sample_type", "chemistry", "log1p_total_counts"]):
    sc.pl.umap(adata, color=colour, ax=ax, show=False, frameon=False,
               legend_loc="on data" if colour == "leiden" else "right margin",
               legend_fontsize=6, size=2)
    ax.set_title(colour)
fig.tight_layout()
fig.savefig(FIGURES / "umap_diagnostics.png", dpi=150)

# %% [markdown]
# ### Batch mixing per cluster, quantified
#
# Normalized entropy of each cluster's batch composition: 1.0 means the cluster's mix
# matches the cohort's overall mix, 0.0 means it is one batch only.
#
# **A low value is not automatically a failure**, which is exactly why this is
# reported per cluster instead of as one number. A patient-private malignant clone
# *should* be low-entropy on `patient_id` — that is the biology, not a defect. Read it
# against `cohort` and `n_genes_ref`, where low entropy has no biological excuse.

# %%
mixing = {}
for key in ("cohort", "n_genes_ref", "patient_id"):
    table = integration.batch_mixing(adata, batch_key=key)
    table.insert(1, "batch_key", key)
    mixing[key] = table
    table.to_csv(OUT / f"batch_mixing_{key}.csv", index=False)

summary = pd.DataFrame({
    key: {
        "median entropy": table["entropy"].median(),
        "min entropy": table["entropy"].min(),
        "clusters below 0.5": int((table["entropy"] < 0.5).sum()),
        "n clusters": len(table),
    }
    for key, table in mixing.items()
}).T
display(summary.round(3))

print("Clusters least mixed by COHORT (the ones to look at — batch has no excuse "
      "here,\nbut a patient-private clone concentrated in one cohort would also "
      "show up):")
display(mixing["cohort"].nsmallest(6, "entropy")[
    ["leiden", "n_cells", "entropy", "dominant", "dominant_pct"]].round(3))

# %% [markdown]
# ### Harmony vs. no correction
#
# The honest comparison: the same clustering diagnostic on the uncorrected PCA. If
# the two were similar, Harmony would be doing nothing and should be dropped.

# %%
# A lightweight shell carrying only `obs` and the uncorrected embedding — copying
# the full object would duplicate a ~6 GB matrix to recluster 50 PCs.
uncorrected = anndata.AnnData(
    X=np.zeros((adata.n_obs, 0), dtype=np.float32),
    obs=adata.obs.copy(),
    obsm={"X_pca": adata.obsm["X_pca"].copy()},
)
integration.cluster_and_embed(uncorrected, use_rep="X_pca", key_added="leiden_raw")
raw_mixing = integration.batch_mixing(uncorrected, batch_key="cohort",
                                      cluster_key="leiden_raw")

comparison = pd.DataFrame({
    "uncorrected (X_pca)": {
        "median cluster entropy by cohort": raw_mixing["entropy"].median(),
        "clusters below 0.5": int((raw_mixing["entropy"] < 0.5).sum()),
        "n clusters": len(raw_mixing),
    },
    "Harmony (X_pca_harmony)": {
        "median cluster entropy by cohort": mixing["cohort"]["entropy"].median(),
        "clusters below 0.5": int((mixing["cohort"]["entropy"] < 0.5).sum()),
        "n clusters": len(mixing["cohort"]),
    },
}).T
display(comparison.round(3))
comparison.to_csv(OUT / "harmony_vs_uncorrected.csv")

fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
sc.pl.umap(uncorrected, color="cohort", ax=axes[0], show=False, frameon=False,
           size=2, title="uncorrected PCA")
sc.pl.umap(adata, color="cohort", ax=axes[1], show=False, frameon=False,
           size=2, title="Harmony")
fig.tight_layout()
fig.savefig(FIGURES / "umap_harmony_vs_uncorrected.png", dpi=150)
del uncorrected

# %% [markdown]
# ### Where the correction worked, and where it did not
#
# The median above hides the most important thing in this stage. Splitting the
# clusters by whether they are plasma-cell-like separates two very different results.

# %%
import scipy.sparse as sp  # noqa: E402

counts = adata.X.tocsc()


def detected_pct(gene: str, mask: np.ndarray) -> float:
    if gene not in adata.var_names:
        return float("nan")
    values = np.asarray(counts[:, adata.var_names.get_loc(gene)].todense()).ravel() > 0
    return 100 * values[mask].mean()


rows = []
cohort_mix = mixing["cohort"].set_index(mixing["cohort"]["leiden"].astype(str))
for cluster in adata.obs["leiden"].cat.categories:
    mask = (adata.obs["leiden"] == cluster).to_numpy()
    rows.append({
        "leiden": cluster,
        "n_cells": int(mask.sum()),
        "n_patients": adata.obs.loc[mask, "patient_id"].nunique(),
        "entropy_cohort": cohort_mix.loc[str(cluster), "entropy"],
        "dominant_cohort": cohort_mix.loc[str(cluster), "dominant"],
        **{gene: round(detected_pct(gene, mask), 1)
           for gene in ("MZB1", "SDC1", "TNFRSF17", "CD3D", "NKG7", "LYZ", "MS4A1")},
    })
profile = pd.DataFrame(rows)
profile["plasma_like"] = profile["MZB1"] > 40
profile.to_csv(OUT / "cluster_profile.csv", index=False)

split = profile.groupby("plasma_like").agg(
    n_clusters=("leiden", "size"),
    n_cells=("n_cells", "sum"),
    median_entropy=("entropy_cohort", "median"),
).rename(index={True: "plasma-cell-like (MZB1 > 40%)", False: "everything else"})
display(split.round(3))

# %% [markdown]
# **Harmony corrected the immune compartment and did not correct the plasma-cell
# compartment** — a ~7× gap in cohort mixing between the two. The three largest
# plasma-cell clusters are one per cohort, each spanning ~30 patients, rather than
# one shared cluster.

# %%
display(profile.query("plasma_like").sort_values("n_cells", ascending=False)[
    ["leiden", "n_cells", "n_patients", "dominant_cohort", "entropy_cohort",
     "MZB1", "SDC1", "TNFRSF17"]].round(3))

# %% [markdown]
# ### Is that a failure? Partly, and it is the stage-04 censoring showing up again
#
# Two readings, and they are not equally likely:
#
# - *It is biology.* The malignant clone is patient-private, so plasma cells should
#   not merge across patients. **This does not explain what is observed**: the split
#   is by **cohort**, and each of the three clusters spans ~30 patients. A
#   patient-private effect would fragment into ~41 clusters, not three.
# - *It is the censoring found at stage 04.* WashU cohorts 1 and 2 were cut at 10,000
#   UMIs before deposit; MMRF was not. Plasma cells are the highest-RNA-content cells
#   in marrow, so **WashU's plasma cells are a truncated subset of the plasma-cell
#   distribution** — the antigen-positive, high-depth tail is simply missing. They are
#   genuinely a different population from MMRF's, and no batch-correction method can
#   restore cells that were never deposited.
#
# The second reading also explains why the effect is specific to this compartment:
# T, NK, myeloid and B cells sit well below 10,000 UMIs in every cohort, so the
# ceiling never touched them, and Harmony mixes them to entropy ~0.75.
#
# ### What follows from it
#
# This is **contained rather than fatal**, and it is contained by decisions made
# before it was observed:
#
# 1. **Per-cell antigen calls are raw counts** (`layers["counts"]`) and never touch
#    this embedding, so stage 08's numbers are unaffected by how well Harmony did.
# 2. **Malignant subclustering at stage 10 is per patient and un-integrated**, so it
#    does not consume this embedding either.
# 3. **Stage 06 annotates at cluster level**, and three cohort-specific plasma-cell
#    clusters all annotate as PlasmaCell. Splitting costs nothing there.
#
# What it does mean: **no cross-cohort comparison of malignant-cell state may be read
# off this embedding**, and stage 08's cohort covariate is not optional. The
# truncate-all-cohorts-at-10,000 sensitivity analysis that stage 04 said was owed is
# now owed twice over — this is the second independent sign of the same problem.

# %% [markdown]
# ## Where the plasma cells are
#
# Not annotation — that is stage 06's job, with three methods compared and thresholds
# declared in advance. This is a sanity check that the object contains a recognisable
# plasma-cell population at all, and that the antigen genes are detected somewhere
# sensible. If BCMA and GPRC5D did not concentrate in an `MZB1`/`SDC1`-high region,
# something upstream would be wrong and stage 06 should not be started.

# %%
fig, axes = plt.subplots(2, 4, figsize=(20, 9))
for ax, gene in zip(axes.ravel(),
                    ["MZB1", "SDC1", "TNFRSF17", "GPRC5D",
                     "CD3D", "NKG7", "LYZ", "HBB"]):
    sc.pl.umap(adata, color=gene, ax=ax, show=False, frameon=False, size=2,
               cmap="viridis")
fig.tight_layout()
fig.savefig(FIGURES / "umap_marker_sanity.png", dpi=150)

# %%
markers = {
    "PlasmaCell": ["SDC1", "CD38", "MZB1", "XBP1", "IRF4"],
    "Antigen": ["TNFRSF17", "GPRC5D", "SLAMF7", "FCRL5"],
    "Bcell": ["MS4A1", "CD79A", "CD19"],
    "Tcell": ["CD3D", "CD3E", "CD8A"],
    "NK": ["NKG7", "GNLY"],
    "Myeloid": ["CD14", "LYZ"],
    "Erythroid": ["HBB", "GYPA"],
    "HSPC": ["CD34", "KIT"],
}
present = {k: [g for g in v if g in adata.var_names] for k, v in markers.items()}
sc.pl.dotplot(adata, present, groupby="leiden", standard_scale="var",
              show=False, figsize=(16, 7))
plt.savefig(FIGURES / "dotplot_markers_by_cluster.png", dpi=150, bbox_inches="tight")

# %% [markdown]
# ## Per-sample cluster composition
#
# Proportions, not counts: sample cell yields vary ~15x across this cohort, so a raw
# count table would read as biology. This is the input to stage 06's composition
# outputs and stage 11's confounder control.

# %%
composition_by_sample = integration.composition_table(adata, by="sample_name")
composition_by_sample.to_csv(OUT / "cluster_composition_by_sample.csv")

fig, ax = plt.subplots(figsize=(14, 6))
order = adata.obs.groupby("sample_name", observed=True)["cohort"].first().sort_values().index
sns.heatmap(composition_by_sample.loc[order], cmap="rocket_r", ax=ax,
            cbar_kws={"label": "fraction of the sample's cells"})
ax.set(xlabel="leiden cluster", ylabel="", title="Cluster composition per sample")
ax.tick_params(axis="y", labelsize=6)
fig.tight_layout()
fig.savefig(FIGURES / "composition_by_sample.png", dpi=150)

# %% [markdown]
# ## What stage 06 receives

# %%
print(f"{adata.n_obs:,} cells x {adata.n_vars:,} genes")
print(f"{adata.obs['sample_name'].nunique()} samples, "
      f"{adata.obs['patient_id'].nunique()} patients, "
      f"{adata.obs['leiden'].nunique()} Leiden clusters")
print(f"\nembeddings : X_pca (uncorrected), X_pca_harmony (use this), X_umap")
print(f"layers     : counts  <- RAW counts; stage 08 reads these and nothing else")
print(f"X          : log1p(CP10K), all {adata.n_vars:,} genes (no .raw — X is not subset)")
print(f"\nwrote {INTEGRATED}")
print("\nREMINDER for stage 10: malignant subclustering is per patient and")
print("un-integrated. Do NOT read obsm['X_pca_harmony'] there.")
