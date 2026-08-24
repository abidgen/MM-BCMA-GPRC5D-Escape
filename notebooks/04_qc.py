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
#     display_name: mm-qc
#     language: python
#     name: mm-qc
# ---

# %% [markdown]
# # 04 — Quality control and doublet detection
#
# **Env: `mm-qc`** (this is the only stage that needs R — `scDblFinder` over `rpy2`).
# Writes `results/04_qc/`.
#
# ## What this stage decides, and why it is delicate
#
# Which cells enter the analysis. That matters more here than in a typical scRNA
# project because the headline metric, `frac_double_negative`, is a **fraction of
# zeros**: any filter that preferentially keeps or drops shallow cells moves it
# directly, and in a direction that reads as biology.
#
# So three things are done differently from the default recipe, all settled in
# `CLAUDE.md` before this notebook was written:
#
# 1. **Thresholds are derived per cohort, never pooled.** The four cohorts differ by
#    ~1.9x in genes detected per cell. A pooled MAD would filter WashU cohort 1 — 23
#    of the 54 myeloma samples — harder than MMRF for a reason that is batch, not
#    quality. This notebook computes the pooled result too, so the cost is shown
#    rather than asserted.
# 2. **`sc-best-practices`' numbers are not copied.** Its 5-MAD counts and 8%
#    mitochondrial cap are defaults for healthy PBMC/BMMC. The *procedure* is adopted;
#    the *values* are recomputed here and written to
#    `results/04_qc/qc_thresholds.csv`.
# 3. **The depositors' own QC is discarded.** Their stated Seurat filter (drop cells
#    above 10,000 UMIs as multiplets) is internally garbled and demonstrably was not
#    applied to what is deposited. Checked below rather than asserted.
#
# **Cells are annotated, not deleted.** Every checkpoint keeps all its barcodes with
# `obs["keep"]` set. "Does this result survive a different QC?" is a question that
# will be asked of a fraction-of-zeros metric, and it can only be answered if the
# filtered cells are still on disk.

# %%
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

REPO = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
sys.path.insert(0, str(REPO / "src"))

from mm_escape import config, io, qc  # noqa: E402

warnings.filterwarnings("ignore", category=FutureWarning)
sns.set_theme(style="whitegrid", context="notebook")

OUT = config.RESULTS_DIR / "04_qc"
CHECKPOINTS = OUT / "samples"
FIGURES = OUT / "figures"
for directory in (OUT, CHECKPOINTS, FIGURES):
    directory.mkdir(parents=True, exist_ok=True)
print("writing to", OUT)

# %% [markdown]
# ## The cohort, and the confounder this stage has to respect
#
# `load_manifest` now carries the S1 clinical metadata as well as the GEO metadata,
# so the depth-vs-cohort structure and the patient mapping are both visible here.

# %%
manifest = io.load_manifest()
print(f"{len(manifest)} samples: "
      f"{(manifest['sample_type'] == 'myeloma').sum()} myeloma, "
      f"{(manifest['sample_type'] == 'normal_bm').sum()} normal BM")

in_paper = manifest.loc[manifest["in_paper_cohort"]]
print(f"in-paper cohort: {len(in_paper)} samples / "
      f"{in_paper['patient_id'].nunique()} patients "
      f"(the paper's 53 / 41)")

manifest.groupby("cohort", observed=True).agg(
    n_samples=("sample_name", "size"),
    chemistry=("chemistry", "first"),
    dead_cell_removal=("dead_cell_removal", "first"),
    n_patients=("patient_id", "nunique"),
)

# %% [markdown]
# ## Pass 1 — per sample: QC metrics and doublets
#
# Per sample, and checkpointed, for two independent reasons:
#
# - **scDblFinder must not see pooled samples.** It simulates doublets by combining
#   cells from the same droplet suspension; pooling lets it manufacture
#   cross-sample "doublets" that no droplet could contain, and loses each sample's
#   own doublet rate — which scales with loading density and therefore with
#   chemistry, the exact axis this cohort varies on.
# - **62 samples with an R call each is long enough that a crash must not cost the
#   whole run.** `run_sample_qc` returns the checkpoint if it exists, so re-running
#   this cell is cheap.

# %%
# Resumable: samples already checkpointed are read back rather than recomputed.
obs, thresholds = qc.run_cohort_qc(
    manifest=manifest,
    checkpoint_dir=CHECKPOINTS,
    group_key="cohort",
)
print(f"{len(obs):,} cells over {obs['sample_name'].nunique()} samples")

# %% [markdown]
# ## The deposit is already filtered — differently in each cohort
#
# **This corrects a claim carried in `CLAUDE.md` until this notebook was run.** That
# claim was that the depositors' stated 10,000-UMI cut "demonstrably was not applied
# to what is deposited", reasoned from a cohort-wide average UMI count. The average
# pooled MMRF with WashU and hid a per-cohort truth. Read per cohort, the boundaries
# are sharp and unmistakable — a `max` of exactly 9,999 or 10.00 is a cutoff, not a
# distribution.

# %%
prefilter = pd.DataFrame({
    "n_cells": obs.groupby("cohort", observed=True).size(),
    "min_UMI": obs.groupby("cohort", observed=True)["total_counts"].min(),
    "max_UMI": obs.groupby("cohort", observed=True)["total_counts"].max(),
    "max_pct_mt": obs.groupby("cohort", observed=True)["pct_counts_mt"].max(),
    "min_genes": obs.groupby("cohort", observed=True)["n_genes_by_counts"].min(),
    "pct_above_10k_UMI": obs.groupby("cohort", observed=True)["total_counts"].apply(
        lambda s: 100 * (s > 10_000).mean()),
})
display(prefilter.round(2))

# %% [markdown]
# Read off the table:
#
# | cohort | UMI ceiling | UMI floor | mt ceiling | gene floor |
# |---|---|---|---|---|
# | WU1, WU2 | **< 10,000** | >= 1,000 | < 20% | >= 200 |
# | MMRF | **none** (max ~269,000) | >= 1,000 | **< 10%** | >= 200 |
# | Donor | none | none | < 20% | >= 200 |
#
# ### Why this is a first-order problem, not bookkeeping
#
# Malignant plasma cells are professional secretors and are the highest-RNA-content
# cells in the marrow. So a 10,000-UMI ceiling does not remove a random slice — it
# removes the cells this project is about. Measured in the *uncensored* cohorts,
# where that population is still visible:

# %%
import anndata  # noqa: E402
import scipy.sparse as sp  # noqa: E402

rows = []
for name in ("MMRF_1695", "MMRF_1424", "BM4"):
    sample = anndata.read_h5ad(CHECKPOINTS / f"{name}.h5ad")
    matrix = sample.X.tocsc() if sp.issparse(sample.X) else sp.csc_matrix(sample.X)
    high = (sample.obs["total_counts"] > 10_000).to_numpy()
    for gene in ("TNFRSF17", "GPRC5D", "SDC1", "CD3D"):
        if gene not in sample.var_names:
            continue
        detected = np.asarray(
            matrix[:, sample.var_names.get_loc(gene)].todense()
        ).ravel() > 0
        rows.append({
            "sample": name, "gene": gene,
            "pct_detected_above_10k": 100 * detected[high].mean(),
            "pct_detected_below_10k": 100 * detected[~high].mean(),
            "enrichment": (detected[high].mean() / max(detected[~high].mean(), 1e-9)),
        })
censoring = pd.DataFrame(rows).round(2)
censoring.to_csv(OUT / "umi_censoring_effect.csv", index=False)
display(censoring)

# %% [markdown]
# The band the WashU deposits had removed is enriched several-fold for `TNFRSF17`
# (BCMA) and **20-40x for `GPRC5D`** — the low-abundance transcript whose dropout is
# already this project's largest measurement worry. So the WashU cohorts had the
# antigen-positive tail of their own tumours cut off before deposit, which will
# **inflate `frac_double_negative`** for those 36 of 54 myeloma samples relative to
# MMRF's 18. That is a bias in the project's own direction of interest.
#
# **What is done about it here: nothing, deliberately.** Two options were considered:
#
# - *Harmonise the censoring* — truncate MMRF and Donor at 10,000 UMIs so every
#   cohort is censored identically. Comparable, but it discards 42% of MMRF's cells
#   including most of its plasma cells, to fix a comparability problem by making
#   every cohort equally damaged.
# - *Carry it forward as a stated, quantified confounder.* Chosen. It goes into
#   stage 08 as a covariate alongside cohort and depth, and stage 08 runs the
#   truncate-everything-at-10k version as a **sensitivity analysis** — where it
#   costs nothing to compute and answers the question directly, rather than costing
#   40% of the best cohort's data here.
#
# QC's job is to remove bad cells, not to destroy good ones to equalise two cohorts.

# %% [markdown]
# ### Our own floors are a safety net, not the filter
#
# `min_genes = 200` sits below what the deposit already enforces, so it flags
# nothing. That is intended: a high floor here would hit the shallow cohorts hardest,
# which is precisely the failure per-cohort thresholds exist to avoid.

# %% [markdown]
# ## The depth confounder, measured
#
# This is the reason for every per-cohort decision below and in stage 05.

# %%
depth = (
    obs.groupby(["cohort", "sample_name"], observed=True)["n_genes_by_counts"]
    .median()
    .reset_index()
)
cohort_depth = depth.groupby("cohort", observed=True)["n_genes_by_counts"].agg(
    ["size", "median", "min", "max"]
)
display(cohort_depth)
print(f"spread across cohorts: "
      f"{cohort_depth['median'].max() / cohort_depth['median'].min():.2f}x")

fig, ax = plt.subplots(figsize=(7, 4))
order = cohort_depth.sort_values("median", ascending=False).index
sns.stripplot(data=depth, x="cohort", y="n_genes_by_counts", order=order,
              size=7, alpha=0.8, ax=ax)
sns.boxplot(data=depth, x="cohort", y="n_genes_by_counts", order=order,
            showcaps=False, boxprops={"facecolor": "none"}, showfliers=False, ax=ax)
ax.set(ylabel="median genes / cell (per sample)", xlabel="",
       title="Sequencing depth tracks cohort — the reason QC is not pooled")
fig.tight_layout()
fig.savefig(FIGURES / "depth_by_cohort.png", dpi=150)

# %% [markdown]
# ## The thresholds this cohort actually produces
#
# Not the tutorial's numbers. One row per cohort per metric, with the interval the
# MAD rule implies and how many cells fall each side of it.
#
# `pct_counts_mt` is **one-sided**: a cell with unusually *few* mitochondrial reads
# is not low quality, so only the upper bound is applied. `mt_bound_by` says whether
# the MAD rule or the 20% cap was the binding constraint.

# %%
thresholds.to_csv(OUT / "qc_thresholds.csv", index=False)
thresholds.round(3)

# %% [markdown]
# ### What pooling would have cost
#
# The comparison the per-cohort decision rests on. If the pooled and per-cohort
# thresholds were interchangeable the decision would be cosmetic; they are not.

# %%
_, pooled = qc.cohort_thresholds(obs, group_key=None)
comparison = pd.concat([
    thresholds.assign(scope="per-cohort"),
    pooled.assign(scope="pooled", group="all"),
])[["scope", "group", "metric", "median", "mad", "lower", "upper"]]
display(comparison.query("metric == 'log1p_total_counts'").round(3))

pooled_flags, _ = qc.cohort_thresholds(obs, group_key=None)
side_by_side = pd.DataFrame({
    "per_cohort": obs.groupby("cohort", observed=True)["outlier"].mean() * 100,
    "pooled": pooled_flags.groupby(obs["cohort"], observed=True)["outlier"].mean() * 100,
})
side_by_side["extra_pct_points_if_pooled"] = (
    side_by_side["pooled"] - side_by_side["per_cohort"]
)
print("% of cells flagged as outliers, per cohort:")
display(side_by_side.round(2))

# %% [markdown]
# ## QC metric distributions
#
# Coloured by cohort so the batch structure is visible in the raw metrics rather
# than only in the summary.

# %%
metrics = ["log1p_total_counts", "log1p_n_genes_by_counts",
           "pct_counts_in_top_20_genes", "pct_counts_mt"]
fig, axes = plt.subplots(2, 2, figsize=(12, 7))
for ax, metric in zip(axes.ravel(), metrics):
    for cohort, block in obs.groupby("cohort", observed=True):
        sns.kdeplot(block[metric], ax=ax, label=cohort, fill=False, linewidth=1.6)
    for _, row in thresholds.query("metric == @metric").iterrows():
        if np.isfinite(row["upper"]):
            ax.axvline(row["upper"], ls=":", lw=0.9, color="0.4")
        if np.isfinite(row["lower"]):
            ax.axvline(row["lower"], ls=":", lw=0.9, color="0.4")
    ax.set(title=metric, ylabel="")
axes[0, 0].legend(title="cohort", fontsize=8)
fig.suptitle("QC metrics by cohort; dotted lines are the per-cohort MAD bounds")
fig.tight_layout()
fig.savefig(FIGURES / "qc_metric_distributions.png", dpi=150)

# %% [markdown]
# ## `pct_counts_in_top_20_genes` is computed, but must not delete cells here
#
# The second thing this cohort forced a departure on. In most tissues a library
# dominated by a handful of transcripts means an empty-ish droplet full of ambient
# soup, and a MAD filter on this metric is standard. `CLAUDE.md` specified it as one
# of the three MAD metrics, and it was implemented that way — then the numbers came
# back:

# %%
_, all_flags = qc.cohort_thresholds(obs, group_key="cohort",
                                    filters=(*qc.DEFAULT_FILTERS, "outlier_top20"))
display(all_flags.query("metric == 'pct_counts_in_top_20_genes'")[
    ["group", "median", "mad", "lower", "upper", "n_below", "n_above", "n_cells"]
].round(2))

top20_cost = (
    obs.assign(flagged=obs["outlier_top20"])
    .groupby("cohort", observed=True)["flagged"].mean().mul(100).round(2)
)
print("% of cells a 5-MAD top-20 filter would delete, per cohort:")
display(top20_cost)

# %% [markdown]
# 17% of MMRF and 15% of WashU 1, against 3% of WashU 2 — far too uneven, and far
# too much, to be catching the same thing in each cohort. So: what *are* those cells?

# %%
sample = anndata.read_h5ad(CHECKPOINTS / "MMRF_1695.h5ad")
matrix = sample.X.tocsr() if sp.issparse(sample.X) else sp.csr_matrix(sample.X)
cut = sample.obs["pct_counts_in_top_20_genes"].quantile(0.90)
high = (sample.obs["pct_counts_in_top_20_genes"] > cut).to_numpy()

block = matrix[high]
share = np.asarray(block.sum(0)).ravel() / block.sum()
top = np.argsort(share)[::-1][:8]
print(f"MMRF_1695, top decile of pct_counts_in_top_20_genes (n={high.sum()}):")
for i in top:
    print(f"   {sample.var_names[i]:10s} {100 * share[i]:5.1f}% of the group's counts")

marker_rows = []
for gene in ("TNFRSF17", "GPRC5D", "SDC1", "MZB1", "CD3D", "LYZ"):
    if gene not in sample.var_names:
        continue
    detected = np.asarray(
        matrix[:, sample.var_names.get_loc(gene)].todense()
    ).ravel() > 0
    marker_rows.append({
        "gene": gene,
        "pct_in_high_top20": 100 * detected[high].mean(),
        "pct_elsewhere": 100 * detected[~high].mean(),
    })
display(pd.DataFrame(marker_rows).round(2))

# %% [markdown]
# Two populations, not one: `IGKC` at ~25% of counts (**plasma cells**) and
# `HBB`/`HBA1`/`HBA2` at ~32% (**erythroid debris**). And the plasma-cell half is the
# project's entire subject — `TNFRSF17` (BCMA) is detected in ~22% of that decile
# against ~1% elsewhere, `SDC1` (CD138) in ~19% against ~0%.
#
# A plasma cell is a professional secretor. A library dominated by immunoglobulin is
# its **normal state**, not a defect. Filtering on this metric would preferentially
# delete antigen-*positive* malignant plasma cells and inflate
# `frac_double_negative` — the same direction as the censoring problem above, and the
# same reason it is unacceptable.
#
# **Decision: `outlier_top20` is computed and reported, but is not in
# `qc.DEFAULT_FILTERS`.** The metric still earns its keep — it is one of the few
# handles on ambient Ig available at all here, since SoupX/DecontX need unfiltered
# matrices this deposit does not have, and stage 08 needs one. It is just not allowed
# to delete cells. The erythroid half is what `pct_counts_hb` is for, and properly
# belongs to stage 06's annotation rather than to a blanket QC cut.
#
# This is a documented deviation from `CLAUDE.md`'s stage-04 specification, made on
# this cohort's own numbers — which is exactly what that specification asked for
# when it said to recompute rather than copy.

# %%
print("filters composing obs['outlier']:", qc.DEFAULT_FILTERS)
print("computed but not filtering:     ",
      tuple(f for f in qc.ALL_FLAGS if f not in qc.DEFAULT_FILTERS))

# %% [markdown]
# ## Sensitivity of the outlier rate to the MAD count
#
# The default is the standard 5 MADs, adopted as a procedure rather than tuned. This
# is what tightening or loosening it would do — recorded so the choice is auditable
# and so a later stage can say how much of a result depends on it.

# %%
sweep = []
for n_mads in (3.0, 4.0, 5.0, 6.0):
    flags, _ = qc.cohort_thresholds(
        obs, group_key="cohort",
        n_mads={metric: n_mads for metric in qc.MAD_METRICS},
    )
    removed = flags["outlier"] | obs["is_doublet"].to_numpy()
    sweep.append({
        "n_mads": n_mads,
        "pct_removed_overall": 100 * removed.mean(),
        **{f"pct_removed_{cohort}": 100 * removed[(obs["cohort"] == cohort).to_numpy()].mean()
           for cohort in sorted(obs["cohort"].unique())},
    })
sweep = pd.DataFrame(sweep)
sweep.to_csv(OUT / "mad_sensitivity.csv", index=False)
sweep.round(2)

# %% [markdown]
# ## Doublets
#
# `scDblFinder` per sample. The expected pattern is that the doublet rate tracks
# loading density, so it should differ by chemistry/cohort rather than being flat —
# a flat rate across a cohort this heterogeneous would suggest the caller is not
# seeing the samples separately.

# %%
doublet_rate = (
    obs.groupby(["cohort", "sample_name"], observed=True)["is_doublet"]
    .mean().mul(100).reset_index()
)
display(doublet_rate.groupby("cohort", observed=True)["is_doublet"]
        .agg(["size", "median", "min", "max"]).round(2))

fig, ax = plt.subplots(figsize=(7, 4))
sns.stripplot(data=doublet_rate, x="cohort", y="is_doublet", size=7, alpha=0.8, ax=ax)
ax.set(ylabel="% doublets (per sample)", xlabel="",
       title="scDblFinder doublet rate by cohort")
fig.tight_layout()
fig.savefig(FIGURES / "doublet_rate_by_cohort.png", dpi=150)

# %% [markdown]
# ## What QC did, per sample
#
# The per-filter columns overlap — a cell can be flagged by several — so they do not
# sum to `n_removed`, which is their union and is the honest total.

# %%
report = qc.qc_report(obs, by="sample_name")
report.to_csv(OUT / "qc_report_by_sample.csv", index=False)

by_cohort = qc.qc_report(obs, by="cohort")
by_cohort.to_csv(OUT / "qc_report_by_cohort.csv", index=False)
display(by_cohort.round(2))

print(f"kept {int(report['n_kept'].sum()):,} of {int(report['n_cells_pre'].sum()):,} "
      f"cells ({100 * report['n_kept'].sum() / report['n_cells_pre'].sum():.1f}%)")
report.sort_values("pct_removed", ascending=False).head(10).round(2)

# %% [markdown]
# ### Did QC flatten the depth gap, or preserve it?
#
# Both answers are informative. Preserved means QC left the biological/batch
# structure alone, which is what per-cohort thresholds are supposed to do. Flattened
# would mean QC itself became the batch correction — and a filter that equalises
# depth across cohorts by discarding the shallow cohort's cells is exactly the
# failure this stage is designed around.

# %%
depth_shift = pd.DataFrame({
    "median_genes_pre": obs.groupby("cohort", observed=True)["n_genes_by_counts"].median(),
    "median_genes_post": obs.loc[obs["keep"]].groupby("cohort", observed=True)[
        "n_genes_by_counts"].median(),
    "pct_removed": obs.groupby("cohort", observed=True)["keep"].apply(
        lambda s: 100 * (~s).mean()),
})
depth_shift["spread_pre"] = (depth_shift["median_genes_pre"].max()
                             / depth_shift["median_genes_pre"].min())
depth_shift["spread_post"] = (depth_shift["median_genes_post"].max()
                              / depth_shift["median_genes_post"].min())
display(depth_shift.round(2))

# %% [markdown]
# ## Normal-donor controls survive QC
#
# The 8 donor samples are stage 07's negative control (polyclonal marrow must yield
# no malignant clone) and stage 09's normal-plasma-cell baseline. They are the
# shallowest cohort after WashU 1, so a pooled QC would have hit them hardest — this
# checks per-cohort thresholds actually protected them.

# %%
donors = obs.loc[obs["sample_type"] == "normal_bm"]
print(f"{len(donors):,} donor cells, {donors['keep'].mean() * 100:.1f}% kept, "
      f"across {donors['sample_name'].nunique()} samples")
display(donors.groupby("sample_name", observed=True).agg(
    n_cells=("keep", "size"),
    pct_kept=("keep", lambda s: 100 * s.mean()),
    median_genes=("n_genes_by_counts", "median"),
).round(1))

# %% [markdown]
# ## Post-QC cohort summary
#
# The numbers stage 05 starts from. `results/04_qc/samples/*.h5ad` holds one
# checkpoint per sample with every barcode retained and `obs["keep"]` set; stage 05
# does the gene-space intersection over those files.

# %%
summary = pd.DataFrame({
    "n_cells_pre": obs.groupby("cohort", observed=True).size(),
    "n_cells_post": obs.loc[obs["keep"]].groupby("cohort", observed=True).size(),
    "n_samples": obs.groupby("cohort", observed=True)["sample_name"].nunique(),
})
summary["pct_kept"] = 100 * summary["n_cells_post"] / summary["n_cells_pre"]
summary.loc["TOTAL"] = [summary["n_cells_pre"].sum(), summary["n_cells_post"].sum(),
                        summary["n_samples"].sum(),
                        100 * summary["n_cells_post"].sum() / summary["n_cells_pre"].sum()]
display(summary.round(1))

# The per-cell QC table is deliberately NOT written out. Every cell's QC metrics and
# flags already live in its sample's checkpoint, and `qc.collect_obs()` rebuilds the
# pooled table from them in ~10 s by reading `obs` alone. A separate 200,040-row copy
# would be a redundant intermediate large enough to matter in git — and the committed
# deliverables of this stage are the CSVs and figures below, not the raw table.
deliverables = sorted(p.name for p in OUT.glob("*.csv"))
print("deliverables in", OUT)
for name in deliverables:
    print("   ", name)
for name in sorted(p.name for p in FIGURES.glob("*.png")):
    print("    figures/" + name)
print(f"\ncheckpoints: {len(list(CHECKPOINTS.glob('*.h5ad')))} files in {CHECKPOINTS}")
print("rebuild the pooled QC table any time with qc.collect_obs()")
