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
#     display_name: mm-integration
#     language: python
#     name: mm-integration
# ---

# %% [markdown]
# # 05b — Integration-method benchmark
#
# **Env: `mm-integration`.** Reads `results/05_integration/integrated.h5ad`
# **read-only**; writes `results/05b_benchmark/`.
#
# ## Why `05b` and not `06`
#
# `CLAUDE.md`'s rule is that number order is execution order with no exceptions. This
# is a side-comparison feeding the stage-05 *choice*, not a new pipeline stage, so it
# takes a letter rather than displacing annotation. It also **must not overwrite the
# incumbent** — candidate embeddings live here until a winner is picked, and the last
# cell asserts `integrated.h5ad` is byte-identical afterwards.
#
# ## Why this exists
#
# Stage 05 used Harmony because it was the obvious default, never because it beat
# anything. `sc-best-practices`' integration chapter recommends running several methods
# and scoring them with scIB rather than assuming. This closes that gap.
#
# ## Why the scoring is deliberately not standard scIB
#
# Stage 04 found the deposit is censored per cohort: WashU 1 and 2 were cut at 10,000
# UMIs before deposit, MMRF and the donors were not. Plasma cells are the
# highest-RNA-content cells in marrow, so MMRF's two largest plasma clusters hold 68%
# and 88% of their cells above that ceiling.
#
# **scIB's batch metrics cannot distinguish "correctly left apart" from "failed to
# merge".** A method that squashes the three plasma islands together scores better on
# kBET/iLISI while manufacturing correspondence where none is recoverable. A naive
# global ranking would structurally reward overcorrection on this dataset.
#
# Note the careful claim: this is **not** an assertion that the cohorts have
# biologically different plasma cells. It is a **non-recoverable sampling/censoring
# asymmetry** — WashU's *observed* distribution is missing its high-RNA portion, so no
# one-to-one correspondence remains for any method to recover. That argument is
# stronger precisely because it does not depend on the difference being biological.
#
# So: **the immune compartment is scored; the plasma compartment is diagnosed.**
#
# ## What this cannot do, stated before any result
#
# **No integration method restores cells that were never deposited.** Whichever arm
# wins, the ascertainment bias remains in the raw counts stage 08 reads, and stage 08
# still owes its truncate-all-cohorts-at-10,000 sensitivity analysis. And the blast
# radius is modest by construction: the embedding feeds only stages 06 and 11, antigen
# calls read `layers["counts"]`, and malignant subclustering is per-patient and
# un-integrated — so **nothing here can move `frac_double_negative`.**

# %%
from __future__ import annotations

import sys
import time
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

from mm_escape import benchmark, config, integration  # noqa: E402

warnings.filterwarnings("ignore", category=FutureWarning)
sns.set_theme(style="whitegrid", context="notebook")
sc.settings.verbosity = 1
SEED = 0

INTEGRATED = config.RESULTS_DIR / "05_integration" / "integrated.h5ad"
OUT = config.RESULTS_DIR / "05b_benchmark"
FIGURES = OUT / "figures"
for directory in (OUT, FIGURES):
    directory.mkdir(parents=True, exist_ok=True)

# Fingerprint the incumbent so the last cell can prove this notebook did not touch it.
INCUMBENT_FINGERPRINT = (INTEGRATED.stat().st_size, INTEGRATED.stat().st_mtime_ns)
print("reading (read-only):", INTEGRATED)
print("writing to:", OUT)

# %%
import torch  # noqa: E402

print("torch", torch.__version__, "| CUDA available:", torch.cuda.is_available())
if torch.cuda.is_available():
    print("       ", torch.cuda.get_device_name(0))
else:
    print("        scVI will run on CPU — slower, but the result is the same.")

# %% [markdown]
# ## Load, and fix what every arm shares
#
# Every arm must see **identical cells, gene space, HVGs, PCs and seed**, or the
# comparison measures five different pipelines rather than five methods.

# %%
adata = anndata.read_h5ad(INTEGRATED)
print(adata.shape, "|", adata.obs["sample_name"].nunique(), "samples,",
      adata.obs["patient_id"].nunique(), "patients,",
      adata.obs["cohort"].nunique(), "cohorts")
print("HVGs:", int(adata.var["highly_variable"].sum()))
assert "counts" in adata.layers, "scVI needs raw counts"

# %% [markdown]
# ## Provisional labels — the non-circular half of the scoring
#
# scIB scores batch removal **against** biological conservation, and the bio half needs
# cell-type labels. Ours come from stage 06 — which consumes the embedding this
# benchmark is choosing. CellTypist breaks the circle: it classifies from
# log-normalized expression and never sees an embedding.
#
# **`majority_voting=False` is load-bearing.** Majority voting smooths predictions over
# an over-clustering, and that clustering comes from a representation — which would
# smuggle an embedding straight back into the labels the embeddings are scored against.
#
# These labels are **provisional and for scoring only.** Stage 06 still runs its full
# three-method comparison; nothing downstream reads these.

# %%
labels = benchmark.provisional_labels(adata)
labels.to_csv(OUT / "provisional_labels.csv")

summary = pd.DataFrame({
    "n_cells": adata.obs.groupby("broad_label", observed=True).size(),
    "compartment": adata.obs.groupby("broad_label", observed=True)["compartment"].first(),
}).sort_values("n_cells", ascending=False)
display(summary)

immune = (adata.obs["compartment"] == "immune").to_numpy()
plasma = (adata.obs["compartment"] == "plasma").to_numpy()
print(f"\nimmune {immune.sum():,} | plasma {plasma.sum():,} | "
      f"other/unmapped {(~immune & ~plasma).sum():,} (excluded from scoring)")

# %% [markdown]
# ## The shared PCA
#
# Computed once. The unintegrated baseline **is** these PCs, and every Harmony arm
# corrects them — so the baseline differs from the corrected arms by batch correction
# alone, not by preprocessing.

# %%
pca = benchmark.shared_pca(adata, n_comps=50, random_state=SEED)
print("shared PCA:", pca.shape)

# %% [markdown]
# ## Run the arms
#
# One common batch key (`sample_name`, the true technical unit) so the comparison is of
# *methods*, plus reference arms for the incumbent's own configuration and for
# `cohort`.
#
# **Why `cohort` is a real arm and not an afterthought:** 42 of 50 patients contribute
# exactly one sample, so `sample_name` and `patient_id` are nearly the same partition —
# "avoid correcting on patient" is a weaker argument than it first appears, and
# correcting on either may erase genuine between-patient immune biology. `cohort` is
# the only batch definition **not** confounded with patient, and it is where the
# demonstrated distortion actually lives.
#
# **BBKNN is absent on purpose** — it yields a corrected neighbour graph and no
# embedding, so `Benchmarker` (which scores `obsm` keys) cannot place it on the same
# footing. Deriving an embedding from its graph would compare a different kind of
# object. It is installed in the env so a graph-only side diagnostic stays possible.

# %%
timings = {}
for arm in benchmark.ARMS:
    key = f"X_{arm.name}"
    if key in adata.obsm:
        print(f"  {arm.name:20s} cached")
        continue
    start = time.time()
    adata.obsm[key] = benchmark.run_arm(
        adata, arm, pca=pca, random_state=SEED
    ).astype(np.float32)
    timings[arm.name] = time.time() - start
    print(f"  {arm.name:20s} {str(arm.batch_keys):45s} "
          f"{adata.obsm[key].shape} {timings[arm.name]:6.0f}s", flush=True)

ARM_KEYS = [f"X_{arm.name}" for arm in benchmark.ARMS]
pd.Series(timings, name="seconds").to_csv(OUT / "arm_runtimes.csv")

# %% [markdown]
# ## Primary scoring — scIB on the immune compartment
#
# Where mixing genuinely should happen, and where the embedding is actually used
# (stages 06 and 11).
#
# `Benchmarker` runs its own standardized clustering for NMI/ARI, so clustering
# parameters are consistent across arms for free; only the seed needs pinning.
#
# **Every arm is scored on `cohort`, whatever it corrected on — and that is deliberate.**
# The arms use three different batch definitions (`sample_name`, the incumbent's three
# covariates, `cohort`), so scoring each against its own key would compare six different
# questions rather than six methods. `cohort` is the common axis: samples nest inside
# it, correcting a finer key should also mix it, and it is where the demonstrated
# distortion — the chemistry, site and 10,000-UMI censoring — actually lives.

# %%
from scib_metrics.benchmark import Benchmarker  # noqa: E402

MAX_CELLS = 40_000  # declared in advance; see the plan's runtime note


def score(mask: np.ndarray, tag: str) -> tuple[pd.DataFrame, int | None]:
    """scIB over the arms, on `mask`, subsampling only if the full set is too large."""
    index = np.flatnonzero(mask)
    subsample = None
    if index.size > MAX_CELLS:
        strata = (
            adata.obs.iloc[index]["sample_name"].astype(str) + "|"
            + adata.obs.iloc[index]["broad_label"].astype(str)
        )
        rng = np.random.default_rng(SEED)
        keep = (
            pd.Series(index)
            .groupby(strata.to_numpy(), observed=True)
            .apply(lambda block: block.sample(
                n=max(1, int(round(len(block) * MAX_CELLS / index.size))),
                random_state=SEED))
        )
        index = np.sort(np.asarray(keep).ravel())
        subsample = index.size
        print(f"  {tag}: subsampled {subsample:,} of {mask.sum():,} "
              f"(stratified by sample x label, seed {SEED})")

    block = adata[index].copy()
    bm = Benchmarker(
        block, batch_key="cohort", label_key="broad_label",
        embedding_obsm_keys=ARM_KEYS, n_jobs=-1,
    )
    bm.prepare()
    bm.benchmark()
    return bm.get_results(min_max_scale=False), subsample


start = time.time()
immune_scores, immune_subsample = score(immune, "immune")
print(f"immune scIB done in {time.time() - start:.0f}s")
immune_scores.to_csv(OUT / "scib_immune.csv")
display(immune_scores.round(3))

# %% [markdown]
# ## Secondary reference only — global scIB
#
# Computed for completeness and **never used for selection.** On this dataset the
# global score rewards whichever method most aggressively merges the censored plasma
# populations, which is the failure mode this design exists to avoid.

# %%
global_scores, global_subsample = score(
    np.ones(adata.n_obs, dtype=bool), "global"
)
global_scores.to_csv(OUT / "scib_global.csv")
display(global_scores.round(3))

# %% [markdown]
# ## Dataset-specific diagnostics
#
# Generic scIB does not know about this deposit's censoring, so two diagnostics are
# added that speak to it directly.
#
# ### Depth association
#
# `R²(log1p(total_counts) ~ latent)`. **The statistic was fixed before running**, and
# the choice is not arbitrary: R² depends only on the *column span* of the embedding, so
# it is **rotation-invariant**. Latent axes are arbitrary and differ between methods, so
# a per-dimension statistic like `max |Spearman|` would rank methods on an accident of
# their parameterisation instead of on how much depth information they carry.
#
# Higher is worse here, because depth is confounded with cohort by the stage-04
# censoring.

# %%
depth = adata.obs["total_counts"].to_numpy()
rows = []
for arm in benchmark.ARMS:
    embedding = adata.obsm[f"X_{arm.name}"]
    row = {
        "arm": arm.name,
        "immune": benchmark.depth_association(embedding[immune], depth[immune]),
        "plasma": benchmark.depth_association(embedding[plasma], depth[plasma]),
    }
    for cls in benchmark.IMMUNE_CLASSES:
        mask = (adata.obs["broad_label"] == cls).to_numpy()
        row[cls] = (benchmark.depth_association(embedding[mask], depth[mask])
                    if mask.sum() > 50 else np.nan)
    rows.append(row)
depth_table = pd.DataFrame(rows).set_index("arm")
depth_table.to_csv(OUT / "depth_association.csv")
display(depth_table.round(3))

# %% [markdown]
# ### Compartment mixing
#
# Cohort-mixing entropy per cluster, split immune vs plasma. Each arm gets its own
# Leiden run on its own embedding, with identical parameters and seed.
#
# **Plasma mixing is reported, never optimized.** A method that merges the three plasma
# islands is flagged for inspection, not credited.

# %%
mixing_rows = []
for arm in benchmark.ARMS:
    shell = anndata.AnnData(
        X=np.zeros((adata.n_obs, 0), dtype=np.float32),
        obs=adata.obs[["cohort", "compartment", "sample_name", "patient_id"]].copy(),
        obsm={"X_emb": adata.obsm[f"X_{arm.name}"].copy()},
    )
    sc.pp.neighbors(shell, n_neighbors=15, use_rep="X_emb", random_state=SEED)
    sc.tl.leiden(shell, resolution=1.0, key_added="cl", flavor="igraph",
                 n_iterations=2, directed=False, random_state=SEED)
    table = benchmark.compartment_mixing(shell, cluster_key="cl")
    for compartment in ("immune", "plasma"):
        block = table[table["compartment"] == compartment]
        mixing_rows.append({
            "arm": arm.name, "compartment": compartment,
            "n_clusters": len(block),
            "median_entropy": block["entropy"].median() if len(block) else np.nan,
        })
    print(f"  {arm.name:20s} {shell.obs['cl'].nunique():3d} clusters", flush=True)

mixing = pd.DataFrame(mixing_rows).pivot(
    index="arm", columns="compartment", values="median_entropy")
mixing.to_csv(OUT / "compartment_mixing.csv")
display(mixing.round(3))

# %% [markdown]
# ## The verdict — computed, not narrated
#
# The rule and its tolerances live in `benchmark.DECISION_TOLERANCES` and were fixed
# before any arm ran, exactly as stage 06 pre-declares its F1 thresholds and stage 10
# pre-registers the γ-secretase hypothesis. `decision.csv` is the source of truth and
# `decision.md` is **rendered from it**, so the write-up cannot drift from the rule.
#
# **Harmony is allowed to win.** "The incumbent survived a comparison it could have
# lost" is a valid, reportable outcome — including against a higher conventional global
# scIB score.

# %%
# `get_results()` appends a `Metric Type` row whose cells are the string
# "Aggregate score", so the numeric columns must be filtered before any float cast.
results = immune_scores.drop(index="Metric Type", errors="ignore")
scores = pd.DataFrame({
    "batch_score": results["Batch correction"].astype(float),
    "bio_score": results["Bio conservation"].astype(float),
})
scores.index = [str(i).removeprefix("X_") for i in scores.index]
scores = scores.loc[[arm.name for arm in benchmark.ARMS]]
scores["depth_r2"] = depth_table["immune"]
scores["plasma_mixing"] = mixing["plasma"]

decision = benchmark.decide(scores)
decision.to_csv(OUT / "decision.csv")
display(decision.round(3))

text = benchmark.render_decision(decision, subsample=immune_subsample)
(OUT / "decision.md").write_text(text)
print(text)

# %% [markdown]
# ## Figures

# %%
fig, axes = plt.subplots(1, 3, figsize=(17, 4.5))
axes[0].scatter(scores["batch_score"], scores["bio_score"], s=60)
for name, row in scores.iterrows():
    axes[0].annotate(name, (row["batch_score"], row["bio_score"]), fontsize=7,
                     xytext=(4, 4), textcoords="offset points")
axes[0].set(xlabel="immune batch correction", ylabel="immune bio conservation",
            title="The trade-off that matters")
mixing.plot.bar(ax=axes[1], rot=30)
axes[1].set(ylabel="median cohort entropy", xlabel="",
            title="Mixing: immune scored, plasma diagnosed")
depth_table[["immune", "plasma"]].plot.bar(ax=axes[2], rot=30)
axes[2].set(ylabel="R²(depth ~ latent)", xlabel="",
            title="Depth association (lower is better)")
for ax in axes[1:]:
    ax.tick_params(axis="x", labelsize=7)
fig.tight_layout()
fig.savefig(FIGURES / "benchmark_summary.png", dpi=150)

# %% [markdown]
# ## The incumbent must be untouched
#
# This notebook is read-only over stage 05. If this assertion ever fails, the benchmark
# has mutated the thing it was supposed to be evaluating.

# %%
after = (INTEGRATED.stat().st_size, INTEGRATED.stat().st_mtime_ns)
assert after == INCUMBENT_FINGERPRINT, (
    f"integrated.h5ad CHANGED during the benchmark: {INCUMBENT_FINGERPRINT} -> {after}"
)
print("integrated.h5ad is byte-identical — the incumbent was not mutated.")

candidates = OUT / "candidate_embeddings.npz"
np.savez_compressed(candidates, **{arm.name: adata.obsm[f"X_{arm.name}"]
                                   for arm in benchmark.ARMS})
print(f"candidate embeddings written to {candidates}")
print("\nIf an arm won, stage 05 is re-run with it and every other choice held fixed.")
print("If Harmony survived, nothing changes and that is recorded as the result.")
