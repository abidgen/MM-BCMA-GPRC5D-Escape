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
# # Stage 10 — is the DN population a coherent transcriptional state?
#
# > Among malignant / clone-compatible plasma cells with observed BCMA/GPRC5D
# > double-negativity, is there evidence of a coherent escape-associated transcriptional
# > state **beyond the two antigen genes themselves**?
#
# Restated as the test actually run: **do DN cells show reproducible, non-antigen
# transcriptional organization that cannot be explained by depth or sample structure
# alone?**
#
# **State ≠ subclone.** Transcriptional coherence licenses "escape-associated
# transcriptional state" and nothing stronger. "Subclone" needs independent genetic
# support, and Stage 07's CNV inference failed its donor negative control — that failed
# method is **not** silently reused here, so every patient carries
# `CNV_SUBCLONE_NOT_EVALUABLE`. An underpowered null is not a negative.
#
# Design frozen in `results/10_dn_coherence/stage10_design.md`.

# %%
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from pathlib import Path

from mm_escape import subclone as SC

OUT = Path("results/10_dn_coherence")
pd.set_option("display.width", 280)
SEED = 20260825

# %% [markdown]
# ## The antigen columns are dropped before anything else happens
#
# Selecting DN cells by antigen status is legitimate and necessary. Letting antigen
# expression into the feature matrix that then "discovers" DN structure is not.
#
# Order matters and is the whole point: normalising first would fold antigen counts into
# every cell's size factor, so perturbing an antigen would shift every other gene's
# normalised value and the "independent" result would not be independent.
# `tests/test_dn_coherence.py` perturbs both genes to zero and to 10,000, singly and
# jointly, and requires **bit-identical** PCs, neighbour graphs and coherence calls.

# %%
A = sc.read_h5ad("results/08_dual_antigen_escape/antigen_states.h5ad")
A.obs["patient_id"] = A.obs["patient_id"].astype(str)
A.obs["sample_name"] = A.obs["sample_name"].astype(str)

genes = list(A.var_names)
counts = sp.csr_matrix(A.layers["counts"])
counts, genes = SC.drop_antigen_features(counts, genes)
obs = A.obs.copy()
del A

assert not (set(genes) & set(SC.ANTIGEN_FEATURES))
print(f"{counts.shape[0]:,} cells x {counts.shape[1]:,} features (antigens dropped)")

is_dn = (obs.observed_state == "double_negative").values
depth = obs.depth_ex_antigen.values.astype(float)

# %% [markdown]
# ## Evaluability, and why provisional tiers are not inputs
#
# ≥20 DN cells, ≥20 antigen-positive comparison cells, n_primary ≥ 100 — **all three are
# pre-existing frozen Stage-08 constants**, none invented here. Comparison is within
# patient only: DN vs antigen-positive clone cells.
#
# Stage 10 runs across **every** evaluable patient, never only the four provisional
# `robust-high` ones — conditioning the biological validation on the outcome it exists to
# test would be circular. No function here branches on a tier label, and a regression test
# asserts that changing tier labels cannot change any coherence output.

# %%
E = pd.read_csv(OUT / "stage10_evaluability.csv")
print(f"evaluable rows {int(E.evaluable.sum())} of {len(E)}")
print(E[~E.evaluable][["patient", "denominator", "n_dn", "n_antigen_pos", "reason"]].to_string(index=False))

# %% [markdown]
# ## Each patient gets a fresh, un-integrated local embedding
#
# Local CP10K + `log1p` → local HVGs → local PCA → local kNN → local Leiden, per patient
# per denominator. **Harmony, the global integrated PCA, the global neighbour graph and
# Stage-05 Leiden structure are never used as biological evidence here.** Because every
# analysis is within patient, patient identity is not an axis at all.

# %% [markdown]
# ## Depth was the dominant finding
#
# Shallow cells are both more likely to be DN *and* more likely to sit together in
# low-dimensional space, so an unconditioned test sees depth structure and calls it
# biology — an artifact pointing in exactly the direction this project hopes to find.

# %%
F = pd.read_csv(OUT / "dn_coherence_final_states.csv")
caught = F[(F.perm_p_unconditioned_primary < 0.05) & (F.perm_p_primary >= 0.05)]
print(caught[["patient", "cohort", "perm_p_unconditioned_primary", "perm_p_primary",
              "morans_i_primary", "depth_ratio_dn_over_pos"]].round(4).to_string(index=False))

# %% [markdown]
# Four MMRF patients are significant under the unconditioned null and **not** under the
# depth-stratified one. Their DN cells really do cluster together — because they are the
# shallow cells. MMRF DN cells run **7–40× shallower** than antigen-positive cells
# (`MMRF_1640`: Moran's I 0.47, unconditioned p = 0.001, stratified p = 0.499).
#
# Without this null, Stage 10 would have reported an escape-associated state across most
# of MMRF.

# %% [markdown]
# ## Repeated samples rejected two otherwise-passing patients

# %%
RS = pd.read_csv(OUT / "repeated_sample_dn_coherence.csv", dtype={"patient": str})
print(RS[RS.patient.isin(["27522", "83942"])]
      [["patient", "sample", "n_cells", "n_dn", "perm_p_depth_stratified",
        "sample_evaluable", "repeated_sample_status"]].round(4).to_string(index=False))

# %% [markdown]
# `27522` and `83942` both reached p ≤ 0.002 under **both** denominators, and both are
# **discordant** across their own samples — the pooled signal in each came from one
# sample. `83942` is the sharpest: one patient, one clone, two protocols, and the
# coherence does not replicate. Both are `NOT_SUPPORTED`.

# %% [markdown]
# ## The cross-patient program criterion is not met
#
# Among the four supported patients, 1,718 of 8,797 genes share DE direction — against
# **1,388 expected by chance** (a gene seen in *k* patients is unanimous with probability
# 2^(1−k), and k is 3–4). Observed/expected = **1.24**. That is not a reproducible escape
# program, and it is reported as a negative rather than dressed up as a gene list.
# Coherence in these four rests on **local neighbourhood organization only**.

# %%
print(F.dn_coherence_state.value_counts().to_string())
print(F[["patient", "cohort", "n_dn_primary", "perm_p_primary", "perm_p_sensitivity",
         "depth_ratio_dn_over_pos", "repeated_sample_status",
         "dn_coherence_state"]].round(4).to_string(index=False))

# %% [markdown]
# ## Provisional tier overlay — added only AFTER the states were frozen

# %%
prov = pd.read_csv("results/08_dual_antigen_escape/risk_tier_provisional/risk_tiers_provisional.csv")
J = F.merge(prov[["patient", "final_tier"]], on="patient")
print(pd.crosstab(J.final_tier, J.dn_coherence_state).to_string())

# %% [markdown]
# **The two sets are completely disjoint.** Not one of the four provisional `robust-high`
# patients shows DN transcriptional coherence, and all four coherence-supported patients
# are provisionally `uncertain`. `MMRF_1267` is the instructive case: strongly enriched
# under the unconditioned null, flatly not under the depth-stratified one.
#
# **This is the most consequential result of the stage.** A measurement-robust elevated DN
# fraction and a coherent DN transcriptional state are, in this cohort, *not the same
# patients* — which is exactly why Stage 10 had to run across all evaluable patients
# rather than only the provisional highs, and why the tiers were relabelled provisional
# before it ran.

# %%
print((OUT / "stage10_dn_coherence_summary.md").read_text())

# %% [markdown]
# ---
# # Level-2 arm — phenotype and programs
#
# **Terminology correction first.** The four Level-1 positives were relabelled
# `DN_STRUCTURE_SUPPORTED`. They rest on non-random local organization only, which
# licenses *"non-random DN organization"* and **not** "escape-associated transcriptional
# state". Every Level-1 number above is preserved unchanged.
#
# Criterion frozen in the `stage10_design.md` Level-2 addendum before any Level-2 number
# existed. Programs are `config.LEVEL2_PROGRAMS`, fixed; the γ-secretase axis is
# pre-registered at exactly five genes and is tested regardless of how the rest turns out.

# %%
from mm_escape import config

print("frozen Level-2 programs:", config.LEVEL2_PROGRAMS)
print("pre-registered gamma-secretase:", config.STATE_PROGRAMS["gamma_secretase"])
print(pd.read_csv(OUT / "program_score_vs_depth.csv").to_string(index=False))

# %% [markdown]
# **The depth screen runs first, by requirement.** Broad activity/ETC-type sets track
# library size as a technical property, so MYC and OXPHOS are reported against depth
# *before* any DN-vs-comparator difference. OXPHOS is the strongest tracker (ρ = 0.274);
# MYC is essentially depth-independent (ρ = 0.076).
#
# Depth matching moves the DN/comparator depth ratio from **0.470 to 0.992** — and it
# changes the answers: `gamma_secretase` looks strongly negative unmatched and is not
# matched, and `myc` **flips sign**. Third time in this project that depth conditioning has
# overturned an apparent result.

# %%
R = pd.read_csv(OUT / "level2_program_cohort_tests.csv")
print(R[R.effect == "matched"].round(4).to_string(index=False))

# %% [markdown]
# Cohort-reproducible under both denominators: **`antigen_presentation`, `oxphos`,
# `interferon` — all DN-lower.** `myc`, `stress`, `upr` and `gamma_secretase` are not.
#
# ## The pre-registered γ-secretase hypothesis is not supported
#
# Depth-matched BH 0.387 (primary) / 0.070 (sensitivity), and the direction is **opposite**
# to the registered prediction of a γ-secretase-high escape phenotype. A clean negative for
# a hypothesis written down before looking.

# %%
print(pd.read_csv(OUT / "gamma_secretase_hypothesis.csv").round(4).to_string(index=False))

# %% [markdown]
# ## Pseudobulk DE — and the caveat that matters most
#
# Paired over patients, raw counts summed over depth-matched cells, one observation per
# patient per group. Per-cell DE was **not** used as a substitute. 190 genes reach
# BH < 0.05 under both denominators, and they are overwhelmingly **ER / secretory-pathway
# and plasma-cell identity genes, down in DN cells** — `SPCS1`, `SPCS2`, `SEC61B`,
# `UBE2J1`, `TMBIM6`, `MZB1`, `B2M`.
#
# BCMA and GPRC5D are themselves secretory-pathway-dependent surface proteins of a
# differentiated plasma cell. So part of "double-negative" may mark a **less
# secretory-differentiated plasma-cell state** rather than antigen-specific escape. That is
# biology rather than artifact, but it is not the biology the escape framing assumes, and
# this data cannot separate them.

# %%
DE = pd.read_csv(OUT / "pseudobulk_de_results.csv")
both = (set(DE[(DE.denominator == "primary") & (DE.p_BH < 0.05)].gene)
        & set(DE[(DE.denominator == "sensitivity") & (DE.p_BH < 0.05)].gene))
print(f"significant under both denominators: {len(both)}")
print(DE[(DE.denominator == "primary") & (DE.gene.isin(both))]
      .nsmallest(15, "p")[["gene", "median_log2FC", "frac_up", "p_BH"]].round(4).to_string(index=False))

# %% [markdown]
# ## LIMITATION — the per-patient Level-2 rule turned out near-vacuous
#
# The predeclared rule supports a patient if **any** cohort-reproducible program runs the
# same direction under both denominators. With three reproducible programs, sign agreement
# is nearly free: **26 of 27 evaluable patients qualify (96%)**.
#
# The rule is **not** retuned after the fact — that is the post-hoc tuning this project
# forbids. But the honest reading is that per-patient Level 2 carries almost no
# discriminative information, and **the real Level-2 result is the cohort-level one.**

# %%
L = pd.read_csv(OUT / "stage10_evidence_levels.csv", dtype={"patient": str})
print(pd.crosstab(L.level1_structure, L.level2_state).to_string())
print("\nlicensed language:")
print(L.licensed_language.value_counts().to_string())

# %% [markdown]
# **Level 3 is `CNV_SUBCLONE_NOT_EVALUABLE` for all 32 patients. No patient is called a
# subclone.** The Level-1 and Level-2 convergent sets are *different* patients
# (`MMRF_1720`/`MMRF_2038`/`MMY34339`/`MMY80649` vs the four measurement-highs), and no
# patient is supported on all three of measurement-high, Level 1 and Level 2. The axes stay
# separate and are never combined into a scalar score.
