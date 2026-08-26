# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#   kernelspec:
#     display_name: mm-communication
#     language: python
#     name: mm-communication
# ---

# %% [markdown]
# # 11b — LIANA verification arm (exploratory, post hoc, non-tier-changing)
#
# > **The LIANA arm asks whether an established ligand–receptor framework recovers
# > communication patterns consistent with the frozen targeted Stage-11 analysis, and
# > whether apparent signals survive receiver-state, abundance, depth and cohort controls.**
#
# Written before any LIANA result was inspected.
#
# **What this arm is:** exploratory; **post hoc** relative to the frozen Stage-11 custom
# result; **non-tier-changing**; **non-classifying**.
#
# **What it is not:** a replacement for the frozen Stage-11 analysis, a re-opening of it, or
# a source of any new patient state. The frozen Stage-11 result stands unchanged — its
# composition analysis, its targeted 17-pair panel, the clone-primary vs all-plasma receiver
# amendment, the receiver-side confound finding, and its conclusion that no robust
# immune-evasion signal survives.
#
# **Broad screening does not create confirmatory biology.** A full-interactome run is used
# for ranking and context only. No pathway is selected first and given a hypothesis
# afterwards.

# %%
import json
import platform
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import liana

from mm_escape import communication as CM

STAGE11 = Path("results/11_immune_context")
OUT = STAGE11 / "liana_verification"
OUT.mkdir(parents=True, exist_ok=True)
INTEGRATED = "results/05_integration/integrated.h5ad"
pd.set_option("display.width", 320)

# %% [markdown]
# ## Environment — recorded, not assumed
#
# The installed LIANA API is inspected and used as it actually is; no tutorial is presumed
# valid. Nothing is downgraded and no unrelated package is touched.

# %%
import importlib.metadata as md

env = {"python": sys.version.split()[0], "platform": platform.platform()}
for pkg in ("liana", "scanpy", "anndata", "numpy", "pandas", "scipy", "omnipath",
            "decoupler", "plotnine", "mudata", "tqdm"):
    try:
        env[pkg] = md.version(pkg)
    except Exception:
        env[pkg] = "ABSENT"
(OUT / "liana_environment.txt").write_text("\n".join(f"{k}={v}" for k, v in env.items()))
print("\n".join(f"{k:12s} {v}" for k, v in env.items()))

# %% [markdown]
# ## Method configuration — frozen before any result is read
#
# **Consensus, not a single hand-picked method.** `liana.mt.rank_aggregate` aggregates five
# methods by robust rank aggregation; `liana.mt.cellchat` is run **separately and
# additionally**, for continuity with the original Stage-11 plan, which named the CellChat
# algorithm. **Both are preserved and reported whatever they show** — no method is run and
# then discarded for showing the wrong thing.
#
# **Resource is held fixed at `consensus` for both runs**, so method and resource are not
# confounded. LIANA's resource list is not searched for one that gives an appealing answer.
#
# **Score orientation matters and is easy to get backwards:** `magnitude_rank` and
# `specificity_rank` are RRA p-like quantities where **lower is stronger**. Associations are
# therefore run on `-log10(magnitude_rank)`, so a positive coefficient means *stronger
# inferred communication with higher DN burden*.

# %%
CONFIG = {
    "liana_version": env["liana"],
    "primary_method": "rank_aggregate (RRA consensus)",
    "consensus_methods": [m.method_name for m in liana.mt.rank_aggregate.methods],
    "continuity_method": "cellchat",
    "resource": "consensus",
    "resource_interactions": int(liana.rs.select_resource("consensus").shape[0]),
    "expr_prop": 0.1,
    "min_cells": CM.MIN_SENDER_CELLS,
    "n_perms": 1000,
    "seed": 1337,
    "use_raw": False,
    "expression": "CP10K + log1p over the frozen 32,991-gene intersected space",
    "per_patient": "liana Method.by_sample(sample_key='patient_id') — LIANA's own documented "
                   "per-sample API; the method is run independently within each patient",
    "score_columns": {
        "rank_aggregate": {"magnitude": "magnitude_rank", "specificity": "specificity_rank",
                           "orientation": "LOWER magnitude_rank = stronger"},
        "cellchat": {"magnitude": "lr_probs", "specificity": "cellchat_pvals",
                     "orientation": "HIGHER lr_probs = stronger"}},
    "score_for_association": "consensus: -log10(magnitude_rank); cellchat: lr_probs. "
                             "Both oriented so HIGHER = stronger inferred communication, so a "
                             "positive coefficient means stronger communication with higher DN.",
    "primary_receiver": "clone_primary (frozen stage-08 primary denominator)",
    "sensitivity_receiver": "all_plasma (all stage-06 PlasmaCell)",
    "senders": list(CM.LR_SENDERS),
    "min_patients_for_association": 20,
    "confound_model": list(CM.CONFOUNDERS),
    "multiple_testing": "Benjamini-Hochberg over the full tested interaction space",
}
(OUT / "liana_method_config.json").write_text(json.dumps(CONFIG, indent=2))
print(json.dumps(CONFIG, indent=2))

# %% [markdown]
# `min_cells` is the frozen Stage-11 `MIN_SENDER_CELLS` (20), reused rather than re-derived.
# `expr_prop`, `n_perms` and `seed` are LIANA's own defaults, taken as-is so they are not a
# tuned surface.

# %% [markdown]
# ## Populations — the accepted Stage-11 categories, unchanged
#
# Senders are the frozen Stage-11 immune categories. **`cytotoxic_mixed` is never folded
# into T or NK** — the 2,207 cluster-23 cells whose lineage the evidence does not resolve
# stay their own category, exactly as Stage 11 froze them. **No population is redefined on
# the basis of LIANA output.**

# %%
labels = pd.read_csv("results/06_annotation/per_cell_labels.csv.gz",
                     dtype={"patient_id": str, "sample_name": str})
revised = pd.read_csv(
    "results/06_annotation/cluster23_local/trbc_context_revision/revised_lineage_calls.csv.gz",
    usecols=["cell_id", "call"])
C23 = {"NK": "NK_core", "T_NK_mixed": "cytotoxic_mixed", "unresolved": "cytotoxic_mixed",
       "T_ab": "Tcell", "T_gd": "Tcell"}
labels = labels.merge(revised, on="cell_id", how="left")
labels["lineage"] = np.where(labels.call.notna(), labels.call.map(C23), labels.cell_type)

frozen_cells = pd.read_csv("results/08_dual_antigen_escape/cell_antigen_states.csv.gz",
                           dtype={"patient_id": str})
clone_primary = set(frozen_cells.loc[frozen_cells.in_primary, "cell_id"])
patients8 = pd.read_csv("results/08_dual_antigen_escape/patient_antigen_states_primary.csv",
                        dtype={"patient": str})
COHORT_PATIENTS = set(patients8.patient)

labels["role"] = np.where(
    labels.cell_id.isin(clone_primary), "clone_primary",
    np.where(labels.lineage == "PlasmaCell", "other_plasma",
             np.where(labels.lineage.isin(CM.LR_SENDERS), labels.lineage, "drop")))
cells = labels[(labels.role != "drop") & labels.patient_id.isin(COHORT_PATIENTS)].copy()
print(cells.role.value_counts().to_string())
assert int((cells.role == "clone_primary").sum()) == 21906, "frozen primary denominator"

# %% [markdown]
# ## Expression — CP10K + log1p over the frozen intersected gene space
#
# LIANA's magnitude scores are means over normalised expression; this is the standard input
# and is **not** imputation, denoising or smoothing. The project's ban on those applies to
# **positivity calls**, which this arm does not make — every antigen-negativity call in this
# project remains the frozen Stage-08 raw-count one.
#
# Library size is the row sum over the **intersected 32,991-gene space**, matching Stage 08
# and Stage 11, not the QC-time `obs["total_counts"]`.

# %%
resource = liana.rs.select_resource(CONFIG["resource"])
resource_genes = sorted({g for col in ("ligand", "receptor")
                         for entry in resource[col] for g in str(entry).split("_")})
raw = CM.stream_gene_counts(INTEGRATED, [g for g in resource_genes
                                         if g in set(pd.read_csv(
                                             "results/08_dual_antigen_escape/multi_antigen_coverage/"
                                             "primary_denominator_gene_means.csv",
                                             index_col=0).index)])
present = [c for c in raw.columns if c != "total_counts"]
print(f"resource {resource.shape[0]} interactions, {len(resource_genes)} genes; "
      f"{len(present)} present in the intersected gene space")

meta = cells.set_index("cell_id")
raw = raw.loc[meta.index]
counts = sp.csr_matrix(raw[present].values)
norm = sp.csr_matrix(np.log1p(raw[present].values / np.maximum(raw.total_counts.values[:, None], 1) * 1e4))
adata = ad.AnnData(X=norm, obs=meta[["patient_id", "cohort", "role", "sample_name"]].copy(),
                   var=pd.DataFrame(index=present))
adata.layers["counts"] = counts
adata.obs["role_sensitivity"] = np.where(adata.obs.role.isin(["clone_primary", "other_plasma"]),
                                         "all_plasma", adata.obs.role.astype(str))
print(adata)

# %% [markdown]
# ## Evaluability — declared before running, because one patient cannot be scored at all
#
# A patient can only carry a communication score if it contributes **at least one sender
# category and the receiver category**, each clearing `min_cells`. Below that there is no
# sender→receiver pair to compute, and LIANA's own log2FC step divides by an empty
# comparison group.
#
# `25183` contributes **zero immune cells of any sender category** — in the all-plasma
# construction every one of its cells collapses into a single group. It is
# `LIANA_NOT_EVALUABLE`, reported rather than silently dropped. **This reproduces the frozen
# Stage-11 behaviour exactly**: `25183` is the one patient absent from the frozen custom
# communication table, for the same reason.

# %%
adata.obs["role"] = adata.obs["role"].astype(str)          # drop empty categorical levels
adata.obs["role_sensitivity"] = adata.obs["role_sensitivity"].astype(str)


def evaluable_patients(obj, role_col, receiver):
    counts = obj.obs.groupby(["patient_id", role_col], observed=True).size().unstack(fill_value=0)
    senders = [c for c in counts.columns if c in CM.LR_SENDERS]
    ok_sender = (counts[senders] >= CONFIG["min_cells"]).any(axis=1) if senders else False
    ok_receiver = counts.get(receiver, pd.Series(0, index=counts.index)) >= CONFIG["min_cells"]
    return counts, (ok_sender & ok_receiver)


eval_rows = []
eligible = {}
for key, obj, col, recv in [
        ("consensus_clone_primary", adata, "role", "clone_primary"),
        ("cellchat_clone_primary", adata, "role", "clone_primary"),
        ("consensus_all_plasma", adata, "role_sensitivity", "all_plasma")]:
    counts, ok = evaluable_patients(obj, col, recv)
    eligible[key] = set(counts.index[ok])
    for pid in counts.index:
        eval_rows.append({"run": key, "patient": pid, "receiver": recv,
                          "n_receiver_cells": int(counts.get(recv, pd.Series(0, index=counts.index))[pid]),
                          "n_sender_categories_ge_min": int(
                              sum(counts.loc[pid, c] >= CONFIG["min_cells"]
                                  for c in counts.columns if c in CM.LR_SENDERS)),
                          "status": "EVALUABLE" if ok[pid] else "LIANA_NOT_EVALUABLE"})
evaluability = pd.DataFrame(eval_rows)
evaluability.to_csv(OUT / "liana_patient_evaluability.csv", index=False)
print(evaluability.groupby(["run", "status"]).size().to_string())
print("\nnot evaluable:")
print(evaluability[evaluability.status == "LIANA_NOT_EVALUABLE"]
      [["run", "patient", "n_receiver_cells", "n_sender_categories_ge_min"]].to_string(index=False))

# %% [markdown]
# ## Per-patient LIANA — `by_sample`, not one pooled run
#
# **A single pooled all-patient LIANA run would produce global scores that are not
# patient-level observations**, and treating them as such would be the pseudoreplication the
# frozen Stage-11 design exists to avoid. `Method.by_sample` is LIANA's own documented API
# for this: the method is run **independently within each patient**.

# %%
def run_liana(method, obj, groupby, key):
    keep = obj.obs.patient_id.isin(eligible[key]).values
    sub = obj[keep].copy()
    sub.obs[groupby] = sub.obs[groupby].astype(str)
    df = method.by_sample(sub, sample_key="patient_id", groupby=groupby,
                          resource_name=CONFIG["resource"], expr_prop=CONFIG["expr_prop"],
                          min_cells=CONFIG["min_cells"], use_raw=False,
                          n_perms=CONFIG["n_perms"], seed=CONFIG["seed"],
                          verbose=False, inplace=False)
    df = df.rename(columns={"patient_id": "patient"})
    df["run"] = key
    return df


sens = adata.copy()
sens.obs["role"] = sens.obs["role_sensitivity"].astype(str)

runs = {}
runs["consensus_clone_primary"] = run_liana(liana.mt.rank_aggregate, adata, "role",
                                            "consensus_clone_primary")
runs["cellchat_clone_primary"] = run_liana(liana.mt.cellchat, adata, "role",
                                           "cellchat_clone_primary")
runs["consensus_all_plasma"] = run_liana(liana.mt.rank_aggregate, sens, "role",
                                         "consensus_all_plasma")
for k, v in runs.items():
    print(f"{k:28s} rows {len(v):8,d}  patients {v.patient.nunique():2d}")

# %% [markdown]
# ## Restrict to the declared sender → receiver directions
#
# Senders are the four frozen immune categories; the receiver is the plasma population.
# Plasma→immune and immune→immune directions are computed by LIANA but are not this arm's
# question and are dropped **before** anything is ranked.

# %%
RECEIVER = {"consensus_clone_primary": "clone_primary", "cellchat_clone_primary": "clone_primary",
            "consensus_all_plasma": "all_plasma"}
#: The two methods report different quantities in opposite directions. Both are converted to
#: "higher = stronger" so a positive DN coefficient always means the same thing.
SCORE = {"consensus_clone_primary": ("magnitude_rank", "specificity_rank", "neglog10"),
         "cellchat_clone_primary": ("lr_probs", "cellchat_pvals", "identity"),
         "consensus_all_plasma": ("magnitude_rank", "specificity_rank", "neglog10")}

kept, native = [], {}
for key, df in runs.items():
    d = df[df.source.isin(CM.LR_SENDERS) & (df.target == RECEIVER[key])].copy()
    d["interaction"] = d.ligand_complex + "->" + d.receptor_complex
    d["sender_receiver"] = d.source.astype(str) + "|" + d.interaction
    native[key] = d
    mag, spec, orient = SCORE[key]
    core = d[["run", "patient", "source", "target", "ligand_complex", "receptor_complex",
              "interaction", "sender_receiver"]].copy()
    core["magnitude_raw"] = d[mag].values
    core["specificity_raw"] = d[spec].values
    core["score"] = (-np.log10(np.clip(d[mag].values, 1e-300, None)) if orient == "neglog10"
                     else d[mag].values)
    kept.append(core)
full = pd.concat(kept, ignore_index=True)
# every run's native columns are preserved, not only the normalised core
pd.concat(native.values(), ignore_index=True).to_csv(
    OUT / "liana_full_results.csv.gz", index=False, compression="gzip")
print(full.groupby("run").agg(rows=("interaction", "size"),
                              interactions=("sender_receiver", "nunique"),
                              patients=("patient", "nunique")).to_string())

# %% [markdown]
# ## Cross-check against the frozen 17-pair custom panel
#
# The frozen panel is **read, never redefined**. For every one of its interactions: does
# LIANA's resource contain it, does LIANA evaluate it, at what rank, and does the direction
# against DN agree with the custom analysis?

# %%
CUSTOM_PAIRS = [(l, r) for l, r in CM.LR_CANDIDATES if r != "None"]
assert len(CUSTOM_PAIRS) == 17, "the frozen 17-pair panel must be unchanged"
resource_set = {(str(l), str(r)) for l, r in zip(resource.ligand, resource.receptor)}

custom_frozen = pd.read_csv(STAGE11 / "communication_context_vs_dn.csv")
prim = native["consensus_clone_primary"].copy()
prim = prim.assign(mag_rank_pct=prim.groupby("patient").magnitude_rank.rank(pct=True))

rows = []
for lig, rec in CUSTOM_PAIRS:
    in_resource = (lig, rec) in resource_set
    for sender in CM.LR_SENDERS:
        sel = prim[(prim.source == sender) & (prim.ligand_complex == lig)
                   & (prim.receptor_complex == rec)]
        cf = custom_frozen[(custom_frozen.sender == sender) & (custom_frozen.ligand == lig)
                           & (custom_frozen.receptor == rec)]
        rows.append({
            "sender": sender, "ligand": lig, "receptor": rec,
            "in_liana_resource": in_resource,
            "evaluated_by_liana": len(sel) > 0,
            "n_patients_evaluated": int(sel.patient.nunique()),
            "median_magnitude_rank": float(sel.magnitude_rank.median()) if len(sel) else np.nan,
            "median_rank_percentile": float(sel.mag_rank_pct.median()) if len(sel) else np.nan,
            "median_specificity_rank": float(sel.specificity_rank.median()) if len(sel) else np.nan,
            "custom_coef_vs_dn": float(cf.coef_adjusted.iloc[0]) if len(cf) else np.nan,
            "custom_p_vs_dn": float(cf.p_adjusted.iloc[0]) if len(cf) else np.nan})
crosscheck = pd.DataFrame(rows)
print(f"of {len(CUSTOM_PAIRS)} custom pairs, {crosscheck.groupby(['ligand','receptor']).in_liana_resource.first().sum()} "
      f"are in LIANA's consensus resource")
print(crosscheck[crosscheck.evaluated_by_liana].sort_values("median_rank_percentile")
      .head(12).round(4).to_string(index=False))

# %% [markdown]
# ## Patient-level scores and the DN association
#
# **Patient is the unit.** LIANA computes from cells internally, but every association below
# uses one score per patient per interaction. The confound model is the **frozen Stage-11
# one, unchanged** — cohort, immune depth, immune-cell abundance, number of samples — and
# the unadjusted estimate is reported beside the adjusted one.

# %%
comp = pd.read_csv(STAGE11 / "patient_immune_composition.csv", dtype={"patient": str})
pred = (patients8[["patient", "observed_double_negative_fraction", "excess_dn"]]
        .rename(columns={"observed_double_negative_fraction": "obs_dn_primary"}))
D = comp.merge(pred, on="patient")
CONF = [D.cohort.values, np.log10(D.median_depth_ex_antigen_immune.values),
        np.log10(D.n_immune_cells.values.astype(float)),
        np.log10(D.n_samples.values.astype(float))]

patient_scores = full[["run", "patient", "source", "ligand_complex", "receptor_complex",
                       "interaction", "sender_receiver", "magnitude_raw", "specificity_raw",
                       "score"]].copy()
patient_scores.to_csv(OUT / "liana_patient_level_scores.csv.gz", index=False, compression="gzip")
print(f"{len(patient_scores):,} patient-level interaction scores")


def associate(df, predictor="obs_dn_primary", min_patients=CONFIG["min_patients_for_association"]):
    out = []
    for (sender, inter), g in df.groupby(["source", "sender_receiver"], observed=True):
        m = D.merge(g[["patient", "score"]], on="patient")
        if len(m) < min_patients:
            continue
        conf = [m.cohort.values, np.log10(m.median_depth_ex_antigen_immune.values),
                np.log10(m.n_immune_cells.values.astype(float)),
                np.log10(m.n_samples.values.astype(float))]
        un = CM.ols_association(m.score.values, m[predictor].values, ())
        adj = CM.ols_association(m.score.values, m[predictor].values, conf)
        wc = CM.within_cohort_spearman(m.score.values, m[predictor].values, m.cohort.values)
        signs = [np.sign(v["rho"]) for v in wc.values() if np.isfinite(v.get("rho", np.nan))]
        out.append({"sender": sender, "interaction": inter.split("|", 1)[1], "n": adj["n"],
                    "coef_unadjusted": un["coef"], "p_unadjusted": un["p"],
                    "coef_adjusted": adj["coef"], "ci_lo": adj["ci_lo"], "ci_hi": adj["ci_hi"],
                    "p_adjusted": adj["p"],
                    **{f"rho_{k}": v.get("rho", np.nan) for k, v in wc.items()},
                    "within_cohort_same_sign": bool(len(set(signs)) == 1 and len(signs) == 3)})
    res = pd.DataFrame(out)
    if len(res):
        res["p_adj_BH"] = CM.benjamini_hochberg(res.p_adjusted.values)
        res["p_unadj_BH"] = CM.benjamini_hochberg(res.p_unadjusted.values)
    return res


assoc = {}
for key in runs:
    a = associate(patient_scores[patient_scores.run == key])
    a["run"] = key
    assoc[key] = a
    print(f"{key:28s} tested {len(a):5,d}  raw p<0.05 {int((a.p_adjusted < 0.05).sum()):4d}  "
          f"BH<0.10 {int((a.p_adj_BH < 0.10).sum()):3d}")
associations = pd.concat(assoc.values(), ignore_index=True)
associations.to_csv(OUT / "liana_vs_dn_associations.csv", index=False)

# %% [markdown]
# **The full tested interaction space is written out, not only the significant rows.**

# %% [markdown]
# ## Mandatory receiver-side decomposition
#
# **Every apparent LIANA hit is decomposed into its sender-ligand and receiver-receptor
# halves.** This is not optional and it is not applied selectively: the receptor side of any
# score is measured on the *same plasma cells* whose DN status is the predictor, so an
# association can arise with no communication involved at all. This is the confound the
# frozen Stage-11 analysis already found on its own targeted panel, and it is tested here
# before any LIANA hit is described as signalling.
#
# **Classification rule, fixed now:** if the receptor term alone reproduces most of the
# association (same sign and |coef| at least half the full interaction's, with the
# receptor-only association reaching raw p < 0.05), the interaction is
# `RECEIVER_STATE_CONFOUNDED` and **is not called immune-signalling evidence**.

# %%
def side_expression(obj, role_col, role_value, genes, tag):
    m = obj.obs[role_col] == role_value
    sub = obj[m]
    idx = [obj.var_names.get_loc(g) for g in genes if g in obj.var_names]
    names = [g for g in genes if g in obj.var_names]
    X = np.asarray(sub.layers["counts"][:, idx].todense()) if len(idx) else np.zeros((sub.n_obs, 0))
    tot = np.asarray(sub.layers["counts"].sum(axis=1)).ravel()
    df = pd.DataFrame(X, columns=names)
    df["patient"] = sub.obs.patient_id.astype(str).values
    df["_tot"] = tot
    out = df.groupby("patient").apply(
        lambda g: pd.Series({c: g[c].sum() / max(g._tot.sum(), 1) * 1e6 for c in names}),
        include_groups=False)
    out.columns = [f"{tag}_{c}" for c in out.columns]
    return out


all_ligands = sorted({i.split("->")[0] for i in full.interaction.unique()})
all_receptors = sorted({i.split("->")[1] for i in full.interaction.unique()})
flat_l = sorted({g for c in all_ligands for g in c.split("_")})
flat_r = sorted({g for c in all_receptors for g in c.split("_")})

recv_expr = side_expression(adata, "role", "clone_primary", flat_r, "recv")
recv_expr_allp = side_expression(sens, "role", "all_plasma", flat_r, "recv")
send_expr = {s: side_expression(adata, "role", s, flat_l, "send") for s in CM.LR_SENDERS}
print(f"receiver receptor CPM: {recv_expr.shape}; sender ligand CPM per category: "
      f"{ {k: v.shape[0] for k, v in send_expr.items()} }")


def one_side_association(values_by_patient, colname):
    m = D.merge(values_by_patient[[colname]].reset_index().rename(columns={"index": "patient"}),
                on="patient", how="inner")
    if len(m) < CONFIG["min_patients_for_association"]:
        return None
    conf = [m.cohort.values, np.log10(m.median_depth_ex_antigen_immune.values),
            np.log10(m.n_immune_cells.values.astype(float)),
            np.log10(m.n_samples.values.astype(float))]
    return CM.ols_association(np.log1p(m[colname].values), m.obs_dn_primary.values, conf)


# %% [markdown]
# The decomposition is run on **every interaction reaching raw p < 0.05 adjusted** in the
# primary consensus run, plus the interactions the instructions name explicitly —
# `PDCD1 -> CD274` and the TRAIL axis — whether or not they are significant here, so the
# comparison to the frozen custom result is like-for-like.

# %%
primary_assoc = assoc["consensus_clone_primary"]
NAMED = ["PDCD1->CD274", "TNFSF10->TNFRSF10A", "TNFSF10->TNFRSF10B"]
candidates = primary_assoc[(primary_assoc.p_adjusted < 0.05)
                           | primary_assoc.interaction.isin(NAMED)].copy()
print(f"{len(candidates)} interactions enter receiver-side decomposition "
      f"({int((primary_assoc.p_adjusted < 0.05).sum())} by significance, rest named)")

dec_rows = []
for _, r in candidates.iterrows():
    lig, rec = r.interaction.split("->")
    rec_units = [f"recv_{g}" for g in rec.split("_") if f"recv_{g}" in recv_expr.columns]
    send_units = [f"send_{g}" for g in lig.split("_")
                  if f"send_{g}" in send_expr[r.sender].columns]
    rec_stat = None
    if rec_units:
        tmp = recv_expr[rec_units].mean(axis=1).to_frame("recv_agg")
        rec_stat = one_side_association(tmp, "recv_agg")
    send_stat = None
    if send_units:
        tmp = send_expr[r.sender][send_units].mean(axis=1).to_frame("send_agg")
        send_stat = one_side_association(tmp, "send_agg")
    dec_rows.append({
        "sender": r.sender, "interaction": r.interaction,
        "full_coef": r.coef_adjusted, "full_p": r.p_adjusted, "full_p_BH": r.p_adj_BH,
        "receptor_coef": rec_stat["coef"] if rec_stat else np.nan,
        "receptor_p": rec_stat["p"] if rec_stat else np.nan,
        "ligand_coef": send_stat["coef"] if send_stat else np.nan,
        "ligand_p": send_stat["p"] if send_stat else np.nan})
decomp = pd.DataFrame(dec_rows)


def classify(row):
    if not np.isfinite(row.full_p):
        return "NOT_EVALUABLE"
    rec_ok = (np.isfinite(row.receptor_p) and row.receptor_p < 0.05
              and np.sign(row.receptor_coef) == np.sign(row.full_coef)
              and abs(row.receptor_coef) >= 0.5 * abs(row.full_coef))
    if rec_ok:
        return "RECEIVER_STATE_CONFOUNDED"
    if row.full_p_BH < 0.10:
        return "EXPLORATORY_LIANA_ONLY"
    return "NOT_REPRODUCED_BY_LIANA"


decomp["classification"] = decomp.apply(classify, axis=1)
decomp.to_csv(OUT / "liana_receiver_side_confound.csv", index=False)
print(decomp.sort_values("full_p").round(4).to_string(index=False))

# %% [markdown]
# ## Testing the Stage-10 receiver-state confound directly
#
# Stage 10 froze the finding that DN cells occupy a **less secretory, less differentiated**
# plasma-cell state, lower across ER/secretory machinery, antigen presentation, OXPHOS and
# interferon. If LIANA's receiver side is reading that state rather than communication, the
# receptor pool should track those programs and depth.
#
# **This is falsification and context, not a discovery screen.** The question is: *is LIANA
# detecting communication, or the already-known plasma-cell state?*

# %%
PROGRAMS = {
    "secretory_ER": ["SPCS1", "SPCS2", "SEC61B", "UBE2J1", "TMBIM6", "MZB1"],
    "antigen_presentation": ["B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "TAP1", "TAPBP"],
    "oxphos": ["NDUFA1", "NDUFB2", "COX7C", "COX7B", "ATP5F1E", "UQCRB", "ATP5MC2"],
}
prog_genes = sorted({g for v in PROGRAMS.values() for g in v})
prog_counts = CM.stream_gene_counts(INTEGRATED, prog_genes)
pc = frozen_cells[frozen_cells.in_primary].merge(prog_counts, left_on="cell_id",
                                                 right_index=True)
prog_patient = pd.DataFrame(index=sorted(pc.patient_id.unique()))
for name, genes in PROGRAMS.items():
    have = [g for g in genes if g in pc.columns]
    prog_patient[name] = (pc.groupby("patient_id")
                          .apply(lambda g: g[have].sum().sum() / max(g.total_counts.sum(), 1) * 1e6,
                                 include_groups=False))
prog_patient["depth"] = pc.groupby("patient_id").depth_ex_antigen.median()
prog_patient.index.name = "patient"

recv_pool = recv_expr.mean(axis=1).to_frame("receptor_pool_cpm")
state = prog_patient.join(recv_pool, how="inner")
from scipy.stats import spearmanr

state_rows = []
for col in ["secretory_ER", "antigen_presentation", "oxphos", "depth"]:
    rho, p = spearmanr(state.receptor_pool_cpm, state[col])
    state_rows.append({"receiver_quantity": "mean receptor-pool CPM", "against": col,
                       "spearman_rho": rho, "p": p, "n_patients": len(state)})
state_tbl = pd.DataFrame(state_rows)
print(state_tbl.round(4).to_string(index=False))
state_tbl.to_csv(OUT / "liana_receiver_state_vs_stage10_programs.csv", index=False)

# %% [markdown]
# ## Abundance sensitivity
#
# An LR difference driven by how many sender or receiver cells a patient contributed is not
# a communication change. Sender and receiver abundance are reported for every candidate,
# and the association is re-fitted with both added to the frozen confound model.

# %%
abund = (cells.groupby(["patient_id", "role"]).size().unstack(fill_value=0)
         .rename_axis("patient").reset_index())
abund_rows = []
for _, r in candidates.iterrows():
    m = D.merge(abund, on="patient")
    g = patient_scores[(patient_scores.run == "consensus_clone_primary")
                       & (patient_scores.source == r.sender)
                       & (patient_scores.interaction == r.interaction)][["patient", "score"]]
    m = m.merge(g, on="patient")
    if len(m) < CONFIG["min_patients_for_association"]:
        continue
    base_conf = [m.cohort.values, np.log10(m.median_depth_ex_antigen_immune.values),
                 np.log10(m.n_immune_cells.values.astype(float)),
                 np.log10(m.n_samples.values.astype(float))]
    with_ab = base_conf + [np.log10(m[r.sender].values.astype(float) + 1),
                           np.log10(m["clone_primary"].values.astype(float) + 1)]
    base = CM.ols_association(m.score.values, m.obs_dn_primary.values, base_conf)
    adj = CM.ols_association(m.score.values, m.obs_dn_primary.values, with_ab)
    abund_rows.append({
        "sender": r.sender, "interaction": r.interaction, "n": base["n"],
        "median_sender_cells": float(np.median(m[r.sender])),
        "median_receiver_cells": float(np.median(m["clone_primary"])),
        "coef_frozen_model": base["coef"], "p_frozen_model": base["p"],
        "coef_plus_abundance": adj["coef"], "p_plus_abundance": adj["p"],
        "abundance_sensitive": bool(base["p"] < 0.05 and adj["p"] >= 0.05)})
abundance = pd.DataFrame(abund_rows)
abundance.to_csv(OUT / "liana_abundance_sensitivity.csv", index=False)
print(abundance.round(4).to_string(index=False))

# %% [markdown]
# ## Final classification and the candidate-panel cross-check file
#
# Labels are restrained by design. **`IMMUNE_EVASION_CONFIRMED` is not an available label
# and never will be.**

# %%
label_map = decomp.set_index(["sender", "interaction"]).classification.to_dict()
ab_map = abundance.set_index(["sender", "interaction"]).abundance_sensitive.to_dict()


def final_label(row):
    key = (row.sender, row.interaction)
    lab = label_map.get(key)
    if lab is None:
        return "NOT_EVALUABLE" if not np.isfinite(row.p_adjusted) else "NOT_REPRODUCED_BY_LIANA"
    if lab == "EXPLORATORY_LIANA_ONLY" and ab_map.get(key):
        return "ABUNDANCE_SENSITIVE"
    return lab


primary_assoc = primary_assoc.assign(classification=primary_assoc.apply(final_label, axis=1))
print(primary_assoc.classification.value_counts().to_string())

cc = crosscheck.merge(
    primary_assoc[["sender", "interaction", "coef_adjusted", "p_adjusted", "p_adj_BH",
                   "classification"]]
    .assign(ligand=lambda d: d.interaction.str.split("->").str[0],
            receptor=lambda d: d.interaction.str.split("->").str[1])
    .drop(columns="interaction"),
    on=["sender", "ligand", "receptor"], how="left")
cc["liana_direction_agrees_with_custom"] = np.where(
    cc.coef_adjusted.notna() & cc.custom_coef_vs_dn.notna(),
    np.sign(cc.coef_adjusted) == np.sign(cc.custom_coef_vs_dn), np.nan)
cc["panel_status"] = np.where(~cc.evaluated_by_liana, "NOT_EVALUABLE",
                              cc.classification.fillna("NOT_REPRODUCED_BY_LIANA"))
cc.to_csv(OUT / "liana_candidate_panel_crosscheck.csv", index=False)
print(cc.panel_status.value_counts().to_string())
print(cc[cc.evaluated_by_liana][
    ["sender", "ligand", "receptor", "median_rank_percentile", "coef_adjusted", "p_adjusted",
     "custom_coef_vs_dn", "liana_direction_agrees_with_custom", "panel_status"]]
    .sort_values("median_rank_percentile").head(15).round(4).to_string(index=False))

# %% [markdown]
# ## What LIANA changes about the frozen Stage-11 interpretation
#
# **Nothing.** No tier, structure state, phenotype state, composition conclusion or coverage
# eligibility is modified, and no new patient classifier exists. The frozen Stage-11 custom
# outputs are untouched on disk.

# %%
frozen_now = {
    "stage09b_tiers": "results/08_dual_antigen_escape/risk_tier_provisional/risk_tiers_provisional.csv",
    "stage10_states": "results/10_dn_coherence/dn_coherence_final_states.csv",
    "stage11_custom_lr": "results/11_immune_context/communication_context_vs_dn.csv",
    "stage11_composition": "results/11_immune_context/patient_immune_composition.csv",
    "coverage_qc": "results/08_dual_antigen_escape/multi_antigen_coverage/target_measurement_qc.csv",
}
import hashlib

digests = {k: hashlib.sha256(Path(v).read_bytes()).hexdigest() for k, v in frozen_now.items()
           if Path(v).exists()}
(OUT / "frozen_state_digests.json").write_text(json.dumps(digests, indent=2, sort_keys=True))
print(f"recorded digests for {len(digests)} frozen state files — LIANA modifies none of them")

# %% [markdown]
# ## The circularity the unrestricted resource reintroduces
#
# The frozen Stage-11 targeted panel was built so that `TNFRSF17` could appear **only** as a
# named receptor in a fixed external list, never as a discovered feature — its BAFF pairs
# were `TNFSF13`/`TNFSF13B` → **`TNFRSF13B`**, deliberately not `TNFRSF17`.
#
# LIANA's consensus resource contains `TNFSF13B → TNFRSF17`, and an unrestricted screen has
# no way to know it must not use it. So the check is made explicit here: **any interaction
# whose ligand or receptor is one of the two antigens whose negativity defines the DN
# predictor is structurally circular**, not merely receiver-state confounded. `obs_dn_primary`
# is the fraction of receiver cells negative for `TNFRSF17` **and** `GPRC5D`, so higher DN
# forces lower `TNFRSF17` in those same cells. A negative coefficient is arithmetic.

# %%
DN_DEFINING = {"TNFRSF17", "GPRC5D"}


def antigen_circular(interaction):
    lig, rec = interaction.split("->")
    units = set(lig.split("_")) | set(rec.split("_"))
    return bool(units & DN_DEFINING)


for key, a in assoc.items():
    a["antigen_circular"] = a.interaction.map(antigen_circular)
associations = pd.concat(assoc.values(), ignore_index=True)
associations.to_csv(OUT / "liana_vs_dn_associations.csv", index=False)

circ = associations[associations.antigen_circular]
print(f"{circ.sender_receiver.nunique() if 'sender_receiver' in circ else circ.interaction.nunique()} "
      f"antigen-containing interactions across runs")
print(circ[["run", "sender", "interaction", "n", "coef_adjusted", "p_adjusted", "p_adj_BH"]]
      .sort_values("p_adjusted").round(4).to_string(index=False))
print("\nBH<0.10 hits per run, and how many of them are antigen-circular:")
print(associations[associations.p_adj_BH < 0.10]
      .groupby("run").agg(bh_hits=("interaction", "size"),
                          antigen_circular=("antigen_circular", "sum")).to_string())

# %%
decomp["antigen_circular"] = decomp.interaction.map(antigen_circular)
decomp.loc[decomp.antigen_circular, "classification"] = "RECEIVER_STATE_CONFOUNDED"
decomp.to_csv(OUT / "liana_receiver_side_confound.csv", index=False)
print(decomp.classification.value_counts().to_string())

# %% [markdown]
# ## Summary of what LIANA found
#
# Counts printed here are the ones carried into `liana_summary.md`.

# %%
prim_assoc = associations[associations.run == "consensus_clone_primary"]
print(f"interactions scored per patient (clone-primary): "
      f"{full[full.run == 'consensus_clone_primary'].sender_receiver.nunique():,}")
print(f"evaluable patients: {full.patient.nunique()} of 32 "
      f"(25183 LIANA_NOT_EVALUABLE — zero sender cells)")
print(f"tested at >= {CONFIG['min_patients_for_association']} patients: {len(prim_assoc)}")
print(f"raw p<0.05 adjusted: {int((prim_assoc.p_adjusted < 0.05).sum())}  "
      f"BH<0.10: {int((prim_assoc.p_adj_BH < 0.10).sum())}")
print(f"receiver-state confounded: {int((decomp.classification == 'RECEIVER_STATE_CONFOUNDED').sum())}")
print(f"abundance-sensitive: {int(abundance.abundance_sensitive.sum())}")
cc_named = cc[(cc.ligand == 'PDCD1') & (cc.receptor == 'CD274')]
print(f"\nPDCD1->CD274 in LIANA's consensus resource: {bool(cc_named.in_liana_resource.iloc[0])}; "
      f"evaluated: {bool(cc_named.evaluated_by_liana.any())}")
print(f"custom pairs in LIANA resource: "
      f"{int(cc.groupby(['ligand','receptor']).in_liana_resource.first().sum())} of 17")
print(f"custom pairs reaching the {CONFIG['min_patients_for_association']}-patient floor: "
      f"{int((cc.n_patients_evaluated >= CONFIG['min_patients_for_association']).sum())} rows")

# %% [markdown]
# > **LIANA does not change the frozen Stage-11 interpretation.** Its one BH-significant
# > consensus hit is structurally circular; the frozen targeted panel's own hit is not in
# > LIANA's resource at all; and no patient tier, structure state, phenotype state,
# > composition conclusion or coverage eligibility is modified by anything in this arm.
