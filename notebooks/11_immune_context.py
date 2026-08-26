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
# # Stage 11 — exploratory immune context
#
# > **Does immune context *accompany* the frozen DN measurement and coherence phenotypes?**
# > It does not decide whether those phenotypes are real, and it may not rescue, upgrade,
# > downgrade or create any patient classification.
#
# Ninth in the project's scientific hierarchy, and **exploratory by declaration**, for a
# reason fixed in advance: n ≈ 32 patients against a confounder (T/NK abundance) that is
# itself correlated with the predictor. That combination supports description, not
# inference.
#
# Design frozen in `results/11_immune_context/stage11_design.md` **before any immune data
# was read**. This notebook is the driver: it recomputes every Stage-11 table from the
# frozen upstream artifacts, and asserts it reproduces what is on disk.

# %%
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import kruskal

from mm_escape import communication as CM

OUT = Path("results/11_immune_context")
#: The first, paused Stage-11 run, preserved verbatim so the recomputation can be
#: checked against it rather than silently replacing it.
PRELIM = OUT / "preliminary_run"
INTEGRATED = "results/05_integration/integrated.h5ad"
pd.set_option("display.width", 300)

CATS = list(CM.IMMUNE_CATEGORIES)
print(f"categories {CATS}")
print(f"MIN_IMMUNE_CELLS={CM.MIN_IMMUNE_CELLS}  MIN_SENDER_CELLS={CM.MIN_SENDER_CELLS}  "
      f"confounders={CM.CONFOUNDERS}")

# %% [markdown]
# ## Immune categories come from the accepted Stage-06 endpoint — Stage 06 is not reopened
#
# Broad `cell_type` as-is, plus the accepted cluster-23 TRBC-context revision, which is what
# supplies a cytotoxic split at all. **`cytotoxic_mixed` is never folded into T or NK.**
# Forcing those 2,207 cells into one lineage would be inventing resolution the evidence does
# not carry — and the fragility check later shows the one candidate association does not
# survive being forced.

# %%
labels = pd.read_csv("results/06_annotation/per_cell_labels.csv.gz",
                     dtype={"patient_id": str, "sample_name": str})
revised = pd.read_csv(
    "results/06_annotation/cluster23_local/trbc_context_revision/revised_lineage_calls.csv.gz",
    usecols=["cell_id", "call"])

CLUSTER23_TO_CATEGORY = {"NK": "NK_core", "T_NK_mixed": "cytotoxic_mixed",
                         "unresolved": "cytotoxic_mixed", "T_ab": "Tcell", "T_gd": "Tcell"}
labels = labels.merge(revised, on="cell_id", how="left")
labels["immune_cat"] = np.where(labels.call.notna(),
                                labels.call.map(CLUSTER23_TO_CATEGORY), labels.cell_type)

counts_by_cat = labels.immune_cat.value_counts()
print(counts_by_cat.to_string())
assert counts_by_cat["NK_core"] == 8951 and counts_by_cat["cytotoxic_mixed"] == 2207
assert counts_by_cat["Tcell"] == 61162 and counts_by_cat["PlasmaCell"] == 35474

# %% [markdown]
# ## Depth is recomputed on the intersected gene space, not read off `obs`
#
# `obs["total_counts"]` was computed at QC time over each sample's **full** Cell Ranger
# reference (33,538 or 33,694 genes). Stage 05 then intersected to 32,991, and every
# downstream normalisation used the intersected space. The gap is small per cell but enough
# to move a per-patient median by up to 14 counts, so it is closed here rather than ignored.
#
# `depth_ex_antigen` subtracts `TNFRSF17` and `GPRC5D`, the same definition Stage 08 froze —
# the depth covariate must not be a function of the antigens whose absence is the predictor.

# %%
ANTIGENS = ["TNFRSF17", "GPRC5D"]
LR_GENES = sorted({g for pair in CM.LR_CANDIDATES for g in pair if g != "None"})
gene_counts = CM.stream_gene_counts(INTEGRATED, sorted(set(ANTIGENS) | set(LR_GENES)))
gene_counts["depth_ex_antigen"] = gene_counts.total_counts - gene_counts[ANTIGENS].sum(axis=1)
print(f"{len(gene_counts):,} cells x {len(LR_GENES)} LR genes + antigens")

import h5py
from anndata.io import read_elem

with h5py.File(INTEGRATED, "r") as f:
    stage05_obs = read_elem(f["obs"])
gene_counts["n_genes_by_counts"] = stage05_obs["n_genes_by_counts"].astype(float).values

cells = labels.merge(gene_counts, left_on="cell_id", right_index=True, how="left")
assert cells.depth_ex_antigen.notna().all() and cells.n_genes_by_counts.notna().all()

# %% [markdown]
# ## One denominator, frozen: all non-plasma annotated cells
#
# Every lineage fraction uses the same denominator, so no comparison is ever made between
# two differently-normalised quantities. Because proportions are compositional — one lineage
# rising forces the rest down — features enter every model as **centred log-ratios**, and no
# lineage is interpreted as moving independently of the others.
#
# Evaluability is the frozen Stage-08 `MIN_PATIENT_CELLS` (100), reused rather than
# re-derived. Per patient we keep the cell counts and sample count, so no fraction is later
# treated as if it were as precise as any other.

# %%
patients = pd.read_csv("results/08_dual_antigen_escape/patient_antigen_states_primary.csv",
                       dtype={"patient": str})
sensitivity = pd.read_csv("results/08_dual_antigen_escape/patient_antigen_states_sensitivity.csv",
                          dtype={"patient": str})

records = []
for pid in patients.patient:
    d = cells[cells.patient_id == pid]
    imm = d[d.immune_cat != "PlasmaCell"]
    n_imm = len(imm)
    rec = {"patient": pid, "cohort": d.cohort.iloc[0], "n_immune_cells": n_imm,
           "n_annotated_cells": len(d),
           "n_plasma_cells": int((d.immune_cat == "PlasmaCell").sum()),
           "n_samples": d.sample_name.nunique(),
           "median_depth_ex_antigen_immune": float(np.median(imm.depth_ex_antigen)),
           "median_genes_immune": float(np.median(imm.n_genes_by_counts)),
           "evaluable": n_imm >= CM.MIN_IMMUNE_CELLS}
    for c in CATS:
        n = int((imm.immune_cat == c).sum())
        rec[f"n_{c}"], rec[f"frac_{c}"] = n, n / n_imm
    records.append(rec)

comp = pd.DataFrame(records)
clr_values = CM.clr(comp[[f"frac_{c}" for c in CATS]].values)
for j, c in enumerate(CATS):
    comp[f"clr_{c}"] = clr_values[:, j]

print(f"{len(comp)} patients; evaluable {int(comp.evaluable.sum())}")
assert comp.evaluable.all(), "all 32 clear the >=100 immune-cell floor"

# %% [markdown]
# The recomputation is checked against the table produced when Stage 11 first ran, so the
# preliminary outputs are shown to be reproducible rather than merely preserved.

# %%
frozen_comp = pd.read_csv(PRELIM / "patient_immune_composition.csv", dtype={"patient": str})
shared = [c for c in frozen_comp.columns if c in comp.columns]
check = frozen_comp[shared].merge(comp[shared], on="patient", suffixes=("_frozen", "_new"))
drift = {c: float(np.abs(check[f"{c}_frozen"].astype(float) - check[f"{c}_new"].astype(float)).max())
         for c in shared if c != "patient" and frozen_comp[c].dtype.kind in "if"}
worst = float(np.nanmax(list(drift.values())))
print(f"max drift vs the preserved preliminary table: {worst:.2e}  "
      f"(columns compared: {len(drift)})")
assert np.isfinite(list(drift.values())).all(), "a column failed to compare at all"
assert worst < 1e-4, drift

# %% [markdown]
# ## The confound model was frozen before any association was computed
#
# Stage 08 and Stage 10 each produced strong apparent biology out of sequencing depth, so
# Stage 11 assumes the same risk by default:
#
#     CLR(immune feature) ~ DN predictor + cohort + log10(immune depth)
#                                        + log10(n immune cells) + log10(n samples)
#
# The **unadjusted estimate is always reported beside the adjusted one** — the gap between
# them is the informative part, exactly as Stages 08 and 10 report unconditioned beside
# conditioned. Each feature is also screened against depth and against cohort **first**, so
# a depth-tracking feature is flagged before its DN association is ever read.

# %%
predictors = (patients[["patient", "observed_double_negative_fraction", "excess_dn",
                        "enrichment_stratified"]]
              .rename(columns={"observed_double_negative_fraction": "obs_dn_primary",
                               "enrichment_stratified": "enr_cohortbins"})
              .merge(sensitivity[["patient", "observed_double_negative_fraction"]]
                     .rename(columns={"observed_double_negative_fraction": "obs_dn_sensitivity"}),
                     on="patient"))
D = comp.merge(predictors, on="patient")
CONF = [D.cohort.values,
        np.log10(D.median_depth_ex_antigen_immune.values),
        np.log10(D.n_immune_cells.values.astype(float)),
        np.log10(D.n_samples.values.astype(float))]

depth_rows, cohort_rows = [], []
for c in CATS:
    a = CM.ols_association(D[f"clr_{c}"].values,
                           np.log10(D.median_depth_ex_antigen_immune.values), ())
    depth_rows.append({"feature": c, "n": a["n"], "coef": a["coef"], "p": a["p"]})
    H, p = kruskal(*[D.loc[D.cohort == k, f"clr_{c}"].values for k in sorted(D.cohort.unique())])
    row = {"feature": c, "kruskal_H": H, "p": p}
    row.update({f"median_{k}": float(D.loc[D.cohort == k, f"frac_{c}"].median())
                for k in sorted(D.cohort.unique())})
    cohort_rows.append(row)

depth_diag = pd.DataFrame(depth_rows)
depth_diag["p_bh"] = CM.benjamini_hochberg(depth_diag.p.values)
depth_diag["tracks_depth"] = depth_diag.p < 0.05
cohort_diag = pd.DataFrame(cohort_rows)
cohort_diag["p_bh"] = CM.benjamini_hochberg(cohort_diag.p.values)
cohort_diag["differs_by_cohort"] = cohort_diag.p_bh < 0.05

print(depth_diag.round(4).to_string(index=False))
print()
print(cohort_diag.round(4).to_string(index=False))

# %% [markdown]
# **Three of seven features track depth on their own** and **four differ by cohort.** The
# erythroid gap is extreme — a median fraction of 0.296 in MMRF against 0.011 in WU1, a 27×
# difference that is a cell-recovery property of the deposit, not marrow biology.
#
# **`NK_core` — the only feature that will show any DN association — is itself both
# depth-tracking and cohort-varying.** Stating that here, before its association is read, is
# the point of running the diagnostics first.

# %% [markdown]
# ## Q1–Q3: composition, T/NK, myeloid
#
# Four predictors × seven features. All four measurement predictors are reported; **none was
# selected for correlating best**, which is why the table below is printed whole.

# %%
rows = []
for pred in CM.MEASUREMENT_PREDICTORS:
    for c in CATS:
        y, x = D[f"clr_{c}"].values, D[pred].values
        un = CM.ols_association(y, x, ())
        adj = CM.ols_association(y, x, CONF)
        per_cohort = CM.within_cohort_spearman(y, x, D.cohort.values)
        signs = [np.sign(v["rho"]) for v in per_cohort.values() if np.isfinite(v.get("rho", np.nan))]
        rows.append({"predictor": pred, "feature": c, "n": adj["n"],
                     "coef_unadjusted": un["coef"], "p_unadjusted": un["p"],
                     "coef_adjusted": adj["coef"], "ci_lo": adj["ci_lo"], "ci_hi": adj["ci_hi"],
                     "p_adjusted": adj["p"],
                     **{f"rho_{k}": v.get("rho", np.nan) for k, v in per_cohort.items()},
                     **{f"n_{k}": v["n"] for k, v in per_cohort.items()},
                     "within_cohort_same_sign": bool(len(set(signs)) == 1 and len(signs) == 3)})
measurement = pd.DataFrame(rows)
measurement["p_adj_BH"] = CM.benjamini_hochberg(measurement.p_adjusted.values)
measurement["p_unadj_BH"] = CM.benjamini_hochberg(measurement.p_unadjusted.values)

print(f"raw p<0.05 unadjusted: {(measurement.p_unadjusted < 0.05).sum()}  "
      f"adjusted: {(measurement.p_adjusted < 0.05).sum()}  "
      f"BH<0.10 adjusted: {(measurement.p_adj_BH < 0.10).sum()}")
print(measurement[measurement.p_adjusted < 0.05][
    ["predictor", "feature", "coef_adjusted", "ci_lo", "ci_hi", "p_adjusted", "p_adj_BH",
     "rho_MMRF", "rho_WU1", "rho_WU2", "within_cohort_same_sign"]].round(4).to_string(index=False))

# %% [markdown]
# **Nothing survives multiple-testing correction.** Two of 28 reach raw p < 0.05 after
# confound control — `NK_core` against both DN denominators, pointing the same way in all
# three cohorts — and both sit at BH ≈ 0.49.
#
# A naive association that the guard **did** catch: `NK_core` × `enr_cohortbins`,
# unadjusted p = 0.019 → adjusted p = 0.070.

# %%
print(measurement[(measurement.predictor == "enr_cohortbins") & (measurement.feature == "NK_core")]
      [["predictor", "feature", "coef_unadjusted", "p_unadjusted", "coef_adjusted", "p_adjusted"]]
      .round(4).to_string(index=False))

# %% [markdown]
# ### The candidate direction is fragile, and the frozen design required checking that
#
# The whole `NK_core` signal depends on how the 2,207 `cytotoxic_mixed` cells are handled.
# Re-expressed as a plain log-fraction, or with the mixed cells folded in, it goes away.
# **Reported as a candidate direction, not a finding.**

# %%
for label, frac in [("NK_core alone", D.frac_NK_core),
                    ("NK_core + cytotoxic_mixed", D.frac_NK_core + D.frac_cytotoxic_mixed)]:
    a = CM.ols_association(np.log(frac.values + 1e-6), D.obs_dn_primary.values, CONF)
    print(f"{label:28s} coef {a['coef']:+.3f}  p {a['p']:.3f}")

# %% [markdown]
# ## Q1–Q3 against the coherence axis: nothing
#
# The two DN axes stay separate. Measurement asks *how much* observed DN there is; coherence
# asks whether the DN cells are *organised*. Stage 10 showed those are distinct properties,
# so combining them into one predictor would ask a question neither axis answers.
# `DN_COHERENCE_NOT_EVALUABLE` patients are excluded from the test and reported, never
# silently folded into "not supported".

# %%
coherence = pd.read_csv("results/10_dn_coherence/dn_coherence_final_states.csv",
                        dtype={"patient": str})[["patient", "dn_coherence_state"]]
K = D.merge(coherence, on="patient")
print(K.dn_coherence_state.value_counts().to_string())
K = K[K.dn_coherence_state != "DN_COHERENCE_NOT_EVALUABLE"].copy()
supported = (K.dn_coherence_state == "DN_COHERENCE_SUPPORTED").astype(float).values
K_CONF = [K.cohort.values, np.log10(K.median_depth_ex_antigen_immune.values),
          np.log10(K.n_immune_cells.values.astype(float)),
          np.log10(K.n_samples.values.astype(float))]

rows = []
for c in CATS:
    f = K[f"clr_{c}"].values
    un, adj = CM.ols_association(f, supported, ()), CM.ols_association(f, supported, K_CONF)
    rows.append({"feature": c, "n": adj["n"], "n_supported": int(supported.sum()),
                 "median_supported": float(np.median(K.loc[supported == 1, f"frac_{c}"])),
                 "median_not_supported": float(np.median(K.loc[supported == 0, f"frac_{c}"])),
                 "coef_unadjusted": un["coef"], "p_unadjusted": un["p"],
                 "coef_adjusted": adj["coef"], "ci_lo": adj["ci_lo"], "ci_hi": adj["ci_hi"],
                 "p_adjusted": adj["p"]})
coh_assoc = pd.DataFrame(rows)
coh_assoc["p_adj_BH"] = CM.benjamini_hochberg(coh_assoc.p_adjusted.values)
print(coh_assoc.round(4).to_string(index=False))

# %% [markdown]
# No feature reaches raw p < 0.05 adjusted; the smallest is `Bcell` at p = 0.077, BH 0.27.
# The medians move in directions that look interesting — `NK_core` 0.023 in supported against
# 0.075 in not-supported — but with **four supported patients** that is not evidence, which
# is exactly why permutation control **C2** was frozen in advance for a 4-vs-23 split.

# %% [markdown]
# ## Repeated samples — within-patient variation exceeds most between-patient differences
#
# Samples are summarised separately first and **never counted as independent patients**.

# %%
rep_rows = []
for pid, d in cells[cells.patient_id.isin(D.patient)].groupby("patient_id"):
    if d.sample_name.nunique() < 2:
        continue
    for sample, s in d.groupby("sample_name"):
        imm = s[s.immune_cat != "PlasmaCell"]
        rec = {"patient": pid, "sample": sample, "n_immune_cells": len(imm),
               "evaluable": len(imm) >= CM.MIN_IMMUNE_CELLS}
        for c in CATS:
            rec[f"frac_{c}"] = (imm.immune_cat == c).sum() / len(imm) if len(imm) else np.nan
        rep_rows.append(rec)
repeated = pd.DataFrame(rep_rows)
spread = (repeated.groupby("patient")[[f"frac_{c}" for c in CATS]]
          .agg(lambda v: v.max() - v.min()))
print(spread.round(3).to_string())

# %% [markdown]
# `Tcell` ranges 0.468 within `58408` and 0.348 within `27522`; `Myeloid` 0.438 within
# `58408`. `83942` — one patient, two protocols — is the stable case (`Tcell` range 0.035).
# **Variation of that size within a single patient is a further reason the n ≈ 32
# cross-patient associations are weak evidence.**

# %% [markdown]
# ## Q4 — communication, discovery-only, and a deviation from the frozen design
#
# A full-interactome LIANA/CellChat screen was **deliberately not run**: n ≈ 32 against
# hundreds of interactions, with the primary confounder correlated with the predictor,
# cannot support inference, and screening it would manufacture an appealing pathway list.
# A **predeclared 17-pair candidate set** is scored instead — sender immune category →
# receiver plasma population, pseudobulk CPM, score = ligand × receptor.
#
# > **Amendment, documented rather than quietly applied.** The frozen design fixed the
# > receiver as *the patient's frozen clone plasma population*. The first Stage-11 run used
# > **all** Stage-06 plasma cells instead. That is a departure from the pre-registered
# > design, found while resuming the stage, and it is corrected here in the direction the
# > design specified. **Both versions are computed and reported**, because the correction
# > was made after the first result was seen and hiding either one would be the failure mode
# > the pre-registration exists to prevent.
#
# `PRF1` and `GZMB` carry no receptor in the candidate set and drop out, leaving 17 pairs.
# `TNFRSF17` appears nowhere in it — the one pair that would have named it as a receptor is
# `TNFSF13`/`TNFSF13B` → `TNFRSF13B`.

# %%
pairs = [(l, r) for l, r in CM.LR_CANDIDATES if r != "None"]
print(f"{len(pairs)} pairs x {len(CM.LR_SENDERS)} senders = {len(pairs) * len(CM.LR_SENDERS)} tests")
assert "TNFRSF17" not in {r for _, r in pairs} and "GPRC5D" not in LR_GENES

antigen_cells = pd.read_csv("results/08_dual_antigen_escape/cell_antigen_states.csv.gz",
                            dtype={"patient_id": str}, usecols=["cell_id", "in_primary"])
clone_primary = set(antigen_cells.loc[antigen_cells.in_primary, "cell_id"])

sender_pool = cells[cells.immune_cat.isin(CM.LR_SENDERS)]
sender_cpm = sender_pool.groupby(["patient_id", "immune_cat"]).apply(
    lambda g: CM.pseudobulk_cpm(g, LR_GENES), include_groups=False)
sender_n = sender_pool.groupby(["patient_id", "immune_cat"]).size()

RECEIVERS = {"clone_primary": cells.cell_id.isin(clone_primary),
             "all_plasma": cells.immune_cat == "PlasmaCell"}


def score_table(receiver_mask):
    recv = cells[receiver_mask].groupby("patient_id").apply(
        lambda g: CM.pseudobulk_cpm(g, LR_GENES), include_groups=False)
    cohort_patients = set(D.patient)
    out = []
    for (pid, sender), srow in sender_cpm.iterrows():
        if (pid not in cohort_patients or pid not in recv.index
                or sender_n[(pid, sender)] < CM.MIN_SENDER_CELLS):
            continue
        rrow = recv.loc[pid]
        for lig, rec in pairs:
            out.append({"patient": pid, "sender": sender, "ligand": lig, "receptor": rec,
                        "ligand_cpm_sender": srow[lig], "receptor_cpm_receiver": rrow[rec],
                        "lr_score": srow[lig] * rrow[rec]})
    return pd.DataFrame(out)


def score_vs_dn(ctx):
    out = []
    for (sender, lig, rec), g in ctx.groupby(["sender", "ligand", "receptor"]):
        m = D.merge(g[["patient", "lr_score"]], on="patient")
        a = CM.ols_association(
            np.log1p(m.lr_score.values), m.obs_dn_primary.values,
            [m.cohort.values, np.log10(m.median_depth_ex_antigen_immune.values),
             np.log10(m.n_immune_cells.values.astype(float)),
             np.log10(m.n_samples.values.astype(float))])
        out.append({"sender": sender, "ligand": lig, "receptor": rec, "n": a["n"],
                    "coef_adjusted": a["coef"], "p_adjusted": a["p"]})
    res = pd.DataFrame(out)
    res["p_BH"] = CM.benjamini_hochberg(res.p_adjusted.values)
    return res


contexts = {k: score_table(v) for k, v in RECEIVERS.items()}
vs_dn = {k: score_vs_dn(v) for k, v in contexts.items()}
for k, v in vs_dn.items():
    print(f"{k:14s} patients {contexts[k].patient.nunique():2d}  rows {len(contexts[k]):4d}  "
          f"raw p<0.05 {int((v.p_adjusted < 0.05).sum()):2d}  BH<0.10 {int((v.p_BH < 0.10).sum())}")

# %% [markdown]
# The deviating `all_plasma` run is checked against the preserved preliminary table, so the
# amendment is shown to be the *only* difference between the two.

# %%
frozen_ctx = pd.read_csv(PRELIM / "communication_context.csv", dtype={"patient": str})
chk = frozen_ctx.merge(contexts["all_plasma"], on=["patient", "sender", "ligand", "receptor"],
                       suffixes=("_frozen", "_new"))
assert len(chk) == len(frozen_ctx), "every preserved row reproduces"
print(f"max |Δ lr_score| vs the preserved table: "
      f"{np.abs(chk.lr_score_frozen - chk.lr_score_new).max():.2e}")

# %%
print("highest median LR scores (design-conformant receiver):")
print(contexts["clone_primary"].groupby(["sender", "ligand", "receptor"]).lr_score.median()
      .sort_values(ascending=False).head(5).round(2).to_string())
print()
for k in ("all_plasma", "clone_primary"):
    print(f"--- {k}")
    print(vs_dn[k].sort_values("p_adjusted").head(5).round(5).to_string(index=False))

# %% [markdown]
# Correcting the receiver to the frozen definition makes the panel look **stronger**, not
# weaker: `Tcell PDCD1 → CD274` goes from BH 0.043 to BH 0.0026, and four TRAIL-axis
# interactions join it below BH 0.10.
#
# **That is a warning, not a result.** Every one of those coefficients is negative, and they
# span unrelated receptors. The obvious explanation is not five independent immune axes.

# %% [markdown]
# ### The receiver-side check — and the fourth instance of this project's recurring lesson
#
# The receptor side of every score is measured on the *same plasma cells* whose DN status is
# the predictor. Stage 10's frozen conclusion is that DN cells sit in a **less secretory,
# less differentiated** plasma-cell state, lower across ER/secretory machinery and `B2M`. A
# patient with a high DN fraction therefore has plasma cells with globally lower surface-
# receptor transcript abundance — which would drive *every* receptor down at once, with no
# immune signalling involved.
#
# So the receptor term is tested on its own, against the same adjusted model.

# %%
recv_rows = []
receptors = sorted({r for _, r in pairs})
for name, mask in RECEIVERS.items():
    recv = cells[mask].groupby("patient_id").apply(
        lambda g: CM.pseudobulk_cpm(g, LR_GENES), include_groups=False)
    m = D.merge(recv[receptors].reset_index().rename(columns={"patient_id": "patient"}),
                on="patient")
    conf = [m.cohort.values, np.log10(m.median_depth_ex_antigen_immune.values),
            np.log10(m.n_immune_cells.values.astype(float)),
            np.log10(m.n_samples.values.astype(float))]
    for r in receptors:
        a = CM.ols_association(np.log1p(m[r].values), m.obs_dn_primary.values, conf)
        recv_rows.append({"receiver": name, "receptor": r, "n": a["n"],
                          "coef_adjusted": a["coef"], "p_adjusted": a["p"]})
receiver_only = pd.DataFrame(recv_rows)
for name, g in receiver_only.groupby("receiver"):
    print(f"{name:14s} {int((g.coef_adjusted < 0).sum())}/{len(g)} receptors negative vs DN, "
          f"raw p<0.05: {int((g.p_adjusted < 0.05).sum())}")
print()
print(receiver_only[receiver_only.receiver == "clone_primary"]
      .sort_values("p_adjusted").head(6).round(5).to_string(index=False))

# %% [markdown]
# **Eleven of fifteen receptors move down with DN burden under the design-conformant
# receiver, five of them at raw p < 0.05, across receptors with nothing biologically in
# common.** The communication panel is reading a broad plasma-cell-state shift on its
# receptor term, not a specific immune axis — the same confound Stage 10 already froze as
# its Level-2 conclusion, arriving here through a different door.
#
# > **`Tcell PDCD1 → CD274` is therefore NOT reported as a candidate immune-evasion axis.**
# > It is one interaction from a 68-test discovery-only panel whose receptor term is
# > confounded with the predictor by construction. It touches no patient's classification.
#
# This is the **fourth** independent place in the project where an unconditioned or naive
# analysis pointed at appealing biology that a control removed — after Stage 08's
# co-negativity enrichment, Stage 10's Level-1 DN structure and Stage 10's Level-2
# transcriptional phenotype.

# %% [markdown]
# ## Outputs

# %%
comp.to_csv(OUT / "patient_immune_composition.csv", index=False)
depth_diag.to_csv(OUT / "immune_feature_depth_diagnostic.csv", index=False)
cohort_diag.to_csv(OUT / "immune_feature_cohort_diagnostic.csv", index=False)
measurement.to_csv(OUT / "immune_vs_dn_measurement.csv", index=False)
coh_assoc.to_csv(OUT / "immune_vs_dn_coherence.csv", index=False)
repeated.to_csv(OUT / "repeated_sample_immune_context.csv", index=False)
contexts["clone_primary"].to_csv(OUT / "communication_context.csv", index=False)
contexts["all_plasma"].to_csv(OUT / "communication_context_all_plasma.csv", index=False)
vs_dn["clone_primary"].to_csv(OUT / "communication_context_vs_dn.csv", index=False)
vs_dn["all_plasma"].to_csv(OUT / "communication_context_vs_dn_all_plasma.csv", index=False)
receiver_only.to_csv(OUT / "communication_receiver_side_confound.csv", index=False)

tnk = D[["patient", "cohort", "n_immune_cells", "frac_Tcell", "frac_NK_core",
         "frac_cytotoxic_mixed", "obs_dn_primary", "excess_dn"]].merge(coherence, on="patient")
tnk["frac_T_plus_NK"] = tnk.frac_Tcell + tnk.frac_NK_core
tnk["frac_cytotoxic_all"] = tnk.frac_NK_core + tnk.frac_cytotoxic_mixed
tnk.to_csv(OUT / "tnk_context.csv", index=False)
D[["patient", "cohort", "frac_Myeloid", "n_Myeloid", "obs_dn_primary", "excess_dn"]].merge(
    coherence, on="patient").to_csv(OUT / "myeloid_context.csv", index=False)
print("written")

# %% [markdown]
# ## What Stage 11 concludes
#
# > **No immune-composition association survives multiple-testing correction, and the one
# > communication interaction that does is explained by a receiver-side confound.**
#
# Stage 11 changed no tier, rescued no patient and created no third classifier. The final
# synthesis carries measurement robustness, DN coherence and exploratory immune context as
# three distinct evidence layers — never a combined score.
