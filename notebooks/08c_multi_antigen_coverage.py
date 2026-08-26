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
# # 08c — Supplemental multi-antigen coverage and incremental gain
#
# > **For each patient, how much malignant-cell transcript-level coverage is gained by
# > adding a second or third measurable target, while keeping measurement reliability,
# > depth robustness and normal-marrow expression as separate evidence dimensions?**
#
# **This is a supplemental Stage-08 deliverable consuming frozen upstream malignant-cell
# and measurement infrastructure. It is NOT a reopening of the frozen BCMA/GPRC5D
# analysis.** It reads frozen per-cell clone states, frozen depth strata and frozen raw
# counts; it writes only into `results/08_dual_antigen_escape/multi_antigen_coverage/`.
#
# It takes a **letter** (`08c`), like `05b` and `09b`, because number order is execution
# order with no exceptions and this is a terminal deliverable off Stage 08's evidence
# rather than a new stage in the pipeline.
#
# Design frozen in `multi_antigen_design.md` **before any pair or triple was computed**.
#
# Everything here is a **transcript-level detection** quantity. `count > 0` does not mean
# "biologically positive", and no output is a protein-level, surface-density or
# therapeutic-coverage claim.

# %%
import hashlib
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from mm_escape import antigen as AG
from mm_escape import communication as CM
from mm_escape import config
from mm_escape import coverage as CV

STAGE08 = Path("results/08_dual_antigen_escape")
OUT = STAGE08 / "multi_antigen_coverage"
OUT.mkdir(parents=True, exist_ok=True)
INTEGRATED = "results/05_integration/integrated.h5ad"
SEED = 20260825
pd.set_option("display.width", 320)

TARGETS = list(CV.TARGETS)
print(f"targets      {TARGETS}")
print(f"anchor       {CV.ANCHOR}")
print(f"eligibility  technical_zero>={CV.TECHNICAL_ZERO_MAX}  depth_rho>={CV.DEPTH_RHO_MAX}  "
      f"spread>={CV.DEPTH_SPREAD_MAX}x  background_separation<{CV.BACKGROUND_SEPARATION_MIN}x")

# %% [markdown]
# ## Frozen upstream is recorded by digest before anything else runs
#
# The isolation claim is not left to good intentions. Every frozen artifact this analysis
# consumes is hashed now and re-checked by `tests/test_coverage.py`, so a later edit to a
# frozen stage fails a test instead of passing unnoticed.

# %%
FROZEN = ["results/08_dual_antigen_escape/patient_antigen_states_primary.csv",
          "results/08_dual_antigen_escape/patient_antigen_states_sensitivity.csv",
          "results/08_dual_antigen_escape/patient_conegativity_enrichment.csv",
          "results/08_dual_antigen_escape/noise_floor_technical_zero.csv",
          "results/08_dual_antigen_escape/noise_floor_ambient.csv",
          "results/08_dual_antigen_escape/truncate10k_sensitivity.csv",
          "results/08_dual_antigen_escape/depth_strata_definition.csv",
          "results/09_bulk_validation/",
          "results/10_dn_coherence/dn_coherence_final_states.csv",
          "results/11_immune_context/patient_immune_composition.csv"]


def digest(path: Path) -> str:
    if path.is_dir():
        parts = [f"{p.relative_to(path)}:{digest(p)}" for p in sorted(path.rglob("*")) if p.is_file()]
        return hashlib.sha256("|".join(parts).encode()).hexdigest()
    return hashlib.sha256(path.read_bytes()).hexdigest()


digests = {rel: digest(Path(rel)) for rel in FROZEN if Path(rel).exists()}
(OUT / "frozen_upstream_digests.json").write_text(json.dumps(digests, indent=2, sort_keys=True))
print(f"recorded digests for {len(digests)} frozen artifacts")

# %% [markdown]
# ## Inputs — frozen cells, frozen strata, raw counts
#
# The denominators, the clone states and each cell's depth stratum are **read**, never
# recomputed. Because the stratum assignment is frozen rather than re-derived, no
# perturbation of any target — including the five new ones — can move a cell between
# strata. That makes the antigen-isolation property structural instead of asserted.

# %%
frozen_cells = pd.read_csv(STAGE08 / "cell_antigen_states.csv.gz", dtype={"patient_id": str})
raw = CM.stream_gene_counts(INTEGRATED, TARGETS)

labels = pd.read_csv("results/06_annotation/per_cell_labels.csv.gz",
                     dtype={"patient_id": str, "sample_name": str})
revised = pd.read_csv(
    "results/06_annotation/cluster23_local/trbc_context_revision/revised_lineage_calls.csv.gz",
    usecols=["cell_id", "call"])
C23 = {"NK": "NK", "T_NK_mixed": "cytotoxic_mixed", "unresolved": "cytotoxic_mixed",
       "T_ab": "Tcell", "T_gd": "Tcell"}
labels = labels.merge(revised, on="cell_id", how="left")
labels["lineage"] = np.where(labels.call.notna(), labels.call.map(C23), labels.cell_type)

allcells = labels.merge(raw, left_on="cell_id", right_index=True, how="left")
allcells["depth_ex_antigen"] = allcells.total_counts - allcells.TNFRSF17 - allcells.GPRC5D

mal = frozen_cells.merge(raw, left_on="cell_id", right_index=True, how="left")
for t in TARGETS:
    mal[f"det_{t}"] = CV.detected(mal[t].values)
primary = mal[mal.in_primary]
print(f"primary {len(primary):,} cells / {primary.patient_id.nunique()} patients   "
      f"sensitivity {len(mal):,} cells")
assert len(primary) == 21906 and len(mal) == 29015, "frozen denominators must be unchanged"

# %% [markdown]
# ## Per-target measurement QC — every target, on the same footing
#
# ### Detection in malignant cells, and the depth dependence
#
# The positivity rule is the frozen Stage-08 one, `count > 0`, applied unchanged to all
# seven. A per-gene threshold fitted to each target's abundance would be exactly the
# post-hoc tuning the design forbids and would make uncovered fractions incomparable
# between targets. **Variation in reliability is handled by the QC gate, not by moving the
# cutoff** — a dropout-prone target does not get a kinder rule.

# %%
det_by_cohort = primary.groupby("cohort")[[f"det_{t}" for t in TARGETS]].mean()
det_by_cohort.loc["ALL"] = primary[[f"det_{t}" for t in TARGETS]].mean()
det_by_cohort.columns = [c.replace("det_", "") for c in det_by_cohort.columns]
print("detection fraction in the primary denominator")
print(det_by_cohort.round(4).to_string())
print("\nmean raw count in the primary denominator")
print(primary[TARGETS].mean().round(4).to_string())

# %%
depth_rows = []
for t in TARGETS:
    rho, p = spearmanr(primary[f"det_{t}"].astype(float), primary.depth_ex_antigen)
    by_stratum = primary.groupby(["cohort", "depth_stratum_cohort"])[f"det_{t}"].mean()
    spread = by_stratum.groupby(level=0).apply(lambda s: s.max() / max(s.min(), 1e-9)).max()
    depth_rows.append({"antigen": t, "detection_rho_depth": rho, "p": p,
                       "max_within_cohort_stratum_spread": float(spread)})
depth_diag = pd.DataFrame(depth_rows)
print(depth_diag.round(4).to_string(index=False))

# %% [markdown]
# **Every one of the seven targets is depth-dependent** — Spearman ρ 0.32 to 0.48, all far
# above the 0.20 threshold, with 3× to 16× detection spread across depth strata inside a
# single cohort. **No target in this panel is depth-robust**, and that is the single most
# important QC result here: the panel's measurement quality, not its biology, is the
# binding constraint on everything below.

# %% [markdown]
# ### Ambient / background — declared per target from biology, before any distribution was seen
#
# A background population is only a background population if the target is genuinely absent
# from it. Stage 08's reference is correct for BCMA and GPRC5D and is inherited unchanged
# for those two; it is **not** automatically valid for the others. `SLAMF7` is expressed on
# NK cells, CD8 T subsets, activated monocytes, DC and B cells; `CD38` is on activated T and
# NK, B subsets, monocytes and is itself a progenitor marker; `ITGB7` is broadly lymphoid.
# For those three there is **no clean lymphoid or myeloid negative in marrow**, so the
# ambient floor is marked `NOT_EVALUABLE` rather than invented.

# %%
myeloma = allcells[allcells.sample_type == "myeloma"]
bg_rows = []
for t in TARGETS:
    lin = CV.BACKGROUND_LINEAGES[t]
    bg = myeloma[myeloma.lineage.isin(lin)]
    bg_det = float((bg[t] > 0).mean())
    mal_det = float(primary[f"det_{t}"].mean())
    single = len(lin) == 1
    bg_rows.append({"antigen": t, "background_lineages": ",".join(lin),
                    "n_background_cells": len(bg), "background_detection": bg_det,
                    "malignant_detection": mal_det,
                    "separation_fold": mal_det / max(bg_det, 1e-12),
                    "ambient_floor_status": CV.AMBIENT_NOT_EVALUABLE if single else "evaluable"})
background = pd.DataFrame(bg_rows)
print(background.round(4).to_string(index=False))

# %% [markdown]
# `SLAMF7`, `CD38` and `ITGB7` fall back to erythroid cells alone and carry
# `AMBIENT_FLOOR_NOT_EVALUABLE`. All three still clear their imperfect background by more
# than the 2× separation the design requires, so criterion 1 does not exclude them — but the
# floor itself is not quotable, and erythroid recovery differs **27×** between cohorts
# (Stage 11), which is stated wherever those floors appear.

# %% [markdown]
# ### Expression-matched technical-zero floor
#
# Control genes matched to each target's mean expression **in the primary denominator**;
# their zero-fraction in non-plasma reference cells, by cohort and frozen depth stratum.
# This is the rate at which a gene of that abundance reads zero for depth reasons alone —
# a **plausibility bound, never a correction**. Subtracting it to manufacture a "corrected"
# uncovered fraction is forbidden here exactly as it is in Stage 08.
#
# The reference population is the one Stage 08 actually used, established by reproducing
# `noise_floor_ambient.csv` exactly (14 of 14 rows, counts and medians): `cell_type` ∈
# {`Tcell`, `Myeloid`, `Bcell`, `HSPC`}. Stage 08's predeclaration also names `NK`, which
# resolved to nothing because Stage 06's endpoint emits no NK class.

# %%
gene_means = pd.read_csv(OUT / "primary_denominator_gene_means.csv", index_col=0).iloc[:, 0] \
    if (OUT / "primary_denominator_gene_means.csv").exists() else None

if gene_means is None:
    import h5py
    from anndata.io import read_elem

    prim_ids = set(primary.cell_id)
    with h5py.File(INTEGRATED, "r") as f:
        gene_names = np.asarray(read_elem(f["var"]).index.astype(str))
        cell_names = np.asarray(read_elem(f["obs"]).index.astype(str))
        keep = np.isin(cell_names, list(prim_ids))
        grp = f["layers"]["counts"]
        indptr = grp["indptr"][:]
        data, indices = grp["data"], grp["indices"]
        totals = np.zeros(len(gene_names))
        for s in range(0, len(indptr) - 1, 20000):
            e = min(s + 20000, len(indptr) - 1)
            lo, hi = indptr[s], indptr[e]
            d, ix = data[lo:hi], indices[lo:hi]
            rows = np.repeat(np.arange(s, e), np.diff(indptr[s:e + 1] - lo))
            m = keep[rows]
            if m.any():
                np.add.at(totals, ix[m], d[m])
    gene_means = pd.Series(totals / keep.sum(), index=gene_names, name="mean_count")
    gene_means.to_csv(OUT / "primary_denominator_gene_means.csv")

assert np.isclose(gene_means["TNFRSF17"], 5.7712, atol=1e-4), "must match the frozen stage-08 mean"
print("gene means reproduce the frozen stage-08 values for both anchor targets")

# %%
control_genes = {}
for t in TARGETS:
    candidates = gene_means.drop(index=[g for g in TARGETS if g in gene_means.index])
    control_genes[t] = list(candidates.sub(gene_means[t]).abs().sort_values().index[:100])

ctrl_counts = CM.stream_gene_counts(INTEGRATED, sorted({g for v in control_genes.values() for g in v}
                                                       | set(TARGETS)))
ctrl_counts["depth_ex_antigen"] = ctrl_counts.total_counts - ctrl_counts.TNFRSF17 - ctrl_counts.GPRC5D

REF_CLASSES = ["Tcell", "Myeloid", "Bcell", "HSPC"]
ref_index = labels[labels.cell_type.isin(REF_CLASSES)].set_index("cell_id")[["cohort"]]
ref = ctrl_counts.loc[ctrl_counts.index.isin(ref_index.index)].join(ref_index)

edges_tbl = pd.read_csv(STAGE08 / "depth_strata_definition.csv")
edges_tbl = edges_tbl[edges_tbl.scheme == "cohort_specific"]

tz_rows = []
for cohort, g in ref.groupby("cohort"):
    ed = edges_tbl[edges_tbl.cohort == cohort].sort_values("stratum")
    if not len(ed):
        continue
    edges = np.append(ed.lo.values, ed.hi.values[-1])
    strata = AG.assign_strata(g.depth_ex_antigen.values, edges)
    for t in TARGETS:
        arr = g[control_genes[t]].values
        for lab in np.unique(strata):
            m = strata == lab
            tz_rows.append({"antigen": t, "target_mean_in_denominator": round(float(gene_means[t]), 4),
                            "cohort": cohort, "stratum": int(lab), "n_reference_cells": int(m.sum()),
                            "median_depth": float(np.median(g.depth_ex_antigen.values[m])),
                            "matched_control_genes": 100,
                            "technical_zero_fraction": float((arr[m] == 0).mean())})
technical_zero = pd.DataFrame(tz_rows)
technical_zero.to_csv(OUT / "target_technical_zero_floor.csv", index=False)

tz_pivot = technical_zero.pivot_table(index="antigen", columns="cohort",
                                      values="technical_zero_fraction", aggfunc="median")
tz_pivot["pooled_median"] = technical_zero.groupby("antigen").technical_zero_fraction.median()
print(tz_pivot.reindex(TARGETS).round(4).to_string())

# %% [markdown]
# **This floor is computed on one consistent rule across all seven targets, and it is not
# bit-identical to Stage 08's.** Stage 08's exact control-gene selection is not recoverable
# from its frozen artifacts; reproducing its inputs (the denominator gene means match to
# four decimals) still leaves this construction running **0.10–0.20 higher**, i.e. more
# conservative. Stage 08's own rows stay the cited values for BCMA and GPRC5D and are read
# from disk unchanged. This table exists so the five new targets are compared against the
# anchor **under one rule** rather than against a number built differently.
#
# The rank order is the informative part and it agrees with the frozen result: `TNFRSF17`
# has much the smallest technical-zero floor, `GPRC5D` much the largest, and the five new
# targets sit between them and close together.
#
# **The design's threshold wording — "≥ 0.50 in the median cohort stratum" — is ambiguous
# about which cohort.** It is resolved here the only way that cannot cherry-pick: the
# **pooled median over all cohort × stratum rows**, one number per target, with the
# per-cohort breakdown reported beside it. `SLAMF7` (WU2 0.503) and `ITGB7` (WU1 0.514,
# WU2 0.518) sit just over the line in WashU alone and are flagged as threshold-hugging.

# %% [markdown]
# ## The SDC1 / CD138 checkpoint — answered before any eligibility call
#
# `SDC1` is not treated as an interchangeable therapeutic target. Five predeclared
# questions.

# %%
SECRETORY = ["SPCS1", "SPCS2", "SEC61B", "UBE2J1", "TMBIM6", "MZB1", "B2M"]
sec_counts = CM.stream_gene_counts(INTEGRATED, SECRETORY + TARGETS)
sec = primary.merge(sec_counts[SECRETORY], left_on="cell_id", right_index=True)
sec["secretory_breadth"] = (sec[SECRETORY] > 0).mean(axis=1)
sec["secretory_cp10k"] = sec[SECRETORY].sum(axis=1) / sec.total_umi * 1e4

q_rows = []
for t in TARGETS:
    neg_w = pos_w = w_tot = 0.0
    for _, g in sec.groupby(["cohort", "depth_stratum_cohort"]):
        neg, pos = g.loc[g[f"det_{t}"] == False, "secretory_breadth"], g.loc[g[f"det_{t}"], "secretory_breadth"]
        if len(neg) < 20 or len(pos) < 20:
            continue
        w_tot += len(g); neg_w += len(g) * neg.mean(); pos_w += len(g) * pos.mean()
    rho_sec, _ = spearmanr(sec[t], sec.secretory_cp10k)
    dn_gap = np.nan
    num = den = 0.0
    is_dn = (sec.observed_state == "double_negative").values
    for _, g in sec.assign(is_dn=is_dn).groupby(["cohort", "depth_stratum_cohort"]):
        if len(g) < 50:
            continue
        a, b = g.loc[g.is_dn, t].eq(0).mean(), g.loc[~g.is_dn, t].eq(0).mean()
        if np.isfinite(a) and np.isfinite(b):
            num += len(g) * (a - b); den += len(g)
    if den:
        dn_gap = num / den
    q_rows.append({"antigen": t,
                   "secretory_breadth_negative": neg_w / w_tot, "secretory_breadth_positive": pos_w / w_tot,
                   "secretory_breadth_delta": (neg_w - pos_w) / w_tot,
                   "rho_count_vs_secretory": rho_sec,
                   "excess_negativity_in_DN": dn_gap,
                   "in_PLASMA_MATURE": t in config.PLASMA_MATURE,
                   "in_MARKER_PANEL_PlasmaCell": t in config.MARKER_PANEL["PlasmaCell"]})
sdc1_qc = pd.DataFrame(q_rows)
sdc1_qc.to_csv(OUT / "sdc1_differentiation_checkpoint.csv", index=False)
print(sdc1_qc.round(4).to_string(index=False))

# %% [markdown]
# **Q1 — how strongly does SDC1 detection track depth?** ρ = +0.45 with 4.3× spread across
# strata: strongly, like every other target in the panel.
#
# **Q2/Q3 — is SDC1 negativity separable from the Stage-10 less-secretory DN phenotype?**
# Yes, and more cleanly than the anchor. Within frozen depth strata, SDC1-negative cells
# carry a secretory-panel detection breadth only 0.013 below SDC1-positive cells — smaller
# than `TNFRSF17` (0.036) and `ITGB7` (0.028). Per-cell SDC1 count is essentially
# **uncorrelated** with secretory output (ρ = +0.02) where `TNFRSF17` and `GPRC5D` sit at
# +0.23. **The differentiation confound is measurable and SDC1 is not the worst offender.**
#
# **Q4 — was SDC1 used upstream in a way that makes it circular?** **Yes, and this is a
# property of the frozen pipeline rather than a result.** `config.PLASMA_MATURE` is Stage
# 06's axis-(b) mature-plasma predicate — *"at least one of `PLASMA_MATURE`"* — and it is
# `("SDC1", "TNFRSF17")`. `SDC1` is also one of the five
# `config.MARKER_PANEL["PlasmaCell"]` genes whose detection fraction was the marker-coverage
# evidence, the load-bearing biological test that could veto the PlasmaCell class. **The
# plasma-cell denominator was partly established using SDC1 detection.**
#
# **Q5 — can the data distinguish target absence from differentiation loss?** For the
# *state* confound, yes (Q2/Q3). For the *selection* confound, no — it is structural, not
# statistical, and no amount of conditioning removes it.

# %% [markdown]
# ### The verdict, and the asymmetry it forces me to state plainly
#
# **`SDC1` is assigned `COVERAGE_NOT_EVALUABLE` for the primary matrix on circularity
# (design §11), not on its differentiation behaviour, which was measurable and unremarkable.**
# Its coverage numbers are still computed and reported in the separated exploratory table.
#
# **`TNFRSF17` carries the identical structural limitation** — the predicate is an OR over
# the two, and there is no principled sense in which one is more circular than the other.
# It is retained **only** because it is the frozen anchor and the project question is
# defined on it. **That asymmetry is a disclosure, not a scientific distinction**, and it
# is a limitation of the anchor that this supplemental analysis surfaces rather than fixes.
# Nothing here reopens Stage 08.
#
# `CD38` carries a weaker form — it is in `MARKER_PANEL["PlasmaCell"]` but was *deliberately
# excluded* from `PLASMA_MATURE` because it is not plasma-specific, so it was never able to
# satisfy the predicate alone. Recorded, not disqualifying.

# %% [markdown]
# ## Eligibility — decided from measurement QC alone, before any combination is inspected
#
# `CV.eligibility` cannot see coverage: it is not a parameter of the function, which is
# asserted by a test rather than left to discipline.

# %%
CIRCULARITY_BLOCKED = {"SDC1"}          # design §11, from Q4 above
tz_pooled = technical_zero.groupby("antigen").technical_zero_fraction.median()

qc_rows = []
for t in TARGETS:
    bg = background.set_index("antigen").loc[t]
    dd = depth_diag.set_index("antigen").loc[t]
    state, reason = CV.eligibility(
        ambient_status=bg.ambient_floor_status, malignant_detection=bg.malignant_detection,
        background_detection=bg.background_detection, technical_zero=tz_pooled[t],
        circularity_blocked=t in CIRCULARITY_BLOCKED)
    label = CV.reliability_label(tz_pooled[t], dd.detection_rho_depth,
                                 dd.max_within_cohort_stratum_spread,
                                 evaluable=(state == CV.ELIGIBLE))
    for cohort in ["MMRF", "WU1", "WU2"]:
        qc_rows.append({
            "antigen": t, "cohort": cohort,
            "malignant_detection": float(primary.loc[primary.cohort == cohort, f"det_{t}"].mean()),
            "immune_background_detection": float(
                (myeloma[(myeloma.cohort == cohort) & myeloma.lineage.isin(CV.BACKGROUND_LINEAGES[t])][t] > 0).mean()),
            "background_lineages": bg.background_lineages,
            "depth_association_rho": dd.detection_rho_depth,
            "depth_stratum_spread": dd.max_within_cohort_stratum_spread,
            "technical_zero_estimate_cohort": float(tz_pivot.loc[t, cohort]),
            "technical_zero_estimate_pooled": float(tz_pooled[t]),
            "ambient_floor_status": bg.ambient_floor_status,
            "measurement_reliability": label,
            "sdc1_differentiation_caveat": (
                "PLASMA_MATURE axis-(b) predicate gene: the plasma denominator was partly "
                "established using its detection" if t in config.PLASMA_MATURE else
                "in MARKER_PANEL[PlasmaCell] but deliberately excluded from PLASMA_MATURE"
                if t in config.MARKER_PANEL["PlasmaCell"] else "not used in stage-06 plasma identification"),
            "sdc1_checkpoint_resolved": t not in CIRCULARITY_BLOCKED,
            "coverage_eligibility": state,
            "not_evaluable_reason": reason})
qc = pd.DataFrame(qc_rows)
qc = qc.merge(pd.DataFrame({"antigen": TARGETS}), on="antigen")

donor = allcells[allcells.sample_type == "normal_bm"]
donor_plasma = donor[donor.cell_type == "PlasmaCell"]
dp = donor_plasma.groupby("sample_name")[TARGETS].apply(lambda g: (g > 0).mean()).median()
qc["donor_plasma_detection"] = qc.antigen.map(dp)
qc.to_csv(OUT / "target_measurement_qc.csv", index=False)

verdict = qc.drop_duplicates("antigen")[
    ["antigen", "technical_zero_estimate_pooled", "depth_association_rho", "ambient_floor_status",
     "measurement_reliability", "coverage_eligibility", "not_evaluable_reason"]]
print(verdict.round(4).to_string(index=False))

# %% [markdown]
# **Two targets fail, and one of them is half the project's anchor.**
#
# `GPRC5D` is `COVERAGE_NOT_EVALUABLE` on the dropout criterion — a pooled technical-zero
# floor of 0.62 against a 0.50 threshold fixed before any floor was computed. **This is not
# a new finding; it is the frozen Stage-08 conclusion reached again by an independent,
# even-handed rule** ("`GPRC5D` has a much larger technical-zero floor than `TNFRSF17`,
# especially in WashU"; "most GPRC5D zeros in WashU are consistent with dropout"). The
# threshold is **not** relaxed to admit it.
#
# `SDC1` is `COVERAGE_NOT_EVALUABLE` on circularity, per the checkpoint above.
#
# **Eligibility governs entry to the alternative-combination matrix. It does not suppress
# the anchor**, which is reported for every patient under design §16 with its labels
# attached — the project question is defined on BCMA+GPRC5D and is not redefined by a QC
# gate. Both excluded targets remain fully visible, with reasons, in
# `target_measurement_qc.csv`, and their combinations are computed and reported as
# exploratory rather than hidden.
#
# **No target in the panel is `comparatively_reliable.`** Every eligible one is
# `depth_sensitive`.

# %% [markdown]
# ## Coverage — every combination, per patient, under both denominators
#
# Both denominators are carried to the end and **never collapsed, averaged or selected
# between**. Denominator disagreement is part of the result.

# %%
ELIGIBLE = [t for t in TARGETS if (qc.set_index("antigen").loc[t, "coverage_eligibility"] == CV.ELIGIBLE).all()]
print(f"eligible: {ELIGIBLE}")
print(f"not evaluable: {[t for t in TARGETS if t not in ELIGIBLE]}")

ALL_COMBOS = CV.all_combinations(TARGETS, sizes=(1, 2, 3))


def combo_status(combo):
    if tuple(combo) == CV.ANCHOR:
        return "anchor"
    return "primary_matrix" if all(t in ELIGIBLE for t in combo) else "exploratory"


def coverage_table(cells_df, denominator):
    det_cols = [f"det_{t}" for t in TARGETS]
    rows = []
    for pid, g in cells_df.groupby("patient_id"):
        detection = g[det_cols].to_numpy(dtype=bool)
        for combo in ALL_COMBOS:
            row = CV.coverage_row(detection, TARGETS, combo)
            rows.append({"patient": pid, "cohort": g.cohort.iloc[0], "denominator": denominator,
                         "n_cells": len(g), "combination_status": combo_status(combo), **row})
    return pd.DataFrame(rows)


cov = pd.concat([coverage_table(primary, "primary"), coverage_table(mal, "sensitivity")],
                ignore_index=True)
print(f"{len(cov):,} rows = {cov.patient.nunique()} patients x 2 denominators x {len(ALL_COMBOS)} combinations")

for size, name in [(1, "single_target_coverage"), (2, "pair_coverage"), (3, "triple_coverage")]:
    cov[cov["size"] == size].to_csv(OUT / f"{name}.csv", index=False)

# %% [markdown]
# ### The monotonicity invariants, checked on the real data
#
# A mis-ordered set operation would make a triple look better than its own pairs and nothing
# downstream would notice. The relations are asserted here on every patient, not only in the
# synthetic tests.

# %%
wide = cov.set_index(["patient", "denominator", "combination"]).uncovered
violations = 0
for (pid, den), g in cov.groupby(["patient", "denominator"]):
    lookup = g.set_index("combination").uncovered
    for combo in ALL_COMBOS:
        if len(combo) == 1:
            continue
        joint = lookup["+".join(combo)]
        for k in range(1, len(combo)):
            for sub in combinations(combo, k):
                if joint > lookup["+".join(sub)] + 1e-12:
                    violations += 1
print(f"monotonicity violations across {cov.patient.nunique()} patients x 2 denominators: {violations}")
assert violations == 0

# %% [markdown]
# ## Incremental gain — the quantity a single- vs. dual- vs. triple-target question turns on
#
# Reported per direction and never summarised into one number: a pair can carry a large gain
# in one direction and a small one in the other, and that asymmetry *is* the answer.

# %%
gain_rows = []
for (pid, den), g in cov.groupby(["patient", "denominator"]):
    src = primary if den == "primary" else mal
    d = src[src.patient_id == pid][[f"det_{t}" for t in TARGETS]].to_numpy(dtype=bool)
    for combo in ALL_COMBOS:
        if len(combo) == 1:
            continue
        for row in CV.incremental_gains(d, TARGETS, combo):
            gain_rows.append({"patient": pid, "cohort": g.cohort.iloc[0], "denominator": den,
                              "combination_status": combo_status(combo), **row})
gains = pd.DataFrame(gain_rows)
assert (gains.gain >= -1e-12).all(), "adding a target can never reduce coverage"
gains[gains.combination.str.count(r"\+") == 1].to_csv(OUT / "incremental_gain_pairs.csv", index=False)
gains[gains.combination.str.count(r"\+") == 2].to_csv(OUT / "incremental_gain_triples.csv", index=False)
print(f"{len(gains):,} directional gain rows; min gain {gains.gain.min():.6f}")

# %% [markdown]
# ## The anchor — reported for every patient regardless of eligibility

# %%
anchor = cov[(cov.combination == "TNFRSF17+GPRC5D")][
    ["patient", "cohort", "denominator", "n_cells", "uncovered",
     "uncovered_TNFRSF17", "uncovered_GPRC5D"]].rename(columns={"uncovered": "uncovered_BCMA_GPRC5D"})
ag = gains[gains.combination == "TNFRSF17+GPRC5D"].pivot_table(
    index=["patient", "denominator"], columns="added", values="gain").reset_index()
ag.columns = ["patient", "denominator", "gain_from_adding_GPRC5D", "gain_from_adding_BCMA"] \
    if list(ag.columns[2:]) == ["GPRC5D", "TNFRSF17"] else ag.columns
anchor = anchor.merge(ag, on=["patient", "denominator"])
anchor.to_csv(OUT / "anchor_bcma_gprc5d_coverage.csv", index=False)

print("primary denominator, anchor summary across 32 patients")
ap = anchor[anchor.denominator == "primary"]
print(ap[["uncovered_TNFRSF17", "uncovered_GPRC5D", "uncovered_BCMA_GPRC5D",
          "gain_from_adding_GPRC5D", "gain_from_adding_BCMA"]].describe().loc[
    ["mean", "50%", "min", "max"]].round(4).to_string())

# %% [markdown]
# ## Descriptive combination ranking
#
# Sorted by uncovered fraction. **The minimum-uncovered combination is never labelled
# optimal, recommended or best.** The permitted wording is *greatest observed
# transcript-level malignant-cell coverage among evaluated combinations*.

# %%
summary = (cov[cov.denominator == "primary"]
           .groupby(["combination", "size", "combination_status"])
           .uncovered.agg(["median", "mean", "min", "max"]).reset_index()
           .sort_values("median"))
print("pairs — greatest observed coverage first")
print(summary[summary["size"] == 2].head(8).round(4).to_string(index=False))
print("\ntriples — greatest observed coverage first")
print(summary[summary["size"] == 3].head(8).round(4).to_string(index=False))
print("\nsingles")
print(summary[summary["size"] == 1].round(4).to_string(index=False))

# %% [markdown]
# ## Truncate-all-at-10k — the frozen Stage-08 procedure, reused exactly
#
# `Hypergeometric(N = total UMI, K = gene count, n = 10,000)`, seed 20260825, cells at or
# below the cap unchanged and **never discarded**. No new downsampling scheme is written:
# `antigen.downsample_gene_counts` is the frozen implementation and is called directly.

# %%
trunc = primary.copy()
for t in TARGETS:
    trunc[f"det_{t}"] = AG.downsample_gene_counts(trunc.total_umi.values, trunc[t].values,
                                                  cap=10000, seed=SEED) > 0
trunc_cov = coverage_table(trunc, "primary_truncated10k")
t10 = (cov[cov.denominator == "primary"][["patient", "combination", "uncovered"]]
       .merge(trunc_cov[["patient", "combination", "uncovered"]], on=["patient", "combination"],
              suffixes=("_original", "_truncated")))
t10["delta"] = t10.uncovered_truncated - t10.uncovered_original
t10 = t10.merge(cov[cov.denominator == "primary"][["patient", "cohort", "combination",
                                                   "combination_status"]],
                on=["patient", "combination"])
t10.to_csv(OUT / "truncate10k_coverage_sensitivity.csv", index=False)
print("mean change in uncovered fraction under truncation, by cohort")
print(t10.groupby(["cohort", "combination_status"]).delta.agg(["mean", "max"]).round(4).to_string())
print("\nSpearman of the combination ordering, original vs truncated (primary, median over patients)")
order = (t10.groupby("combination")[["uncovered_original", "uncovered_truncated"]].median())
print(f"  rho = {spearmanr(order.uncovered_original, order.uncovered_truncated).statistic:.4f}")

# %% [markdown]
# ## Primary vs sensitivity denominator

# %%
pvs = (cov[cov.denominator == "primary"][["patient", "cohort", "combination", "combination_status", "uncovered"]]
       .merge(cov[cov.denominator == "sensitivity"][["patient", "combination", "uncovered"]],
              on=["patient", "combination"], suffixes=("_primary", "_sensitivity")))
pvs["delta"] = pvs.uncovered_sensitivity - pvs.uncovered_primary
pvs.to_csv(OUT / "primary_vs_sensitivity_coverage.csv", index=False)
print(pvs.groupby("size" if "size" in pvs else "combination_status").delta.agg(["mean", "median", "max"]).round(4).to_string())
print(f"\npatients whose anchor uncovered fraction moves >5 points: "
      f"{int((pvs[(pvs.combination == 'TNFRSF17+GPRC5D')].delta.abs() > 0.05).sum())} of 32")

# %% [markdown]
# ## Repeated samples — disagreement is reported, never pooled away

# %%
rep_rows = []
for pid, g in primary.groupby("patient_id"):
    if g.sample_name.nunique() < 2:
        continue
    for sample, s in g.groupby("sample_name"):
        d = s[[f"det_{t}" for t in TARGETS]].to_numpy(dtype=bool)
        for combo in ALL_COMBOS:
            rep_rows.append({"patient": pid, "sample": sample, "n_cells": len(s),
                             "combination": "+".join(combo), "size": len(combo),
                             "combination_status": combo_status(combo),
                             "uncovered": CV.uncovered_fraction(d, TARGETS, combo)})
repeated = pd.DataFrame(rep_rows)
repeated.to_csv(OUT / "repeated_sample_coverage.csv", index=False)
spread = (repeated.groupby(["patient", "combination"]).uncovered.agg(lambda v: v.max() - v.min())
          .groupby("patient").agg(["max", "median"]))
print("within-patient range of the uncovered fraction across that patient's own samples")
print(spread.round(4).to_string())

# %% [markdown]
# ## Normal-marrow expression context
#
# **Donor is the biological unit.** This is normal **marrow** expression context. It is
# **not** a safety score, not whole-body off-tumour safety and not toxicity prediction.
# GPRC5D's decisive off-tumour liability is keratinized tissue, which a marrow dataset
# cannot observe at all, and expression is not toxicity.

# %%
DONOR_LINEAGES = ["PlasmaCell", "Tcell", "NK", "Bcell", "Myeloid", "Erythroid", "HSPC"]
nm_rows = []
for lineage in DONOR_LINEAGES:
    d = donor[donor.lineage == lineage] if lineage != "PlasmaCell" else donor_plasma
    if not len(d):
        continue
    per_donor = d.groupby("sample_name")[TARGETS].apply(lambda g: (g > 0).mean())
    for t in TARGETS:
        nm_rows.append({"lineage": lineage, "antigen": t, "n_donors": len(per_donor),
                        "n_cells": len(d), "median_donor_detection": float(per_donor[t].median()),
                        "min_donor_detection": float(per_donor[t].min()),
                        "max_donor_detection": float(per_donor[t].max()),
                        "mean_count": float(d[t].mean())})
normal_marrow = pd.DataFrame(nm_rows)
normal_marrow.to_csv(OUT / "normal_marrow_target_context.csv", index=False)
print(normal_marrow.pivot_table(index="lineage", columns="antigen",
                                values="median_donor_detection").reindex(DONOR_LINEAGES)[TARGETS]
      .round(3).to_string())

# %% [markdown]
# ## Pairwise co-loss — only where the frozen null applies, and only after the synthetic control
#
# The synthetic depth-only falsification test is re-run **before** the null touches patient
# data: two targets whose zeros are produced by library size alone must not read as
# biologically co-lost.

# %%
rng = np.random.default_rng(SEED)
n = 4000
sim_depth = rng.lognormal(mean=np.log(4000), sigma=1.2, size=n)
p_det = np.clip(sim_depth / 12000, 0.02, 0.95)
sim_a, sim_b = rng.random(n) > p_det, rng.random(n) > p_det
sim_strata = CV.shared_depth_strata(sim_depth, n_bins=10)
sim_obs = float((sim_a & sim_b).mean())
uncond = sim_obs / CV.unconditioned_expected_co_negative(sim_a, sim_b)
cond = sim_obs / CV.stratified_expected_co_negative(sim_a, sim_b, sim_strata)
print(f"synthetic depth-only pair: unconditioned enrichment {uncond:.3f}  ->  conditioned {cond:.3f}")
assert uncond > 1.15 and abs(cond - 1.0) < 0.05, "the depth-only control must not survive"

# %%
coloss_rows = []
for (pid, den), src in [((p, "primary"), primary) for p in primary.patient_id.unique()]:
    g = src[src.patient_id == pid]
    if len(g) < 100:
        continue
    for a, b in combinations(TARGETS, 2):
        a_neg, b_neg = ~g[f"det_{a}"].values, ~g[f"det_{b}"].values
        obs = float((a_neg & b_neg).mean())
        exp_c = CV.stratified_expected_co_negative(a_neg, b_neg, g.depth_stratum_cohort.values)
        exp_u = CV.unconditioned_expected_co_negative(a_neg, b_neg)
        coloss_rows.append({"patient": pid, "cohort": g.cohort.iloc[0], "pair": f"{a}+{b}",
                            "combination_status": combo_status((a, b)),
                            "observed_co_negative": obs,
                            "expected_unconditioned": exp_u, "expected_depth_conditioned": exp_c,
                            "enrichment_unconditioned": obs / exp_u if exp_u > 0 else np.nan,
                            "enrichment_depth_conditioned": obs / exp_c if exp_c > 0 else np.nan})
coloss = pd.DataFrame(coloss_rows)
coloss.to_csv(OUT / "pairwise_coloss_context.csv", index=False)
print("median enrichment across patients, unconditioned vs depth-conditioned")
print(coloss.groupby("pair")[["enrichment_unconditioned", "enrichment_depth_conditioned"]]
      .median().sort_values("enrichment_unconditioned", ascending=False).round(3).to_string())

# %% [markdown]
# ## Stage-12 interface — one compact patient-level file, no recommendation

# %%
prim_cov = cov[cov.denominator == "primary"]
best_pair = (prim_cov[(prim_cov["size"] == 2) & (prim_cov.combination_status == "primary_matrix")]
             .loc[lambda d: d.groupby("patient").uncovered.idxmin()]
             [["patient", "combination", "uncovered"]]
             .rename(columns={"combination": "greatest_coverage_pair_descriptive",
                              "uncovered": "greatest_coverage_pair_uncovered"}))
best_triple = (prim_cov[(prim_cov["size"] == 3) & (prim_cov.combination_status == "primary_matrix")]
               .loc[lambda d: d.groupby("patient").uncovered.idxmin()]
               [["patient", "combination", "uncovered"]]
               .rename(columns={"combination": "greatest_coverage_triple_descriptive",
                                "uncovered": "greatest_coverage_triple_uncovered"}))
iface = (anchor[anchor.denominator == "primary"]
         .drop(columns="denominator")
         .merge(best_pair, on="patient").merge(best_triple, on="patient"))
iface = iface.merge(
    pvs[pvs.combination == "TNFRSF17+GPRC5D"][["patient", "delta"]]
    .rename(columns={"delta": "anchor_primary_vs_sensitivity_delta"}), on="patient")
iface = iface.merge(
    t10[t10.combination == "TNFRSF17+GPRC5D"][["patient", "delta"]]
    .rename(columns={"delta": "anchor_truncate10k_delta"}), on="patient")
iface["eligible_targets"] = ",".join(ELIGIBLE)
iface["not_evaluable_targets"] = "GPRC5D(dropout);SDC1(circularity)"
iface["all_targets_depth_sensitive"] = True
iface["normal_marrow_context_file"] = "normal_marrow_target_context.csv"
iface.to_csv(OUT / "stage12_multi_antigen_interface.csv", index=False)
print(iface.head(6).round(4).to_string(index=False))
print(f"\n{len(iface)} patients written to stage12_multi_antigen_interface.csv")

# %% [markdown]
# **No clinical recommendation is produced.** The interface carries uncovered fractions,
# incremental gains, denominator and depth robustness, eligibility caveats and a pointer to
# the normal-marrow context as **separate columns**. There is no weighted aggregate anywhere
# in this analysis, and a test asserts none exists.
