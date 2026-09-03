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
#     display_name: mm-annotation
#     language: python
#     name: mm-annotation
# ---

# %% [markdown]
# # 06 — Cell type annotation
#
# **Env: `mm-annotation`.** Reads `results/05_integration/integrated.h5ad`, writes
# `results/06_annotation/`.
#
# ## What this stage owes downstream, and nothing else
#
# One column: `obs["cell_type"]`, seven coarse classes plus `Ambiguous`. Stages 07-12
# read it and never branch on annotation logic. That decoupling is what lets this
# comparison be redone or reversed later without touching anything else.
#
# Three methods run separately and the winner is chosen **per class** against bars
# declared before any result was looked at. The methods are expected to fail on
# *different* classes, so one verdict for the whole stage would throw away good labels
# to punish an unrelated weakness.
#
# ## The two things that decide this stage
#
# **Concordance is not accuracy.** The manual labels are a third opinion from the same
# expression matrix, not ground truth, so F1-against-manual measures *agreement*. Two
# automated methods agreeing is the strongest agreement available and still only
# agreement, between references that share marker-biology priors and may share their
# blind spots.
#
# **So marker coverage is a veto.** A class whose assigned cells do not express that
# class's markers is rejected regardless of concordance. `MARKER_COVERAGE_MIN = 0.30`
# was recorded in the main project document before any coverage number was computed, for
# the same
# reason the F1 bars were.
#
# ## What the plasma-cell boundary needs
#
# Stage 07's malignant calling — and therefore `frac_double_negative`'s denominator —
# rests on it. Stage 05 split plasma cells into three cohort-shaped clusters because
# WashU was cut at 10,000 UMIs before deposit, so the plausible failure mode here is
# good annotation on donor marrow plus systematic under-recognition of plasma cells in
# a censored myeloma cohort. Marker coverage is therefore reported by cohort and by
# stage-05 cluster, not only overall, and **cross-cohort plasma geometry is not read as
# malignant-state biology** — it is an ascertainment artefact.

# %%
from __future__ import annotations

import gc
import json
import resource
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

REPO = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
sys.path.insert(0, str(REPO / "src"))

from mm_escape import annotation as ann          # noqa: E402
from mm_escape import config                     # noqa: E402

IN_H5AD = REPO / "results" / "05_integration" / "integrated.h5ad"
OUT = REPO / "results" / "06_annotation"
OUT.mkdir(parents=True, exist_ok=True)

sc.settings.verbosity = 1
CLUSTER_KEY = "leiden"


def rss(tag: str) -> None:
    """Peak RSS so far. This stage loads the same object stage 05 peaked ~20 GB on,
    and the box has 30 GB — memory is a real constraint here, not a footnote."""
    gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2
    print(f"[mem] peak RSS {gb:5.1f} GB  after {tag}", flush=True)

print("scanpy", sc.__version__)
print("classes  ", ann.BROAD_CLASSES)
print("F1 bars  ", config.CONCORDANCE_THRESHOLDS)
print("coverage ", config.MARKER_COVERAGE_MIN, "(declared before any result)")

# %%
# Stage 06 ran twice. v1 is preserved at `results/06_annotation_v1/` and was NOT
# accepted — see the main project document and docs/decisions-archive.md. v2 revises the
# manual
# reference and adds the lineage-exclusivity veto; the acceptance thresholds are
# unchanged.
#
# CellTypist and SingleR are reused from v1 rather than recomputed. That is legitimate
# here and only here: their inputs (the stage-05 expression matrix), their models and
# references, and the prediction procedure are all untouched by this revision, which
# changes only the MANUAL panel and the validation framework. The barcode identity is
# asserted below, so a silently different input fails loudly instead of being reused.
V_PREV = REPO / "results" / "06_annotation_c2d_accepted" / "annotated.h5ad"
REUSED = ("cell_type_celltypist", "cell_type_singler_nov", "cell_type_singler_hpca",
          "celltypist_raw", "celltypist_percell", "celltypist_conf",
          "singler_nov_raw", "singler_hpca_raw")

if V_PREV.exists():
    adata = sc.read_h5ad(V_PREV)
    missing = [c for c in REUSED if c not in adata.obs]
    assert not missing, f"previous predictions incomplete: {missing}"
    ref_obs = sc.read_h5ad(IN_H5AD, backed="r").obs_names
    assert list(adata.obs_names) == list(ref_obs), "barcodes differ from stage 05 — do not reuse"
    # Drop everything v2 recomputes, so no v1 result can leak into a v2 table.
    for c in ["cell_type_manual", "cell_type", "cell_type_fine", "annotation_source",
              "annotation_conf", *[c for c in adata.obs.columns
                                   if c.startswith(("score_", "program_"))]]:
        if c in adata.obs:
            del adata.obs[c]
    REUSE = True
    print(f"reusing prior automated predictions for {len(REUSED)} columns")
else:
    adata = sc.read_h5ad(IN_H5AD)
    REUSE = False

print(adata.shape)
print(adata.obs[CLUSTER_KEY].value_counts().sort_index().to_string())
rss("load")

# %% [markdown]
# ## Method A — manual marker panel, at cluster level
#
# Cluster level, never per cell. At **1,162 median genes per cell** a per-cell marker
# call on a gene that dropped out is a wrong call, not a missing one; clustering is
# what absorbs that. A cluster with no clear winner is recorded `Ambiguous` rather than
# forced into a class.

# %%
with warnings.catch_warnings():
    warnings.simplefilter("always")
    score_cols = ann.score_marker_panel(adata)
print("marker scores:", score_cols)
rss("marker scoring")

# v3: identity is adjudicated on DETECTION FRACTIONS, not on cross-panel score_genes
# argmax. Those module scores subtract a control set drawn from each gene's own
# expression bin, so each panel carries its own baseline offset — measured here as
# 0.2036 favouring NK over T, larger than the 0.068/0.196 margins by which two T-cell
# clusters were called NK in v1/v2. The scores above stay as descriptive within-program
# quantities and no longer decide anything.
detection = ann.marker_detection_by_cluster(adata, CLUSTER_KEY)
detection.to_csv(OUT / "marker_detection_by_cluster.csv")

manual_tbl = ann.adjudicate_clusters(adata, CLUSTER_KEY, out_key="cell_type_manual")
manual_tbl.to_csv(OUT / "manual_cluster_adjudication.csv")
print(manual_tbl[["winner", "reason", "lead_positive", "runner_up_positive",
                  "supported", "survivors"]].to_string())

# %% [markdown]
# ### The DE pass is not optional
#
# It is the only step that can surface a population the seven-class panel does not
# cover at all — pDC, erythroid progenitors, a doublet-driven cluster.

# %%
sc.tl.rank_genes_groups(adata, CLUSTER_KEY, method="wilcoxon", key_added="rank_manual")
de_top = pd.DataFrame(
    {c: [g for g in adata.uns["rank_manual"]["names"][c][:10]]
     for c in adata.obs[CLUSTER_KEY].cat.categories}
).T
de_top.columns = [f"top{i+1}" for i in range(de_top.shape[1])]
de_top.to_csv(OUT / "cluster_de_top_genes.csv")
rss("wilcoxon DE")
print(de_top.iloc[:, :6].to_string())

# %% [markdown]
# ## Method B — CellTypist
#
# Run on **expression, never the Harmony embedding**, with majority voting over the
# stage-05 Leiden partition so the three methods share one partition and differ only in
# labelling.
#
# Note this is deliberately *not* the configuration used at stage 05b. There,
# `majority_voting=False` was load-bearing — the benchmark was selecting an embedding,
# so its labels had to be embedding-independent or the comparison would be circular.
# Here the embedding is already fixed, so voting over it is the right choice.

# %%
import anndata as ad                             # noqa: E402

if not REUSE:
    import celltypist                                        # noqa: E402
    from celltypist import models                            # noqa: E402
    # Everything that needs `layers["counts"]` happens here, then the layer is dropped.
    # A full `adata.copy()` here is what killed the first run: it duplicates a
    # 172,940 x 32,991 object on a 30 GB box that stage 05 already peaked ~20 GB on.
    # Build a lean object that carries only what CellTypist reads instead.
    pseudo = sc.get.aggregate(adata, by=CLUSTER_KEY, func="mean", layer="counts")
    pseudo_df = pd.DataFrame(
        pseudo.layers["mean"], index=pseudo.obs_names, columns=pseudo.var_names
    ).T
    print("cluster pseudobulk for SingleR:", pseudo_df.shape)

    ct_input = ad.AnnData(
        X=adata.layers["counts"].copy(),
        obs=pd.DataFrame(index=adata.obs_names),
        var=pd.DataFrame(index=adata.var_names),
    )
    del adata.layers["counts"]
    gc.collect()
    rss("counts layer freed")

    sc.pp.normalize_total(ct_input, target_sum=1e4)
    sc.pp.log1p(ct_input)

    models.download_models(model=["Immune_All_Low.pkl", "Immune_All_High.pkl"], force_update=False)
    # CellTypist is run in CELL CHUNKS with majority_voting=False, and the vote is applied
    # afterwards by `ann.majority_vote`. This is not a methodological change — it is
    # arithmetically the same assignment — it is a memory one. CellTypist densifies to
    # scale, and 172,940 cells x 5,951 model features is ~8 GB on top of an already
    # resident object; running it whole OOM-killed this stage twice. A regression test
    # pins chunked == unchunked.
    CHUNK = 20_000
    per_cell_labels, per_cell_conf = [], []
    for start in range(0, ct_input.n_obs, CHUNK):
        part = ct_input[start:start + CHUNK].copy()
        res = celltypist.annotate(part, model="Immune_All_Low.pkl", majority_voting=False)
        per_cell_labels.append(res.predicted_labels["predicted_labels"].to_numpy())
        per_cell_conf.append(res.probability_matrix.max(axis=1).to_numpy())
        del part, res
        gc.collect()
        print(f"  celltypist {min(start + CHUNK, ct_input.n_obs):>7,}/{ct_input.n_obs:,}", flush=True)

    adata.obs["celltypist_percell"] = np.concatenate(per_cell_labels)
    adata.obs["celltypist_conf"] = np.concatenate(per_cell_conf)
    adata.obs["celltypist_raw"] = ann.majority_vote(
        adata.obs["celltypist_percell"], adata.obs[CLUSTER_KEY].astype(str)
    ).to_numpy()
    del ct_input, per_cell_labels, per_cell_conf
    gc.collect()
    rss("celltypist")
else:
    pseudo_df = None      # SingleR is reused too; no pseudobulk needed
    if "counts" in adata.layers:
        del adata.layers["counts"]
        gc.collect()
    rss("reused celltypist")


adata.obs["cell_type_celltypist"] = ann.collapse_labels(
    adata.obs["celltypist_raw"], ann.CELLTYPIST_TO_BROAD
).to_numpy()

unmapped = sorted(set(adata.obs.loc[adata.obs["cell_type_celltypist"] == config.AMBIGUOUS_LABEL,
                                    "celltypist_raw"].astype(str)))
print("CellTypist native labels:", adata.obs["celltypist_raw"].nunique())
print("unmapped ->", unmapped or "none")
print(adata.obs["cell_type_celltypist"].value_counts().to_string())

# %% [markdown]
# ## Method C — SingleR
#
# Chosen to cover CellTypist's predictable blind spot rather than duplicate its
# strengths: `Immune_All_*` is immune-only, so erythroid and HSPC are where automated
# annotation is most likely to fail here.
#
# **Both references are run and both are reported; neither is silently substituted.**
# `NovershternHematopoieticData` is the documented primary. It has **no plasma-cell
# label at any level** — its B-lineage stops at "Mature B cells class switched" — so it
# *cannot* return PlasmaCell and will push plasma cells into B cells. That is a property
# of the reference, not of the cells, and it is exactly why the plasma-cell boundary is
# judged on marker coverage rather than on automated agreement.
# `HumanPrimaryCellAtlasData` carries `B_cell:Plasma_cell`, but only in `label.fine`.

# %%
if not REUSE:
    # rpy2 resolves R via $R_HOME, which a bare `python notebooks/06_annotation.py` does
    # not have (only an activated env or the registered kernel does). Set it from this
    # interpreter's own prefix so the stage runs identically as a script and as a notebook.
    import os                                        # noqa: E402

    os.environ.setdefault("R_HOME", str(Path(sys.executable).parent.parent / "lib" / "R"))

    import rpy2.robjects as ro                       # noqa: E402
    from rpy2.robjects import numpy2ri               # noqa: E402
    from rpy2.robjects.conversion import localconverter  # noqa: E402

    ro.r("suppressPackageStartupMessages({library(SingleR); library(celldex)})")
    with localconverter(ro.default_converter + numpy2ri.converter):
        ro.globalenv["mat"] = pseudo_df.to_numpy()
    ro.globalenv["rn"] = ro.StrVector(pseudo_df.index.tolist())
    ro.globalenv["cn"] = ro.StrVector(pseudo_df.columns.tolist())
    ro.r("rownames(mat) <- rn; colnames(mat) <- cn; mat <- log1p(mat)")

    singler_out = {}
    for refname, level in [("NovershternHematopoieticData", "label.main"),
                           ("HumanPrimaryCellAtlasData", "label.fine")]:
        try:
            ro.r(f'ref <- celldex::{refname}()')
            ro.r(f'pred <- SingleR(test=mat, ref=ref, labels=ref${level})')
            labs = list(ro.r("as.character(pred$labels)"))
            pruned = list(ro.r('as.character(ifelse(is.na(pred$pruned.labels), "PRUNED", pred$pruned.labels))'))
            singler_out[refname] = pd.DataFrame(
                {"label": labs, "pruned": pruned}, index=pseudo_df.columns
            )
            print(f"{refname} ({level}): OK, {len(set(labs))} distinct labels")
        except Exception as exc:                                    # noqa: BLE001
            singler_out[refname] = None
            print(f"{refname} ({level}): FAILED -> {exc}")

    # Primary = Novershtern; HPCA is the documented fallback and also the only SingleR route
    # to a plasma-cell call. Both are carried as separate label columns, never merged.
    for refname, short in [("NovershternHematopoieticData", "singler_nov"),
                           ("HumanPrimaryCellAtlasData", "singler_hpca")]:
        tbl = singler_out[refname]
        if tbl is None:
            print(f"{short}: not available, column not written")
            continue
        per_cell_raw = adata.obs[CLUSTER_KEY].astype(str).map(tbl["label"].to_dict())
        adata.obs[f"{short}_raw"] = per_cell_raw.to_numpy()
        adata.obs[f"cell_type_{short}"] = ann.collapse_labels(
            per_cell_raw, ann.SINGLER_TO_BROAD
        ).to_numpy()
        print(f"\n{short}:")
        print(adata.obs[f"cell_type_{short}"].value_counts().to_string())
else:
    singler_out = {"NovershternHematopoieticData": "reused from v1",
                   "HumanPrimaryCellAtlasData": "reused from v1"}
    for short in ("singler_nov", "singler_hpca"):
        print(f"{short}: reused")
        print(adata.obs[f"cell_type_{short}"].value_counts().to_string())

# %% [markdown]
# ## Marker coverage — the decisive, biological test
#
# For each method, mean scaled expression of a class's own markers within the cells
# that method assigned to that class. Scaled across labels (the
# `dotplot(standard_scale="var")` transform), so a gene that is middling everywhere
# scores ~0 and only genuine enrichment scores high.

# %%
METHODS = {"manual": "cell_type_manual", "celltypist": "cell_type_celltypist"}
if "cell_type_singler_nov" in adata.obs:
    METHODS["singler_nov"] = "cell_type_singler_nov"
if "cell_type_singler_hpca" in adata.obs:
    METHODS["singler_hpca"] = "cell_type_singler_hpca"

coverage = {m: ann.marker_coverage(adata, k) for m, k in METHODS.items()}
cov_wide = pd.DataFrame({m: t["coverage"] for m, t in coverage.items()})
cov_wide.to_csv(OUT / "marker_coverage_by_method.csv")
print(cov_wide.round(3).to_string())
print(f"\nveto line: {config.MARKER_COVERAGE_MIN}")
rss("marker coverage")

# %% [markdown]
# ### Lineage exclusivity — the second veto (v2)
#
# Coverage asks whether cells labelled X express X's markers, which is
# precision-like and therefore blind to a class that has *swallowed* another lineage.
# This asks the complementary question: do those same cells carry strong positive
# evidence for a lineage incompatible with X? Detection-based, so dropout can only
# hide a contradiction and never invent one.

# %%
contradiction = {m: ann.contradiction_rate(adata, k) for m, k in METHODS.items()}
contra_wide = pd.DataFrame({m: t["contradiction_rate"] for m, t in contradiction.items()})
contra_wide.to_csv(OUT / "lineage_contradiction_by_method.csv")
for m, t in contradiction.items():
    t.to_csv(OUT / f"lineage_contradiction_{m}.csv")
print(contra_wide.round(3).to_string())
print(f"\nveto line: {config.CONTRADICTION_MAX_RATE}")
rss("lineage exclusivity")

# %% [markdown]
# ### Coverage by cohort, by diagnosis, and by stage-05 plasma cluster
#
# This is the check the stage turns on. Good annotation on donor marrow plus systematic
# under-recognition of plasma cells in one censored myeloma cohort would pass every
# global number here and still break stage 07.

# %%
strata = []
for label, col in [("cohort", "cohort"), ("diagnosis", "sample_type")]:
    for val in sorted(adata.obs[col].astype(str).unique()):
        sub = adata[adata.obs[col].astype(str) == val]
        for m, k in METHODS.items():
            c = ann.marker_coverage(sub, k)["coverage"]
            strata.append(pd.Series(c, name=(label, val, m)))
cov_strata = pd.DataFrame(strata)
cov_strata.index = pd.MultiIndex.from_tuples(cov_strata.index, names=["stratum", "value", "method"])
cov_strata.to_csv(OUT / "marker_coverage_by_stratum.csv")
print(cov_strata["PlasmaCell"].round(3).to_string())

# %%
# The three cohort-shaped plasma clusters from stage 05, individually.
prof = pd.read_csv(REPO / "results" / "05_integration" / "cluster_profile.csv")
plasma_clusters = prof.loc[prof["plasma_like"], "leiden"].astype(str).tolist()
rows = []
for cl in plasma_clusters:
    sub = adata[adata.obs[CLUSTER_KEY].astype(str) == cl]
    pcs = config.MARKER_PANEL["PlasmaCell"]
    present = [g for g in pcs if g in sub.var_names]
    X = sub[:, present].X
    X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
    rows.append({
        "leiden": cl, "n_cells": sub.n_obs,
        "dominant_cohort": sub.obs["cohort"].astype(str).mode().iat[0],
        "median_umi": float(np.median(sub.obs["total_counts"])),
        **{f"pct_detect_{g}": float((X[:, i] > 0).mean() * 100) for i, g in enumerate(present)},
        "manual": sub.obs["cell_type_manual"].astype(str).mode().iat[0],
        "celltypist": sub.obs["cell_type_celltypist"].astype(str).mode().iat[0],
    })
plasma_tbl = pd.DataFrame(rows).set_index("leiden")
plasma_tbl.to_csv(OUT / "plasma_cluster_marker_support.csv")
print(plasma_tbl.to_string())

# %% [markdown]
# ## Concordance — agreement, not accuracy
#
# Labelled as concordance everywhere. It is the weakest of the three evidence types and
# is presented first for that reason, not because it is the strongest.

# %%
conc = {
    m: ann.per_class_concordance(adata.obs["cell_type_manual"], adata.obs[k])
    for m, k in METHODS.items() if m != "manual"
}
f1_wide = pd.DataFrame({m: t["f1"] for m, t in conc.items()})
f1_wide["threshold"] = pd.Series(config.CONCORDANCE_THRESHOLDS)
f1_wide.to_csv(OUT / "concordance_f1.csv")
for m, t in conc.items():
    t.to_csv(OUT / f"concordance_{m}.csv")
print(f1_wide.round(3).to_string())

# %% [markdown]
# ## The pre-declared decision
#
# Coverage veto first, then the concordance bar, then highest F1, then fallback to
# manual — and if manual also fails the veto, the class is unresolved and says so.
# Thresholds are **not** revisited now that results exist.

# %%
decision = ann.decide_per_class(conc, coverage, contradiction)
decision.to_csv(OUT / "annotation_decision.csv")
print(decision[["chosen_method", "reason", "f1_chosen", "coverage_chosen",
                "contradiction_chosen", "vetoed_by_coverage",
                "vetoed_by_contradiction", "not_evaluable"]].to_string())

# Explicit NOT_EVALUABLE audit across all seven classes — no combination may be
# silently read as a pass or a veto.
state_cols = [c for c in decision.columns if c.startswith(("cov_state_", "contra_state_"))]
states = decision[state_cols]
states.to_csv(OUT / "veto_states.csv")
print("\n--- NOT_EVALUABLE combinations (all seven classes) ---")
ne = [(cls, c) for cls in decision.index for c in state_cols
      if decision.loc[cls, c] == ann.NOT_EVALUABLE]
for cls, c in ne:
    print(f"  {cls:11s} {c}")
print(f"  total: {len(ne)}")

# %%
summary = ann.assemble_final_labels(
    adata, decision, METHODS,
    conf_keys={"celltypist": "celltypist_conf"},
)
summary.to_csv(OUT / "final_cell_type_counts.csv")
print(summary.round(2).to_string())

# %% [markdown]
# ## State programs — continuous, non-exclusive, never identity
#
# A cycling plasma cell is PlasmaCell **plus** a high cell-cycle score, not a "Cycling"
# cell type. These are covariates for stages 10-12 and never enter `cell_type`.

# %%
with warnings.catch_warnings():
    warnings.simplefilter("always")
    prog_cols = ann.score_state_programs(adata)
print("programs:", prog_cols)

prog_summary = adata.obs.groupby("cell_type", observed=True)[prog_cols].mean()
prog_summary.to_csv(OUT / "state_program_means_by_cell_type.csv")
adata.obs[["cell_type", *prog_cols]].to_csv(OUT / "state_program_scores.csv.gz",
                                            compression="gzip")
print(prog_summary.round(3).to_string())

# %% [markdown]
# ## Outputs

# %%
per_cell = adata.obs[
    [c for c in ["sample_name", "patient_id", "cohort", "sample_type", CLUSTER_KEY,
                 "cell_type_manual", "cell_type_celltypist", "cell_type_singler_nov",
                 "cell_type_singler_hpca", "cell_type", "annotation_source",
                 "annotation_conf"] if c in adata.obs]
]
per_cell.to_csv(OUT / "per_cell_labels.csv.gz", compression="gzip")

for key, name in [("patient_id", "composition_by_patient"), ("sample_name", "composition_by_sample")]:
    ann.composition_by(adata, key).to_csv(OUT / f"{name}.csv")

adata.obs["cell_type_fine"] = adata.obs["celltypist_raw"].astype(str)
adata.write_h5ad(OUT / "annotated.h5ad", compression="gzip")

meta = {
    "n_cells": int(adata.n_obs),
    "methods": list(METHODS),
    "singler_references": {k: ("OK" if v is not None else "FAILED") for k, v in singler_out.items()},
    "coverage_min": config.MARKER_COVERAGE_MIN,
    "f1_thresholds": config.CONCORDANCE_THRESHOLDS,
    "decision": decision["chosen_method"].to_dict(),
    "contradiction_max": config.CONTRADICTION_MAX_RATE,
    "contradiction_min_genes": config.CONTRADICTION_MIN_GENES,
    "reused_v2_predictions": bool(REUSE),
    "manual_adjudication": {"detect_min": config.MANUAL_MARKER_DETECT_MIN,
                            "positive_min": config.MANUAL_POSITIVE_MIN,
                            "margin": config.MANUAL_DECISION_MARGIN},
}
(OUT / "run_metadata.json").write_text(json.dumps(meta, indent=2))
print(json.dumps(meta, indent=2))
