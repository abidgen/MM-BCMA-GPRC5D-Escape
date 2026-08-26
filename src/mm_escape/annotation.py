"""
Stage 06 — cell type annotation.

WHAT THIS STAGE PRODUCES
------------------------
`obs["cell_type"]` — seven coarse classes plus `Ambiguous` — and nothing else that is
load-bearing. Stages 07-12 read that column and never branch on annotation logic. The
provenance columns (`annotation_source`, `annotation_conf`) exist so a mixed-provenance
label is traceable; `cell_type_fine` is stage 11's convenience and is never load-bearing.

THREE METHODS, COMPARED, CHOSEN PER CLASS
-----------------------------------------
Manual marker scoring, CellTypist and SingleR are run separately and the winner is
picked **per class** against bars declared before any result was looked at
(`config.CONCORDANCE_THRESHOLDS`, `config.MARKER_COVERAGE_MIN`). The methods are
expected to fail on *different* classes, so a single verdict for the whole stage would
throw away good labels to punish an unrelated weakness.

CONCORDANCE IS NOT ACCURACY, AND MARKER EVIDENCE OUTRANKS IT
------------------------------------------------------------
The manual labels are a third opinion derived from the same expression matrix, not
ground truth, so F1-against-manual measures *agreement*. Two automated methods agreeing
with each other is the strongest agreement available — and still only agreement, between
references that share canonical marker biology and may share its blind spots.

So marker coverage is a **veto**: a class whose assigned cells do not actually express
that class's markers is rejected regardless of how well the methods agree. High
agreement on a biologically unsupported label is agreement on an error. `decide_per_class`
enforces this ordering and there is no way to pass a class on concordance alone.

A REFERENCE CAN BE STRUCTURALLY INCAPABLE OF A CLASS
-----------------------------------------------------
`NovershternHematopoieticData` has no plasma-cell label at any level — its B-lineage
stops at "Mature B cells class switched". SingleR against it therefore *cannot* return
PlasmaCell, and will push plasma cells into B cells. That is a property of the
reference, not a failure of the cell. `HumanPrimaryCellAtlasData` carries
`B_cell:Plasma_cell`, but only in `label.fine`. Both are run and both are reported;
neither is silently substituted for the other.

IDENTITY AND STATE ARE DIFFERENT AXES
--------------------------------------
`score_state_programs` writes continuous float columns and never touches `cell_type`.
A cycling plasma cell is PlasmaCell **plus** a high cell-cycle score, not a "Cycling"
cell type. Collapsing state into identity discards exactly the intermediate cells
stage 10 needs.
"""

from __future__ import annotations

import re
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from . import config

__all__ = [
    "BROAD_CLASSES",
    "CELLTYPIST_TO_BROAD",
    "SINGLER_TO_BROAD",
    "collapse_labels",
    "majority_vote",
    "score_marker_panel",
    "manual_labels_from_clusters",
    "marker_detection_by_cluster",
    "adjudicate_clusters",
    "class_support",
    "PASS",
    "FAIL",
    "NOT_EVALUABLE",
    "evaluate_veto",
    "marker_coverage",
    "lineage_evidence",
    "contradiction_rate",
    "contradiction_concentration",
    "cytotoxic_lineage_calls",
    "per_class_concordance",
    "decide_per_class",
    "assemble_final_labels",
    "score_state_programs",
    "composition_by",
]

#: The seven project classes, in `config.MARKER_PANEL` order. `Ambiguous` is not one of
#: them — it is the explicit "no class won" outcome, counted and reported, never a class
#: a method can be judged on.
BROAD_CLASSES: tuple[str, ...] = tuple(config.MARKER_PANEL)

AMBIGUOUS = config.AMBIGUOUS_LABEL

#: CellTypist `Immune_All_*` labels -> project classes.
#:
#: `ILC` maps to NK deliberately and this was checked rather than assumed at stage 05b:
#: in this cohort the ILC calls are 98.7% NKG7+, 92.2% GNLY+, 85.9% KLRD1+ and 1.3%
#: MS4A1+. They are NK cells.
#:
#: Anything not listed collapses to `Ambiguous` rather than being guessed at — that is
#: what `collapse_labels` does with an unmapped label, and the unmapped rate is a
#: reported quantity.
CELLTYPIST_TO_BROAD: dict[str, str] = {
    # plasma / B lineage
    "Plasma cells": "PlasmaCell",
    "Plasmablasts": "PlasmaCell",
    "B cells": "Bcell",
    "Memory B cells": "Bcell",
    "Naive B cells": "Bcell",
    "Follicular B cells": "Bcell",
    "Germinal center B cells": "Bcell",
    "Large pre-B cells": "Bcell",
    "Small pre-B cells": "Bcell",
    "Pro-B cells": "Bcell",
    "Proliferative germinal center B cells": "Bcell",
    # T
    "T cells": "Tcell",
    "Tcm/Naive helper T cells": "Tcell",
    "Tcm/Naive cytotoxic T cells": "Tcell",
    "Tem/Effector helper T cells": "Tcell",
    "Tem/Temra cytotoxic T cells": "Tcell",
    "Tem/Trm cytotoxic T cells": "Tcell",
    "Regulatory T cells": "Tcell",
    "Type 1 helper T cells": "Tcell",
    "Type 17 helper T cells": "Tcell",
    "Follicular helper T cells": "Tcell",
    "Cycling T cells": "Tcell",
    "Double-negative thymocytes": "Tcell",
    "Double-positive thymocytes": "Tcell",
    "MAIT cells": "Tcell",
    "gamma-delta T cells": "Tcell",
    # NK / ILC
    "NK cells": "NK",
    "ILC": "NK",
    "ILC1": "NK",
    "ILC2": "NK",
    "ILC3": "NK",
    "CD16+ NK cells": "NK",
    "CD16- NK cells": "NK",
    "Cycling NK cells": "NK",
    # myeloid
    "Monocytes": "Myeloid",
    "Classical monocytes": "Myeloid",
    "Non-classical monocytes": "Myeloid",
    "Intermediate macrophages": "Myeloid",
    "Macrophages": "Myeloid",
    "DC": "Myeloid",
    "DC1": "Myeloid",
    "DC2": "Myeloid",
    "pDC": "Myeloid",
    "pDC precursor": "Myeloid",
    "Migratory DC": "Myeloid",
    "Mast cells": "Myeloid",
    "Neutrophils": "Myeloid",
    "Granulocytes": "Myeloid",
    "Myelocytes": "Myeloid",
    "Promyelocytes": "Myeloid",
    "Monocyte precursor": "Myeloid",
    "Cycling monocytes": "Myeloid",
    "Cycling DCs": "Myeloid",
    # erythroid / megakaryocyte
    "Erythroid": "Erythroid",
    "Erythrocytes": "Erythroid",
    "Early erythroid": "Erythroid",
    "Mid erythroid": "Erythroid",
    "Late erythroid": "Erythroid",
    "Megakaryocyte precursor": "Erythroid",
    "Megakaryocytes/platelets": "Erythroid",
    "MEMP": "Erythroid",
    # progenitors
    "HSC/MPP": "HSPC",
    "HSC": "HSPC",
    "MPP": "HSPC",
    "CMP": "HSPC",
    "GMP": "HSPC",
    "MEP": "HSPC",
    "ELP": "HSPC",
    "CLP": "HSPC",
    "Early MK": "HSPC",
    "Cycling HSC/MPP": "HSPC",
}

#: SingleR labels -> project classes, covering **both** references.
#:
#: `NovershternHematopoieticData` (label.main) and `HumanPrimaryCellAtlasData`
#: (label.main and label.fine) use different vocabularies; one map serves both because
#: the keys do not collide.
#:
#: Note what is deliberately ABSENT: Novershtern has no plasma-cell key, because the
#: reference has no such class. Its plasma cells arrive labelled "B cells" and will map
#: to Bcell — which is why PlasmaCell must be judged on marker coverage, not on
#: SingleR/Novershtern agreement.
SINGLER_TO_BROAD: dict[str, str] = {
    # --- Novershtern label.main ---
    "B cells": "Bcell",
    "CD4+ T cells": "Tcell",
    "CD8+ T cells": "Tcell",
    "NK T cells": "Tcell",
    "NK cells": "NK",
    "Monocytes": "Myeloid",
    "Dendritic cells": "Myeloid",
    "Granulocytes": "Myeloid",
    "Basophils": "Myeloid",
    "Eosinophils": "Myeloid",
    "Erythroid cells": "Erythroid",
    "Megakaryocytes": "Erythroid",
    "HSCs": "HSPC",
    "CMPs": "HSPC",
    "GMPs": "HSPC",
    "MEPs": "HSPC",
    # --- HPCA label.main ---
    "B_cell": "Bcell",
    "T_cells": "Tcell",
    "NK_cell": "NK",
    "Monocyte": "Myeloid",
    "Macrophage": "Myeloid",
    "DC": "Myeloid",
    "Neutrophils": "Myeloid",
    "Myelocyte": "Myeloid",
    "Pro-Myelocyte": "Myeloid",
    "Platelets": "Erythroid",
    "Erythroblast": "Erythroid",
    "HSC_CD34+": "HSPC",
    "HSC_-G-CSF": "HSPC",
    "CMP": "HSPC",
    "GMP": "HSPC",
    "MEP": "HSPC",
    "Pre-B_cell_CD34-": "Bcell",
    "Pro-B_cell_CD34+": "Bcell",
    "BM": AMBIGUOUS,
    "BM & Prog.": AMBIGUOUS,
    # --- HPCA label.fine: the only route to a plasma-cell call from SingleR ---
    "B_cell:Plasma_cell": "PlasmaCell",
    "B_cell:Naive": "Bcell",
    "B_cell:Memory": "Bcell",
    "B_cell:immature": "Bcell",
    "B_cell:Germinal_center": "Bcell",
    "B_cell:CXCR4-_centrocyte": "Bcell",
    "B_cell:CXCR4+_centroblast": "Bcell",
}


def collapse_labels(
    labels: Iterable[str],
    mapping: Mapping[str, str],
    *,
    unmapped: str = AMBIGUOUS,
) -> pd.Series:
    """Collapse a method's native labels onto the seven project classes.

    An unmapped label becomes `unmapped` (default `Ambiguous`) rather than being
    guessed at or dropped. The unmapped rate is a reported diagnostic — a method whose
    labels largely fail to map is telling you something, and silently coercing them
    would hide it.

    HPCA `label.fine` keys like ``B_cell:Plasma_cell`` are matched exactly first; a
    fine label with no exact key falls back to its ``label.main`` prefix (the part
    before ``:``), so unlisted fine variants still land on the right broad class.
    """
    out = []
    for raw in labels:
        if raw is None or (isinstance(raw, float) and np.isnan(raw)):
            out.append(unmapped)
            continue
        key = str(raw).strip()
        if key in mapping:
            out.append(mapping[key])
        elif ":" in key and key.split(":", 1)[0] in mapping:
            out.append(mapping[key.split(":", 1)[0]])
        else:
            out.append(unmapped)
    return pd.Series(out, dtype=object)


def majority_vote(
    labels: Iterable[str],
    clusters: Iterable[str],
    *,
    min_prop: float = 0.0,
    undecided: str = "Heterogeneous",
) -> pd.Series:
    """Assign every cell in an over-cluster its cluster's most frequent label.

    This is what `celltypist.annotate(majority_voting=True)` does internally, done
    explicitly here for one reason: CellTypist **densifies** to scale, and at
    172,940 cells x 5,951 model features that is ~8 GB on top of an already-resident
    object, which OOM-killed the run twice. Predicting in cell chunks with
    `majority_voting=False` and voting afterwards is arithmetically the same and fits
    in memory.

    `min_prop` is CellTypist's own threshold: a cluster whose winning label does not
    reach that share of its cells is `undecided` rather than being assigned a plurality
    label. Default 0.0 matches CellTypist's default (plain mode).
    """
    df = pd.DataFrame({"label": list(map(str, labels)), "cluster": list(map(str, clusters))})
    if len(df) == 0:
        return pd.Series([], dtype=object)

    winners = {}
    for cl, grp in df.groupby("cluster", observed=True):
        counts = grp["label"].value_counts()
        top, n = counts.index[0], counts.iloc[0]
        winners[cl] = top if (n / len(grp)) >= min_prop else undecided
    return df["cluster"].map(winners).astype(object)


def score_marker_panel(
    adata,
    panel: Mapping[str, Sequence[str]] | None = None,
    *,
    prefix: str = "score_",
) -> list[str]:
    """Score each class's marker panel with `scanpy.tl.score_genes`.

    Writes one float column per class and returns the column names. Genes absent from
    the harmonized space are skipped with a warning rather than failing the run —
    unlike `REQUIRED_GENES`, these panels are indicative, not a contract.
    """
    import scanpy as sc
    import warnings

    panel = config.MARKER_PANEL if panel is None else panel
    written = []
    for cls, genes in panel.items():
        present = [g for g in genes if g in adata.var_names]
        missing = sorted(set(genes) - set(present))
        if missing:
            warnings.warn(f"{cls}: markers absent from gene space: {missing}", stacklevel=2)
        if not present:
            warnings.warn(f"{cls}: no markers present, score not written", stacklevel=2)
            continue
        col = f"{prefix}{cls}"
        sc.tl.score_genes(adata, present, score_name=col)
        written.append(col)
    return written


def manual_labels_from_clusters(
    adata,
    cluster_key: str,
    *,
    prefix: str = "score_",
    margin: float = 0.0,
    out_key: str = "cell_type_manual",
) -> pd.DataFrame:
    """Assign one label per *cluster* from mean marker scores, then broadcast to cells.

    Cluster level, never per cell: at 1,162 median genes per cell a per-cell marker
    call on a dropped-out gene is a wrong call, not a missing one, and clustering is
    what absorbs that.

    A cluster whose best and second-best class means differ by less than `margin` is
    recorded as `Ambiguous` rather than forced into the winner. Returns the per-cluster
    table (means, winner, runner-up, margin) as assignment evidence.
    """
    cols = [c for c in adata.obs.columns if c.startswith(prefix)]
    if not cols:
        raise ValueError(f"no marker-score columns with prefix {prefix!r}; run score_marker_panel first")

    means = adata.obs.groupby(cluster_key, observed=True)[cols].mean()
    classes = [c[len(prefix):] for c in cols]
    means.columns = classes

    order = np.argsort(-means.to_numpy(), axis=1)
    best = [classes[i] for i in order[:, 0]]
    second = [classes[i] for i in order[:, 1]] if len(classes) > 1 else [None] * len(means)
    top = np.take_along_axis(means.to_numpy(), order[:, :1], axis=1).ravel()
    runner = (
        np.take_along_axis(means.to_numpy(), order[:, 1:2], axis=1).ravel()
        if len(classes) > 1
        else np.full(len(means), -np.inf)
    )
    gap = top - runner
    winner = [b if g >= margin else AMBIGUOUS for b, g in zip(best, gap)]

    table = pd.DataFrame(
        {"winner": winner, "best": best, "second": second, "top_score": top,
         "runner_up_score": runner, "margin": gap},
        index=means.index,
    ).join(means.add_prefix("mean_"))

    adata.obs[out_key] = (
        adata.obs[cluster_key].map(dict(zip(table.index, table["winner"]))).astype("category")
    )
    return table


#: Three-state veto outcome. NaN is a THIRD state, never coerced to a boolean.
#:
#: v3 shipped a real bug from conflating them: `decide_per_class` read a NaN
#: contradiction as "pass" and a NaN coverage as "fail", so NK's decision reported
#: "fallback also failed the coverage veto" when the truth was that the manual
#: classifier produced no NK population at all and there was nothing to evaluate.
#: A class that cannot be assessed has not passed and has not failed, and the
#: distinction changes what you do about it.
PASS = "PASS"
FAIL = "FAIL"
NOT_EVALUABLE = "NOT_EVALUABLE"


def evaluate_veto(value: float, threshold: float, *, direction: str) -> str:
    """Three-state veto evaluation. `direction` is "min" (>= passes) or "max" (<= passes).

    A non-finite value is `NOT_EVALUABLE` — it can neither make a method a winner nor
    record it as vetoed.
    """
    if direction not in ("min", "max"):
        raise ValueError(f"direction must be 'min' or 'max', got {direction!r}")
    if value is None or not np.isfinite(value):
        return NOT_EVALUABLE
    if direction == "min":
        return PASS if value >= threshold else FAIL
    return PASS if value <= threshold else FAIL


def marker_detection_by_cluster(
    adata,
    cluster_key: str,
    panel: Mapping[str, Sequence[str]] | None = None,
) -> pd.DataFrame:
    """Detection fraction of every panel gene, per cluster.

    Detection fractions live on a common [0, 1] scale with no per-panel normalisation,
    which is precisely what `scanpy.tl.score_genes` does not offer and why its outputs
    cannot be compared across panels. This is the evidence `adjudicate_clusters` reads.
    """
    panel = config.MARKER_PANEL if panel is None else panel
    # Must cover every gene any support rule reads, not just MARKER_PANEL — the
    # Myeloid subprograms, the HSPC core and the mature-plasma predicate all consult
    # this table, and a gene missing from it silently reads as "not present" rather
    # than "not detected".
    genes = [g for gs in panel.values() for g in gs]
    genes += [g for gs in config.MYELOID_SUBPROGRAMS.values() for g in gs]
    genes += list(config.MYELOID_MONO_ANCHORS) + list(config.MYELOID_MONO_CONTEXT)
    genes += list(config.MYELOID_DC_ANCHORS) + list(config.MYELOID_DC_CONTEXT)
    genes += list(config.MYELOID_PDC_CORE)
    genes += list(config.MYELOID_SUPPORTING) + list(config.PDC_CONTEXTUAL)
    genes += list(config.HSPC_CORE)
    genes += [g for gs in config.HSPC_CONTEXT.values() for g in gs]
    genes += list(config.PLASMA_SECRETORY) + list(config.PLASMA_MATURE)
    genes = [g for g in dict.fromkeys(genes) if g in adata.var_names]
    X = adata[:, genes].X
    X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
    det = pd.DataFrame(X > 0, index=adata.obs_names, columns=genes)
    clusters = pd.Series(np.asarray(adata.obs[cluster_key]), index=adata.obs_names).astype(str)
    out = det.groupby(clusters, observed=True).mean()
    out.index.name = cluster_key
    return out


def class_support(det_row: pd.Series, *, detect_min: float, positive_min: float) -> dict:
    """Positive support for each class in one cluster, from its own biology.

    Not one rule for seven classes. `Myeloid` is an umbrella over distinct programs and
    is supported by ANY coherent subprogram; `HSPC` is primitive identity kept separate
    from lineage priming; `PlasmaCell` needs concordant evidence on two independent
    axes. The rest use their `MARKER_PANEL` fraction as before.
    """
    def frac(genes) -> float:
        present = [g for g in genes if g in det_row.index]
        return float((det_row[present] >= detect_min).mean()) if present else 0.0

    def hits(genes) -> list[str]:
        return [g for g in genes if g in det_row.index and det_row[g] >= detect_min]

    out: dict[str, dict] = {}

    # Three independent routes. Monocyte is a plain program; conventional DC needs a
    # myeloid/DC-restricted ANCHOR *and* APC context, so shared MHC-II biology can
    # corroborate but never identify; pDC has its own core and never passes through
    # the MHC-II route.
    sub = {name: frac(genes) for name, genes in config.MYELOID_SUBPROGRAMS.items()}
    mono_anchor = frac(config.MYELOID_MONO_ANCHORS)
    mono_context = frac(config.MYELOID_MONO_CONTEXT)
    sub["monocyte"] = min(mono_anchor, mono_context)
    dc_anchor = frac(config.MYELOID_DC_ANCHORS)
    dc_context = frac(config.MYELOID_DC_CONTEXT)
    sub["dc"] = min(dc_anchor, dc_context)      # a route is only as strong as its weaker axis
    sub["pdc"] = frac(config.MYELOID_PDC_CORE)
    passing = [n for n, f in sub.items() if f >= positive_min]
    out["Myeloid"] = {"fraction": max(sub.values()) if sub else 0.0, "subprograms": sub,
                      "subprograms_passing": passing,
                      "supporting": hits(config.MYELOID_SUPPORTING),
                      "mono_anchor": mono_anchor, "mono_context": mono_context,
                      "mono_anchor_hits": hits(config.MYELOID_MONO_ANCHORS),
                      "dc_anchor": dc_anchor, "dc_context": dc_context,
                      "dc_anchor_hits": hits(config.MYELOID_DC_ANCHORS),
                      "pdc_contextual": hits(config.PDC_CONTEXTUAL),
                      "supported": bool(passing)}

    core = frac(config.HSPC_CORE)
    out["HSPC"] = {"fraction": core, "core_hits": hits(config.HSPC_CORE),
                   "priming": {k: hits(v) for k, v in config.HSPC_CONTEXT.items()},
                   "supported": core >= positive_min}

    sec_hits = hits(config.PLASMA_SECRETORY)
    mat_hits = hits(config.PLASMA_MATURE)
    n_sec = len([g for g in config.PLASMA_SECRETORY if g in det_row.index])
    secretory_ok = n_sec > 0 and len(sec_hits) == n_sec
    mature_ok = len(mat_hits) >= 1
    out["PlasmaCell"] = {"fraction": frac(config.MARKER_PANEL["PlasmaCell"]),
                         "secretory_hits": sec_hits, "mature_hits": mat_hits,
                         "secretory_ok": secretory_ok, "mature_ok": mature_ok,
                         "supported": bool(secretory_ok and mature_ok)}

    for cls, genes in config.MARKER_PANEL.items():
        if cls in out:
            continue
        f = frac(genes)
        out[cls] = {"fraction": f, "hits": hits(genes), "supported": f >= positive_min}
    return out


def adjudicate_clusters(
    adata,
    cluster_key: str,
    *,
    panel: Mapping[str, Sequence[str]] | None = None,
    programs: Mapping[str, Sequence[str]] | None = None,
    pairs: Mapping[str, Sequence[str]] | None = None,
    detect_min: float | None = None,
    positive_min: float | None = None,
    margin: float | None = None,
    contradiction_max: float | None = None,
    min_genes: int | None = None,
    out_key: str = "cell_type_manual",
) -> pd.DataFrame:
    """Assign one label per cluster from positive evidence plus exclusion evidence.

    Replaces the v1/v2 cross-panel `score_genes` argmax, which was not calibrated for
    comparison across panels normalised against different control sets. The rule is:

        manual annotation = positive lineage evidence + specificity/exclusion evidence

    applied in this order:

    1. **Support.** A class is supported when at least `positive_min` of its
       `MARKER_PANEL` genes are detected in at least `detect_min` of the cluster's
       cells. A majority of the panel, because single markers are shared between
       lineages — `NKG7`/`GNLY` are a cytotoxic-granule program, not an NK identity.
    2. **Exclusion.** A supported class is disqualified when the cluster's
       contradiction rate against that class's incompatible lineages exceeds
       `contradiction_max`. This is what stops a cytotoxic T-cell cluster being called
       NK on shared effector genes: NK is contradicted by T-lineage evidence, and a
       cluster that is 87% CD3D-positive carries it overwhelmingly.
    3. **Adjudication.** Exactly one survivor wins. None → `Ambiguous`. Several → the
       leader must beat the runner-up by `margin` on positive evidence, else
       `Ambiguous`.

    `Ambiguous` is a real outcome here, not a failure path: a cluster with genuinely
    balanced evidence for two lineages is unresolved at this clustering resolution, and
    saying so is more useful than a tie-break that manufactures a decision.
    """
    panel = config.MARKER_PANEL if panel is None else panel
    pairs = config.CONTRADICTION_PAIRS if pairs is None else pairs
    detect_min = config.MANUAL_MARKER_DETECT_MIN if detect_min is None else detect_min
    positive_min = config.MANUAL_POSITIVE_MIN if positive_min is None else positive_min
    margin = config.MANUAL_DECISION_MARGIN if margin is None else margin
    contradiction_max = (
        config.CONTRADICTION_MAX_RATE if contradiction_max is None else contradiction_max
    )

    det = marker_detection_by_cluster(adata, cluster_key, panel)
    ev = lineage_evidence(adata, programs, min_genes=min_genes)
    clusters = pd.Series(np.asarray(adata.obs[cluster_key]), index=adata.obs_names).astype(str)
    contra_by_cluster = ev.groupby(clusters, observed=True).mean()

    rows = {}
    for cl in det.index:
        ev = class_support(det.loc[cl], detect_min=detect_min, positive_min=positive_min)
        positive = {c: ev[c]["fraction"] for c in panel}
        excluded = {}
        for cls in panel:
            worst = 0.0
            for lin in pairs.get(cls, ()):
                if lin in contra_by_cluster.columns:
                    worst = max(worst, float(contra_by_cluster.loc[cl, lin]))
            excluded[cls] = worst

        # The strongest SUPPORTED hypothesis is adjudicated. If it is disqualified by
        # contradictory lineage evidence the cluster is UNRESOLVED — the rule does not
        # cascade to a weaker runner-up, even one that would independently pass its own
        # exclusion check.
        #
        #     disqualification of the strongest supported biological hypothesis
        #         = unresolved evidence
        #     NOT = discard it and promote the next available label
        #
        # v3 cascaded, and Leiden 23 — NK evidence 0.80, T 0.50, NK disqualified at
        # 0.528 — was absorbed into Tcell despite NK being the better-supported
        # hypothesis. That is how a mixed NK/gamma-delta-T population disappears into a
        # clean-looking label.
        supported = [c for c in panel if ev[c]["supported"]]
        survivors = [c for c in supported if excluded[c] <= contradiction_max]

        if not supported:
            winner, reason = AMBIGUOUS, "no class supported"
            lead = second = float("nan")
        else:
            ranked = sorted(supported, key=lambda c: positive[c], reverse=True)
            lead = positive[ranked[0]]
            second = positive[ranked[1]] if len(ranked) > 1 else 0.0
            if lead - second < margin:
                winner, reason = AMBIGUOUS, (
                    f"competing lineages within margin ({lead - second:.2f} < {margin})"
                )
            elif excluded[ranked[0]] > contradiction_max:
                winner, reason = AMBIGUOUS, (
                    f"strongest supported hypothesis ({ranked[0]}) disqualified by "
                    f"contradictory lineage evidence ({excluded[ranked[0]]:.3f} > "
                    f"{contradiction_max}); not cascaded to runner-up"
                )
            else:
                winner, reason = ranked[0], "positive evidence + exclusion"

        rows[cl] = {
            "winner": winner, "reason": reason,
            "lead_positive": lead, "runner_up_positive": second,
            **{f"positive_{c}": positive[c] for c in panel},
            **{f"excluded_{c}": excluded[c] for c in panel},
            "supported": ",".join(sorted(supported)),
            "myeloid_subprograms_passing": ",".join(ev["Myeloid"]["subprograms_passing"]),
            "myeloid_sub_monocyte": ev["Myeloid"]["subprograms"].get("monocyte", 0.0),
            "myeloid_sub_dc": ev["Myeloid"]["subprograms"].get("dc", 0.0),
            "myeloid_sub_pdc": ev["Myeloid"]["subprograms"].get("pdc", 0.0),
            "hspc_core_hits": ",".join(ev["HSPC"]["core_hits"]),
            "plasma_secretory_ok": ev["PlasmaCell"]["secretory_ok"],
            "plasma_mature_ok": ev["PlasmaCell"]["mature_ok"],
            "plasma_mature_hits": ",".join(ev["PlasmaCell"]["mature_hits"]),
            "survivors": ",".join(sorted(survivors)),
            "strongest_supported": (sorted(supported, key=lambda c: positive[c],
                                           reverse=True)[0] if supported else ""),
        }

    table = pd.DataFrame(rows).T
    table.index.name = cluster_key
    adata.obs[out_key] = (
        clusters.map(table["winner"].to_dict()).astype("category")
    )
    return table


def marker_coverage(
    adata,
    label_key: str,
    panel: Mapping[str, Sequence[str]] | None = None,
    *,
    groups: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Mean scaled expression of each class's own markers within cells labelled that class.

    This is the biological-validity test, and the one that can veto a class. For each
    marker gene the per-label mean expression is min-max scaled across labels — the same
    transform `scanpy.pl.dotplot(standard_scale="var")` shows — so the result is on
    [0, 1] and comparable to `config.MARKER_COVERAGE_MIN`.

    Scaling *across labels* is what makes this a test rather than a restatement: a gene
    that is uniformly middling everywhere scores ~0 for every class, and only a gene
    genuinely enriched in its own class scores high there.

    Returns one row per class present in `label_key`, with the per-marker scaled values
    and the mean over them (`coverage`). Classes absent from the data get NaN coverage
    and are reported as not evaluable rather than as failures.
    """
    panel = config.MARKER_PANEL if panel is None else panel
    labels = pd.Series(np.asarray(adata.obs[label_key]), index=adata.obs_names).astype(str)
    if groups is None:
        groups = sorted(set(labels) - {AMBIGUOUS})

    genes = [g for gs in panel.values() for g in gs]
    genes = [g for g in dict.fromkeys(genes) if g in adata.var_names]
    X = adata[:, genes].X
    X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
    expr = pd.DataFrame(X, index=adata.obs_names, columns=genes)

    per_group = expr.groupby(labels, observed=True).mean().reindex(groups)
    lo, hi = per_group.min(axis=0), per_group.max(axis=0)
    scaled = (per_group - lo) / (hi - lo).replace(0, np.nan)

    rows = {}
    for cls, markers in panel.items():
        present = [g for g in markers if g in scaled.columns]
        if cls not in scaled.index or not present:
            rows[cls] = {"coverage": np.nan, "n_markers": len(present),
                         "n_cells": int((labels == cls).sum())}
            continue
        vals = scaled.loc[cls, present]
        rows[cls] = {"coverage": float(np.nanmean(vals)), "n_markers": len(present),
                     "n_cells": int((labels == cls).sum()),
                     **{f"scaled_{g}": float(vals[g]) for g in present}}
    out = pd.DataFrame(rows).T
    out.index.name = "cell_class"
    return out


def lineage_evidence(
    adata,
    programs: Mapping[str, Sequence[str]] | None = None,
    *,
    min_genes: int | None = None,
) -> pd.DataFrame:
    """Per cell, does it carry strong POSITIVE evidence for each lineage?

    Strong evidence = at least `min_genes` of that lineage's genes detected
    (count > 0). Detection, never absence: dropout can only hide evidence, so this
    under-calls and cannot manufacture a contradiction out of a zero. That is what
    makes it safe to challenge an NK call with T-lineage evidence without ever
    requiring a true NK cell to be literally TCR-negative.

    Returns a boolean frame, cells x lineages.
    """
    programs = config.LINEAGE_PROGRAMS if programs is None else programs
    min_genes = config.CONTRADICTION_MIN_GENES if min_genes is None else min_genes

    out = {}
    for name, genes in programs.items():
        present = [g for g in genes if g in adata.var_names]
        if not present:
            out[name] = np.zeros(adata.n_obs, dtype=bool)
            continue
        X = adata[:, present].X
        X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
        out[name] = (X > 0).sum(axis=1) >= min_genes
    return pd.DataFrame(out, index=adata.obs_names)


def contradiction_rate(
    adata,
    label_key: str,
    *,
    programs: Mapping[str, Sequence[str]] | None = None,
    pairs: Mapping[str, Sequence[str]] | None = None,
    min_genes: int | None = None,
    classes: Sequence[str] = BROAD_CLASSES,
) -> pd.DataFrame:
    """Share of each class's cells carrying evidence for an INCOMPATIBLE lineage.

    This is the complement to `marker_coverage`, and it exists because coverage asks
    a precision-like question — *do the cells I labelled X express X's markers?* — and
    therefore cannot see a class that has SWALLOWED cells from another lineage. In
    v1 the manual NK call scored coverage 1.00 while 63% of its cells were CD3D+:
    those cells really are NKG7/GNLY-high, because cytotoxic T cells are. Coverage had
    no way to object. This does.

    Per-lineage rates are reported alongside the maximum, so a failure names the
    lineage that caused it rather than just failing.
    """
    pairs = config.CONTRADICTION_PAIRS if pairs is None else pairs
    ev = lineage_evidence(adata, programs, min_genes=min_genes)
    labels = pd.Series(np.asarray(adata.obs[label_key]), index=adata.obs_names).astype(str)

    rows = {}
    for cls in classes:
        mask = (labels == cls).to_numpy()
        n = int(mask.sum())
        incompatible = list(pairs.get(cls, ()))
        row = {"n_cells": n, "n_incompatible_lineages": len(incompatible)}
        if n == 0 or not incompatible:
            # No cells, or a class with no contradictions by design (HSPC).
            # Not evaluable and zero are different outcomes; NaN says which.
            row["contradiction_rate"] = 0.0 if (n and not incompatible) else np.nan
            rows[cls] = row
            continue
        any_hit = np.zeros(n, dtype=bool)
        for lin in incompatible:
            if lin not in ev.columns:
                continue
            hit = ev[lin].to_numpy()[mask]
            row[f"rate_{lin}"] = float(hit.mean())
            any_hit |= hit
        row["contradiction_rate"] = float(any_hit.mean())
        rows[cls] = row
    out = pd.DataFrame(rows).T
    out.index.name = "cell_class"
    return out


def contradiction_concentration(
    adata,
    label_key: str,
    cluster_key: str,
    *,
    programs: Mapping[str, Sequence[str]] | None = None,
    pairs: Mapping[str, Sequence[str]] | None = None,
    min_genes: int | None = None,
    classes: Sequence[str] = BROAD_CLASSES,
) -> pd.DataFrame:
    """Where inside a class do its contradictory cells sit?

    The class-level rate from `contradiction_rate` cannot distinguish **mild diffuse
    noise** — a few percent scattered evenly, which is what residual doublets and
    ambient look like — from **one bad cluster** carrying almost all of it, which is a
    mislabelled population hiding inside an otherwise clean class. Those need different
    responses, and the single number hides which one you have.

    Returns one row per (class, cluster) with that cluster's own contradiction rate and
    `share_of_class_contradictions` — the fraction of the class's total contradictory
    cells living in that cluster. A class whose top one or two clusters hold most of its
    contradictions is concentrated; a flat profile is diffuse.

    Reporting only. This is not part of the veto and cannot change a decision.
    """
    pairs = config.CONTRADICTION_PAIRS if pairs is None else pairs
    ev = lineage_evidence(adata, programs, min_genes=min_genes)
    labels = pd.Series(np.asarray(adata.obs[label_key]), index=adata.obs_names).astype(str)
    clusters = pd.Series(np.asarray(adata.obs[cluster_key]), index=adata.obs_names).astype(str)

    rows = []
    for cls in classes:
        incompatible = [l for l in pairs.get(cls, ()) if l in ev.columns]
        mask = (labels == cls).to_numpy()
        if not incompatible or not mask.any():
            continue
        hit = np.zeros(mask.sum(), dtype=bool)
        for lin in incompatible:
            hit |= ev[lin].to_numpy()[mask]
        if not hit.any():
            continue                      # nothing to concentrate; omitted by design
        sub = pd.DataFrame({"cluster": clusters[mask].to_numpy(), "hit": hit})
        g = sub.groupby("cluster", observed=True)["hit"].agg(["size", "sum"])
        total = int(g["sum"].sum())
        for cl, r in g.iterrows():
            if r["sum"] == 0:
                continue
            rows.append({
                "cell_class": cls, "cluster": cl,
                "n_cells_in_class_and_cluster": int(r["size"]),
                "n_contradictory": int(r["sum"]),
                "cluster_contradiction_rate": float(r["sum"] / r["size"]),
                "share_of_class_contradictions": float(r["sum"] / total),
            })
    out = pd.DataFrame(rows)
    if len(out):
        out = out.sort_values(
            ["cell_class", "share_of_class_contradictions"], ascending=[True, False]
        ).reset_index(drop=True)
    return out


def cytotoxic_lineage_calls(adata, *, min_genes: int | None = None) -> pd.DataFrame:
    """Per-cell lineage calls for the cytotoxic-lymphocyte compartment.

    Identity is established only by coordinated lineage-defining machinery:
    `config.T_IDENTITY_ANCHORS` for T, `config.NK_IDENTITY` for NK. `config.T_CONTEXT`
    (TRBC1/TRBC2) and `config.CYTOTOXIC_STATE` (NKG7/GNLY/PRF1/GZMB/GZMA/CTSW) are
    measured and returned but can never make a cell positive for a lineage.

    The rule, stated once: **a transcript that frequently occurs without a lineage's
    defining machinery may provide context but cannot independently establish that
    lineage.** This is the same anchor/context principle already used by the Myeloid
    routes; it is applied here to T because the Leiden-23 diagnostic showed isolated
    TRBC expression is insufficient evidence of T-lineage commitment.

    `min_genes` defaults to `config.CONTRADICTION_MIN_GENES` — the project's existing
    positive-evidence rule. **No new numeric parameter is introduced.**
    """
    min_genes = config.CONTRADICTION_MIN_GENES if min_genes is None else min_genes
    layer = "counts" if "counts" in adata.layers else None

    def hits(genes):
        gs = [g for g in genes if g in adata.var_names]
        if not gs:
            return np.zeros(adata.n_obs, dtype=int)
        X = adata[:, gs].layers[layer] if layer else adata[:, gs].X
        X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
        return (X > 0).sum(axis=1)

    n_t = hits(config.T_IDENTITY_ANCHORS)
    n_ctx = hits(config.T_CONTEXT)
    n_nk = hits(config.NK_IDENTITY)
    n_gd = hits(config.GD_IDENTITY)
    n_cyt = hits(config.CYTOTOXIC_STATE)

    t_pos, nk_pos, gd_pos = n_t >= min_genes, n_nk >= min_genes, n_gd >= min_genes
    calls = []
    for t, nk, gd in zip(t_pos, nk_pos, gd_pos):
        if t and not nk:
            calls.append("T_gd" if gd else "T_ab")
        elif nk and not t:
            calls.append("NK")
        elif t and nk:
            calls.append("T_NK_mixed")
        else:
            calls.append("unresolved")
    return pd.DataFrame(
        {"n_T_identity": n_t, "n_T_context": n_ctx, "n_NK_identity": n_nk,
         "n_gd_identity": n_gd, "n_cytotoxic_state": n_cyt,
         "T_identity_pos": t_pos, "NK_identity_pos": nk_pos, "gd_pos": gd_pos,
         "T_context_pos": n_ctx > 0, "call": calls},
        index=adata.obs_names,
    )


def per_class_concordance(
    reference: Iterable[str],
    prediction: Iterable[str],
    *,
    classes: Sequence[str] = BROAD_CLASSES,
) -> pd.DataFrame:
    """Per-class F1 and Jaccard of `prediction` against `reference`, plus support.

    **This measures agreement, not accuracy** — `reference` here is the manual label
    set, which is a third opinion derived from the same matrix, not ground truth. Every
    caller must present the result as concordance.
    """
    ref = pd.Series(list(map(str, reference)), dtype=object).reset_index(drop=True)
    pred = pd.Series(list(map(str, prediction)), dtype=object).reset_index(drop=True)
    if len(ref) != len(pred):
        raise ValueError(f"length mismatch: {len(ref)} vs {len(pred)}")

    rows = {}
    for cls in classes:
        r, p = (ref == cls), (pred == cls)
        tp, fp, fn = int((r & p).sum()), int((~r & p).sum()), int((r & ~p).sum())
        prec = tp / (tp + fp) if tp + fp else 0.0
        rec = tp / (tp + fn) if tp + fn else 0.0
        f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
        union = int((r | p).sum())
        rows[cls] = {"f1": f1, "precision": prec, "recall": rec,
                     "jaccard": tp / union if union else 0.0,
                     "n_reference": int(r.sum()), "n_prediction": int(p.sum())}
    out = pd.DataFrame(rows).T
    out.index.name = "cell_class"
    return out


def decide_per_class(
    concordance: Mapping[str, pd.DataFrame],
    coverage: Mapping[str, pd.DataFrame],
    contradiction: Mapping[str, pd.DataFrame] | None = None,
    *,
    thresholds: Mapping[str, float] | None = None,
    coverage_min: float | None = None,
    contradiction_max: float | None = None,
    classes: Sequence[str] = BROAD_CLASSES,
    fallback: str = "manual",
) -> pd.DataFrame:
    """Apply the pre-declared per-class rule and return the decision table.

    `concordance` and `coverage` are keyed by method name; `concordance` holds only the
    automated methods (there is nothing to compare the manual labels against but
    themselves), `coverage` should also carry the fallback method so its biological
    support is on the record.

    The rule, in this order and no other:

    1. **Two biological vetoes, applied before any concordance is consulted.**
       A method must clear BOTH for a class or it cannot win that class, whatever its
       concordance — this is what stops two automated methods agreeing their way past
       biology.
       - **Marker coverage** (`coverage_min`): do the assigned cells express the
         class's own markers? Catches a class with no marker support.
       - **Lineage exclusivity** (`contradiction_max`): do they also carry strong
         evidence for an incompatible lineage? Catches a class that has swallowed
         another one — which coverage structurally cannot see, and which is exactly
         how v1 shipped 33,556 "NK" cells that were 63% CD3D+.
    2. Of the methods that survive the veto, one must also clear the class's
       concordance bar.
    3. Where several qualify, the highest F1 wins and the rest are recorded as also
       having passed.
    4. Where none qualifies, the class falls back to `fallback` — **but only if the
       fallback itself passes the coverage veto**. If it does not, the class is
       `Ambiguous` and that is reported as an unresolved class, not papered over.

    Thresholds are read from `config` and are **never** adjusted after results are
    seen; that is the entire point of declaring them in advance.
    """
    thresholds = config.CONCORDANCE_THRESHOLDS if thresholds is None else thresholds
    coverage_min = config.MARKER_COVERAGE_MIN if coverage_min is None else coverage_min
    contradiction_max = (
        config.CONTRADICTION_MAX_RATE if contradiction_max is None else contradiction_max
    )
    # A veto that was NOT COMPUTED is not in force; a veto that was computed and came
    # back NaN for a class is NOT_EVALUABLE. Those are different things and conflating
    # them is the bug this revision exists to fix, so the distinction is explicit.
    contradiction_in_force = bool(contradiction)
    contradiction = {} if contradiction is None else contradiction

    def cov(method: str, cls: str) -> float:
        tbl = coverage.get(method)
        if tbl is None or cls not in tbl.index:
            return float("nan")
        return float(tbl.loc[cls, "coverage"])

    def contra(method: str, cls: str) -> float:
        tbl = contradiction.get(method)
        if tbl is None or cls not in tbl.index:
            return float("nan")
        return float(tbl.loc[cls, "contradiction_rate"])

    def cov_state(method: str, cls: str) -> str:
        return evaluate_veto(cov(method, cls), coverage_min, direction="min")

    def contra_state(method: str, cls: str) -> str:
        if not contradiction_in_force:
            return PASS
        return evaluate_veto(contra(method, cls), contradiction_max, direction="max")

    rows = {}
    for cls in classes:
        bar = float(thresholds.get(cls, 0.0))
        cand = []
        for method, tbl in concordance.items():
            f1 = float(tbl.loc[cls, "f1"]) if cls in tbl.index else float("nan")
            c = cov(method, cls)
            x = contra(method, cls)
            cand.append({"method": method, "f1": f1, "coverage": c, "contradiction": x,
                         "coverage_state": cov_state(method, cls),
                         "contradiction_state": contra_state(method, cls),
                         "f1_ok": bool(np.isfinite(f1) and f1 >= bar)})
        # Only an explicit PASS on BOTH vetoes makes a method eligible. NOT_EVALUABLE is
        # neither a pass nor a veto: it cannot win a class and is not recorded as vetoed.
        eligible = [d for d in cand if d["coverage_state"] == PASS
                    and d["contradiction_state"] == PASS and d["f1_ok"]]
        vetoed = [d["method"] for d in cand
                  if d["f1_ok"] and d["coverage_state"] == FAIL]
        contra_vetoed = [d["method"] for d in cand
                         if d["f1_ok"] and d["coverage_state"] == PASS
                         and d["contradiction_state"] == FAIL]
        not_evaluable = [f"{d['method']}:cov" for d in cand
                         if d["coverage_state"] == NOT_EVALUABLE]
        not_evaluable += [f"{d['method']}:contra" for d in cand
                          if d["contradiction_state"] == NOT_EVALUABLE]

        if eligible:
            win = max(eligible, key=lambda d: d["f1"])
            chosen, reason = win["method"], "passed concordance + coverage"
            passed = ",".join(sorted(d["method"] for d in eligible))
            f1_used, cov_used = win["f1"], win["coverage"]
        else:
            fb_cov_state, fb_contra_state = cov_state(fallback, cls), contra_state(fallback, cls)
            fb_cov = cov(fallback, cls)
            if fb_cov_state == PASS and fb_contra_state == PASS:
                chosen, reason = fallback, "no automated method qualified; fell back"
                f1_used, cov_used, passed = float("nan"), fb_cov, ""
            else:
                chosen = AMBIGUOUS
                if NOT_EVALUABLE in (fb_cov_state, fb_contra_state):
                    which = "coverage" if fb_cov_state == NOT_EVALUABLE else "contradiction"
                    reason = (f"no method passed; fallback {fallback} is NOT_EVALUABLE "
                              f"for {which} (it assigns no cells to this class)")
                else:
                    why = "coverage" if fb_cov_state == FAIL else "lineage-exclusivity"
                    reason = f"no method passed; fallback also failed the {why} veto"
                f1_used, cov_used, passed = float("nan"), fb_cov, ""

        rows[cls] = {
            "chosen_method": chosen, "reason": reason,
            "f1_threshold": bar, "coverage_min": coverage_min,
            "f1_chosen": f1_used, "coverage_chosen": cov_used,
            "also_passed": passed,
            "vetoed_by_coverage": ",".join(sorted(vetoed)),
            "vetoed_by_contradiction": ",".join(sorted(contra_vetoed)),
            "not_evaluable": ",".join(sorted(not_evaluable)),
            **{f"cov_state_{d['method']}": d["coverage_state"] for d in cand},
            **{f"contra_state_{d['method']}": d["contradiction_state"] for d in cand},
            f"cov_state_{fallback}": cov_state(fallback, cls),
            f"contra_state_{fallback}": contra_state(fallback, cls),
            "contradiction_max": contradiction_max,
            "contradiction_chosen": contra(str(chosen), cls),
            **{f"f1_{d['method']}": d["f1"] for d in cand},
            **{f"coverage_{d['method']}": d["coverage"] for d in cand},
            **{f"contradiction_{d['method']}": d["contradiction"] for d in cand},
            f"coverage_{fallback}": cov(fallback, cls),
            f"contradiction_{fallback}": contra(fallback, cls),
        }
    out = pd.DataFrame(rows).T
    out.index.name = "cell_class"
    return out


def assemble_final_labels(
    adata,
    decision: pd.DataFrame,
    label_keys: Mapping[str, str],
    *,
    conf_keys: Mapping[str, str] | None = None,
    out_key: str = "cell_type",
    source_key: str = "annotation_source",
    conf_key: str = "annotation_conf",
) -> pd.DataFrame:
    """Resolve per-method labels into one `cell_type`, deterministically.

    THIS LAYER RESOLVES COMPETING OUTPUTS OF DECISIONS ALREADY MADE. It reads only the
    per-method labels and the class-level decision table. It invents no marker
    evidence, consults no expression, and alters no class-level threshold.

    v4's rule was: class *C* may only be assigned by the method that won *C*. That is
    bookkeeping, not biology, and it lost cells. Where CellTypist, SingleR/Novershtern
    and SingleR/HPCA all said Myeloid but `manual` had won the Myeloid decision and
    called the cluster Ambiguous, nothing could claim those cells and 6,533 unambiguous
    monocytes became ownerless. It also gave `PlasmaCell` silent precedence in any
    multi-claim tie, purely because it comes first in `BROAD_CLASSES`.

    The replacement, applied per cell in this order:

    1. **Authoritative claims.** Class *C* is authoritatively claimed when the method
       that won *C* labels this cell *C*. Exactly one such claim wins outright — the
       class-level decision named that method the adjudicator for that class, so its
       claim is the strongest evidence assembly has.
    2. **Conflict between authoritative claims** — two winners each claiming their own
       class for the same cell — is a genuine disagreement, not a precedence question.
       It falls through to the consensus test rather than being resolved by ordering.
    3. **Consensus.** With no single authoritative claim, count the classes assigned by
       every available method. A class is taken only on a **strict majority** of the
       methods that expressed an opinion (`Ambiguous` is not an opinion). A 2-2 or
       2-1-1 split is not consensus and yields `Ambiguous`.
    4. **A class-level `Ambiguous` verdict is never overridden.** If the decision table
       could not resolve class *C*, assembly cannot assign *C* however many methods
       vote for it — that verdict is a biological finding and this layer must not
       reverse it.

    `annotation_source` records `"consensus"` for cells resolved at step 3, so
    provenance stays traceable and the two paths are distinguishable downstream.
    """
    n = adata.n_obs
    methods = [m for m in label_keys if label_keys[m] in adata.obs]
    labels = {m: np.asarray(adata.obs[label_keys[m]]).astype(str) for m in methods}

    resolvable = {
        cls: str(decision.loc[cls, "chosen_method"])
        for cls in decision.index
        if str(decision.loc[cls, "chosen_method"]) not in (AMBIGUOUS, "nan")
    }
    # class -> its adjudicating method, for classes that were resolved at class level
    winners = {c: m for c, m in resolvable.items() if m in labels}

    final = np.full(n, AMBIGUOUS, dtype=object)
    source = np.full(n, "", dtype=object)
    conf = np.full(n, np.nan, dtype=float)

    for i in range(n):
        claims = [c for c, m in winners.items() if labels[m][i] == c]
        if len(claims) == 1:
            final[i] = claims[0]
            source[i] = winners[claims[0]]
            continue
        # 0 claims, or >1 conflicting authoritative claims -> consensus
        votes: dict[str, int] = {}
        for m in methods:
            lab = labels[m][i]
            # Neither an Ambiguous verdict nor an unevaluated method is an OPINION.
            # Counting either would let abstentions dilute a real majority.
            if lab and lab not in (AMBIGUOUS, NOT_EVALUABLE, "nan", "None"):
                votes[lab] = votes.get(lab, 0) + 1
        total = sum(votes.values())
        if not total:
            continue
        top = max(votes, key=lambda c: votes[c])
        if votes[top] * 2 > total and top in winners:
            final[i] = top
            source[i] = "consensus"

    if conf_keys:
        for m, key in conf_keys.items():
            if key not in adata.obs:
                continue
            vals = np.asarray(adata.obs[key], dtype=float)
            take = source == m
            conf[take] = vals[take]

    adata.obs[out_key] = pd.Categorical(final, categories=list(BROAD_CLASSES) + [AMBIGUOUS])
    adata.obs[source_key] = pd.Categorical(source)
    adata.obs[conf_key] = conf

    summary = (
        pd.Series(final).value_counts().rename("n_cells").to_frame()
        .assign(pct=lambda d: 100 * d["n_cells"] / n)
    )
    summary.index.name = "cell_type"
    return summary


def score_state_programs(
    adata,
    programs: Mapping[str, Sequence[str]] | None = None,
    *,
    prefix: str = "program_",
) -> list[str]:
    """Score orthogonal cell-state programs as **continuous, non-exclusive** covariates.

    These never become cell types and never enter `cell_type`. A cycling plasma cell is
    PlasmaCell plus a high `program_cell_cycle`, not a "Cycling" identity.
    """
    import scanpy as sc
    import warnings

    programs = config.STATE_PROGRAMS if programs is None else programs
    written = []
    for name, genes in programs.items():
        present = [g for g in genes if g in adata.var_names]
        if not present:
            warnings.warn(f"program {name}: no genes present, not scored", stacklevel=2)
            continue
        if len(present) < len(genes):
            warnings.warn(
                f"program {name}: {len(genes) - len(present)} gene(s) absent, scoring on {len(present)}",
                stacklevel=2,
            )
        col = f"{prefix}{name}"
        sc.tl.score_genes(adata, present, score_name=col)
        written.append(col)
    return written


def composition_by(adata, group_key: str, *, label_key: str = "cell_type") -> pd.DataFrame:
    """Cell-type composition per group, as counts and within-group fractions.

    Proportions are compositional data — one type rising mechanically pushes the others
    down — so any *comparison* between groups needs `scCODA`, not a per-type test on
    these fractions. This function reports; it does not test.
    """
    counts = (
        adata.obs.groupby([group_key, label_key], observed=True)
        .size().unstack(fill_value=0).sort_index()
    )
    frac = counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0)
    return pd.concat({"n_cells": counts, "fraction": frac}, axis=1)
