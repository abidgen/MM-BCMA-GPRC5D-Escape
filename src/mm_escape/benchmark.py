"""
Stage 05b — integration-method benchmark.

WHY THIS EXISTS
---------------
Stage 05 used Harmony because it was the obvious default, never because it was compared
against anything. `sc-best-practices`' integration chapter recommends running several
methods and scoring them with scIB rather than assuming one wins. This module is the
machinery for doing that on *this* dataset, where the standard scoring would give the
wrong answer.

WHY THE SCORING IS NOT STANDARD scIB
------------------------------------
Stage 04 found the deposit is censored per cohort: WashU cohorts 1 and 2 were cut at
10,000 UMIs before deposit, MMRF and the donors were not. Plasma cells are the
highest-RNA-content cells in marrow, so MMRF's two largest plasma clusters hold 68% and
88% of their cells above that ceiling — cells WashU cannot contain.

The consequence for benchmarking is specific and severe. **scIB's batch metrics cannot
distinguish "correctly left apart" from "failed to merge."** A method that forces the
three plasma islands together scores better on kBET/iLISI while manufacturing
correspondence between populations where none is recoverable. A naive global ranking
would therefore structurally reward overcorrection.

Note the careful claim: this is **not** an assertion that the cohorts have biologically
different plasma cells. It is a **non-recoverable sampling/censoring asymmetry** —
WashU's *observed* distribution is missing its high-RNA portion, so no one-to-one
correspondence remains for any method to recover. That argument is stronger precisely
because it does not depend on the difference being biological.

So: **the immune compartment is scored, the plasma compartment is diagnosed.** Global
scIB is computed and reported as a secondary reference, never used for selection.

WHAT THIS BENCHMARK CANNOT DO
-----------------------------
**No integration method restores cells that were never deposited.** A beautifully mixed
latent space has not undone the ascertainment bias in the raw counts that stage 08
reads. Whichever arm wins, stage 08 still owes its truncate-all-cohorts-at-10,000
sensitivity analysis. Selecting scVI must not create the impression the censoring was
"handled".

And the blast radius is modest by construction: the embedding feeds only stages 06 and
11. Antigen calls read `layers["counts"]` and malignant subclustering is per-patient
un-integrated, so **no outcome here can move `frac_double_negative`.**
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from . import integration

__all__ = [
    "ArmSpec",
    "CELLTYPIST_TO_BROAD",
    "provisional_labels",
    "ARMS",
    "shared_pca",
    "run_unintegrated",
    "run_harmony",
    "run_scvi",
    "run_scanorama",
    "run_arm",
    "depth_association",
    "compartment_mixing",
    "decide",
    "render_decision",
    "IMMUNE_CLASSES",
    "DECISION_TOLERANCES",
]

#: Broad classes treated as the immune compartment for primary scoring. `plasma` is
#: deliberately absent — it is diagnosed, not scored. Erythroid/HSPC are included
#: because they are marrow populations the embedding must still keep distinct.
IMMUNE_CLASSES = ("T", "NK", "B", "myeloid", "erythroid", "HSPC")

#: The incumbent's name, used as the comparison baseline throughout.
INCUMBENT = "harmony_stage05"

#: CellTypist `Immune_All_High` label -> broad benchmark class. Anything absent maps to
#: `other` and is excluded from scoring.
#:
#: **Two things were checked rather than assumed** when this was written, against the
#: real predictions on all 172,940 cells:
#:
#: 1. **`ILC` is NK.** It is 8% of the marrow, which would be implausible for genuine
#:    innate lymphoid cells. Marker check: NKG7 98.7%, GNLY 92.2%, KLRD1 85.9%,
#:    MS4A1 1.3%. `Immune_All_High` folds NK into its ILC class.
#: 2. **`Immune_All_High` does cover erythroid and HSPC here**, contrary to the
#:    expectation recorded in the project plan that an immune-only reference would be
#:    blind to
#:    them — Erythroid 14,103 cells at HBB 99.7%, HSC/MPP 2,625 at CD34 58.4%. So no
#:    hand-set marker override is needed, which keeps arbitrary thresholds out of the
#:    benchmark entirely.
#:
#: Every other retained class was marker-verified the same way: T cells CD3D 85.8%,
#: B cells MS4A1 86.4%, Plasma cells MZB1 96.0%, Monocytes LYZ 99.1%.
CELLTYPIST_TO_BROAD: dict[str, str] = {
    "T cells": "T",
    "Double-positive thymocytes": "T",
    "Double-negative thymocytes": "T",
    "ETP": "T",
    "ILC": "NK",
    "ILC precursor": "NK",
    "B cells": "B",
    "B-cell lineage": "B",
    "Plasma cells": "plasma",
    "Monocytes": "myeloid",
    "Macrophages": "myeloid",
    "Mono-mac": "myeloid",
    "Monocyte precursor": "myeloid",
    "DC": "myeloid",
    "pDC": "myeloid",
    "DC precursor": "myeloid",
    "Promyelocytes": "myeloid",
    "Myelocytes": "myeloid",
    "Granulocytes": "myeloid",
    "Mast cells": "myeloid",
    "Erythroid": "erythroid",
    "Erythrocytes": "erythroid",
    "HSC/MPP": "HSPC",
    # Megakaryocytic and non-haematopoietic calls are a few hundred cells between them
    # and belong to none of the seven project classes; they fall through to `other`
    # and are excluded from scoring rather than forced into a class they do not fit.
}

#: Tolerances for the decision rule, **declared before running** — this mirrors how
#: stage 06 pre-declares its F1 thresholds and stage 10 pre-registers the γ-secretase
#: hypothesis. Fixing them in advance is what stops "choose the best" from becoming
#: post-hoc rationalisation of whichever arm looks tidier.
DECISION_TOLERANCES = {
    # A candidate must actually beat the incumbent on batch removal, not merely tie.
    "batch_improvement_min": 0.01,
    # Biology may not be sacrificed for mixing. A small tolerance absorbs run-to-run
    # noise without licensing a real loss.
    "bio_loss_max": 0.02,
    # Depth predictability (R^2 of total_counts on the latent space) must not rise
    # materially above the incumbent's.
    "depth_r2_increase_max": 0.05,
    # A jump in plasma mixing this large triggers the overcorrection check.
    "plasma_mixing_jump": 0.15,
}


@dataclass(frozen=True)
class ArmSpec:
    """One benchmark arm: a method plus the batch definition it corrects on."""

    name: str
    method: str
    batch_keys: tuple[str, ...]
    note: str = ""


#: The arms. One common batch key (`sample_name`, the true technical unit) so the
#: comparison is of *methods*; plus two reference arms so the incumbent and the
#: only-non-confounded batch definition are both represented.
#:
#: `sample_name` and `patient_id` are nearly the same partition here — **42 of 50
#: patients contribute exactly one sample** — so "avoid correcting on patient" is a
#: weaker argument than it looks, and correcting on either may erase genuine
#: between-patient immune biology. `cohort` (4 levels) is the only batch definition NOT
#: confounded with patient, and is where the demonstrated distortion actually lives,
#: which is why it is a real arm rather than an afterthought.
#:
#: BBKNN is absent on purpose: it produces a corrected neighbour graph and no
#: embedding, so `scib_metrics.Benchmarker` — which scores `obsm` keys — cannot place
#: it on the same footing. Deriving an embedding from its graph would compare a
#: different kind of object.
ARMS: tuple[ArmSpec, ...] = (
    ArmSpec("unintegrated", "none", (), "baseline: same HVGs/scaling/PCs, no correction"),
    ArmSpec("harmony_sample", "harmony", ("sample_name",)),
    ArmSpec("scvi_sample", "scvi", ("sample_name",)),
    ArmSpec("scanorama_sample", "scanorama", ("sample_name",)),
    ArmSpec(INCUMBENT, "harmony", ("patient_id", "n_genes_ref", "cohort"),
            "the stage-05 configuration, as shipped"),
    ArmSpec("harmony_cohort", "harmony", ("cohort",), "the only non-patient-confounded key"),
    ArmSpec("scvi_cohort", "scvi", ("cohort",), "as above"),
)


# ---------------------------------------------------------------------------
# Provisional labels — the non-circular half of the scoring
# ---------------------------------------------------------------------------

def provisional_labels(adata: AnnData, *, model: str = "Immune_All_High.pkl",
                       inplace: bool = True) -> pd.DataFrame:
    """CellTypist labels for scoring, deliberately embedding-independent.

    scIB scores batch removal *against* biological conservation, and the bio half needs
    cell-type labels. Ours come from stage 06 — which consumes the embedding this
    benchmark is choosing. CellTypist breaks that circle because it classifies from
    log-normalized expression and never sees an embedding.

    **`majority_voting=False` is load-bearing, not a default.** Majority voting
    smooths predictions over an over-clustering, and that clustering is computed from a
    representation — which would smuggle an embedding straight back into the labels the
    embeddings are then scored against.

    Requires `adata.X` to be log1p(CP10K), which is what stage 05 leaves behind and
    what CellTypist documents as its input.

    Returns a frame with `celltypist_label` (the raw call), `broad_label` (collapsed via
    `CELLTYPIST_TO_BROAD`) and `compartment` (immune / plasma / other). These are
    **provisional and for scoring only** — stage 06 still runs its full three-method
    comparison, and nothing downstream should read them.
    """
    import celltypist

    celltypist.models.download_models(model=[model], force_update=False)
    prediction = celltypist.annotate(adata, model=model, majority_voting=False)
    raw = prediction.predicted_labels["predicted_labels"].astype(str)
    raw.index = adata.obs_names

    broad = raw.map(CELLTYPIST_TO_BROAD).fillna("other")
    compartment = np.where(
        broad.isin(IMMUNE_CLASSES), "immune",
        np.where(broad == "plasma", "plasma", "other"),
    )
    frame = pd.DataFrame(
        {"celltypist_label": raw, "broad_label": broad, "compartment": compartment},
        index=adata.obs_names,
    )
    if inplace:
        for column in frame.columns:
            adata.obs[column] = pd.Categorical(frame[column])
    return frame


# ---------------------------------------------------------------------------
# Runners — every arm sees identical cells, genes, HVGs and seeds
# ---------------------------------------------------------------------------

def shared_pca(adata: AnnData, *, n_comps: int = 50, random_state: int = 0) -> np.ndarray:
    """Scale HVGs on a throwaway copy and return PCs, computed ONCE for all arms.

    Computing this once is what makes the comparison about integration rather than
    about preprocessing: the unintegrated baseline and every Harmony arm start from the
    same coordinates, and scVI/Scanorama start from the same HVG set.
    """
    subset = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(subset, max_value=10)
    sc.tl.pca(subset, n_comps=n_comps, svd_solver="arpack", random_state=random_state)
    return np.asarray(subset.obsm["X_pca"])


def run_unintegrated(adata: AnnData, batch_keys: Sequence[str], *, pca: np.ndarray,
                     **_) -> np.ndarray:
    """The baseline: the shared PCs, uncorrected.

    Deliberately PCA space rather than raw or log expression, so the contrast against
    the corrected arms isolates batch correction alone.
    """
    return pca


def run_harmony(adata: AnnData, batch_keys: Sequence[str], *, pca: np.ndarray,
                max_iter: int = 20, random_state: int = 0, **_) -> np.ndarray:
    """Harmony on the shared PCs. Accepts one or several covariates."""
    from harmonypy import run_harmony as _run

    meta = adata.obs[list(batch_keys)].astype(str)
    result = _run(pca, meta, list(batch_keys), max_iter_harmony=max_iter,
                  random_state=random_state)
    corrected = np.asarray(result.Z_corr)
    if corrected.shape != pca.shape:
        corrected = corrected.T
    if corrected.shape != pca.shape:
        raise ValueError(
            f"Harmony returned {np.asarray(result.Z_corr).shape}, neither {pca.shape} "
            f"nor its transpose."
        )
    return corrected


def run_scvi(adata: AnnData, batch_keys: Sequence[str], *, n_latent: int = 30,
             max_epochs: int | None = None, random_state: int = 0, **_) -> np.ndarray:
    """scVI on raw counts over the shared HVGs.

    Counts, not log-normalized values: scVI models the count distribution and the
    library size explicitly, which is the whole reason it is a principled candidate on
    a dataset whose confound *is* depth.

    Extra batch keys beyond the first are passed as categorical covariates rather than
    concatenated, so the incumbent's multi-covariate configuration has a fair analogue.
    """
    import scvi

    scvi.settings.seed = random_state
    subset = adata[:, adata.var["highly_variable"]].copy()
    subset.X = subset.layers["counts"].copy()

    primary, *extra = list(batch_keys)
    scvi.model.SCVI.setup_anndata(
        subset, batch_key=primary,
        categorical_covariate_keys=[k for k in extra] or None,
    )
    model = scvi.model.SCVI(subset, n_latent=n_latent)
    model.train(max_epochs=max_epochs)
    return np.asarray(model.get_latent_representation())


def run_scanorama(adata: AnnData, batch_keys: Sequence[str], *, pca: np.ndarray,
                  n_comps: int = 50, **_) -> np.ndarray:
    """Scanorama over the shared HVGs, batches ordered deterministically.

    Scanorama returns per-batch blocks, so the rows are reassembled back into the
    caller's cell order — getting that wrong yields an embedding that clusters happily
    and means nothing.
    """
    import scanorama

    if len(batch_keys) != 1:
        raise ValueError("Scanorama takes a single batch key.")
    key = batch_keys[0]
    subset = adata[:, adata.var["highly_variable"]]
    labels = adata.obs[key].astype(str)

    order = sorted(labels.unique())
    blocks, indices = [], []
    for value in order:
        mask = np.flatnonzero((labels == value).to_numpy())
        indices.append(mask)
        block = subset[mask].X
        blocks.append(np.asarray(block.todense()) if hasattr(block, "todense")
                      else np.asarray(block))
    genes = [list(subset.var_names)] * len(blocks)

    corrected, _ = scanorama.integrate(blocks, genes, dimred=n_comps)

    out = np.empty((adata.n_obs, corrected[0].shape[1]), dtype=np.float32)
    for mask, block in zip(indices, corrected):
        out[mask] = np.asarray(block, dtype=np.float32)
    return out


RUNNERS: dict[str, Callable[..., np.ndarray]] = {
    "none": run_unintegrated,
    "harmony": run_harmony,
    "scvi": run_scvi,
    "scanorama": run_scanorama,
}


def run_arm(adata: AnnData, arm: ArmSpec, *, pca: np.ndarray, **kwargs) -> np.ndarray:
    """Dispatch one arm and sanity-check the shape it returns."""
    if arm.method not in RUNNERS:
        raise ValueError(f"unknown method {arm.method!r}; known: {sorted(RUNNERS)}")
    embedding = RUNNERS[arm.method](adata, arm.batch_keys, pca=pca, **kwargs)
    embedding = np.asarray(embedding)
    if embedding.shape[0] != adata.n_obs:
        raise ValueError(
            f"arm {arm.name!r} returned {embedding.shape[0]} rows for {adata.n_obs} "
            f"cells — the batch blocks were probably not reassembled in cell order."
        )
    return embedding


# ---------------------------------------------------------------------------
# Dataset-specific diagnostics
# ---------------------------------------------------------------------------

def depth_association(embedding: np.ndarray, depth: np.ndarray) -> float:
    """R² of sequencing depth regressed on the latent space.

    **This statistic is fixed before the benchmark runs, and the choice is not
    arbitrary.** R² from a linear regression depends only on the *column span* of the
    embedding, so it is **rotation-invariant** — which matters because latent axes are
    arbitrary and differ between methods. A per-dimension statistic such as
    `max |Spearman(depth, dim_k)|` is not rotation-invariant and would rank methods on
    an accident of their parameterisation rather than on how much depth information
    their representation carries.

    `depth` is log1p-transformed first: library size is roughly log-normal, and
    regressing the raw counts would let a handful of 250,000-UMI plasma cells dominate.

    Higher means the representation encodes depth more strongly. For this dataset that
    is a warning sign, because depth is confounded with cohort by the stage-04
    censoring.
    """
    x = np.asarray(embedding, dtype=np.float64)
    y = np.log1p(np.asarray(depth, dtype=np.float64))
    if x.shape[0] != y.shape[0]:
        raise ValueError(f"embedding has {x.shape[0]} rows, depth has {y.shape[0]}")
    if y.size < 2 or np.allclose(y, y[0]):
        return 0.0

    design = np.column_stack([np.ones(x.shape[0]), x - x.mean(axis=0)])
    coefficients, *_ = np.linalg.lstsq(design, y, rcond=None)
    residual = y - design @ coefficients
    total = float(((y - y.mean()) ** 2).sum())
    if total == 0:
        return 0.0
    return float(1.0 - (residual ** 2).sum() / total)


def compartment_mixing(adata: AnnData, *, cluster_key: str, batch_key: str = "cohort",
                       compartment_key: str = "compartment") -> pd.DataFrame:
    """Cohort-mixing entropy per cluster, split by compartment.

    Thin wrapper over `integration.batch_mixing` so there is one implementation of the
    entropy, not two. Returns a per-compartment median alongside the per-cluster table.
    """
    table = integration.batch_mixing(adata, batch_key=batch_key, cluster_key=cluster_key)
    compartment = (
        adata.obs.groupby(cluster_key, observed=True)[compartment_key]
        .agg(lambda values: values.astype(str).mode().iat[0])
    )
    table[compartment_key] = table[cluster_key].map(compartment).to_numpy()
    return table


# ---------------------------------------------------------------------------
# The decision rule, as code
# ---------------------------------------------------------------------------

def decide(scores: pd.DataFrame, *, incumbent: str = INCUMBENT,
           tolerances: dict[str, float] | None = None) -> pd.DataFrame:
    """Apply the pre-declared rule and return one row per arm.

    `scores` needs one row per arm indexed by arm name, with columns
    `batch_score`, `bio_score`, `depth_r2`, `plasma_mixing`.

    The verdict is **computed, not narrated** — `decision.md` is rendered from this
    frame, so the write-up cannot drift from the rule. A method replaces the incumbent
    only if all four criteria hold:

    1. `batch_improved`      immune batch removal actually improves;
    2. `bio_preserved`       broad cell identity is not sacrificed for that mixing;
    3. `depth_ok`            depth predictability does not rise materially;
    4. `overcorrection_ok`   plasma mixing never *earns* a method anything. A large
                             jump alone only flags it for inspection, since a method
                             could improve plasma geometry legitimately — but a jump
                             **together with** rising depth association is
                             disqualifying, because that pair is the signature of the
                             censoring being smoothed over rather than respected.
    """
    tolerances = {**DECISION_TOLERANCES, **(tolerances or {})}
    required = {"batch_score", "bio_score", "depth_r2", "plasma_mixing"}
    missing = required - set(scores.columns)
    if missing:
        raise ValueError(f"scores is missing {sorted(missing)}")
    if incumbent not in scores.index:
        raise ValueError(f"incumbent {incumbent!r} not among arms {list(scores.index)}")

    base = scores.loc[incumbent]
    out = scores.copy()

    out["batch_improved"] = (
        out["batch_score"] > base["batch_score"] + tolerances["batch_improvement_min"]
    )
    out["bio_preserved"] = out["bio_score"] >= base["bio_score"] - tolerances["bio_loss_max"]
    depth_rise = out["depth_r2"] - base["depth_r2"]
    out["depth_ok"] = depth_rise <= tolerances["depth_r2_increase_max"]

    plasma_jump = out["plasma_mixing"] - base["plasma_mixing"]
    out["plasma_flagged"] = plasma_jump > tolerances["plasma_mixing_jump"]
    # A flag alone is not fatal; a flag with rising depth association is.
    out["overcorrection_ok"] = ~(out["plasma_flagged"] & (depth_rise > 0))

    out["eligible"] = (
        out["batch_improved"] & out["bio_preserved"]
        & out["depth_ok"] & out["overcorrection_ok"]
    )
    out.loc[incumbent, "eligible"] = False  # the incumbent cannot replace itself
    return out


def _markdown_table(frame: pd.DataFrame, *, digits: int = 4) -> str:
    """Render a DataFrame as a markdown table.

    Hand-rolled rather than `DataFrame.to_markdown()`, which needs `tabulate` —
    `mm-core` does not carry it, and adding a dependency to `env-core.yml` for one
    formatting call is exactly the kind of casual change that broke that env once.
    """
    header = ["arm", *(str(column) for column in frame.columns)]
    rows = [header, ["---"] * len(header)]
    for name, row in frame.iterrows():
        cells = [str(name)]
        for value in row:
            if isinstance(value, (bool, np.bool_)):
                cells.append("yes" if value else "no")
            elif isinstance(value, (int, float, np.floating, np.integer)):
                cells.append(f"{float(value):.{digits}g}")
            else:
                cells.append(str(value))
        rows.append(cells)
    return "\n".join("| " + " | ".join(row) + " |" for row in rows)


def render_decision(decision: pd.DataFrame, *, incumbent: str = INCUMBENT,
                    subsample: int | None = None, notes: str = "") -> str:
    """Render `decision.md` FROM the decision table, so prose cannot drift from it."""
    eligible = decision.index[decision["eligible"]].tolist()
    if eligible:
        winner = decision.loc[eligible, "batch_score"].idxmax()
        verdict = (
            f"**{winner}** replaces the incumbent: it is the eligible arm with the "
            f"strongest immune batch removal."
        )
    else:
        winner = incumbent
        verdict = (
            f"**No arm qualified. `{incumbent}` stays.** That is a real result, not a "
            f"failure to find one — the incumbent survived a comparison it could have "
            f"lost."
        )

    lines = [
        "# Stage 05b — integration benchmark decision",
        "",
        "## The rule, declared before running",
        "",
        "A method replaces the incumbent only if **all four** hold. Thresholds are in",
        "`benchmark.DECISION_TOLERANCES` and were fixed before any arm was run.",
        "",
        "1. `batch_improved` — immune batch removal actually improves.",
        "2. `bio_preserved` — broad cell identity is not sacrificed for that mixing.",
        "3. `depth_ok` — R²(depth ~ latent) does not rise materially. R² is used",
        "   because it depends only on the embedding's column span and is therefore",
        "   rotation-invariant; latent axes are arbitrary across methods.",
        "4. `overcorrection_ok` — plasma-cell mixing never earns a method anything. A",
        "   jump alone flags it; a jump **with** rising depth association disqualifies",
        "   it, that pair being the signature of the censoring being smoothed over.",
        "",
        "Scoring is on the **immune compartment**. Global scIB is reported as a",
        "secondary reference and is never used for selection.",
        "",
        "## Verdict",
        "",
        verdict,
        "",
        "## Per-arm results",
        "",
        _markdown_table(decision),
        "",
        "## What this benchmark cannot do",
        "",
        "**No integration method restores cells that were never deposited.** WashU",
        "cohorts 1 and 2 were cut at 10,000 UMIs before deposit; the high-RNA portion",
        "of their plasma-cell population is absent from the counts, whichever arm",
        "wins. A well-mixed latent space has not undone that ascertainment bias.",
        "**Stage 08 still owes its truncate-all-cohorts-at-10,000 sensitivity",
        "analysis, regardless of the outcome here.**",
        "",
        "Nor can this move the headline metric: the embedding feeds only stages 06 and",
        "11, antigen calls read `layers['counts']`, and malignant subclustering is",
        "per-patient and un-integrated.",
    ]
    if subsample:
        lines += ["", f"Metrics computed on a stratified subsample of {subsample:,} cells."]
    if notes:
        lines += ["", notes]
    return "\n".join(lines) + "\n"
