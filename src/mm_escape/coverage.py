"""Supplemental multi-antigen coverage — uncovered fractions and incremental gain.

Design frozen in
``results/08_dual_antigen_escape/multi_antigen_coverage/multi_antigen_design.md``.

**This consumes frozen Stage-07/08 malignant-cell and measurement infrastructure; it does
not reopen the frozen BCMA/GPRC5D analysis.** Nothing here recomputes a Stage-08 number, and
nothing here may write into a frozen stage's namespace.

Two rules are enforced structurally rather than by convention:

* **Depth strata are read, never re-derived.** The primary analysis uses the frozen
  per-cell ``depth_stratum_cohort`` from Stage 08, and any stratification that genuinely
  must be recomputed calls :mod:`mm_escape.antigen`'s utility. There is no binning code in
  this module — ``test_no_local_depth_binning`` asserts that by inspecting the source.
* **No utility score.** Uncovered fraction, incremental gain, measurement reliability,
  depth robustness and normal-marrow expression stay separate. There is no principled
  weight to combine them with, and a composite would hide the inputs a reader could
  otherwise disagree with.
"""
from __future__ import annotations

from itertools import combinations

import numpy as np

from . import antigen as _antigen

__all__ = [
    "TARGETS", "ANCHOR", "BACKGROUND_LINEAGES", "AMBIENT_NOT_EVALUABLE",
    "ELIGIBLE", "NOT_EVALUABLE", "RELIABILITY_LABELS",
    "TECHNICAL_ZERO_MAX", "DEPTH_RHO_MAX", "DEPTH_SPREAD_MAX", "BACKGROUND_SEPARATION_MIN",
    "detected", "uncovered_fraction", "coverage_row", "all_combinations",
    "incremental_gains", "reliability_label", "eligibility",
    "shared_depth_strata", "stratified_expected_co_negative", "unconditioned_expected_co_negative",
]

#: The candidate panel from the analysis plan. Order is fixed so combination ids are stable.
TARGETS = ("TNFRSF17", "GPRC5D", "SLAMF7", "FCRL5", "CD38", "SDC1", "ITGB7")

#: The project question. Alternatives are compared *to* this pair, never substituted for it.
ANCHOR = ("TNFRSF17", "GPRC5D")

#: Declared per target from biology in the frozen design §5, before any background
#: distribution was inspected. A background population is only a background population if
#: the target is genuinely absent from it — several of these are broadly expressed across
#: marrow lineages, and reusing BCMA/GPRC5D's reference wholesale would invent a floor.
BACKGROUND_LINEAGES: dict[str, tuple[str, ...]] = {
    # frozen stage-08 reference, inherited unchanged
    "TNFRSF17": ("Tcell", "Myeloid", "Bcell", "HSPC", "NK"),
    "GPRC5D": ("Tcell", "Myeloid", "Bcell", "HSPC", "NK"),
    # B-lineage restricted: B cells are NOT a negative for FCRL5
    "FCRL5": ("Tcell", "Myeloid", "Erythroid", "HSPC", "NK"),
    # CD138 is absent from mature haematopoietic non-plasma lineages
    "SDC1": ("Tcell", "Myeloid", "Bcell", "Erythroid", "HSPC", "NK"),
    # CS1 is on NK, CD8 T subsets, activated monocytes, DC and B cells — no clean
    # lymphoid or myeloid negative exists in marrow
    "SLAMF7": ("Erythroid",),
    # CD38 is on activated T/NK, B subsets, monocytes and progenitors — it is itself an
    # HSPC marker, so HSPC is emphatically not a negative
    "CD38": ("Erythroid",),
    # integrin beta-7 is broadly lymphoid (gut-homing T and B) and on NK
    "ITGB7": ("Myeloid", "Erythroid"),
}

AMBIENT_NOT_EVALUABLE = "AMBIENT_FLOOR_NOT_EVALUABLE"
ELIGIBLE = "COVERAGE_ELIGIBLE"
NOT_EVALUABLE = "COVERAGE_NOT_EVALUABLE"

RELIABILITY_LABELS = ("comparatively_reliable", "depth_sensitive",
                      "strongly_dropout_prone", "not_evaluable")

#: Frozen in design §7/§9 before any target's floor was computed.
TECHNICAL_ZERO_MAX = 0.50        # above this, a gene of that abundance reads zero for depth alone
DEPTH_RHO_MAX = 0.20             # detection-vs-depth Spearman above this is depth-sensitive
DEPTH_SPREAD_MAX = 2.0           # >= 2-fold detection spread across strata within a cohort
BACKGROUND_SEPARATION_MIN = 2.0  # malignant detection must clear background by this factor


def detected(counts):
    """The frozen Stage-08 positivity rule, applied unchanged to every target: ``count > 0``.

    Deliberately not tuned per target. A per-gene threshold fitted to each target's
    abundance would be the post-hoc tuning the design forbids, and would make uncovered
    fractions incomparable between targets. Variation in reliability is handled by the QC
    gate, not by moving the cutoff — a dropout-prone target does not get a kinder rule, it
    gets excluded or explicitly qualified.

    A pure function of raw integer counts. Nothing normalised, imputed, denoised or
    smoothed may reach it.
    """
    c = np.asarray(counts)
    if c.dtype.kind == "f" and not np.all(np.isfinite(c) & (c == np.floor(c))):
        raise ValueError("positivity must be called on raw integer counts, not "
                         "normalised/imputed values")
    return c > 0


def uncovered_fraction(detection, targets, combination):
    """Fraction of cells negative for **every** target in `combination`.

    `detection` is a (n_cells, n_targets) boolean array; `targets` names its columns.
    """
    detection = np.asarray(detection, dtype=bool)
    if detection.ndim != 2:
        raise ValueError("detection must be 2-D (cells x targets)")
    if detection.shape[0] == 0:
        return np.nan
    idx = [list(targets).index(t) for t in combination]
    return float((~detection[:, idx]).all(axis=1).mean())


def all_combinations(eligible, sizes=(1, 2, 3)):
    """Every combination of the given sizes, in the fixed `TARGETS` order."""
    ordered = [t for t in TARGETS if t in set(eligible)]
    out = []
    for k in sizes:
        out.extend(combinations(ordered, k))
    return out


def coverage_row(detection, targets, combination):
    """Uncovered fraction plus the per-member marginals, for one combination."""
    row = {"combination": "+".join(combination), "size": len(combination),
           "uncovered": uncovered_fraction(detection, targets, combination)}
    for t in combination:
        row[f"uncovered_{t}"] = uncovered_fraction(detection, targets, (t,))
    return row


def incremental_gains(detection, targets, combination):
    """Directional gain from adding each member to the rest of the combination.

    ``gain_X_given_rest = P(rest all negative) - P(rest and X all negative)``.

    Reported per direction and never summarised into one number: a pair can carry a large
    gain in one direction and a small one in the other, and that asymmetry is the answer to
    a single- vs. dual-target question.
    """
    out = []
    for t in combination:
        rest = tuple(x for x in combination if x != t)
        if not rest:
            continue
        base = uncovered_fraction(detection, targets, rest)
        joint = uncovered_fraction(detection, targets, combination)
        out.append({"combination": "+".join(combination), "added": t,
                    "given": "+".join(rest), "uncovered_given": base,
                    "uncovered_combination": joint, "gain": base - joint})
    return out


def reliability_label(technical_zero, depth_rho, depth_spread, evaluable=True):
    """Descriptive measurement-reliability label. **Never a weight.**

    Order matters: dropout dominance is checked before depth sensitivity, because a gene
    that reads zero half the time for depth reasons alone is not merely "depth-sensitive".
    """
    if not evaluable:
        return "not_evaluable"
    if np.isfinite(technical_zero) and technical_zero >= TECHNICAL_ZERO_MAX:
        return "strongly_dropout_prone"
    if (np.isfinite(depth_rho) and depth_rho >= DEPTH_RHO_MAX) or (
            np.isfinite(depth_spread) and depth_spread >= DEPTH_SPREAD_MAX):
        return "depth_sensitive"
    return "comparatively_reliable"


def eligibility(*, ambient_status, malignant_detection, background_detection,
                technical_zero, circularity_blocked=False):
    """`COVERAGE_ELIGIBLE` / `COVERAGE_NOT_EVALUABLE` from the frozen design §9 rules only.

    Never a function of how good the target's coverage looks — that information is not even
    passed in, which is the point.
    """
    reasons = []
    if ambient_status == AMBIENT_NOT_EVALUABLE:
        sep = (malignant_detection / background_detection
               if np.isfinite(background_detection) and background_detection > 0 else np.inf)
        if not np.isfinite(sep) or sep < BACKGROUND_SEPARATION_MIN:
            reasons.append("no clean background population and malignant detection is not "
                           f"{BACKGROUND_SEPARATION_MIN:g}x its background")
    if np.isfinite(technical_zero) and technical_zero >= TECHNICAL_ZERO_MAX:
        reasons.append(f"expression-matched technical-zero fraction {technical_zero:.2f} "
                       f">= {TECHNICAL_ZERO_MAX:g}")
    if circularity_blocked:
        reasons.append("fails the differentiation/circularity checkpoint (design §11)")
    return (NOT_EVALUABLE if reasons else ELIGIBLE), "; ".join(reasons)


# ---------------------------------------------------------------------------
# Shared depth machinery. These are thin delegations on purpose: this module must
# contain no binning logic of its own, and the tests check the source for it.
# ---------------------------------------------------------------------------

def shared_depth_strata(depth, n_bins=5, min_cells=20):
    """Depth strata via the validated Stage-08 utility — never a local reimplementation.

    Used only where a stratification genuinely must be recomputed (the truncate-10k
    sensitivity). The primary analysis reads Stage 08's frozen per-cell assignment instead,
    so no perturbation of any target can move a cell between strata.
    """
    edges = _antigen.quantile_edges(depth, n_bins=n_bins)
    return _antigen.merge_sparse_strata(_antigen.assign_strata(depth, edges),
                                        min_cells=min_cells)


def stratified_expected_co_negative(a_neg, b_neg, strata):
    """Depth-conditioned independence baseline, via the frozen Stage-08 closed form.

    A *technical baseline* the observed value is compared against — never a correction.
    Multiplying the marginals assumes exactly the independence a co-loss test exists to
    interrogate.
    """
    return _antigen.stratified_expected_dn(a_neg, b_neg, strata)


def unconditioned_expected_co_negative(a_neg, b_neg):
    """Patient-wide ``P(A-) * P(B-)``, reported beside the conditioned value.

    The gap between the two is the depth artifact, quantified.
    """
    return _antigen.unconditioned_expected_dn(a_neg, b_neg)
