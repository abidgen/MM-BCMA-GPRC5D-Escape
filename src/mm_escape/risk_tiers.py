"""Final risk tiering — the frozen rule from `risk_tier_final/risk_tier_policy.md`.

The whole design point: **`robust-high` means the high-DN conclusion survives not just the
technical sensitivity analyses, but also the reasonable uncertainty over what numerical DN
fraction should count as "high."**

Two things this module deliberately does not do. It never reads a cohort label, so cohort
cannot enter the rule. And it never reads a bulk/Stage-09 field, so marginal validation
cannot upgrade or downgrade a patient — Stage 09 is interpretation context, not a scoring
axis, and it does not override joint-state uncertainty.
"""
from __future__ import annotations

import numpy as np

__all__ = [
    "TAU_LOW", "TAU_HIGH_SET", "MIN_N_PRIMARY", "MIN_CELLS_PER_SAMPLE",
    "ROBUST_HIGH", "ROBUST_LOW", "UNCERTAIN", "THRESHOLD_ROBUST", "THRESHOLD_DEPENDENT",
    "UNCERTAINTY_REASONS", "side", "sample_agreement", "tier_at", "final_tier",
]

#: Predeclared in stage08_predeclaration.md §9 as the smallest DN fraction the project
#: claims to resolve. Not altered because `robust-low` may be empty.
TAU_LOW = 0.05

#: Transparent project conventions, NOT clinically validated boundaries. Stage 09 found no
#: biologically validated joint-DN threshold and no cutoff is supported by bulk.
TAU_HIGH_SET = (0.20, 0.25, 0.33)

MIN_N_PRIMARY = 100          # frozen stage-08 minimum-n rule
MIN_CELLS_PER_SAMPLE = 20    # pre-existing stage-08 per-stratum floor

ROBUST_HIGH = "robust-high"
ROBUST_LOW = "robust-low"
UNCERTAIN = "uncertain"
THRESHOLD_ROBUST = "THRESHOLD_ROBUST"
THRESHOLD_DEPENDENT = "THRESHOLD_DEPENDENT"

UNCERTAINTY_REASONS = (
    "uncertain_low_n", "uncertain_denominator", "uncertain_depth",
    "uncertain_repeated_sample", "uncertain_null_scheme",
    "uncertain_dropout_compatible", "uncertain_intermediate_dn", "uncertain_threshold",
)


def side(x, tau_high, tau_low=TAU_LOW):
    """The only numeric comparison in the rule. Everything else is a flip test on it."""
    if x is None or not np.isfinite(x):
        return "not_evaluable"
    if x >= tau_high:
        return "high"
    if x <= tau_low:
        return "low"
    return "intermediate"


def sample_agreement(per_sample, tau_high, tau_low=TAU_LOW,
                     min_cells=MIN_CELLS_PER_SAMPLE):
    """Do a patient's assessable samples fall on the same side of the boundary?

    Returns "agree", "disagree", or "not_assessable_for_agreement". A patient with fewer
    than two assessable samples **cannot disagree and is never credited with agreeing** —
    the distinction matters, because most patients here contribute one sample and their
    stability is unobserved rather than demonstrated.

    `per_sample` is an iterable of (sample_id, n_cells, dn_fraction).
    """
    usable = [(s, n, v) for s, n, v in per_sample if n >= min_cells and np.isfinite(v)]
    if len(usable) < 2:
        return "not_assessable_for_agreement"
    sides = {side(v, tau_high, tau_low) for _, _, v in usable}
    return "agree" if len(sides) == 1 else "disagree"


def tier_at(ev, tau_high, tau_low=TAU_LOW):
    """Apply the frozen rule at one `TAU_HIGH`. Returns (tier, reasons dict).

    `ev` carries only Stage-07/08 evidence: n_primary, dn_primary, dn_sensitivity,
    dn_trunc10k, enr_ci_lo, enr_cohortbins, enr_globalbins, per_sample. No cohort, no bulk.
    """
    reasons = {r: False for r in UNCERTAINTY_REASONS}

    s_prim = side(ev["dn_primary"], tau_high, tau_low)
    s_sens = side(ev["dn_sensitivity"], tau_high, tau_low)
    s_trunc = side(ev["dn_trunc10k"], tau_high, tau_low)
    agreement = sample_agreement(ev.get("per_sample", []), tau_high, tau_low)

    reasons["uncertain_low_n"] = bool(ev["n_primary"] < MIN_N_PRIMARY)
    reasons["uncertain_denominator"] = bool(s_sens != s_prim)
    reasons["uncertain_depth"] = bool(s_trunc != s_prim)
    reasons["uncertain_repeated_sample"] = bool(agreement == "disagree")

    # qualitative sign flip across the null itself, not a tuned delta
    a, b = ev.get("enr_cohortbins", np.nan), ev.get("enr_globalbins", np.nan)
    reasons["uncertain_null_scheme"] = bool(
        np.isfinite(a) and np.isfinite(b) and (a - 1.0) * (b - 1.0) < 0)

    reasons["uncertain_intermediate_dn"] = bool(s_prim == "intermediate")

    # "beyond depth-conditioned marginal dropout": the enrichment CI excludes 1 from
    # above. NOT the permutation p-value, which never gates.
    enr_lo = ev.get("enr_ci_lo", np.nan)
    beyond_dropout = bool(np.isfinite(enr_lo) and enr_lo > 1.0)

    blocking = any(reasons[r] for r in (
        "uncertain_low_n", "uncertain_denominator", "uncertain_depth",
        "uncertain_repeated_sample", "uncertain_null_scheme"))

    if s_prim == "high" and not blocking and beyond_dropout:
        return ROBUST_HIGH, reasons
    if s_prim == "high" and not blocking and not beyond_dropout:
        reasons["uncertain_dropout_compatible"] = True
        return UNCERTAIN, reasons
    if s_prim == "low" and not blocking and not beyond_dropout:
        return ROBUST_LOW, reasons
    if s_prim == "low" and not blocking and beyond_dropout:
        reasons["uncertain_dropout_compatible"] = True
        return UNCERTAIN, reasons
    return UNCERTAIN, reasons


def final_tier(ev, tau_high_set=TAU_HIGH_SET, tau_low=TAU_LOW):
    """Collapse the three runs. `robust-high` requires it under **all three** thresholds.

    Qualifying only at the most permissive cutoff is not enough: the high-risk reading
    must not depend on whether "meaningfully elevated" means 20%, 25% or 33%.
    """
    per_tau = {t: tier_at(ev, t, tau_low) for t in tau_high_set}
    tiers = {t: v[0] for t, v in per_tau.items()}
    reasons = {r: any(v[1][r] for v in per_tau.values()) for r in UNCERTAINTY_REASONS}

    robustness = (THRESHOLD_ROBUST if len(set(tiers.values())) == 1
                  else THRESHOLD_DEPENDENT)
    if robustness == THRESHOLD_DEPENDENT:
        reasons["uncertain_threshold"] = True

    if all(v == ROBUST_HIGH for v in tiers.values()):
        final = ROBUST_HIGH
    elif all(v == ROBUST_LOW for v in tiers.values()):
        final = ROBUST_LOW
    else:
        final = UNCERTAIN
    return final, robustness, tiers, reasons
