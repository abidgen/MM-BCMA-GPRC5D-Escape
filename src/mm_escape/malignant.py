"""
Stage 07 — dominant-clone membership from immunoglobulin evidence.

SINGLE AXIS. SAY SO EVERY TIME.
--------------------------------
CNV was attempted as an independent malignancy axis and rejected on its healthy-donor
negative control, before any disease CNV was inspected. It is frozen `NOT_EVALUABLE`
cohort-wide. What remains is kappa/lambda class plus V-gene usage — the same molecule,
the same biological event. `CLONE_SUPPORTED` therefore means *high-confidence membership
in the dominant plasma-cell clone*, *not* malignancy proven by two orthogonal assays.

If stage 08 uses it as the denominator it is **a high-specificity dominant-clone
denominator**, never "exhaustive recovery of all malignant plasma cells".

V ABSENCE IS NOT NEGATIVE EVIDENCE
-----------------------------------
10x 3' chemistry captures V segments at ~1 UMI and J segments essentially not at all.
A missing dominant-V call may be dropout, shallow depth, a true non-clone cell, or
something else. So:

    positive V detection  -> may establish clone support
    absent V detection    -> may NEVER establish incompatibility

Incompatibility requires its own positive evidence: minority-chain class, or a coherent
alternative V at `config.ALT_V_MIN_UMI`.
"""

from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np
import pandas as pd

from . import config

__all__ = [
    "CLONE_SUPPORTED", "CLONE_COMPATIBLE_V_UNOBSERVED", "CLONE_INCOMPATIBLE",
    "CLONE_UNCERTAIN", "V_EVALUABLE", "V_PARTIAL", "V_NOT_EVALUABLE",
    "light_chain_class", "patient_v_evaluability", "clone_membership",
]

CLONE_SUPPORTED = "CLONE_SUPPORTED"
CLONE_COMPATIBLE_V_UNOBSERVED = "CLONE_COMPATIBLE_V_UNOBSERVED"
CLONE_INCOMPATIBLE = "CLONE_INCOMPATIBLE"
CLONE_UNCERTAIN = "CLONE_UNCERTAIN"

V_EVALUABLE = "V_EVALUABLE"
V_PARTIAL = "V_PARTIAL"
V_NOT_EVALUABLE = "V_NOT_EVALUABLE"


def light_chain_class(kappa_umi, lambda_umi) -> np.ndarray:
    """Per-cell kappa/lambda call from RAW light-chain UMI. Ratio, never presence.

    Immunoglobulin is the most ambient-contaminated transcript family in this tissue,
    so a presence call is far noisier than it looks; a ratio is robust to a shared
    additive background.
    """
    k = np.asarray(kappa_umi, dtype=float)
    l = np.asarray(lambda_umi, dtype=float)
    tot = k + l
    frac_k = np.divide(k, np.maximum(tot, 1))
    out = np.full(len(k), "insufficient", dtype=object)
    ok = tot >= config.LC_CLASS_MIN_UMI
    out[ok & (frac_k >= config.LC_CLASS_MIN_FRAC)] = "kappa"
    out[ok & ((1 - frac_k) >= config.LC_CLASS_MIN_FRAC)] = "lambda"
    out[ok & (out == "insufficient")] = "ambiguous"
    return out


def patient_v_evaluability(
    n_v_positive: int, pct_lcv: float, top_v_frac: float, enrichment: float
) -> str:
    """Three-state V evaluability, from the frozen constants.

    `top_v_frac` sits inside the empirical gap between healthy donors (max 0.378) and
    evaluable disease patients (min 0.562). A patient below it is not "negative" — it
    is simply not usable for V-based clone membership, and its cells fall through to
    `CLONE_UNCERTAIN` rather than being called normal.
    """
    if not np.isfinite(top_v_frac) or top_v_frac < config.DOMINANT_V_MIN_FRAC:
        return V_NOT_EVALUABLE
    if not np.isfinite(enrichment) or enrichment < config.DOMINANT_V_MIN_ENRICHMENT:
        return V_NOT_EVALUABLE
    if n_v_positive >= config.V_EVALUABLE_MIN_CELLS and pct_lcv >= config.V_EVALUABLE_MIN_PCT:
        return V_EVALUABLE
    if n_v_positive >= config.V_PARTIAL_MIN_CELLS and pct_lcv >= config.V_PARTIAL_MIN_PCT:
        return V_PARTIAL
    return V_NOT_EVALUABLE


def clone_membership(
    patient_clonality: str,
    v_state: str,
    lc_class: str,
    dominant_class: str,
    dominant_v_detected: bool,
    alt_v_detected: bool,
) -> str:
    """The per-cell state machine. Every combination is explicit; nothing is implicit.

    Order matters and is deliberate:

    1. **Positive incompatibility first.** Minority-chain class, or a coherent
       alternative V, is positive evidence against membership and outranks everything.
       Absence of the dominant V is *not* in this list.
    2. **Support requires all four**: a confidently clonal patient, a V-evaluable
       patient, dominant-class compatibility, and positive dominant-V detection.
    3. **Compatible-but-unobserved** is its own state, never folded into either
       support or incompatibility — it is exactly the population whose status the data
       cannot settle.
    4. **Everything else is uncertain**, including every cell of a patient that is not
       confidently clonal or not V-evaluable.
    """
    # 1. positive incompatibility
    if dominant_class in ("kappa", "lambda") and lc_class in ("kappa", "lambda") \
            and lc_class != dominant_class:
        return CLONE_INCOMPATIBLE
    if alt_v_detected and not dominant_v_detected:
        return CLONE_INCOMPATIBLE

    # patient-level preconditions for any positive membership claim
    if patient_clonality != "CLONAL_STRONG" or v_state != V_EVALUABLE:
        return CLONE_UNCERTAIN

    if lc_class != dominant_class:          # ambiguous / insufficient light chain
        return CLONE_UNCERTAIN

    # 2 / 3
    return CLONE_SUPPORTED if dominant_v_detected else CLONE_COMPATIBLE_V_UNOBSERVED
