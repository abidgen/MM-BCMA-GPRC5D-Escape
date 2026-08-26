"""Stage 07 clone-membership rules. Data-free.

The invariant this file exists to protect: **V absence is not negative evidence**, and
the classifier is completely independent of TNFRSF17/GPRC5D.
"""
from __future__ import annotations

import numpy as np
import pytest

from mm_escape import config
from mm_escape import malignant as M

STRONG, EVAL = "CLONAL_STRONG", M.V_EVALUABLE


def call(**kw):
    base = dict(patient_clonality=STRONG, v_state=EVAL, lc_class="kappa",
                dominant_class="kappa", dominant_v_detected=True, alt_v_detected=False)
    base.update(kw)
    return M.clone_membership(**base)


# --- A. positive support ---------------------------------------------------
def test_A_dominant_class_plus_dominant_v_is_clone_supported():
    assert call() == M.CLONE_SUPPORTED


# --- B. V dropout ----------------------------------------------------------
def test_B_dominant_class_without_v_is_compatible_v_unobserved():
    assert call(dominant_v_detected=False) == M.CLONE_COMPATIBLE_V_UNOBSERVED


# --- C. minority chain -----------------------------------------------------
def test_C_minority_chain_cannot_be_clone_supported():
    assert call(lc_class="lambda") == M.CLONE_INCOMPATIBLE
    assert call(lc_class="lambda", dominant_v_detected=True) != M.CLONE_SUPPORTED


# --- D. V absence alone ----------------------------------------------------
def test_D_v_absence_alone_cannot_create_incompatibility():
    """THE hard invariant. 10x 3' captures V at ~1 UMI; absence is uninformative."""
    assert call(dominant_v_detected=False) != M.CLONE_INCOMPATIBLE
    assert call(dominant_v_detected=False, lc_class="ambiguous") != M.CLONE_INCOMPATIBLE
    assert call(dominant_v_detected=False, lc_class="insufficient") != M.CLONE_INCOMPATIBLE


# --- E. V-not-evaluable patient -------------------------------------------
@pytest.mark.parametrize("state", [M.V_PARTIAL, M.V_NOT_EVALUABLE])
def test_E_non_evaluable_patient_cannot_yield_clone_supported(state):
    assert call(v_state=state) == M.CLONE_UNCERTAIN


def test_E_non_clonal_patient_cannot_yield_clone_supported():
    for c in ("CLONAL_WEAK", "NO_RESTRICTION", "NOT_EVALUABLE"):
        assert call(patient_clonality=c) == M.CLONE_UNCERTAIN


# --- F. alternative coherent V --------------------------------------------
def test_F_alternative_v_is_not_folded_into_the_dominant_clone():
    assert call(dominant_v_detected=False, alt_v_detected=True) == M.CLONE_INCOMPATIBLE


def test_F_alt_v_requires_two_umi_not_one():
    """1 UMI is the typical V detection level and must not establish incompatibility."""
    assert config.ALT_V_MIN_UMI >= 2


# --- ambiguous / insufficient light chain ---------------------------------
@pytest.mark.parametrize("lc", ["ambiguous", "insufficient"])
def test_uncallable_light_chain_is_uncertain_not_supported(lc):
    assert call(lc_class=lc) == M.CLONE_UNCERTAIN
    assert call(lc_class=lc, dominant_v_detected=False) == M.CLONE_UNCERTAIN


def test_every_combination_is_explicit():
    """No input combination may fall through to something undefined."""
    valid = {M.CLONE_SUPPORTED, M.CLONE_COMPATIBLE_V_UNOBSERVED,
             M.CLONE_INCOMPATIBLE, M.CLONE_UNCERTAIN}
    for pc in ("CLONAL_STRONG", "CLONAL_WEAK", "NO_RESTRICTION", "NOT_EVALUABLE"):
        for vs in (M.V_EVALUABLE, M.V_PARTIAL, M.V_NOT_EVALUABLE):
            for lc in ("kappa", "lambda", "ambiguous", "insufficient"):
                for dv in (True, False):
                    for av in (True, False):
                        assert M.clone_membership(pc, vs, lc, "kappa", dv, av) in valid


# --- light-chain classing --------------------------------------------------
def test_light_chain_class_is_a_ratio_not_a_presence_call():
    # 10 kappa / 1 lambda -> kappa; 5/5 -> ambiguous; 1/1 -> below the UMI floor
    out = M.light_chain_class([10, 5, 1, 0], [1, 5, 1, 0])
    assert list(out) == ["kappa", "ambiguous", "insufficient", "insufficient"]


def test_light_chain_floor_and_fraction_are_the_frozen_constants():
    assert config.LC_CLASS_MIN_UMI == 3 and config.LC_CLASS_MIN_FRAC == 0.80


# --- evaluability gate -----------------------------------------------------
def test_dominant_v_threshold_sits_inside_the_donor_disease_gap():
    """Donors max 0.378; evaluable disease begins 0.562."""
    assert 0.378 < config.DOMINANT_V_MIN_FRAC < 0.562


def test_top_v_below_threshold_is_never_evaluable():
    assert M.patient_v_evaluability(5000, 0.99, 0.40, 50.0) == M.V_NOT_EVALUABLE


def test_low_enrichment_v_gene_is_not_evaluable():
    """A V gene common across the cohort carries no membership information."""
    assert M.patient_v_evaluability(5000, 0.99, 0.99, 1.5) == M.V_NOT_EVALUABLE


def test_evaluability_tiers():
    assert M.patient_v_evaluability(60, 0.6, 0.9, 20) == M.V_EVALUABLE
    assert M.patient_v_evaluability(25, 0.3, 0.9, 20) == M.V_PARTIAL
    assert M.patient_v_evaluability(5, 0.05, 0.9, 20) == M.V_NOT_EVALUABLE


# --- G / H / I: antigen invariance, CNV isolation, rule immutability -------
def test_G_no_antigen_gene_appears_anywhere_in_stage07_rules():
    import inspect
    src = inspect.getsource(M)
    for g in ("TNFRSF17", "GPRC5D"):
        assert g not in src, f"{g} appears in the stage-07 classifier source"


def test_G_clone_membership_signature_admits_no_antigen_input():
    import inspect
    params = set(inspect.signature(M.clone_membership).parameters)
    assert params == {"patient_clonality", "v_state", "lc_class", "dominant_class",
                      "dominant_v_detected", "alt_v_detected"}


def test_H_cnv_takes_no_part_in_stage07_clone_calls():
    """CNV is frozen NOT_EVALUABLE cohort-wide; no CNV input may reach these rules."""
    import inspect
    src = inspect.getsource(M).lower()
    assert "cnv" not in inspect.signature(M.clone_membership).parameters
    assert "infercnv" not in src


def test_I_stage07_constants_are_the_frozen_values():
    """An antigen-bias diagnostic must never feed back into these."""
    assert config.DOMINANT_V_MIN_FRAC == 0.50
    assert config.V_EVALUABLE_MIN_CELLS == 50 and config.V_PARTIAL_MIN_CELLS == 20
    assert config.V_EVALUABLE_MIN_PCT == 0.50 and config.V_PARTIAL_MIN_PCT == 0.20
    assert config.DOMINANT_V_MIN_ENRICHMENT == 3.0
    assert config.ALT_V_MIN_UMI == 2
