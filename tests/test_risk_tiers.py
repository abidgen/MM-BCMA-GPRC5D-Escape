"""Frozen risk-tier rule, invariants A-I. Written before the rule was applied to data."""
import inspect

import numpy as np
import pytest

from mm_escape import risk_tiers as RT


def ev(**kw):
    """A patient that is cleanly robust-high unless a field is overridden."""
    base = dict(n_primary=500, dn_primary=0.40, dn_sensitivity=0.42, dn_trunc10k=0.41,
                enr_ci_lo=1.05, enr_cohortbins=1.10, enr_globalbins=1.12,
                per_sample=[("s1", 300, 0.40), ("s2", 200, 0.41)])
    base.update(kw)
    return base


# ------------------------------------------------------ A. significance not required
def test_robust_high_without_any_p_value_present():
    """No p-value is supplied at all; the patient must still qualify.

    Significance would encode sequencing power as biology — the four significant patients
    all sit in the deepest cohort — so it may support a high call but never gate it.
    """
    e = ev()
    assert "perm_p" not in e
    final, robustness, tiers, _ = RT.final_tier(e)
    assert final == RT.ROBUST_HIGH and robustness == RT.THRESHOLD_ROBUST


def test_rule_never_reads_a_p_value():
    src = inspect.getsource(RT)
    code = "\n".join(l for l in src.splitlines() if not l.strip().startswith("#"))
    assert "perm_p" not in code and "pvalue" not in code and "p_value" not in code


# ----------------------------------------------------- B. significance not sufficient
@pytest.mark.parametrize("override,reason", [
    (dict(n_primary=50), "uncertain_low_n"),
    (dict(dn_sensitivity=0.10), "uncertain_denominator"),
    (dict(dn_trunc10k=0.10), "uncertain_depth"),
    (dict(per_sample=[("s1", 300, 0.40), ("s2", 200, 0.02)]), "uncertain_repeated_sample"),
    (dict(enr_globalbins=0.90), "uncertain_null_scheme"),
])
def test_strong_enrichment_cannot_override_a_robustness_failure(override, reason):
    e = ev(enr_ci_lo=3.0, enr_cohortbins=3.5, **override)
    final, _, _, reasons = RT.final_tier(e)
    assert final == RT.UNCERTAIN
    assert reasons[reason] is True


# --------------------------------------------------------------- C. threshold invariance
def test_only_tau_high_varies_across_the_three_runs():
    """Same evidence, same code path; only the cutoff differs."""
    e = ev(dn_primary=0.30, dn_sensitivity=0.30, dn_trunc10k=0.30,
           per_sample=[("s1", 300, 0.30)])
    out = {t: RT.tier_at(e, t)[0] for t in RT.TAU_HIGH_SET}
    assert out[0.20] == RT.ROBUST_HIGH and out[0.25] == RT.ROBUST_HIGH
    assert out[0.33] == RT.UNCERTAIN          # 0.30 < 0.33 -> intermediate
    assert RT.tier_at(e, 0.33)[1]["uncertain_intermediate_dn"] is True


# --------------------------------------------------------- D. threshold-dependent high
def test_high_at_020_but_not_at_033_is_finally_uncertain():
    e = ev(dn_primary=0.22, dn_sensitivity=0.23, dn_trunc10k=0.22,
           per_sample=[("s1", 300, 0.22)])
    final, robustness, tiers, reasons = RT.final_tier(e)
    assert tiers[0.20] == RT.ROBUST_HIGH and tiers[0.33] == RT.UNCERTAIN
    assert final == RT.UNCERTAIN
    assert robustness == RT.THRESHOLD_DEPENDENT and reasons["uncertain_threshold"] is True


# ------------------------------------------------------------------- E. robust-high
def test_high_under_all_three_thresholds_with_no_failures_is_robust_high():
    e = ev(dn_primary=0.55, dn_sensitivity=0.57, dn_trunc10k=0.56,
           per_sample=[("s1", 400, 0.55), ("s2", 100, 0.58)])
    final, robustness, tiers, reasons = RT.final_tier(e)
    assert final == RT.ROBUST_HIGH and robustness == RT.THRESHOLD_ROBUST
    assert set(tiers.values()) == {RT.ROBUST_HIGH}
    assert not any(reasons.values())


def test_large_dn_with_no_evidence_beyond_dropout_is_not_robust_high():
    """The case the criterion exists for: burden alone must not promote a patient."""
    e = ev(dn_primary=0.60, dn_sensitivity=0.62, dn_trunc10k=0.61, enr_ci_lo=0.99,
           enr_cohortbins=1.00, enr_globalbins=1.00, per_sample=[("s1", 800, 0.60)])
    final, _, _, reasons = RT.final_tier(e)
    assert final == RT.UNCERTAIN and reasons["uncertain_dropout_compatible"] is True


# ------------------------------------------------------------ F. empty robust-low OK
def test_suite_does_not_require_any_robust_low_patient_to_exist():
    """An empty robust-low tier is an acceptable scientific outcome, so the rule is
    tested on synthetic evidence and never asserted against the real cohort."""
    e = ev(dn_primary=0.02, dn_sensitivity=0.03, dn_trunc10k=0.02, enr_ci_lo=0.95,
           enr_cohortbins=0.98, enr_globalbins=0.99, per_sample=[("s1", 500, 0.02)])
    assert RT.final_tier(e)[0] == RT.ROBUST_LOW


def test_low_dn_with_co_negativity_excess_is_not_robust_low():
    e = ev(dn_primary=0.02, dn_sensitivity=0.03, dn_trunc10k=0.02, enr_ci_lo=1.15,
           enr_cohortbins=1.40, enr_globalbins=1.20, per_sample=[("s1", 500, 0.02)])
    final, _, _, reasons = RT.final_tier(e)
    assert final == RT.UNCERTAIN and reasons["uncertain_dropout_compatible"] is True


# --------------------------------------------------------------- G. Stage-09 isolation
def test_bulk_fields_cannot_change_a_tier():
    base = RT.final_tier(ev())
    for extra in ({}, {"bulk_TNFRSF17_tpm": 500.0, "bulk_GPRC5D_tpm": 0.0},
                  {"bulk_available": False}, {"bulk_rho": -0.9}):
        assert RT.final_tier(ev(**extra))[:3] == base[:3]


class _RecordingEv(dict):
    """Records every evidence key the rule actually reads."""

    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        self.read = set()

    def __getitem__(self, key):
        self.read.add(key)
        return super().__getitem__(key)

    def get(self, key, default=None):
        self.read.add(key)
        return super().get(key, default)


ALLOWED_EVIDENCE_KEYS = {
    "n_primary", "dn_primary", "dn_sensitivity", "dn_trunc10k",
    "enr_ci_lo", "enr_cohortbins", "enr_globalbins", "per_sample",
}


def test_rule_reads_only_whitelisted_stage07_08_evidence():
    """Stronger than text matching: record what the rule actually touches.

    Nothing bulk-derived, nothing cohort-derived, and no p-value can influence a tier if
    the rule never reads such a key.
    """
    rec = _RecordingEv(ev(cohort="MMRF", bulk_GPRC5D_tpm=500.0, perm_p=1e-9,
                          bulk_available=True))
    RT.final_tier(rec)
    assert rec.read <= ALLOWED_EVIDENCE_KEYS, f"rule read {rec.read - ALLOWED_EVIDENCE_KEYS}"
    for forbidden in ("cohort", "bulk_GPRC5D_tpm", "perm_p", "bulk_available"):
        assert forbidden not in rec.read


# ------------------------------------------------------------------ H. cohort isolation
def test_cohort_label_cannot_change_a_tier():
    base = RT.final_tier(ev())
    for c in ("MMRF", "WU1", "WU2", None):
        assert RT.final_tier(ev(cohort=c))[:3] == base[:3]


def test_rule_reads_no_cohort_key_even_when_one_is_supplied():
    rec = _RecordingEv(ev(cohort="WU2"))
    RT.final_tier(rec)
    assert "cohort" not in rec.read


# ------------------------------------------------- I. multiple uncertainty reasons kept
def test_a_patient_can_carry_several_uncertainty_reasons_at_once():
    e = ev(n_primary=40, dn_sensitivity=0.05, dn_trunc10k=0.04, enr_globalbins=0.8,
           per_sample=[("s1", 30, 0.40), ("s2", 25, 0.02)])
    _, _, _, reasons = RT.final_tier(e)
    fired = [r for r, v in reasons.items() if v]
    assert len(fired) >= 4
    for r in ("uncertain_low_n", "uncertain_denominator", "uncertain_depth",
              "uncertain_repeated_sample", "uncertain_null_scheme"):
        assert reasons[r] is True


# ------------------------------------------ single-sample patients are never "agreeing"
def test_single_sample_patient_is_not_credited_with_agreement():
    assert RT.sample_agreement([("s1", 900, 0.4)], 0.25) == "not_assessable_for_agreement"


def test_sample_below_the_cell_floor_does_not_count_toward_agreement():
    assert RT.sample_agreement([("s1", 900, 0.4), ("s2", 5, 0.01)],
                               0.25) == "not_assessable_for_agreement"


def test_frozen_constants_are_unchanged():
    assert RT.TAU_LOW == 0.05
    assert RT.TAU_HIGH_SET == (0.20, 0.25, 0.33)
    assert RT.MIN_N_PRIMARY == 100 and RT.MIN_CELLS_PER_SAMPLE == 20
