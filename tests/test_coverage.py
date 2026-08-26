"""Supplemental multi-antigen coverage — invariants, frozen BEFORE any coverage result.

The arithmetic here is simple enough to look obviously correct and is exactly the kind
that goes quietly wrong: a mis-ordered set operation makes a triple look better than its
own pairs, and nothing downstream would notice. Every relation the analysis relies on is
asserted rather than assumed.

The structural guards matter more than the arithmetic ones. This project has been misled
four separate times by naive depth handling, so the rule is that there is exactly one
depth-stratification implementation and this module calls it.
"""
from __future__ import annotations

import hashlib
import inspect
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pytest

from mm_escape import antigen as AG
from mm_escape import coverage as CV
from mm_escape import subclone as SC

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "results" / "08_dual_antigen_escape" / "multi_antigen_coverage"

requires_coverage = pytest.mark.skipif(
    not (OUT / "single_target_coverage.csv").exists(),
    reason="needs the supplemental coverage outputs (run notebooks/08c_multi_antigen_coverage.py)",
)


def random_detection(n_cells=400, n_targets=7, seed=0):
    rng = np.random.default_rng(seed)
    # heterogeneous per-target detection rates, as the real panel has
    rates = rng.uniform(0.05, 0.9, n_targets)
    return rng.random((n_cells, n_targets)) < rates


# --------------------------------------------------------------- A. monotonicity
def test_a_pair_is_never_less_covered_than_either_constituent_single():
    det = random_detection(seed=1)
    for a, b in combinations(CV.TARGETS, 2):
        pair = CV.uncovered_fraction(det, CV.TARGETS, (a, b))
        for single in (a, b):
            assert pair <= CV.uncovered_fraction(det, CV.TARGETS, (single,)) + 1e-12


def test_a_triple_is_never_less_covered_than_any_contained_pair():
    det = random_detection(seed=2)
    for combo in combinations(CV.TARGETS, 3):
        triple = CV.uncovered_fraction(det, CV.TARGETS, combo)
        for pair in combinations(combo, 2):
            assert triple <= CV.uncovered_fraction(det, CV.TARGETS, pair) + 1e-12


def test_a_monotonicity_holds_on_a_degenerate_all_negative_panel():
    """The edge case that would hide an off-by-one: nothing is ever detected, so every
    combination is fully uncovered and no inequality may be violated."""
    det = np.zeros((50, len(CV.TARGETS)), dtype=bool)
    for k in (1, 2, 3):
        for combo in combinations(CV.TARGETS, k):
            assert CV.uncovered_fraction(det, CV.TARGETS, combo) == 1.0


# ------------------------------------------------------------- B. incremental gain
def test_b_adding_a_target_can_never_reduce_coverage():
    det = random_detection(seed=3)
    for combo in list(combinations(CV.TARGETS, 2)) + list(combinations(CV.TARGETS, 3)):
        for row in CV.incremental_gains(det, CV.TARGETS, combo):
            assert row["gain"] >= -1e-12, row


def test_b_gain_is_directional_and_both_directions_are_reported():
    det = random_detection(seed=4)
    rows = CV.incremental_gains(det, CV.TARGETS, CV.ANCHOR)
    assert {r["added"] for r in rows} == set(CV.ANCHOR)
    # a pair with very different marginals must show very different directional gains
    det2 = np.zeros((100, len(CV.TARGETS)), dtype=bool)
    det2[:90, 0] = True          # TNFRSF17 detected in 90%
    det2[:10, 1] = True          # GPRC5D detected in 10%
    g = {r["added"]: r["gain"] for r in CV.incremental_gains(det2, CV.TARGETS, CV.ANCHOR)}
    assert g["GPRC5D"] != pytest.approx(g["TNFRSF17"])


def test_b_gain_equals_the_definition_in_the_frozen_design():
    det = random_detection(seed=5)
    a, b = CV.ANCHOR
    p_a = CV.uncovered_fraction(det, CV.TARGETS, (a,))
    p_ab = CV.uncovered_fraction(det, CV.TARGETS, (a, b))
    rows = {r["added"]: r for r in CV.incremental_gains(det, CV.TARGETS, CV.ANCHOR)}
    assert rows["GPRC5D"]["gain"] == pytest.approx(p_a - p_ab)


# ------------------------------------------------------------- D. antigen isolation
def test_d_perturbing_one_target_changes_only_combinations_containing_it():
    det = random_detection(seed=6)
    perturbed = det.copy()
    j = CV.TARGETS.index("SLAMF7")
    perturbed[:, j] = True                      # drive SLAMF7 fully positive
    for k in (1, 2, 3):
        for combo in combinations(CV.TARGETS, k):
            before = CV.uncovered_fraction(det, CV.TARGETS, combo)
            after = CV.uncovered_fraction(perturbed, CV.TARGETS, combo)
            if "SLAMF7" in combo:
                assert after == 0.0
            else:
                assert after == pytest.approx(before), combo


def test_d_the_anchor_is_untouched_by_any_new_target():
    """The frozen BCMA/GPRC5D result must be invariant to everything added here."""
    det = random_detection(seed=7)
    anchor_before = CV.uncovered_fraction(det, CV.TARGETS, CV.ANCHOR)
    for name in ("SLAMF7", "FCRL5", "CD38", "SDC1", "ITGB7"):
        p = det.copy()
        p[:, CV.TARGETS.index(name)] = ~p[:, CV.TARGETS.index(name)]
        assert CV.uncovered_fraction(p, CV.TARGETS, CV.ANCHOR) == anchor_before


# ------------------------------------------------------------------ E. eligibility
def test_e_a_not_evaluable_target_cannot_be_offered_as_a_combination():
    eligible = [t for t in CV.TARGETS if t != "CD38"]
    combos = CV.all_combinations(eligible)
    assert all("CD38" not in c for c in combos)
    assert CV.ANCHOR in combos


def test_e_eligibility_reads_only_measurement_qc_and_never_coverage():
    """Coverage is not even a parameter — the rule cannot see how good a target looks."""
    params = set(inspect.signature(CV.eligibility).parameters)
    assert params == {"ambient_status", "malignant_detection", "background_detection",
                      "technical_zero", "circularity_blocked"}
    assert not (params & {"uncovered", "gain", "coverage", "rank"})


def test_e_dropout_prone_target_is_excluded_with_a_stated_reason():
    state, reason = CV.eligibility(ambient_status="ok", malignant_detection=0.4,
                                   background_detection=0.01, technical_zero=0.72)
    assert state == CV.NOT_EVALUABLE and "technical-zero" in reason


def test_e_unclean_background_only_blocks_when_separation_also_fails():
    """No clean negative population is survivable if the target is plainly detected far
    above its imperfect background; it is not if the two are comparable."""
    ok, _ = CV.eligibility(ambient_status=CV.AMBIENT_NOT_EVALUABLE,
                           malignant_detection=0.80, background_detection=0.05,
                           technical_zero=0.2)
    bad, reason = CV.eligibility(ambient_status=CV.AMBIENT_NOT_EVALUABLE,
                                 malignant_detection=0.30, background_detection=0.25,
                                 technical_zero=0.2)
    assert ok == CV.ELIGIBLE
    assert bad == CV.NOT_EVALUABLE and "background" in reason


# ------------------------------------------------------------------- F. raw counts
def test_f_positivity_rejects_normalised_or_imputed_input():
    """Whether a transcript is genuinely absent is the entire question; smoothing over the
    zeros erases the measurement."""
    with pytest.raises(ValueError):
        CV.detected(np.array([0.0, 0.37, 1.82]))          # CP10K-like
    assert CV.detected(np.array([0, 1, 5])).tolist() == [False, True, True]


def test_f_positivity_is_the_frozen_stage08_rule_for_every_target():
    counts = np.array([0, 1, 2, 9])
    assert CV.detected(counts).tolist() == [False, True, True, True]
    # count >= 2 is the declared stage-08 *sensitivity* rule and is not the primary here
    assert CV.detected(np.array([1]))[0]


# --------------------------------------------------------------- G. no utility score
def test_g_no_weighted_aggregate_exists_anywhere_in_the_module():
    """Checked against executable code with docstrings stripped — the module *discusses*
    why there is no utility score, and that prose must not be what makes this pass."""
    import ast

    tree = ast.parse(inspect.getsource(CV))
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef)):
            body = node.body
            if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant) \
                    and isinstance(body[0].value.value, str):
                node.body = body[1:] or [ast.Pass()]
    code = ast.unparse(ast.fix_missing_locations(tree))
    code = "\n".join(l for l in code.splitlines() if not l.strip().startswith("#"))
    for token in ("lambda_", "utility", "weighted_score", "composite", "risk_score"):
        assert token not in code, token

    public = [n for n in dir(CV) if not n.startswith("_")]
    assert not [n for n in public if "utility" in n.lower() or "composite" in n.lower()]
    assert not [n for n in public if n.lower().endswith("score")]


@requires_coverage
def test_g_no_output_column_combines_the_evidence_dimensions():
    import pandas as pd

    for f in OUT.glob("*.csv"):
        cols = [c.lower() for c in pd.read_csv(f, nrows=0).columns]
        for c in cols:
            assert "utility" not in c and "composite" not in c, (f.name, c)


# --------------------------------------------------------------------- H. patient unit
@requires_coverage
def test_h_primary_output_is_one_row_per_patient_denominator_combination():
    import pandas as pd

    p = pd.read_csv(OUT / "pair_coverage.csv", dtype={"patient": str})
    assert not p.duplicated(["patient", "denominator", "combination"]).any()
    assert set(p.denominator) == {"primary", "sensitivity"}


# --------------------------------------------------------- J. shared depth utility
def test_j_module_contains_no_local_depth_binning():
    """One implementation, called from everywhere. Four separate wrappers around the same
    idea is how a fifth subtly-different one gets written."""
    src = inspect.getsource(CV)
    code = "\n".join(l for l in src.splitlines() if not l.strip().startswith("#"))
    for banned in ("np.quantile", "np.percentile", "searchsorted", "np.histogram",
                   "pd.qcut", "pd.cut", "np.digitize", "linspace"):
        assert banned not in code, f"local binning primitive {banned!r} in coverage.py"


def test_j_strata_helper_delegates_to_the_stage08_utility():
    depth = np.concatenate([np.linspace(500, 3000, 150), np.linspace(3000, 40000, 150)])
    via_helper = CV.shared_depth_strata(depth)
    edges = AG.quantile_edges(depth, n_bins=5)
    direct = AG.merge_sparse_strata(AG.assign_strata(depth, edges), min_cells=20)
    assert np.array_equal(via_helper, direct)


def test_j_expected_co_negative_delegates_to_the_frozen_closed_form():
    rng = np.random.default_rng(11)
    a = rng.random(500) < 0.4
    b = rng.random(500) < 0.3
    s = CV.shared_depth_strata(rng.gamma(2, 3000, 500))
    assert CV.stratified_expected_co_negative(a, b, s) == AG.stratified_expected_dn(a, b, s)
    assert CV.unconditioned_expected_co_negative(a, b) == AG.unconditioned_expected_dn(a, b)


def test_j_the_stage10_adaptive_wrapper_is_still_the_only_other_one():
    """If a third binning wrapper appears, this fails and forces consolidation."""
    assert callable(SC.adaptive_depth_bins) and callable(AG.assign_strata)
    assert not hasattr(CV, "depth_bins") and not hasattr(CV, "quantile_edges")


# --------------------------------------------------- K. synthetic depth-only null
def _depth_only_pair(n=3000, seed=42):
    """Two targets whose zeros are produced by library size alone.

    Detection probability for each is a function of the cell's own depth and nothing else,
    so the pair is co-negative purely because shallow cells read zero for both. There is no
    biological co-loss to find, and a conditioned null must say so.
    """
    rng = np.random.default_rng(seed)
    depth = rng.lognormal(mean=np.log(4000), sigma=1.2, size=n)
    p = np.clip(depth / 12000, 0.02, 0.95)
    return (rng.random(n) > p), (rng.random(n) > p), depth


def test_k_depth_only_pair_looks_enriched_unconditioned():
    """The failure mode this test exists to catch must actually occur, or the guard below
    is vacuous."""
    a, b, _ = _depth_only_pair()
    observed = float((a & b).mean())
    expected = CV.unconditioned_expected_co_negative(a, b)
    assert observed / expected > 1.15, "synthetic scenario is not producing the artifact"


def test_k_depth_only_pair_does_not_survive_the_conditioned_null():
    a, b, depth = _depth_only_pair()
    strata = CV.shared_depth_strata(depth, n_bins=10)
    observed = float((a & b).mean())
    conditioned = CV.stratified_expected_co_negative(a, b, strata)
    assert abs(observed / conditioned - 1.0) < 0.05, (
        "a purely depth-generated pair must not read as biologically co-lost")


def test_k_conditioned_null_still_detects_a_genuine_co_loss():
    """The control must not have neutered the statistic."""
    rng = np.random.default_rng(7)
    n = 3000
    depth = rng.lognormal(mean=np.log(4000), sigma=0.9, size=n)
    shared = rng.random(n) < 0.25                     # a real co-negative subpopulation
    a = shared | (rng.random(n) < 0.10)
    b = shared | (rng.random(n) < 0.10)
    strata = CV.shared_depth_strata(depth, n_bins=10)
    observed = float((a & b).mean())
    assert observed / CV.stratified_expected_co_negative(a, b, strata) > 1.5


def test_k_permutation_null_is_the_frozen_stage08_one():
    a, b, depth = _depth_only_pair(n=800, seed=3)
    strata = CV.shared_depth_strata(depth, n_bins=10)
    _, p = AG.permutation_null_dn(a, b, strata, n_perm=200, seed=20260825)
    assert 0.0 < p <= 1.0


# ------------------------------------------------------------------ L. SDC1 safeguard
def test_l_sdc1_cannot_be_eligible_while_the_circularity_flag_is_set():
    state, reason = CV.eligibility(ambient_status="ok", malignant_detection=0.5,
                                   background_detection=0.01, technical_zero=0.1,
                                   circularity_blocked=True)
    assert state == CV.NOT_EVALUABLE and "circularity" in reason


def test_l_sdc1_is_recorded_as_a_stage06_plasma_identity_gene():
    """The checkpoint exists because this is true of the frozen pipeline, not because SDC1
    looked bad. `PLASMA_MATURE` is stage 06's axis-(b) mature-plasma predicate."""
    from mm_escape import config

    assert "SDC1" in config.PLASMA_MATURE
    assert "SDC1" in config.MARKER_PANEL["PlasmaCell"]
    assert "CD38" in config.MARKER_PANEL["PlasmaCell"]
    assert "CD38" not in config.PLASMA_MATURE      # deliberately excluded from axis (b)


@requires_coverage
def test_l_sdc1_qc_questions_were_answered_before_any_eligibility_call():
    import pandas as pd

    qc = pd.read_csv(OUT / "target_measurement_qc.csv")
    row = qc[qc.antigen == "SDC1"]
    assert len(row), "SDC1 must remain visible in the QC output whatever its verdict"
    assert row.coverage_eligibility.notna().all()
    assert row.sdc1_differentiation_caveat.notna().all()
    if (row.coverage_eligibility == CV.ELIGIBLE).any():
        assert row.sdc1_checkpoint_resolved.all(), (
            "SDC1 may only be eligible if its differentiation/depth checkpoint resolved")


# ---------------------------------------------------- I. frozen-stage isolation
FROZEN_ARTIFACTS = [
    "results/08_dual_antigen_escape/patient_antigen_states_primary.csv",
    "results/08_dual_antigen_escape/patient_antigen_states_sensitivity.csv",
    "results/08_dual_antigen_escape/patient_conegativity_enrichment.csv",
    "results/08_dual_antigen_escape/noise_floor_technical_zero.csv",
    "results/08_dual_antigen_escape/noise_floor_ambient.csv",
    "results/08_dual_antigen_escape/truncate10k_sensitivity.csv",
    "results/08_dual_antigen_escape/depth_strata_definition.csv",
    "results/09_bulk_validation/",
    "results/10_dn_coherence/dn_coherence_final_states.csv",
    "results/11_immune_context/patient_immune_composition.csv",
]


def _digest(path: Path) -> str:
    if path.is_dir():
        parts = [f"{p.relative_to(path)}:{_digest(p)}" for p in sorted(path.rglob("*"))
                 if p.is_file()]
        return hashlib.sha256("|".join(parts).encode()).hexdigest()
    return hashlib.sha256(path.read_bytes()).hexdigest()


@requires_coverage
def test_i_supplemental_analysis_wrote_nothing_into_a_frozen_namespace():
    """Every supplemental output lives under `multi_antigen_coverage/`, and the frozen
    artifacts carry the digests recorded when this analysis ran."""
    manifest_path = OUT / "frozen_upstream_digests.json"
    assert manifest_path.exists(), "the run must record the frozen digests it consumed"
    recorded = json.loads(manifest_path.read_text())
    for rel, digest in recorded.items():
        p = REPO / rel
        if not p.exists():
            pytest.skip(f"{rel} absent in this tree")
        assert _digest(p) == digest, f"FROZEN ARTIFACT CHANGED: {rel}"


@requires_coverage
def test_i_all_supplemental_outputs_are_in_the_separated_namespace():
    stage08 = REPO / "results" / "08_dual_antigen_escape"
    produced = {"single_target_coverage.csv", "pair_coverage.csv", "triple_coverage.csv",
                "incremental_gain_pairs.csv", "incremental_gain_triples.csv"}
    for name in produced:
        assert (OUT / name).exists()
        assert not (stage08 / name).exists(), f"{name} leaked into the stage-08 namespace"
