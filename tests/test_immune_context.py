"""Stage-11 falsification tests, frozen and passing BEFORE real associations were read.

Stages 08 and 10 both manufactured strong apparent biology out of sequencing depth. These
tests assert the Stage-11 model does not repeat that.
"""
import inspect

import numpy as np
import pytest

from mm_escape import communication as CM

N_PAT = 32
COHORTS = np.array(["MMRF"] * 10 + ["WU1"] * 11 + ["WU2"] * 11)


def scenario(seed=0):
    """A cohort/depth structure resembling the real one: MMRF far deeper."""
    rng = np.random.default_rng(seed)
    depth = np.where(COHORTS == "MMRF", rng.normal(25000, 3000, N_PAT),
                     rng.normal(5800, 800, N_PAT))
    n_immune = rng.integers(200, 6000, N_PAT).astype(float)
    n_samples = rng.integers(1, 3, N_PAT).astype(float)
    # DN burden also tracks cohort, as it really does -> predictor and confound entangled
    dn = np.where(COHORTS == "MMRF", rng.normal(0.22, 0.06, N_PAT),
                  rng.normal(0.33, 0.08, N_PAT)).clip(0.01, 0.95)
    conf = [COHORTS, np.log10(depth), np.log10(n_immune), np.log10(n_samples)]
    return dn, conf, depth, rng


# ------------------------------------------------------- A. synthetic cohort confound
def test_a_feature_driven_only_by_cohort_is_not_reported_as_a_dn_association():
    """The predictor is entangled with cohort, so an unadjusted model is fooled and the
    adjusted one must not be."""
    dn, conf, depth, rng = scenario(1)
    feature = np.where(COHORTS == "MMRF", 2.0, -1.0) + rng.normal(0, 0.05, N_PAT)

    unadjusted = CM.ols_association(feature, dn, ())
    adjusted = CM.ols_association(feature, dn, conf)

    assert unadjusted["p"] < 0.05, "unadjusted model should be fooled by cohort"
    assert adjusted["p"] > 0.05, "cohort-adjusted model must not report it"


# -------------------------------------------------------- B. synthetic depth confound
def test_b_feature_driven_only_by_depth_does_not_survive_depth_adjustment():
    dn, conf, depth, rng = scenario(2)
    feature = 3.0 * np.log10(depth) + rng.normal(0, 0.02, N_PAT)

    unadjusted = CM.ols_association(feature, dn, ())
    adjusted = CM.ols_association(feature, dn, conf)

    assert unadjusted["p"] < 0.05
    assert adjusted["p"] > 0.05


def test_b2_depth_only_feature_is_flagged_by_the_depth_pre_test():
    """Each feature is screened against depth alone before its DN association is read."""
    _, _, depth, rng = scenario(3)
    feature = 3.0 * np.log10(depth) + rng.normal(0, 0.02, N_PAT)
    assert CM.ols_association(feature, np.log10(depth), ())["p"] < 1e-6


# ------------------------------------------------------------ C. permuted patient labels
def test_c_permuting_dn_labels_keeps_the_false_positive_rate_at_nominal():
    dn, conf, depth, rng = scenario(4)
    feature = rng.normal(0, 1, N_PAT)
    hits = 0
    trials = 400
    for i in range(trials):
        permuted = np.random.default_rng(1000 + i).permutation(dn)
        if CM.ols_association(feature, permuted, conf)["p"] < 0.05:
            hits += 1
    assert 0.01 < hits / trials < 0.12, f"false-positive rate {hits / trials:.3f}"


def test_c2_permutation_also_controlled_for_a_categorical_coherence_predictor():
    dn, conf, depth, rng = scenario(5)
    feature = rng.normal(0, 1, N_PAT)
    coherence = np.zeros(N_PAT)
    coherence[:4] = 1.0                       # 4 supported, as observed
    hits = sum(
        CM.ols_association(feature, np.random.default_rng(2000 + i).permutation(coherence),
                           conf)["p"] < 0.05
        for i in range(400))
    assert hits / 400 < 0.12


# ------------------------------------------------------------- real signal is detectable
def test_the_adjusted_model_still_detects_a_genuine_dn_association():
    """The controls must not have neutered the model."""
    dn, conf, depth, rng = scenario(6)
    feature = 8.0 * dn + rng.normal(0, 0.2, N_PAT)
    out = CM.ols_association(feature, dn, conf)
    assert out["p"] < 0.01 and out["ci_lo"] > 0


# ------------------------------------------------------------------ tier isolation
def test_changing_provisional_tier_labels_cannot_alter_a_stage11_result():
    dn, conf, depth, rng = scenario(7)
    feature = rng.normal(0, 1, N_PAT)
    base = CM.ols_association(feature, dn, conf)
    for labels in (np.array(["robust-high"] * N_PAT), np.array(["uncertain"] * N_PAT),
                   np.array(["robust-low"] * N_PAT)):
        assert CM.ols_association(feature, dn, conf) == base
        assert labels.shape == (N_PAT,)          # labels exist but never enter the model


def test_module_never_references_a_provisional_tier_or_an_antigen_gene():
    code = "\n".join(l for l in inspect.getsource(CM).splitlines()
                     if not l.strip().startswith(("#", '"', "'")))
    for token in ("robust-high", "robust_high", "robust-low", "TAU_HIGH"):
        assert token not in code
    assert "GPRC5D" not in code
    # TNFRSF17 may appear once, only inside the fixed external LR receptor list
    assert code.count("TNFRSF17") <= 1


def test_antigens_are_not_immune_features():
    assert not (set(CM.IMMUNE_CATEGORIES) & {"TNFRSF17", "GPRC5D"})


# --------------------------------------------------------------------- compositional
def test_clr_removes_the_sum_to_one_constraint():
    frac = np.array([0.5, 0.3, 0.15, 0.05])
    out = CM.clr(frac)
    assert abs(out.sum()) < 1e-9
    doubled = CM.clr(frac * 2 / (frac * 2).sum())
    assert np.allclose(out, doubled)           # scale-invariant, as a composition must be


def test_clr_handles_a_zero_lineage_without_blowing_up():
    assert np.all(np.isfinite(CM.clr(np.array([0.9, 0.1, 0.0, 0.0]))))


# ---------------------------------------------------------------- reporting discipline
def test_within_cohort_reports_not_evaluable_rather_than_no_relationship():
    f = np.arange(8.0)
    p = np.arange(8.0)
    coh = np.array(["A"] * 6 + ["B"] * 2)
    out = CM.within_cohort_spearman(f, p, coh)
    assert out["A"]["status"] == "evaluable"
    assert out["B"]["status"].startswith(CM.NOT_EVALUABLE)
    assert np.isnan(out["B"]["rho"])


def test_benjamini_hochberg_is_monotone_and_bounded():
    p = np.array([0.001, 0.01, 0.02, 0.3, 0.9])
    adj = CM.benjamini_hochberg(p)
    assert np.all(np.diff(adj) >= -1e-12) and adj.max() <= 1.0 and np.all(adj >= p - 1e-12)


def test_ols_is_not_evaluable_at_tiny_n_rather_than_returning_a_number():
    out = CM.ols_association([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], ())
    assert out["status"] == CM.NOT_EVALUABLE


def test_frozen_stage11_constants():
    assert CM.MIN_IMMUNE_CELLS == 100
    assert CM.CONFOUNDERS == ("cohort", "log_depth", "log_n_immune", "log_n_samples")
    assert len(CM.MEASUREMENT_PREDICTORS) == 4
    assert "cytotoxic_mixed" in CM.IMMUNE_CATEGORIES and "NK_core" in CM.IMMUNE_CATEGORIES


# ============================================================================
# Stage-11 resume (2026-08-26): the extraction helpers, and the invariants that
# make the recomputed stage-11 tables checkable against the preserved first run.
# ============================================================================

from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
STAGE11 = REPO / "results" / "11_immune_context"
INTEGRATED = REPO / "results" / "05_integration" / "integrated.h5ad"

requires_stage11 = pytest.mark.skipif(
    not (STAGE11 / "patient_immune_composition.csv").exists(),
    reason="needs the stage-11 outputs (run notebooks/11_immune_context.py)",
)
requires_integrated = pytest.mark.skipif(
    not INTEGRATED.exists(), reason="needs results/05_integration/integrated.h5ad",
)


def test_pseudobulk_is_pooled_not_a_mean_of_per_cell_rates():
    """A mean over cells weights a 300-UMI cell like a 20,000-UMI one.

    In this cohort that means weighting a WashU cell like an MMRF cell, which is the
    depth confound the whole stage is built to avoid.
    """
    counts = pd.DataFrame({"G": [1.0, 1.0], "total_counts": [100.0, 9900.0]})
    pooled = CM.pseudobulk_cpm(counts, ["G"])["G"]
    per_cell_mean = float(np.mean(counts.G / counts.total_counts * 1e6))
    assert pooled == pytest.approx(2 / 10000 * 1e6)
    assert pooled < per_cell_mean / 2, "pooling must not reduce to the per-cell mean"


def test_pseudobulk_returns_nan_for_an_empty_group_not_zero():
    """Absent is not the same as zero — a patient with no cells of a lineage has no
    estimate, and must not be entered into a regression as a measured zero."""
    empty = pd.DataFrame({"G": [], "total_counts": []})
    assert np.isnan(CM.pseudobulk_cpm(empty, ["G"])["G"])


def test_lr_candidate_set_never_names_the_antigens_as_discovered_features():
    receptors = {r for _, r in CM.LR_CANDIDATES if r != "None"}
    ligands = {l for l, _ in CM.LR_CANDIDATES}
    assert "GPRC5D" not in receptors | ligands
    # TNFRSF17 may appear only as a fixed external receptor name, never as a ligand
    assert "TNFRSF17" not in ligands


def test_sender_floor_is_the_frozen_stage08_group_minimum():
    """Reused, not re-derived — a new constant here would be a new free parameter."""
    from mm_escape import subclone

    assert CM.MIN_SENDER_CELLS == subclone.MIN_GROUP_CELLS == 20
    assert CM.MIN_IMMUNE_CELLS == subclone.MIN_PATIENT_CELLS == 100


@requires_integrated
def test_stream_gene_counts_totals_match_the_stage08_depth_definition():
    """The row sum must be over the *intersected* gene space.

    `obs["total_counts"]` was computed at QC time over each sample's full Cell Ranger
    reference, before stage 05 intersected to 32,991 genes, so it runs a few counts
    high. Using it shifts per-patient median depth by up to 14 counts — small, but the
    depth covariate is doing most of the work in this stage's models.
    """
    frame = CM.stream_gene_counts(INTEGRATED, ["TNFRSF17", "GPRC5D"])
    assert len(frame) == 172_940
    assert (frame.total_counts > 0).all()
    assert frame[["TNFRSF17", "GPRC5D"]].to_numpy().min() >= 0
    assert (frame.total_counts >= frame[["TNFRSF17", "GPRC5D"]].sum(axis=1)).all()


@requires_integrated
def test_stream_gene_counts_raises_on_a_missing_gene_rather_than_returning_zeros():
    """A silently-zero column would read as biological absence — the exact confusion
    this project exists to avoid."""
    with pytest.raises(KeyError):
        CM.stream_gene_counts(INTEGRATED, ["TNFRSF17", "NOT_A_REAL_GENE"])


@requires_stage11
def test_stage11_recomputation_matches_the_preserved_preliminary_run():
    """The paused first run is preserved, and the notebook's recomputation reproduces it.

    This is what turns 'preliminary outputs on disk' into 'a reproducible stage'.
    """
    prelim = STAGE11 / "preliminary_run" / "patient_immune_composition.csv"
    if not prelim.exists():
        pytest.skip("preliminary run not preserved in this tree")
    a = pd.read_csv(prelim, dtype={"patient": str}).set_index("patient").sort_index()
    b = pd.read_csv(STAGE11 / "patient_immune_composition.csv",
                    dtype={"patient": str}).set_index("patient").sort_index()
    assert list(a.index) == list(b.index)
    numeric = [c for c in a.columns if a[c].dtype.kind in "if" and c in b.columns]
    assert len(numeric) >= 20
    for col in numeric:
        assert np.allclose(a[col], b[col], atol=1e-4), col


@requires_stage11
def test_every_stage11_patient_clears_the_immune_evaluability_floor():
    comp = pd.read_csv(STAGE11 / "patient_immune_composition.csv", dtype={"patient": str})
    assert len(comp) == 32
    assert (comp.n_immune_cells >= CM.MIN_IMMUNE_CELLS).all()
    fractions = comp[[f"frac_{c}" for c in CM.IMMUNE_CATEGORIES]].sum(axis=1)
    assert np.allclose(fractions, 1.0), "one denominator, all lineages"


@requires_stage11
def test_no_composition_association_is_reported_as_surviving_correction():
    """The stage's headline is a negative one, and a future edit must not quietly
    turn it positive."""
    m = pd.read_csv(STAGE11 / "immune_vs_dn_measurement.csv")
    assert len(m) == 28
    assert (m.p_adj_BH >= 0.10).all(), "stage 11 reports nothing at BH < 0.10"


@requires_stage11
def test_the_communication_hit_is_matched_by_a_receiver_side_confound():
    """`Tcell PDCD1 -> CD274` is not a candidate immune axis: the receptor term is
    measured on the same plasma cells whose DN status is the predictor, and most
    receptors move down together."""
    confound = pd.read_csv(STAGE11 / "communication_receiver_side_confound.csv")
    clone = confound[confound.receiver == "clone_primary"]
    assert (clone.coef_adjusted < 0).sum() >= 10, "broad receptor-side shift"
    assert "CD274" in set(clone[clone.p_adjusted < 0.05].receptor)


@requires_stage11
def test_both_receiver_definitions_are_kept_on_disk():
    """The design-conformant run and the deviating first run are both preserved, so the
    amendment can be audited rather than taken on trust."""
    for name in ("communication_context.csv", "communication_context_all_plasma.csv",
                 "communication_context_vs_dn.csv",
                 "communication_context_vs_dn_all_plasma.csv"):
        assert (STAGE11 / name).exists(), name
