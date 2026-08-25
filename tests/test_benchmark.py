"""Stage 05b — the benchmark's arithmetic and its decision rule.

Almost all of this is data-free by design. The parts most worth protecting are the
depth statistic (which must be rotation-invariant, or it ranks methods on an accident
of their parameterisation) and the decision rule (which must be able to reject every
candidate, or "benchmark" means "pick something new"). Neither needs a matrix.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from mm_escape import benchmark


# ---------------------------------------------------------------------------
# Arm configuration
# ---------------------------------------------------------------------------

def test_arms_share_one_common_batch_key_plus_reference_arms():
    by_name = {arm.name: arm for arm in benchmark.ARMS}
    common = [a for a in benchmark.ARMS if a.batch_keys == ("sample_name",)]
    assert {a.method for a in common} == {"harmony", "scvi", "scanorama"}
    assert by_name["unintegrated"].batch_keys == ()
    # The incumbent keeps its own multi-covariate configuration, as shipped.
    assert by_name[benchmark.INCUMBENT].batch_keys == (
        "patient_id", "n_genes_ref", "cohort"
    )
    # cohort is the only batch definition not confounded with patient (42 of 50
    # patients contribute exactly one sample), so it gets real arms.
    assert {a.name for a in benchmark.ARMS if a.batch_keys == ("cohort",)} == {
        "harmony_cohort", "scvi_cohort"
    }


def test_bbknn_is_absent_on_purpose():
    """It yields a graph, not an embedding, so Benchmarker cannot score it."""
    assert "bbknn" not in {arm.method for arm in benchmark.ARMS}
    assert "bbknn" not in benchmark.RUNNERS


def test_plasma_is_not_in_the_scored_immune_compartment():
    """Plasma cells are diagnosed, never optimized — see the module docstring."""
    assert "plasma" not in benchmark.IMMUNE_CLASSES
    assert set(benchmark.IMMUNE_CLASSES) == {
        "T", "NK", "B", "myeloid", "erythroid", "HSPC"
    }


def test_run_arm_rejects_an_unknown_method():
    from anndata import AnnData

    adata = AnnData(X=np.zeros((5, 0), dtype=np.float32))
    arm = benchmark.ArmSpec("bogus", "not_a_method", ("cohort",))
    with pytest.raises(ValueError, match="unknown method"):
        benchmark.run_arm(adata, arm, pca=np.zeros((5, 2)))


# ---------------------------------------------------------------------------
# depth_association — the statistic fixed before running
# ---------------------------------------------------------------------------

def test_depth_association_is_one_when_depth_lies_in_the_span():
    rng = np.random.default_rng(0)
    embedding = rng.normal(size=(500, 5))
    depth = np.expm1(embedding[:, 0] * 0.5 + 3.0)
    assert benchmark.depth_association(embedding, depth) == pytest.approx(1.0, abs=1e-6)


def test_depth_association_is_near_zero_for_an_unrelated_embedding():
    rng = np.random.default_rng(1)
    embedding = rng.normal(size=(2000, 5))
    depth = np.expm1(rng.normal(size=2000) + 5.0)
    assert benchmark.depth_association(embedding, depth) < 0.05


def test_depth_association_is_rotation_invariant():
    """The whole reason R² was chosen over a per-dimension correlation.

    Latent axes are arbitrary and differ between methods, so a statistic that changes
    under rotation would rank methods on their parameterisation rather than on how much
    depth information they carry.
    """
    rng = np.random.default_rng(2)
    embedding = rng.normal(size=(800, 6))
    depth = np.expm1(embedding @ rng.normal(size=6) * 0.3 + 4.0)

    q, _ = np.linalg.qr(rng.normal(size=(6, 6)))
    rotated = embedding @ q

    assert benchmark.depth_association(embedding, depth) == pytest.approx(
        benchmark.depth_association(rotated, depth), abs=1e-8
    )

    # A per-dimension statistic is NOT invariant — this is what we avoided.
    def max_abs_corr(x, y):
        return max(abs(np.corrcoef(x[:, k], np.log1p(y))[0, 1]) for k in range(x.shape[1]))

    assert max_abs_corr(embedding, depth) != pytest.approx(
        max_abs_corr(rotated, depth), abs=1e-3
    )


def test_depth_association_handles_degenerate_input():
    rng = np.random.default_rng(3)
    embedding = rng.normal(size=(50, 3))
    assert benchmark.depth_association(embedding, np.full(50, 1000.0)) == 0.0


def test_depth_association_rejects_mismatched_lengths():
    with pytest.raises(ValueError, match="rows"):
        benchmark.depth_association(np.zeros((10, 2)), np.zeros(9))


# ---------------------------------------------------------------------------
# The decision rule
# ---------------------------------------------------------------------------

def _scores(**overrides) -> pd.DataFrame:
    base = {
        benchmark.INCUMBENT: dict(batch_score=0.50, bio_score=0.70,
                                  depth_r2=0.20, plasma_mixing=0.11),
        "candidate": dict(batch_score=0.60, bio_score=0.70,
                          depth_r2=0.20, plasma_mixing=0.11),
    }
    base["candidate"].update(overrides)
    return pd.DataFrame(base).T


def test_a_clean_improvement_is_eligible():
    decision = benchmark.decide(_scores())
    assert decision.loc["candidate", "eligible"]


def test_the_incumbent_can_never_replace_itself():
    decision = benchmark.decide(_scores())
    assert not decision.loc[benchmark.INCUMBENT, "eligible"]


def test_better_mixing_bought_with_biology_is_rejected():
    decision = benchmark.decide(_scores(bio_score=0.60))
    assert decision.loc["candidate", "batch_improved"]
    assert not decision.loc["candidate", "bio_preserved"]
    assert not decision.loc["candidate", "eligible"]


def test_better_mixing_bought_with_depth_encoding_is_rejected():
    decision = benchmark.decide(_scores(depth_r2=0.40))
    assert not decision.loc["candidate", "depth_ok"]
    assert not decision.loc["candidate", "eligible"]


def test_merely_matching_the_incumbent_is_not_an_improvement():
    decision = benchmark.decide(_scores(batch_score=0.50))
    assert not decision.loc["candidate", "batch_improved"]
    assert not decision.loc["candidate", "eligible"]


def test_a_plasma_mixing_jump_alone_only_flags():
    """A method could improve plasma geometry legitimately, so a jump is not fatal."""
    decision = benchmark.decide(_scores(plasma_mixing=0.80))
    assert decision.loc["candidate", "plasma_flagged"]
    assert decision.loc["candidate", "overcorrection_ok"]
    assert decision.loc["candidate", "eligible"]


def test_a_plasma_jump_with_rising_depth_association_is_disqualifying():
    """The signature of the censoring being smoothed over rather than respected."""
    decision = benchmark.decide(_scores(plasma_mixing=0.80, depth_r2=0.22))
    assert decision.loc["candidate", "plasma_flagged"]
    assert decision.loc["candidate", "depth_ok"]        # rise is within tolerance
    assert not decision.loc["candidate", "overcorrection_ok"]
    assert not decision.loc["candidate", "eligible"]


def test_decide_validates_its_inputs():
    with pytest.raises(ValueError, match="missing"):
        benchmark.decide(pd.DataFrame({"batch_score": [1.0]}, index=["a"]))
    frame = _scores().drop(index=benchmark.INCUMBENT)
    with pytest.raises(ValueError, match="incumbent"):
        benchmark.decide(frame)


# ---------------------------------------------------------------------------
# Rendering — the write-up is generated from the table, never hand-written
# ---------------------------------------------------------------------------

def test_render_reports_the_incumbent_surviving_as_a_real_result():
    decision = benchmark.decide(_scores(batch_score=0.50))
    text = benchmark.render_decision(decision)
    assert "No arm qualified" in text
    assert benchmark.INCUMBENT in text
    # The caveat that must never be dropped.
    assert "truncate-all-cohorts-at-10,000" in text
    assert "restores cells that were never deposited" in text


def test_render_names_the_winner_when_one_qualifies():
    scores = _scores()
    scores.loc["other", :] = [0.55, 0.70, 0.20, 0.11]
    decision = benchmark.decide(scores)
    text = benchmark.render_decision(decision)
    # Both qualify; the stronger batch score wins.
    assert "**candidate** replaces the incumbent" in text


def test_render_states_a_subsample_when_one_was_used():
    decision = benchmark.decide(_scores())
    assert "40,000" in benchmark.render_decision(decision, subsample=40_000)


# ---------------------------------------------------------------------------
# The provisional-label mapping — verified against real predictions when written
# ---------------------------------------------------------------------------

def test_ilc_maps_to_nk():
    """`Immune_All_High` folds NK into its ILC class.

    Checked on all 172,940 cells rather than assumed: the ILC call is 8% of the marrow
    (implausible for genuine innate lymphoid cells) at NKG7 98.7%, GNLY 92.2%,
    KLRD1 85.9%, MS4A1 1.3%.
    """
    assert benchmark.CELLTYPIST_TO_BROAD["ILC"] == "NK"
    assert "NK cells" not in benchmark.CELLTYPIST_TO_BROAD


def test_mapping_covers_every_broad_class_including_the_expected_blind_spots():
    """CLAUDE.md expected an immune-only reference to miss erythroid and HSPC.

    On this cohort it does not — Erythroid 14,103 cells at HBB 99.7%, HSC/MPP 2,625 at
    CD34 58.4% — so no hand-set marker override is needed and no arbitrary threshold
    enters the benchmark.
    """
    produced = set(benchmark.CELLTYPIST_TO_BROAD.values())
    assert produced == {*benchmark.IMMUNE_CLASSES, "plasma"}
    assert benchmark.CELLTYPIST_TO_BROAD["Erythroid"] == "erythroid"
    assert benchmark.CELLTYPIST_TO_BROAD["HSC/MPP"] == "HSPC"


def test_unmapped_labels_are_excluded_rather_than_forced_into_a_class():
    """Megakaryocytic and non-haematopoietic calls fit none of the seven classes."""
    for label in ("Megakaryocytes/platelets", "Fibroblasts", "Endothelial cells",
                  "Epithelial cells", "Cycling cells"):
        assert label not in benchmark.CELLTYPIST_TO_BROAD


def test_plasma_never_lands_in_the_immune_compartment():
    """The compartment split is what keeps plasma cells out of the scored set."""
    assert benchmark.CELLTYPIST_TO_BROAD["Plasma cells"] == "plasma"
    assert "plasma" not in benchmark.IMMUNE_CLASSES
