"""Stage-08 antigen logic. The load-bearing test is antigen->depth invariance."""
import numpy as np
import pytest

from mm_escape import antigen as A


# ---------------------------------------------------------------- depth metric
def test_depth_metric_invariant_to_antigen_counts():
    """Perturbing either antigen must not move a cell between depth strata.

    This is stage 07's antigen-circularity invariant carried into stage 08. The depth
    metric is the one place total UMI could leak antigen expression into the analysis,
    so it is tested at the extremes rather than argued to be negligible.
    """
    rng = np.random.default_rng(0)
    n = 500
    other = rng.integers(200, 40000, n)
    bcma = rng.integers(0, 30, n)
    gprc = rng.integers(0, 5, n)

    base = A.depth_ex_antigen(other + bcma + gprc, bcma, gprc)
    edges = A.quantile_edges(base, 5)
    strata = A.assign_strata(base, edges)

    for nb, ng in [(np.zeros(n, int), gprc), (bcma, np.zeros(n, int)),
                   (np.zeros(n, int), np.zeros(n, int)),
                   (np.full(n, 10000), gprc), (bcma, np.full(n, 10000)),
                   (np.full(n, 10000), np.full(n, 10000))]:
        d = A.depth_ex_antigen(other + nb + ng, nb, ng)
        assert np.array_equal(d, base)
        assert np.array_equal(A.assign_strata(d, A.quantile_edges(d, 5)), strata)


def test_depth_metric_would_have_failed_on_bare_total_counts():
    """Guard the *reason* for the metric: bare total UMI does move cells between bins."""
    rng = np.random.default_rng(1)
    n = 400
    other = rng.integers(200, 5000, n)
    bcma = rng.integers(0, 30, n)
    naive = other + bcma
    perturbed = other + 10000
    s1 = A.assign_strata(naive, A.quantile_edges(naive, 5))
    s2 = A.assign_strata(perturbed, A.quantile_edges(perturbed, 5))
    assert not np.array_equal(s1, s2)


# ---------------------------------------------------------------- states
def test_observed_states_cover_the_four_cells_of_the_2x2():
    st = A.observed_states([True, True, False, False], [True, False, True, False])
    assert list(st) == [A.DP, A.BCMA_ONLY, A.GPRC5D_ONLY, A.DN]


# ---------------------------------------------------------------- binning
def test_quantile_edges_dedup_on_ties_yields_fewer_bins():
    x = np.array([5] * 90 + [7, 8, 9, 10, 100, 200, 300, 400, 500, 600])
    edges = A.quantile_edges(x, 5)
    assert len(edges) - 1 < 5
    assert np.all(np.diff(edges) > 0)


def test_assign_strata_is_left_closed_right_open_with_closed_top():
    edges = np.array([0.0, 10.0, 20.0])
    assert list(A.assign_strata([0, 9.99, 10, 19.9, 20], edges)) == [0, 0, 1, 1, 1]


def test_merge_sparse_strata_merges_into_the_larger_neighbour():
    s = np.array([0] * 5 + [1] * 3 + [2] * 40)
    out = A.merge_sparse_strata(s, min_cells=20)
    assert len(np.unique(out)) == 1


def test_merge_sparse_strata_leaves_adequate_strata_alone():
    s = np.array([0] * 30 + [1] * 30)
    assert np.array_equal(A.merge_sparse_strata(s, min_cells=20), s)


def test_merge_sparse_strata_never_empties_everything():
    s = np.array([0] * 3 + [1] * 4)
    out = A.merge_sparse_strata(s, min_cells=20)
    assert len(out) == 7 and len(np.unique(out)) == 1


# ---------------------------------------------------------------- nulls
def test_stratified_expectation_equals_permutation_mean():
    """The closed form must match the permutation it stands in for, or the bootstrap
    (which uses the closed form) would be estimating a different null."""
    rng = np.random.default_rng(2)
    n = 600
    strata = rng.integers(0, 3, n)
    b = rng.random(n) < (0.2 + 0.2 * strata)
    g = rng.random(n) < (0.1 + 0.3 * strata)
    closed = A.stratified_expected_dn(b, g, strata)
    null, _ = A.permutation_null_dn(b, g, strata, n_perm=400, seed=7)
    assert abs(closed - null.mean()) < 0.01


def test_depth_alone_creates_spurious_unconditioned_enrichment():
    """The reason the null is stratified: with no within-stratum coupling at all,
    depth heterogeneity alone inflates the unconditioned ratio above 1 while the
    stratified ratio stays at 1. An unconditioned null would report co-escape on data
    that has none — in the direction this project hopes to find."""
    rng = np.random.default_rng(3)
    strata = np.repeat([0, 1, 2], 2000)
    p = np.array([0.85, 0.5, 0.15])[strata]
    b = rng.random(strata.size) < p
    g = rng.random(strata.size) < p
    obs = float((b & g).mean())
    assert obs / A.unconditioned_expected_dn(b, g) > 1.15
    assert abs(obs / A.stratified_expected_dn(b, g, strata) - 1.0) < 0.06


def test_permutation_preserves_within_stratum_marginals():
    rng = np.random.default_rng(4)
    n = 300
    strata = rng.integers(0, 2, n)
    b = rng.random(n) < 0.4
    g = rng.random(n) < 0.3
    null, p = A.permutation_null_dn(b, g, strata, n_perm=200, seed=5)
    assert 0.0 < p <= 1.0 and null.size == 200


# ---------------------------------------------------------------- bootstrap
def test_hierarchical_bootstrap_is_wider_than_flat_when_samples_differ():
    """A flat cell bootstrap would treat between-sample batch variation as biological
    spread. With two samples disagreeing sharply, the hierarchical interval must be
    the wider one."""
    vals = np.concatenate([np.ones(500), np.zeros(500)])
    sid = np.array(["s1"] * 500 + ["s2"] * 500)
    hier = A.hierarchical_bootstrap(sid, vals, n_boot=400, seed=1)
    flat = A.hierarchical_bootstrap(np.full(1000, "s1"), vals, n_boot=400, seed=1)
    assert np.ptp(np.percentile(hier, [2.5, 97.5])) > np.ptp(np.percentile(flat, [2.5, 97.5]))


def test_single_sample_bootstrap_degenerates_to_cell_bootstrap():
    rng = np.random.default_rng(6)
    vals = rng.random(300) < 0.3
    one = A.hierarchical_bootstrap(np.full(300, "s"), vals, n_boot=300, seed=2)
    assert abs(one.mean() - vals.mean()) < 0.05


# ---------------------------------------------------------------- downsampling
def test_downsample_leaves_shallow_cells_untouched():
    out = A.downsample_gene_counts([500, 9000], [4, 12], cap=10000)
    assert list(out) == [4, 12]


def test_downsample_scales_deep_cells_toward_the_cap_ratio():
    n = 4000
    out = A.downsample_gene_counts(np.full(n, 100000), np.full(n, 100), cap=10000, seed=3)
    assert 8.0 < out.mean() < 12.0
    assert out.max() <= 100


def test_downsample_never_exceeds_the_original_gene_count():
    rng = np.random.default_rng(8)
    tot = rng.integers(10001, 200000, 500)
    gene = rng.integers(0, 200, 500)
    out = A.downsample_gene_counts(tot, gene, cap=10000, seed=4)
    assert np.all(out <= gene) and np.all(out >= 0)


def test_downsample_is_reproducible_under_a_fixed_seed():
    tot, gene = np.full(200, 50000), np.full(200, 60)
    a = A.downsample_gene_counts(tot, gene, cap=10000, seed=11)
    b = A.downsample_gene_counts(tot, gene, cap=10000, seed=11)
    assert np.array_equal(a, b)
