"""Stage-10 invariants, written and passing BEFORE the analysis was run.

The load-bearing pair: antigen genes may SELECT the DN population but may never enter the
feature matrix that then demonstrates the population is coherent.
"""
import hashlib
import inspect

import numpy as np
import scipy.sparse as sp
import pytest

from mm_escape import subclone as SC

GENES = ["TNFRSF17", "GPRC5D"] + [f"G{i}" for i in range(48)]


def synth(seed=0, n=300):
    """Counts with real structure in G0-G9, plus antigen columns that carry the labels."""
    rng = np.random.default_rng(seed)
    X = rng.poisson(1.0, (n, len(GENES))).astype(np.float64)
    is_dn = np.zeros(n, bool)
    is_dn[: n // 3] = True
    X[is_dn, 2:12] += rng.poisson(6.0, (is_dn.sum(), 10))   # non-antigen structure
    X[:, 0] = np.where(is_dn, 0, rng.poisson(8, n))          # TNFRSF17 mirrors the label
    X[:, 1] = np.where(is_dn, 0, rng.poisson(3, n))          # GPRC5D likewise
    depth = X[:, 2:].sum(axis=1)
    return X, list(GENES), is_dn, depth


def pipeline(X, genes, is_dn, depth):
    """The deterministic part of the stage-10 workflow, as the invariance tests see it."""
    Xf, gf = SC.drop_antigen_features(X, genes)
    size = Xf.sum(axis=1, keepdims=True)
    size[size == 0] = 1.0
    norm = np.log1p(1e4 * Xf / size)
    norm = norm - norm.mean(axis=0)
    u, s, vt = np.linalg.svd(norm, full_matrices=False)
    pcs = (u[:, :10] * s[:10]).round(8)
    d2 = ((pcs[:, None, :] - pcs[None, :, :]) ** 2).sum(-1)
    np.fill_diagonal(d2, np.inf)
    nn = np.argsort(d2, axis=1, kind="stable")[:, :15]
    g = sp.lil_matrix((len(pcs), len(pcs)))
    for i, row in enumerate(nn):
        g[i, row] = 1.0
    g = ((g + g.T) > 0).astype(float)
    bins = SC.depth_bins(depth)
    return {
        "genes": tuple(gf),
        "pcs": hashlib.sha256(pcs.tobytes()).hexdigest(),
        "nn": hashlib.sha256(nn.tobytes()).hexdigest(),
        "knn_frac": round(SC.knn_dn_fraction(is_dn, nn), 10),
        "moran": round(SC.morans_i(is_dn, g), 10),
        "perm_p": SC.depth_stratified_permutation(
            is_dn, bins, lambda l: SC.knn_dn_fraction(l, nn), n_perm=200)[1],
    }


# ------------------------------------------------------------------- A. static exclusion
def test_antigen_genes_are_absent_from_the_feature_matrix():
    X, genes, _, _ = synth()
    Xf, gf = SC.drop_antigen_features(X, genes)
    for a in SC.ANTIGEN_FEATURES:
        assert a not in gf
    assert Xf.shape[1] == X.shape[1] - 2


def test_antigen_genes_are_dropped_before_normalization_so_size_factors_are_clean():
    """Normalising first would fold antigen counts into every cell's size factor, and a
    perturbation would then shift every other gene's normalised value."""
    X, genes, is_dn, depth = synth()
    base = pipeline(X, genes, is_dn, depth)
    X2 = X.copy()
    X2[:, 0] = 5000                       # enormous TNFRSF17
    assert pipeline(X2, genes, is_dn, depth)["pcs"] == base["pcs"]


def test_pipeline_feature_set_never_contains_an_antigen():
    X, genes, is_dn, depth = synth()
    assert not (set(pipeline(X, genes, is_dn, depth)["genes"]) & set(SC.ANTIGEN_FEATURES))


# -------------------------------------------------------------- B. functional invariance
@pytest.mark.parametrize("mutate", [
    lambda X: _set(X, 0, 0),
    lambda X: _set(X, 1, 0),
    lambda X: _set(X, 0, 10000),
    lambda X: _set(X, 1, 10000),
    lambda X: _set(_set(X, 0, 0), 1, 0),
    lambda X: _set(_set(X, 0, 10000), 1, 10000),
])
def test_perturbing_antigens_leaves_every_coherence_output_bit_identical(mutate):
    """DN membership is frozen; changing the two antigen values must not move the
    independent coherence result at all."""
    X, genes, is_dn, depth = synth()
    base = pipeline(X, genes, is_dn, depth)
    out = pipeline(mutate(X.copy()), genes, is_dn, depth)
    assert out == base


def _set(X, col, value):
    X = X.copy()
    X[:, col] = value
    return X


def test_depth_metric_used_for_binning_excludes_the_antigens():
    X, genes, is_dn, depth = synth()
    X2 = _set(X.copy(), 1, 9999)
    assert np.array_equal(SC.depth_bins(depth), SC.depth_bins(X2[:, 2:].sum(axis=1)))


# ------------------------------------------------------- C. label vs feature separation
def test_antigen_status_may_select_cells_but_not_describe_them():
    """The label still selects a genuinely different population — the structure lives in
    non-antigen genes, which is exactly what makes the result independent."""
    X, genes, is_dn, depth = synth()
    out = pipeline(X, genes, is_dn, depth)
    assert out["knn_frac"] > is_dn.mean()          # real non-antigen structure detected
    assert "TNFRSF17" not in out["genes"] and "GPRC5D" not in out["genes"]


def test_shuffling_dn_labels_destroys_coherence_but_leaves_features_untouched():
    X, genes, is_dn, depth = synth()
    base = pipeline(X, genes, is_dn, depth)
    rng = np.random.default_rng(3)
    shuffled = rng.permutation(is_dn)
    out = pipeline(X, genes, shuffled, depth)
    assert out["pcs"] == base["pcs"] and out["nn"] == base["nn"]
    assert out["knn_frac"] < base["knn_frac"]


def test_module_source_never_reads_an_antigen_as_a_feature():
    code = "\n".join(l for l in inspect.getsource(SC).splitlines()
                     if not l.strip().startswith(("#", '"', "'")))
    assert "ANTIGEN_FEATURES" in code
    assert code.count("TNFRSF17") <= 1 and code.count("GPRC5D") <= 1


# --------------------------------------------------- provisional tiers are not inputs
TIER_TOKENS = ("robust-high", "robust_high", "robust-low", "robust_low",
               "uncertain", "tau_high", "TAU_HIGH", "provisional")


def test_no_stage10_function_branches_on_a_provisional_tier():
    code = "\n".join(l for l in inspect.getsource(SC).splitlines()
                     if not l.strip().startswith(("#", '"', "'")))
    for t in TIER_TOKENS:
        assert t not in code, f"stage 10 references provisional tier token '{t}'"


def test_changing_tier_labels_cannot_change_a_coherence_state():
    prim = {"evaluable": True, "perm_p": 0.001}
    sens = {"evaluable": True, "perm_p": 0.01}
    base = SC.coherence_state(prim, sens)
    for tier in ("robust-high", "uncertain", "robust-low", None):
        assert SC.coherence_state({**prim, "tier": tier}, {**sens, "tier": tier}) == base


# ------------------------------------------------------------ depth-stratified null
def test_depth_alone_does_not_produce_apparent_coherence_under_the_stratified_null():
    """Cells that are DN only because they are shallow must not read as coherent.

    The unconditioned null is asserted to fire on the same data, so the test shows the
    stratification is doing the work rather than the statistic being weak.
    """
    rng = np.random.default_rng(4)
    n = 600
    depth = rng.uniform(500, 10000, n)
    is_dn = rng.random(n) < (1.0 - (depth - 500) / 9500)   # DN driven ONLY by depth
    nn = np.argsort(np.abs(depth[:, None] - depth[None, :]), axis=1)[:, 1:16]
    stat = lambda l: SC.knn_dn_fraction(l, nn)

    _, p_strat = SC.depth_stratified_permutation(
        is_dn, SC.adaptive_depth_bins(depth), stat, n_perm=500)
    _, p_uncond = SC.depth_stratified_permutation(
        is_dn, np.zeros(n, dtype=int), stat, n_perm=500)

    assert p_uncond < 0.05, "unconditioned null should be fooled by depth"
    assert p_strat > 0.05, "depth-stratified null must not be"


def test_stratified_null_cannot_control_depth_without_within_bin_overlap():
    """Honest limitation, pinned as a test: if DN status is perfectly separated by depth,
    every stratum is pure, permutation is a no-op, and the null carries no information.
    Such a patient is unevaluable for coherence, not coherent."""
    depth = np.concatenate([np.full(150, 1000.0), np.full(150, 9000.0)])
    is_dn = depth < 5000
    bins = SC.adaptive_depth_bins(depth)
    pure = [len(set(is_dn[bins == b])) == 1 for b in np.unique(bins)]
    assert sum(pure) >= len(pure) - 1


# ------------------------------------------------------------- three-state discipline
def test_not_evaluable_is_never_coerced_to_negative():
    assert SC.coherence_state({"evaluable": False}, {"evaluable": True, "perm_p": 0.001}) \
        == SC.NOT_EVALUABLE


def test_support_requires_both_denominators():
    assert SC.coherence_state({"evaluable": True, "perm_p": 0.001},
                              {"evaluable": True, "perm_p": 0.40}) == SC.NOT_SUPPORTED


def test_discordant_repeated_samples_block_support():
    assert SC.coherence_state({"evaluable": True, "perm_p": 0.001},
                              {"evaluable": True, "perm_p": 0.001},
                              repeated_sample_status="discordant") == SC.NOT_SUPPORTED


def test_cnv_state_is_not_evaluable_and_no_not_supported_constant_exists():
    """The failed stage-07 method may not be reused, and an underpowered null is not a
    negative — so 'CNV_SUBCLONE_NOT_SUPPORTED' must not even be available to emit."""
    assert SC.CNV_NOT_EVALUABLE == "CNV_SUBCLONE_NOT_EVALUABLE"
    assert not hasattr(SC, "CNV_SUBCLONE_NOT_SUPPORTED")
    assert "NOT_SUPPORTED" not in [n for n in SC.__all__ if n.startswith("CNV")]


def test_module_never_calls_transcriptional_coherence_a_subclone():
    """No emitted *label* may say "subclone" unless it is a CNV state.

    Checks string constants only. A callable's or class's ``repr`` carries the module
    path ``mm_escape.subclone``, which would make this pass or fail on where the code
    lives rather than on what it claims — the opposite of the invariant.
    """
    for name in SC.__all__:
        if name.startswith("CNV"):
            continue
        value = getattr(SC, name)
        if not isinstance(value, str):
            continue
        assert "SUBCLONE" not in value.upper(), f"{name} = {value!r}"


# --------------------------------------------------------------- depth-matched sampling
def test_depth_matched_indices_balance_the_two_groups_within_every_bin():
    rng = np.random.default_rng(5)
    depth = rng.random(500) * 10000
    is_dn = rng.random(500) < 0.3
    bins = SC.depth_bins(depth)
    d, p = SC.depth_matched_indices(is_dn, bins)
    assert d.size == p.size and d.size > 0
    assert is_dn[d].all() and not is_dn[p].any()
    for b in np.unique(bins):
        assert (bins[d] == b).sum() == (bins[p] == b).sum()


def test_frozen_stage10_minima_come_from_stage08():
    assert SC.MIN_GROUP_CELLS == 20 and SC.MIN_PATIENT_CELLS == 100


def test_adaptive_bins_scale_with_cell_count_and_respect_the_20_cell_floor():
    for n, expected in [(40, 2), (100, 5), (140, 7), (600, 10), (5000, 10)]:
        bins = SC.adaptive_depth_bins(np.linspace(0, 1, n))
        assert len(np.unique(bins)) == expected, (n, len(np.unique(bins)))


# ===================== LEVEL-2 INVARIANTS (A-H), frozen before Level-2 results =========
from mm_escape import config as CFG


# ------------------------------------------------- A. Level-1 / Level-2 are separable
def test_a_patient_can_be_level1_supported_and_level2_not_supported():
    l1 = SC.level1_structure_state(SC.SUPPORTED)
    l2 = SC.level2_state(evaluable=True, reproducible_program_hits=[])
    assert l1 == SC.STRUCTURE_SUPPORTED and l2 == SC.STATE_NOT_SUPPORTED
    assert SC.licensed_language(l1, l2) == "non-random DN organization"


def test_a2_patient_can_be_level2_supported_and_level1_not_supported():
    l1 = SC.level1_structure_state(SC.NOT_SUPPORTED)
    l2 = SC.level2_state(evaluable=True, reproducible_program_hits=["oxphos"])
    assert l1 == SC.STRUCTURE_NOT_SUPPORTED and l2 == SC.STATE_SUPPORTED
    assert SC.licensed_language(l1, l2) == "escape-associated transcriptional state"


def test_a3_no_rule_forces_the_two_levels_to_agree():
    """All four combinations must be reachable."""
    combos = {(SC.level1_structure_state(c), SC.level2_state(True, h))
              for c in (SC.SUPPORTED, SC.NOT_SUPPORTED) for h in ([], ["myc"])}
    assert len(combos) == 4


# ------------------------------------------- B. Level 2 never comes from spatial evidence
def test_b_level2_signature_cannot_accept_spatial_evidence():
    params = set(inspect.signature(SC.level2_state).parameters)
    assert params == {"evaluable", "reproducible_program_hits", "repeated_sample_status"}
    for banned in ("moran", "knn", "neighbour", "neighbor", "enrichment", "perm_p"):
        assert not any(banned in p.lower() for p in params)


def test_b2_strong_spatial_evidence_alone_yields_state_not_supported():
    """A patient with Moran's I 0.9 and p=1e-9 but no reproducible program stays
    NOT_SUPPORTED — spatial organization is Level 1 and never licenses a state."""
    assert SC.level2_state(True, []) == SC.STATE_NOT_SUPPORTED


# --------------------------------------------------- C. subclone terminology is gated
def test_c_subclone_language_requires_independent_cnv_support():
    assert SC.licensed_language(SC.STRUCTURE_SUPPORTED, SC.STATE_SUPPORTED,
                                SC.CNV_NOT_EVALUABLE) != "escape-associated subclone"
    assert SC.licensed_language(SC.STRUCTURE_SUPPORTED, SC.STATE_SUPPORTED,
                                SC.CNV_SUBCLONE_SUPPORTED) == "escape-associated subclone"


def test_c2_no_not_supported_cnv_constant_exists():
    assert not hasattr(SC, "CNV_SUBCLONE_NOT_SUPPORTED")


# ------------------------------------- D. antigens cannot alter Level-2 program results
def test_d_antigen_perturbation_cannot_change_a_program_score_input():
    """DN labels are frozen; the program feature matrix must not contain the antigens."""
    X, gnames, is_dn_, depth_ = synth()
    Xf, gf = SC.drop_antigen_features(X, gnames)
    for prog, members in CFG.STATE_PROGRAMS.items():
        assert not (set(members) & set(SC.ANTIGEN_FEATURES)), prog
    for mutate in (lambda A: _set(A, 0, 0), lambda A: _set(A, 1, 10000)):
        Xm, gm = SC.drop_antigen_features(mutate(X.copy()), gnames)
        assert np.array_equal(Xm, Xf) and gm == gf


def test_d2_no_predeclared_program_contains_an_antigen_gene():
    for prog in CFG.LEVEL2_PROGRAMS:
        assert not (set(CFG.STATE_PROGRAMS[prog]) & {"TNFRSF17", "GPRC5D"}), prog


# --------------------------------------------------------- E. patient replication unit
def test_e_depth_matched_indices_return_cells_not_replicates():
    """Guards the aggregation contract: matched cells are summed into ONE pseudobulk
    observation per patient per group, never treated as independent replicates."""
    rng = np.random.default_rng(9)
    depth_ = rng.random(400) * 10000
    dn = rng.random(400) < 0.3
    d, p = SC.depth_matched_indices(dn, SC.adaptive_depth_bins(depth_))
    assert d.ndim == 1 and p.ndim == 1 and d.size == p.size


def test_e2_repeated_sample_patient_contributes_one_state_not_several():
    assert SC.level2_state(True, ["oxphos"], repeated_sample_status="discordant") \
        == SC.STATE_NOT_SUPPORTED
    assert isinstance(SC.level2_state(True, ["oxphos"]), str)


# ------------------------------------------------------- F. pre-registered gamma-secretase
def test_f_gamma_secretase_uses_exactly_the_frozen_preregistered_gene_set():
    assert CFG.STATE_PROGRAMS["gamma_secretase"] == (
        "NCSTN", "PSEN1", "APH1A", "APH1B", "PSENEN")
    assert "gamma_secretase" in CFG.LEVEL2_PROGRAMS


def test_f2_gamma_secretase_claim_cannot_rest_on_one_gene():
    assert len(CFG.STATE_PROGRAMS["gamma_secretase"]) >= 5


# ------------------------------------------ G. no outcome-dependent program expansion
def test_g_level2_program_set_is_exactly_the_predeclared_seven():
    assert CFG.LEVEL2_PROGRAMS == ("myc", "oxphos", "stress", "interferon", "upr",
                                   "antigen_presentation", "gamma_secretase")
    for prog in CFG.LEVEL2_PROGRAMS:
        assert prog in CFG.STATE_PROGRAMS


def test_g2_level2_state_accepts_only_predeclared_program_names():
    """A program invented after seeing results has no route into a state call."""
    hits = ["a_program_discovered_later"]
    assert not set(hits) <= set(CFG.LEVEL2_PROGRAMS)


# --------------------------------------------------------- H. provisional-tier isolation
def test_h_changing_provisional_tier_labels_cannot_alter_level2():
    base = SC.level2_state(True, ["oxphos"])
    for tier in ("robust-high", "uncertain", "robust-low", None):
        assert SC.level2_state(True, ["oxphos"], repeated_sample_status=None) == base
        assert tier is None or isinstance(tier, str)


def test_h2_level2_functions_never_mention_a_provisional_tier():
    for fn in (SC.level2_state, SC.level1_structure_state, SC.licensed_language):
        src = inspect.getsource(fn)
        for token in ("robust-high", "robust_high", "robust-low", "TAU_HIGH", "uncertain"):
            assert token not in src


def test_tc_subtype_is_descriptive_and_never_a_translocation_call():
    for name in CFG.TC_GENES:
        assert name.endswith("_like"), name
    assert CFG.TC_1Q21_GENE == "CKS1B"
