"""End-to-end synthetic tests of the frozen Stage-10 pseudobulk design.

The Codex audit's MEDIUM finding #3: `test_e_depth_matched_indices_return_cells_not_
replicates` asserts an array shape, not that the production DE builds exactly one
pseudobulk per patient per group. The production code did not exist in the repo at the
time, so that test could not have been written. It exists now, recovered verbatim at
`production/stage10/s10e_pseudobulk_de_decoupler_tc.py`, so the design is testable.

These tests run the **recovered loop shape** over synthetic counts, using the same
committed primitives the driver calls (`SC.adaptive_depth_bins`,
`SC.depth_matched_indices`). They assert behaviour — replication unit, matching order,
antigen exclusion — never array shape alone. No frozen artifact is read or written.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import numpy as np
import pytest
import scipy.sparse as sp

from mm_escape import subclone as SC

REPO = Path(__file__).resolve().parents[1]
DRIVER = REPO / "production" / "stage10" / "s10e_pseudobulk_de_decoupler_tc.py"
SEED = 20260825


# --------------------------------------------------------------- synthetic cohort
def _synthetic_cohort(n_patients=6, cells_per_patient=200, n_genes=40, seed=0):
    """A cohort where DN cells are systematically shallower than comparator cells.

    That is the confound the frozen design exists to defeat, so the fixture reproduces
    it rather than a clean case.
    """
    rng = np.random.default_rng(seed)
    patient, sample, is_dn, depth, rows = [], [], [], [], []
    for p in range(n_patients):
        n_samples = 2 if p < 2 else 1          # two patients contribute repeated samples
        for c in range(cells_per_patient):
            dn = c % 4 == 0
            d = rng.integers(300, 1500) if dn else rng.integers(1200, 6000)
            patient.append(f"P{p}")
            sample.append(f"P{p}_s{(c // 3) % n_samples}")   # DN cells span both samples
            is_dn.append(dn)
            depth.append(float(d))
            base = rng.poisson(d / n_genes, n_genes).astype(float)
            if dn:
                base[:5] *= 0.5                # a genuine DN-lower program
            rows.append(base)
    return (np.array(patient), np.array(sample), np.array(is_dn, dtype=bool),
            np.array(depth), sp.csr_matrix(np.vstack(rows)),
            [f"G{i}" for i in range(n_genes)])


def _build_pseudobulks(patient, is_dn, depth, counts, patients, seed=SEED):
    """The recovered driver's loop, verbatim in shape.

    for each patient: adaptive depth bins -> depth-matched indices -> ONE summed
    pseudobulk vector per group. Matching happens on cells, before any summation.
    """
    dn_mat, pos_mat, used = [], [], {}
    for p in patients:
        idx = np.flatnonzero(patient == p)
        bins = SC.adaptive_depth_bins(depth[idx])
        di, pi = SC.depth_matched_indices(is_dn[idx], bins, seed)
        used[p] = (idx[di], idx[pi])
        dn_mat.append(np.asarray(counts[idx[di]].sum(axis=0)).ravel())
        pos_mat.append(np.asarray(counts[idx[pi]].sum(axis=0)).ravel())
    return np.vstack(dn_mat), np.vstack(pos_mat), used


@pytest.fixture(scope="module")
def cohort():
    return _synthetic_cohort()


# =================================================== A. one pseudobulk per patient/group
def test_a1_each_patient_contributes_exactly_one_pseudobulk_per_group(cohort):
    patient, _sample, is_dn, depth, counts, _genes = cohort
    pats = sorted(set(patient))
    dn_mat, pos_mat, _ = _build_pseudobulks(patient, is_dn, depth, counts, pats)
    assert dn_mat.shape[0] == len(pats)
    assert pos_mat.shape[0] == len(pats)
    assert dn_mat.shape[0] == pos_mat.shape[0], "groups must stay paired over patients"


def test_a2_the_de_replication_unit_is_the_patient_not_the_cell(cohort):
    """n for the paired test is the patient count — thousands of cells must not inflate it."""
    patient, _sample, is_dn, depth, counts, _genes = cohort
    pats = sorted(set(patient))
    dn_mat, pos_mat, used = _build_pseudobulks(patient, is_dn, depth, counts, pats)
    n_matched_cells = sum(len(d) + len(p) for d, p in used.values())
    assert n_matched_cells > 10 * len(pats), "fixture is too small to distinguish the two"
    # the matrix the test statistic sees has one row per patient, not per cell
    assert dn_mat.shape[0] == len(pats) < n_matched_cells


def test_a3_adding_cells_to_one_patient_does_not_add_a_replicate(cohort):
    patient, _sample, is_dn, depth, counts, _genes = cohort
    pats = sorted(set(patient))
    before, _, _ = _build_pseudobulks(patient, is_dn, depth, counts, pats)
    # duplicate every cell of P0 — more cells, same patient
    extra = np.flatnonzero(patient == "P0")
    patient2 = np.concatenate([patient, patient[extra]])
    is_dn2 = np.concatenate([is_dn, is_dn[extra]])
    depth2 = np.concatenate([depth, depth[extra]])
    counts2 = sp.vstack([counts, counts[extra]]).tocsr()
    after, _, _ = _build_pseudobulks(patient2, is_dn2, depth2, counts2, pats)
    assert after.shape[0] == before.shape[0] == len(pats)


# ============================================ B. repeated samples are not extra patients
def test_b1_repeated_samples_collapse_into_their_patient(cohort):
    """Two patients here carry two samples each. The design must still see 6 patients."""
    patient, sample, is_dn, depth, counts, _genes = cohort
    assert len(set(sample)) > len(set(patient)), "fixture must contain repeated samples"
    pats = sorted(set(patient))
    dn_mat, _, _ = _build_pseudobulks(patient, is_dn, depth, counts, pats)
    assert dn_mat.shape[0] == len(set(patient)) == 6
    assert dn_mat.shape[0] < len(set(sample))


def test_b2_a_multi_sample_patient_pools_its_samples_into_one_row(cohort):
    patient, sample, is_dn, depth, counts, _genes = cohort
    pats = sorted(set(patient))
    _, _, used = _build_pseudobulks(patient, is_dn, depth, counts, pats)
    dn_idx, _ = used["P0"]
    assert len(set(sample[dn_idx])) > 1, "P0's pseudobulk must draw from both its samples"


def test_b3_the_driver_iterates_patients_not_samples():
    src = DRIVER.read_text()
    assert "for p in pats" in src
    assert "obs.patient_id==p" in src.replace(" ", "")
    assert "sample_name" not in src.split("pseudobulk DE")[-1].split("decoupler")[0]


# ================================================= C. depth matching precedes pseudobulk
def test_c1_matching_happens_on_cells_before_any_summation(cohort):
    """If summation came first, the two groups' depths could not be equalised at all."""
    patient, _sample, is_dn, depth, counts, _genes = cohort
    pats = sorted(set(patient))
    _, _, used = _build_pseudobulks(patient, is_dn, depth, counts, pats)
    for p in pats:
        dn_idx, pos_idx = used[p]
        assert len(dn_idx) == len(pos_idx), f"{p}: groups not matched 1:1"
        assert set(dn_idx).isdisjoint(set(pos_idx))


def test_c2_matching_removes_the_depth_gap_that_exists_before_it(cohort):
    patient, _sample, is_dn, depth, counts, _genes = cohort
    pats = sorted(set(patient))
    _, _, used = _build_pseudobulks(patient, is_dn, depth, counts, pats)
    raw = np.median(depth[is_dn]) / np.median(depth[~is_dn])
    dn_all = np.concatenate([used[p][0] for p in pats])
    pos_all = np.concatenate([used[p][1] for p in pats])
    matched = np.median(depth[dn_all]) / np.median(depth[pos_all])
    assert raw < 0.8, "fixture must start confounded"
    assert abs(matched - 1.0) < abs(raw - 1.0), "matching must shrink the depth gap"
    assert 0.7 < matched < 1.4


def test_c3_the_driver_calls_bins_then_match_then_sum_in_that_order():
    src = DRIVER.read_text()
    i_bins = src.index("SC.adaptive_depth_bins")
    i_match = src.index("SC.depth_matched_indices")
    i_sum = src.index("sum(axis=0)")
    assert i_bins < i_match < i_sum, "order is bins -> match -> pseudobulk"


def test_c4_a_patient_with_no_matchable_bin_contributes_nothing_rather_than_a_biased_row():
    """Zero matched cells must not become a silently unmatched pseudobulk."""
    is_dn = np.ones(30, dtype=bool)            # no comparator cells at all
    bins = SC.adaptive_depth_bins(np.linspace(100, 5000, 30))
    di, pi = SC.depth_matched_indices(is_dn, bins, SEED)
    assert di.size == 0 and pi.size == 0


# ================================================ D. antigens never become DE features
def test_d1_antigen_columns_are_dropped_before_normalization(cohort):
    _p, _s, _dn, _d, counts, genes = cohort
    genes = list(genes)
    genes[3], genes[7] = "TNFRSF17", "GPRC5D"
    mat, kept = SC.drop_antigen_features(counts, genes)
    assert "TNFRSF17" not in kept and "GPRC5D" not in kept
    assert mat.shape[1] == counts.shape[1] - 2
    assert len(kept) == mat.shape[1]


def test_d2_dropping_happens_before_the_size_factor_is_formed(cohort):
    """Order is the whole point: normalise first and the antigen leaks into every gene."""
    _p, _s, _dn, _d, counts, genes = cohort
    genes = list(genes)
    genes[3] = "TNFRSF17"
    dense = counts.toarray().astype(float)
    spiked = dense.copy()
    spiked[:, 3] += 5000.0                      # a huge antigen perturbation

    def cpm_after_drop(mat):
        m, _ = SC.drop_antigen_features(sp.csr_matrix(mat), genes)
        m = m.toarray()
        return 1e6 * m / m.sum(axis=1, keepdims=True)

    np.testing.assert_allclose(cpm_after_drop(dense), cpm_after_drop(spiked))

    def cpm_before_drop(mat):
        cpm = 1e6 * mat / mat.sum(axis=1, keepdims=True)
        m, _ = SC.drop_antigen_features(sp.csr_matrix(cpm), genes)
        return m.toarray()

    assert not np.allclose(cpm_before_drop(dense), cpm_before_drop(spiked)), (
        "normalising first must visibly contaminate every gene — that is why order matters")


def test_d3_no_predeclared_level2_program_contains_an_antigen():
    from mm_escape import config as CFG

    for name in CFG.LEVEL2_PROGRAMS:
        genes = set(CFG.STATE_PROGRAMS[name])
        assert not genes & set(SC.ANTIGEN_FEATURES), f"{name} contains an antigen"


def test_d4_the_stage10_chain_drops_antigens_before_anything_else():
    for step in ("s10a_level1_structure_and_per_patient_de.py",
                 "s10c_program_scores_and_cohort_tests.py"):
        src = (REPO / "production" / "stage10" / step).read_text()
        assert "antigens gone" in src.lower(), step


# ================================================= E. determinism of the recovered design
def test_e1_the_same_seed_gives_the_same_pseudobulks(cohort):
    patient, _s, is_dn, depth, counts, _g = cohort
    pats = sorted(set(patient))
    a, b, _ = _build_pseudobulks(patient, is_dn, depth, counts, pats, seed=SEED)
    c, d, _ = _build_pseudobulks(patient, is_dn, depth, counts, pats, seed=SEED)
    np.testing.assert_array_equal(a, c)
    np.testing.assert_array_equal(b, d)


def test_e2_the_driver_fixes_its_seed():
    assert re.search(r"SEED\s*=\s*\d+", DRIVER.read_text())
    assert "SC.depth_matched_indices(dn,bins,SEED)" in DRIVER.read_text().replace(" ", "")


def test_e3_the_driver_still_parses():
    ast.parse(DRIVER.read_text(), filename=str(DRIVER))
