"""Stage-09 invariants A-G. The load-bearing ones are A (joint-DN isolation) and
B (frozen Stage-08 outputs cannot be touched)."""
import hashlib
import inspect
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mm_escape import bulk as BK

S09 = Path("results/09_bulk_validation")
S08 = Path("results/08_dual_antigen_escape")
FORBIDDEN = ("double_negative", "frac_double_negative", "conegativity", "co_negativity",
             "enrichment", "dual_negative")

requires_stage09 = pytest.mark.skipif(
    not (S09 / "sc_bulk_sample_mapping.csv").exists(), reason="stage 09 outputs absent")


# --------------------------------------------------------------- A. joint-DN isolation
def test_no_stage09_function_takes_a_joint_dn_argument():
    """Bulk cannot validate the joint state, so no entry point may even accept it."""
    for name in BK.__all__:
        obj = getattr(BK, name)
        if not callable(obj):
            continue
        for param in inspect.signature(obj).parameters:
            assert not any(f in param.lower() for f in FORBIDDEN), (
                f"{name}() takes a joint-DN parameter '{param}'")


def test_bulk_module_source_never_references_joint_dn_quantities():
    src = inspect.getsource(BK)
    code = "\n".join(l for l in src.splitlines()
                     if not l.strip().startswith(("#", '"', "'")))
    for f in ("frac_double_negative", "co_negativity", "conegativity"):
        assert f not in code, f"bulk.py computes on '{f}'"


@requires_stage09
def test_stage09_outputs_contain_no_joint_dn_column():
    for csv in S09.glob("*.csv"):
        cols = " ".join(pd.read_csv(csv, nrows=0).columns).lower()
        for f in ("double_negative", "dual_negative", "conegativity"):
            assert f not in cols, f"{csv.name} carries a joint-DN column"


# ------------------------------------------------------- B. frozen Stage-08 is untouched
@requires_stage09
def test_stage09_functions_cannot_modify_frozen_stage08_outputs():
    def digest():
        return {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                for p in sorted(S08.glob("*.csv"))}

    before = digest()
    rng = np.random.default_rng(0)
    g, t = rng.integers(0, 5, 200), rng.integers(500, 5000, 200)
    BK.pseudobulk_cpm(g, t, np.ones(200, bool))
    BK.spearman_or_not_evaluable(rng.random(10), rng.random(10))
    assert digest() == before, "a stage-09 call changed a frozen stage-08 file"


# --------------------------------------------------------- C. both denominators present
@requires_stage09
def test_both_denominators_are_reported_independently():
    assert (S09 / "sc_marginal_antigen_primary.csv").exists()
    assert (S09 / "sc_marginal_antigen_sensitivity.csv").exists()
    r = pd.read_csv(S09 / "bulk_vs_sc_by_cohort.csv")
    assert {"primary", "sensitivity"} <= set(r.denominator)
    for antigen in ("TNFRSF17", "GPRC5D"):
        sub = r[(r.antigen == antigen) & (r.stratum == "pooled")]
        assert {"primary", "sensitivity"} <= set(sub.denominator), antigen


# ------------------------------------------------------------- D. antigens kept separate
@requires_stage09
def test_antigens_are_analysed_separately_and_never_combined():
    a = pd.read_csv(S09 / "bulk_vs_sc_tnfrsf17.csv")
    b = pd.read_csv(S09 / "bulk_vs_sc_gprc5d.csv")
    assert set(a.antigen) == {"TNFRSF17"} and set(b.antigen) == {"GPRC5D"}
    cols = " ".join(pd.read_csv(S09 / "sc_marginal_antigen_primary.csv", nrows=0).columns)
    assert "combined" not in cols.lower() and "score" not in cols.lower()


# ------------------------------------------------------- E. no bulk -> cell imputation
def test_pseudobulk_cpm_returns_a_scalar_not_per_cell_values():
    rng = np.random.default_rng(1)
    g, t = rng.integers(0, 4, 500), rng.integers(1000, 9000, 500)
    out = BK.pseudobulk_cpm(g, t, np.ones(500, bool))
    assert isinstance(out, float) and not isinstance(out, np.ndarray)


def test_no_stage09_function_returns_an_array_sized_like_cells():
    """Nothing may hand back a per-cell vector that a caller could write onto obs."""
    rng = np.random.default_rng(2)
    g, t = rng.integers(0, 4, 300), rng.integers(1000, 9000, 300)
    for out in (BK.pseudobulk_cpm(g, t, np.ones(300, bool)),
                BK.spearman_or_not_evaluable(rng.random(20), rng.random(20))):
        assert not isinstance(out, np.ndarray)


# --------------------------------------------------------------- F. patient-level unit
def test_select_one_pair_per_patient_yields_exactly_one_row_per_patient():
    m = pd.DataFrame({"sc_patient": ["A", "A", "A", "B"],
                      "sc_sample": ["A_1", "A_2", "A_4", "B_1"],
                      "timepoint": ["1", "2", "4", "1"]})
    counts = {"A_1": 500, "A_2": 3, "A_4": 90, "B_1": 400}
    out = BK.select_one_pair_per_patient(m, counts)
    assert len(out) == out.sc_patient.nunique() == 2
    # earliest timepoint clearing the floor, not the largest and not the best-correlating
    assert out.set_index("sc_patient").loc["A", "sc_sample"] == "A_1"


def test_select_one_pair_per_patient_drops_samples_below_the_cell_floor():
    m = pd.DataFrame({"sc_patient": ["A", "A"], "sc_sample": ["A_1", "A_2"],
                      "timepoint": ["1", "2"]})
    out = BK.select_one_pair_per_patient(m, {"A_1": 5, "A_2": 900})
    assert list(out.sc_sample) == ["A_2"]


@requires_stage09
def test_correlation_input_has_one_observation_per_patient():
    d = pd.read_csv(S09 / "sc_marginal_antigen_primary.csv")
    assert d.sc_patient.is_unique, "a patient contributes more than one observation"


@requires_stage09
def test_ambiguous_matches_are_not_evaluable_rather_than_inferred():
    m = pd.read_csv(S09 / "sc_bulk_sample_mapping.csv")
    for bulk_id in ("47499", "98433", "59114_2"):
        row = m[m.bulk_sample == bulk_id]
        assert len(row) == 1 and row.match_status.iloc[0] == BK.NOT_EVALUABLE
        assert pd.isna(row.sc_patient.iloc[0]), f"{bulk_id} was silently matched"


# --------------------------------------------------------------- G. normal-donor unit
@requires_stage09
def test_normal_marrow_summary_preserves_donor_identity():
    n = pd.read_csv(S09 / "normal_marrow_antigen_context.csv")
    assert "donor" in n.columns and n.donor.nunique() > 1
    assert n.groupby(["donor", "population"]).size().max() == 1
    assert (n.n_cells > 0).all()


# --------------------------------------------------------------- small-n discipline
def test_spearman_reports_not_evaluable_rather_than_no_relationship():
    out = BK.spearman_or_not_evaluable([1, 2, 3], [3, 2, 1])
    assert out["status"].startswith(BK.NOT_EVALUABLE) and np.isnan(out["rho"])


def test_multi_run_bulk_is_averaged_not_summed():
    """MMRF_1686 stacks two runs; TPM is already per-run normalised."""
    src = inspect.getsource(BK.read_mmrf_bulk)
    assert ".mean()" in src and "per[g].sum()" not in src
