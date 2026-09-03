"""Stage-12 synthesis guards.

Stage 12 is a **synthesis stage**: it consumes 29 frozen artifacts, writes only into
`results/12_final_synthesis/`, and may not create a score, a ranking, a classifier or a new
statistical test. These tests enforce that, plus the frozen counts every Stage-12 output
inherits.

Data-backed tests skip cleanly when Stage 12 has not been run, so the suite stays runnable
on a fresh clone.
"""
from __future__ import annotations

import ast
import csv
import hashlib
import re
from pathlib import Path

import pandas as pd
import pytest

from mm_escape import subclone as SC
from mm_escape import synthesis as SY

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "results" / "12_final_synthesis"
MANIFEST = REPO / "provenance" / "frozen_artifacts_pre_stage12.tsv"
PRODUCER = REPO / "notebooks" / "12_final_synthesis.py"
MODULE = REPO / "src" / "mm_escape" / "synthesis.py"

requires_stage12 = pytest.mark.skipif(
    not (OUT / "stage12_patient_evidence_matrix.csv").exists(),
    reason="Stage 12 has not been executed")

TIER_EXPECT = {"robust-high": 4, "uncertain": 28}
L1_EXPECT = {"SUPPORTED": 4, "NOT_SUPPORTED": 23, "NOT_EVALUABLE": 5}
L2_EXPECT = {"SUPPORTED": 26, "NOT_SUPPORTED": 1, "NOT_EVALUABLE": 5}
JOINT_EXPECT = {
    ("uncertain", "NOT_SUPPORTED", "SUPPORTED"): 18,
    ("uncertain", "NOT_EVALUABLE", "NOT_EVALUABLE"): 5,
    ("robust-high", "NOT_SUPPORTED", "SUPPORTED"): 4,
    ("uncertain", "SUPPORTED", "SUPPORTED"): 4,
    ("uncertain", "NOT_SUPPORTED", "NOT_SUPPORTED"): 1,
}


def _sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 22), b""):
            h.update(c)
    return h.hexdigest()


def _manifest():
    with open(MANIFEST) as fh:
        return {r["path"]: r for r in csv.DictReader(fh, delimiter="\t")}


def _matrix():
    d = pd.read_csv(OUT / "stage12_patient_evidence_matrix.csv", dtype={"patient_id": str})
    d["L1"] = d.level1_state.str.replace("DN_STRUCTURE_", "", regex=False)
    d["L2"] = d.level2_state.str.replace("DN_STATE_", "", regex=False)
    return d


# =============================================================== A. provenance
@requires_stage12
def test_a1_input_manifest_lists_29_inputs():
    m = pd.read_csv(OUT / "stage12_input_manifest.tsv", sep="\t")
    assert len(m) == 29
    for col in ("relative_path", "upstream_stage", "sha256", "file_size", "producer",
                "producer_commit", "environment", "synthesis_role"):
        assert col in m.columns


@requires_stage12
def test_a2_every_recorded_input_hash_matches_the_committed_manifest():
    man = _manifest()
    m = pd.read_csv(OUT / "stage12_input_manifest.tsv", sep="\t")
    for _, r in m.iterrows():
        assert r.relative_path in man, r.relative_path
        assert r.sha256 == man[r.relative_path]["sha256"], r.relative_path
        assert str(r.file_size) == man[r.relative_path]["bytes"], r.relative_path


@requires_stage12
def test_a3_recorded_input_hashes_still_match_disk():
    m = pd.read_csv(OUT / "stage12_input_manifest.tsv", sep="\t")
    for _, r in m.iterrows():
        p = REPO / r.relative_path
        assert p.exists(), r.relative_path
        assert _sha(p) == r.sha256, f"{r.relative_path} mutated since Stage 12 ran"


@requires_stage12
def test_a4_design_snapshot_records_the_binding_design_hash():
    txt = (OUT / "stage12_design_snapshot.md").read_text()
    assert _sha(REPO / "docs" / "stage12_design.md") in txt
    assert "notebooks/12_final_synthesis.py" in txt


# ======================================================= B. namespace isolation
@requires_stage12
def test_b1_stage12_wrote_only_into_its_own_namespace():
    producer = PRODUCER.read_text()
    for m in re.finditer(r"""to_csv\(\s*([^)]*)\)""", producer):
        arg = m.group(1)
        if "results" in arg:
            assert "OUT" in arg, f"write outside the Stage-12 namespace: {arg[:80]}"
    assert 'OUT = REPO / "results" / "12_final_synthesis"' in producer


@requires_stage12
def test_b2_no_upstream_artifact_changed():
    """The whole 393-row freeze must still verify after Stage 12."""
    changed = [p for p, r in _manifest().items()
               if not (REPO / p).exists() or _sha(REPO / p) != r["sha256"]]
    assert not changed, f"upstream artifacts mutated: {changed[:5]}"


def test_b3_producer_never_writes_into_an_upstream_stage_directory():
    src = PRODUCER.read_text()
    for stage in ("04_qc", "05_integration", "06_annotation", "07_malignant_plasma",
                  "09_bulk_validation", "10_dn_coherence", "11_immune_context"):
        for m in re.finditer(rf"""(to_csv|write_text|savefig|open)\([^)]*{stage}""", src):
            raise AssertionError(f"producer writes into {stage}: {m.group(0)}")


# ============================================================= C. frozen counts
@requires_stage12
def test_c1_patient_rows():
    d = _matrix()
    assert len(d) == 32
    assert d.patient_id.is_unique
    # the invariant is that IDs stayed STRINGS -- silent int64 coercion once produced an
    # empty join. pandas 3 gives StringDtype for dtype=str, pandas 2 gives object, so
    # assert the behaviour rather than the dtype name.
    assert not pd.api.types.is_numeric_dtype(d.patient_id)
    assert all(isinstance(v, str) for v in d.patient_id)
    assert "25183" in set(d.patient_id) and "MMRF_1267" in set(d.patient_id)


@requires_stage12
def test_c2_tier_counts_unchanged():
    assert _matrix().provisional_measurement_tier.value_counts().to_dict() == TIER_EXPECT


@requires_stage12
def test_c3_level1_counts_unchanged():
    assert _matrix().L1.value_counts().to_dict() == L1_EXPECT


@requires_stage12
def test_c4_level2_counts_unchanged():
    assert _matrix().L2.value_counts().to_dict() == L2_EXPECT


@requires_stage12
def test_c5_level3_all_not_evaluable():
    d = _matrix()
    assert (d.level3_state == "CNV_SUBCLONE_NOT_EVALUABLE").all()
    assert (~d.cnv_evaluable.astype(bool)).all()


@requires_stage12
def test_c6_repeated_primary_patients_are_seven_not_eight():
    """60359 has zero primary-denominator cells and is not a repeated patient here."""
    r = pd.read_csv(OUT / "stage12_repeated_patient_summary.csv", dtype={"patient_id": str})
    assert len(r) == 7
    assert "60359" not in set(r.patient_id)


# ================================================================== D. cross-tab
@requires_stage12
def test_d1_joint_contingency_matches_the_design():
    d = _matrix()
    got = d.groupby(["provisional_measurement_tier", "L1", "L2"]).size().to_dict()
    assert got == JOINT_EXPECT
    assert len(got) == 5, "only 5 of 18 joint cells may be occupied"


@requires_stage12
def test_d2_level1_and_level2_evaluability_are_perfectly_coupled():
    d = _matrix()
    assert int(((d.L1 == "NOT_EVALUABLE") != (d.L2 == "NOT_EVALUABLE")).sum()) == 0


@requires_stage12
def test_d3_measurement_high_and_structure_supported_are_disjoint():
    d = _matrix()
    rh = set(d[d.provisional_measurement_tier == "robust-high"].patient_id)
    l1 = set(d[d.L1 == "SUPPORTED"].patient_id)
    assert rh & l1 == set()


@requires_stage12
def test_d4_no_patient_has_convergent_evidence():
    d = _matrix()
    conv = d[(d.provisional_measurement_tier == "robust-high") &
             (d.L1 == "SUPPORTED") & (d.L2 == "SUPPORTED")]
    assert len(conv) == 0


@requires_stage12
def test_d5_concordance_file_reproduces_the_crosstabs():
    cm = pd.read_csv(OUT / "stage12_concordance_matrix.csv")
    assert set(cm.table) == {"tier_x_level1", "tier_x_level2", "level1_x_level2",
                             "joint_tier_level1_level2"}
    assert int(cm[cm.table == "tier_x_level1"].n.sum()) == 32


# ============================================= E. no classifier, score or ranking
BANNED_CALLS = ("PCA", "TSNE", "UMAP", "KMeans", "AgglomerativeClustering",
                "LogisticRegression", "RandomForest", "fit_transform", "sklearn",
                "argsort", "rank(")


def test_e1_producer_contains_no_classifier_or_reduction():
    src = PRODUCER.read_text()
    for bad in BANNED_CALLS:
        assert bad not in src, f"Stage-12 producer references {bad}"


def test_e2_no_weighted_aggregation_across_evidence_columns():
    """AST scan with docstrings stripped: no dot products, no axis-1 sums over evidence."""
    for path in (PRODUCER, MODULE):
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.BinOp) and isinstance(node.op, ast.MatMult):
                raise AssertionError(f"{path.name}: matrix-multiply aggregation")
            if isinstance(node, ast.Call):
                fn = getattr(node.func, "attr", "")
                if fn in ("dot", "average", "wsum"):
                    raise AssertionError(f"{path.name}: {fn}() aggregation")
                if fn == "sum":
                    for kw in node.keywords:
                        if kw.arg == "axis" and getattr(kw.value, "value", None) == 1:
                            raise AssertionError(f"{path.name}: row-wise sum over columns")


def test_e3_no_score_or_rank_identifiers():
    banned = re.compile(r"\b(risk_score|utility_score|composite_score|evidence_score|"
                        r"weighted_score|patient_rank|risk_rank|points_total)\b")
    for path in (PRODUCER, MODULE):
        assert not banned.search(path.read_text()), path.name


@requires_stage12
def test_e4_no_output_column_is_a_score_or_rank():
    d = _matrix()
    for c in d.columns:
        low = c.lower()
        assert "score" not in low or low == "n_uncertainty_flags", c
        assert "rank" not in low and "risk_class" not in low, c


@requires_stage12
def test_e5_n_uncertainty_flags_is_a_count_not_a_score():
    d = _matrix()
    flags = ["low_n", "denominator_unstable", "depth_sensitive", "repeated_sample_unstable",
             "null_scheme_sensitive", "dropout_compatible", "intermediate_dn"]
    recomputed = d[flags].astype(bool).sum(axis=1)
    assert (d.n_uncertainty_flags >= recomputed).all()
    assert d.n_uncertainty_flags.dtype.kind in "iu"


@requires_stage12
def test_e6_no_patient_synthesis_category_was_created():
    """The design concluded: evidence matrix only, no category."""
    d = _matrix()
    for c in d.columns:
        assert "category" not in c.lower() and "synthesis_class" not in c.lower(), c
    if "evidence_profile" in d.columns:
        # permitted only as a lossless concatenation of four frozen states
        for v in d.evidence_profile:
            assert v.startswith("tier=") and "|L1=" in v and "|L2=" in v and "|L3=" in v
        assert d.evidence_profile.nunique() == 5


# ============================================================ F. Level-2 method
def test_f1_stage10_method_wording_is_wilcoxon_not_pydeseq2():
    text = SY.STAGE10_DE_METHOD.lower()
    assert "wilcoxon signed-rank" in text and "benjamini" in text
    for bad in ("pydeseq2", "deseq2", "negative-binomial", "~ patient + group"):
        assert bad not in text


@requires_stage12
def test_f2_no_stage12_output_attributes_the_de_to_pydeseq2():
    for p in list(OUT.glob("*.md")) + list(OUT.glob("*.csv")):
        if p.name == "stage12_design_snapshot.md":
            continue          # the snapshot quotes the design, which discusses the mislabel
        low = p.read_text().lower()
        for bad in ("pydeseq2", "deseq2", "negative-binomial glm"):
            assert bad not in low, f"{p.name} attributes Stage-10 DE to {bad}"


@requires_stage12
def test_f3_gamma_secretase_is_reported_as_not_supported():
    cl = pd.read_csv(OUT / "stage12_claim_ladder.csv")
    txt = (OUT / "stage12_summary.md").read_text().lower()
    assert "not supported" in txt and "gamma-secretase" in txt.replace("γ", "gamma")
    assert "near-significant" not in txt


@requires_stage12
def test_f4_single_denominator_bh_for_gamma_secretase_is_not_quoted():
    """0.0694 is the sensitivity-only value; quoting it alone misrepresents the result."""
    txt = (OUT / "stage12_summary.md").read_text()
    for bad in ("0.0694", "0.069"):
        assert bad not in txt


@requires_stage12
def test_f5_only_the_three_both_denominator_programmes_are_called_reproducible():
    txt = (OUT / "stage12_summary.md").read_text().lower()
    assert "antigen presentation" in txt and "oxphos" in txt and "interferon" in txt
    assert "both-denominators rule" in txt or "both denominators" in txt


# ========================================================= G. genomic semantics
def test_g1_cnv_not_supported_is_never_emitted():
    for path in (PRODUCER, MODULE):
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and isinstance(node.value, str):
                assert "CNV_SUBCLONE_NOT_SUPPORTED" not in node.value, path.name
    assert not hasattr(SC, "CNV_SUBCLONE_NOT_SUPPORTED")


@requires_stage12
def test_g2_not_evaluable_is_never_coerced_to_negative():
    d = _matrix()
    ne = d[d.L1 == "NOT_EVALUABLE"]
    assert len(ne) == 5
    assert (ne.level1_state == "DN_STRUCTURE_NOT_EVALUABLE").all()
    assert (ne.level2_state == "DN_STATE_NOT_EVALUABLE").all()
    assert (ne.main_uncertainty == "not_evaluable").all()
    for txt in ne.biological_interpretation:
        assert "absence of evidence, not negative evidence" in txt


# ========================================================== H. immune semantics
@requires_stage12
def test_h1_no_patient_is_marked_immune_supported():
    d = _matrix()
    assert (~d.composition_supported.astype(bool)).all()
    assert (~d.communication_supported.astype(bool)).all()
    assert (~d.liana_independent_support.astype(bool)).all()


@requires_stage12
def test_h2_liana_circular_hit_is_not_independent_support():
    txt = (OUT / "stage12_summary.md").read_text().lower()
    assert "circular" in txt
    cs = pd.read_csv(OUT / "stage12_cohort_summary.csv")
    assert cs[cs.domain == "immune"].strength_of_evidence.iloc[0] == "NOT_SUPPORTED"


# ======================================================== I. coverage semantics
@requires_stage12
def test_i1_gprc5d_and_sdc1_remain_not_evaluable():
    d = _matrix()
    for v in d.not_evaluable_targets:
        assert "GPRC5D" in v and "SDC1" in v


@requires_stage12
def test_i2_no_therapeutic_ranking_terms_in_coverage_columns():
    d = _matrix()
    cov_cols = [c for c in d.columns if "cover" in c.lower() or "pair" in c.lower()
                or "triple" in c.lower()]
    for c in cov_cols:
        assert "best" not in c.lower() and "optimal" not in c.lower()
        assert "recommend" not in c.lower()
    cs = pd.read_csv(OUT / "stage12_cohort_summary.csv")
    assert cs[cs.domain == "coverage"].strength_of_evidence.iloc[0] == "EXPLORATORY"


@requires_stage12
def test_i3_coverage_is_flagged_depth_sensitive_for_every_patient():
    assert _matrix().coverage_depth_sensitive.astype(bool).all()


# ============================================================= J. claim guards
#: Phrases Stage 12 may never assert. They are permitted ONLY in the columns whose
#: purpose is to record what must not be said.
PROHIBITED = ("confirmed escape subclone", "confirmed subclone",
              "pre-existing resistant clone", "high-risk patient", "low-risk patient",
              "immune escape mechanism", "immune-evasion axis", "optimal target pair",
              "best therapeutic", "gprc5d is redundant", "clinically safe",
              "normal tissue safe", "cnv-negative", "no subclone exists",
              "level-2 predicts patient risk", "bulk validates double-negative")
PROHIBITION_COLUMNS = {"prohibited_wording", "prohibited_claim", "prohibited_interpretation"}


@requires_stage12
def test_j1_prohibited_claims_appear_only_in_prohibition_columns():
    for p in list(OUT.glob("*.csv")) + list(OUT.glob("*.tsv")):
        d = pd.read_csv(p, sep="\t" if p.suffix == ".tsv" else ",", dtype=str).fillna("")
        for c in d.columns:
            if c in PROHIBITION_COLUMNS:
                continue
            joined = " ".join(d[c].astype(str)).lower()
            for bad in PROHIBITED:
                assert bad not in joined, f"{p.name}:{c} asserts '{bad}'"


@requires_stage12
def test_j2_summary_uses_prohibited_phrases_only_when_rejecting_them():
    txt = (OUT / "stage12_summary.md").read_text().lower()
    ctx = re.compile(r"prohibit|must not|never|rejected|not an? |unavailable|cannot say|"
                     r"is not\b|no combination|explicitly rejected|\| `not_")
    for bad in PROHIBITED:
        for m in re.finditer(re.escape(bad), txt):
            window = txt[max(0, m.start() - 220):m.end() + 120]
            assert ctx.search(window), f"summary asserts '{bad}'"


@requires_stage12
def test_j3_allowed_claim_never_exceeds_the_licensed_language():
    d = _matrix()
    for _, r in d.iterrows():
        lic = SC.licensed_language(r.level1_state, r.level2_state, SC.CNV_NOT_EVALUABLE)
        assert r.allowed_claim == SY.allowed_claim(lic)
        assert "subclone" not in r.allowed_claim.lower()


@requires_stage12
def test_j4_every_row_prohibits_the_two_universal_claims():
    for v in _matrix().prohibited_claim:
        assert "no genomic subclone claim" in v
        assert "no clinical recommendation" in v


# ============================================================== K. determinism
def test_k1_synthesis_text_functions_are_pure_and_deterministic():
    flags = {"uncertain_low_n": True, "uncertain_dropout_compatible": True}
    a = SY.main_uncertainty("DN_STRUCTURE_NOT_SUPPORTED", "DN_STATE_SUPPORTED", flags)
    b = SY.main_uncertainty("DN_STRUCTURE_NOT_SUPPORTED", "DN_STATE_SUPPORTED", flags)
    assert a == b == "low_n"
    assert SY.measurement_interpretation("robust-high", 0, False) == \
        SY.measurement_interpretation("robust-high", 0, False)


def test_k2_not_evaluable_outranks_every_measurement_flag():
    flags = {c: True for c in SY.FLAG_COLUMNS}
    assert SY.main_uncertainty("DN_STRUCTURE_NOT_EVALUABLE", "DN_STATE_NOT_EVALUABLE",
                               flags) == "not_evaluable"


def test_k3_uncertainty_priority_is_the_frozen_order():
    assert SY.UNCERTAINTY_PRIORITY == (
        "not_evaluable", "low_n", "dropout_compatible", "denominator",
        "repeated_sample", "depth", "null_scheme", "intermediate_dn", "none")


def test_k4_direction_sign_treats_missing_as_not_evaluable():
    assert SY.direction_sign(float("nan")) == "not_evaluable"
    assert SY.direction_sign(None) == "not_evaluable"
    assert SY.direction_sign(-1.0) == "lower_in_DN"
    assert SY.direction_sign(1.0) == "higher_in_DN"


# ================================================= L. no new statistics were run
def test_l1_producer_runs_no_new_statistical_test():
    src = PRODUCER.read_text()
    for bad in ("ttest", "mannwhitney", "wilcoxon(", "fisher_exact", "chi2_contingency",
                "pearsonr", "spearmanr", "linregress", "ols(", "multipletests",
                "benjamini_hochberg("):
        assert bad not in src.lower(), f"Stage-12 producer runs a new test: {bad}"


@requires_stage12
def test_l2_no_new_pvalue_column_was_created():
    d = _matrix()
    pcols = [c for c in d.columns if c.lower().startswith("p_") or c.lower().endswith("_p")
             or "pval" in c.lower()]
    #: only frozen Stage-08/10 p-values may be carried through
    allowed = {"perm_p_depth_stratified", "depth_stratified_p_primary",
               "depth_stratified_p_sensitivity", "unconditioned_p_primary"}
    assert set(pcols) <= allowed, f"new p-value columns: {set(pcols) - allowed}"


# ================================================================= M. cohort docs
@requires_stage12
def test_m1_cohort_summary_has_seven_domains_and_frozen_vocabulary():
    cs = pd.read_csv(OUT / "stage12_cohort_summary.csv")
    assert list(cs.domain) == ["measurement", "co-loss", "structure", "phenotype",
                               "genomic", "immune", "coverage"]
    assert cs.strength_of_evidence.isin(SY.STRENGTH_VOCABULARY).all()


@requires_stage12
def test_m2_claim_ladder_covers_the_twelve_claims_with_frozen_vocabulary():
    cl = pd.read_csv(OUT / "stage12_claim_ladder.csv")
    assert len(cl) == 12
    assert cl.current_status.isin(SY.STRENGTH_VOCABULARY).all()
    for col in ("claim", "evidence_required", "frozen_evidence", "current_status",
                "allowed_wording", "prohibited_wording"):
        assert col in cl.columns


@requires_stage12
def test_m3_uncertainty_register_keeps_three_separate_classes():
    ur = pd.read_csv(OUT / "stage12_uncertainty_register.csv")
    assert set(ur.uncertainty_class) == {"measurement", "biological", "external_validity"}
    for c in ur.columns:
        assert "score" not in c.lower() and "total" not in c.lower()
