"""Stage-11b LIANA verification arm — guards, not a re-analysis.

LIANA is an *additional exploratory arm*. The frozen Stage-11 custom result is not
reopened, not replaced and not edited, and LIANA may not create a patient state. Most of
what follows checks exactly that, because the failure mode here is not a wrong number — it
is a well-dressed screen quietly acquiring authority it was never given.
"""
from __future__ import annotations

import hashlib
import inspect
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mm_escape import communication as CM

REPO = Path(__file__).resolve().parents[1]
STAGE11 = REPO / "results" / "11_immune_context"
LIANA = STAGE11 / "liana_verification"
NOTEBOOK = REPO / "notebooks" / "11b_liana_verification.py"

requires_liana_run = pytest.mark.skipif(
    not (LIANA / "liana_method_config.json").exists(),
    reason="needs the LIANA arm outputs (run notebooks/11b_liana_verification.py)")


def config():
    return json.loads((LIANA / "liana_method_config.json").read_text())


# ------------------------------------------------------------------ A / B. receivers
@requires_liana_run
def test_a_primary_receiver_is_the_frozen_clone_primary_population():
    cfg = config()
    assert cfg["primary_receiver"].startswith("clone_primary")
    scores = pd.read_csv(LIANA / "liana_patient_level_scores.csv.gz", dtype={"patient": str})
    assert "consensus_clone_primary" in set(scores.run)


@requires_liana_run
def test_b_all_plasma_is_a_separate_sensitivity_receiver_not_a_replacement():
    """The two receivers must both exist and stay distinguishable — the frozen Stage-11
    amendment turned on exactly this distinction."""
    cfg = config()
    assert cfg["sensitivity_receiver"].startswith("all_plasma")
    scores = pd.read_csv(LIANA / "liana_patient_level_scores.csv.gz", dtype={"patient": str})
    runs = set(scores.run)
    assert {"consensus_clone_primary", "consensus_all_plasma"} <= runs
    a = scores[scores.run == "consensus_clone_primary"]
    b = scores[scores.run == "consensus_all_plasma"]
    assert len(a) and len(b)
    merged = a.merge(b, on=["patient", "source", "interaction"], suffixes=("_p", "_s"))
    assert len(merged), "the two receivers must be comparable on shared interactions"
    assert not np.allclose(merged.score_p, merged.score_s), (
        "clone-primary and all-plasma receivers must not be the same computation")


# --------------------------------------------------------- C. no pseudoreplication
@requires_liana_run
def test_c_dn_association_uses_patient_level_scores_not_cells():
    """n can never exceed the patient count — cells are not replicates."""
    assoc = pd.read_csv(LIANA / "liana_vs_dn_associations.csv")
    n_patients = pd.read_csv(
        REPO / "results/08_dual_antigen_escape/patient_antigen_states_primary.csv",
        dtype={"patient": str}).patient.nunique()
    assert assoc.n.max() <= n_patients, "an association saw more units than there are patients"
    assert assoc.n.min() >= config()["min_patients_for_association"]


@requires_liana_run
def test_c_per_patient_scores_come_from_a_per_patient_computation():
    """A pooled run's global scores must not be presented as patient observations."""
    assert "by_sample" in config()["per_patient"]
    src = NOTEBOOK.read_text()
    assert "by_sample(" in src and "sample_key=\"patient_id\"" in src
    scores = pd.read_csv(LIANA / "liana_patient_level_scores.csv.gz", dtype={"patient": str})
    per_pair = scores.groupby(["run", "source", "interaction"]).score.nunique()
    assert per_pair.max() > 1, "scores are identical across patients — this is a pooled run"


# ------------------------------------------------------- D. the frozen panel is frozen
def test_d_frozen_17_pair_candidate_panel_is_unchanged():
    pairs = [(l, r) for l, r in CM.LR_CANDIDATES if r != "None"]
    assert len(pairs) == 17
    assert ("PDCD1", "CD274") in pairs and ("TNFSF10", "TNFRSF10A") in pairs
    assert len(CM.LR_CANDIDATES) == 19          # PRF1/GZMB carry no receptor
    assert CM.LR_SENDERS == ("Tcell", "NK_core", "cytotoxic_mixed", "Myeloid")


@requires_liana_run
def test_d_crosscheck_covers_every_frozen_pair_and_sender():
    cc = pd.read_csv(LIANA / "liana_candidate_panel_crosscheck.csv")
    pairs = {(l, r) for l, r in CM.LR_CANDIDATES if r != "None"}
    assert {(l, r) for l, r in zip(cc.ligand, cc.receptor)} == pairs
    assert len(cc) == len(pairs) * len(CM.LR_SENDERS)
    for col in ("in_liana_resource", "evaluated_by_liana", "panel_status"):
        assert col in cc.columns


# --------------------------------------------- E. receiver-only falsification exists
@requires_liana_run
def test_e_receiver_side_decomposition_exists_for_every_highlighted_interaction():
    dec = pd.read_csv(LIANA / "liana_receiver_side_confound.csv")
    assoc = pd.read_csv(LIANA / "liana_vs_dn_associations.csv")
    primary = assoc[assoc.run == "consensus_clone_primary"]
    highlighted = {(r.sender, r.interaction) for _, r in primary[primary.p_adjusted < 0.05].iterrows()}
    covered = {(r.sender, r.interaction) for _, r in dec.iterrows()}
    assert highlighted <= covered, "every raw-significant interaction must be decomposed"
    assert {"receptor_coef", "receptor_p", "ligand_coef", "ligand_p"} <= set(dec.columns)

    # The named interactions must be *accounted for*: either decomposed, or explicitly
    # recorded as unevaluable. Silence about them is what this forbids.
    cc = pd.read_csv(LIANA / "liana_candidate_panel_crosscheck.csv")
    for lig, rec in (("PDCD1", "CD274"), ("TNFSF10", "TNFRSF10A"), ("TNFSF10", "TNFRSF10B")):
        name = f"{lig}->{rec}"
        if name in set(dec.interaction):
            continue
        rows = cc[(cc.ligand == lig) & (cc.receptor == rec)]
        assert len(rows), f"{name} is neither decomposed nor present in the cross-check"
        floor = config()["min_patients_for_association"]
        # Legitimate reasons not to decompose: LIANA's resource does not contain it, or it
        # never reached the patient floor so no association exists to decompose. Being
        # testable and silently skipped is what this forbids.
        untestable = (~rows.in_liana_resource) | (rows.n_patients_evaluated < floor)
        assert untestable.all(), (
            f"{name} reached the {floor}-patient floor but was never decomposed:\n"
            f"{rows[['sender', 'n_patients_evaluated', 'panel_status']].to_string(index=False)}")


@requires_liana_run
def test_e_receiver_confounded_interactions_are_labelled_not_called_signalling():
    dec = pd.read_csv(LIANA / "liana_receiver_side_confound.csv")
    allowed = {"CONSISTENT_WITH_TARGETED_PANEL", "NOT_REPRODUCED_BY_LIANA",
               "RECEIVER_STATE_CONFOUNDED", "ABUNDANCE_SENSITIVE",
               "EXPLORATORY_LIANA_ONLY", "NOT_EVALUABLE"}
    assert set(dec.classification) <= allowed


# ------------------------------------------------------------- F. the frozen model
@requires_liana_run
def test_f_confound_model_is_the_frozen_stage11_one():
    assert tuple(config()["confound_model"]) == CM.CONFOUNDERS
    assoc = pd.read_csv(LIANA / "liana_vs_dn_associations.csv")
    for col in ("coef_unadjusted", "p_unadjusted", "coef_adjusted", "p_adjusted",
                "ci_lo", "ci_hi", "p_adj_BH"):
        assert col in assoc.columns, col
    assert {"rho_MMRF", "rho_WU1", "rho_WU2"} <= set(assoc.columns)


def test_f_no_liana_specific_model_was_invented():
    """The notebook must call the shared stage-11 estimator, not a local one."""
    src = NOTEBOOK.read_text()
    assert "CM.ols_association" in src
    assert "CM.benjamini_hochberg" in src
    for banned in ("statsmodels", "mixedlm", "GLM(", "logit("):
        assert banned not in src, f"a separate modelling stack appeared: {banned}"


# ------------------------------------------------- G. the full space is retained
@requires_liana_run
def test_g_full_interaction_space_is_written_not_only_significant_rows():
    assoc = pd.read_csv(LIANA / "liana_vs_dn_associations.csv")
    assert (assoc.p_adjusted >= 0.05).sum() > 0, "non-significant rows must be retained"
    assert len(assoc) > 100, "the tested space must be reported whole"
    assert (LIANA / "liana_full_results.csv.gz").exists()


@requires_liana_run
def test_g_only_one_resource_was_used_so_databases_were_not_shopped():
    cfg = config()
    assert cfg["resource"] == "consensus"
    src = NOTEBOOK.read_text()
    for other in ("cellphonedb'", "cellchatdb'", "italk", "ramilowski2015", "connectomedb2020"):
        assert f"select_resource('{other}" not in src and f'select_resource("{other}' not in src


@requires_liana_run
def test_g_both_methods_are_preserved_not_only_the_stronger():
    scores = pd.read_csv(LIANA / "liana_patient_level_scores.csv.gz", dtype={"patient": str})
    assert {"consensus_clone_primary", "cellchat_clone_primary"} <= set(scores.run)


# ------------------------------------------ H / I / J. LIANA changes nothing
FROZEN_STATE = {
    "stage09b_tiers": "results/08_dual_antigen_escape/risk_tier_provisional/risk_tiers_provisional.csv",
    "stage10_states": "results/10_dn_coherence/dn_coherence_final_states.csv",
    "stage11_custom_lr": "results/11_immune_context/communication_context_vs_dn.csv",
    "stage11_composition": "results/11_immune_context/patient_immune_composition.csv",
    "coverage_qc": "results/08_dual_antigen_escape/multi_antigen_coverage/target_measurement_qc.csv",
}


@requires_liana_run
def test_hi_liana_did_not_modify_any_frozen_tier_state_or_custom_lr_output():
    recorded = json.loads((LIANA / "frozen_state_digests.json").read_text())
    for key, rel in FROZEN_STATE.items():
        p = REPO / rel
        if not p.exists() or key not in recorded:
            pytest.skip(f"{rel} absent")
        assert hashlib.sha256(p.read_bytes()).hexdigest() == recorded[key], (
            f"LIANA arm changed a frozen artifact: {rel}")


@requires_liana_run
def test_h_liana_outputs_live_only_in_their_own_namespace():
    produced = {p.name for p in LIANA.glob("*")}
    assert produced, "the LIANA arm produced nothing"
    for name in produced:
        assert not (STAGE11 / name).exists() or (STAGE11 / name).is_dir(), (
            f"{name} leaked into the frozen stage-11 namespace")


@requires_liana_run
def test_j_no_liana_patient_classifier_was_created():
    """Interaction labels are allowed; a per-patient state is not."""
    banned_cols = {"liana_tier", "liana_state", "liana_class", "liana_risk",
                   "immune_evasion", "patient_classification"}
    for f in LIANA.glob("*.csv"):
        cols = {c.lower() for c in pd.read_csv(f, nrows=0).columns}
        assert not (cols & banned_cols), (f.name, cols & banned_cols)
    scores = pd.read_csv(LIANA / "liana_patient_level_scores.csv.gz", nrows=5)
    assert not ({"tier", "state", "classification"} & {c.lower() for c in scores.columns})


def test_j_forbidden_label_is_never_emitted_as_a_value():
    """`IMMUNE_EVASION_CONFIRMED` may be *named as forbidden* in prose, but must never be
    emitted as a label. The distinction matters: the notebook and summary say the label is
    unavailable, and a naive substring scan would flag that prohibition as a violation."""
    if not LIANA.exists():
        pytest.skip("arm not run")
    for f in LIANA.glob("*.csv"):
        df = pd.read_csv(f)
        for col in df.columns:
            if df[col].dtype == object:
                assert "IMMUNE_EVASION_CONFIRMED" not in set(df[col].astype(str)), (f.name, col)
    allowed = {"CONSISTENT_WITH_TARGETED_PANEL", "NOT_REPRODUCED_BY_LIANA",
               "RECEIVER_STATE_CONFOUNDED", "ABUNDANCE_SENSITIVE",
               "EXPLORATORY_LIANA_ONLY", "NOT_EVALUABLE"}
    dec = pd.read_csv(LIANA / "liana_receiver_side_confound.csv")
    assert set(dec.classification) <= allowed


@requires_liana_run
def test_j_arm_is_declared_exploratory_and_non_tier_changing():
    summary = LIANA / "liana_summary.md"
    if not summary.exists():
        pytest.skip("summary not yet written")
    text = summary.read_text().lower()
    assert "exploratory" in text
    assert "non-tier-changing" in text or "no tier" in text


# ============================================================================
# Frozen 2026-08-26. These guard the accepted LIANA interpretation itself, not
# just the mechanics that produced it.
# ============================================================================

#: The two antigens whose negativity defines `obs_dn_primary`. An interaction touching
#: either is structurally circular for the DN-vs-communication question.
DN_DEFINING = {"TNFRSF17", "GPRC5D"}


def _is_antigen_circular(interaction: str) -> bool:
    lig, rec = interaction.split("->")
    return bool((set(lig.split("_")) | set(rec.split("_"))) & DN_DEFINING)


@requires_liana_run
def test_antigen_circularity_flag_is_present_and_correct():
    """Permanent. The screen's top hit is circular by construction, and the flag that says
    so must not quietly disappear from a future rerun."""
    a = pd.read_csv(LIANA / "liana_vs_dn_associations.csv")
    assert "antigen_circular" in a.columns, "the antigen-circularity flag is mandatory"
    recomputed = a.interaction.map(_is_antigen_circular)
    assert (a.antigen_circular.astype(bool) == recomputed).all()


@requires_liana_run
def test_the_circular_row_is_preserved_and_flagged_never_excluded():
    """Deleting the row would hide what the screen ranked first. It is flagged, not dropped."""
    a = pd.read_csv(LIANA / "liana_vs_dn_associations.csv")
    hit = a[(a.interaction == "TNFSF13B->TNFRSF17") & (a.sender == "Myeloid")]
    assert len(hit) >= 1, "the top consensus hit must remain in the output"
    assert hit.antigen_circular.all()
    dec = pd.read_csv(LIANA / "liana_receiver_side_confound.csv")
    row = dec[dec.interaction == "TNFSF13B->TNFRSF17"]
    assert len(row) and (row.classification == "RECEIVER_STATE_CONFOUNDED").all()


@requires_liana_run
def test_every_bh_significant_consensus_hit_is_accounted_for():
    """A BH hit is only allowed to stand if it is circular, receiver-confounded or
    abundance-sensitive. An unexplained survivor would contradict the frozen conclusion and
    must fail loudly rather than be discovered later in prose."""
    a = pd.read_csv(LIANA / "liana_vs_dn_associations.csv")
    dec = pd.read_csv(LIANA / "liana_receiver_side_confound.csv")
    ab = pd.read_csv(LIANA / "liana_abundance_sensitivity.csv")
    confounded = {(r.sender, r.interaction) for _, r in
                  dec[dec.classification == "RECEIVER_STATE_CONFOUNDED"].iterrows()}
    abundance = {(r.sender, r.interaction) for _, r in ab[ab.abundance_sensitive].iterrows()}
    hits = a[(a.run == "consensus_clone_primary") & (a.p_adj_BH < 0.10)]
    for _, r in hits.iterrows():
        key = (r.sender, r.interaction)
        assert r.antigen_circular or key in confounded or key in abundance, (
            f"unexplained BH-significant LIANA hit: {key}")


@requires_liana_run
def test_25183_is_a_declared_evaluability_exclusion_not_a_software_failure():
    """It contributes zero immune sender cells, so no sender->receiver pair exists. That is a
    biological/evaluability exclusion and is consistent with the frozen Stage-11 design,
    where 25183 is likewise the one patient absent from the custom communication table."""
    ev = pd.read_csv(LIANA / "liana_patient_evaluability.csv", dtype={"patient": str})
    rows = ev[ev.patient == "25183"]
    assert len(rows) == 3, "25183 must be recorded for every run, not omitted"
    assert (rows.status == "LIANA_NOT_EVALUABLE").all()
    assert (rows.n_sender_categories_ge_min == 0).all()
    assert (ev[ev.status == "EVALUABLE"].groupby("run").size() == 31).all()

    frozen_custom = pd.read_csv(STAGE11 / "communication_context.csv", dtype={"patient": str})
    assert "25183" not in set(frozen_custom.patient), (
        "the frozen custom arm also excludes 25183 — the two must agree")


@requires_liana_run
def test_pdcd1_cd274_is_absent_from_the_resource_not_disproved():
    """`NOT_EVALUABLE_BY_LIANA_RESOURCE`. LIANA neither supports nor refutes it."""
    cc = pd.read_csv(LIANA / "liana_candidate_panel_crosscheck.csv")
    rows = cc[(cc.ligand == "PDCD1") & (cc.receptor == "CD274")]
    assert len(rows) == len(CM.LR_SENDERS)
    assert not rows.in_liana_resource.any()
    assert not rows.evaluated_by_liana.any()
    assert (rows.panel_status == "NOT_EVALUABLE").all()
    assert rows.coef_adjusted.isna().all(), "an untested interaction cannot carry an estimate"


@requires_liana_run
def test_no_interaction_is_labelled_consistent_with_the_targeted_panel():
    """Frozen result: LIANA is a partial verification arm, not a reproduction."""
    cc = pd.read_csv(LIANA / "liana_candidate_panel_crosscheck.csv")
    assert (cc.panel_status == "CONSISTENT_WITH_TARGETED_PANEL").sum() == 0
    in_resource = cc.groupby(["ligand", "receptor"]).in_liana_resource.first()
    assert int(in_resource.sum()) == 8, "8 of 17 frozen pairs are in LIANA's resource"


@requires_liana_run
def test_method_disagreement_stays_visible():
    """CellChat produces more BH hits than the consensus. That disagreement is a reason for
    caution and must remain inspectable, not be collapsed to the more exciting arm."""
    a = pd.read_csv(LIANA / "liana_vs_dn_associations.csv")
    runs = set(a.run)
    assert {"consensus_clone_primary", "cellchat_clone_primary", "consensus_all_plasma"} <= runs
    per_run = a[a.p_adj_BH < 0.10].groupby("run").size()
    assert per_run.get("cellchat_clone_primary", 0) > per_run.get("consensus_clone_primary", 0)
