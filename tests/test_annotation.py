"""Stage 06 annotation logic. Data-free: everything here runs on synthetic AnnData."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from mm_escape import annotation as ann
from mm_escape import config

AMB = config.AMBIGUOUS_LABEL


# --------------------------------------------------------------------------
# label collapse
# --------------------------------------------------------------------------

def test_celltypist_map_targets_are_project_classes():
    bad = set(ann.CELLTYPIST_TO_BROAD.values()) - set(ann.BROAD_CLASSES) - {AMB}
    assert not bad, f"CellTypist map emits non-project classes: {bad}"


def test_singler_map_targets_are_project_classes():
    bad = set(ann.SINGLER_TO_BROAD.values()) - set(ann.BROAD_CLASSES) - {AMB}
    assert not bad, f"SingleR map emits non-project classes: {bad}"


def test_ilc_maps_to_nk():
    """Checked against markers at stage 05b, not assumed: ILC here is NK."""
    assert ann.CELLTYPIST_TO_BROAD["ILC"] == "NK"


def test_novershtern_has_no_plasma_route():
    """The reference has no plasma-cell class, so the map must not invent one.

    Novershtern's B-lineage stops at 'Mature B cells class switched'. Its plasma cells
    arrive as 'B cells'. If someone later adds a Novershtern->PlasmaCell key, that is a
    silent substitution and this test is what catches it.
    """
    assert ann.SINGLER_TO_BROAD["B cells"] == "Bcell"
    novershtern_keys = {"B cells", "CD4+ T cells", "CD8+ T cells", "NK cells", "HSCs",
                        "CMPs", "GMPs", "MEPs", "Erythroid cells", "Monocytes",
                        "Dendritic cells", "Granulocytes", "Basophils", "Eosinophils",
                        "Megakaryocytes", "NK T cells"}
    assert not any(ann.SINGLER_TO_BROAD[k] == "PlasmaCell" for k in novershtern_keys)


def test_hpca_fine_is_the_only_singler_plasma_route():
    assert ann.SINGLER_TO_BROAD["B_cell:Plasma_cell"] == "PlasmaCell"


def test_collapse_unmapped_becomes_ambiguous_not_guessed():
    out = ann.collapse_labels(["T cells", "Wat", None], ann.CELLTYPIST_TO_BROAD)
    assert list(out) == ["Tcell", AMB, AMB]


def test_collapse_fine_label_falls_back_to_main_prefix():
    """An unlisted HPCA fine variant still lands on its main class, not Ambiguous."""
    out = ann.collapse_labels(["B_cell:some_new_subtype"], ann.SINGLER_TO_BROAD)
    assert list(out) == ["Bcell"]


def test_collapse_fine_exact_key_beats_prefix():
    """B_cell:Plasma_cell must NOT be collapsed to Bcell by the prefix rule."""
    out = ann.collapse_labels(["B_cell:Plasma_cell"], ann.SINGLER_TO_BROAD)
    assert list(out) == ["PlasmaCell"]


# --------------------------------------------------------------------------
# fixtures
# --------------------------------------------------------------------------

def _adata(labels, expr: dict[str, list[float]], *, pad: int = 0):
    """Synthetic AnnData. `pad` adds filler genes, which `scanpy.tl.score_genes`
    needs as a control pool — with only the panel genes present it raises
    "No control genes found in any cut"."""
    import anndata as ad

    expr = dict(expr)
    n = len(next(iter(expr.values())))
    rng = np.random.default_rng(0)
    for i in range(pad):
        expr[f"FILLER{i}"] = list(rng.uniform(0.0, 4.0, size=n))
    genes = list(expr)
    X = np.array([expr[g] for g in genes], dtype=float).T
    a = ad.AnnData(X)
    a.var_names = genes
    a.obs_names = [f"c{i}" for i in range(a.n_obs)]
    a.obs["lab"] = pd.Categorical(labels)
    return a


# --------------------------------------------------------------------------
# marker scoring / cluster assignment
# --------------------------------------------------------------------------

def test_score_marker_panel_warns_and_skips_absent_genes():
    a = _adata(["x"] * 6, {"CD3D": [5, 5, 5, 0, 0, 0], "MS4A1": [0, 0, 0, 5, 5, 5]}, pad=60)
    with pytest.warns(UserWarning):
        written = ann.score_marker_panel(a, {"Tcell": ("CD3D", "NOPE"), "HSPC": ("ABSENT",)})
    assert written == ["score_Tcell"]
    assert "score_HSPC" not in a.obs


def test_manual_labels_are_per_cluster_not_per_cell():
    a = _adata(["x"] * 4, {"g": [1.0, 1, 1, 1]})
    a.obs["clust"] = pd.Categorical(["0", "0", "1", "1"])
    a.obs["score_Tcell"] = [9.0, -9.0, 0.0, 0.0]   # cluster 0 mean 0, cell-level split
    a.obs["score_NK"] = [0.0, 0.0, 5.0, 5.0]
    tbl = ann.manual_labels_from_clusters(a, "clust")
    assert tbl.loc["1", "winner"] == "NK"
    # both cells of cluster 0 get the SAME label despite opposite per-cell scores
    assert a.obs["cell_type_manual"].iloc[0] == a.obs["cell_type_manual"].iloc[1]


def test_manual_margin_records_ambiguous_rather_than_forcing():
    a = _adata(["x"] * 2, {"g": [1.0, 1]})
    a.obs["clust"] = pd.Categorical(["0", "0"])
    a.obs["score_Tcell"] = [1.0, 1.0]
    a.obs["score_NK"] = [0.99, 0.99]
    tbl = ann.manual_labels_from_clusters(a, "clust", margin=0.5)
    assert tbl.loc["0", "winner"] == AMB
    assert set(a.obs["cell_type_manual"]) == {AMB}


def test_manual_requires_scores_first():
    a = _adata(["x"], {"g": [1.0]})
    a.obs["clust"] = pd.Categorical(["0"])
    with pytest.raises(ValueError, match="score_marker_panel"):
        ann.manual_labels_from_clusters(a, "clust")


# --------------------------------------------------------------------------
# marker coverage
# --------------------------------------------------------------------------

def test_marker_coverage_high_when_class_expresses_its_own_markers():
    labels = ["Tcell"] * 3 + ["NK"] * 3
    a = _adata(labels, {"CD3D": [9, 9, 9, 0, 0, 0], "NKG7": [0, 0, 0, 9, 9, 9]})
    cov = ann.marker_coverage(a, "lab", {"Tcell": ("CD3D",), "NK": ("NKG7",)})
    assert cov.loc["Tcell", "coverage"] == pytest.approx(1.0)
    assert cov.loc["NK", "coverage"] == pytest.approx(1.0)


def test_marker_coverage_low_when_label_lacks_marker_support():
    """The failure mode the veto exists for: cells labelled Tcell that are not T cells."""
    labels = ["Tcell"] * 3 + ["NK"] * 3
    a = _adata(labels, {"CD3D": [0, 0, 0, 9, 9, 9], "NKG7": [0, 0, 0, 9, 9, 9]})
    cov = ann.marker_coverage(a, "lab", {"Tcell": ("CD3D",), "NK": ("NKG7",)})
    assert cov.loc["Tcell", "coverage"] == pytest.approx(0.0)
    assert cov.loc["Tcell", "coverage"] < config.MARKER_COVERAGE_MIN


def test_marker_coverage_absent_class_is_nan_not_zero():
    """Not evaluable and failed are different outcomes and must not be conflated."""
    a = _adata(["Tcell"] * 3, {"CD3D": [9, 9, 9], "HBB": [0, 0, 0]})
    cov = ann.marker_coverage(a, "lab", {"Tcell": ("CD3D",), "Erythroid": ("HBB",)})
    assert np.isnan(cov.loc["Erythroid", "coverage"])
    assert cov.loc["Erythroid", "n_cells"] == 0


# --------------------------------------------------------------------------
# concordance
# --------------------------------------------------------------------------

def test_perfect_agreement_gives_f1_one():
    ref = ["Tcell", "NK", "PlasmaCell"]
    out = ann.per_class_concordance(ref, ref)
    for cls in ("Tcell", "NK", "PlasmaCell"):
        assert out.loc[cls, "f1"] == pytest.approx(1.0)


def test_total_disagreement_gives_f1_zero():
    out = ann.per_class_concordance(["Tcell", "Tcell"], ["NK", "NK"])
    assert out.loc["Tcell", "f1"] == 0.0


def test_concordance_length_mismatch_raises():
    with pytest.raises(ValueError, match="length mismatch"):
        ann.per_class_concordance(["Tcell"], ["Tcell", "NK"])


# --------------------------------------------------------------------------
# the decision rule
# --------------------------------------------------------------------------

def _tbl(d, col):
    return pd.DataFrame({col: d}).rename_axis("cell_class")


def test_highest_f1_wins_when_several_qualify():
    conc = {"celltypist": _tbl({"Tcell": 0.97}, "f1"), "singler": _tbl({"Tcell": 0.93}, "f1")}
    cov = {"celltypist": _tbl({"Tcell": 0.9}, "coverage"),
           "singler": _tbl({"Tcell": 0.9}, "coverage"),
           "manual": _tbl({"Tcell": 0.9}, "coverage")}
    d = ann.decide_per_class(conc, cov, classes=("Tcell",))
    assert d.loc["Tcell", "chosen_method"] == "celltypist"
    assert "singler" in d.loc["Tcell", "also_passed"]


def test_class_below_its_bar_falls_back_to_manual():
    conc = {"celltypist": _tbl({"PlasmaCell": 0.80}, "f1")}   # bar is 0.95
    cov = {"celltypist": _tbl({"PlasmaCell": 0.9}, "coverage"),
           "manual": _tbl({"PlasmaCell": 0.9}, "coverage")}
    d = ann.decide_per_class(conc, cov, classes=("PlasmaCell",))
    assert d.loc["PlasmaCell", "chosen_method"] == "manual"
    assert "fell back" in d.loc["PlasmaCell", "reason"]


def test_marker_coverage_vetoes_a_class_both_methods_agree_on():
    """THE load-bearing case: perfect automated agreement, no biological support.

    Both methods clear the concordance bar comfortably. Neither has marker support.
    The class must NOT go to either of them — agreement on a biologically unsupported
    label is agreement on an error.
    """
    conc = {"celltypist": _tbl({"Tcell": 0.99}, "f1"), "singler": _tbl({"Tcell": 0.98}, "f1")}
    cov = {"celltypist": _tbl({"Tcell": 0.05}, "coverage"),
           "singler": _tbl({"Tcell": 0.04}, "coverage"),
           "manual": _tbl({"Tcell": 0.85}, "coverage")}
    d = ann.decide_per_class(conc, cov, classes=("Tcell",))
    assert d.loc["Tcell", "chosen_method"] == "manual"
    assert set(d.loc["Tcell", "vetoed_by_coverage"].split(",")) == {"celltypist", "singler"}


def test_class_unresolved_when_even_the_fallback_fails_coverage():
    conc = {"celltypist": _tbl({"HSPC": 0.99}, "f1")}
    cov = {"celltypist": _tbl({"HSPC": 0.01}, "coverage"),
           "manual": _tbl({"HSPC": 0.02}, "coverage")}
    d = ann.decide_per_class(conc, cov, classes=("HSPC",))
    assert d.loc["HSPC", "chosen_method"] == AMB
    assert "coverage veto" in d.loc["HSPC", "reason"]


def test_thresholds_come_from_config_and_are_the_declared_ones():
    assert config.CONCORDANCE_THRESHOLDS["PlasmaCell"] == 0.95
    assert all(config.CONCORDANCE_THRESHOLDS[c] == 0.90 for c in ("Tcell", "NK", "Myeloid"))
    assert all(config.CONCORDANCE_THRESHOLDS[c] == 0.85 for c in ("Bcell", "Erythroid", "HSPC"))
    assert config.MARKER_COVERAGE_MIN == 0.30


# --------------------------------------------------------------------------
# final labels + state programs
# --------------------------------------------------------------------------

def test_unclaimed_cells_are_ambiguous_not_absorbed():
    a = _adata(["x"] * 3, {"g": [1.0, 1, 1]})
    a.obs["ct"] = ["Tcell", "Nonsense", "NK"]
    dec = pd.DataFrame({"chosen_method": {"Tcell": "celltypist", "NK": "celltypist"}}).rename_axis("cell_class")
    summary = ann.assemble_final_labels(a, dec, {"celltypist": "ct"})
    assert list(a.obs["cell_type"]) == ["Tcell", AMB, "NK"]
    assert summary.loc[AMB, "n_cells"] == 1


def test_final_labels_record_provenance():
    a = _adata(["x"] * 2, {"g": [1.0, 1]})
    a.obs["ct"] = ["Tcell", "Tcell"]
    a.obs["conf"] = [0.8, 0.4]
    dec = pd.DataFrame({"chosen_method": {"Tcell": "celltypist"}}).rename_axis("cell_class")
    ann.assemble_final_labels(a, dec, {"celltypist": "ct"}, conf_keys={"celltypist": "conf"})
    assert set(a.obs["annotation_source"]) == {"celltypist"}
    assert list(a.obs["annotation_conf"]) == [0.8, 0.4]


def test_state_programs_are_floats_and_never_enter_cell_type():
    a = _adata(["x"] * 4, {"MKI67": [5.0, 0, 5, 0], "TOP2A": [5.0, 0, 5, 0]}, pad=60)
    a.obs["cell_type"] = pd.Categorical(["PlasmaCell"] * 4)
    written = ann.score_state_programs(a, {"cell_cycle": ("MKI67", "TOP2A")})
    assert written == ["program_cell_cycle"]
    assert a.obs["program_cell_cycle"].dtype.kind == "f"
    assert set(a.obs["cell_type"]) == {"PlasmaCell"}


def test_programs_are_non_exclusive():
    """A cell may carry several programs at once; they are covariates, not categories."""
    a = _adata(["x"] * 4, {"MKI67": [5.0, 0, 5, 0], "ISG15": [5.0, 5, 0, 0]}, pad=60)
    ann.score_state_programs(a, {"cell_cycle": ("MKI67",), "interferon": ("ISG15",)})
    both = (a.obs["program_cell_cycle"] > 0) & (a.obs["program_interferon"] > 0)
    assert both.any()


# --------------------------------------------------------------------------
# majority voting (done in-house because CellTypist densifies to scale)
# --------------------------------------------------------------------------

def test_majority_vote_assigns_cluster_mode_to_every_cell():
    labels = ["Tcell", "Tcell", "NK", "Myeloid", "Myeloid", "Myeloid"]
    clusters = ["0", "0", "0", "1", "1", "1"]
    out = ann.majority_vote(labels, clusters)
    assert list(out) == ["Tcell"] * 3 + ["Myeloid"] * 3


def test_majority_vote_min_prop_marks_undecided_rather_than_taking_a_plurality():
    labels = ["Tcell", "NK", "Myeloid"]          # 1/3 each, no majority
    out = ann.majority_vote(labels, ["0"] * 3, min_prop=0.5)
    assert set(out) == {"Heterogeneous"}


def test_majority_vote_default_min_prop_matches_celltypist_plain_mode():
    labels = ["Tcell", "NK", "NK"]
    out = ann.majority_vote(labels, ["0"] * 3)
    assert set(out) == {"NK"}


def test_majority_vote_is_stable_across_chunking():
    """The point of the reimplementation: chunked prediction must give the same
    answer as one pass, or the memory fix has changed the result."""
    rng = np.random.default_rng(0)
    labels = list(rng.choice(["Tcell", "NK", "Myeloid"], size=300))
    clusters = list(rng.choice(["0", "1", "2"], size=300))
    whole = ann.majority_vote(labels, clusters)
    chunked = ann.majority_vote(
        [l for c in (labels[:100], labels[100:200], labels[200:]) for l in c],
        [c for g in (clusters[:100], clusters[100:200], clusters[200:]) for c in g],
    )
    assert list(whole) == list(chunked)


def test_majority_vote_empty_input():
    assert len(ann.majority_vote([], [])) == 0


# --------------------------------------------------------------------------
# v2 — revised lineage definitions
# --------------------------------------------------------------------------

def test_erythroid_panel_is_not_globin_driven():
    """v1's Erythroid was ("HBB","GYPA"); stage 04 showed HBB is broadly ambient."""
    ery = config.MARKER_PANEL["Erythroid"]
    assert "HBB" not in ery, "HBB must not drive an erythroid assignment"
    specific = {"GYPA", "AHSP", "ALAS2", "CA1"}
    assert len(specific & set(ery)) >= 3, "panel must be led by lineage-specific evidence"
    assert len(ery) >= 5, "a 2-gene panel is what made Erythroid runner-up in 18/30 clusters"


def test_erythroid_contradiction_program_excludes_every_globin():
    """Globins are the dominant ambient species here; they may support an ID but must
    never be able to ACCUSE a cell of being erythroid."""
    prog = set(config.LINEAGE_PROGRAMS["erythroid"])
    assert not (prog & {"HBB", "HBA1", "HBA2", "HBD"})


def test_myeloid_contradiction_program_excludes_lyz():
    """LYZ is the other big ambient species in marrow; same argument."""
    assert "LYZ" not in config.LINEAGE_PROGRAMS["myeloid"]
    assert "LYZ" in config.MARKER_PANEL["Myeloid"]      # still fine for positive ID


def test_t_panel_is_the_tcr_complex():
    t = set(config.MARKER_PANEL["Tcell"])
    assert {"CD3D", "CD3E", "TRAC"} <= t
    assert not (t & {"CD4", "CD8A"}), "CD4 is on monocytes, CD8A on NK/DC subsets"


def test_nk_panel_adds_an_nk_restricted_gene_beyond_the_cytotoxic_program():
    """NKG7/GNLY are a cytotoxic-granule program shared with cytotoxic T cells."""
    nk = set(config.MARKER_PANEL["NK"])
    assert nk & {"KLRF1", "NCAM1", "KLRD1"}, "need NK-restricted evidence, not just cytotoxic"


def test_contradiction_pairs_respect_b_plasma_lineage_relationship():
    """Plasma cells ARE B-lineage — the plasmablast continuum is real biology, not a
    contradiction, and flagging it would be an error."""
    assert "B" not in config.CONTRADICTION_PAIRS["PlasmaCell"]
    assert "plasma" not in config.CONTRADICTION_PAIRS["Bcell"]


def test_hspc_has_no_contradictions_by_design():
    """Progenitors legitimately co-express lineage-priming programs."""
    assert config.CONTRADICTION_PAIRS["HSPC"] == ()


def test_nk_and_erythroid_are_contradicted_by_t_lineage():
    """The two pairs this whole revision exists for."""
    assert "T" in config.CONTRADICTION_PAIRS["NK"]
    assert "T" in config.CONTRADICTION_PAIRS["Erythroid"]


# --------------------------------------------------------------------------
# v2 — exclusivity statistic
# --------------------------------------------------------------------------

def _lineage_adata(labels, per_cell_genes: list[dict[str, float]]):
    """Build an AnnData over every LINEAGE_PROGRAMS gene; unlisted genes are 0."""
    import anndata as ad

    genes = sorted({g for gs in config.LINEAGE_PROGRAMS.values() for g in gs})
    X = np.array([[cell.get(g, 0.0) for g in genes] for cell in per_cell_genes], dtype=float)
    a = ad.AnnData(X)
    a.var_names = genes
    a.obs_names = [f"c{i}" for i in range(a.n_obs)]
    a.obs["lab"] = pd.Categorical(labels)
    return a


def test_cytotoxic_t_cells_are_not_accepted_as_nk_on_cytotoxic_markers_alone():
    """THE case this revision exists for.

    Cells that are NKG7/GNLY-high AND carry strong CD3/TRAC evidence are cytotoxic
    T cells. They pass NK marker coverage (they really are cytotoxic), so only the
    contradiction check can object — and it must.
    """
    cells = [{"NKG7": 9, "GNLY": 9, "CD3D": 3, "TRAC": 3} for _ in range(10)]
    a = _lineage_adata(["NK"] * 10, cells)
    contra = ann.contradiction_rate(a, "lab", classes=("NK",))
    assert contra.loc["NK", "contradiction_rate"] == pytest.approx(1.0)
    assert contra.loc["NK", "contradiction_rate"] > config.CONTRADICTION_MAX_RATE


def test_genuine_nk_cells_pass_without_requiring_tcr_absence():
    """Dropout makes absence unreliable, so the rule is positive-evidence only.
    A true NK cell with ONE stray TCR transcript must still pass."""
    cells = [{"NKG7": 9, "GNLY": 9, "KLRD1": 5, "CD3D": 1} for _ in range(10)]
    a = _lineage_adata(["NK"] * 10, cells)
    contra = ann.contradiction_rate(a, "lab", classes=("NK",))
    assert contra.loc["NK", "contradiction_rate"] == pytest.approx(0.0)


def test_broad_ambient_hbb_alone_cannot_create_an_erythroid_call():
    """Every cell carries ambient haemoglobin. That must produce no erythroid
    evidence at all, because globins are excluded from the erythroid program."""
    cells = [{"HBB": 5, "HBA1": 5, "HBA2": 5} for _ in range(10)]
    a = _lineage_adata(["Tcell"] * 10, cells)
    ev = ann.lineage_evidence(a)
    assert not ev["erythroid"].any(), "globins must not constitute erythroid evidence"


def test_true_erythroid_evidence_is_detected():
    cells = [{"GYPA": 4, "AHSP": 4, "ALAS2": 3} for _ in range(5)]
    a = _lineage_adata(["Tcell"] * 5, cells)
    assert ann.lineage_evidence(a)["erythroid"].all()


def test_min_genes_requires_two_independent_genes():
    """One transcript can be ambient or mismapped; two of a complex is much harder
    to explain that way."""
    one = _lineage_adata(["NK"], [{"CD3D": 9}])
    two = _lineage_adata(["NK"], [{"CD3D": 9, "TRAC": 1}])
    assert not ann.lineage_evidence(one)["T"].iloc[0]
    assert ann.lineage_evidence(two)["T"].iloc[0]


def test_hspc_never_flagged_however_mixed_its_programs():
    cells = [{"CD3D": 5, "TRAC": 5, "GYPA": 5, "AHSP": 5, "CD14": 5, "FCN1": 5}] * 5
    a = _lineage_adata(["HSPC"] * 5, cells)
    contra = ann.contradiction_rate(a, "lab", classes=("HSPC",))
    assert contra.loc["HSPC", "contradiction_rate"] == 0.0


# --------------------------------------------------------------------------
# v2 — the second veto in the decision rule
# --------------------------------------------------------------------------

def test_class_passing_coverage_can_still_fail_contradictory_lineage_validation():
    """Exactly v1's NK: coverage 1.00, concordance clears the bar, but 63% of the
    cells carry T-lineage evidence. Coverage cannot see this; exclusivity must."""
    conc = {"celltypist": _tbl({"NK": 0.99}, "f1")}
    cov = {"celltypist": _tbl({"NK": 1.00}, "coverage"), "manual": _tbl({"NK": 1.00}, "coverage")}
    contra = {"celltypist": _tbl({"NK": 0.63}, "contradiction_rate"),
              "manual": _tbl({"NK": 0.63}, "contradiction_rate")}
    d = ann.decide_per_class(conc, cov, contra, classes=("NK",))
    assert d.loc["NK", "chosen_method"] == AMB
    assert "lineage-exclusivity" in d.loc["NK", "reason"]
    assert d.loc["NK", "vetoed_by_contradiction"] == "celltypist"


def test_new_validation_leaves_classes_without_contradictory_evidence_untouched():
    """A clean class must decide identically with and without the new check."""
    conc = {"celltypist": _tbl({"PlasmaCell": 0.98}, "f1")}
    cov = {"celltypist": _tbl({"PlasmaCell": 1.0}, "coverage"),
           "manual": _tbl({"PlasmaCell": 1.0}, "coverage")}
    contra = {"celltypist": _tbl({"PlasmaCell": 0.02}, "contradiction_rate"),
              "manual": _tbl({"PlasmaCell": 0.02}, "contradiction_rate")}
    before = ann.decide_per_class(conc, cov, classes=("PlasmaCell",))
    after = ann.decide_per_class(conc, cov, contra, classes=("PlasmaCell",))
    assert before.loc["PlasmaCell", "chosen_method"] == after.loc["PlasmaCell", "chosen_method"] == "celltypist"


def test_not_evaluable_contradiction_is_neither_pass_nor_veto():
    """NaN is NOT_EVALUABLE: it cannot win the class, and it is not recorded as vetoed."""
    conc = {"celltypist": _tbl({"HSPC": 0.99}, "f1")}
    cov = {"celltypist": _tbl({"HSPC": 1.0}, "coverage"), "manual": _tbl({"HSPC": 1.0}, "coverage")}
    contra = {"celltypist": _tbl({"HSPC": float("nan")}, "contradiction_rate"),
              "manual": _tbl({"HSPC": 0.0}, "contradiction_rate")}
    d = ann.decide_per_class(conc, cov, contra, classes=("HSPC",))
    assert d.loc["HSPC", "contra_state_celltypist"] == ann.NOT_EVALUABLE
    assert d.loc["HSPC", "chosen_method"] != "celltypist"
    assert "celltypist" not in d.loc["HSPC", "vetoed_by_contradiction"]


def test_fallback_must_also_clear_the_exclusivity_veto():
    """v1's structural flaw: manual was both the concordance reference and the
    fallback, so a wrong manual labelling could not be escaped. It must now clear
    the same vetoes as anyone else."""
    conc = {"celltypist": _tbl({"Erythroid": 0.62}, "f1")}      # below its 0.85 bar
    cov = {"celltypist": _tbl({"Erythroid": 1.0}, "coverage"),
           "manual": _tbl({"Erythroid": 1.0}, "coverage")}
    contra = {"celltypist": _tbl({"Erythroid": 0.05}, "contradiction_rate"),
              "manual": _tbl({"Erythroid": 0.47}, "contradiction_rate")}
    d = ann.decide_per_class(conc, cov, contra, classes=("Erythroid",))
    assert d.loc["Erythroid", "chosen_method"] == AMB
    assert "fallback also failed" in d.loc["Erythroid", "reason"]


def test_acceptance_thresholds_are_untouched_by_the_v2_revision():
    """This revision fixes the reference and the validation, NOT the bars."""
    assert config.CONCORDANCE_THRESHOLDS["PlasmaCell"] == 0.95
    assert all(config.CONCORDANCE_THRESHOLDS[c] == 0.90 for c in ("Tcell", "NK", "Myeloid"))
    assert all(config.CONCORDANCE_THRESHOLDS[c] == 0.85 for c in ("Bcell", "Erythroid", "HSPC"))
    assert config.MARKER_COVERAGE_MIN == 0.30


# --------------------------------------------------------------------------
# contradiction concentration (reporting only — cannot change a decision)
# --------------------------------------------------------------------------

def test_concentration_distinguishes_one_bad_cluster_from_diffuse_noise():
    """The whole point: the class-level rate is identical in both cases here."""
    tc = {"CD3D": 5, "TRAC": 5}
    clean = {"NKG7": 5}
    # concentrated: all 4 contradictory cells in cluster "1"
    conc_cells = [clean] * 8 + [tc] * 4
    conc_clust = ["0"] * 6 + ["1"] * 6
    # diffuse: the same 4 spread across both clusters
    diff_cells = [clean] * 4 + [tc] * 2 + [clean] * 4 + [tc] * 2
    diff_clust = ["0"] * 6 + ["1"] * 6

    for cells, clust, expect_top in ((conc_cells, conc_clust, 1.0), (diff_cells, diff_clust, 0.5)):
        a = _lineage_adata(["NK"] * 12, cells)
        a.obs["clust"] = pd.Categorical(clust)
        rate = ann.contradiction_rate(a, "lab", classes=("NK",)).loc["NK", "contradiction_rate"]
        assert rate == pytest.approx(4 / 12)          # identical class-level rate
        top = ann.contradiction_concentration(a, "lab", "clust", classes=("NK",))
        assert top["share_of_class_contradictions"].max() == pytest.approx(expect_top)


def test_concentration_omits_classes_with_no_contradictions():
    a = _lineage_adata(["NK"] * 4, [{"NKG7": 5}] * 4)
    a.obs["clust"] = pd.Categorical(["0"] * 4)
    assert len(ann.contradiction_concentration(a, "lab", "clust", classes=("NK",))) == 0


def test_concentration_shares_sum_to_one_per_class():
    a = _lineage_adata(["NK"] * 6, [{"CD3D": 5, "TRAC": 5}] * 3 + [{"NKG7": 5}] * 3)
    a.obs["clust"] = pd.Categorical(["0", "1", "2", "0", "1", "2"])
    out = ann.contradiction_concentration(a, "lab", "clust", classes=("NK",))
    assert out["share_of_class_contradictions"].sum() == pytest.approx(1.0)


# --------------------------------------------------------------------------
# v3 — evidence-based adjudication replaces cross-panel score_genes argmax
# --------------------------------------------------------------------------

def _cluster_adata(spec: dict[str, dict[str, float]], n_per: int = 200):
    """One synthetic cluster per key. `spec[cluster][gene] = detection fraction`.

    Genes not named are absent from every cell. Detection is deterministic: the first
    `frac * n_per` cells of the cluster express the gene, so the tests assert on exact
    fractions rather than on sampling noise.
    """
    import anndata as ad

    genes = sorted({g for gs in config.MARKER_PANEL.values() for g in gs}
                   | {g for gs in config.LINEAGE_PROGRAMS.values() for g in gs}
                   | {g for gs in config.MYELOID_SUBPROGRAMS.values() for g in gs}
                   | set(config.MYELOID_MONO_ANCHORS) | set(config.MYELOID_MONO_CONTEXT)
                   | set(config.MYELOID_DC_ANCHORS) | set(config.MYELOID_DC_CONTEXT)
                   | set(config.MYELOID_PDC_CORE)
                   | set(config.MYELOID_SUPPORTING) | set(config.PDC_CONTEXTUAL)
                   | set(config.HSPC_CORE) | set(config.PLASMA_SECRETORY)
                   | set(config.PLASMA_MATURE)
                   | {g for gs in config.HSPC_CONTEXT.values() for g in gs})
    blocks, labels = [], []
    for cl, fracs in spec.items():
        M = np.zeros((n_per, len(genes)))
        for gi, g in enumerate(genes):
            k = int(round(fracs.get(g, 0.0) * n_per))
            M[:k, gi] = 5.0
        blocks.append(M)
        labels += [cl] * n_per
    a = ad.AnnData(np.vstack(blocks))
    a.var_names = genes
    a.obs_names = [f"c{i}" for i in range(a.n_obs)]
    a.obs["clust"] = pd.Categorical(labels)
    return a


T_STRONG = {"CD3D": 0.87, "CD3E": 0.72, "CD3G": 0.55, "TRAC": 0.78,
            "TRBC1": 0.50, "TRBC2": 0.66}
NK_STRONG = {"KLRD1": 0.90, "KLRF1": 0.75, "NCAM1": 0.40, "NKG7": 0.99, "GNLY": 0.94}
CYTOTOXIC = {"NKG7": 0.95, "GNLY": 0.60, "KLRD1": 0.50}


def test_cross_panel_score_genes_magnitude_cannot_determine_identity():
    """The v1/v2 defect: identity must not be decided by module-score magnitude.

    `adjudicate_clusters` reads detection fractions only. Planting wildly divergent
    score_genes columns must not change a single assignment.
    """
    a = _cluster_adata({"x": {**T_STRONG}})
    before = ann.adjudicate_clusters(a, "clust")["winner"].to_dict()
    a.obs["score_NK"] = 99.0        # would have won any argmax
    a.obs["score_Tcell"] = -99.0
    after = ann.adjudicate_clusters(a, "clust")["winner"].to_dict()
    assert before == after == {"x": "Tcell"}


def test_strong_cd3_supports_t_identity_despite_high_cytotoxic_genes():
    """Clusters 3/12 in kind: full TCR complex AND a high cytotoxic program."""
    a = _cluster_adata({"x": {**T_STRONG, **CYTOTOXIC}})
    tbl = ann.adjudicate_clusters(a, "clust")
    assert tbl.loc["x", "winner"] == "Tcell"
    # NK was disqualified by T-lineage evidence, not merely out-scored
    assert "NK" not in tbl.loc["x", "survivors"]


def test_cytotoxic_genes_alone_cannot_make_a_t_population_nk():
    a = _cluster_adata({"x": {**T_STRONG, "NKG7": 0.99, "GNLY": 0.99}})
    assert ann.adjudicate_clusters(a, "clust").loc["x", "winner"] != "NK"


def test_genuine_nk_supported_without_requiring_absence_of_t_transcripts():
    """A true NK cluster carrying scattered TCR transcripts must still be called NK."""
    a = _cluster_adata({"x": {**NK_STRONG, "CD3D": 0.10, "TRAC": 0.08}})
    tbl = ann.adjudicate_clusters(a, "clust")
    assert tbl.loc["x", "winner"] == "NK"
    assert tbl.loc["x", "excluded_NK"] <= config.CONTRADICTION_MAX_RATE


def test_cytotoxic_genes_alone_cannot_make_a_t_population_nk():
    a = _cluster_adata({"x": {**T_STRONG, "NKG7": 0.99, "GNLY": 0.99}})
    assert ann.adjudicate_clusters(a, "clust").loc["x", "winner"] != "NK"


def test_genuine_nk_supported_without_requiring_absence_of_t_transcripts():
    """A true NK cluster carrying scattered TCR transcripts must still be called NK."""
    a = _cluster_adata({"x": {**NK_STRONG, "CD3D": 0.10, "TRAC": 0.08}})
    tbl = ann.adjudicate_clusters(a, "clust")
    assert tbl.loc["x", "winner"] == "NK"
    assert tbl.loc["x", "excluded_NK"] <= config.CONTRADICTION_MAX_RATE


def test_mutually_contradictory_lineages_both_disqualified():
    """A cluster that is both B and myeloid is a doublet cluster, not a cell type.
    Each lineage contradicts the other, both are disqualified, and the result is
    Ambiguous with the reason stated."""
    # Built from genes the myeloid CONTRADICTION program actually reads
    # (CD14/FCN1/MNDA/ITGAM), so both lineages genuinely accuse each other. Note the
    # contradiction program is deliberately narrower than the identification
    # subprograms — it excludes ambient-prone genes — so a monocyte program lacking
    # these four raises no myeloid contradiction. That asymmetry is intentional and
    # is reported with the Checkpoint-2 results.
    a = _cluster_adata({"x": {"MS4A1": 0.80, "CD79A": 0.80, "CD19": 0.80,
                              "CD14": 0.80, "FCN1": 0.80, "MNDA": 0.80,
                              "ITGAM": 0.80, "LST1": 0.80, "FCER1G": 0.80,
                              "AIF1": 0.80, "CTSS": 0.80, "TYROBP": 0.80,
                              "CSF1R": 0.80, "VCAN": 0.80, "MS4A7": 0.80,
                              "SERPINA1": 0.80}})
    tbl = ann.adjudicate_clusters(a, "clust")
    assert tbl.loc["x", "positive_Bcell"] == pytest.approx(1.0)
    assert tbl.loc["x", "myeloid_subprograms_passing"]      # a real monocyte program
    assert tbl.loc["x", "winner"] == AMB
    # Either route -- margin, or the strongest hypothesis being disqualified by the
    # other lineage -- is a correct refusal to call a doublet cluster a cell type.
    assert ("within margin" in tbl.loc["x", "reason"]
            or "disqualified" in tbl.loc["x", "reason"])


def test_balanced_competing_evidence_returns_ambiguous_not_a_tie_break():
    """The margin route, constructed independently of the cluster-3/12 failure.

    PlasmaCell and Bcell are DELIBERATELY not contradictory (plasma cells are B-lineage;
    the plasmablast continuum is real), so neither can disqualify the other and the
    decision falls to the margin alone. With equal positive evidence the answer must be
    Ambiguous — not a win for whichever panel has more markers, higher raw detection, or
    comes first in iteration order. This is the test that the replacement has not simply
    swapped one calibration artifact for another.
    """
    a = _cluster_adata({"x": {"SDC1": 0.80, "CD38": 0.80, "MZB1": 0.80,
                              "XBP1": 0.80, "IRF4": 0.80,
                              "MS4A1": 0.80, "CD79A": 0.80, "CD19": 0.80}})
    tbl = ann.adjudicate_clusters(a, "clust")
    assert tbl.loc["x", "positive_PlasmaCell"] == pytest.approx(1.0)
    assert tbl.loc["x", "positive_Bcell"] == pytest.approx(1.0)
    assert "PlasmaCell" in tbl.loc["x", "survivors"] and "Bcell" in tbl.loc["x", "survivors"]
    assert tbl.loc["x", "winner"] == AMB
    assert "within margin" in tbl.loc["x", "reason"]


def test_ambiguous_is_not_awarded_to_the_panel_with_more_markers():
    """Tcell has 6 markers, Bcell 3. Equal fractional support must not favour Tcell."""
    a = _cluster_adata({"x": {**T_STRONG, "MS4A1": 0.80, "CD79A": 0.80, "CD19": 0.80}})
    tbl = ann.adjudicate_clusters(a, "clust")
    assert tbl.loc["x", "positive_Tcell"] == pytest.approx(1.0)
    assert tbl.loc["x", "positive_Bcell"] == pytest.approx(1.0)
    assert tbl.loc["x", "winner"] == AMB


def test_unsupported_cluster_is_ambiguous_with_a_stated_reason():
    a = _cluster_adata({"x": {"NKG7": 0.99}})       # one shared gene, no program
    tbl = ann.adjudicate_clusters(a, "clust")
    assert tbl.loc["x", "winner"] == AMB
    assert tbl.loc["x", "reason"] == "no class supported"


def test_detection_below_the_declared_floor_does_not_count_as_expressed():
    weak = {g: config.MANUAL_MARKER_DETECT_MIN - 0.05 for g in config.MARKER_PANEL["Tcell"]}
    a = _cluster_adata({"x": weak})
    assert ann.adjudicate_clusters(a, "clust").loc["x", "positive_Tcell"] == 0.0


def test_v3_leaves_every_v2_acceptance_threshold_untouched():
    assert config.CONCORDANCE_THRESHOLDS["PlasmaCell"] == 0.95
    assert config.MARKER_COVERAGE_MIN == 0.30
    assert config.CONTRADICTION_MIN_GENES == 2
    assert config.CONTRADICTION_MAX_RATE == 0.25
    assert config.MARKER_PANEL["NK"] == ("NCAM1", "NKG7", "GNLY", "KLRD1", "KLRF1")


# --------------------------------------------------------------------------
# D — three-state vetoes: NOT_EVALUABLE is never coerced
# --------------------------------------------------------------------------

def test_evaluate_veto_three_states():
    assert ann.evaluate_veto(0.9, 0.3, direction="min") == ann.PASS
    assert ann.evaluate_veto(0.1, 0.3, direction="min") == ann.FAIL
    assert ann.evaluate_veto(float("nan"), 0.3, direction="min") == ann.NOT_EVALUABLE
    assert ann.evaluate_veto(0.1, 0.25, direction="max") == ann.PASS
    assert ann.evaluate_veto(0.9, 0.25, direction="max") == ann.FAIL
    assert ann.evaluate_veto(float("nan"), 0.25, direction="max") == ann.NOT_EVALUABLE
    assert ann.evaluate_veto(None, 0.3, direction="min") == ann.NOT_EVALUABLE


def test_not_evaluable_coverage_cannot_make_a_winner():
    """Novershtern/PlasmaCell shape: the reference has no such class."""
    conc = {"singler_nov": _tbl({"PlasmaCell": 0.99}, "f1")}
    cov = {"singler_nov": _tbl({"PlasmaCell": float("nan")}, "coverage"),
           "manual": _tbl({"PlasmaCell": 1.0}, "coverage")}
    contra = {"singler_nov": _tbl({"PlasmaCell": float("nan")}, "contradiction_rate"),
              "manual": _tbl({"PlasmaCell": 0.05}, "contradiction_rate")}
    d = ann.decide_per_class(conc, cov, contra, classes=("PlasmaCell",))
    assert d.loc["PlasmaCell", "chosen_method"] == "manual"
    assert d.loc["PlasmaCell", "cov_state_singler_nov"] == ann.NOT_EVALUABLE
    assert "singler_nov" not in d.loc["PlasmaCell", "vetoed_by_coverage"]


def test_not_evaluable_fallback_is_reported_as_such_not_as_a_failure():
    """manual/NK shape: the manual classifier produces no NK population."""
    conc = {"celltypist": _tbl({"NK": 0.50}, "f1")}
    cov = {"celltypist": _tbl({"NK": 1.0}, "coverage"),
           "manual": _tbl({"NK": float("nan")}, "coverage")}
    contra = {"celltypist": _tbl({"NK": 0.56}, "contradiction_rate"),
              "manual": _tbl({"NK": float("nan")}, "contradiction_rate")}
    d = ann.decide_per_class(conc, cov, contra, classes=("NK",))
    assert d.loc["NK", "chosen_method"] == AMB
    assert "NOT_EVALUABLE" in d.loc["NK", "reason"]
    assert "assigns no cells" in d.loc["NK", "reason"]


def test_not_evaluable_is_listed_explicitly_in_the_decision_table():
    conc = {"m1": _tbl({"NK": 0.99}, "f1")}
    cov = {"m1": _tbl({"NK": float("nan")}, "coverage"), "manual": _tbl({"NK": 1.0}, "coverage")}
    contra = {"m1": _tbl({"NK": 0.01}, "contradiction_rate"),
              "manual": _tbl({"NK": 0.01}, "contradiction_rate")}
    d = ann.decide_per_class(conc, cov, contra, classes=("NK",), fallback="manual")
    assert "m1:cov" in d.loc["NK", "not_evaluable"]


def test_not_evaluable_contradiction_no_longer_silently_passes():
    """v3 read NaN contradiction as pass. It must now be neither."""
    conc = {"m1": _tbl({"NK": 0.99}, "f1")}
    cov = {"m1": _tbl({"NK": 1.0}, "coverage"), "manual": _tbl({"NK": float("nan")}, "coverage")}
    contra = {"m1": _tbl({"NK": float("nan")}, "contradiction_rate")}
    d = ann.decide_per_class(conc, cov, contra, classes=("NK",))
    assert d.loc["NK", "contra_state_m1"] == ann.NOT_EVALUABLE
    assert d.loc["NK", "chosen_method"] != "m1"      # cannot win on unevaluated evidence


def test_A_cluster_3_12_shape_strongest_is_t_and_t_passes_exclusion():
    """T strongest and uncontradicted, NK carrying substantial shared-cytotoxic
    support: resolves to Tcell. This is clusters 3 and 12."""
    a = _cluster_adata({"x": {**T_STRONG, **CYTOTOXIC}})
    tbl = ann.adjudicate_clusters(a, "clust")
    assert tbl.loc["x", "strongest_supported"] == "Tcell"
    assert tbl.loc["x", "winner"] == "Tcell"
    assert "NK" in tbl.loc["x", "supported"]
    assert "NK" not in tbl.loc["x", "survivors"]


def test_B_cluster_23_shape_disqualified_strongest_does_not_cascade():
    """NK strongest but contradicted; T has real, uncontradicted support.

    The result must be Ambiguous. This is the specific regression: v3 cascaded to T
    and absorbed a mixed NK/gamma-delta population into a clean label.
    """
    a = _cluster_adata({"x": {**NK_STRONG, "CD3E": 0.30, "TRBC1": 0.61, "TRBC2": 0.47}})
    tbl = ann.adjudicate_clusters(a, "clust")
    assert tbl.loc["x", "strongest_supported"] == "NK"
    assert tbl.loc["x", "positive_NK"] > tbl.loc["x", "positive_Tcell"]
    assert tbl.loc["x", "excluded_Tcell"] <= config.CONTRADICTION_MAX_RATE   # T itself is fine
    assert tbl.loc["x", "winner"] == AMB
    assert "not cascaded to runner-up" in tbl.loc["x", "reason"]


def test_C_no_supported_lineage_is_ambiguous_via_insufficient_support():
    a = _cluster_adata({"x": {"NKG7": 0.99}})
    tbl = ann.adjudicate_clusters(a, "clust")
    assert tbl.loc["x", "winner"] == AMB
    assert tbl.loc["x", "reason"] == "no class supported"


def test_nk_markers_on_a_tcr_positive_cluster_with_equal_support_is_ambiguous():
    """Full TCR AND full NK panel at high detection is genuinely mixed evidence;
    under the no-cascade rule it is unresolved rather than adjudicated to either."""
    a = _cluster_adata({"x": {**NK_STRONG, **T_STRONG}})
    assert ann.adjudicate_clusters(a, "clust").loc["x", "winner"] == AMB


# --------------------------------------------------------------------------
# Checkpoint 1 — final-label assembly
#
# The assembly layer resolves competing outputs from decisions ALREADY MADE. It never
# invents marker evidence and never alters a class-level threshold. Its only inputs are
# the per-method labels and the class-level decision table.
# --------------------------------------------------------------------------

def _assembly_adata(labels: dict[str, str], n: int = 10):
    """One synthetic cluster; `labels` maps method -> the class that method assigns."""
    import anndata as ad

    a = ad.AnnData(np.zeros((n, 2)))
    a.var_names = ["g1", "g2"]
    a.obs_names = [f"c{i}" for i in range(n)]
    for method, lab in labels.items():
        a.obs[f"cell_type_{method}"] = pd.Categorical([lab] * n)
    return a


def _decision(winners: dict[str, str]):
    return pd.DataFrame({"chosen_method": winners}).rename_axis("cell_class")


_KEYS = {m: f"cell_type_{m}" for m in ("manual", "celltypist", "singler_nov", "singler_hpca")}


def test_assembly_case1_broad_consensus_is_not_lost_to_class_winner_bookkeeping():
    """Clusters 22/25 shape: three methods agree Myeloid, manual (the class winner
    for Myeloid) says Ambiguous. The cluster must not go ownerless purely because of
    which method happened to win the class-level decision."""
    a = _assembly_adata({"manual": AMB, "celltypist": "Myeloid",
                         "singler_nov": "Myeloid", "singler_hpca": "Myeloid"})
    dec = _decision({"Myeloid": "manual", "PlasmaCell": "celltypist",
                     "HSPC": "singler_nov", "NK": AMB})
    ann.assemble_final_labels(a, dec, _KEYS)
    assert set(a.obs["cell_type"]) == {"Myeloid"}
    assert set(a.obs["annotation_source"]) == {"consensus"}


def test_assembly_case2_cross_claim_has_a_deterministic_documented_path():
    """Cluster-26 shape: the winner of one class points elsewhere and the winner of
    another points back. Resolution must be explicit, not accidental."""
    a = _assembly_adata({"manual": "PlasmaCell", "celltypist": "Myeloid",
                         "singler_nov": "Myeloid", "singler_hpca": "Myeloid"})
    dec = _decision({"Myeloid": "manual", "PlasmaCell": "celltypist"})
    ann.assemble_final_labels(a, dec, _KEYS)
    # 3 of 4 methods say Myeloid -> strict majority of expressed opinions
    assert set(a.obs["cell_type"]) == {"Myeloid"}


def test_assembly_case3_genuine_cross_claim_can_still_return_ambiguous():
    """Two comparably legitimate claims: assembly must NOT force ownership."""
    a = _assembly_adata({"manual": "Bcell", "celltypist": "Myeloid",
                         "singler_nov": "Myeloid", "singler_hpca": "Bcell"})
    dec = _decision({"Myeloid": "singler_nov", "Bcell": "manual"})
    ann.assemble_final_labels(a, dec, _KEYS)
    assert set(a.obs["cell_type"]) == {AMB}          # 2-2 is not a strict majority


def test_assembly_never_overrides_a_class_level_ambiguous_verdict():
    """Cluster 23 shape: every automated method says NK, but NK was ruled unresolved
    at class level. Assembly must not resurrect it — that would smuggle a biological
    decision into the bookkeeping layer."""
    a = _assembly_adata({"manual": AMB, "celltypist": "NK",
                         "singler_nov": "NK", "singler_hpca": "NK"})
    dec = _decision({"NK": AMB, "Myeloid": "celltypist"})
    ann.assemble_final_labels(a, dec, _KEYS)
    assert set(a.obs["cell_type"]) == {AMB}


def test_assembly_conflicting_authoritative_claims_fall_through_to_consensus():
    """Two winners each claiming their own class is a genuine disagreement, so it is
    NOT resolved by precedence — it falls through to the consensus test. Here consensus
    is unambiguous (3 of 4 say Bcell) and resolves it; where consensus is unclear the
    result is Ambiguous, which the 2-2 and iteration-order tests cover."""
    a = _assembly_adata({"manual": "Myeloid", "celltypist": "Bcell",
                         "singler_nov": "Bcell", "singler_hpca": "Bcell"})
    dec = _decision({"Myeloid": "manual", "Bcell": "singler_nov"})
    ann.assemble_final_labels(a, dec, _KEYS)
    assert set(a.obs["cell_type"]) == {"Bcell"}
    assert set(a.obs["annotation_source"]) == {"consensus"}


def test_assembly_single_authoritative_claim_is_assigned_with_provenance():
    a = _assembly_adata({"manual": "Myeloid", "celltypist": AMB,
                         "singler_nov": AMB, "singler_hpca": AMB})
    dec = _decision({"Myeloid": "manual"})
    ann.assemble_final_labels(a, dec, _KEYS)
    assert set(a.obs["cell_type"]) == {"Myeloid"}
    assert set(a.obs["annotation_source"]) == {"manual"}


def test_assembly_does_not_resolve_ties_by_broad_classes_iteration_order():
    """v4 gave PlasmaCell precedence simply because it is first in BROAD_CLASSES.
    Order must not decide anything."""
    a = _assembly_adata({"manual": "PlasmaCell", "celltypist": "HSPC",
                         "singler_nov": "PlasmaCell", "singler_hpca": "HSPC"})
    dec = _decision({"PlasmaCell": "manual", "HSPC": "celltypist"})
    ann.assemble_final_labels(a, dec, _KEYS)
    assert set(a.obs["cell_type"]) == {AMB}


def test_assembly_invents_no_marker_evidence():
    """The assembly layer reads labels and the decision table only. Changing expression
    must not change its output."""
    a = _assembly_adata({"manual": AMB, "celltypist": "Myeloid",
                         "singler_nov": "Myeloid", "singler_hpca": "Myeloid"})
    dec = _decision({"Myeloid": "manual"})
    ann.assemble_final_labels(a, dec, _KEYS)
    first = list(a.obs["cell_type"])
    a.X = np.ones_like(a.X) * 999.0
    ann.assemble_final_labels(a, dec, _KEYS)
    assert list(a.obs["cell_type"]) == first


def test_assembly_ambiguous_is_not_a_consensus_opinion():
    a = _assembly_adata({"manual": AMB, "celltypist": AMB,
                         "singler_nov": "Myeloid", "singler_hpca": AMB})
    dec = _decision({"Myeloid": "manual"})
    ann.assemble_final_labels(a, dec, _KEYS)
    # 1 of 1 expressed opinion is a strict majority; the three Ambiguous do not dilute it
    assert set(a.obs["cell_type"]) == {"Myeloid"}


def test_assembly_not_evaluable_is_not_a_consensus_opinion():
    a = _assembly_adata({"manual": ann.NOT_EVALUABLE, "celltypist": "Myeloid",
                         "singler_nov": ann.NOT_EVALUABLE, "singler_hpca": ann.NOT_EVALUABLE})
    dec = _decision({"Myeloid": "manual"})
    ann.assemble_final_labels(a, dec, _KEYS)
    assert set(a.obs["cell_type"]) == {"Myeloid"}


def test_assembly_two_agreeing_opinions_with_noise_resolve():
    """T, T, NOT_EVALUABLE, Ambiguous -> T on 2 of 2 expressed opinions."""
    a = _assembly_adata({"manual": "Tcell", "celltypist": "Tcell",
                         "singler_nov": ann.NOT_EVALUABLE, "singler_hpca": AMB})
    dec = _decision({"Tcell": "singler_nov"})     # winner abstains -> consensus route
    ann.assemble_final_labels(a, dec, _KEYS)
    assert set(a.obs["cell_type"]) == {"Tcell"}
    assert set(a.obs["annotation_source"]) == {"consensus"}


def test_assembly_one_one_split_with_noise_is_ambiguous():
    """T, NK, NOT_EVALUABLE, Ambiguous -> 1-1 is not a strict majority."""
    a = _assembly_adata({"manual": "Tcell", "celltypist": "NK",
                         "singler_nov": ann.NOT_EVALUABLE, "singler_hpca": AMB})
    dec = _decision({"Tcell": "singler_nov", "NK": "singler_hpca"})
    ann.assemble_final_labels(a, dec, _KEYS)
    assert set(a.obs["cell_type"]) == {AMB}


# --------------------------------------------------------------------------
# Checkpoint 2 — structured biological reference
# --------------------------------------------------------------------------

def _support(det: dict[str, float]):
    genes = sorted(set(det) | {g for gs in config.MYELOID_SUBPROGRAMS.values() for g in gs}
                   | set(config.MYELOID_MONO_ANCHORS) | set(config.MYELOID_MONO_CONTEXT)
                   | set(config.MYELOID_DC_ANCHORS) | set(config.MYELOID_DC_CONTEXT)
                   | set(config.MYELOID_PDC_CORE) | set(config.PDC_CONTEXTUAL)
                   | set(config.MYELOID_SUPPORTING) | set(config.HSPC_CORE)
                   | set(config.PLASMA_SECRETORY) | set(config.PLASMA_MATURE)
                   | {g for gs in config.MARKER_PANEL.values() for g in gs}
                   | {"IGLL1", "MPO", "GZMB"})
    row = pd.Series({g: det.get(g, 0.0) for g in genes})
    return ann.class_support(row, detect_min=config.MANUAL_MARKER_DETECT_MIN,
                             positive_min=config.MANUAL_POSITIVE_MIN)


HI = 0.80


def test_ontology_places_pdc_under_broad_myeloid_consistently():
    """Declared, not accidental: the CellTypist collapse map has done this since v1."""
    assert ann.CELLTYPIST_TO_BROAD["pDC"] == "Myeloid"
    assert config.MYELOID_PDC_CORE          # pDC is an independent Myeloid route


# --- Myeloid positive -------------------------------------------------------

def test_monocyte_with_low_cd14_is_myeloid():
    """Original intent preserved: CD14 is not required for a monocyte call. Updated
    for the anchor+context architecture — broad context alone is no longer enough, so
    the committed anchors are supplied and CD14 stays below the floor."""
    ev = _support({"LST1": HI, "FCER1G": HI, "CTSS": HI, "AIF1": HI, "TYROBP": HI,
                   "CSF1R": HI, "FCN1": HI, "VCAN": HI, "MS4A7": HI,
                   "CD14": 0.21, "LYZ": HI})
    assert ev["Myeloid"]["supported"]
    assert "monocyte" in ev["Myeloid"]["subprograms_passing"]


def test_mhcII_alone_no_longer_supports_myeloid():
    """SUPERSEDED Checkpoint-2 test, kept inverted as a regression guard.

    It originally asserted that an MHC-II-high, CD14-low cell IS Myeloid — which is
    exactly the case a B cell, a T cell and a progenitor also satisfy, and it is why
    the C2 suite passed while the reference destroyed Bcell and HSPC. Shared
    professional-APC biology must not identify myeloid lineage."""
    ev = _support({"CST3": HI, "CTSS": HI, "HLA-DRA": HI, "HLA-DRB1": HI,
                   "HLA-DPA1": HI, "HLA-DPB1": HI, "CD74": HI, "CD14": 0.10})
    assert not ev["Myeloid"]["supported"]


def test_pdc_program_supports_broad_myeloid_without_monocyte_markers():
    ev = _support({"TCF4": HI, "IRF7": HI, "IRF8": HI, "LILRA4": HI, "PLD4": HI})
    assert ev["Myeloid"]["supported"]
    assert ev["Myeloid"]["subprograms_passing"] == ["pdc"]
    assert ev["Myeloid"]["subprograms"]["monocyte"] < config.MANUAL_POSITIVE_MIN


# --- Myeloid negative -------------------------------------------------------

def test_plasma_cell_with_ambient_lyz_only_is_not_myeloid():
    ev = _support({"LYZ": HI, "MZB1": HI, "XBP1": HI, "SDC1": HI})
    assert not ev["Myeloid"]["supported"]
    assert ev["PlasmaCell"]["supported"]


def test_t_cell_with_ambient_lyz_only_is_not_myeloid():
    ev = _support({"LYZ": HI, "CD3D": HI, "CD3E": HI, "TRAC": HI, "TRBC2": HI})
    assert not ev["Myeloid"]["supported"]


def test_ambient_lyz_heavy_background_cannot_create_myeloid():
    """LYZ detected everywhere, CD14/ITGAM low, no coherent subprogram."""
    ev = _support({"LYZ": 0.99, "CD14": 0.05, "ITGAM": 0.05})
    assert not ev["Myeloid"]["supported"]
    assert ev["Myeloid"]["supporting"] == ["LYZ"]


def test_gzmb_alone_does_not_make_a_population_pdc_or_myeloid():
    ev = _support({"GZMB": HI, "NKG7": HI, "GNLY": HI})
    assert not ev["Myeloid"]["supported"]


# --- HSPC / progenitor ------------------------------------------------------

def test_lymphoid_progenitor_without_mpo_is_hspc():
    ev = _support({"CD34": HI, "SOX4": HI, "HLF": HI, "IGLL1": HI})
    assert ev["HSPC"]["supported"]
    assert ev["HSPC"]["priming"]["lymphoid_primed"] == ["IGLL1"]
    assert ev["HSPC"]["priming"]["myeloid_primed"] == []


def test_myeloid_primed_progenitor_without_igll1_is_hspc():
    ev = _support({"CD34": HI, "GATA2": HI, "SPINK2": HI, "MPO": HI})
    assert ev["HSPC"]["supported"]
    assert ev["HSPC"]["priming"]["myeloid_primed"] == ["MPO"]


def test_mature_b_cell_without_progenitor_program_is_not_hspc():
    ev = _support({"MS4A1": HI, "CD79A": HI, "CD19": HI})
    assert not ev["HSPC"]["supported"]
    assert ev["Bcell"]["supported"]


# --- the mature-plasma predicate -------------------------------------------

def test_progenitor_with_mzb1_xbp1_is_not_stolen_by_plasmacell():
    """Cluster-24 shape. Secretory axis present, mature axis absent."""
    ev = _support({"CD34": HI, "SOX4": HI, "HLF": HI, "IGLL1": HI,
                   "MZB1": HI, "XBP1": HI, "SDC1": 0.02, "TNFRSF17": 0.01})
    assert ev["PlasmaCell"]["secretory_ok"]
    assert not ev["PlasmaCell"]["mature_ok"]
    assert not ev["PlasmaCell"]["supported"]
    assert ev["HSPC"]["supported"]


def test_true_plasma_cell_satisfies_both_axes():
    ev = _support({"MZB1": HI, "XBP1": HI, "SDC1": HI})
    assert ev["PlasmaCell"]["secretory_ok"] and ev["PlasmaCell"]["mature_ok"]
    assert ev["PlasmaCell"]["supported"]


def test_cd38_alone_cannot_satisfy_the_mature_axis():
    """CD38 is on activated T, NK, pro-B and progenitors — not plasma-specific."""
    assert "CD38" not in config.PLASMA_MATURE
    ev = _support({"MZB1": HI, "XBP1": HI, "CD38": HI, "SDC1": 0.0, "TNFRSF17": 0.0})
    assert not ev["PlasmaCell"]["supported"]


def test_secretory_axis_alone_cannot_satisfy_the_predicate():
    ev = _support({"MZB1": HI, "XBP1": HI})
    assert ev["PlasmaCell"]["secretory_ok"]
    assert not ev["PlasmaCell"]["supported"]


def test_mature_axis_alone_cannot_satisfy_the_predicate():
    ev = _support({"SDC1": HI, "MZB1": 0.05, "XBP1": 0.05})
    assert ev["PlasmaCell"]["mature_ok"]
    assert not ev["PlasmaCell"]["secretory_ok"]
    assert not ev["PlasmaCell"]["supported"]


def test_mixed_progenitor_plasma_without_clean_separation_is_ambiguous():
    """Both fully supported and within the margin -> unresolved, not forced."""
    a = _cluster_adata({"x": {"CD34": HI, "SOX4": HI, "HLF": HI,
                              "MZB1": HI, "XBP1": HI, "SDC1": HI}})
    tbl = ann.adjudicate_clusters(a, "clust")
    assert tbl.loc["x", "winner"] == AMB


def test_checkpoint2_changes_no_threshold():
    assert config.MANUAL_POSITIVE_MIN == 0.5
    assert config.MANUAL_MARKER_DETECT_MIN == 0.25
    assert config.MANUAL_DECISION_MARGIN == 0.15
    assert config.MARKER_COVERAGE_MIN == 0.30
    assert config.CONTRADICTION_MIN_GENES == 2
    assert config.CONTRADICTION_MAX_RATE == 0.25


# --------------------------------------------------------------------------
# Checkpoint 2b — conventional DC: anchor + context
# --------------------------------------------------------------------------

MHCII = {"HLA-DRA": HI, "HLA-DRB1": HI, "HLA-DPA1": HI, "HLA-DPB1": HI,
         "CD74": HI, "CST3": HI, "CTSS": HI}
DC_ANCHOR = {"FCER1G": HI, "TYROBP": HI, "LST1": HI, "AIF1": HI}


def test_mhcII_high_b_cell_without_anchor_is_not_myeloid():
    """Clusters 10/21: B cells are professional APCs and scored 0.86 on the flat list."""
    ev = _support({**MHCII, "MS4A1": HI, "CD79A": HI, "CD19": HI})
    assert not ev["Myeloid"]["supported"]
    assert ev["Bcell"]["supported"]


def test_mhcII_high_progenitor_without_anchor_is_not_myeloid():
    ev = _support({**MHCII, "CD34": HI, "SOX4": HI, "HLF": HI})
    assert not ev["Myeloid"]["supported"]
    assert ev["HSPC"]["supported"]


def test_mhcII_high_t_cell_without_anchor_is_not_myeloid():
    """Cluster 3: the exact regression this redesign exists to undo."""
    ev = _support({**MHCII, "CD3D": HI, "CD3E": HI, "CD3G": HI,
                   "TRAC": HI, "TRBC1": HI, "TRBC2": HI})
    assert not ev["Myeloid"]["supported"]
    assert ev["Tcell"]["supported"]


def test_conventional_dc_with_anchor_and_context_is_myeloid():
    ev = _support({**MHCII, **DC_ANCHOR, "CD1C": HI, "FCER1A": HI})
    assert ev["Myeloid"]["supported"]
    assert "dc" in ev["Myeloid"]["subprograms_passing"]


def test_ambient_lyz_plus_mhcII_without_anchor_is_insufficient():
    ev = _support({**MHCII, "LYZ": 0.99, "CD14": 0.10})
    assert not ev["Myeloid"]["supported"]


def test_pdc_does_not_pass_through_the_mhcII_dc_route():
    """pDCs are MHC-II-low; the pDC route must stand alone."""
    ev = _support({"TCF4": HI, "LILRA4": HI, "CLEC4C": HI, "IL3RA": HI, "PLD4": HI})
    assert ev["Myeloid"]["supported"]
    assert "pdc" in ev["Myeloid"]["subprograms_passing"]
    assert "dc" not in ev["Myeloid"]["subprograms_passing"]


def test_irf7_irf8_gzmb_are_contextual_and_cannot_carry_pdc():
    ev = _support({"IRF7": HI, "IRF8": HI, "GZMB": HI})
    assert not ev["Myeloid"]["supported"]


def test_dc_route_is_the_weaker_of_its_two_axes():
    anchor_only = _support(DC_ANCHOR)
    assert anchor_only["Myeloid"]["subprograms"]["dc"] == 0.0


# --------------------------------------------------------------------------
# Checkpoint 2c — monocyte route: committed anchor + broad context
# --------------------------------------------------------------------------

MONO_CONTEXT = {"FCER1G": HI, "TYROBP": HI, "LST1": HI, "AIF1": HI, "CTSS": HI, "LYZ": HI}
MONO_ANCHOR = {"CSF1R": HI, "FCN1": HI, "VCAN": HI, "MS4A7": HI, "SERPINA1": HI}
PROGENITOR = {"CD34": HI, "SOX4": HI, "SPINK2": HI}


def test_progenitor_with_broad_innate_context_is_not_myeloid():
    """Leiden-24 shape: real progenitor evidence plus broad innate machinery, but no
    committed-monocyte anchor. Broad innate biology is not lineage commitment."""
    ev = _support({**PROGENITOR, **MONO_CONTEXT})
    assert not ev["Myeloid"]["supported"]
    assert ev["HSPC"]["supported"]


def test_committed_monocyte_with_anchor_and_context_is_myeloid():
    ev = _support({**MONO_ANCHOR, **MONO_CONTEXT})
    assert ev["Myeloid"]["supported"]
    assert "monocyte" in ev["Myeloid"]["subprograms_passing"]


def test_broad_myeloid_context_alone_is_insufficient():
    ev = _support(MONO_CONTEXT)
    assert not ev["Myeloid"]["supported"]
    assert ev["Myeloid"]["subprograms"]["monocyte"] == 0.0


def test_ambient_lyz_plus_broad_context_alone_is_insufficient():
    ev = _support({"LYZ": 0.99, "FCER1G": HI, "TYROBP": HI, "CTSS": HI})
    assert not ev["Myeloid"]["supported"]


def test_committed_anchor_alone_without_context_is_insufficient():
    """The route is the weaker of its two axes, symmetric with DC."""
    ev = _support(MONO_ANCHOR)
    assert ev["Myeloid"]["subprograms"]["monocyte"] == 0.0


def test_dc_route_still_works_after_the_monocyte_redesign():
    ev = _support({**MHCII, **DC_ANCHOR, "CD1C": HI, "FCER1A": HI})
    assert ev["Myeloid"]["supported"]
    assert "dc" in ev["Myeloid"]["subprograms_passing"]


def test_pdc_route_still_works_after_the_monocyte_redesign():
    ev = _support({"TCF4": HI, "LILRA4": HI, "CLEC4C": HI, "IL3RA": HI, "PLD4": HI})
    assert ev["Myeloid"]["supported"]
    assert "pdc" in ev["Myeloid"]["subprograms_passing"]


def test_cluster24_shaped_case_resolves_away_from_myeloid():
    """Strong progenitor + moderate broad context + weak committed anchor."""
    a = _cluster_adata({"x": {**PROGENITOR, "FCER1G": 0.39, "TYROBP": 0.46,
                              "LST1": 0.64, "AIF1": 0.77, "CTSS": 0.61,
                              "FCN1": 0.26, "FCGR3A": 0.13, "CSF1R": 0.10}})
    tbl = ann.adjudicate_clusters(a, "clust")
    assert tbl.loc["x", "winner"] in ("HSPC", AMB)
    assert tbl.loc["x", "winner"] != "Myeloid"


def test_lilrb1_is_not_a_committed_monocyte_anchor():
    """ILT2 is on B cells, NK cells and T subsets — dropped from the anchor set."""
    assert "LILRB1" not in config.MYELOID_MONO_ANCHORS


def test_checkpoint2c_changes_no_threshold_and_no_frozen_component():
    assert config.MANUAL_POSITIVE_MIN == 0.5
    assert config.MANUAL_MARKER_DETECT_MIN == 0.25
    assert config.MANUAL_DECISION_MARGIN == 0.15
    assert config.MARKER_COVERAGE_MIN == 0.30
    assert config.CONTRADICTION_MIN_GENES == 2
    assert config.CONTRADICTION_MAX_RATE == 0.25
    assert config.HSPC_CORE == ("CD34", "HLF", "SPINK2", "GATA2", "MEIS1", "SOX4")
    assert config.PLASMA_SECRETORY == ("MZB1", "XBP1")
    assert config.PLASMA_MATURE == ("SDC1", "TNFRSF17")
    assert config.LINEAGE_PROGRAMS["myeloid"] == ("CD14", "FCN1", "MNDA", "ITGAM")


# --------------------------------------------------------------------------
# Checkpoint 2d — conventional DC: cDC-specific anchor + APC/broad context
# --------------------------------------------------------------------------

BROAD_INNATE = {"FCER1G": HI, "TYROBP": HI, "LST1": HI, "AIF1": HI, "CTSS": HI}
CDC_ANCHOR = {"FCER1A": HI, "CD1C": HI, "CLEC10A": HI, "CD1E": HI}


def test_cross_route_anchor_context_consistency():
    """No gene may be a lineage-restricted anchor in one Myeloid route while serving
    as generic context in another. This is the C2d defect, checked exhaustively."""
    roles: dict[str, set[str]] = {}
    for genes, role in ((config.MYELOID_MONO_ANCHORS, "ANCHOR"),
                        (config.MYELOID_MONO_CONTEXT, "context"),
                        (config.MYELOID_DC_ANCHORS, "ANCHOR"),
                        (config.MYELOID_DC_CONTEXT, "context"),
                        (config.MYELOID_PDC_CORE, "ANCHOR"),
                        (config.PDC_CONTEXTUAL, "context")):
        for g in genes:
            roles.setdefault(g, set()).add(role)
    conflicts = {g for g, r in roles.items() if len(r) > 1}
    assert not conflicts, f"anchor/context contradiction across routes: {sorted(conflicts)}"


def test_cdc_anchor_set_needs_two_independent_anchors():
    """Four genes, so the 0.5 axis cannot be satisfied by one accidental marker."""
    assert len(config.MYELOID_DC_ANCHORS) >= 4
    one = _support({**MHCII, **BROAD_INNATE, "CD1C": HI})
    assert one["Myeloid"]["subprograms"]["dc"] < config.MANUAL_POSITIVE_MIN


def test_broad_innate_context_cannot_establish_cdc():
    ev = _support({**MHCII, **BROAD_INNATE})
    assert not ev["Myeloid"]["supported"]
    assert ev["Myeloid"]["subprograms"]["dc"] == 0.0


def test_mhcII_plus_broad_innate_without_cdc_anchor_is_not_cdc():
    ev = _support({**MHCII, **BROAD_INNATE, "LYZ": HI})
    assert not ev["Myeloid"]["supported"]


def test_progenitor_with_broad_innate_and_apc_context_is_not_myeloid_via_dc():
    """Leiden-24 shape under C2d: real progenitor evidence, broad innate + APC
    context, no cDC-specific anchor."""
    ev = _support({**PROGENITOR, **MHCII, **BROAD_INNATE})
    assert not ev["Myeloid"]["supported"]
    assert ev["HSPC"]["supported"]


def test_b_cell_with_hlaII_and_cd74_without_cdc_anchor_is_not_myeloid():
    ev = _support({**MHCII, "MS4A1": HI, "CD79A": HI, "CD19": HI})
    assert not ev["Myeloid"]["supported"]
    assert ev["Bcell"]["supported"]


def test_coherent_cdc_anchors_plus_apc_context_is_myeloid():
    ev = _support({**CDC_ANCHOR, **MHCII, **BROAD_INNATE})
    assert ev["Myeloid"]["supported"]
    assert "dc" in ev["Myeloid"]["subprograms_passing"]


def test_monocyte_route_remains_independent_of_cdc_anchors():
    ev = _support({**MONO_ANCHOR, **MONO_CONTEXT})
    assert ev["Myeloid"]["supported"]
    assert "monocyte" in ev["Myeloid"]["subprograms_passing"]
    assert "dc" not in ev["Myeloid"]["subprograms_passing"]


def test_pdc_route_remains_independent_of_cdc_anchors():
    ev = _support({"TCF4": HI, "LILRA4": HI, "CLEC4C": HI, "IL3RA": HI, "PLD4": HI})
    assert ev["Myeloid"]["supported"]
    assert ev["Myeloid"]["subprograms_passing"] == ["pdc"]


def test_flt3_is_not_a_cdc_anchor():
    """FLT3 is on hematopoietic progenitors; using it would reintroduce the
    progenitor false positive this revision removes."""
    assert "FLT3" not in config.MYELOID_DC_ANCHORS


def test_c2d_changes_no_threshold_and_not_the_contradiction_panel():
    assert config.MANUAL_POSITIVE_MIN == 0.5
    assert config.MANUAL_MARKER_DETECT_MIN == 0.25
    assert config.MARKER_COVERAGE_MIN == 0.30
    assert config.CONTRADICTION_MIN_GENES == 2
    assert config.CONTRADICTION_MAX_RATE == 0.25
    assert config.LINEAGE_PROGRAMS["myeloid"] == ("CD14", "FCN1", "MNDA", "ITGAM")
    assert config.MYELOID_MONO_ANCHORS == ("CSF1R", "FCN1", "VCAN", "MS4A7",
                                           "SERPINA1", "CD300E", "FCGR3A")
    assert config.MYELOID_PDC_CORE == ("TCF4", "LILRA4", "CLEC4C", "IL3RA", "PLD4")


# --------------------------------------------------------------------------
# T-lineage revision — TRBC1/2 move from identity to context
#
# Operational conclusion the revision rests on: isolated TRBC1/TRBC2 expression is
# insufficient evidence of T-lineage commitment, because it frequently occurs without
# coordinated CD3/TRAC expression in cells carrying strong NK-lineage evidence.
# No numeric threshold changes; `min_genes` reuses CONTRADICTION_MIN_GENES.
# --------------------------------------------------------------------------

def _cyto(counts: dict[str, int], n: int = 5):
    import anndata as ad
    genes = sorted(set(config.T_IDENTITY_ANCHORS) | set(config.T_CONTEXT)
                   | set(config.NK_IDENTITY) | set(config.GD_IDENTITY)
                   | set(config.CYTOTOXIC_STATE))
    X = np.tile([[float(counts.get(g, 0)) for g in genes]], (n, 1))
    a = ad.AnnData(X)
    a.var_names = genes
    a.obs_names = [f"c{i}" for i in range(n)]
    a.layers["counts"] = a.X.copy()
    return a


NK_ANCH = {"KLRD1": 4, "KLRF1": 3, "FCGR3A": 2, "KLRC1": 2}
CYT_HI = {"NKG7": 9, "GNLY": 8, "PRF1": 5, "GZMB": 6}


def test_A_genuine_t_cell_is_t_lineage():
    r = ann.cytotoxic_lineage_calls(_cyto({"CD3D": 3, "CD3E": 4, "TRAC": 5,
                                           "TRBC1": 6, "TRBC2": 4}))
    assert set(r["call"]) == {"T_ab"}
    assert r["T_identity_pos"].all()


def test_B_nk_with_trbc_context_is_not_t_lineage():
    """THE case the revision exists for: strong NK identity, TRBC present, no CD3/TRAC."""
    r = ann.cytotoxic_lineage_calls(_cyto({**NK_ANCH, "TRBC1": 4, "TRBC2": 3}))
    assert set(r["call"]) == {"NK"}
    assert not r["T_identity_pos"].any()
    assert r["T_context_pos"].all()          # context still measured and reported


def test_C_cytotoxic_t_cell_keeps_t_identity():
    """Shared cytotoxic genes must not cost a real T cell its lineage."""
    r = ann.cytotoxic_lineage_calls(_cyto({"CD3D": 4, "CD3E": 3, "TRAC": 4, **CYT_HI}))
    assert set(r["call"]) == {"T_ab"}


def test_D_trbc_only_cell_is_not_t_lineage():
    r = ann.cytotoxic_lineage_calls(_cyto({"TRBC1": 5, "TRBC2": 4}))
    assert set(r["call"]) == {"unresolved"}
    assert not r["T_identity_pos"].any()


def test_E_genuine_t_and_nk_evidence_stays_mixed():
    """Coordinated T anchors AND NK anchors: preserved as mixed, no forced winner."""
    r = ann.cytotoxic_lineage_calls(_cyto({"CD3D": 4, "CD3E": 3, "TRAC": 4, **NK_ANCH}))
    assert set(r["call"]) == {"T_NK_mixed"}


def test_F_trbc_output_is_preserved_not_removed():
    r = ann.cytotoxic_lineage_calls(_cyto({**NK_ANCH, "TRBC1": 7, "TRBC2": 2}))
    assert "n_T_context" in r.columns and "T_context_pos" in r.columns
    assert int(r["n_T_context"].iloc[0]) == 2      # both TRBC genes still counted
    assert config.T_CONTEXT == ("TRBC1", "TRBC2")
    assert set(config.T_CONTEXT) <= set(config.MARKER_PANEL["Tcell"])   # still in the panel


def test_G_cytotoxic_state_alone_establishes_no_lineage():
    r = ann.cytotoxic_lineage_calls(_cyto(CYT_HI))
    assert set(r["call"]) == {"unresolved"}
    assert not r["T_identity_pos"].any() and not r["NK_identity_pos"].any()


def test_trdc_alone_is_not_gamma_delta_identity():
    r = ann.cytotoxic_lineage_calls(_cyto({**NK_ANCH, "TRDC": 6}))
    assert set(r["call"]) == {"NK"}
    assert not r["gd_pos"].any()


def test_gamma_delta_needs_coordinated_trdc_and_trgc():
    r = ann.cytotoxic_lineage_calls(_cyto({"CD3D": 3, "CD3E": 3, "TRDC": 4, "TRGC1": 3}))
    assert set(r["call"]) == {"T_gd"}


def test_revision_introduces_no_new_numeric_parameter():
    a = _cyto({"CD3D": 3, "CD3E": 3})
    assert (ann.cytotoxic_lineage_calls(a)["call"]
            == ann.cytotoxic_lineage_calls(a, min_genes=config.CONTRADICTION_MIN_GENES)["call"]).all()


def test_global_panels_are_untouched_by_this_revision():
    """The revision is validated on cluster 23 only; global C2d behaviour must not move."""
    assert config.MARKER_PANEL["Tcell"] == ("CD3D","CD3E","CD3G","TRAC","TRBC1","TRBC2")
    assert config.LINEAGE_PROGRAMS["T"] == ("CD3D","CD3E","CD3G","TRAC","TRBC1","TRBC2")
    assert config.MANUAL_POSITIVE_MIN == 0.5
    assert config.CONTRADICTION_MIN_GENES == 2
    assert config.CONTRADICTION_MAX_RATE == 0.25
    assert config.MARKER_COVERAGE_MIN == 0.30
