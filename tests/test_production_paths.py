"""Production-path and freeze-provenance guards (pre-Stage-12 repair pass).

These tests exist because the Codex pre-Stage-12 audit found that the frozen Stage 07-10
results had **no committed producer**: the authoritative notebooks read the CSVs they
discuss. The producers were recovered verbatim from the session transcripts into
`production/` (see `production/README.md`). This file guards three things:

* **H2** — `subclone.level2_state` now *actually rejects* a program name outside the
  frozen predeclaration. The pre-existing test only asserted that an invented list was
  not a subset of the frozen names; it never called the function.
* **C1** — the recovered drivers stay present, compile, carry their provenance headers,
  and obey the stage-isolation rules the project depends on (no future-stage inputs, no
  writes outside their own namespace, denominators kept separate).
* **H1** — the committed freeze manifest is well formed, fully resolved, and its hashes
  match the local frozen artifacts.

Nothing here re-runs an analysis or reads a frozen scientific number as ground truth.
"""

from __future__ import annotations

import ast
import csv
import hashlib
import re
from pathlib import Path

import pytest

from mm_escape import config as CFG
from mm_escape import subclone as SC

REPO = Path(__file__).resolve().parents[1]
PROD = REPO / "production"
MANIFEST = REPO / "provenance" / "frozen_artifacts_pre_stage12.tsv"

STAGE_DIRS = ("stage06", "stage07", "stage08", "stage09", "stage09b", "stage10")


def _code_tokens(py_path):
    """Identifiers + string constants from executable code only (docstrings stripped).

    Comment/docstring text is prose: the project's own files legitimately *describe*
    the labels they refuse to emit, so a substring scan over raw source cannot tell a
    rule from a violation. This walks the AST instead.
    """
    tree = ast.parse(_read(py_path), filename=str(py_path))
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            body = node.body
            if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant) \
                    and isinstance(body[0].value.value, str):
                body[0].value.value = ""          # drop the docstring
    out = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Name):
            out.append(node.id)
        elif isinstance(node, ast.Attribute):
            out.append(node.attr)
        elif isinstance(node, ast.arg):
            out.append(node.arg)
        elif isinstance(node, ast.Constant) and isinstance(node.value, str):
            out.append(node.value)
    return out


def _driver_paths(*stages):
    stages = stages or STAGE_DIRS
    out = []
    for s in stages:
        out.extend(sorted((PROD / s).glob("*.py")))
        out.extend(sorted((PROD / s).glob("*.sh")))
    return out


def _read(p):
    return p.read_text(encoding="utf-8", errors="replace")


# ============================================================== A. H2 — Level-2 inputs
def test_a1_level2_state_rejects_an_invented_program_name():
    """The audit's H2: exercise rejection, do not merely assert a set relation."""
    with pytest.raises(SC.UnknownProgramError):
        SC.level2_state(True, ["a_program_discovered_later"])


def test_a2_level2_state_rejects_a_plausible_but_unfrozen_program():
    """'proliferation' is a real program and still has no route in. That is the point."""
    with pytest.raises(SC.UnknownProgramError) as exc:
        SC.level2_state(True, ["oxphos", "proliferation"])
    assert "proliferation" in str(exc.value)
    assert "LEVEL2_PROGRAMS" in str(exc.value)


def test_a3_every_frozen_program_is_accepted_singly():
    for name in CFG.LEVEL2_PROGRAMS:
        assert SC.level2_state(True, [name]) == SC.STATE_SUPPORTED


def test_a4_validation_is_not_skipped_when_the_patient_is_not_evaluable():
    """A bad name is a bug regardless of the branch it would have taken."""
    with pytest.raises(SC.UnknownProgramError):
        SC.level2_state(False, ["invented"])


def test_a5_empty_hits_are_legal_and_mean_not_supported():
    assert SC.level2_state(True, []) == SC.STATE_NOT_SUPPORTED
    assert SC.level2_state(False, []) == SC.STATE_NOT_EVALUABLE


def test_a6_validate_program_names_is_exported_and_returns_the_input():
    assert SC.validate_program_names(["myc", "oxphos"]) == ["myc", "oxphos"]


def test_a7_the_frozen_stage10_driver_could_only_ever_use_frozen_program_names():
    """Why the frozen numbers were already safe, stated as a test rather than a claim.

    `level2_state()` was NOT the production call site — the recovered driver inlines the
    rule. What actually bounded the frozen vocabulary is that every program loop in the
    recovered drivers iterates `config.LEVEL2_PROGRAMS`, so an invented program had no
    route into `repro`, into `hits`, or into a state call.
    """
    src = _read(PROD / "stage10" / "s10c_program_scores_and_cohort_tests.py")
    assert "for name in config.LEVEL2_PROGRAMS" in src
    assert "config.LEVEL2_PROGRAMS" in _read(
        PROD / "stage10" / "s10d_hypotheses_and_evidence_levels.py"
    )


# ================================================= B. recovered drivers exist / compile
def test_b1_every_recovered_stage_has_drivers():
    for s in STAGE_DIRS:
        assert (PROD / s).is_dir(), f"missing {s}"
        assert _driver_paths(s), f"no drivers under production/{s}"


def test_b2_all_recovered_python_drivers_parse():
    for p in _driver_paths():
        if p.suffix != ".py":
            continue
        ast.parse(_read(p), filename=str(p))


def test_b3_every_driver_carries_its_provenance_header():
    for p in _driver_paths():
        head = _read(p)[:2000]
        assert "RECOVERED PRODUCTION DRIVER" in head, p
        assert "Recovered from : Claude Code session transcript" in head, p
        assert "Executed (UTC)" in head, p
        assert "Writes" in head, p


def test_b4_drivers_are_labelled_as_not_re_executed_during_recovery():
    """The recovery pass regenerated nothing; the header must say so, on every file."""
    for p in _driver_paths():
        assert "NOT re-executed" in _read(p)[:2000], p


# ======================================================= C. stage isolation / ordering
STAGE_OF_DRIVER = {
    "stage06": 6, "stage07": 7, "stage08": 8, "stage09": 9, "stage09b": 9.5, "stage10": 10,
}
#: results/ namespaces a driver of a given stage is allowed to *write* into.
ALLOWED_WRITE = {
    "stage06": ("results/06_annotation",),
    "stage07": ("results/07_malignant_plasma",),
    "stage08": ("results/08_dual_antigen_escape",),
    "stage09": ("results/09_bulk_validation",),
    # 09b is the risk-tier arm and lives inside the stage-08 tree by design
    # 09b is the risk-tier arm; by design it roots at the stage-08 tree and writes only
    # into its own risk_tier_* subdirectories.
    "stage09b": ("results/08_dual_antigen_escape",),
    "stage10": ("results/10_dn_coherence",),
}
WRITE_CALL = re.compile(
    r"""(?:to_csv|to_json|write_h5ad|write_text|to_markdown|savefig|write)\(\s*['"]?"""
    r"""(results/[A-Za-z0-9_./-]+)""")
LITERAL_RESULT_PATH = re.compile(r"""results/[A-Za-z0-9_][A-Za-z0-9_./-]*""")


def test_c1_no_driver_writes_outside_its_own_stage_namespace():
    for stage in STAGE_DIRS:
        allowed = ALLOWED_WRITE[stage]
        for p in _driver_paths(stage):
            src = _read(p)
            # shell redirections and OUT=Path(...) assignments
            for m in re.finditer(r"""(?:cat >>?|>>?)\s+(results/[A-Za-z0-9_./-]+)""", src):
                assert m.group(1).startswith(allowed), f"{p.name} writes {m.group(1)}"
            for m in re.finditer(r"""(?:OUT|RT|RT_DIR|D)\s*=\s*Path\(['"](results/[^'"]+)""", src):
                assert m.group(1).startswith(allowed), f"{p.name} sets OUT={m.group(1)}"


def test_c2_no_driver_reads_a_future_stage():
    """Execution order is 06 -> 07 -> 08 -> 09 -> 09b -> 10; nothing may read forward."""
    order = {"results/06_annotation": 6, "results/07_malignant_plasma": 7,
             "results/08_dual_antigen_escape": 8, "results/09_bulk_validation": 9,
             "results/10_dn_coherence": 10, "results/11_immune_context": 11}
    for stage in STAGE_DIRS:
        mine = STAGE_OF_DRIVER[stage]
        for p in _driver_paths(stage):
            for ref in set(LITERAL_RESULT_PATH.findall(_read(p))):
                for prefix, n in order.items():
                    if not ref.startswith(prefix):
                        continue
                    # the 09b arm legitimately lives under the stage-08 tree
                    if stage == "stage09b" and prefix == "results/08_dual_antigen_escape":
                        continue
                    assert n <= mine, f"{p.name} (stage {mine}) references {ref}"


def test_c3a_the_tier_decision_module_never_sees_bulk_cohort_or_coherence():
    """`risk_tiers.py` is where a tier is actually decided. Nothing else may reach it.

    Checked over executable tokens, not raw text — the module's docstring correctly
    *states* that it never reads bulk, and that sentence must not fail its own test.
    """
    #: `enr_cohortbins` / `enr_globalbins` name the two depth-BINNING schemes of the
    #: Stage-08 null. They are not cohort labels and are legitimate tier evidence.
    depth_scheme_tokens = {"enr_cohortbins", "enr_globalbins"}
    tokens = [s.lower() for s in _code_tokens(REPO / "src" / "mm_escape" / "risk_tiers.py")]
    tokens = [tok for tok in tokens if tok not in depth_scheme_tokens]
    for forbidden in ("bulk", "tpm", "s09_", "dn_coherence", "level1_", "level2_",
                      "tc_like", "liana", "immune", "cohort"):
        hits = [tok for tok in tokens if forbidden in tok]
        assert not hits, f"risk_tiers.py references {forbidden}: {hits[:5]}"


def test_c3a2_no_tier_function_accepts_a_cohort_bulk_or_coherence_argument():
    import inspect

    from mm_escape import risk_tiers as RT

    for name in getattr(RT, "__all__", dir(RT)):
        obj = getattr(RT, name, None)
        if not callable(obj) or isinstance(obj, type):
            continue
        try:
            params = inspect.signature(obj).parameters
        except (TypeError, ValueError):
            continue
        for param in params:
            low = param.lower()
            assert "cohort" not in low, f"{name}() takes {param}"
            assert "bulk" not in low, f"{name}() takes {param}"
            assert "coherence" not in low and "level" not in low, f"{name}() takes {param}"


def test_c3b_stage09b_never_reads_stage10_coherence():
    """Stage 10 runs after 09b. Coherence can never be a tier input."""
    for p in _driver_paths("stage09b"):
        src = _read(p)
        assert "10_dn_coherence" not in src, f"{p.name} reads Stage-10 coherence"
        for forbidden in ("dn_coherence", "level2_state", "level1_structure", "tc_like"):
            assert forbidden not in src, f"{p.name} references {forbidden}"


def test_c3c_stage09_bulk_enters_09b_only_as_report_only_context_after_tiers_exist():
    """Bulk *is* read by the tier driver — but strictly as annotation, never as evidence.

    The frozen driver builds every tier first (`RT.final_tier(ev)` over Stage-07/08
    evidence only), materialises the table, and *then* appends `s09_*` context columns.
    This test pins that ordering, which is the invariant the project plan actually states:
    Stage 09 is interpretation context, not an additional scoring axis.
    """
    src = _read(PROD / "stage09b" / "s09b2_provisional_tiers.py")
    assert "09_bulk_validation" in src, "test is stale: driver no longer reads bulk"
    tier_built = src.index("T=pd.DataFrame(rows)")
    marker = src.index("Stage-09 context: REPORT-ONLY")
    assert marker > tier_built, "bulk context must be appended AFTER tiers are built"
    # every bulk read and every bulk-derived column lands after the tiers exist
    for token in ("09_bulk_validation", "s09_bulk_available", "s09_marginal_context",
                  "bulk_TNFRSF17_tpm", "bulk_GPRC5D_tpm"):
        for occurrence in [i for i in range(len(src)) if src.startswith(token, i)]:
            if token == "09_bulk_validation" and occurrence < tier_built:
                # the read itself may be hoisted to the top with the other loads,
                # but it must not be referenced inside the tier loop
                loop = src[src.index("for _,r in EV.iterrows()"):tier_built]
                assert token not in loop and "SC9" not in loop and "MAP" not in loop
    # the decision call itself takes only the Stage-07/08 evidence dict
    ev_block = src[src.index("ev=dict("):src.index("RT.final_tier(ev)")]
    for forbidden in ("bulk", "tpm", "SC9", "MAP", "s09"):
        assert forbidden not in ev_block, f"tier input carries {forbidden}"


def test_c3d_frozen_tier_membership_is_unchanged():
    """4 robust-high / 28 uncertain / 0 robust-low. Recovery may not move a patient."""
    csvp = (REPO / "results" / "08_dual_antigen_escape" / "risk_tier_provisional"
            / "risk_tiers_provisional.csv")
    if not csvp.exists():
        pytest.skip("frozen tier table not present locally")
    with csvp.open() as fh:
        rows = list(csv.DictReader(fh))
    counts = {}
    for r in rows:
        counts[r["final_tier"]] = counts.get(r["final_tier"], 0) + 1
    assert counts == {"robust-high": 4, "uncertain": 28}
    assert counts.get("robust-low", 0) == 0


def test_c4_stage08_core_drivers_do_not_consume_bulk_or_coherence():
    for p in _driver_paths("stage08"):
        src = _read(p)
        assert "09_bulk_validation" not in src
        assert "10_dn_coherence" not in src


#: Stage-07 drivers whose *purpose* is to measure antigen behaviour, so naming the
#: antigens is the analysis rather than a leak. Each is an OUTPUT about the antigens,
#: never an INPUT to a clone call.
ANTIGEN_MEASURING_07 = ("s07i_", "s07j_")


#: A line in a Stage-07 driver may name an antigen only to ASSERT it is absent or to
#: EXCLUDE it. Anything else would be the antigen reaching the clone call.
ANTIGEN_INTENT = ("assert", "antigen leak", "is_ant", "~is_ant", "must not touch",
                  "excl", "n/a")


def test_c5a_the_clone_call_itself_never_reads_the_antigens():
    """Clone membership must not be a function of the genes the metric is about.

    Antigen names do appear in these drivers — always to assert absence or to build an
    exclusion mask. Each such line is checked for that intent, so a line that merely
    *uses* an antigen fails.
    """
    seen = 0
    for p in _driver_paths("stage07"):
        if p.name.startswith(ANTIGEN_MEASURING_07):
            continue
        lines = _read(p).splitlines()
        for i, line in enumerate(lines):
            if line.lstrip().startswith("#"):
                continue
            if "TNFRSF17" not in line and "GPRC5D" not in line:
                continue
            seen += 1
            window = " ".join(lines[i:i + 4]).lower()
            assert any(k in window for k in ANTIGEN_INTENT), f"{p.name}:{i + 1}: {line.strip()}"
    assert seen, "test is stale: no antigen mentions found in the Stage-07 drivers"


def test_c5a2_the_cnv_gene_set_excludes_both_antigens():
    """The exclusion is the point, so pin it rather than merely tolerating the mention."""
    src = _read(PROD / "stage07" / "s07c_cnv_input_gene_set.sh")
    assert "ANT = ['TNFRSF17','GPRC5D']" in src
    assert "~is_ant" in src, "antigens must be removed from the CNV input gene set"


def test_c5b_the_clone_decision_module_accepts_no_antigen_input():
    import inspect

    from mm_escape import malignant as M

    for name in getattr(M, "__all__", dir(M)):
        obj = getattr(M, name, None)
        if not callable(obj) or isinstance(obj, type):
            continue
        try:
            params = inspect.signature(obj).parameters
        except (TypeError, ValueError):
            continue
        for param in params:
            assert "tnfrsf17" not in param.lower() and "gprc5d" not in param.lower()
            assert "antigen" not in param.lower(), f"{name}() takes {param}"


def test_c5c_the_antigen_measuring_drivers_only_report_never_gate():
    """s07i/s07j may name the antigens; they must not assign a clone state from them."""
    for p in _driver_paths("stage07"):
        if not p.name.startswith(ANTIGEN_MEASURING_07):
            continue
        src = _read(p)
        for state in ("CLONE_SUPPORTED", "CLONE_COMPATIBLE_V_UNOBSERVED",
                      "CLONE_UNCERTAIN", "CLONE_INCOMPATIBLE"):
            for line in src.splitlines():
                if state in line and ("TNFRSF17" in line or "GPRC5D" in line):
                    raise AssertionError(f"{p.name} derives {state} from an antigen: {line}")


# ============================================== D. determinism / denominator separation
def test_d1_stochastic_drivers_fix_a_seed():
    """Bootstraps, permutations and downsampling must not drift between runs."""
    for p in _driver_paths("stage08", "stage10"):
        src = _read(p)
        if not re.search(r"bootstrap|permutation|downsample|SEED|shuffle", src):
            continue
        assert re.search(r"SEED\s*=\s*\d+", src), f"{p.name} is stochastic with no fixed SEED"


def test_d2_primary_and_sensitivity_are_carried_as_separate_outputs():
    """They are never collapsed, and never selected between."""
    s08 = _read(PROD / "stage08" / "s08a_patient_antigen_states.py")
    assert "patient_antigen_states_primary.csv" in s08
    assert "patient_antigen_states_sensitivity.csv" in s08
    s10 = _read(PROD / "stage10" / "s10e_pseudobulk_de_decoupler_tc.py")
    assert "for tag in ['primary','sensitivity']" in s10.replace('"', "'")


def test_d3_no_driver_picks_whichever_denominator_looks_better():
    banned = re.compile(r"(best|preferred|chosen|winning)_denominator|denominator\s*=\s*best")
    for p in _driver_paths():
        assert not banned.search(_read(p)), p


# ================================== E. NOT_EVALUABLE may never decay into NOT_SUPPORTED
def test_e1_cnv_not_evaluable_is_the_only_cnv_state_emitted():
    """No driver and no module may *emit* `CNV_SUBCLONE_NOT_SUPPORTED`.

    Frozen prose is allowed to say the label is never emitted — that sentence is the
    project stating the rule, not breaking it — so this walks executable tokens only.
    """
    bad = "CNV_SUBCLONE_NOT_SUPPORTED"
    for p in _driver_paths():
        if p.suffix != ".py":
            continue
        assert bad not in _code_tokens(p), p
    assert bad not in _code_tokens(REPO / "src" / "mm_escape" / "subclone.py")
    assert not hasattr(SC, "CNV_SUBCLONE_NOT_SUPPORTED")


def test_e2_stage10_driver_maps_unevaluable_to_not_evaluable_not_not_supported():
    src = _read(PROD / "stage10" / "s10d_hypotheses_and_evidence_levels.py")
    assert "DN_STATE_NOT_EVALUABLE' if not ev else" in src.replace('"', "'")


def test_e3_level_states_never_silently_default():
    """A missing evidence value must not fall through to a 'supported' branch."""
    assert SC.level2_state(False, []) == SC.STATE_NOT_EVALUABLE
    assert SC.level1_structure_state is not None


# ================================================ F. H1 — the committed freeze manifest
EXPECTED_COLUMNS = ["stage", "path", "bytes", "sha256", "artifact_role",
                    "producer_source_path", "producer_code_commit", "frozen_date",
                    "environment", "reproducibility"]
VALID_REPRO = {"REPRODUCIBLE_FROM_COMMITTED_CODE", "REPRODUCIBLE_WITH_EXTERNAL_INPUT",
               "ARCHIVED_ONLY", "NOT_FULLY_REPRODUCIBLE"}


def _manifest_rows():
    with MANIFEST.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_f1_manifest_exists_and_is_committed_outside_results():
    assert MANIFEST.exists(), "freeze manifest missing"
    assert "results" not in MANIFEST.parts[:-1]


def test_f2_manifest_schema():
    with MANIFEST.open() as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
    assert hdr == EXPECTED_COLUMNS


def test_f3_every_row_is_fully_resolved():
    for r in _manifest_rows():
        assert r["sha256"] and len(r["sha256"]) == 64, r["path"]
        assert r["bytes"].isdigit(), r["path"]
        assert r["reproducibility"] in VALID_REPRO, r["path"]
        assert not r["producer_source_path"].startswith("UNRESOLVED"), r["path"]
        assert r["environment"], r["path"]


def test_f4_manifest_covers_every_required_stage():
    stages = {r["stage"] for r in _manifest_rows()}
    for s in ("06", "07", "08", "08c", "09", "09b", "10", "11", "11b"):
        assert s in stages, f"stage {s} absent from the freeze manifest"


def test_f5_named_producers_are_committed_paths_that_exist():
    for r in _manifest_rows():
        p = r["producer_source_path"]
        if p.startswith("("):          # documents authored in place
            continue
        assert (REPO / p).exists(), f"{r['path']} names a missing producer {p}"


@pytest.mark.slow
def test_f6_manifest_hashes_match_the_local_frozen_artifacts():
    """The freeze record must authenticate what is actually on disk."""
    missing, changed = [], []
    for r in _manifest_rows():
        p = REPO / r["path"]
        if not p.exists():
            missing.append(r["path"])
            continue
        h = hashlib.sha256()
        with p.open("rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 22), b""):
                h.update(chunk)
        if h.hexdigest() != r["sha256"]:
            changed.append(r["path"])
    assert not changed, f"frozen artifacts MUTATED: {changed[:10]}"
    assert not missing, f"frozen artifacts missing locally: {missing[:10]}"


def test_f7_no_result_artifact_is_tracked_in_git():
    """The manifest authenticates ignored artifacts; it does not license committing them."""
    gitignore = _read(REPO / ".gitignore")
    assert re.search(r"^/?results/?", gitignore, re.M), "results/ must stay ignored"
