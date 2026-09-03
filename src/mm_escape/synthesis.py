"""Stage 12 — deterministic synthesis text and derivations.

Stage 12 is a **synthesis stage**: every quantity it reports already exists in a frozen
artifact. This module holds only the small amount of logic that turns frozen states into
the predeclared text fields and derived flags of
``results/12_final_synthesis/stage12_patient_evidence_matrix.csv``.

Three rules bind everything here, and each is enforced by ``tests/test_stage12_synthesis.py``:

* **No scoring.** Nothing in this module may sum, weight, average or rank an evidence axis.
  ``n_uncertainty_flags`` is a *count of flags*, never an evidence score.
* **No freeform prose.** Every text field is a pure function of frozen inputs, selected from
  a fixed template set declared in ``docs/stage12_design.md``. Running twice must give
  byte-identical output.
* **``NOT_EVALUABLE`` is never negative.** It carries no evidential weight in either
  direction and may not be coerced to ``False``, ``0`` or "not supported".
"""
from __future__ import annotations

__all__ = [
    "UNCERTAINTY_PRIORITY", "FLAG_COLUMNS", "STRENGTH_VOCABULARY",
    "STAGE10_DE_METHOD", "PHENOTYPE_COMPATIBILITY_NOTE", "CNV_NOT_EVALUABLE_REASON",
    "IMMUNE_CONTEXT_SUMMARY", "ANCHOR_VS_ALTERNATIVE_NOTE", "COVERAGE_QC_NOTE",
    "main_uncertainty", "measurement_interpretation", "biological_interpretation",
    "allowed_claim", "prohibited_claim", "evidence_profile", "direction_sign",
]

#: Display priority for the ``main_uncertainty`` text field. Declared in the design before
#: execution and never reordered afterwards. This selects *which caveat to name first* on a
#: patient's row. It is NOT an ordering of patients and NOT a weighting of evidence.
UNCERTAINTY_PRIORITY = (
    "not_evaluable", "low_n", "dropout_compatible", "denominator",
    "repeated_sample", "depth", "null_scheme", "intermediate_dn", "none",
)

#: Frozen Stage-08/09b flag columns, in their upstream spelling.
FLAG_COLUMNS = (
    "uncertain_low_n", "uncertain_denominator", "uncertain_depth",
    "uncertain_repeated_sample", "uncertain_null_scheme",
    "uncertain_dropout_compatible", "uncertain_intermediate_dn", "uncertain_threshold",
)

#: The only evidence-strength values Stage 12 may emit. Semantics in the design.
STRENGTH_VOCABULARY = (
    "STRONG", "SUPPORTED_WITH_CAVEATS", "EXPLORATORY", "NOT_SUPPORTED", "NOT_EVALUABLE",
)

#: The recovered Stage-10 producer never imports pydeseq2. This wording is mandatory in
#: every Stage-12 output; the historical ``~ patient + group`` design string in
#: ``pseudobulk_de_evaluability.csv`` describes intent, not a fitted model.
STAGE10_DE_METHOD = (
    "Stage-10 differential expression used depth-matched patient pseudobulks. For each "
    "gene, paired patient-level DN-versus-comparator log-fold changes were tested with a "
    "two-sided Wilcoxon signed-rank test, followed by Benjamini-Hochberg correction. "
    "Patient is the biological unit."
)

#: Carried on every row whose Level-2 state is supported.
PHENOTYPE_COMPATIBILITY_NOTE = (
    "weakly discriminative: 26 of 27 evaluable patients satisfy the frozen per-patient "
    "Level-2 rule; indicates compatibility with the cohort-level DN-associated programme, "
    "not patient-specific evidence of a distinct escape state"
)

CNV_NOT_EVALUABLE_REASON = (
    "infercnvpy donor negative control failed (donor plasma-cell false-positive rate "
    "0.0-50.6% at z>3); the method was rejected before any disease CNV was inspected, so "
    "it contributes no evidence in either direction"
)

IMMUNE_CONTEXT_SUMMARY = (
    "no robust independent immune association (cohort-level): 0 of 28 composition tests at "
    "BH<0.10; the targeted ligand-receptor panel is receiver-state confounded; the single "
    "LIANA consensus hit is antigen-circular"
)

ANCHOR_VS_ALTERNATIVE_NOTE = (
    "alternative eligible combinations show a lower observed transcript-level uncovered "
    "fraction than the BCMA+GPRC5D anchor, but this reflects differences in transcript "
    "detection rate (1.8-2.8x) rather than therapeutic superiority; no combination is "
    "optimal, recommended or best"
)

COVERAGE_QC_NOTE = (
    "no target in the evaluated panel is depth-robust (detection-vs-depth rho 0.32-0.48); "
    "GPRC5D is COVERAGE_NOT_EVALUABLE on technical-zero burden (0.62 vs a 0.50 gate fixed "
    "beforehand) while remaining the frozen Stage-08 anchor; SDC1 is COVERAGE_NOT_EVALUABLE "
    "on circularity; TNFRSF17 carries the identical selection-dependence as a disclosure"
)

_NOT_EVALUABLE = "NOT_EVALUABLE"


def _is_not_evaluable(state: str) -> bool:
    return str(state).endswith(_NOT_EVALUABLE)


def direction_sign(value) -> str:
    """Sign of a per-patient matched programme delta, as a descriptive label.

    These are **descriptive signs, not tests**. The frozen design tests programmes at
    cohort level; per-patient deltas feed only the weakly-discriminative Level-2 rule, so
    no p-value is ever attached to a patient for a programme direction.
    """
    if value is None:
        return "not_evaluable"
    try:
        v = float(value)
    except (TypeError, ValueError):
        return "not_evaluable"
    if v != v:                       # NaN
        return "not_evaluable"
    if v < 0:
        return "lower_in_DN"
    if v > 0:
        return "higher_in_DN"
    return "no_difference"


def main_uncertainty(level1_state, level2_state, flags: dict) -> str:
    """Name the highest-priority caveat for this patient, by fixed display priority.

    Non-evaluability outranks every measurement flag: a patient whose Level-1/2 evidence
    could not be assessed has an evaluability caveat, not merely a noisy estimate.
    """
    if _is_not_evaluable(level1_state) or _is_not_evaluable(level2_state):
        return "not_evaluable"
    lookup = {
        "low_n": "uncertain_low_n",
        "dropout_compatible": "uncertain_dropout_compatible",
        "denominator": "uncertain_denominator",
        "repeated_sample": "uncertain_repeated_sample",
        "depth": "uncertain_depth",
        "null_scheme": "uncertain_null_scheme",
        "intermediate_dn": "uncertain_intermediate_dn",
    }
    for name in UNCERTAINTY_PRIORITY:
        if name in ("not_evaluable", "none"):
            continue
        if bool(flags.get(lookup[name], False)):
            return name
    return "none"


def measurement_interpretation(tier: str, n_flags: int, low_n: bool) -> str:
    """Template selected by (tier, flag load). Measurement statement only."""
    if tier == "robust-high":
        return ("observed DN estimate survived every frozen sensitivity analysis "
                "(denominator, depth, truncate-10k, repeated-sample, null-scheme, "
                "threshold); measurement-robust only, not a biological classification")
    if tier == "robust-low":
        return ("observed DN estimate consistently low across every frozen sensitivity "
                "assumption; measurement statement only")
    if low_n:
        return ("observed DN estimate is measurement-uncertain, primarily because the "
                "malignant-cell denominator is small")
    if n_flags >= 3:
        return (f"observed DN estimate is measurement-uncertain under {n_flags} frozen "
                "sensitivity flags")
    return ("observed DN estimate is measurement-uncertain under "
            f"{n_flags} frozen sensitivity flag(s)")


def biological_interpretation(level1_state: str, level2_state: str) -> str:
    """Template selected by (Level-1, Level-2). Level-3 is always appended."""
    genomic = "genomic subclone evidence not evaluable"
    if _is_not_evaluable(level1_state) and _is_not_evaluable(level2_state):
        return ("DN structure and DN-associated phenotype were both not evaluable for this "
                f"patient; this is an absence of evidence, not negative evidence; {genomic}")
    l1_sup = level1_state.endswith("SUPPORTED") and not _is_not_evaluable(level1_state)
    l2_sup = level2_state.endswith("SUPPORTED") and not _is_not_evaluable(level2_state)
    if l1_sup and l2_sup:
        body = ("DN cells are non-randomly organised beyond the depth-stratified null and "
                "are compatible with the cohort-level DN-associated programme; "
                "compatibility is weakly discriminative (26 of 27 evaluable patients)")
    elif l1_sup:
        body = ("DN cells are non-randomly organised beyond the depth-stratified null; no "
                "DN-associated phenotype support")
    elif l2_sup:
        body = ("no support for non-random DN organisation beyond the depth-stratified "
                "null; DN cells are compatible with the cohort-level DN-associated "
                "programme, which is weakly discriminative (26 of 27 evaluable patients)")
    else:
        body = ("no support for non-random DN organisation and no support for the "
                "cohort-level DN-associated programme")
    return f"{body}; {genomic}"


def allowed_claim(licensed_language: str) -> str:
    """Wrap the frozen ``subclone.licensed_language`` output as permitted wording."""
    mapping = {
        "escape-associated subclone": "escape-associated subclone",
        "escape-associated transcriptional state": (
            "escape-associated transcriptional state (cohort-level phenotype; "
            "patient-level compatibility is weakly discriminative)"),
        "non-random DN organization": "non-random DN organization",
        "none": ("observed transcript-level double-negative fraction only"),
    }
    return mapping.get(str(licensed_language), "observed transcript-level "
                       "double-negative fraction only")


def prohibited_claim(level1_state: str, level2_state: str) -> str:
    """Constant per (Level-1, Level-2), always naming the two universal prohibitions."""
    parts = ["no genomic subclone claim", "no clinical recommendation"]
    if not (level1_state.endswith("SUPPORTED") and not _is_not_evaluable(level1_state)):
        parts.append("no claim of non-random DN organization")
    if _is_not_evaluable(level1_state) or _is_not_evaluable(level2_state):
        parts.append("not-evaluable must not be read as negative evidence")
    parts.append("no antigen-specific escape mechanism")
    return "; ".join(parts)


def evidence_profile(tier: str, level1_state: str, level2_state: str,
                     level3_state: str) -> str:
    """Lossless concatenation of four frozen states. NOT a category and NOT ordinal.

    Takes exactly the five values enumerated in the design's concordance section. Its only
    purpose is grouping rows; it introduces no new vocabulary and cannot be ordered.
    """
    return (f"tier={tier}|L1={level1_state}|L2={level2_state}|L3={level3_state}")
