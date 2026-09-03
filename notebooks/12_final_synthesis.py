# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#   kernelspec:
#     display_name: mm-core
#     language: python
#     name: mm-core
# ---

# %% [markdown]
# # Stage 12 — final synthesis / decision packet
#
# **Binding specification: `docs/stage12_design.md`.** This notebook implements that design
# and nothing else. It is a **synthesis stage**: it runs no new statistical test, fits no
# model, creates no threshold, and produces no score.
#
# Every quantity reported here already exists in a frozen artifact. Stage 12 lays six frozen
# evidence axes side by side, per patient and per cohort, in a form that makes their
# **disagreement** legible, and states precisely which claims that evidence licenses.
#
# **Reads** 29 frozen artifacts, read-only, all verified against
# `provenance/frozen_artifacts_pre_stage12.tsv`.
# **Writes** only into `results/12_final_synthesis/`.
#
# Where planning-era text conflicts with the recovered production code under `production/`,
# **the recovered producer is authoritative** — most importantly, Stage-10 differential
# expression is a paired patient-level Wilcoxon signed-rank test, never pydeseq2.

# %%
from __future__ import annotations

import csv
import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

sys.path.insert(0, str(Path.cwd() / "src"))
from mm_escape import subclone as SC          # noqa: E402
from mm_escape import synthesis as SY         # noqa: E402

REPO = Path.cwd()
OUT = REPO / "results" / "12_final_synthesis"
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

MANIFEST = REPO / "provenance" / "frozen_artifacts_pre_stage12.tsv"
DESIGN = REPO / "docs" / "stage12_design.md"

pd.set_option("display.width", 240)
STR = {"patient": str, "patient_id": str}     # int64 coercion once produced an empty join

# %% [markdown]
# ## Step 1 — Verify provenance before reading any scientific input
#
# Abort on the first mismatch. Never regenerate an upstream file and never update the
# committed manifest to agree with disk.

# %%
def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


MANIFEST_ROWS = {}
with open(MANIFEST) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        MANIFEST_ROWS[row["path"]] = row
print(f"committed manifest rows: {len(MANIFEST_ROWS)}")

#: (relative_path, upstream_stage, synthesis_role) — the 29 inputs named in the design.
INPUTS = [
    ("results/08_dual_antigen_escape/patient_antigen_states_primary.csv", "08",
     "DN fraction, marginal detection, enrichment, bootstrap, perm p (primary)"),
    ("results/08_dual_antigen_escape/patient_antigen_states_sensitivity.csv", "08",
     "same quantities under the sensitivity denominator"),
    ("results/08_dual_antigen_escape/patient_bootstrap_intervals.csv", "08",
     "bootstrap CI bounds, low_n, single_sample"),
    ("results/08_dual_antigen_escape/patient_conegativity_enrichment.csv", "08",
     "conditioned/unconditioned co-negativity enrichment"),
    ("results/08_dual_antigen_escape/patient_evidence_states.csv", "08",
     "per-patient uncertainty flags"),
    ("results/08_dual_antigen_escape/repeated_sample_antigen_consistency.csv", "08",
     "per-sample DN within repeated patients"),
    ("results/08_dual_antigen_escape/truncate10k_sensitivity.csv", "08",
     "WashU 10k censoring sensitivity"),
    ("results/08_dual_antigen_escape/primary_vs_sensitivity_denominator_comparison.csv", "08",
     "denominator sensitivity"),
    ("results/08_dual_antigen_escape/noise_floor_technical_zero.csv", "08",
     "expression-matched technical-zero (dropout) burden"),
    ("results/08_dual_antigen_escape/noise_floor_ambient.csv", "08",
     "ambient floor, the opposite-signed bias"),
    ("results/08_dual_antigen_escape/depth_stratified_null.csv", "08",
     "null-scheme sensitivity (cohort vs global bins)"),
    ("results/08_dual_antigen_escape/risk_tier_provisional/risk_tiers_provisional.csv", "09b",
     "provisional measurement tier and uncertainty flags"),
    ("results/08_dual_antigen_escape/risk_tier_design/patient_evidence_matrix.csv", "09b",
     "tier evidence inputs"),
    ("results/10_dn_coherence/stage10_evidence_levels.csv", "10",
     "Level-1/2/3 states and licensed language"),
    ("results/10_dn_coherence/stage10_evaluability.csv", "10",
     "per-patient L1/L2 statistics, depth ratios, DE counts"),
    ("results/10_dn_coherence/dn_local_structure_by_patient.csv", "10",
     "Level-1 statistics per denominator"),
    ("results/10_dn_coherence/dn_program_scores_by_patient.csv", "10",
     "per-patient programme deltas, raw and depth-matched"),
    ("results/10_dn_coherence/level2_program_cohort_tests.csv", "10",
     "cohort-level programme results (the interpretable Level-2 result)"),
    ("results/10_dn_coherence/pseudobulk_de_results.csv", "10",
     "pseudobulk DE, both denominators"),
    ("results/10_dn_coherence/gamma_secretase_hypothesis.csv", "10",
     "the pre-registered gamma-secretase negative"),
    ("results/10_dn_coherence/repeated_sample_dn_coherence.csv", "10",
     "Level-1 repeated-sample consistency"),
    ("results/10_dn_coherence/tc_like_subtype.csv", "10",
     "descriptive TC-like expression proxy"),
    ("results/09_bulk_validation/bulk_vs_sc_by_cohort.csv", "09",
     "marginal bulk validation, cohort-split"),
    ("results/09_bulk_validation/normal_marrow_antigen_context.csv", "09",
     "normal-marrow expression context"),
    ("results/11_immune_context/immune_vs_dn_measurement.csv", "11",
     "immune composition associations"),
    ("results/11_immune_context/liana_verification/liana_vs_dn_associations.csv", "11b",
     "LIANA associations and antigen-circularity flag"),
    ("results/08_dual_antigen_escape/multi_antigen_coverage/stage12_multi_antigen_interface.csv",
     "08c", "per-patient multi-antigen coverage interface"),
    ("results/08_dual_antigen_escape/multi_antigen_coverage/target_measurement_qc.csv", "08c",
     "target eligibility, technical-zero burden, depth sensitivity"),
    ("results/04_qc/umi_censoring_effect.csv", "04", "bias-direction table"),
]
assert len(INPUTS) == 29, f"design names 29 inputs, got {len(INPUTS)}"

verified = []
failures = []
for rel, stage, role in INPUTS:
    p = REPO / rel
    if rel not in MANIFEST_ROWS:
        failures.append((rel, "NOT IN COMMITTED MANIFEST")); continue
    if not p.exists():
        failures.append((rel, "MISSING ON DISK")); continue
    digest, size = sha256_of(p), p.stat().st_size
    exp = MANIFEST_ROWS[rel]
    if digest != exp["sha256"]:
        failures.append((rel, "HASH MISMATCH")); continue
    if str(size) != exp["bytes"]:
        failures.append((rel, "SIZE MISMATCH")); continue
    verified.append(dict(relative_path=rel, upstream_stage=stage, sha256=digest,
                         file_size=size, producer=exp["producer_source_path"],
                         producer_commit=exp["producer_code_commit"],
                         environment=exp["environment"], synthesis_role=role))

if failures:
    for f in failures:
        print("PROVENANCE FAILURE:", f)
    raise SystemExit("STOP: upstream provenance verification failed. Stage 12 not executed.")
print(f"provenance verified: {len(verified)}/29 inputs, hash and size")

# %% [markdown]
# ## Step 2 — Pin execution provenance and snapshot the binding design

# %%
def _git(*args):
    return subprocess.run(["/usr/bin/git", *args], capture_output=True, text=True,
                          cwd=REPO).stdout.strip()


EXEC = dict(
    commit=_git("rev-parse", "HEAD"),
    commit_short=_git("rev-parse", "--short", "HEAD"),
    nearest_tag=_git("describe", "--tags", "--abbrev=0"),
    executed_utc=datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    environment="mm-core",
    python=platform.python_version(),
    pandas=pd.__version__, numpy=np.__version__, matplotlib=matplotlib.__version__,
    design_path=str(DESIGN.relative_to(REPO)), design_sha256=sha256_of(DESIGN),
    design_bytes=DESIGN.stat().st_size,
)
for k, v in EXEC.items():
    print(f"  {k:16s} {v}")

pd.DataFrame(verified).to_csv(OUT / "stage12_input_manifest.tsv", sep="\t", index=False)
print(f"\nwrote {OUT/'stage12_input_manifest.tsv'} ({len(verified)} rows)")

snapshot = [
    "# Stage-12 design snapshot (execution record)", "",
    "Stage 12 was executed against the frozen design below. This snapshot is a verbatim",
    "copy taken at execution time; it is not edited afterwards except for these provenance",
    "fields.", "",
    "| field | value |", "|---|---|",
    f"| binding design | `{EXEC['design_path']}` |",
    f"| design SHA256 | `{EXEC['design_sha256']}` |",
    f"| design bytes | {EXEC['design_bytes']} |",
    f"| execution commit | `{EXEC['commit']}` |",
    f"| nearest tag | `{EXEC['nearest_tag']}` |",
    f"| executed (UTC) | {EXEC['executed_utc']} |",
    f"| environment | {EXEC['environment']} (python {EXEC['python']}) |",
    f"| producer | `notebooks/12_final_synthesis.py` |",
    f"| upstream inputs verified | {len(verified)}/29 (SHA256 + size) |", "",
    "---", "", DESIGN.read_text(),
]
(OUT / "stage12_design_snapshot.md").write_text("\n".join(snapshot))
print(f"wrote {OUT/'stage12_design_snapshot.md'}")

# %% [markdown]
# ## Step 3 — Load frozen inputs (read-only) and assert the frozen invariants
#
# Assert immediately, before doing any work. Abort on any deviation.

# %%
def rd(rel, **kw):
    return pd.read_csv(REPO / rel, **kw)


s08p = rd(INPUTS[0][0], dtype=STR)
s08s = rd(INPUTS[1][0], dtype=STR)
boot = rd(INPUTS[2][0], dtype=STR)
rep08 = rd(INPUTS[5][0], dtype=STR)
trunc = rd(INPUTS[6][0], dtype=STR)
pvs = rd(INPUTS[7][0], dtype=STR)
tiers = rd(INPUTS[11][0], dtype=STR)
lvl = rd(INPUTS[13][0], dtype=STR)
evalu = rd(INPUTS[14][0], dtype=STR)
struct = rd(INPUTS[15][0], dtype=STR)
progs = rd(INPUTS[16][0], dtype=STR)
cohort_tests = rd(INPUTS[17][0])
de = rd(INPUTS[18][0])
rep10 = rd(INPUTS[20][0], dtype=STR)
liana = rd(INPUTS[25][0])
cov = rd(INPUTS[26][0], dtype=STR)
tqc = rd(INPUTS[27][0])
immune = rd(INPUTS[24][0])

TIER_EXPECT = {"robust-high": 4, "uncertain": 28}
L1_EXPECT = {"DN_STRUCTURE_SUPPORTED": 4, "DN_STRUCTURE_NOT_SUPPORTED": 23,
             "DN_STRUCTURE_NOT_EVALUABLE": 5}
L2_EXPECT = {"DN_STATE_SUPPORTED": 26, "DN_STATE_NOT_SUPPORTED": 1,
             "DN_STATE_NOT_EVALUABLE": 5}

assert len(s08p) == 32, f"expected 32 Stage-08 primary patients, got {len(s08p)}"
assert s08p.patient.is_unique, "duplicate patient in Stage-08 primary"
got_tier = tiers.final_tier.value_counts().to_dict()
assert got_tier == TIER_EXPECT, f"tier counts changed: {got_tier}"
assert got_tier.get("robust-low", 0) == 0, "robust-low appeared"
got_l1 = lvl.level1_structure.value_counts().to_dict()
assert got_l1 == L1_EXPECT, f"Level-1 counts changed: {got_l1}"
got_l2 = lvl.level2_state.value_counts().to_dict()
assert got_l2 == L2_EXPECT, f"Level-2 counts changed: {got_l2}"
assert (lvl.level3_cnv == "CNV_SUBCLONE_NOT_EVALUABLE").all(), "Level-3 changed"
assert len(lvl) == 32 and set(lvl.patient) == set(s08p.patient), "patient set changed"
print("frozen invariants hold:")
print(f"  tiers   {got_tier}")
print(f"  Level-1 {got_l1}")
print(f"  Level-2 {got_l2}")
print(f"  Level-3 CNV_SUBCLONE_NOT_EVALUABLE x {int((lvl.level3_cnv=='CNV_SUBCLONE_NOT_EVALUABLE').sum())}")
print(f"  patients {len(s08p)}")

# %% [markdown]
# ## Step 4 — Assemble the patient evidence matrix
#
# Axes A–F kept structurally separate. Every field is a verbatim copy from a frozen
# artifact or a predeclared derivation. No transformation beyond renaming and sign
# extraction; no aggregation across axes.

# %%
sp = evalu[evalu.denominator == "primary"].set_index("patient")
ss = evalu[evalu.denominator == "sensitivity"].set_index("patient")
pp = progs[progs.denominator == "primary"].set_index("patient")
b = boot.set_index("patient")
t = tiers.set_index("patient")
lv = lvl.set_index("patient")
cv = cov.set_index("patient")
s8p = s08p.set_index("patient")
s8s = s08s.set_index("patient")
r10 = rep10.drop_duplicates("patient").set_index("patient")

rows = []
for pid in sorted(s08p.patient):
    P, S, T, L, C = s8p.loc[pid], s8s.loc[pid], t.loc[pid], lv.loc[pid], cv.loc[pid]
    E, F = sp.loc[pid], (ss.loc[pid] if pid in ss.index else None)
    G = pp.loc[pid] if pid in pp.index else None
    flags = {c: bool(T[c]) for c in SY.FLAG_COLUMNS if c in T.index}
    n_flags = int(sum(flags.values()))
    l1s, l2s, l3s = L.level1_structure, L.level2_state, L.level3_cnv
    ap = SY.direction_sign(G["antigen_presentation_delta_matched"]) if G is not None else "not_evaluable"
    ox = SY.direction_sign(G["oxphos_delta_matched"]) if G is not None else "not_evaluable"
    ifn = SY.direction_sign(G["interferon_delta_matched"]) if G is not None else "not_evaluable"
    rec = dict(
        # --- identity -------------------------------------------------------
        patient_id=pid, cohort=P.cohort, n_samples=int(P.n_samples),
        repeated_sample_flag=bool(int(P.n_samples) > 1),
        in_paper_cohort=(pid != "25183"),
        # --- A. measurement -------------------------------------------------
        n_primary_cells=int(P.n_cells),
        n_sensitivity_cells=int(S.n_cells) if pid in s8s.index else np.nan,
        observed_dn_primary=float(P.observed_double_negative_fraction),
        observed_dn_sensitivity=float(S.observed_double_negative_fraction) if pid in s8s.index else np.nan,
        bcma_detection=float(P.bcma_detect), gprc5d_detection=float(P.gprc5d_detect),
        provisional_measurement_tier=T.final_tier,
        bootstrap_ci_lower=float(b.loc[pid, "dn_ci_lo"]),
        bootstrap_ci_upper=float(b.loc[pid, "dn_ci_hi"]),
        enrichment_depth_conditioned=float(P.enrichment_stratified),
        enrichment_unconditioned=float(P.enrichment_unconditioned),
        perm_p_depth_stratified=float(P.perm_p),
        denominator_unstable=bool(T.uncertain_denominator),
        depth_sensitive=bool(T.uncertain_depth), low_n=bool(T.uncertain_low_n),
        repeated_sample_unstable=bool(T.uncertain_repeated_sample),
        null_scheme_sensitive=bool(T.uncertain_null_scheme),
        dropout_compatible=bool(T.uncertain_dropout_compatible),
        intermediate_dn=bool(T.uncertain_intermediate_dn),
        n_uncertainty_flags=n_flags,          # a COUNT of flags, never an evidence score
        # --- B. structure ---------------------------------------------------
        level1_state=l1s, level1_evaluable=bool(E.evaluable),
        level1_not_evaluable_reason=("" if bool(E.evaluable) else str(E.reason)),
        morans_i_primary=float(E.morans_i) if pd.notna(E.morans_i) else np.nan,
        morans_i_sensitivity=float(F.morans_i) if F is not None and pd.notna(F.morans_i) else np.nan,
        knn_dn_frac_primary=float(E.knn_dn_frac) if pd.notna(E.knn_dn_frac) else np.nan,
        best_cluster_enr_primary=float(E.best_cluster_enr) if pd.notna(E.best_cluster_enr) else np.nan,
        depth_stratified_p_primary=float(E.perm_p_depth_stratified) if pd.notna(E.perm_p_depth_stratified) else np.nan,
        depth_stratified_p_sensitivity=float(F.perm_p_depth_stratified) if F is not None and pd.notna(F.perm_p_depth_stratified) else np.nan,
        unconditioned_p_primary=float(E.perm_p_unconditioned) if pd.notna(E.perm_p_unconditioned) else np.nan,
        n_depth_bins=int(E.n_depth_bins) if pd.notna(E.n_depth_bins) else np.nan,
        repeated_structure_consistent=(str(r10.loc[pid, "repeated_sample_status"])
                                       if pid in r10.index else "not_assessable"),
        # --- C. phenotype ---------------------------------------------------
        level2_state=l2s, level2_evaluable=bool(E.evaluable),
        phenotype_compatibility_note=(SY.PHENOTYPE_COMPATIBILITY_NOTE
                                      if l2s == "DN_STATE_SUPPORTED" else ""),
        antigen_presentation_direction=ap, oxphos_direction=ox, interferon_direction=ifn,
        secretory_direction=ap,   # secretory read-out is the antigen-presentation programme
        n_de_padj05_primary=int(E.n_de_padj05) if pd.notna(E.n_de_padj05) else np.nan,
        depth_ratio_matched=float(G["depth_ratio_matched"]) if G is not None and pd.notna(G["depth_ratio_matched"]) else np.nan,
        # --- D. genomic -----------------------------------------------------
        level3_state=l3s, cnv_evaluable=False,
        cnv_not_evaluable_reason=SY.CNV_NOT_EVALUABLE_REASON,
        # --- E. immune (cohort-level negative, recorded per row) -------------
        immune_context_summary=SY.IMMUNE_CONTEXT_SUMMARY,
        composition_supported=False, communication_supported=False,
        liana_independent_support=False,
        immune_confound_flag="receiver_state_confounded_and_antigen_circular",
        liana_evaluable=(pid != "25183"),
        # --- F. coverage ----------------------------------------------------
        anchor_pair_uncovered=float(C.uncovered_BCMA_GPRC5D),
        lowest_observed_uncovered_eligible_pair=C.greatest_coverage_pair_descriptive,
        lowest_observed_uncovered_eligible_pair_value=float(C.greatest_coverage_pair_uncovered),
        lowest_observed_uncovered_eligible_triple=C.greatest_coverage_triple_descriptive,
        lowest_observed_uncovered_eligible_triple_value=float(C.greatest_coverage_triple_uncovered),
        uncovered_bcma_alone=float(C.uncovered_TNFRSF17),
        uncovered_gprc5d_alone=float(C.uncovered_GPRC5D),
        gain_from_adding_gprc5d=float(C.gain_from_adding_GPRC5D),
        anchor_vs_alternative_note=SY.ANCHOR_VS_ALTERNATIVE_NOTE,
        coverage_depth_sensitive=bool(C.all_targets_depth_sensitive),
        coverage_repeat_sensitive=bool(int(P.n_samples) > 1),
        coverage_qc_note=SY.COVERAGE_QC_NOTE,
        eligible_targets=C.eligible_targets, not_evaluable_targets=C.not_evaluable_targets,
    )
    rec["measurement_interpretation"] = SY.measurement_interpretation(
        T.final_tier, n_flags, bool(T.uncertain_low_n))
    rec["biological_interpretation"] = SY.biological_interpretation(l1s, l2s)
    rec["main_uncertainty"] = SY.main_uncertainty(l1s, l2s, flags)
    rec["allowed_claim"] = SY.allowed_claim(
        SC.licensed_language(l1s, l2s, SC.CNV_NOT_EVALUABLE))
    rec["prohibited_claim"] = SY.prohibited_claim(l1s, l2s)
    rec["evidence_profile"] = SY.evidence_profile(T.final_tier, l1s, l2s, l3s)
    rows.append(rec)

M = pd.DataFrame(rows)
assert len(M) == 32 and M.patient_id.is_unique, "patient row invariant violated"
M.to_csv(OUT / "stage12_patient_evidence_matrix.csv", index=False)
print(f"wrote stage12_patient_evidence_matrix.csv — {M.shape[0]} rows x {M.shape[1]} columns")
print(M.provisional_measurement_tier.value_counts().to_dict())

# %% [markdown]
# ## Step 5 — Concordance / discordance
#
# The analytical heart of Stage 12. Counting only — no association test is run: with 4 and 4
# positives among 32 a formal test would be badly underpowered, and the contingency table is
# the honest presentation. Disjointness is reported as an observed fact, not a significant one.

# %%
M["L1"] = M.level1_state.str.replace("DN_STRUCTURE_", "", regex=False)
M["L2"] = M.level2_state.str.replace("DN_STATE_", "", regex=False)

ct_l1 = pd.crosstab(M.provisional_measurement_tier, M.L1)
ct_l2 = pd.crosstab(M.provisional_measurement_tier, M.L2)
ct_ll = pd.crosstab(M.L1, M.L2)
joint = (M.groupby(["provisional_measurement_tier", "L1", "L2"]).size()
         .reset_index(name="n").sort_values("n", ascending=False))
print("=== tier x Level-1 ==="); print(ct_l1)
print("\n=== tier x Level-2 ==="); print(ct_l2)
print("\n=== Level-1 x Level-2 ==="); print(ct_ll)
print("\n=== occupied joint profiles ==="); print(joint.to_string(index=False))

long = []
for name, ct in (("tier_x_level1", ct_l1), ("tier_x_level2", ct_l2), ("level1_x_level2", ct_ll)):
    for r in ct.index:
        for c in ct.columns:
            long.append(dict(table=name, row=r, column=c, n=int(ct.loc[r, c])))
for _, r in joint.iterrows():
    long.append(dict(table="joint_tier_level1_level2",
                     row=f"{r.provisional_measurement_tier}|{r.L1}", column=r.L2, n=int(r.n)))
CM = pd.DataFrame(long)
CM.to_csv(OUT / "stage12_concordance_matrix.csv", index=False)

n_occupied = len(joint)
ne_coupled = int(((M.L1 == "NOT_EVALUABLE") != (M.L2 == "NOT_EVALUABLE")).sum())
rh = set(M[M.provisional_measurement_tier == "robust-high"].patient_id)
l1s_ = set(M[M.L1 == "SUPPORTED"].patient_id)
converge = M[(M.provisional_measurement_tier == "robust-high") &
             (M.L1 == "SUPPORTED") & (M.L2 == "SUPPORTED")]
print(f"\noccupied cells: {n_occupied} of 18")
print(f"patients NOT_EVALUABLE on exactly one of L1/L2: {ne_coupled} (0 = perfectly coupled)")
print(f"robust-high n L1-supported overlap: {len(rh & l1s_)} (0 = disjoint)")
print(f"convergent robust-high + L1 + L2: {len(converge)}")
assert n_occupied == 5 and ne_coupled == 0 and len(rh & l1s_) == 0 and len(converge) == 0
print("wrote stage12_concordance_matrix.csv")

# %% [markdown]
# ## Step 6 — Repeated-patient summary
#
# Ranges, never means. `60359` has no primary-denominator cells and is not a repeated
# patient in this analysis, so 7 patients qualify, not 8.

# %%
rep_pat = sorted(M[M.repeated_sample_flag].patient_id)
print("repeated patients in the primary denominator:", rep_pat, f"(n={len(rep_pat)})")
assert len(rep_pat) == 7, f"expected 7 repeated-primary patients, got {len(rep_pat)}"

rr = []
for pid in rep_pat:
    s = rep08[rep08.patient == pid].copy()
    s["obs_dn"] = pd.to_numeric(s.obs_dn, errors="coerce")
    s["n_cells"] = pd.to_numeric(s.n_cells, errors="coerce")
    row = M[M.patient_id == pid].iloc[0]
    tau = [c for c in ("tier_tau020", "tier_tau025", "tier_tau033") if c in tiers.columns]
    tv = tiers[tiers.patient == pid][tau].iloc[0].tolist() if tau else []
    rr.append(dict(
        patient_id=pid, n_samples=int(row.n_samples),
        total_primary_cells=int(row.n_primary_cells),
        per_sample_cell_counts=";".join(f"{a}:{int(n)}" for a, n in zip(s["sample"], s.n_cells.fillna(0))),
        dn_min=float(s.obs_dn.min()), dn_max=float(s.obs_dn.max()),
        dn_range=float(s.obs_dn.max() - s.obs_dn.min()),
        measurement_tier=row.provisional_measurement_tier,
        measurement_tier_stability=("stable across TAU_HIGH 0.20/0.25/0.33"
                                    if len(set(tv)) <= 1 else "varies across TAU_HIGH"),
        level1_repeated_consistency=row.repeated_structure_consistent,
        coverage_min=float(row.lowest_observed_uncovered_eligible_pair_value),
        coverage_max=float(row.anchor_pair_uncovered),
        coverage_range=float(row.anchor_pair_uncovered - row.lowest_observed_uncovered_eligible_pair_value),
        small_denominator_flag=bool((s.n_cells < 20).any()),
        small_denominator_note=("one or more samples below the frozen 20-cell floor; "
                                "per-sample DN values are not interpretable in isolation"
                                if (s.n_cells < 20).any() else ""),
    ))
RP = pd.DataFrame(rr)
RP.to_csv(OUT / "stage12_repeated_patient_summary.csv", index=False)
print(RP[["patient_id", "n_samples", "total_primary_cells", "dn_min", "dn_max",
          "dn_range", "small_denominator_flag"]].round(3).to_string(index=False))
print("\nNOTE: 60359 has zero primary-denominator cells and is excluded by that fact, "
      "not by a filter.")

# %% [markdown]
# ## Step 7 — Cohort summary, claim ladder, uncertainty register
#
# Pre-populated from the design; every value was verified against frozen artifacts during
# the design pass and is re-asserted below where cheap.

# %%
V = SY.STRENGTH_VOCABULARY
cohort = [
 dict(domain="measurement",
      question="Is the observed DN measurement technically credible?",
      result=("observed DN common; median 0.335 (0.017-0.783) over 32 patients and "
              "21,906 primary-denominator cells"),
      strength_of_evidence="STRONG",
      major_control="denominator, depth, truncate-10k, repeated-sample and threshold sensitivity",
      major_limitation="GPRC5D pooled technical-zero 0.620; WashU 10,000-UMI deposit censoring",
      allowed_interpretation="observed transcript-level DN burden",
      prohibited_interpretation="true target-negative clone prevalence"),
 dict(domain="co-loss",
      question="Is BCMA/GPRC5D co-negativity enriched beyond technical expectation?",
      result=("unconditioned enrichment median 1.052 (max 4.606) collapses to 1.009 "
              "(max 1.750) under the cohort-specific depth-stratified null; 4 of 32 patients "
              "retain significance, all in the deepest cohort"),
      strength_of_evidence="NOT_SUPPORTED",
      major_control="cohort-specific depth-stratified permutation null",
      major_limitation="significance tracks depth and power, not necessarily biology",
      allowed_interpretation="most apparent co-negativity enrichment is attributable to sequencing depth",
      prohibited_interpretation="coordinated biological antigen co-loss"),
 dict(domain="structure",
      question="Are DN cells non-randomly organised beyond the depth null?",
      result="4 supported / 23 not supported / 5 not evaluable",
      strength_of_evidence="SUPPORTED_WITH_CAVEATS",
      major_control="depth-stratified within-patient permutation with adaptive depth bins",
      major_limitation=("conservative by design; 5 patients not evaluable; MMRF_1640 shows "
                        "Moran's I 0.470 with unconditioned p 0.001 and depth-stratified p 0.499"),
      allowed_interpretation="non-random DN organization in a minority of patients",
      prohibited_interpretation="subclone; pre-existing resistant clone"),
 dict(domain="phenotype",
      question="Do DN cells share a transcriptional phenotype?",
      result=("cohort-level DN-lower antigen presentation, OXPHOS and interferon under the "
              "frozen both-denominators rule; 190 DE genes significant under both; "
              "per-patient support 26 of 27 evaluable"),
      strength_of_evidence="SUPPORTED_WITH_CAVEATS",
      major_control="depth-matched cells before pseudobulk; both-denominators reproducibility rule",
      major_limitation=("confounded with a broad secretory/differentiation shift; BCMA and "
                        "GPRC5D are themselves secretory-pathway-dependent surface proteins"),
      allowed_interpretation="cohort-level DN-associated transcriptional state",
      prohibited_interpretation="antigen-specific escape mechanism; individualized phenotype risk"),
 dict(domain="genomic",
      question="Is there genomic evidence for a distinct DN subclone?",
      result="32 of 32 patients CNV_SUBCLONE_NOT_EVALUABLE",
      strength_of_evidence="NOT_EVALUABLE",
      major_control="healthy-donor negative control, run before any disease interpretation",
      major_limitation="donor false-positive rate spans 0.0-50.6% at z>3",
      allowed_interpretation="genomic subclone evidence is not evaluable",
      prohibited_interpretation="absence of a subclone; any use of the word subclone"),
 dict(domain="immune",
      question="Does immune context explain or track DN burden?",
      result=("0 of 28 composition tests at BH<0.10 (smallest 0.49); targeted LR panel is "
              "receiver-state confounded; LIANA 1 of 87 at BH<0.10 and antigen-circular"),
      strength_of_evidence="NOT_SUPPORTED",
      major_control="receptor-only falsification; permanent antigen-circularity test",
      major_limitation="n~32 with a confounder correlated with the predictor",
      allowed_interpretation="no robust independent immune association",
      prohibited_interpretation="immune-evasion mechanism; PDCD1-CD274 axis"),
 dict(domain="coverage",
      question="How does observed uncovered fraction change under alternative target sets?",
      result=("eligible alternatives lower the observed uncovered fraction in 32 of 32 "
              "patients (median advantage 0.098); no target is depth-robust; 2 of 7 targets "
              "are COVERAGE_NOT_EVALUABLE"),
      strength_of_evidence="EXPLORATORY",
      major_control="target-specific expression-matched technical-zero floors; predeclared eligibility gate",
      major_limitation=("detection rates differ 1.8-2.8x; GPRC5D dropout; SDC1 and TNFRSF17 "
                        "selection dependence"),
      allowed_interpretation=("greatest observed transcript-level malignant-cell coverage "
                              "among evaluated combinations"),
      prohibited_interpretation="optimal / recommended / best pair; GPRC5D is redundant; any safety claim"),
]
CS = pd.DataFrame(cohort)
assert list(CS.domain) == ["measurement", "co-loss", "structure", "phenotype", "genomic",
                           "immune", "coverage"]
assert CS.strength_of_evidence.isin(V).all()
CS.to_csv(OUT / "stage12_cohort_summary.csv", index=False)
print(CS[["domain", "strength_of_evidence"]].to_string(index=False))

# %%
ladder = [
 ("observed DN burden exists at baseline",
  "per-patient DN over a defended, antigen-circularity-controlled denominator",
  "median 0.335 (0.017-0.783); 32 patients; 21,906 primary cells", "STRONG",
  "observed transcript-level double-negative fraction",
  "true antigen-negative clone prevalence"),
 ("the DN measurement is technically robust",
  "survival of every frozen sensitivity analysis",
  "4 of 32 patients carry zero uncertainty flags; 0 robust-low", "SUPPORTED_WITH_CAVEATS",
  "measurement-robust for 4 patients (measurement statement only)",
  "biologically high-risk patient"),
 ("DN co-loss is enriched beyond depth",
  "conditioned enrichment above 1 under a depth-stratified null",
  "median conditioned enrichment 1.009; 4 of 32 significant, all in the deepest cohort",
  "NOT_SUPPORTED",
  "4 patients retain significance in the cohort where the test has power",
  "coordinated biological antigen co-loss across the cohort"),
 ("DN cells are non-randomly organised",
  "Level-1 support under the depth-stratified null, both denominators",
  "4 supported / 23 not supported / 5 not evaluable", "SUPPORTED_WITH_CAVEATS",
  "non-random DN organization in 4 patients",
  "escape subclone; pre-existing resistant clone"),
 ("a DN-associated transcriptional phenotype exists",
  "cohort-level programmes reproducible under both denominators with consistent sign",
  "antigen presentation, OXPHOS and interferon DN-lower; 190 DE genes under both denominators",
  "SUPPORTED_WITH_CAVEATS",
  "cohort-level DN-associated transcriptional state",
  "individualized phenotype-based risk; Level-2 predicts patient risk"),
 ("DN reflects an antigen-specific escape state",
  "separation of the phenotype from a general secretory/differentiation shift",
  "not separable: BCMA and GPRC5D are secretory-pathway-dependent surface proteins",
  "NOT_SUPPORTED", "(no wording permitted)", "antigen-specific escape mechanism"),
 ("a genomic DN subclone exists",
  "CNV evidence from a calibrated assay",
  "32 of 32 CNV_SUBCLONE_NOT_EVALUABLE; donor negative control failed", "NOT_EVALUABLE",
  "genomic subclone evidence is not evaluable",
  "confirmed escape subclone; no subclone exists; CNV-negative"),
 ("an immune-evasion mechanism operates",
  "independent immune association surviving correction and falsification",
  "0 of 28 composition tests at BH<0.10; LR receiver-state confounded; LIANA 1 of 87 and circular",
  "NOT_SUPPORTED", "no robust independent immune association",
  "immune escape mechanism; PDCD1-CD274 immune-evasion axis"),
 ("alternative target sets reduce the observed uncovered fraction",
  "per-patient uncovered fractions over eligible targets",
  "lower in 32 of 32 patients; median advantage 0.098", "EXPLORATORY",
  "greatest observed transcript-level malignant-cell coverage among evaluated combinations",
  "therapeutically superior combination"),
 ("an alternative target combination is therapeutically superior",
  "comparable detection plus efficacy, protein-level density and toxicity evidence",
  "detection differs 1.8-2.8x; no target is depth-robust; 2 of 7 not evaluable",
  "NOT_SUPPORTED", "(no wording permitted)",
  "optimal target pair; best therapeutic combination; GPRC5D is redundant"),
 ("normal-marrow data establish safety",
  "normal-tissue atlas covering the relevant off-tumour sites",
  "8 donors, 647 plasma cells, marrow only; GPRC5D keratinized-tissue liability unobservable here",
  "NOT_EVALUABLE", "normal marrow expression context",
  "clinically safe; normal tissue safe"),
 ("a final patient risk classifier is justified",
  "convergent, commensurable evidence across axes",
  ("no patient is simultaneously measurement-robust-high, Level-1 supported and Level-2 "
   "supported; only 5 of 18 joint cells occupied; Level-2 supported in 26 of 27 evaluable"),
  "NOT_SUPPORTED", "separate evidence axes reported side by side",
  "high-risk patient; low-risk patient; composite risk score"),
]
CL = pd.DataFrame(ladder, columns=["claim", "evidence_required", "frozen_evidence",
                                   "current_status", "allowed_wording", "prohibited_wording"])
assert CL.current_status.isin(V).all() and len(CL) == 12
CL.to_csv(OUT / "stage12_claim_ladder.csv", index=False)
print(CL[["claim", "current_status"]].to_string(index=False))

# %%
unc = [
 ("measurement", "dropout / technical-zero burden",
  "GPRC5D pooled technical-zero 0.620 vs TNFRSF17 0.276; 15 of 32 patients flagged dropout_compatible"),
 ("measurement", "sequencing depth",
  "GPRC5D detection spans 9.4x across cohorts and tracks depth monotonically; DN moves the opposite way"),
 ("measurement", "WashU 10,000-UMI deposit censoring",
  "truncate-10k leaves WashU exactly unchanged and raises MMRF DN by a mean 0.059"),
 ("measurement", "denominator sensitivity",
  "primary to sensitivity median DN shift +0.032; 12 of 32 patients move more than 5 points"),
 ("measurement", "low malignant-cell count",
  "4 patients flagged low_n; cell counts vary widely across patients"),
 ("measurement", "repeated-sample instability",
  "7 repeated patients; 2 discordant; within-patient DN range up to 0.571 on very small samples"),
 ("measurement", "null-scheme sensitivity",
  "cohort vs global depth bins: median |delta| 0.0011, and MMRF_1413 crosses 1.0 between schemes"),
 ("biological", "absence of Level-1 structure",
  "23 of 32 patients show no non-random DN organization; 5 not evaluable"),
 ("biological", "weakly discriminative Level 2",
  "26 of 27 evaluable patients satisfy the frozen per-patient rule"),
 ("biological", "phenotype specificity",
  "the DN-associated phenotype cannot be separated from a general secretory/differentiation shift"),
 ("biological", "CNV non-evaluability",
  "32 of 32 CNV_SUBCLONE_NOT_EVALUABLE; carries no evidential weight in either direction"),
 ("external_validity", "RNA versus protein",
  "CAR-T binds surface protein; BCMA is actively shed by gamma-secretase and GPRC5D transcript "
  "correlates imperfectly with surface density"),
 ("external_validity", "baseline versus post-treatment selection",
  "this cohort is untreated baseline; no therapy selection has occurred and none can be inferred"),
 ("external_validity", "marrow versus whole body",
  "GPRC5D keratinized-tissue liability is structurally unobservable in a marrow dataset"),
 ("external_validity", "cohort and protocol confounding",
  "depth, site, chemistry and deposit censoring move together and cannot be separated in this design"),
]
UR = pd.DataFrame(unc, columns=["uncertainty_class", "source", "frozen_evidence"])
assert set(UR.uncertainty_class) == {"measurement", "biological", "external_validity"}
UR.to_csv(OUT / "stage12_uncertainty_register.csv", index=False)
print(UR.uncertainty_class.value_counts().to_dict(),
      "\n(three separate registers; no combined uncertainty score)")

# %% [markdown]
# ## Step 8 — Figures
#
# Declared sort/order rules only. No clustering, no seriation, no optimisation of visual
# pattern, and no ordinal colour ramp on a categorical axis.

# %%
CAT = {"SUPPORTED": "#2c7fb8", "NOT_SUPPORTED": "#f0f0f0", "NOT_EVALUABLE": "#bdbdbd",
       "robust-high": "#238b45", "uncertain": "#f0f0f0", "flag": "#d94801",
       "no_flag": "#f0f0f0", "agree": "#2c7fb8", "disagree": "#d94801",
       "not_assessable_for_agreement": "#bdbdbd", "not_assessable": "#bdbdbd"}

F2 = M.sort_values(["provisional_measurement_tier", "patient_id"],
                   ascending=[True, True]).reset_index(drop=True)  # tier, then patient_id
cols = [("measurement tier", F2.provisional_measurement_tier),
        ("Level 1 structure", F2.L1), ("Level 2 phenotype", F2.L2),
        ("Level 3 genomic", pd.Series(["NOT_EVALUABLE"] * len(F2))),
        ("repeated-sample", F2.repeated_structure_consistent)]
flagcols = [("low n", F2.low_n), ("dropout-compatible", F2.dropout_compatible),
            ("denominator", F2.denominator_unstable), ("depth", F2.depth_sensitive),
            ("repeated-sample flag", F2.repeated_sample_unstable),
            ("null scheme", F2.null_scheme_sensitive), ("intermediate DN", F2.intermediate_dn)]
labels = [c[0] for c in cols] + [c[0] for c in flagcols]
grid = []
for _, s in cols:
    grid.append([CAT.get(str(v), "#f0f0f0") for v in s])
for _, s in flagcols:
    grid.append([CAT["flag"] if bool(v) else CAT["no_flag"] for v in s])

fig, ax = plt.subplots(figsize=(9.5, 11))
for j, colcolours in enumerate(grid):
    for i, c in enumerate(colcolours):
        ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor=c, edgecolor="white", lw=0.7))
ax.set_xlim(0, len(grid)); ax.set_ylim(0, len(F2)); ax.invert_yaxis()
ax.set_xticks(np.arange(len(labels)) + 0.5)
ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
ax.set_yticks(np.arange(len(F2)) + 0.5)
ax.set_yticklabels(F2.patient_id, fontsize=7)
ax.axvline(len(cols), color="black", lw=1.2)
for s in ax.spines.values():
    s.set_visible(False)
ax.set_title("Figure 2 — Patient evidence matrix (categorical)\n"
             "sorted by provisional measurement tier, then patient_id; colours encode "
             "category identity only", fontsize=9)
ax.legend(handles=[Patch(facecolor=CAT["SUPPORTED"], label="supported / agree"),
                   Patch(facecolor=CAT["NOT_SUPPORTED"], label="not supported / absent"),
                   Patch(facecolor=CAT["NOT_EVALUABLE"], label="not evaluable"),
                   Patch(facecolor=CAT["robust-high"], label="robust-high (measurement)"),
                   Patch(facecolor=CAT["flag"], label="uncertainty flag present")],
          loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=7, frameon=False)
fig.tight_layout(); fig.savefig(FIG / "figure2_patient_evidence_matrix.png", dpi=200,
                                bbox_inches="tight"); plt.close(fig)
print("wrote figure2_patient_evidence_matrix.png")

# %%
mt = cohort_tests[cohort_tests.effect == "matched"]
raw = cohort_tests[cohort_tests.effect == "raw"]
order = ["antigen_presentation", "oxphos", "interferon", "upr", "myc", "stress",
         "gamma_secretase"]
fig, ax = plt.subplots(figsize=(9, 5.2))
for i, prog in enumerate(order):
    for den, off, mk in (("primary", -0.16, "o"), ("sensitivity", 0.16, "s")):
        r = mt[(mt.program == prog) & (mt.denominator == den)].iloc[0]
        rr_ = raw[(raw.program == prog) & (raw.denominator == den)].iloc[0]
        ax.scatter(rr_.median_delta, i + off, s=34, facecolor="none",
                   edgecolor="#bdbdbd", marker=mk, zorder=2)
        repro = (mt[(mt.program == prog) & (mt.denominator == "primary")].p_BH.iloc[0] < 0.10
                 and mt[(mt.program == prog) & (mt.denominator == "sensitivity")].p_BH.iloc[0] < 0.10)
        ax.scatter(r.median_delta, i + off, s=58,
                   color="#2c7fb8" if repro else "#737373", marker=mk, zorder=3)
        ax.text(r.median_delta, i + off, f"  BH {r.p_BH:.3f}", va="center", fontsize=6.5)
ax.axvline(0, color="black", lw=0.8)
ax.set_yticks(range(len(order)))
ax.set_yticklabels([p.replace("_", " ") for p in order], fontsize=9)
ax.invert_yaxis(); ax.set_xlabel("median per-patient delta (DN minus comparator)", fontsize=9)
ax.set_title("Figure 4 — Stage-10 cohort phenotype: seven predeclared programmes\n"
             "filled = depth-matched, open grey = raw (unmatched); circle = primary, "
             "square = sensitivity\nblue = reproducible under the both-denominators rule "
             "(BH<0.10 in BOTH, consistent sign)", fontsize=8.5)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(FIG / "figure4_stage10_cohort_phenotype.png", dpi=200)
plt.close(fig)
print("wrote figure4_stage10_cohort_phenotype.png")

# %%
combos = sorted(set(M.lowest_observed_uncovered_eligible_pair))       # alphabetical, never by value
data = [("TNFRSF17+GPRC5D (anchor)", M.anchor_pair_uncovered.values)]
for c in combos:
    data.append((c, M.loc[M.lowest_observed_uncovered_eligible_pair == c,
                          "lowest_observed_uncovered_eligible_pair_value"].values))
fig, ax = plt.subplots(figsize=(9, 5))
bp = ax.boxplot([d[1] for d in data], vert=True, widths=0.55, patch_artist=True)
for i, patch in enumerate(bp["boxes"]):
    patch.set_facecolor("#fdd0a2" if i == 0 else "#deebf7"); patch.set_edgecolor("#525252")
for med in bp["medians"]:
    med.set_color("black")
ax.set_xticklabels([d[0] for d in data], rotation=20, ha="right", fontsize=8)
ax.set_ylabel("observed transcript-level uncovered fraction", fontsize=9)
ax.set_title("Figure 5 — Multi-antigen coverage (descriptive; combinations ordered "
             "alphabetically, NOT by value)\n"
             "GPRC5D: COVERAGE_NOT_EVALUABLE (technical-zero 0.62 vs a 0.50 gate fixed "
             "beforehand) - remains the frozen Stage-08 anchor\n"
             "SDC1: excluded on circularity | TNFRSF17: identical selection-dependence, "
             "retained as anchor (disclosure) | no target is depth-robust", fontsize=7.8)
ax.text(0.5, -0.30, "Lower uncovered fraction for alternatives is a detection-rate artifact "
        "(targets detected 1.8-2.8x more often than GPRC5D).\nNo combination is optimal, "
        "recommended or best.", transform=ax.transAxes, ha="center", fontsize=7.5, color="#525252")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(FIG / "figure5_multi_antigen_coverage.png", dpi=200,
                                bbox_inches="tight"); plt.close(fig)
print("wrote figure5_multi_antigen_coverage.png")

# %%
fig, ax = plt.subplots(figsize=(10, 4.4))
lanes = [("A. Measurement\nrobustness", "4 robust-high / 28 uncertain / 0 robust-low", "#2c7fb8"),
         ("B. DN structure\n(Level 1)", "4 / 23 / 5", "#2c7fb8"),
         ("C. DN phenotype\n(Level 2)", "26 / 1 / 5\nweakly discriminative", "#2c7fb8"),
         ("D. Genomic\n(Level 3)", "32 NOT_EVALUABLE", "#bdbdbd")]
ctx = [("E. Immune context\n(Stages 11 / 11b)", "no robust independent signal", "#f0f0f0"),
       ("F. Multi-antigen coverage\n(Stage 08c)", "descriptive; no target depth-robust", "#f0f0f0")]
for i, (name, val, col) in enumerate(lanes):
    ax.add_patch(plt.Rectangle((i * 2.5, 2.0), 2.2, 1.3, facecolor=col, alpha=0.35,
                               edgecolor="#525252"))
    ax.text(i * 2.5 + 1.1, 2.95, name, ha="center", fontsize=8.5, weight="bold")
    ax.text(i * 2.5 + 1.1, 2.35, val, ha="center", fontsize=7.2)
for i, (name, val, col) in enumerate(ctx):
    ax.add_patch(plt.Rectangle((i * 5 + 1.2, 0.2), 4.6, 1.0, facecolor=col,
                               edgecolor="#525252", linestyle="--"))
    ax.text(i * 5 + 3.5, 0.85, name, ha="center", fontsize=8.5, weight="bold")
    ax.text(i * 5 + 3.5, 0.42, val, ha="center", fontsize=7.2)
ax.text(5.0, 3.75, "Six parallel evidence axes — NOT a linear escalation and never summed",
        ha="center", fontsize=9.5, weight="bold")
ax.text(5.0, 1.55, "no patient is simultaneously measurement-robust-high, Level-1 supported "
        "and Level-2 supported  |  only 5 of 18 joint cells occupied",
        ha="center", fontsize=7.6, color="#a63603")
ax.set_xlim(-0.4, 10.4); ax.set_ylim(0, 4.1); ax.axis("off")
ax.set_title("Figure 1 — Evidence architecture", fontsize=10)
fig.tight_layout(); fig.savefig(FIG / "figure1_evidence_architecture.png", dpi=200)
plt.close(fig)
print("wrote figure1_evidence_architecture.png")

# %%
fig, ax = plt.subplots(figsize=(8, 5.4))
mk = {"SUPPORTED": "o", "NOT_SUPPORTED": "s", "NOT_EVALUABLE": "X"}
for state, g in M.groupby("L1"):
    for tier, gg in g.groupby("provisional_measurement_tier"):
        ax.scatter(gg.observed_dn_primary, gg.morans_i_primary, marker=mk[state], s=78,
                   facecolor=("#238b45" if tier == "robust-high" else "none"),
                   edgecolor=("#238b45" if tier == "robust-high" else "#525252"),
                   linewidths=1.3, label=f"L1 {state} / {tier}", zorder=3)
r = M[M.patient_id == "MMRF_1640"]
if len(r):
    ax.annotate("MMRF_1640\nMoran's I 0.470\nunconditioned p 0.001\ndepth-stratified p 0.499",
                (float(r.observed_dn_primary.iloc[0]), float(r.morans_i_primary.iloc[0])),
                textcoords="offset points", xytext=(14, -8), fontsize=7,
                arrowprops=dict(arrowstyle="->", lw=0.8, color="#a63603"), color="#a63603")
ax.set_xlabel("observed DN fraction (primary denominator)", fontsize=9)
ax.set_ylabel("Moran's I on the DN label (primary)", fontsize=9)
ax.set_title("Figure 3 — Measurement versus Level-1 structure: the axes are discordant\n"
             "all 4 measurement robust-high patients are Level-1 NOT_SUPPORTED", fontsize=9)
ax.legend(fontsize=6.5, frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(FIG / "figure3_measurement_vs_structure.png", dpi=200,
                                bbox_inches="tight"); plt.close(fig)
print("wrote figure3_measurement_vs_structure.png")

# %% [markdown]
# ## Step 9 — Narrative synthesis
#
# Sections 9 and 10 are generated from `stage12_claim_ladder.csv`, not written freehand, so
# the narrative cannot drift from the table. The central result is **discordance**; the text
# must not smooth it over.

# %%
def fmt(x, n=3):
    return f"{float(x):.{n}f}"


n_rh = int((M.provisional_measurement_tier == "robust-high").sum())
n_unc = int((M.provisional_measurement_tier == "uncertain").sum())
n_l1s = int((M.L1 == "SUPPORTED").sum())
n_l1ns = int((M.L1 == "NOT_SUPPORTED").sum())
n_l1ne = int((M.L1 == "NOT_EVALUABLE").sum())
n_l2s = int((M.L2 == "SUPPORTED").sum())
n_l2ns = int((M.L2 == "NOT_SUPPORTED").sum())
n_l2ne = int((M.L2 == "NOT_EVALUABLE").sum())
med_dn = M.observed_dn_primary.median()
gs_p = cohort_tests[(cohort_tests.program == "gamma_secretase") &
                    (cohort_tests.effect == "matched")]

lines = []
A = lines.append
A("# Stage 12 — final synthesis")
A("")
A(f"**Executed:** {EXEC['executed_utc']} · commit `{EXEC['commit_short']}` · "
  f"tag `{EXEC['nearest_tag']}` · environment `{EXEC['environment']}`")
A(f"**Binding design:** `docs/stage12_design.md` (SHA256 `{EXEC['design_sha256'][:16]}…`)")
A(f"**Upstream inputs:** 29 frozen artifacts, all verified against the committed manifest "
  f"by SHA256 and size.")
A("")
A("> **Stage 12 is synthesis only.** It ran no new statistical test, fitted no model, "
  "created no threshold, and produced no score. Every quantity below already existed in a "
  "frozen artifact.")
A("")
A("---")
A("")
A("## 1. What was measured")
A("")
A("For each of **32 patients**, the fraction of malignant plasma cells in the frozen clone "
  "denominator whose **raw UMI count is zero for both `TNFRSF17` (BCMA) and `GPRC5D`** — "
  "the observed double-negative (DN) fraction. Raw counts, `count > 0`, no normalisation "
  "and no imputation: whether a transcript is genuinely absent is the question, so any "
  "method that borrows expression from neighbouring cells would erase the measurement.")
A("")
A(f"The primary denominator is **{int(M.n_primary_cells.sum()):,} cells** "
  "(`CLONE_SUPPORTED`); the sensitivity denominator adds clone-compatible cells whose "
  f"dominant V transcript was not observed, giving **{int(M.n_sensitivity_cells.sum()):,} "
  "cells**. Both are reported for every patient and are never averaged or selected between.")
A("")
A(f"**Observed DN fraction: median {fmt(med_dn)}, range {fmt(M.observed_dn_primary.min())}–"
  f"{fmt(M.observed_dn_primary.max())}.** This is a measurement, not an escape probability.")
A("")
A("## 2. How robust the measurement is")
A("")
A(f"Provisional measurement tiers: **{n_rh} robust-high, {n_unc} uncertain, 0 robust-low**. "
  "These are *measurement* tiers — a `robust-high` label means the observed DN fraction "
  "survived the denominator, depth, repeated-sample, null-scheme and threshold sensitivity "
  "analyses, and nothing more.")
A("")
A("**Zero `robust-low` patients is itself a result.** At this depth, in this cohort, no "
  "patient can be called confidently low-escape.")
A("")
A("Three technical facts bound every number above:")
A("")
A("- **GPRC5D is measured far less reliably than BCMA.** Detection spans **9.4×** across "
  "cohorts (MMRF 0.364 / WU1 0.121 / WU2 0.039), tracks depth monotonically, and DN moves "
  "the opposite way. The pooled expression-matched technical-zero floor is **0.620** for "
  "GPRC5D against **0.276** for TNFRSF17. Most GPRC5D zeros in WashU are consistent with "
  "dropout.")
A("- **The WashU deposit was censored at 10,000 UMIs before deposit.** Truncating every "
  "cohort at 10k leaves WashU *exactly* unchanged and raises MMRF DN by a mean **+0.059** "
  "while lowering GPRC5D detection by 0.186 — the deposit's own bias, reproduced on demand.")
A("- **Denominator choice moves DN directionally.** Primary → sensitivity shifts DN by a "
  "median **+0.032**, with 12 of 32 patients moving more than 5 points. Adding the shallower "
  "V-unobserved cells *raises* DN, which is what a depth-driven metric does.")
A("")
A("**Dropout is bounded, not corrected.** No dropout-corrected DN estimate is produced or "
  "claimed, and the technical-zero floor may not be subtracted or divided to manufacture one.")
A("")
A("## 3. Whether DN cells are structured")
A("")
A(f"Level-1 (non-random DN organisation beyond the depth-stratified null): "
  f"**{n_l1s} supported, {n_l1ns} not supported, {n_l1ne} not evaluable**.")
A("")
A("Supported: " + ", ".join(f"`{p}`" for p in sorted(M[M.L1 == 'SUPPORTED'].patient_id)) + ".")
A("")
A("Level 1 licenses **non-random DN organization and nothing more** — not a transcriptional "
  "state, not a subclone, and not antigen-specificity.")
A("")
A("The worked example that shows why the depth null is load-bearing is **`MMRF_1640`**: "
  "Moran's I **0.470** on the DN label with an unconditioned permutation p of **0.001**, "
  "and a depth-stratified p of **0.499**. By every descriptive measure it looks like a "
  "structured DN population; conditioned on depth it is indistinguishable from noise. It is "
  "`DN_STRUCTURE_NOT_SUPPORTED`.")
A("")
A("Apparent co-loss enrichment collapses the same way: unconditioned enrichment median "
  "**1.052** (max 4.606) becomes **1.009** (max 1.750) under the cohort-specific "
  "depth-stratified null, and significant patients fall from 8 to **4 — all in the deepest "
  "cohort**, i.e. where the test has power rather than necessarily where the biology is.")
A("")
A("## 4. Whether DN cells share a phenotype")
A("")
A(f"Level-2 (compatibility with the cohort-level DN-associated programme): "
  f"**{n_l2s} supported, {n_l2ns} not supported, {n_l2ne} not evaluable**.")
A("")
A("> **The per-patient Level-2 state is weakly discriminative: 26 of 27 evaluable patients "
  "satisfy the frozen rule.** `DN_STATE_SUPPORTED` means *compatible with the cohort-level "
  "programme*, never strong patient-specific evidence of a distinct escape state. The rule "
  "was not retuned after this was observed.")
A("")
A("**The scientifically informative Level-2 result is the cohort-level phenotype.** Under "
  "the frozen both-denominators rule (BH < 0.10 under *both* the primary and sensitivity "
  "denominators, with consistent sign), three programmes are reproducibly **lower in DN "
  "cells**:")
A("")
A("| programme | BH (primary) | BH (sensitivity) | direction |")
A("|---|---:|---:|---|")
for prog in ("antigen_presentation", "oxphos", "interferon"):
    p_ = mt[(mt.program == prog) & (mt.denominator == "primary")].iloc[0]
    s_ = mt[(mt.program == prog) & (mt.denominator == "sensitivity")].iloc[0]
    A(f"| {prog.replace('_',' ')} | {p_.p_BH:.4f} | {s_.p_BH:.4f} | lower in DN |")
A("")
A("Pseudobulk differential expression gives **190 genes significant under both "
  "denominators**, dominated on the DN-lower side by ER/secretory machinery and mature "
  "plasma-cell identity genes (`SPCS1`, `SPCS2`, `SEC61B`, `UBE2J1`, `TMBIM6`, `MZB1`, "
  "`ITM2C`, `B2M`).")
A("")
A(f"> **Method.** {SY.STAGE10_DE_METHOD}")
A("")
A("**The DN phenotype is compatible with a less secretory, less differentiated plasma-cell "
  "state.** It does **not** establish an antigen-specific escape mechanism, and the two "
  "cannot be separated in this data: BCMA and GPRC5D are themselves "
  "secretory-pathway-dependent surface proteins, so a cell with a less active secretory "
  "programme would be expected to carry less of both.")
A("")
A("**The pre-registered γ-secretase hypothesis is not supported.** The five-gene programme "
  "(`NCSTN`, `PSEN1`, `APH1A`, `APH1B`, `PSENEN`) was frozen before any result and no gene "
  "was added afterwards. It **fails the frozen both-denominators rule**, and its direction "
  "is *negative* — opposite to the pre-registered γ-secretase-high prediction. Depth "
  "matching is what removed it: before matching it appeared significant and negative, which "
  "would have been a confident, wrong-signed claim.")
A("")
A("## 5. Whether a genomic subclone is demonstrated")
A("")
A(f"**No. All {len(M)} patients are `CNV_SUBCLONE_NOT_EVALUABLE`.**")
A("")
A("Stage 07 calibrated `infercnvpy` on eight healthy donors *before* inspecting any disease "
  "sample. Donor plasma-cell false-positive rates span **0.0%–50.6%** at z > 3, with one "
  "donor at a median z of **+3.03**. The method was rejected on its own negative control, so "
  "it contributes no evidence in either direction.")
A("")
A("> This is **not** `NOT_SUPPORTED`. The project can say *\"genomic subclone evidence is "
  "not evaluable\"* and cannot say *\"there is no subclone\"*. Only Level 3 licenses the word "
  "*subclone*, so that word is unavailable for every patient in this cohort.")
A("")
A("## 6. Whether immune context adds support")
A("")
A("**It does not.** Stage 11 and its Stage-11b LIANA verification arm are exploratory by "
  "declaration and non-tier-changing.")
A("")
A("- **Composition: 0 of 28 tests reach BH < 0.10** (smallest BH ≈ 0.49). Two `NK_core` "
  "associations reach raw significance with the same sign in all three cohorts, but both are "
  "fragile to re-parameterisation and `NK_core` was pre-flagged as depth-tracking.")
A("- **The targeted ligand–receptor panel is receiver-state confounded.** Tested alone, "
  "**11 of 15 receptors move down with DN burden**, across receptors with nothing "
  "biologically in common — the panel reads the Stage-10 plasma-cell-state shift through its "
  "receptor term.")
A("- **LIANA does not rescue it.** Of 87 interactions tested, **1** reaches BH < 0.10: "
  "`Myeloid TNFSF13B → TNFRSF17`. Its receptor is one of the two antigens that *define* the "
  "predictor, so the negative coefficient is **arithmetic**. Decomposition of all 12 "
  "raw-significant interactions gives 3 receiver-state-confounded, 9 not reproduced, and "
  "**0 surviving as LIANA-only exploratory findings**.")
A("")
A("**No immune-evasion mechanism is established**, and `PDCD1 → CD274` is not an "
  "immune-evasion axis. LIANA is never described as validating immune evasion.")
A("")
A("## 7. What multi-antigen transcript coverage shows")
A("")
A("Descriptive transcript coverage over seven candidate targets. **Measurement quality, not "
  "biology, is the binding constraint.**")
A("")
A(f"- Anchor BCMA+GPRC5D uncovered fraction: median **{fmt(M.anchor_pair_uncovered.median())}**; "
  f"BCMA alone **{fmt(M.uncovered_bcma_alone.median())}**; GPRC5D alone "
  f"**{fmt(M.uncovered_gprc5d_alone.median())}**.")
A(f"- Median gain from adding GPRC5D to BCMA: **{fmt(M.gain_from_adding_gprc5d.median())}**.")
A("- **Two of seven targets are `COVERAGE_NOT_EVALUABLE`**: `GPRC5D` on technical-zero "
  "burden (0.62 against a 0.50 gate fixed before any floor existed) and `SDC1` on "
  "circularity. **No target in the panel is depth-robust.**")
A(f"- Eligible alternative combinations show a lower observed uncovered fraction in "
  f"**{int((M.anchor_pair_uncovered > M.lowest_observed_uncovered_eligible_pair_value).sum())} "
  f"of {len(M)} patients**.")
A("")
A("> **This is a detection-rate artifact, not a therapeutic finding.** The eligible "
  "alternatives are detected 1.8–2.8× more often than `GPRC5D`. A target that reads zero "
  "most of the time for depth reasons cannot add apparent coverage — which is exactly why "
  "`GPRC5D` failed the comparative QC gate **before** any gain was computed.")
A("")
A("Two distinctions that must not blur: `GPRC5D` **remains the frozen Stage-08 anchor** and "
  "fails only this supplemental comparative criterion; `TNFRSF17` carries the **identical** "
  "selection-dependence as `SDC1` and is retained only because it is the anchor — a "
  "disclosure, not a scientific distinction. **Normal marrow is expression context, never "
  "safety**: GPRC5D's decisive keratinized-tissue liability is structurally unobservable in "
  "a marrow dataset.")
A("")
A("**The only permitted wording is *greatest observed transcript-level malignant-cell "
  "coverage among evaluated combinations*.** No combination is optimal, recommended or best.")
A("")
A("## 8. What remains unresolved")
A("")
A("**The evidence axes are substantially discordant, and that is the central result.**")
A("")
A(f"- Only **{n_occupied} of 18** possible (tier × Level-1 × Level-2) cells are occupied.")
A("- **Level-1 and Level-2 evaluability are perfectly coupled**: the same 5 patients are "
  "`NOT_EVALUABLE` on both, and zero patients are non-evaluable on exactly one.")
A(f"- **The measurement-robust-high and Level-1-supported sets are disjoint** "
  f"(overlap = {len(rh & l1s_)}). All {n_rh} `robust-high` patients are Level-1 "
  "`NOT_SUPPORTED`; all {0} Level-1-supported patients are measurement-`uncertain`."
  .replace("{0}", str(n_l1s)))
A(f"- **No patient is simultaneously measurement-robust-high, Level-1 supported and "
  f"Level-2 supported** (n = {len(converge)}).")
A("")
A("The occupied profiles:")
A("")
A("| provisional tier | Level 1 | Level 2 | n |")
A("|---|---|---|---:|")
for _, r in joint.iterrows():
    A(f"| {r.provisional_measurement_tier} | {r.L1} | {r.L2} | {int(r.n)} |")
A("")
A("Uncertainty is reported in **three separate registers — measurement, biological and "
  "external-validity — and is never combined into one score** "
  "(`stage12_uncertainty_register.csv`). External-validity limits in particular are not "
  "addressable by any analysis of this dataset: transcript versus surface protein, untreated "
  "baseline versus post-therapy selection, marrow versus whole body, and cohort/protocol "
  "confounding in which depth, site, chemistry and deposit censoring all move together.")
A("")
A("Seven patients contribute repeated samples to the primary denominator "
  "(`stage12_repeated_patient_summary.csv`), reported as **ranges, never means**. "
  "`27522`'s apparent DN trajectory across six timepoints rests on samples of 1, 10, 7 and "
  "45 cells, so per-sample values are not interpretable in isolation. `60359` has **zero "
  "primary-denominator cells** and is excluded by that fact rather than by a filter.")
A("")
A("## 9. Claims supported")
A("")
A("Generated from `stage12_claim_ladder.csv`.")
A("")
A("| claim | status | permitted wording |")
A("|---|---|---|")
for _, r in CL.iterrows():
    if r.current_status in ("STRONG", "SUPPORTED_WITH_CAVEATS", "EXPLORATORY"):
        A(f"| {r.claim} | `{r.current_status}` | {r.allowed_wording} |")
A("")
A("## 10. Claims explicitly rejected")
A("")
A("| claim | status | prohibited wording |")
A("|---|---|---|")
for _, r in CL.iterrows():
    if r.current_status in ("NOT_SUPPORTED", "NOT_EVALUABLE"):
        A(f"| {r.claim} | `{r.current_status}` | {r.prohibited_wording} |")
A("")
A("**No patient synthesis category was created.** The design concluded *evidence matrix "
  "only*, after inspecting the joint cross-tabulation: with 5 of 18 cells occupied, a "
  "singleton at `MMY67868`, Level-2 supported in 26 of 27 evaluable patients, and the "
  "measurement and structure axes disjoint, any category scheme would either relabel three "
  "frozen columns losslessly or merge cells arbitrarily. **No composite score, no ranking, "
  "and no risk classifier exists anywhere in Stage 12.**")
A("")
A("---")
A("")
A("*Stage 12 consumed 29 frozen artifacts, wrote only into "
  "`results/12_final_synthesis/`, and changed no upstream scientific result.*")

(OUT / "stage12_summary.md").write_text("\n".join(lines) + "\n")
print(f"wrote stage12_summary.md ({len(lines)} lines)")

# %% [markdown]
# ## Step 10 — Re-verify the upstream freeze after execution
#
# Stage 12 must not have touched anything upstream. Re-hash all 393 manifested artifacts.

# %%
changed = []
for path, row in MANIFEST_ROWS.items():
    p = REPO / path
    if not p.exists():
        changed.append((path, "MISSING")); continue
    if sha256_of(p) != row["sha256"]:
        changed.append((path, "HASH CHANGED"))
if changed:
    for c in changed:
        print("UPSTREAM MUTATION:", c)
    raise SystemExit("STOP: Stage-12 integrity failure — an upstream artifact changed.")
print(f"upstream freeze re-verified: {len(MANIFEST_ROWS)}/{len(MANIFEST_ROWS)} unchanged")

written = sorted(p.relative_to(REPO).as_posix()
                 for p in OUT.rglob("*") if p.is_file())
print(f"\nStage-12 outputs written ({len(written)}):")
for w in written:
    print("  ", w)
print("\nSTAGE 12 EXECUTED AS SYNTHESIS ONLY")
