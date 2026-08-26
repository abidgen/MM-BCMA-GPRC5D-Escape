# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#   kernelspec:
#     display_name: mm-core
#     language: python
#     name: mm-core
# ---

# %% [markdown]
# # Stage 09 — orthogonal **marginal** antigen validation
#
# Do independent bulk RNA-seq (GSE223061) and normal marrow support the patient-level
# marginal `TNFRSF17` and `GPRC5D` abundance seen in single-cell data?
#
# **Hard invariant.** Bulk RNA can validate marginal antigen abundance. It **cannot**
# determine whether the same individual cells are simultaneously BCMA-negative and
# GPRC5D-negative — a tumour that is 50% BCMA⁺GPRC5D⁻ plus 50% BCMA⁻GPRC5D⁺ shows healthy
# bulk expression of *both* genes while containing zero dual-positive cells. Bulk destroys
# the joint distribution over cells by construction, so nothing here validates
# `frac_double_negative`, co-negativity enrichment, or any cell-level escape population.
#
# **Stage 08 is frozen.** Nothing in this notebook modifies, recalibrates or reruns any
# Stage-08 measurement. Where Stage 09 disagrees with Stage 08, the disagreement is
# reported — Stage 08 is not changed to make them agree. Guardrails:
# `results/09_bulk_validation/stage09_design.md`, written before any bulk file was opened.

# %%
import numpy as np
import pandas as pd
from pathlib import Path

from mm_escape import bulk as BK

OUT = Path("results/09_bulk_validation")
RAW = Path("raw/unpacked_bulk")
pd.set_option("display.width", 260)

# %% [markdown]
# ## The deposit holds two different assays in two different formats
#
# WashU-1 is kallisto **transcript**-level TPM inside per-sample tarballs; MMRF is
# **gene**-level TPM on Ensembl IDs. They are never pooled: **MMRF is CD138⁺ sorted** and
# pairs with malignant-cell pseudobulk, **WashU 1 is unsorted BMMC** and pairs with
# whole-sample pseudobulk. Pooling would make a large share of comparisons measure tumour
# burden instead of antigen abundance.

# %%
t2g = BK.antigen_transcript_map("raw/gtf/Homo_sapiens.GRCh38.93.gtf.gz")
print(f"antigen transcripts in Ensembl 93: {len(t2g)}")

records = []
for f in sorted(RAW.glob("GSM*")):
    gsm, rest = f.name.split("_", 1)
    if f.name.endswith(".tar.gz"):
        records.append({"gsm": gsm, "bulk_sample": rest[:-7], "bulk_cohort": "WU1",
                        "specimen": "unsorted BMMC", "usable": True,
                        **BK.read_washu_bulk(f, t2g)})
    else:
        records.append({"gsm": gsm, "bulk_sample": rest[:-len("_tpm.tsv.gz")],
                        "bulk_cohort": "MMRF", "specimen": "CD138+ sorted",
                        **BK.read_mmrf_bulk(f)})
B = pd.DataFrame(records)
print(f"bulk files {len(B)}; usable {int(B.usable.sum())}; "
      f"empty stubs {sorted(B[~B.usable].bulk_sample)}")

# %% [markdown]
# Three forensic findings, each of which would have corrupted the analysis silently:
#
# * **2 header-only 114-byte stubs** (`MMRF_1505`, `MMRF_2259`) → `NOT_EVALUABLE`.
# * **`MMRF_1686` stacks two sequencing runs in one file** (80,590 rows, TPM summing to
#   2×10⁶). TPM is already per-run normalised, so `read_mmrf_bulk` **averages** the runs —
#   summing would have inflated its `TNFRSF17` from 145 to 289.
# * **Three bulk IDs have no scRNA counterpart**: `47499` (scRNA has `47491_1`/`47491_2`),
#   `59114_2` (scRNA has `_1` and `_4`), and `98433`. **No near-miss is resolved by
#   inference** — every one is `NOT_EVALUABLE`.

# %%
print(B[B.n_runs > 1][["bulk_sample", "n_runs", "TNFRSF17_tpm", "GPRC5D_tpm"]].to_string(index=False))

# %% [markdown]
# ## One observation per biological patient
#
# The selection rule was frozen before any correlation ran: **the earliest matched
# timepoint clearing the pre-existing ≥20-cell floor**, never the timepoint that agrees
# better with bulk, and multiple samples from one patient are never counted as independent
# patients. `56203_2` (17 cells) falls out; `27522` contributes its first timepoint only.

# %%
mapping = pd.read_csv(OUT / "sc_bulk_sample_mapping.csv")
print(mapping.match_status.value_counts().to_string())
matched = mapping[mapping.match_status == BK.MATCHED_EXACT]
print(f"\nexact matches {len(matched)} over {matched.sc_patient.nunique()} patients; "
      f"one-to-many: {sorted(matched.sc_patient.value_counts()[lambda x: x > 1].index)}")

# %% [markdown]
# ## Frozen single-cell metric
#
# **Primary: pseudobulk CPM over the frozen denominator cells** — it measures *abundance*,
# the quantity bulk also measures, aggregates sparsity rather than thresholding it, and
# needs no cutoff. **Secondary, descriptive: detection fraction.** No third summary was
# computed, and the primary metric was not reselected after seeing which correlates better.
# Every comparison runs independently under **both** frozen Stage-07 denominators.

# %%
R = pd.read_csv(OUT / "bulk_vs_sc_by_cohort.csv")
print(R[["comparison", "n", "rho", "p", "status"]].to_string(index=False))

# %% [markdown]
# **Where the assay pairing is correct — CD138⁺-sorted bulk against malignant-cell
# pseudobulk — both antigens reproduce at ρ = 0.933 (n = 9).** That is the strongest
# external support this project has for its marginal antigen measurements.
#
# WU1's GPRC5D ρ = −0.60 is **not** evidence of an inverse relationship: n = 5, p = 0.28,
# and the pairing is deliberately mismatched (unsorted bulk vs sorted-clone pseudobulk).
# Under the composition-appropriate whole-sample pairing it is ρ = 0.10 — uninformative
# rather than inverse. Both WU1 results are underpowered and support no conclusion.
# **Pooled rows mix two assays** and are reported for completeness only.

# %% [markdown]
# ## The WashU GPRC5D question
#
# WashU bulk carries **more** GPRC5D than MMRF bulk (median 116.2 vs 74.4 TPM) while WashU
# single-cell detects it at **one quarter** the rate (0.074 vs 0.298). The sharpest case is
# `77570`: bulk GPRC5D 116.6 — comparable to the highest MMRF values — against sc detection
# 0.024 and 3.9 CPM.
#
# This is **single-cell dropout / capture limitation strongly implicated**, independently
# corroborating Stage 08's technical-zero floor. **It changes nothing**: not the floor, not
# the GPRC5D zero calls, not the depth-conditioned null, not the DN enrichment.
#
# The converse also appears. `59114` has near-zero bulk for *both* antigens with normal
# single-cell CPM and only 22 denominator cells — the tumour-purity caveat, not a
# single-cell failure. **The deposit carries no purity column**, so disagreement may reflect
# specimen composition as well as dropout, and not every disagreement is a technical failure.

# %%
print(pd.read_csv(OUT / "bulk_vs_sc_truncate10k_sensitivity.csv").to_string(index=False))

# %% [markdown]
# Equalising the depth ceiling barely moves bulk agreement. **A validation sensitivity, not
# a reason to change Stage 08** — the frozen truncation procedure was reused unchanged.

# %% [markdown]
# ## Normal marrow — donor is the biological unit
#
# 647 donor plasma cells across 8 donors are **never treated as 647 replicates**.

# %%
N = pd.read_csv(OUT / "normal_marrow_antigen_context.csv")
print(N[N.population == "PlasmaCell"].to_string(index=False))

# %% [markdown]
# Normal marrow plasma cells carry substantial BCMA (detection median 0.377) and **almost
# no detectable GPRC5D** (median 0.009, and 0.000 in every donor's T and myeloid cells);
# malignant cells exceed both (0.768 / 0.154). Across donors GPRC5D detection tracks depth
# at Spearman 0.756 — the same depth dependence seen everywhere else in this project.
#
# **This is marrow expression context only.** It cannot address whole-body safety,
# on-target/off-tumour toxicity in other organs, malignant clone identity, or any patient's
# joint-DN probability. GPRC5D's decisive off-tumour liability is keratinized tissue, which
# a marrow dataset cannot observe at all.

# %%
print((OUT / "stage09_bulk_validation_summary.md").read_text())
