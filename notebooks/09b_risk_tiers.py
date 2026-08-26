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
# # 09b — Provisional risk tiers (measurement-robust only)
#
# > **PROVISIONAL — Stage-08/09b measurement-robust classification only. Stage 10
# > biological coherence assessment is required before final biological interpretation.**
#
# **`robust-high` means the high-DN conclusion survives not just the technical sensitivity
# analyses, but also the reasonable uncertainty over what numerical DN fraction should
# count as "high."** It does **not** mean the DN population is biologically coherent —
# that is Stage 10's question, and it sits between this provisional tiering and any final
# biological risk classification.
#
# The rule is frozen in `results/08_dual_antigen_escape/risk_tier_provisional/risk_tier_policy.md`
# and was written before it was applied to any patient. This notebook applies it; it does
# not define it, and it does not modify Stage 08 or Stage 09.
#
# It takes a **letter** (`09b`) for the same reason `05b` did: number order is execution
# order with no exceptions, and this is a terminal step off Stage 08's evidence rather than
# a new stage in the pipeline.

# %%
import numpy as np
import pandas as pd
from pathlib import Path

from mm_escape import risk_tiers as RT

RT_DIR = Path("results/08_dual_antigen_escape/risk_tier_provisional")
EV = pd.read_csv("results/08_dual_antigen_escape/risk_tier_design/patient_evidence_matrix.csv")
pd.set_option("display.width", 300)

print(f"TAU_LOW={RT.TAU_LOW}  TAU_HIGH_SET={RT.TAU_HIGH_SET}  "
      f"min n={RT.MIN_N_PRIMARY}  min cells/sample={RT.MIN_CELLS_PER_SAMPLE}")

# %% [markdown]
# ## Thresholds are conventions, not biology
#
# `TAU_LOW = 0.05` was predeclared in `stage08_predeclaration.md` §9 as the smallest DN
# fraction the project claims to resolve. It is **not** relaxed because `robust-low` may
# turn out empty — an empty tier is an acceptable scientific outcome.
#
# `TAU_HIGH ∈ {0.20, 0.25, 0.33}` are **transparent project conventions and sensitivity
# thresholds, not clinically validated boundaries**. Stage 09 found no biologically
# validated joint-DN threshold, and no cutoff is supported by bulk. They were not derived
# from natural breaks, and patient membership was not inspected before the policy froze.

# %%
def parse_per_sample(s):
    if not isinstance(s, str) or not s:
        return []
    return [(a, int(n), float(v)) for a, n, v in (tok.split(":") for tok in s.split("|"))]


records = []
for _, r in EV.iterrows():
    # Only Stage-07/08 evidence is assembled. No cohort, no bulk, no p-value.
    ev = dict(n_primary=int(r.n_primary), dn_primary=r.obs_dn_primary,
              dn_sensitivity=r.obs_dn_sensitivity, dn_trunc10k=r.obs_dn_trunc10k_primary,
              enr_ci_lo=r.enr_lo_primary, enr_cohortbins=r.enr_primary,
              enr_globalbins=r.enr_globalbins_primary,
              per_sample=parse_per_sample(r.per_sample_dn_primary))
    final, robustness, tiers, reasons = RT.final_tier(ev)
    records.append({"patient": r.patient, "final_tier": final,
                    "threshold_robustness": robustness,
                    **{f"tier_tau{int(t*100):03d}": v for t, v in tiers.items()},
                    **reasons})
T = pd.DataFrame(records)
print(T.final_tier.value_counts().to_string())
print(T.threshold_robustness.value_counts().to_string())

# %% [markdown]
# ## The same rule at all three thresholds
#
# Only `TAU_HIGH` varies between runs; every other criterion is identical. `robust-high`
# requires qualifying under **all three** — qualifying at 0.20 alone is not enough, so the
# high-risk reading cannot depend on whether "meaningfully elevated" means a fifth, a
# quarter or a third of the clone.

# %%
print(pd.read_csv(RT_DIR / "risk_tier_threshold_comparison.csv").to_string(index=False))

# %% [markdown]
# **All 32 patients are `THRESHOLD_ROBUST`** — not one changes tier across the set, so
# `uncertain_threshold` never fires. The binding constraints are dropout compatibility and
# intermediate DN, not the cutoff.

# %% [markdown]
# ## Significance did not decide the tiers
#
# The four patients with significant depth-conditioned enrichment and the four
# `robust-high` patients **overlap in only one**. Three `robust-high` patients have
# non-significant permutation p-values; three significant patients are `uncertain`, blocked
# by low DN, intermediate DN or denominator disagreement. A significance gate would have
# encoded sequencing power as biology, since the significant patients all sit in the
# deepest cohort.

# %%
F = pd.read_csv(RT_DIR / "risk_tiers_provisional.csv")
print(F[["patient", "cohort", "n_primary", "dn_primary", "enr_ci_lo",
         "perm_p_reported_not_used", "final_tier", "uncertainty_reasons"]]
      .round(4).to_string(index=False))

# %% [markdown]
# ## The limitation that matters most
#
# The enrichment CI lower bounds for all four `robust-high` patients are **1.0077, 1.0062,
# 1.0056 and 1.0009** — membership hinges on a bootstrap CI lower bound in the *third
# decimal place*, and `MMY74196` would flip on bootstrap noise alone. This is stated rather
# than fixed: re-specifying the criterion after seeing which patients it admits is exactly
# the post-hoc tuning the design exists to prevent. **Read the tier as provisional.**

# %% [markdown]
# ## Cohort diagnostic — after assignment, and diagnostic only

# %%
print(pd.read_csv(RT_DIR / "risk_tier_by_cohort_diagnostic.csv").to_string(index=False))

# %% [markdown]
# Skew remains and points **opposite** to the one that was feared: with significance
# blocked, the deepest cohort (MMRF) yields the fewest `robust-high` calls and the
# shallowest-detecting (WU2) the most — because MMRF is highest on every instability
# dimension and its patients are routed to `uncertain` by robustness failures.
#
# Since Stage 09 showed WU2-type GPRC5D negativity is the most dropout-compatible in the
# cohort, this concentrates `robust-high` where the negative calls are least reliable.
# Reported as a limitation. **No tier was adjusted on the basis of it**; any methodological
# response needs a separate reviewed experiment.
#
# The data-availability channel runs the other way too: `robust-high` patients have *less*
# data (0.000 multi-sample, 0.250 matched bulk) than `uncertain` ones (0.214, 0.464). So
# **all four passed the repeated-sample criterion vacuously** — their within-patient
# stability is unobserved, not demonstrated, which is why the rule never credits a
# single-sample patient with "agreeing".

# %%
print((RT_DIR / "risk_tier_provisional_summary.md").read_text())
