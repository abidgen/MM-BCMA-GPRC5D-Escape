# RESUME HERE

**Last session: 2026-08-26. Stage 08c and the Stage-11 LIANA verification arm (11b) are both
COMPLETE / ACCEPTED / FROZEN. Stage 8 is closed. Stop point is below.**

## Current stopping point

> **Stages 01–11 are complete through the Stage-11 LIANA verification arm. The supplemental
> multi-antigen coverage analysis is complete and frozen. Stage 12 has not started.**
>
> **Before Stage 12, the next planned work is an independent whole-codebase audit with Codex
> and independent project-process reconstruction documents from Claude Code and Codex.**
>
> **Do not perform additional biological analyses before that checkpoint unless a concrete
> reproducibility defect is discovered.**

## Freeze status — all sources agree

| stage | status |
|---|---|
| 01–05b · 06 · 07 | **FROZEN** |
| **08 core** dual-antigen escape | **FROZEN** — do not reopen |
| **08c supplemental** multi-antigen coverage | **COMPLETE / ACCEPTED / FROZEN** — Stage 8 closed |
| 09 bulk validation · 09b provisional tiers | **FROZEN** |
| 10 DN coherence | **FROZEN** |
| **11 custom immune context** | **FROZEN** |
| **11b LIANA verification arm** | **COMPLETE / ACCEPTED / FROZEN** |
| **12 decision packet** | **NOT STARTED** |

## Where the project stands

**Stages 01–11 are complete and FROZEN. The supplemental 08c coverage deliverable is
COMPLETE / ACCEPTED / FROZEN, and the 11b LIANA verification arm is complete.**
**507 pass, 1 skip** (`pytest -m "not slow"`, `mm-core`).

Nothing is staged or committed.

## What this session did

### 1. Stage 11 — resumed, corrected, accepted, FROZEN

Stage 11 had tables on disk but **no driver**. It now has one
(`notebooks/11_immune_context.py`), which recomputes everything from frozen upstream and
asserts agreement with the preserved first run in `results/11_immune_context/preliminary_run/`
— max drift **1.78e-15**. Three dated amendments in `stage11_design.md`: the receiver
definition was corrected back to the frozen design (both versions kept); the resulting
"stronger" result was traced to a **receiver-side confound** (11 of 15 receptors move down
with DN burden); and the depth covariate was corrected to the intersected gene space.

> **No immune association survives correction. `Tcell PDCD1 → CD274` is NOT an
> immune-evasion axis** and must never be written up as one.

### 2. Stage 08c — supplemental multi-antigen coverage — **FROZEN**

`notebooks/08c_multi_antigen_coverage.py` → `results/08_dual_antigen_escape/multi_antigen_coverage/`.
Design frozen **before** any pair or triple was computed. **No frozen artifact was touched** —
verified by timestamp and by `frozen_upstream_digests.json`, which `tests/test_coverage.py`
re-checks on every run.

> **Measurement quality, not biology, is the binding constraint on this panel.**

- **No target is depth-robust.** Detection-vs-depth ρ = 0.32–0.48 for all seven.
- **`GPRC5D` is `COVERAGE_NOT_EVALUABLE`** (technical-zero 0.62 ≥ 0.50) — the frozen
  Stage-08 conclusion reached again by an independent rule. **Threshold not relaxed.**
- **`SDC1` is `COVERAGE_NOT_EVALUABLE` on circularity**, not differentiation.
  `config.PLASMA_MATURE = ("SDC1", "TNFRSF17")` is Stage 06's axis-(b) plasma predicate.
  **`TNFRSF17` carries the identical limitation and is kept only as the frozen anchor —
  a disclosure, not a distinction.**
- **Anchor:** BCMA uncovered 0.353, GPRC5D 0.899, pair 0.335, **GPRC5D gain 0.011**.
  **Never read as GPRC5D being clinically redundant.**
- **Alternative pairs beat the anchor in 32/32 patients — a detection-rate artifact.**
  No combination is optimal, recommended or best.
- **All 21 pairs collapse to ~1.0 co-loss under depth conditioning.**
- Truncate-10k: WashU exactly unchanged, MMRF +0.05–0.07, ordering ρ 0.996.
- Two patients (`27522` 0.571, `59114` 0.547) vary more within themselves than most
  patients differ from each other.

### 3. Stage 11b — LIANA verification arm — **FROZEN**

`notebooks/11b_liana_verification.py` → `results/11_immune_context/liana_verification/`.
**Exploratory · post hoc · non-tier-changing · non-classifying.** The frozen Stage-11 custom
analysis was **not** reopened, replaced or rewritten.

`liana 1.8.1`; `rank_aggregate` RRA consensus (5 methods) **plus** `cellchat` for continuity,
both on the `consensus` resource held fixed; per-patient via LIANA's own `Method.by_sample`;
`min_cells` 20 (frozen `MIN_SENDER_CELLS`); frozen Stage-11 confound model unchanged.

- **31 of 32 patients evaluable** — `25183` has zero sender cells, matching frozen Stage 11.
- 1,050 interactions/patient; **87 tested at ≥20 patients; 12 raw p<0.05; 1 at BH<0.10.**
- **The one consensus BH hit is `Myeloid TNFSF13B → TNFRSF17` — structurally circular.**
  Its receptor is `TNFRSF17`, half of what defines `obs_dn_primary`, so the negative
  coefficient is arithmetic. Reproduced by CellChat and the all-plasma receiver.
  **100% of consensus BH hits are antigen-circular.**
- **3 `RECEIVER_STATE_CONFOUNDED`, 9 `NOT_REPRODUCED_BY_LIANA`, 0 `EXPLORATORY_LIANA_ONLY`,
  2 `ABUNDANCE_SENSITIVE`, 0 `CONSISTENT_WITH_TARGETED_PANEL`.**
- **`PDCD1 → CD274` is not in LIANA's resource → `NOT_EVALUABLE`.** LIANA can neither
  reproduce nor contradict it; Stage 11 had already classified it as receiver-side confounded.
- Only **8 of 17** frozen custom pairs are in LIANA's resource, and only one row cleared the
  patient floor — LIANA's `expr_prop` filter removes the sparse TRAIL/granzyme interactions.
- **No tier, state, composition conclusion or coverage eligibility changed. No classifier.**

**Accepted interpretation (final wording, not to be strengthened):** Stage 11 found no robust
independent evidence that immune composition or ligand–receptor communication explains the
observed DN phenotype. The targeted LR analysis was receiver-state confounded, and LIANA
verification did not rescue that interpretation. The strongest LIANA consensus association was
structurally circular because its receptor was `TNFRSF17`, one of the antigens defining the DN
predictor. **LIANA is not a validation of immune evasion.**

Kept explicitly separate from the custom Stage-11 analysis, which **was not performed with
LIANA**. `25183` is a **biological/evaluability exclusion, not a software failure**.
`PDCD1 → CD274` = `NOT_EVALUABLE_BY_LIANA_RESOURCE` — **neither disproved nor supported**.
LIANA is a **partial methodological verification arm, not a full reproduction** of the 17-pair
panel; missing resource coverage is **not** a negative biological result.


## Frozen state, by axis — never combined into a scalar score

| axis | states | where |
|---|---|---|
| measurement (08/09b) | 4 robust-high · 28 uncertain · **0 robust-low** | `results/08_dual_antigen_escape/risk_tier_provisional/` |
| Level-1 DN structure | 4 · 23 · 5 | `results/10_dn_coherence/` |
| Level-2 DN phenotype | 26 · 1 · 5 | `results/10_dn_coherence/` |
| Level-3 genomic | 32 `CNV_SUBCLONE_NOT_EVALUABLE` | — |
| immune context (11) | **nothing survives correction** | `results/11_immune_context/` |
| multi-antigen coverage (08c) | 5 eligible · **2 not evaluable** | `.../multi_antigen_coverage/` |
| LIANA verification (11b) | **no credible LIANA-only interaction** | `results/11_immune_context/liana_verification/` |

**No patient is simultaneously measurement-robust-high, Level-1 supported and Level-2
supported.** `DN_STATE_SUPPORTED` is weakly discriminative (26/27). Measurement tiers are
**provisional** — measurement-robust only, never a biological classification.

## Stage 12 — what it consumes and what it may not do

Inputs, as **separate columns, never a composite**:
- `stage12_multi_antigen_interface.csv` (08c) — one row per patient.
- Provisional measurement tiers (09b), DN coherence levels (10), immune context (11),
  the bias-direction table, matched-bulk validation (09).

Hard constraints carried forward:
- **Risk tiers, never a rank ordering.** Caterpillar plot with CIs, never a bar chart.
- **No composite risk score anywhere**, and no `coverage − λ · exposure` utility.
- **No combination may be called optimal, recommended or best.** Permitted wording:
  *greatest observed transcript-level malignant-cell coverage among evaluated combinations*.
- **Normal marrow is expression context, never safety.** GPRC5D's keratinized-tissue
  liability is unobservable here and stays a cited external caveat.
- **The mRNA-vs-protein limitation is stated explicitly and mechanistically.**
- Decision rules stated in advance, not fitted after seeing the ranking.

## Reproducibility status

- **`.py` is authoritative** for every notebook; `.ipynb` generated explicitly with
  `jupytext --to notebook`. **`jupytext --sync` was never used.**
- `notebooks/08c_multi_antigen_coverage.py` and `notebooks/11b_liana_verification.py` were
  both **executed end-to-end** via `nbclient` from the repo root.
- **Do not rerun the expensive LIANA computations** (~20 min) unless verifying a concrete
  defect.

## Open items

- **Codex whole-codebase audit + independent process-reconstruction documents — next.**
- **Stage 12 decision packet — after that checkpoint.** Everything upstream is frozen.
- Phase 2 (GSE117156) strictly after Phase 1. Never merged with GSE223060 (MARS-seq vs 10x);
  separate `phase2_`-prefixed pipeline.

## Rules that bit during these sessions

- **Never `jupytext --sync`** — it treated the `.ipynb` as authoritative and silently dropped
  cells appended to the `.py`, three times. Use
  `jupytext --to notebook notebooks/NN_x.py -o notebooks/NN_x.ipynb`, then execute and confirm
  the expected numbers appear.
- **`jupyter nbconvert --execute` chdir's to the notebook's directory**, breaking every
  relative path. Execute with `nbclient` and `resources={"metadata": {"path": "."}}` from the
  repo root. (A stray empty `notebooks/results/` tree from an earlier session was removed.)
- **`obs["total_counts"]` is not the depth Stage 08 normalised against** — it is a QC-time sum
  over each sample's full Cell Ranger reference, taken before the stage-05 intersection to
  32,991 genes. Recompute the row sum on the intersected space
  (`communication.stream_gene_counts` returns it).
- Read a CSV of patient IDs with `dtype={"patient": str}` — silent int64 coercion produced an
  empty join in notebook 10.
- **Stage 08's ambient reference is `cell_type` ∈ {Tcell, Myeloid, Bcell, HSPC}**, established
  by reproducing `noise_floor_ambient.csv` exactly. Its predeclaration also names `NK`, which
  resolved to nothing because Stage 06 emits no NK class.
