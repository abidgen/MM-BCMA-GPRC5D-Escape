# Marker Panels and Malignant-Cell Identification — Reference

**MM Dual-Antigen BCMA/GPRC5D project · stages 06 and 07**

Q&A backup for the J&J seminar. Everything here is read out of
`projects/mm-dual-antigen/src/mm_escape/config.py`, `malignant.py`, and the frozen
stage-06/07 results — nothing is recalled from memory. Line references are to
`config.py` unless stated.

Two separate questions, answered by two separate mechanisms:

| Question | Stage | Mechanism |
|---|---|---|
| What cell type is this? | 06 | Marker panels + three-method adjudication |
| Is this plasma cell malignant? | 07 | Light-chain restriction, refined by V-gene usage |

---

# Part A — Cell-type annotation (stage 06)

## A1. The seven-class marker panel

`MARKER_PANEL` (line 248). Seven classes; no "Cycling" or "Activated" class exists —
those are states, not identities, and are carried as continuous scores.

| Class | Markers |
|---|---|
| **PlasmaCell** | `SDC1`, `CD38`, `MZB1`, `XBP1`, `IRF4` |
| **Bcell** | `MS4A1`, `CD79A`, `CD19` |
| **Tcell** | `CD3D`, `CD3E`, `CD3G`, `TRAC`, `TRBC1`, `TRBC2` |
| **NK** | `NCAM1`, `NKG7`, `GNLY`, `KLRD1`, `KLRF1` |
| **Myeloid** | `CD14`, `LYZ`, `ITGAM` |
| **Erythroid** | `GYPA`, `AHSP`, `ALAS2`, `CA1`, `HBA1`, `HBA2` |
| **HSPC** | `CD34`, `KIT` |

Unresolved clusters stay `Ambiguous` (`AMBIGUOUS_LABEL`) rather than being forced into
a class.

## A2. The plasma-cell call needs two independent axes

This is the part worth knowing, because PlasmaCell sets the denominator of the whole
project's headline metric.

A secretory program alone is **not** sufficient. Cluster 24 exposed why: a lymphoid
progenitor expressing `MZB1`/`XBP1` was pulled toward PlasmaCell by secretory genes
alone. So a PlasmaCell claim requires concordant evidence on two axes (line 547):

| Axis | Constant | Genes | Rule |
|---|---|---|---|
| (a) secretory program | `PLASMA_SECRETORY` | `MZB1`, `XBP1` | **all of** |
| (b) mature plasma identity | `PLASMA_MATURE` | `SDC1`, `TNFRSF17` | **at least one of** |

**`CD38` is deliberately excluded from axis (b)** — it is expressed on activated T cells,
NK cells, pro-B cells and progenitors, and is not plasma-specific. It stays in
`MARKER_PANEL` for descriptive scoring but cannot satisfy the predicate alone. `MZB1`/
`XBP1` alone — axis (a) without (b) — also cannot.

> **Known circularity, disclosed in stage 08c:** both `PLASMA_MATURE` genes are also
> antigens the project measures. `SDC1` was ruled `COVERAGE_NOT_EVALUABLE` on exactly
> this basis. `TNFRSF17` carries the identical selection-dependence and is retained only
> because it is the frozen anchor — **a disclosure, not a scientific distinction**. The
> predicate operated at cluster level and no cell was removed for its own `TNFRSF17` zero.

## A3. Contradiction programs are a *different* gene list

`LINEAGE_PROGRAMS` (line 733) is used **only** to detect that a cell contradicts its
assigned class. It is deliberately not `MARKER_PANEL`:

> A panel used to *identify* a class can afford an ambient-prone gene among several. A
> program used to *accuse* a cell of another lineage cannot — a false accusation is
> manufactured out of ambient RNA, whereas a missed one is merely hidden by dropout.

| Lineage | Contradiction program | Excluded, and why |
|---|---|---|
| T | `CD3D`, `CD3E`, `CD3G`, `TRAC`, `TRBC1`, `TRBC2` | — |
| B | `MS4A1`, `CD79A`, `CD79B`, `CD19` | — |
| erythroid | `GYPA`, `AHSP`, `ALAS2`, `CA1` | **globins** (`HBB`/`HBA1`/`HBA2`) — dominant ambient species |
| myeloid | `CD14`, `FCN1`, `MNDA`, `ITGAM` | **`LYZ`** — dominant ambient species |

There is **no plasma/B contradiction program at all**, because immunoglobulin is the
third big ambient species in marrow.

`CONTRADICTION_PAIRS` (line 749) — which lineages are incompatible with each class:

| Class | Contradicted by |
|---|---|
| PlasmaCell | T, erythroid, myeloid |
| Bcell | T, erythroid, myeloid |
| Tcell | B, erythroid, myeloid |
| NK | T, B, erythroid, myeloid |
| Myeloid | T, B, erythroid |
| Erythroid | T, B, myeloid |
| HSPC | *(none)* |

**What is absent is as load-bearing as what is present.** PlasmaCell is not contradicted
by B and vice versa — plasma cells *are* B-lineage and the plasmablast continuum is real
biology. HSPC has no contradictions at all: progenitors legitimately co-express
lineage-priming programs, and flagging that would be a biology error, not a QC finding.

**Contradiction is detected, never inferred from absence.** Dropout can only hide
evidence, so a detection-based rule under-calls and can never manufacture a contradiction
from a zero.

## A4. Identity versus state — the T/NK distinction

The original annotation over-called NK from `NKG7`/`GNLY` alone; 22,132 cells were
corrected. The fix was to separate identity anchors from shared effector state:

| Constant | Genes | Status |
|---|---|---|
| `T_IDENTITY_ANCHORS` (701) | `CD3D`, `CD3E`, `CD3G`, `TRAC` | T-lineage **identity** |
| `T_CONTEXT` (704) | `TRBC1`, `TRBC2` | measured and reported, **never sufficient alone** |
| `NK_IDENTITY` (707) | `KLRD1`, `KLRF1`, `NCAM1`, `FCGR3A`, `KLRC1` | NK identity |
| `GD_IDENTITY` (712) | `TRDC`, `TRGC1`, `TRGC2` | γδ axis |
| `CYTOTOXIC_STATE` (716) | `NKG7`, `GNLY`, `PRF1`, `GZMB`, `GZMA`, `CTSW` | **state — establishes no lineage** |

Frozen operational conclusion, quoted verbatim:

> Isolated `TRBC1`/`TRBC2` expression is insufficient evidence of T-lineage commitment
> because it frequently occurs without coordinated `CD3`/`TRAC` expression in cells with
> strong NK-lineage evidence.

`TRDC` alone is likewise not promoted to identity: in Leiden 23, `TRGC` detection was
near-independent of `TRDC` (28.2% vs 21.6%), and 92.6% of `TRDC`+ cells carried strong NK
evidence.

## A5. How a class was actually decided

Three methods compared per class — **manual marker adjudication, `celltypist`, `SingleR`** —
against bars declared before any result was computed and identical across all three
revisions.

| Gate | Value | Note |
|---|---|---|
| Concordance (F1) — PlasmaCell | **0.95** | sets the stage-07 denominator |
| Concordance — T / NK / Myeloid | **0.90** | defines stage-08's ambient noise floor |
| Concordance — B / Erythroid / HSPC | **0.85** | nothing downstream is load-bearing |
| `MARKER_COVERAGE_MIN` | **0.30** | min mean scaled expression of a class's own markers |
| `CONTRADICTION_MIN_GENES` | 2 | genes detected to call a contradiction |
| `CONTRADICTION_MAX_RATE` | 0.25 | max contradicting fraction tolerated |

**F1 here is concordance, not accuracy.** Manual labels are a third opinion from the same
matrix, not ground truth. The biological evidence is the marker-coverage test, and it can
**veto** a class regardless of how well the three methods agree — high agreement on an
unsupported label is agreement on an error.

Evidence is pooled at **cluster** level, not per cell: at 1,162 median genes/cell, a
per-cell marker call on a dropped-out gene is a *wrong* call rather than a missing one.

## A6. Result

| Class | n | % |
|---|---:|---:|
| Tcell | 60,896 | 35.2 |
| PlasmaCell | 35,474 | 20.5 |
| Myeloid | 33,817 | 19.6 |
| Erythroid | 16,224 | 9.4 |
| Bcell | 12,226 | 7.1 |
| Ambiguous (Leiden 23) | 11,424 | 6.6 |
| HSPC | 2,879 | 1.7 |
| **Total** | **172,940** | |

**The plasma boundary agrees with the source paper's own annotation on 32,307 of 32,337
cells.** Downstream stages read `obs["cell_type"]` and nothing else.

---

# Part B — Malignant-cell identification (stage 07)

Malignancy is **not** called from a marker. Plasma cells are plasma cells; the question
is whether a given plasma cell belongs to the patient's dominant clone.

## B1. Axis 1 — light-chain restriction (ratio, never presence)

Genes (`REQUIRED_GENES["light_chain"]`, line 167):

```
kappa   IGKC
lambda  IGLC1  IGLC2  IGLC3  IGLC4  IGLC5  IGLC6  IGLC7
```

All `IGLC` members matter because the call is a ratio.

**Per-cell rule** (`malignant.light_chain_class`):

```
total = kappa_UMI + lambda_UMI
frac_k = kappa_UMI / max(total, 1)

total < LC_CLASS_MIN_UMI (3)            -> "insufficient"
frac_k      >= LC_CLASS_MIN_FRAC (0.80) -> "kappa"
(1 - frac_k) >= 0.80                    -> "lambda"
otherwise                               -> "ambiguous"
```

> Immunoglobulin is the most ambient-contaminated transcript family in this tissue, so a
> presence call is far noisier than it looks; a ratio is robust to a shared additive
> background.

**Per-patient clonality**, from the dominance statistic D, with cut points placed against
the observed donor ceiling (donor max **0.824**):

| State | Rule | Rationale |
|---|---|---|
| `CLONAL_STRONG` | D ≥ **0.95** | above every donor with wide margin; high-specificity because this sets the stage-08 denominator |
| `CLONAL_WEAK` | 0.85 ≤ D < 0.95 | above the donor ceiling but short of the specificity bar |
| `NO_RESTRICTION` | D < 0.85 | encloses all five qualifying donors including the 0.824 outlier |

## B2. Axis 2 — patient-specific V-gene usage

**This is a higher-specificity refinement of the same immunoglobulin clonality axis — not
an independent second axis.** J segments are uncapturable at 3′ (`IGKJ` 0.00%, `IGHJ`
0.04%), so there is no clonotype, only patient-specific V usage.

| Constant | Value | Meaning |
|---|---|---|
| `DOMINANT_V_MIN_FRAC` | **0.50** | top-V fraction required |
| `DOMINANT_V_MIN_ENRICHMENT` | **3.0** | fold enrichment over other patients |
| `V_EVALUABLE_MIN_CELLS` | 50 | V-positive cells |
| `V_EVALUABLE_MIN_PCT` | 0.50 | of plasma cells LCV-positive |
| `V_PARTIAL_MIN_CELLS` | 20 | |
| `V_PARTIAL_MIN_PCT` | 0.20 | |

The 0.50 threshold sits inside the **empirical gap** between assessable healthy donors
(max **0.378**) and evaluable disease patients (min **0.562**) — 1.32× the donor maximum,
and deliberately *not* placed at the disease minimum.

> Three donors are `V_NOT_EVALUABLE` on 12–47 plasma cells; one of them reads 0.625 on
> **12 cells**. They are excluded from the donor reference for that reason, not to make
> the gap look better.

## B3. The per-cell state machine

`malignant.clone_membership`. Order is deliberate.

1. **Positive incompatibility first** — a minority light-chain class, or a coherent
   alternative V at ≥2 UMI, is positive evidence *against* membership and outranks
   everything. **Absence of the dominant V is not in this list.**
2. **Support requires all four** — `CLONAL_STRONG` patient, `V_EVALUABLE` patient,
   dominant-class compatibility, and positive dominant-V detection.
3. **Compatible-but-unobserved is its own state**, never folded into support or
   incompatibility — it is exactly the population the data cannot settle.
4. **Everything else is uncertain.**

| State | Meaning |
|---|---|
| `CLONE_SUPPORTED` | **21,906 cells / 32 patients** — the primary denominator |
| `CLONE_COMPATIBLE_V_UNOBSERVED` | **7,109 cells** — clone-compatible, dominant V not observed |
| `CLONE_INCOMPATIBLE` | positive evidence against membership |
| `CLONE_UNCERTAIN` | everything else |

Primary denominator **21,906**; sensitivity denominator **29,015** (primary + V-unobserved).
Both are reported for every patient and **never collapsed or selected between**.

**V absence is never negative evidence.** V detectability is strongly depth-dependent:
supported-vs-unobserved median UMI ratio is 1.79× pooled but **17.6× in MMRF**. The
decisive evidence that this is technical rather than biological — patient `83942`/
`MMY83942` is one patient, one clone, one dominant V, and splits **0.746 vs 0.351** across
two protocols.

## B4. Axis 3 — attempted and rejected

CNV was attempted as a genuinely independent malignancy axis and **rejected before any
disease sample was inspected**, because it failed its healthy-donor negative control.

`infercnvpy`, window 100, step 10, same-sample T + myeloid reference, 2,928 windows, IG
loci and the antigens excluded:

| Donor | n plasma | plasma z median | frac z > 3 |
|---|---:|---:|---:|
| BM5 | 98 | −0.89 | 0.0% |
| BM2 | 149 | +0.78 | 3.4% |
| BM6 | 156 | +0.43 | 9.0% |
| BM4 | 90 | +0.39 | **24.4%** |
| ND_090617 | 81 | **+3.03** | **50.6%** |

- False-positive rate spans **0.0% → 50.6%**
- One healthy donor has a **median** plasma z above 3
- A second independent donor also fails, so it is not one outlier
- z > 5 still leaves 7.8% (BM4) and 9.9% (ND_090617)
- Reference imbalance does not explain it — BM4 is most skewed (0.958) but ND_090617
  fails worst at a balanced 0.678
- Leave-one-donor-out is flat

**Consequence, accepted rather than worked around:** the rule was *not* relaxed to
compensate. All 32 patients are `CNV_SUBCLONE_NOT_EVALUABLE`, and because only Level 3
licenses the word, **"subclone" is unavailable for every patient in this cohort.**

---

# What each mechanism does and does not license

| Mechanism | Licenses | Does **not** license |
|---|---|---|
| Marker panel | "this cluster is plasma / T / myeloid …" | any statement about malignancy |
| Light-chain restriction | "this patient has a dominant clone" | that an individual cell is a clone member |
| + V-gene usage | "this cell is a high-specificity dominant-clone member" | that it is genetically distinct |
| CNV | *nothing — not evaluable* | "there is no subclone"; "CNV-negative" |

`CLONE_SUPPORTED` is a **high-specificity dominant-clone core, not an exhaustive
malignant-cell set.**

