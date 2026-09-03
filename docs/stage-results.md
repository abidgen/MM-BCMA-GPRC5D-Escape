# Stage results — run output for stages 01 through 11b (incl. 05b, 08c, 09b)

Split out of the main project document on 2026-08-24. That document keeps, per completed
stage, a short block saying what it produced and what binds downstream; this file is the
full record.

**The frozen `results/` artifacts outrank the tables reproduced here — but only as
authenticated by the committed freeze manifest `provenance/frozen_artifacts_pre_stage12.tsv`.**
`results/` is gitignored, so an on-disk file is durable evidence only when its SHA256
matches that manifest; the committed producers under `production/` and `notebooks/` are
what define how each artifact was derived.
Where a table below duplicates a committed artifact the path is named; if the two ever
disagree, the CSV is right and this file is stale.

**Corrected 2026-08-24 — median genes/cell.** Earlier revisions of the main project
document quoted
"~2,044 median genes/cell" throughout. That figure came from the **R build** (61
samples, 181,336 cells, fixed thresholds) and was carried forward into the Python
rebuild without re-deriving it. Recomputed from
`results/05_integration/integrated.h5ad`:

| population | cells | median genes/cell |
|---|---:|---:|
| whole post-QC cohort | 172,940 | **1,162** |
| myeloma samples | 135,669 | 1,211 |
| donor marrows | 37,271 | 1,074 |
| plasma-like compartment (11 clusters) | 39,893 | **1,521** |
| everything else | 133,047 | 1,111 |

The real cohort is **shallower than the number the project had been arguing from** —
roughly half. Every argument that depended on it (dropout dominates for `GPRC5D`,
annotate at cluster level not per cell, within-clone CNV resolution is underpowered,
TC class is per-patient not per-cell) is *strengthened*, not weakened. Use **1,162**
for whole-cohort statements and **1,521** for statements about the malignant/plasma
compartment specifically.

---

**01-03 — Data acquisition** (`notebooks/01_download_data.ipynb`,
`02_check_files.ipynb`, `03_build_manifest.ipynb`; env: `mm-qc`). Notebooks wrapping the
verified `scripts/01-03`, which remain as a CLI fallback. Status: **all three written,
executed and passing (01/02 added 2026-08-21)** — 62/62 `triplet-ok`, manifest
byte-identical via both paths, all assertions green.

Scripts 01 and 02 are **bash**, so their notebooks `subprocess` out to them rather than
importing a function as notebook 03 does — that is what "wrap, don't duplicate" means
for a shell script. Notebook 02 parses the script's classification table into a
DataFrame for assertions, but never re-derives the classification itself.

Each adds what a notebook is for. **01**: a pre-flight stating whether the run will
download or skip (so a 1 GB transfer is never started by accident), post-conditions on
the 62 sample dirs and their `.extraction_complete` markers, and the bulk inventory
assertion that caught the 12-vs-13 WashU error above. **02**: assertions on the
`triplet-ok` classification, the reference split, a manifest-vs-disk cross-check, and
the **checksum finding** — 62 `genes.tsv` files, only **3 distinct checksums**, which is
the fact that makes the Ensembl-ID reconstruction in `gene_space.py` possible at all. It
asserts exactly 3, so a new reference build fails loudly instead of being merged
silently. **03** reports the reference split and runs the required-gene assertions
(which is how the `NSD2`/`WHSC1` symbol drift was found).

Note notebook 03's kernelspec said `mm-dual-antigen` (an R-build leftover) and was
repointed to `mm-qc` on 2026-08-21, matching the env this stage actually declares.

**04 — QC + doublets** (`notebooks/04_qc.ipynb`, `src/mm_escape/qc.py`; env:
`mm-qc`). Custom loader (handles `counts.mtx`/`genes.tsv` naming). QC metrics via
`scanpy.pp.calculate_qc_metrics` including `pct_counts_in_top_20_genes`. MAD-based
outlier filtering, re-derived per this cohort. `scDblFinder` via the `rpy2`
bridge. **MAD thresholds are derived per cohort, not pooled** — MMRF cells carry ~1.9x
the genes/cell of WU1's, so a pooled MAD would flag much of WashU cohort 1 as
low-quality for a batch reason (see the GEO metadata section). `56203_1` is loaded with
its repaired gene column, not excluded. Checkpoint each sample's post-QC AnnData
individually (mirrors the R build's resumable-per-sample design).

**RUN AND COMPLETE 2026-08-24.** 204,040 pre-QC cells over all 62 samples ->
**172,940 kept (84.8%)**. `results/04_qc/` holds the thresholds, the per-sample and
per-cohort reports, the MAD sensitivity sweep, three figures, and one checkpoint per
sample under `samples/` with **every barcode retained** and `obs["keep"]` set — cells
are annotated, never deleted, so stage 08 can ask what QC cost. The stage runs in two
passes: per sample (metrics + `scDblFinder`, ~210 s for the cohort) then per cohort
(MAD thresholds, ~11 s, off `obs` alone — it never concatenates the matrices).

Two things the data forced, both departures from what this section originally
specified, and both made on this cohort's numbers rather than on preference:

1. **`pct_counts_in_top_20_genes` is computed and reported but does NOT filter.**
   A 5-MAD band on it flags 17% of MMRF and 15% of WU1 against 3% of WU2 — too uneven
   to be catching one thing. Inspecting those cells (MMRF_1695, top decile) shows two
   populations: `IGKC` at ~25% of counts (plasma cells) and `HBB`/`HBA1/2` at ~32%
   (erythroid debris). The plasma-cell half is the project's subject — `TNFRSF17`
   detected in **21.8%** of that decile vs **0.8%** elsewhere, `SDC1` **18.8%** vs
   **0.0%**. A plasma cell is a professional secretor; an Ig-dominated library is its
   normal state, not a defect. Filtering on it would preferentially delete
   antigen-**positive** malignant cells and inflate `frac_double_negative`. The metric
   is kept because it is one of the few ambient-Ig handles available at all (SoupX
   needs unfiltered matrices this deposit lacks) — it is just not allowed to delete
   cells. See `qc.DEFAULT_FILTERS` vs `qc.ALL_FLAGS`; opting it back in is one
   explicit argument, for a sensitivity re-run.
2. **The deposit is already filtered, differently per cohort — see the correction
   below.**

Thresholds this cohort actually produced (`results/04_qc/qc_thresholds.csv` has all of
them; `pct_counts_mt` is one-sided because a cell with unusually *few* mitochondrial
reads is not low quality):

| cohort | log1p_total_counts band | pct_counts_mt cap | % removed |
|---|---|---|---|
| MMRF | [4.24, 13.60] (MAD 0.94) | 12.6% (never binds — see below) | 4.7% |
| WU1 | [6.10, 10.24] (MAD 0.41) | 8.6% | 14.1% |
| WU2 | [5.83, 10.37] (MAD 0.45) | 12.9% | 18.0% |
| Donor | [6.37, 10.08] (MAD 0.37) | 11.7% | 22.6% |

MMRF's band is ~2x wider than the others' because its own depth distribution is much
more dispersed, so at 5 MADs it removes nothing but doublets. The removal rate is
stable from 4 to 6 MADs (16.5% -> 14.6% overall), so nothing downstream hangs on the
exact count. `min_genes = 200` flags nothing — deliberately, it is a safety net, and a
higher floor would hit the shallow cohorts hardest.

**CORRECTION — the deposit IS pre-filtered, and differently in each cohort.** An
earlier note in this project held that the depositors' stated Seurat filter (drop
cells above 10,000 UMIs) "was not applied to what is deposited", reasoning from a
cohort-wide average UMI count. That average pooled MMRF with WashU and hid a
per-cohort truth. Read per cohort the boundaries are unmistakable — a `max` of exactly
9,999 or 10.00 is a cutoff, not a distribution:

| cohort | UMI ceiling | UMI floor | pct_mt ceiling | gene floor |
|---|---|---|---|---|
| WU1, WU2 | **< 10,000** | >= 1,000 | < 20% | >= 200 |
| MMRF | none (max ~269,000) | >= 1,000 | **< 10%** | >= 200 |
| Donor | none (max ~119,000) | none | < 20% | >= 200 |

**This is a first-order problem for stage 08, not bookkeeping.** Malignant plasma
cells are the highest-RNA-content cells in marrow, so a 10,000-UMI ceiling does not
remove a random slice. Measured in the uncensored cohorts, where the band is still
visible, cells above 10,000 UMIs are enriched **3-21x for `TNFRSF17`** and **20-70x
for `GPRC5D`** (`results/04_qc/umi_censoring_effect.csv`). So **36 of the 54 myeloma
samples had the antigen-positive tail of their own tumours removed before deposit**,
which inflates `frac_double_negative` for WU1/WU2 relative to MMRF — a bias in the
project's own direction of interest, and one that is baked in and cannot be undone.

Deliberately **not** corrected at stage 04. Truncating MMRF and Donor to match would
discard 42% of MMRF's cells to make every cohort equally damaged; QC's job is to
remove bad cells, not destroy good ones to equalise two cohorts. It is carried forward
as a quantified confounder, and **stage 08 must run the truncate-everything-at-10k
version as a sensitivity analysis**, where it costs nothing and answers the question
directly. Note this also means the "1.9x cohort depth gap"
(`docs/dataset-ground-truth.md`, GEO series metadata) is partly
*censoring*, not only chemistry.

**05 — Gene-space intersection, integration, clustering**
(`notebooks/05_integration_clustering.ipynb`, `src/mm_escape/gene_space.py` +
`integration.py`; env: `mm-core`). **Canonicalize gene symbols first, then** intersect
gene sets across retained samples (hard-fail with specifics if required genes don't
survive). The harmonization step is not optional — see the Data section's symbol-drift
table; without it `NSD2` is dropped and stage 10 loses t(4;14) entirely. Report how many
genes the harmonization recovers, so a regression here is visible rather than silent.
`anndata.concat(join="inner")`. Normalize, HVG, PCA, `harmonypy` keyed on
`patient_id` **with `n_genes_ref` AND `cohort` as additional covariates** — the build
and the cohort are different axes and neither substitutes for the other (two WU1
samples sit on the 33538 build, the four `ND_*` donors on 33694), Leiden clustering,
UMAP. Diagnostic UMAP colored by reference version (`n_genes_ref`) to confirm the
intersection actually neutralized processing batch. **This three-covariate Harmony
configuration was benchmarked against six alternatives at stage 05b (2026-08-24) and
retained** — it is now a tested choice rather than a default.

**RUN AND COMPLETE 2026-08-24.** 172,940 cells x **32,991 genes**, 30 Leiden
clusters, `results/05_integration/integrated.h5ad` (1.3 GB gzipped). Whole pipeline
~190 s, **peak ~20 GB RAM** — the one stage that actually concatenates the matrices,
so it is the machine-size constraint for the project. Gene space came out exactly as
`gene_space.py` predicted: 22,164 on raw symbols -> **32,991 on Ensembl IDs
(+10,827)**, 11,140 drifted symbols joined correctly, all required genes present with
`NSD2` resolving against `WHSC1` in the older build.

**Harmony works on the immune compartment and does NOT work on the plasma-cell
compartment — and that is the stage-04 censoring showing up a second time.** Median
per-cluster cohort-mixing entropy:

| | clusters | cells | median entropy by cohort |
|---|---|---|---|
| plasma-cell-like (`MZB1` > 40%) | 11 | 39,893 | **0.105** |
| everything else | 19 | 133,047 | **0.751** |

(Uncorrected PCA for reference: 0.341 median over 54 clusters, 34 of them below 0.5.
Harmony: 0.621 over 30, 13 below 0.5. The correction is doing real work — it is just
doing it unevenly across compartments.)

The three largest plasma-cell clusters are **one per cohort**, each spanning ~30
patients. That rules out the benign explanation: a patient-private clone would
fragment into ~41 clusters, not three cohort-shaped ones. The likely cause is the
stage-04 finding — WashU was cut at 10,000 UMIs and MMRF was not, plasma cells are
the highest-RNA-content cells in marrow, so **WashU's plasma cells are a truncated
subset of the plasma-cell distribution** and no batch-correction method can restore
cells that were never deposited. It is compartment-specific for the same reason:
T/NK/myeloid/B cells sit well below 10,000 UMIs everywhere, so the ceiling never
touched them.

**Contained, not fatal** — and contained by decisions made before it was observed:
per-cell antigen calls are raw counts and never touch this embedding; stage 10's
malignant subclustering is per-patient and un-integrated; stage 06 annotates at
cluster level, where three cohort-specific plasma-cell clusters all annotate as
PlasmaCell at no cost. **What it does forbid: reading any cross-cohort comparison of
malignant-cell state off this embedding.** And it makes stage 08's cohort covariate
and the truncate-all-at-10k sensitivity analysis mandatory rather than advisable —
this is the second independent sign of the same problem.

**Preprocessing, and two findings from reporting it.** CP10K + `log1p` with raw
integers kept in `layers["counts"]`; 2,000 HVGs of 32,991 selected with
`batch_key="patient_id"` (median HVG is variable in 15 of ~50 patients, min 10, so the
selection is not one cohort's artefact); scale + PCA to 50 components on a throwaway
HVG-subset copy, so `X` stays unscaled and complete. 50 PCs is 28.4% of variance and
PCs 31-50 add only 2.7 points, so the count is generous rather than load-bearing.

1. **`GPRC5D` is not a highly variable gene** — mean **0.061** against `TNFRSF17`'s
   **0.492**, an 8x gap, and HVG in only 6 patients. This affects nothing: the
   embedding does not need it to find plasma cells, and stage 08 reads
   `layers["counts"]`. **Its value is evidential.** This document argues repeatedly
   that dropout matters more for GPRC5D than for BCMA because GPRC5D is a low-abundance
   GPCR transcript; that has been an assertion from the literature, and this is the
   first number from *this cohort* supporting it. A materially higher share of
   "GPRC5D-negative" calls will be technical zeros, so **GPRC5D-negative calls warrant
   more scepticism than BCMA-negative ones** — which is what stage 08's
   expression-matched false-negative floor exists to quantify. The antigen panel is
   deliberately **not** forced into the HVG set: that would bias the embedding toward
   the genes under study for no gain.

2. **The plasma-cell integration failure IS the stage-04 censoring, now measured
   rather than inferred.** Median UMIs per cell, split by compartment:

   | compartment | MMRF | WU1 | WU2 | MMRF/WU1 |
   |---|---|---|---|---|
   | non-plasma | 5,829 | 3,273 | 2,879 | **1.8x** |
   | plasma-like | 22,477 | 5,036 | 4,888 | **4.5x** |

   MMRF's two largest plasma-cell clusters have **68%** and **88%** of their cells
   above 10,000 UMIs — cells the WashU deposits cannot contain, because WashU was cut
   at that ceiling. WashU's plasma clusters instead press up against it (7.5-8.3% of
   cells in the 9,000-9,999 band, against 1.2% for MMRF's cluster 7). So Harmony is not
   failing on plasma cells: what separates them is a **non-recoverable
   sampling/censoring asymmetry**, not established biology. WashU's *observed*
   plasma-cell distribution is missing its high-RNA portion, so **no one-to-one
   population correspondence remains for any method to recover**, and no correction
   restores cells that were never deposited. (An earlier wording here said "the
   populations genuinely differ", which implied biological divergence between cohorts
   and claimed more than the data supports.) The compartment-specificity follows: T/NK/myeloid/B sit well
   below 10,000 UMIs everywhere, so the ceiling never touched them, and Harmony mixes
   them to 0.75. `results/05_integration/depth_by_compartment.csv` and
   `plasma_cluster_depth.csv`; two regression tests pin the asymmetry.

One defect found and fixed: `gene_space.to_canonical_symbols` named the `var` index
`canonical_symbol` while also keeping a `canonical_symbol` column holding the
*unsuffixed* symbol, which differs for the 9 collision-resolved genes. AnnData
refuses to write such an index, so this failed only at `write_h5ad` — after every
in-memory test had passed. The index is now named `symbol`, and
`tests/test_integration.py` round-trips through `.h5ad` so a serialization-only bug
cannot hide again.

**Constrain integration's blast radius (new).** Harmony keyed on `patient_id` is
correct for the immune compartment but carries a real risk for the tumor: the
malignant clone is *patient-private by definition*, so forcing patients together
can blend genuinely distinct clones into one blob and erase the heterogeneity this
project exists to measure. Therefore:
- Integrated embedding is used for **immune-compartment annotation and clustering
  only** (stages 06, 11).
- **All malignant subclustering is done per patient, un-integrated** (stage 10
  depends on this and must not read the Harmony embedding).
- State explicitly in the notebook that **per-cell antigen calls are raw counts and
  are therefore integration-independent** — this is what contains the risk, and it
  is the answer to "did Harmony distort your escape fractions?" (it cannot; the
  calls never touch the embedding).

**05b — Integration-method benchmark** (`notebooks/05b_integration_benchmark.ipynb`,
`src/mm_escape/benchmark.py`; env: **`mm-integration`**). Added 2026-08-24.
`results/05b_benchmark/`.

**Why `05b` and not a number.** Number order is execution order with no exceptions.
This is a side-comparison feeding the stage-05 *choice*, not a new pipeline stage, so
it takes a letter rather than displacing annotation. It reads
`results/05_integration/integrated.h5ad` **read-only** and asserts the file is
byte-identical afterwards; candidate embeddings live in `05b_benchmark/` until a winner
is picked.

**Why it exists.** Stage 05 used Harmony because it was the obvious default, never
because it beat anything. `sc-best-practices`' integration chapter recommends running
several methods and scoring them with scIB rather than assuming one wins.

**Why the scoring is deliberately NOT standard scIB — the design point that matters
more than the leaderboard.** scIB's batch metrics cannot distinguish *"correctly left
apart"* from *"failed to merge"*. Given stage 04's censoring, a method that squashes
the three plasma islands together scores better on kBET/iLISI while manufacturing
correspondence where none is recoverable. **A naive global ranking would structurally
reward overcorrection on this dataset.** Therefore:

- **The immune compartment is scored; the plasma compartment is diagnosed.** Global
  scIB is computed and reported as a **secondary reference, never used for selection**.
- Plasma mixing **never contributes positively**. A jump alone flags an arm for
  inspection; a jump **together with** rising depth association disqualifies it, that
  pair being the signature of the censoring being smoothed over.

**Arms — one common batch key, plus reference arms.** All scored on `cohort` regardless
of what they corrected on: with three different batch definitions in play, scoring each
against its own key would compare seven different *questions* rather than seven
*methods*. `cohort` is the common axis, samples nest inside it, and it is where the
demonstrated distortion lives.

| arm | corrects on |
|---|---|
| `unintegrated` | — (`X_pca`, same HVGs/scaling/PCs) |
| `harmony_sample`, `scvi_sample`, `scanorama_sample` | `sample_name` |
| `harmony_stage05` *(incumbent)* | `patient_id` + `n_genes_ref` + `cohort` |
| `harmony_cohort`, `scvi_cohort` | `cohort` |

`sample_name` is the true technical unit, but **42 of 50 patients contribute exactly
one sample**, so `sample_name` and `patient_id` are nearly the same partition and
"avoid correcting on patient" is weaker than it looks. **`cohort` is the only batch
definition not confounded with patient**, which is why it is a real arm.

**Labels are provisional CellTypist, deliberately embedding-independent.** scIB scores
batch removal *against* bio conservation, and the bio half needs labels — which come
from stage 06, which consumes the embedding under selection. CellTypist breaks that
circle by classifying from expression. **`majority_voting=False` is load-bearing**:
majority voting smooths over an over-clustering, which would smuggle an embedding back
into the labels. Two things were checked rather than assumed:
- **`Immune_All_High`'s `ILC` class is NK** (8% of marrow; NKG7 98.7%, GNLY 92.2%,
  KLRD1 85.9%, MS4A1 1.3%).
- **It does cover erythroid and HSPC here**, contrary to the expectation recorded
  elsewhere in this document that an immune-only reference would be blind to them
  (Erythroid 14,103 cells at HBB 99.7%; HSC/MPP 2,625 at CD34 58.4%). So **no hand-set
  marker thresholds enter the benchmark at all.**

**The decision rule is declared before running** (`benchmark.DECISION_TOLERANCES`),
exactly as stage 06 pre-declares its F1 thresholds. A method replaces the incumbent only
if all four hold: `batch_improved`, `bio_preserved`, `depth_ok`, `overcorrection_ok`.
The verdict is **computed, not narrated** — `decision.csv` is the source of truth and
`decision.md` is rendered from it. **Harmony is allowed to win**, including against a
higher conventional global scIB score.

`depth_ok` uses **`R²(log1p(total_counts) ~ latent)`**, fixed in advance. R² depends
only on the embedding's **column span** and is therefore **rotation-invariant** — latent
axes are arbitrary across methods, so a per-dimension `max |Spearman|` would rank
methods on an accident of their parameterisation.

**PRIMARY RESULT (immune compartment, scored on `cohort`):**

| arm | batch correction | bio conservation | total |
|---|---|---|---|
| unintegrated | 0.450 | 0.691 | 0.595 |
| **`harmony_sample`** | **0.615** | **0.718** | **0.677** |
| `scvi_sample` | 0.570 | 0.701 | 0.649 |
| `scanorama_sample` | 0.450 | 0.723 | 0.614 |
| **`harmony_stage05` (incumbent)** | **0.427** | 0.700 | 0.591 |
| `harmony_cohort` | 0.591 | 0.706 | 0.660 |
| `scvi_cohort` | 0.492 | 0.690 | 0.611 |

**The headline is about configuration, not method.** Harmony wins — but the stage-05
*configuration* is the **worst batch corrector of all seven arms (0.427), below even no
integration at all (0.450)**, while plain Harmony on `sample_name` leads on batch
removal *and* is second only to Scanorama on bio conservation. Correcting on
`patient_id` + `n_genes_ref` + `cohort` apparently split the correction across
covariates rather than strengthening it. **This finding exists only because the
incumbent was entered as its own arm rather than assumed**, and it is the argument for
doing the benchmark at all.

**VERDICT: no arm qualified — `harmony_stage05` stays.** Which is the outcome the rule
was explicitly written to allow, and the reason for writing it in advance.

The full table (`results/05b_benchmark/decision.csv`, rendered into `decision.md`):

| arm | batch | bio | depth R² | plasma mixing | vs incumbent | eligible |
|---|---|---|---|---|---|---|
| unintegrated | 0.450 | 0.691 | 0.660 | 0.014 | depth +0.291 | no |
| `harmony_sample` | **0.615** | **0.718** | 0.509 | 0.515 | depth +0.140, plasma **13.5x** | no |
| `scvi_sample` | 0.570 | 0.701 | 0.541 | 0.452 | depth +0.172, plasma **11.8x** | no |
| `scanorama_sample` | 0.450 | 0.723 | 0.690 | 0.161 | depth +0.321 | no |
| **`harmony_stage05`** | 0.427 | 0.700 | **0.369** | **0.038** | — | (incumbent) |
| `harmony_cohort` | 0.591 | 0.706 | 0.607 | 0.771 | depth +0.238, plasma **20.2x** | no |
| `scvi_cohort` | 0.492 | 0.690 | 0.576 | 0.017 | depth +0.207 | no |

**The arms that win on conventional scIB are precisely the arms that merge the censored
plasma populations.** `harmony_sample` posts the best batch *and* bio scores — and
mixes the plasma compartment **13.5x** harder than the incumbent while encoding more
depth. `harmony_cohort` reaches **20.2x**. Meanwhile the two arms that leave the plasma
populations apart (`unintegrated` 0.014, `scvi_cohort` 0.017) are the ones with no
meaningful batch gain. **A standard global scIB benchmark would have selected
`harmony_sample`**, which buys its higher score substantially by fusing populations that
cannot be fused. That is the failure mode this design was built to catch, and it
occurred.

The incumbent is simultaneously the **worst batch corrector (0.427, below unintegrated's
0.450)** and by a wide margin the **least depth-encoding (R² 0.369 against 0.51-0.69)**
and the **least plasma-merging (0.038)**. Correcting on `patient_id` + `n_genes_ref` +
`cohort` evidently trades cohort mixing for exactly the properties this dataset needs.

**Stated honestly: `depth_ok` did all the gating.** Every non-incumbent arm failed it,
so the +0.05 tolerance — fixed before the spread (0.37-0.69) was known — is what
produced a clean sweep. Two things keep that from being a threshold artifact:
`harmony_sample`, `scvi_sample` and `harmony_cohort` **independently fail
`overcorrection_ok`**, so they lose even with the depth criterion removed entirely; and
relaxing depth enough to admit anything admits only `scanorama_sample` and
`scvi_cohort`, the two arms with the weakest batch gains. The tolerance is **not**
re-tuned after the fact — doing so is precisely the post-hoc rationalisation
pre-declaration exists to prevent.

**Also recorded: scVI encodes depth in plasma cells, badly.** Plasma-compartment
R²(depth ~ latent) is 0.793 (`scvi_sample`) and 0.850 (`scvi_cohort`) against the
incumbent's 0.528 and `harmony_sample`'s 0.319. Its explicit library-size model appears
to put depth *into* the latent space for the compartment where depth is the
confound — the opposite of the reason it was the principled candidate.

Runtimes (RTX 5070): Harmony ~15 s per arm, Scanorama ~119 s, scVI ~190 s.

**What the benchmark cannot do, and this must be propagated to stage 08.** **No
integration method restores cells that were never deposited.** A well-mixed latent
space has not undone the ascertainment bias in the raw counts stage 08 reads. Whichever
arm wins, **stage 08 still owes its truncate-all-cohorts-at-10,000 sensitivity
analysis** — selecting a fancier method must not create the impression the censoring was
"handled". And nothing here can move `frac_double_negative` at all: the embedding feeds
only stages 06 and 11, antigen calls read `layers["counts"]`, and malignant
subclustering is per-patient un-integrated.



---

## Module validation logs

One-time verification records moved here from the main project document's status section
on 2026-08-24. Not standing instructions.

**`io.py` — validated against the real files (2026-08-24).** All **62** samples load in
~2 s; the four-sample failure-mode set (`MMRF_1695` = 33538 build, `27522_1` = 33694
build, `BM4` = donor control, `56203_1` = truncated deposit) round-trips through
`gene_space.py` to 32,991 genes with all 65 required genes present. Verified: the
genes × cells → cells × genes transpose (against random raw `.mtx` triplets plus exact
total-count and nnz equality), deposited gene order preserved untouched for
`attach_ensembl_ids`'s positional join, `<sample_name>_<barcode>` obs_names unique
across the concat, the `56203_1` prefix repair plus its two failure modes, and the GEO
metadata join. Cohort **pre-QC** totals: **204,040 cells over 62 samples**
(155,890 myeloma / 48,150 donor).

`load_manifest()` emits `cohort`, `chemistry`, `dead_cell_removal` and `diagnosis` from
`resources/sample_metadata/`, and every cell carries them — stage 04 needs `cohort`
before it can derive thresholds.

**`qc.py` + `notebooks/04_qc.ipynb` — run 2026-08-24.** 204,040 → 172,940 cells;
per-cohort thresholds derived and documented; `scDblFinder` run on all 62 samples; one
checkpoint per sample under `results/04_qc/samples/` with every barcode retained.

**`integration.py` + `notebooks/05_integration_clustering.ipynb` — run 2026-08-24.**
172,940 × 32,991, 30 clusters, embedding at `results/05_integration/integrated.h5ad`.

**`benchmark.py` + `notebooks/05b_integration_benchmark.ipynb` — run 2026-08-24.**
`envs/env-integration.yml` built, `mm-integration` kernel registered,
`tests/test_benchmark.py` (24 data-free tests) committed. The incumbent survived;
nothing about stage 05's output changed as a result.

## Test suite counts (regression baselines)

| condition | result |
|---|---|
| fresh clone, no `raw/` | 100 pass, 57 skip |
| with the deposit present | 155 pass, 2 skip, ~27 s |
| **2026-08-26, stages 01-11 complete** | **457 pass, 1 skip, 4 deselected, ~32 s** (`pytest -m "not slow"`, `mm-core`) |
| **2026-08-26, + stage 08c coverage** | **488 pass, 1 skip, 4 deselected, ~33 s** |
| **2026-08-26, + stage 11b LIANA arm** | **507 pass, 1 skip, 4 deselected, ~32 s** |
| **2026-09-03, post-Stage-12, full suite** | **620 passed, 2 skipped, 6 warnings** (`pytest -q`, `mm-core`) |
| **2026-09-03, post-Stage-12, not-slow** | **616 passed, 1 skipped, 5 deselected** (`pytest -m "not slow"`) |

Gates: `conftest.requires_data` (needs `raw/`), `requires_s1` (Supplementary Table S1
is a journal file, not part of the GEO deposit), `requires_r` (scDblFinder lives only
in `mm-qc`; the suite's home `mm-core` carries no R). `pytest -m "not slow"` skips the
two full-cohort passes.

---

## Stage 06 — annotation: full methodology and revision history

Moved out of the main project document 2026-08-25 when stage 06 was frozen. This is the
record of *how* the accepted C2d labels were produced and why three revisions were
needed; the rules that still bind downstream work stay in that document. **No rule was
changed in the move.**

**06 — Annotation** (`notebooks/06_annotation.ipynb`, `src/mm_escape/annotation.py`;
env: **`mm-annotation`**). Follows `sc-best-practices.org`'s annotation chapter.
**Three methods are run and compared, and the choice is made per class against
thresholds declared in advance** — not "`celltypist` and/or marker scoring", which
would leave a load-bearing decision to be settled implicitly by whatever ran first.

**What is load-bearing.** This stage feeds three things and only three:

| Downstream need | Labels that matter | Cost of getting it wrong |
|---|---|---|
| Stage 07 malignant calling | PlasmaCell vs. everything else | wrong plasma-cell set → wrong denominator for `frac_double_negative` |
| Stage 08 ambient noise floor | T / NK / myeloid purity | a plasma cell leaking into the "confidently antigen-negative" population inflates the floor and biases every antigen call |
| Stage 11 confounder control | T/NK abundance, ideally with subsets | composition artifact misread as immune evasion |

Fine-grained subtypes (CD4 Tcm vs. Tem) benefit stage 11 only and **must never be the
reason a method is chosen**.

**Method A — manual, marker-based**, at cluster level (not per cell; clustering absorbs
dropout, which matters at 1,162 median genes/cell). `scanpy.tl.score_genes` over the
project panel: PlasmaCell (SDC1/CD38/MZB1/XBP1/IRF4), Bcell (MS4A1/CD79A/CD19), Tcell
(CD3D/CD3E/CD8A/CD4), NK (NCAM1/NKG7/GNLY), Myeloid (CD14/LYZ/ITGAM), Erythroid
(HBB/GYPA), HSPC (CD34/KIT). Then `dotplot(..., standard_scale="var")` as saved
assignment evidence, and `rank_genes_groups(method="wilcoxon")` +
`filter_rank_genes_groups`. **The DE step is not optional** — it is the only thing that
can surface a population the seven-class panel does not cover (pDC, erythroid
progenitors, a doublet-driven cluster). Record ambiguous clusters as ambiguous.

**Method B — `celltypist`.** Normalize to 10,000 counts/cell then `log1p` (its stated
requirement); run on **expression, never the Harmony embedding**. Enumerate models with
`models_description()` rather than assuming names; default `Immune_All_Low.pkl` +
`Immune_All_High.pkl`, plus a healthy-BM model if one ships.
`annotate(..., majority_voting=True, over_clustering=<stage-05 leiden key>)` — passing
the existing Leiden key is what makes the methods comparable (same partition, different
labelings). Retain `conf_score`.

**Method C — `SingleR`**, chosen to cover CellTypist's blind spot rather than duplicate
its strengths: `Immune_All_*` is immune-only, so its predictable failure is erythroid and
HSPC, which `celldex`'s `NovershternHematopoieticData` covers
(`HumanPrimaryCellAtlasData` as fallback; verify both exist in the installed `celldex`).
Run with `clusters=<leiden key>`; retain `pruned.labels` and the delta/score matrix as
SingleR's own low-confidence signal.

**Two caveats about plasma cells, pointing opposite ways.** The references contain
*normal* plasma cells, so malignant PCs will be labelled "Plasma cells" — correct and
sufficient, since malignant-vs-normal is stage 07's job; do not count a method down for
failing to find the tumor. But conversely, because no malignant class exists in the
reference, a heavily aneuploid clone may be labelled something else or split across
labels. **Run the plasma-cell marker-coverage test on myeloma marrows specifically, not
only the donors** — otherwise a systematic failure on exactly the cells this project
measures passes unnoticed. Stage 05's three cohort-specific plasma clusters cost
annotation nothing (all three annotate as PlasmaCell) but make this check more important.

**The comparison** (artifacts to `results/06_annotation/`):
1. Confusion matrices — manual × CellTypist, manual × SingleR, **CellTypist × SingleR** —
   at cluster and cell level, plus ARI and per-class F1/Jaccard. ARI alone is
   insufficient: a method can score well overall while failing on plasma cells, the one
   class that must not be wrong. The two automated methods agreeing with *each other* is
   the strongest evidence available, being trained on different references.
2. **The marker-coverage test — the decisive one.** Dotplot the *manual* panel grouped by
   each *automated* method's labels. If CellTypist's T cells are CD3D/CD3E-high, its
   plasma cells MZB1/SDC1-high, and so on, the automated labels already encode what the
   manual panel encodes.
3. Confidence/coverage: `conf_score` per cluster, SingleR pruned-`NA` rate and deltas,
   the fraction of cells unassigned or labelled outside the panel.

**These are concordance scores, not validation accuracy.** The manual annotation is a
third opinion from the same expression matrix, not ground truth — "F1 against manual"
measures *agreement*, and the two automated methods agreeing is agreement between
references that share marker-biology priors and may share their blind spots. Label it
all **concordance**. The load-bearing evidence of *biological* validity is the
marker-coverage test: a label set can be perfectly self-consistent and biologically
wrong. Headline figure, in this order — **label concordance → marker-expression validity
→ uncertainty/unassigned** — concordance first because it is the weakest of the three.

**Benchmarked against the original paper's annotation (2026-08-25, diagnostic only).**
The source publication (Cancer Res 2023;83:1214-1233, PMC10102848) annotated by manual
marker review at cluster level — no SingleR, no atlas, no exclusion criteria. C2d agrees
with it on **24 of 30 Leiden clusters (78.9% of cells)**.
- **The plasma-cell boundary is essentially concordant: 32,307 of 32,337
  plasma-compartment cells agree** (the only disagreement is 30 cells in Leiden 11). So
  the denominator of `frac_double_negative` is not in dispute, and **stages 07 and 08 are
  unaffected by any of this**.
- **The one consequential disagreement is the cytotoxic-lymphocyte compartment.** The
  paper defines NK as `NKG7` + `GNLY` with no requirement that CD3/TCR be absent, which
  cannot separate conventional NK from cytotoxic alpha-beta T or gamma-delta T. Leiden 3
  and 12 (22,132 cells) carry a full TCR complex with `KLRF1` at 8-13% and `NCAM1` at
  1-3%; Leiden 23 is NK-enriched (`KLRF1` 75%) but carries `TRDC` at 66.4%.
- **Local analysis of Leiden 23 is therefore motivated by the stage-11 T/NK objective,
  not by any problem with the plasma boundary.** Full benchmark:
  `results/06_annotation/paper_annotation_benchmark/`.
- Where the paper was *better*: its cDC markers (`FCER1A`/`CLEC10A`) never had the MHC-II
  exposure that broke our C2 attempt, and C2d's cDC anchors converged on those same two
  genes. Its three-gene plasma definition already carried the two-axis logic our
  mature-plasma predicate formalises.

**T-lineage identity requires coordinated machinery (2026-08-25).**
**Isolated `TRBC1`/`TRBC2` expression is insufficient evidence of T-lineage
commitment, because it frequently occurs without coordinated `CD3`/`TRAC` expression in
cells with strong NK-lineage evidence.** Validated on Leiden 23: among 5,788 cells with
strong NK evidence and no CD3/TRAC, 74.8% were TRBC-positive; among TRBC-positive mixed
cells only 12.9% carried both CD3 and TRAC. Candidate mechanisms — germline/unrearranged
TRB transcription, ambient spillover, residual multiplets — are **hypotheses**; the
diagnostic did not distinguish them (TRBC intensity exceeds a one-molecule ambient
profile, TRBC positivity does *not* rise with sample T-cell abundance, r = -0.163,
p = 0.3, and the doublet-score gap is too small to explain the compartment). Do not
state any mechanism as established.

**This is not a threshold change.** Of the 4,080 cells that move mixed -> NK under the
revision, only 30.2% carry a single TRBC UMI and 100% remain TRBC-positive at a median
of 2 UMI. The transition is driven by the requirement for coordination with
lineage-defining machinery, not by thresholding away weak signal. `config.T_IDENTITY_ANCHORS`
(CD3D/E/G, TRAC) and `config.T_CONTEXT` (TRBC1/2) are predeclared; `MARKER_PANEL`,
`LINEAGE_PROGRAMS` and every numeric parameter are unchanged, and the revision is so far
**validated on cluster 23 only** — global C2d is untouched.

**Leiden clusters are evidence-aggregation units, not assumed ground-truth cell types.**
This is the premise the whole stage rests on and it is deliberate. At 1,162 median
genes/cell a per-cell marker call on a dropped-out gene is a *wrong* call rather than a
missing one, so evidence is pooled at cluster level to stabilise it — the cluster is a
denominator for detection fractions, not a claim that its cells are one type. Two
consequences follow. A cluster whose evidence is genuinely mixed **stays `Ambiguous`**
rather than being forced into the nearest broad label; and `Ambiguous` is then a
*signal*, marking exactly where targeted local subclustering is warranted (Leiden 23 is
the live example — NK-enriched, `TRDC` 66.4%, unresolved on purpose). Never resolve a
mixed cluster by loosening a global rule; resolve it locally or leave it unresolved.

**Known limitation — the Myeloid contradiction programme is not route-symmetric
(2026-08-25, documented, not yet fixed).** Broad Myeloid identification now has three
route-specific evidence sets (monocyte, cDC, pDC), each anchor+context. The exclusion
side is still the single conservative panel `CD14`/`FCN1`/`MNDA`/`ITGAM`, which
represents the monocyte route and neither of the others — a cDC or pDC population can
satisfy identification while raising almost no Myeloid contradiction elsewhere. This
**favours specificity over sensitivity**: contradictions are under-called, never
invented, which is the safe direction and the same principle as `CONTRADICTION_MIN_GENES`.
A route-aware revision (add `CSF1R`, `VCAN`, `MS4A7`, `CD300E`; exclude `FCGR3A` as
NK-borne; never add `LST1`/`LYZ`/`CTSS`/MHC-II) is proposed but **deliberately deferred
to a separate experiment** — one variable at a time.

**Stage 06 has run three times. v1 (`results/06_annotation_v1/`) and v2
(`results/06_annotation_v2/`) are preserved and were not accepted; v3 is the live
result.** Each revision corrected a *definition* or a *method* that measurement showed
to be broken — never an acceptance bar. The F1 bars, `MARKER_COVERAGE_MIN`,
`CONTRADICTION_MIN_GENES` and `CONTRADICTION_MAX_RATE` are identical across all three.

**v3 — manual identity is adjudicated on evidence, not module-score magnitude.**
v1/v2 assigned each cluster by argmax over `scanpy.tl.score_genes` from different panels.
Those scores subtract a control set drawn from each gene's own expression bin, so every
panel carries its own baseline offset and **they are not comparable across panels**. On
this dataset the T panel sits 0.2036 below the NK panel purely from that effect — larger
than the 0.068 and 0.196 margins by which two T-cell clusters were called NK. Clusters
are now adjudicated on **detection fractions**, which share a common [0, 1] scale:
`manual annotation = positive lineage evidence + specificity/exclusion evidence`
(`MANUAL_MARKER_DETECT_MIN` 0.25, `MANUAL_POSITIVE_MIN` 0.5, `MANUAL_DECISION_MARGIN`
0.15, exclusion reusing `CONTRADICTION_PAIRS`/`CONTRADICTION_MAX_RATE`). `score_genes`
remains a descriptive within-program quantity and no longer decides identity.
**`Ambiguous` is a real outcome**: a cluster with balanced evidence for two lineages is
unresolved at this resolution, and saying so beats a tie-break that manufactures a call.

 v1 executed the pre-declared rules faithfully and
without retuning, and still produced NK = 33,556 against Tcell = 19,133 in bone marrow,
with Erythroid the runner-up class in 18 of 30 clusters. The diagnosis was two
*definitional* defects, not threshold ones: the manual Erythroid panel was `HBB` + `GYPA`
and stage 04 had already shown haemoglobin to be the dominant ambient species (HBB is
detected in 61-85% of cells of every class), and the NK panel rested on `NKG7`/`GNLY`,
which are a cytotoxic-granule program shared with cytotoxic T cells. **The correction was
to the manual reference specification and the validation framework. The F1 bars and
`MARKER_COVERAGE_MIN` were not changed** — see `docs/decisions-archive.md`.

**The decision rule — per class, declared before looking**, so "choose the best" cannot
become post-hoc rationalization. A class goes to an automated method when its
marker-coverage test passes, its own confidence signal is not flagging the cluster, and
concordance with manual clears a pre-set F1 bar:
- **PlasmaCell: F1 ≥ 0.95** — strictest, because it sets the metric's denominator.
- **T / NK / Myeloid: F1 ≥ 0.90** — these define stage 08's noise floor.
- **Bcell / Erythroid / HSPC: F1 ≥ 0.85** — nothing downstream is load-bearing on these.
- **Marker coverage ≥ 0.30** (`config.MARKER_COVERAGE_MIN`), for every class and every
  method — the minimum mean scaled expression of a class's own `MARKER_PANEL` markers,
  within the cells a method assigned to that class. Declared 2026-08-25, **before any
  coverage result was computed**, for the same reason the F1 bars are. This is the
  first veto: a class failing it is rejected regardless of concordance, because high
  agreement on a biologically unsupported label is agreement on an error.
- **Lineage-contradiction rate ≤ 0.25** (`config.CONTRADICTION_MAX_RATE`) — the second
  veto, added for v2 and declared before any v2 result existed. See below.

**Lineage exclusivity — the second veto (added v2, 2026-08-25).** Marker coverage asks
a *precision*-like question: do the cells labelled X express X's markers? It therefore
**cannot see a class that has swallowed another lineage**, which is a recall failure.
Stage 06 v1 shipped 33,556 "NK" cells at coverage 1.00 of which 63% were CD3D⁺ — they
genuinely were NKG7/GNLY-high, because cytotoxic T cells are, and coverage had no way
to object. So a complementary check runs alongside it:

| constant | value | meaning |
|---|---|---|
| `LINEAGE_PROGRAMS` | 4 programs | ambient-robust evidence per lineage |
| `CONTRADICTION_PAIRS` | per class | which lineages are incompatible with each class |
| `CONTRADICTION_MIN_GENES` | 2 | genes **detected** for "strong evidence" |
| `CONTRADICTION_MAX_RATE` | 0.25 | share of a class's cells that may carry incompatible evidence |

Four design points, all fixed in advance and none derived from an observed cluster:
- **Detection, never absence.** Dropout can only *hide* evidence, so a detection-based
  rule under-calls contradictions and can never manufacture one from a zero. This is
  what lets an NK call be challenged by T-lineage evidence **without** requiring a true
  NK cell to be literally TCR-negative, which at 1,162 genes/cell would be unreliable.
- **The contradiction programs are not the identification panels.** A panel that
  *identifies* a class can carry an ambient-prone gene among several; a program that
  *accuses* a cell of another lineage cannot, because a false accusation is created out
  of ambient rather than hidden by dropout. So the globins are excluded from
  `erythroid` and `LYZ` from `myeloid` — the two dominant ambient species here.
- **What has no contradiction is load-bearing too.** PlasmaCell is not contradicted by
  B and vice versa (plasma cells *are* B-lineage; the plasmablast continuum is real),
  and **HSPC has no contradictions at all**, because progenitors legitimately co-express
  lineage-priming programs and flagging that would be a biology error.
- **The fallback must clear both vetoes.** v1's structural flaw was that manual was
  simultaneously the concordance reference and the fallback, so a wrong manual labelling
  could not be escaped: the automated methods "failed" by disagreeing with a bad
  reference, and the class then fell back to that same bad reference. A class whose
  fallback also fails a veto is now reported **unresolved**, not quietly assigned.

Where both qualify, take the higher concordance and record that both passed. Where
neither does, the class falls back to the manual cluster label. **A failed
marker-coverage test vetoes a class regardless of concordance** — high agreement on a
biologically unsupported label is agreement on an error. Expected but **not** assumed:
immune classes and plasma cells from CellTypist, erythroid/HSPC from SingleR or manual.
**The numbers decide.** Outcome to `results/06_annotation/annotation_decision.md`.

**Interface contract — downstream stages read `cell_type` and nothing else:**
- `obs["cell_type"]` — the canonical load-bearing coarse label, seven project classes.
- `obs["cell_type_fine"]` — CellTypist fine label where available, else `NA`. Stage 11
  only; never load-bearing.
- `obs["annotation_source"]` — per cell: `celltypist` | `singler` | `manual`. Without it
  a mixed-provenance label column is untraceable.
- `obs["annotation_conf"]` — the winning method's confidence for that cell.
- `config.ANNOTATION_DECISION` — the per-class method map, so no downstream module
  branches on annotation logic.

This decoupling is deliberate: the comparison can be redone or reversed later without
touching stages 07-12.

**Orthogonal cell-state programs — continuous, never categorical.** A cell has one
`cell_type` but can carry several active programs at once. Score with
`scanpy.tl.score_genes`, store as float `obs` columns, **never collapse into
`cell_type`**:

| Program | Why it matters here |
|---|---|
| Cell cycle (`MKI67`, `TOP2A`, `PCNA`) | a proliferative escape subclone is a different risk from a quiescent one — feeds stage 10 |
| Interferon (`ISG15`, `IFI6`, `STAT1`, `MX1`) | immune-pressure marker; feeds stage 11 |
| Antigen presentation (`B2M`, HLA I/II) | `B2M` loss is a documented myeloma immune-escape route. CAR-T is MHC-independent so it does not affect the escape metric, but it is a *competing* evasion mechanism and belongs in the stage 11/12 interpretation |
| UPR (`XBP1`, `ATF4`, `HSPA5`, `DDIT3`) | plasma cells are professional secretors; UPR tone is core PC biology |
| Hypoxia / stress | standard confounder — cheap to score, expensive to discover late |

A cycling plasma cell is `cell_type == "PlasmaCell"` **plus** a high cell-cycle score,
not a "Cycling" identity. If any method emits a proliferation label as an identity,
remap it to PlasmaCell + score.

**Per-patient composition is a first-class output**, not a by-product: malignant-PC
fraction of the marrow (tumor burden, context for stage 12) and T/NK/myeloid abundance
(stage 11's primary confounder). Any composition *comparison* between groups uses
**`scCODA`**, not a raw proportion test — proportions are compositional (they sum to 1,
so one type rising forces others down) and naive per-type tests are anticonservative.
`scCODA` needs its own env.

---

## Risk tiers are PROVISIONAL (recorded 2026-08-25)

**The Stage-08/09b risk tiers are PROVISIONAL, and the word "final" has been retired from
them (2026-08-25).** They live at `results/08_dual_antigen_escape/risk_tier_provisional/`
and are **measurement-robust provisional tiers under the frozen Stage-08/09b rule** — a
`robust-high` label means the observed DN fraction survived the denominator, depth,
repeated-sample, null-scheme and threshold sensitivity analyses, and **nothing more**.
> **Stage 10 sits between provisional measurement-robust tiering and any final biological
> risk classification.** No provisional tier may be cited as a final biological escape
> classification until stage 10's coherence states are frozen.
Provisional outcome: 4 `robust-high` (`MMRF_1267`, `MMY18273`, `MMY74196`, `MMY98423`),
28 `uncertain`, **0 `robust-low`**, all 32 `THRESHOLD_ROBUST` across
`TAU_HIGH` ∈ {0.20, 0.25, 0.33}. The relabelling changed **no** patient membership, no
threshold, no measurement and not the algorithm — provenance and terminology only.

---

## Stage 10 frozen; the two axes stay separate (2026-08-25)

**Stage 10 is RUN, ACCEPTED and FROZEN exactly as generated (2026-08-25).**
4 `DN_COHERENCE_SUPPORTED` (`MMRF_1720`, `MMRF_2038`, `MMY34339`, `MMY80649`),
23 `DN_COHERENCE_NOT_SUPPORTED`, 5 `DN_COHERENCE_NOT_EVALUABLE`, and
**`CNV_SUBCLONE_NOT_EVALUABLE` for all 32** — the failed stage-07 CNV method was not
reused, and `CNV_SUBCLONE_NOT_SUPPORTED` is never emitted because an underpowered null is
not a negative. Coherence licenses **"escape-associated transcriptional state"**, never
"subclone". Full results: `results/10_dn_coherence/`.

### The two frozen axes stay separate — no composite score, ever

> **Elevated DN burden and coherent DN-associated transcriptional organization are
> separable properties in this cohort.**

| axis | states |
|---|---|
| **measurement** (stage 08/09b) | `measurement_robust_high` · `measurement_uncertain` · `measurement_robust_low` |
| **biological coherence** (stage 10) | `DN_COHERENCE_SUPPORTED` · `DN_COHERENCE_NOT_SUPPORTED` · `DN_COHERENCE_NOT_EVALUABLE` |
| **genetic** | `CNV_SUBCLONE_NOT_EVALUABLE` (all patients) |

**Do not invent a post-hoc rule collapsing these into one scalar escape-risk score.** The
observed cross-tab **is itself a result and must not be "fixed"**: 4 measurement-robust-high
patients, **all** `DN_COHERENCE_NOT_SUPPORTED`; 4 `DN_COHERENCE_SUPPORTED` patients, **all**
measurement-uncertain; **zero patients with convergent high-measurement + coherence-supported
evidence.** The final synthesis displays measurement robustness, DN coherence and
exploratory immune context as **three distinct evidence layers**.

### `MMRF_1640` — the project's primary methodological example

Moran's I **0.47**, unconditioned permutation p **0.001**, depth-stratified p **0.499**.
Use it to illustrate the project-level lesson, not as a patient anecdote:

> **A visually and statistically strong DN-coherence signal can be generated by sequencing
> depth alone, and conditioning on depth can eliminate the apparent biology.**

Its force is that **the same technical confound reproduces across two different
statistical frameworks** — stage 08's count-based co-negativity enrichment (unconditioned
median 1.052 → depth-conditioned 1.009; `MMRF_1505` 4.61 → 1.42) and stage 10's
neighbourhood coherence. Same artifact, different mathematics, same direction — toward the
project's own hypothesis. **No stage-10 result is altered because of this example.**

---

## Stage 10 COMPLETE and FROZEN — three-level hierarchy (2026-08-25)

**Stage 10 is COMPLETE and FROZEN (2026-08-25), through the full three-level evidence
hierarchy.** Results: `results/10_dn_coherence/`.

### The three levels are independent axes and are never combined

| level | states | licenses |
|---|---|---|
| **1 — DN structure** | `DN_STRUCTURE_SUPPORTED` **4** · `NOT_SUPPORTED` **23** · `NOT_EVALUABLE` **5** | *non-random DN organization* — **nothing more** |
| **2 — DN phenotype** | `DN_STATE_SUPPORTED` **26** · `NOT_SUPPORTED` **1** (`MMY67868`) · `NOT_EVALUABLE` **5** | *escape-associated transcriptional state* |
| **3 — genomic** | `CNV_SUBCLONE_NOT_EVALUABLE` **32 (all)** | nothing |

**Level 1 alone is NOT an escape-associated transcriptional state.**
**`CNV_SUBCLONE_NOT_SUPPORTED` is never emitted or implied, and no patient may be called a
genetic escape subclone.**

> **The predeclared per-patient Level-2 criterion proved weakly discriminative, with 26/27
> evaluable patients satisfying it. `DN_STATE_SUPPORTED` therefore indicates compatibility
> with the cohort-level DN-associated program and should not be interpreted as strong
> patient-specific evidence of a distinct escape state.** The rule was **not** retuned
> after seeing that result. **The scientifically useful Level-2 result is the cohort-level
> transcriptional phenotype.**

### The frozen Level-2 biological conclusion

> DN cells exhibit a reproducible **cohort-level** transcriptional shift characterised by
> reduced plasma-cell secretory/differentiation, antigen-presentation, OXPHOS and
> interferon programs. **This supports a biological state associated with observed DN
> status, but does not establish antigen-specific escape.** Patient-level Level-2
> classification is weakly discriminative, and independent genomic subclone evidence is not
> evaluable.

**Keep "biology associated with DN status" and "antigen-specific escape mechanism"
distinct. The latter is not claimed.**

Reproducible DN-**lower** programs (BH < 0.10 under **both** denominators, consistent
sign): **antigen presentation · OXPHOS · interferon**. Not reproducible: **MYC · stress ·
UPR · γ-secretase**. Do not search additional programs.

**Pseudobulk DE — 190 genes significant under both denominators**, paired over patients on
depth-matched cells. Many DN-lower genes are ER/secretory-machinery and mature plasma-cell
identity genes (`SPCS1`, `SPCS2`, `SEC61B`, `UBE2J1`, `TMBIM6`, `MZB1`, `B2M`).
**The DN phenotype is compatible with a less secretory / less differentiated plasma-cell
state.** Never write that DN cells "specifically evolved an antigen-escape program".
**The present data cannot separate a broad plasma-cell-state shift from an antigen-specific
escape mechanism**, and BCMA/GPRC5D are themselves secretory-pathway-dependent surface
proteins, which is exactly why the two are confounded here.

### The pre-registered γ-secretase result is a clean NEGATIVE

Frozen five-gene hypothesis, **no gene added after seeing results**: `NCSTN`, `PSEN1`,
`APH1A`, `APH1B`, `PSENEN`. **Not supported** under the frozen both-denominators rule
(depth-matched BH 0.387 primary / 0.070 sensitivity). **The direction is opposite to the
pre-registered γ-secretase-high prediction.** This stays visible in the stage-10 summary.

### The depth lesson — now the THIRD independent instance

Depth conditioning altered the biological interpretation again at Level 2: matching moved
the DN/comparator depth ratio from **~0.470 to ~0.992**; **γ-secretase looked strongly
negative before matching and not after**; **MYC changed sign**; OXPHOS was the strongest
depth-tracking program tested (ρ 0.274).

> **Three independent places where naive/unconditioned analysis was materially misleading:
> (1) stage 08 co-negativity enrichment, (2) stage 10 Level-1 DN structure, (3) stage 10
> Level-2 transcriptional phenotype.**

`MMRF_1640` remains the project-level example — **Moran's I 0.47, unconditioned p 0.001,
depth-stratified p 0.499**. This is part of the frozen methodological result, **not** a
prompt for another correction or rerun.

**decoupler** (Hallmark/PROGENy/CollecTRI, full tested space retained): E2F targets / G2M /
MAPK up, p53 / OXPHOS / JAK-STAT down. **Limitation, retained: decoupler ULM significance
is based on the fitted contrast vector and is not an independent patient-replicated
hypothesis test.** These p-values may not create patient states or upgrade Level-2 evidence.

**TC-like subtype** stays per-patient, a transcriptional proxy, **descriptive only**, never
a translocation diagnosis, and never an input to structure/state calls. No association
testing.

### Cross-axis — four separate axes, no combined score, ever

Measurement **4 / 28 / 0** · Level-1 **4 / 23 / 5** · Level-2 **26 / 1 / 5** ·
Level-3 **32 not evaluable**.

> **No patient is simultaneously measurement-robust-high, Level-1 structure-supported and
> Level-2 state-supported.**
>
> **The apparent measurement-high × Level-2 overlap of four patients is not strongly
> discriminative, because Level-2 support occurs in 26 of 27 evaluable patients.**

**Do not invent a combination rule or a scalar risk score.**


---

## Stage 11 — exploratory immune context: RUN, ACCEPTED and FROZEN (2026-08-26)

Full frozen record in the main project document; design and its three dated amendments in
`results/11_immune_context/stage11_design.md`; narrative in `stage11_immune_summary.md`.
This section carries the run output.

### Reproducibility — the reason this stage was reopened at all

The first (paused) Stage-11 run left tables on disk and **no driver**. On resume,
`notebooks/11_immune_context.py` was written as that driver, recomputing every table from the
frozen upstream artifacts. The paused run is preserved verbatim in
`results/11_immune_context/preliminary_run/` and the notebook asserts agreement with it:

| check | result |
|---|---|
| `patient_immune_composition.csv`, 27 numeric columns × 32 patients | max drift **1.78e-15** |
| `communication_context.csv`, all 1,768 rows | every row reproduces, max \|Δ lr_score\| 2.9e-3 (CSV rounding) |
| immune category counts | Tcell 61,162 · Myeloid 33,817 · Erythroid 16,224 · Bcell 12,226 · NK_core 8,951 · HSPC 2,879 · cytotoxic_mixed 2,207 |

### Feature diagnostics — run before any DN association was read

| feature | depth coef | depth p | depth BH | cohort BH | MMRF | WU1 | WU2 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Tcell | −0.374 | 0.724 | 0.845 | 0.113 | 0.302 | 0.513 | 0.422 |
| **NK_core** | **−1.987** | **0.027** | 0.093 | **0.040** | 0.032 | 0.091 | 0.058 |
| **cytotoxic_mixed** | **−4.752** | **0.047** | 0.109 | **0.011** | 0.005 | 0.020 | 0.010 |
| Myeloid | +1.681 | 0.100 | 0.175 | 0.989 | 0.236 | 0.213 | 0.198 |
| Bcell | −0.139 | 0.896 | 0.896 | 0.802 | 0.060 | 0.041 | 0.083 |
| **HSPC** | **+3.262** | **0.0009** | **0.007** | **0.004** | 0.050 | 0.009 | 0.032 |
| **Erythroid** | +2.307 | 0.405 | 0.568 | **0.040** | **0.296** | **0.011** | 0.026 |

Three of seven track depth; four differ by cohort; the erythroid gap is **27×** and is a
cell-recovery property of the deposit, not marrow biology.

### Q1–Q3 — 4 predictors × 7 features = 28 tests

**0 at BH < 0.10.** Two reach raw p < 0.05 adjusted, both `NK_core`, both BH ≈ 0.49:

| predictor | coef | 95% CI | p | BH | ρ MMRF / WU1 / WU2 |
|---|---:|---|---:|---:|---|
| `obs_dn_primary` | +1.702 | 0.045 – 3.359 | 0.045 | 0.491 | +0.39 / +0.45 / +0.51 |
| `obs_dn_sensitivity` | +2.007 | 0.104 – 3.910 | 0.040 | 0.491 | +0.14 / +0.35 / +0.63 |

**Fragility check (required by the frozen design):** as a plain log-fraction, `NK_core` alone
gives coef +2.485 p = **0.164**, and `NK_core + cytotoxic_mixed` +2.246 p = **0.211**. The
association does not survive a change in how the mixed cytotoxic cells are handled.

**Guard that worked:** `NK_core` × `enr_cohortbins`, unadjusted p = 0.019 → adjusted p = 0.070.

**Coherence axis (27 patients, 4 supported): nothing.** Smallest is `Bcell` p = 0.077,
BH 0.27. Medians move (`NK_core` 0.023 supported vs 0.075 not) but 4 vs 23 is not evidence.

### Repeated samples — within-patient range of each lineage fraction

| patient | Tcell | NK_core | cyto_mixed | Myeloid | Bcell | HSPC | Erythroid |
|---|---:|---:|---:|---:|---:|---:|---:|
| 27522 | 0.348 | 0.047 | 0.010 | 0.315 | 0.127 | 0.019 | **0.504** |
| 47491 | 0.298 | 0.044 | 0.012 | 0.249 | 0.011 | 0.010 | 0.006 |
| 56203 | 0.216 | 0.023 | 0.036 | 0.093 | 0.068 | 0.011 | 0.079 |
| 58408 | **0.468** | 0.003 | 0.000 | **0.438** | 0.030 | 0.001 | 0.001 |
| 59114 | 0.034 | 0.003 | 0.001 | 0.187 | 0.136 | 0.015 | 0.005 |
| 81012 | 0.151 | 0.016 | 0.007 | 0.103 | 0.020 | 0.027 | 0.008 |
| 83942 | 0.035 | 0.008 | 0.045 | 0.043 | 0.018 | 0.008 | 0.016 |

`83942` — one patient, two protocols — is the stable case. Elsewhere within-patient variation
exceeds most between-patient differences.

### Q4 communication — both receivers, and the confound that explains both

17 predeclared pairs × 4 senders = 68 tests; senders gated at `MIN_SENDER_CELLS` = 20;
31 patients, 1,768 rows under each receiver definition.

Highest median LR scores (design-conformant receiver): `Myeloid TNFSF13B → TNFRSF13B` 13,414 ·
`NK_core GZMA → F2R` 4,052 · `cytotoxic_mixed GZMA → F2R` 2,034 · `NK_core IFNG → IFNGR2` 1,171.

| receiver | raw p < 0.05 | BH < 0.10 | `Tcell PDCD1 → CD274` |
|---|---:|---:|---|
| all plasma (paused first run, preserved) | 3 | 1 | coef −3.041, p 0.00066, BH 0.043 |
| **clone plasma (design-conformant)** | 11 | **5** | coef −4.348, p 0.00004, **BH 0.0026** |

The five BH < 0.10 hits under the design-conformant receiver are `PDCD1→CD274`,
`Tcell TNFSF10→TNFRSF10A` (−7.08), `Tcell TNFSF10→TNFRSF10B` (−5.45),
`Myeloid TNFSF10→TNFRSF10B` (−6.05), `Myeloid TNFSF10→TNFRSF10A` (−7.70) — **all negative,
across unrelated receptors.**

**Receiver-side test (Amendment 2), each receptor alone under the identical adjusted model:**

| receiver | receptors negative vs DN | raw p < 0.05 |
|---|---:|---:|
| all plasma | 8 / 15 | 2 |
| **clone plasma** | **11 / 15** | **5** |

`MICA` p = 0.00057 · `CD274` p = 0.0020 · `MICB` p = 0.0023 · `TNFRSF10A` p = 0.014 ·
`CD86` p = 0.035.

> **The communication panel reads a broad plasma-cell-state shift on its receptor term, not a
> specific immune axis.** `Tcell PDCD1 → CD274` is **not** reported as a candidate
> immune-evasion axis.

### Module additions

`communication.stream_gene_counts` block-reads named columns out of a CSR layer without
materialising the 172,940 × 32,991 matrix (~20 GB); `communication.pseudobulk_cpm` pools rather
than averaging per-cell rates, because a mean would weight a 300-UMI WashU cell like a
20,000-UMI MMRF one. `LR_SENDERS`, `MIN_SENDER_CELLS` (= the frozen `MIN_GROUP_CELLS`, 20).
`mm_escape` is now editable-installed in `mm-communication`.


---

## Stage 08c — supplemental multi-antigen coverage: COMPLETE / ACCEPTED / FROZEN (2026-08-26)

Full frozen record in the main project document; design in
`results/08_dual_antigen_escape/multi_antigen_coverage/multi_antigen_design.md`; narrative in
`multi_antigen_coverage_summary.md`. This section carries the run output.

**Supplemental Stage-08 deliverable consuming frozen upstream infrastructure.** No frozen
artifact was modified — verified by filesystem timestamp and by
`frozen_upstream_digests.json`, which `tests/test_coverage.py` re-checks.

### Per-target measurement QC

Detection fraction in the primary denominator (21,906 cells, 32 patients):

| target | MMRF | WU1 | WU2 | ALL | mean count | depth ρ | stratum spread | technical-zero (pooled) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `TNFRSF17` | 0.835 | 0.831 | 0.667 | 0.772 | 5.77 | +0.446 | 3.1× | **0.276** |
| `GPRC5D` | 0.589 | 0.129 | 0.062 | 0.216 | 0.63 | +0.480 | **15.8×** | **0.620** |
| `SLAMF7` | 0.712 | 0.427 | 0.531 | 0.534 | 1.52 | +0.436 | 6.4× | 0.429 |
| `FCRL5` | 0.697 | 0.575 | 0.549 | 0.595 | 2.06 | +0.419 | 6.9× | 0.425 |
| `CD38` | 0.757 | 0.383 | 0.411 | 0.484 | 1.66 | +0.461 | 6.0× | 0.416 |
| `SDC1` | 0.799 | 0.487 | 0.516 | 0.574 | 2.41 | +0.448 | 4.3× | 0.409 |
| `ITGB7` | 0.567 | 0.245 | 0.401 | 0.380 | 1.55 | +0.318 | **12.8×** | 0.449 |

**No target is `comparatively_reliable`; every eligible one is `depth_sensitive`.**

Background separation: `TNFRSF17` 27.9× · `GPRC5D` 151× · `FCRL5` 50.7× · `SDC1` 42.2× ·
`SLAMF7` 22.9×* · `CD38` 8.2×* · `ITGB7` 4.6×  (*ambient floor `NOT_EVALUABLE` — erythroid
only, no clean marrow negative exists for these targets).

**The technical-zero floor is computed on one consistent rule across all seven and is not
bit-identical to Stage 08's.** Denominator gene means reproduce exactly (`TNFRSF17` 5.7712,
`GPRC5D` 0.6287) but this construction runs 0.10–0.20 higher (more conservative); Stage 08's
own rows stay the cited values for the anchor pair. The rank order agrees with the frozen
result. Threshold applied as the **pooled median over all cohort × stratum rows** because the
design's "median cohort stratum" wording did not name a cohort; `SLAMF7` (WU2 0.503) and
`ITGB7` (WU1 0.514 / WU2 0.518) are flagged threshold-hugging.

### The SDC1 checkpoint

| target | secretory-breadth delta (neg − pos) | ρ(count, secretory) | in `PLASMA_MATURE` | in `MARKER_PANEL[PlasmaCell]` |
|---|---:|---:|---|---|
| `TNFRSF17` | **−0.036** | **+0.228** | **yes** | no |
| `GPRC5D` | −0.016 | +0.227 | no | no |
| `SLAMF7` | −0.015 | +0.028 | no | no |
| `FCRL5` | −0.004 | −0.106 | no | no |
| `CD38` | −0.013 | +0.157 | no | **yes** |
| `SDC1` | −0.013 | **+0.020** | **yes** | **yes** |
| `ITGB7` | −0.028 | +0.215 | no | no |

**SDC1's differentiation confound is measurable and it is not the worst offender** — the
anchor gene `TNFRSF17` is. **SDC1 is excluded on circularity instead:** `PLASMA_MATURE` is
`("SDC1", "TNFRSF17")`, Stage 06's axis-(b) predicate, so the plasma denominator was partly
established using SDC1 detection. **`TNFRSF17` carries the identical limitation and is
retained only as the frozen anchor — a disclosure, not a distinction.**

### Eligibility

`COVERAGE_ELIGIBLE`: `TNFRSF17`, `SLAMF7`, `FCRL5`, `CD38`, `ITGB7`.
`COVERAGE_NOT_EVALUABLE`: **`GPRC5D`** (technical-zero 0.62 ≥ 0.50) and **`SDC1`**
(circularity). Both remain visible in `target_measurement_qc.csv` with reasons, and their
combinations are computed and reported as `exploratory`.

### Coverage — 4,032 rows (32 patients × 2 denominators × 63 combinations), 0 monotonicity violations

**Anchor, primary denominator:**

| quantity | median | mean | min | max |
|---|---:|---:|---:|---:|
| uncovered by BCMA | 0.353 | 0.327 | 0.024 | 0.800 |
| uncovered by GPRC5D | 0.899 | 0.821 | 0.155 | 1.000 |
| uncovered by BCMA+GPRC5D | **0.335** | 0.308 | 0.017 | 0.783 |
| gain from adding GPRC5D | **0.011** | 0.019 | 0.000 | 0.152 |
| gain from adding BCMA | 0.547 | 0.513 | 0.138 | 0.858 |

24 of 32 patients show a GPRC5D gain below 2 points. **Not evidence of clinical redundancy**
— GPRC5D is detected in 21.6% of cells with a 0.62 technical-zero floor, which is why it
failed QC before any gain was computed.

**Median uncovered fraction (primary), greatest observed coverage first.**
Singles: `TNFRSF17` 0.353 · `SDC1` 0.454* · `FCRL5` 0.468 · `SLAMF7` 0.518 · `CD38` 0.575 ·
`ITGB7` 0.762 · `GPRC5D` 0.899*.
Pairs: `TNFRSF17+SDC1` 0.173* · `TNFRSF17+SLAMF7` 0.214 · `TNFRSF17+FCRL5` 0.217 ·
`TNFRSF17+ITGB7` 0.242 · `TNFRSF17+CD38` 0.250 · … `TNFRSF17+GPRC5D` 0.335*.
Triples: `TNFRSF17+FCRL5+SDC1` 0.101* · `TNFRSF17+SLAMF7+SDC1` 0.110* ·
`TNFRSF17+SLAMF7+FCRL5` 0.129 · `TNFRSF17+FCRL5+ITGB7` 0.141.  (*exploratory)

Every eligible alternative pair beats the anchor in **32/32 patients** (median advantage
0.098, max 0.361); per-patient greatest-coverage eligible pair is `TNFRSF17+SLAMF7` (11),
`TNFRSF17+FCRL5` (11), `TNFRSF17+ITGB7` (5), `TNFRSF17+CD38` (5). **A detection-rate
artifact, not a therapeutic finding.** Best eligible pair 0.183 → best eligible triple 0.119.

### Robustness

- **Truncate-10k:** WashU **exactly 0.000** change (deposit ceiling confirmed on the full
  panel); MMRF mean +0.052 anchor / +0.070 primary matrix, max +0.357. Ordering ρ = **0.996**.
- **Primary vs sensitivity:** anchor median +0.032, mean +0.045, max +0.256; **12/32 patients
  move >5 points**. Never collapsed.
- **Repeated samples:** anchor uncovered range within patient — `27522` **0.571** (6
  timepoints), `59114` **0.547**, then 0.098, 0.059, 0.045, 0.016, 0.015. Two patients vary
  more within themselves than most patients differ from each other.

### Pairwise co-loss — the depth lesson generalises

Synthetic depth-only control re-run first: **1.323 unconditioned → 0.996 conditioned.**
On real data **all 21 pairs collapse to ~1.0** under depth conditioning (unconditioned
1.01–1.21 → conditioned 1.00–1.06); `TNFRSF17+GPRC5D` 1.066 → 1.009, reproducing the frozen
Stage-08 value. **The Stage-08 negative co-negativity result holds across the entire panel.**

### Normal-marrow expression context (8 donors, donor as unit, median detection)

| lineage | TNFRSF17 | GPRC5D | SLAMF7 | FCRL5 | CD38 | SDC1 | ITGB7 |
|---|---:|---:|---:|---:|---:|---:|---:|
| PlasmaCell | 0.377 | 0.009 | 0.165 | 0.335 | **0.656** | 0.308 | 0.091 |
| Tcell | 0.002 | 0.000 | 0.026 | 0.002 | 0.017 | 0.000 | 0.166 |
| NK | 0.000 | 0.000 | 0.180 | 0.001 | 0.168 | 0.000 | 0.166 |
| Bcell | 0.028 | 0.000 | 0.005 | 0.132 | 0.133 | 0.000 | 0.054 |
| Myeloid | 0.011 | 0.000 | 0.050 | 0.010 | 0.043 | 0.000 | 0.076 |
| Erythroid | 0.001 | 0.000 | 0.004 | 0.000 | 0.006 | 0.000 | 0.007 |
| HSPC | 0.000 | 0.000 | 0.000 | 0.000 | **0.353** | 0.000 | 0.017 |

**Expression context only — never safety.** GPRC5D's decisive liability is keratinized
tissue, unobservable here, and this table would misleadingly make it look safest.

### Module and tests

`src/mm_escape/coverage.py` — coverage algebra, per-target QC labelling and eligibility, and
**thin delegations** to `antigen.py` for every depth operation. It contains **no binning code**;
`test_j_module_contains_no_local_depth_binning` scans its source for `np.quantile`,
`searchsorted`, `pd.cut` and friends. The primary analysis reads Stage 08's frozen per-cell
`depth_stratum_cohort`, so no perturbation of any target can move a cell between strata.

`tests/test_coverage.py` — 30 tests covering monotonicity (A), incremental gain (B),
denominator isolation (C), antigen isolation (D), eligibility (E), raw counts (F), no utility
score (G, AST-based with docstrings stripped), patient unit (H), frozen-stage isolation (I,
digest-based), shared depth utility (J), the synthetic depth-only null (K, including a
positive control that a genuine co-loss is still detected) and the SDC1 safeguard (L).


---

## Stage 11b — LIANA verification arm: COMPLETE / ACCEPTED / FROZEN (2026-08-26)

**Exploratory · post hoc · non-tier-changing · non-classifying.** A dated addendum to
Stage 11, kept **explicitly separate** from the original custom Stage-11 communication
analysis, which is not reopened, replaced or rewritten — and which **was not performed with
LIANA**. Full narrative: `results/11_immune_context/liana_verification/liana_summary.md`.

> **Accepted interpretation (final wording, not to be strengthened).** Stage 11 found no
> robust independent evidence that immune composition or ligand–receptor communication
> explains the observed DN phenotype. The targeted LR analysis was receiver-state confounded,
> and LIANA verification did not rescue that interpretation. The strongest LIANA consensus
> association was structurally circular because its receptor was `TNFRSF17`, one of the
> antigens defining the DN predictor. **LIANA is not a validation of immune evasion.**

### Environment and configuration

`liana 1.8.1` · Python 3.12.13 · scanpy 1.12.3 · anndata 0.12.19 · numpy 2.0.2 ·
pandas 2.3.3 · scipy 1.18.0 · omnipath 1.0.12. LIANA was already installed; nothing was
downgraded and no unrelated package touched. The **installed** API was inspected rather than
a tutorial assumed.

| item | value |
|---|---|
| primary | `rank_aggregate` — RRA consensus of CellPhoneDB, Connectome, log2FC, NATMI, SingleCellSignalR |
| continuity | `cellchat`, run separately and preserved |
| resource | `consensus` (4,624 interactions), fixed for both methods |
| per-patient | `Method.by_sample(sample_key="patient_id")` |
| params | `expr_prop` 0.1, `n_perms` 1000, `seed` 1337 (LIANA defaults); `min_cells` 20 (frozen `MIN_SENDER_CELLS`) |
| orientation | consensus `magnitude_rank` lower = stronger → `-log10`; cellchat `lr_probs` higher = stronger → identity |

### Evaluability

≥1 sender category **and** the receiver, each ≥20 cells. **31 of 32 patients evaluable**;
`25183` is `LIANA_NOT_EVALUABLE` in all three runs because it contributes **zero immune
sender cells**, so no sender→receiver pair exists. **This is a biological/evaluability
exclusion, not a software failure**, and it is consistent with the frozen Stage-11
communication design — `25183` is likewise the one patient absent from the frozen custom
communication table.

### Tested space

| run | interactions/patient | patients | tested ≥20 pts | raw p<0.05 | BH<0.10 |
|---|---:|---:|---:|---:|---:|
| consensus, clone-primary | 1,050 | 31 | 87 | 12 | **1** |
| cellchat, clone-primary | 1,050 | 31 | 87 | 10 | **5** |
| consensus, all-plasma | 1,001 | 31 | 94 | 10 | **1** |

Frozen Stage-11 confound model used unchanged; unadjusted reported beside adjusted; BH over
the full space; the whole tested space written out, not only significant rows.

### The one consensus hit is structurally circular

`Myeloid TNFSF13B → TNFRSF17` — consensus coef **−1.284**, p 0.0001, **BH 0.0054**;
CellChat BH 0.0104; all-plasma BH 0.0271. **Its receptor is `TNFRSF17`**, one of the two
antigens whose negativity defines `obs_dn_primary`, so higher DN forces lower `TNFRSF17` in
the same cells — the negative coefficient is arithmetic. Receptor-side alone reproduces 94%
of it (coef −1.203, p 0.020).

**This association is not interpretable as evidence of immune communication driving the DN
phenotype.** 100% of consensus BH<0.10 hits are antigen-circular. The frozen targeted panel
had deliberately paired BAFF with `TNFRSF13B`, never `TNFRSF17`; LIANA's resource contains
`TNFSF13B → TNFRSF17` and an unrestricted screen cannot know it must not use it.

**The row is flagged and preserved in the raw LIANA output, never excluded**, and the
**antigen-circularity test is permanent** — any interaction touching `TNFRSF17` or `GPRC5D` is
structurally circular for the DN-vs-communication question, asserted by
`tests/test_liana_verification.py`.

### Mandatory receiver-side decomposition

All 12 raw-significant interactions decomposed into ligand and receptor halves:
**3 `RECEIVER_STATE_CONFOUNDED`** (`TNFSF13B→TNFRSF17`, `CALM1→PTPRA` receptor −6.95
p<1e-4, `VCAN→ITGA4` receptor −3.69 p 0.017), **9 `NOT_REPRODUCED_BY_LIANA`**,
**0 `EXPLORATORY_LIANA_ONLY`**. **2 `ABUNDANCE_SENSITIVE`** (`TOR2A→ATP5F1B` p 0.021→0.060;
`Tcell HMGB1→CXCR4` p 0.034→0.057).

Pooled Stage-10 receiver-state test, reported honestly as weak: receptor pool vs antigen
presentation ρ −0.349 (p 0.050); secretory/ER −0.157, OXPHOS +0.122, depth +0.171, all n.s.

### Frozen 17-pair cross-check

**8 of 17 pairs are in LIANA's consensus resource.** Absent: `PDCD1→CD274`, `CTLA4→CD86`,
`GZMK→F2R`, `IFNG→IFNGR1/2`, `IL6→IL6R`, `KLRK1→MICA/MICB/ULBP2`. Of the 8 present, only
**one row** (`Myeloid TNFSF13B→TNFRSF13B`, 26 patients) reached the 20-patient floor — the
rest were evaluated in 1–6 patients because LIANA's own `expr_prop` filter drops the sparsely
expressed TRAIL and granzyme interactions. That row: LIANA −0.948 p 0.0052 BH 0.150 vs custom
−1.436 p 0.549 — **directions agree, significance does not**, conservatively labelled
`NOT_REPRODUCED_BY_LIANA`. **0 `CONSISTENT_WITH_TARGETED_PANEL`.**

> **`PDCD1 → CD274` = `NOT_EVALUABLE_BY_LIANA_RESOURCE`** — absent from the `consensus`
> resource (`in_liana_resource = False`, `panel_status = NOT_EVALUABLE`). **LIANA did not
> disprove it and did not support it.** The frozen custom Stage-11 conclusion remains
> operative: receiver-state confounded, not accepted as an immune-evasion axis.

**LIANA is a partial methodological verification arm, not a full reproduction of the targeted
17-pair panel**, and missing resource coverage is **not** a negative biological result.

### Method agreement

Consensus and CellChat share only 3 raw-significant interactions of 12 and 10. CellChat's
`MIF→TNFRSF14` reaches BH<0.10 in three senders but **flips sign between cohorts** in all
three, as does `TOR2A→ATP5F1B`. Only the antigen-circular hit is cohort-consistent.
**CellChat-specific hits are not promoted merely because there are more of them; the method
disagreement stays visible and is itself a reason for caution.**

### Conclusion

**No LIANA-only interaction remains credible after the controls, and LIANA does not change
the frozen Stage-11 interpretation.** It reaches the same negative conclusion by a different
route and adds one new observation: an unrestricted LR resource reintroduces an
antigen-circular interaction and ranks it first. **No patient tier or state changed; no
LIANA-specific classifier exists.** `tests/test_liana_verification.py` — 18 tests (A–J),
including digest checks on five frozen state files.

---

## Stage 12 — final synthesis: RUN, COMPLETE (2026-08-27)

This file's scope is deliberately stages 01 through 11b, per its title, and its content
above is unchanged by what follows. **Stage 12 subsequently consumed every result
recorded in this file** — 29 frozen artifacts, hash-verified against
`provenance/frozen_artifacts_pre_stage12.tsv` before use — and closed Phase 1. It added
no new statistical test, fitted no model, and changed nothing recorded above. Full
outcome: `results/12_final_synthesis/stage12_summary.md`; design it executed against:
`docs/stage12_design.md`; condensed record: the main project document's Stage 12 block.

Headline: observed transcript-level BCMA/GPRC5D double-negativity is common at baseline
(median 0.335 across 32 patients), but is dominated by measurement limitations rather
than demonstrable biology, and the measurement-robust, structurally-supported and
phenotype-supported patient sets are disjoint — no patient is convergent across all
three. Stage 12 therefore produced a six-axis per-patient evidence matrix rather than a
risk tier or ranking, which is itself part of the result.
