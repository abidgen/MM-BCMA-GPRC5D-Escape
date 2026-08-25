# Stage results — run output for stages 01 through 05b

Split out of `CLAUDE.md` on 2026-08-24. `CLAUDE.md` keeps, per completed stage, a
short block saying what it produced and what binds downstream; this file is the full
record.

**The `results/*.csv` files are the source of truth, not the tables reproduced here.**
Where a table below duplicates a committed artifact the path is named; if the two ever
disagree, the CSV is right and this file is stale.

**Corrected 2026-08-24 — median genes/cell.** Earlier revisions of `CLAUDE.md` quoted
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

One-time verification records moved here from the `CLAUDE.md` status section on
2026-08-24. Not standing instructions.

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

Gates: `conftest.requires_data` (needs `raw/`), `requires_s1` (Supplementary Table S1
is a journal file, not part of the GEO deposit), `requires_r` (scDblFinder lives only
in `mm-qc`; the suite's home `mm-core` carries no R). `pytest -m "not slow"` skips the
two full-cohort passes.
