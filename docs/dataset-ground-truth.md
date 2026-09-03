# Dataset ground truth — GSE223060 / GSE223061

Split out of the main project document on 2026-08-24. That document keeps the short list
of rules that bind code; this file is the evidence behind them — archive forensics, the
Ensembl-ID reconstruction, the GEO series metadata, the bulk inventory, and the
Supplementary Table S1 patient mapping.

**Everything here was verified against the real files.** Do not re-guess it, and do
not re-derive it from the paper text — where the deposit's own metadata disagrees
with the files, the files win (see "The deposit's own processing metadata is
unreliable").

---


**Source:** GSE223060 (scRNA-seq) / GSE223061 (matched bulk RNA-seq — **no longer
"optional"; it is now the stage 09 orthogonal validation set**), NCBI GEO. Paper:
"Single-Cell Discovery and Multiomic Characterization of Therapeutic Targets in
Multiple Myeloma" (Cancer Research, 2023). 53 BM samples from 41 patients per the
paper's main analysis (MMRF Immune Atlas Pilot + 2 WashU cohorts); the GEO archive
itself bundles 62 sample entries (includes normal/control BM beyond the 53 disease
samples).

**No raw FASTQ/SRA exists for this series.** Confirmed directly against NCBI. Only
processed Cell Ranger *filtered* output is public — not raw, not unfiltered
(matters for ambient-RNA correction, see below). Do not attempt to find or process
raw reads for this dataset; this is a closed question.

**Archive structure (confirmed by direct extraction):**
```
GSE223060_RAW.tar                              # outer archive, ~970 MB
  |-- <GSM_ID>_<sample_name>.tar.gz             # one per-sample archive (62 total)
        |-- (extract to) raw/samples/<GSM_ID>_<sample_name>/
              |-- <sample_name>/                # ONE EXTRA NESTING LEVEL
                    |-- barcodes.tsv             # uncompressed, no .gz
                    |-- genes.tsv                # single column, gene symbols only
                    |-- counts.mtx               # NOT "matrix.mtx"
```

**Critical implications for the Python loader:**
- Filenames are old Cell Ranger v2-style (`counts.mtx`, `genes.tsv`), not the
  `matrix.mtx`/`features.tsv` that `scanpy.read_10x_mtx()` expects — that function
  will silently fail to find these files. Load explicitly with
  `scanpy.read_mtx()` on `counts.mtx`, then attach `.obs`/`.var` from
  `barcodes.tsv`/`genes.tsv` manually, rather than relying on any 10x-directory
  auto-detection helper.
- `genes.tsv` is a **single column of gene symbols**, not the usual 2-3 column
  Ensembl+symbol file — don't write code that assumes a `gene_id` column exists.
- Data sits one directory level deeper than the GSM-named folder.

Sample naming observed: mix of `MMRF_####` (MMRF Immune Atlas Pilot cohort), bare
numeric IDs with `_N` suffixes (e.g. `27522_1` through `27522_6` — likely multiple
timepoints/fractions per patient), `MMY#####` IDs, `ND_######`, and `BM#` (likely
normal bone marrow controls). **Disease stage/cytogenetic metadata per sample is
NOT in the filenames** — pull it from the paper's Supplementary Table S1 before
building the analysis manifest.

**CRITICAL — mixed Cell Ranger references across samples (confirmed):** The 62
samples were processed against **two** references, distinguishable by row count in
`genes.tsv`:
- 33538 genes — 37 samples
- 33694 genes — 25 samples (including `56203_1`, see below)

**`56203_1` is NOT a third reference — it is a truncated deposit, and it is repaired
rather than excluded (corrected 2026-08-24).** The earlier reading (a "22184-gene
build missing `TNFRSF17`") was a misdiagnosis, and the exclusion decision rested on
it. What the files actually show:

```
counts.mtx header:      33694 1837 2135520      <- a normal 33694-build matrix
genes.tsv:              22185 rows, ends 'KBTBD', NO trailing newline
33694 reference row 22185:      'KBTBD7'
rows 1..22184 vs the reference: identical, a strict prefix
```

The gene-file write failed part-way through. `TNFRSF17` (canonical row 25539) and
`IGLC1/2/3` (rows 32548-32552) were never absent from a reference — they were past
the cut. `GPRC5D` (row 20472), `SLAMF7`, `FCRL5`, `SDC1` and `CD38` were all present
even in the truncated file.

The "22184 genes" figure was itself a **`wc -l` artifact**: the file has no trailing
newline, so `wc -l` undercounts by one. There are 22185 written rows.

**Decision: repair and retain.** `io.read_sample` substitutes the canonical column for
the declared build — taken from the committed, position-verified gene map, not from
another sample's file — after asserting that every written row matches the reference
and that the final partial row is a prefix of the symbol it was cut from. If either
assertion fails the load raises rather than substituting, so this is a provable repair
and not a guess. Recovers 1,837 cells and a second sample for patient 56203, which
matters if the `_N` suffixes are serial timepoints (see below). `config.EXCLUDED_SAMPLES`
is now empty; the mechanism lives in `config.TRUNCATED_GENE_FILES`.

**CRITICAL — the two references also use different HGNC symbol vintages, so a naive
symbol intersection silently drops genes present in BOTH builds.** Confirmed
2026-08-20 by direct comparison, and caught by the stage-03 notebook's required-gene
assertions rather than by inspection — which is the argument for keeping those
assertions:

| 33538 (newer symbols) | 33694 (older symbols) | consequence if unharmonized |
|---|---|---|
| `NSD2`    | `WHSC1`   | **t(4;14) becomes uncallable** — the highest-risk MM translocation, and the `MS` class of the TC subgrouping at stage 10 |
| `TENT5C`  | `FAM46C`  | loses a recurrently-deleted MM tumour suppressor (1p12) |
| `NSD3`    | `WHSC1L1` | a **different gene** from NSD2 — do not conflate the two |
| `ATP5F1A` | `ATP5A1`  | OXPHOS program member (stage 10) |

**Decision (superseded 2026-08-21 — see "Ensembl-ID reconstruction" below): join on
Ensembl ID, not on symbols at all.** Symbol canonicalization recovered these four
(22,164 → 22,168 genes); the ID join recovers 32,991 and makes the alias map a
regression assertion rather than the mechanism. The table above still matters — it is
what the assertions check and why they exist.
The four above are the ones this project's gene lists depend on; a fuller HGNC-based
reconciliation is worth doing in `gene_space.py` since the same drift certainly affects
genes outside the required set. **Never assume a missing required gene is biologically
absent — check for a legacy symbol first.**

Note also that the ~11.4k/11.5k symbols unique to each build are dominated by
**annotation-version noise** — versioned clone identifiers such as `AC000032.1` vs
`AC000032.2` — not genuinely absent genes. So 22,164 *understates* the recoverable gene
space. That is tolerable for lncRNA/clone entries (HVG selection never reaches them) but
was not tolerable for `NSD2`, which is exactly why the assertions exist.

### Ensembl-ID reconstruction — SOLVED and verified 2026-08-21

The deposit has no ID column: `genes.tsv` is symbol-only and contains **zero `ENSG`
strings across all 62 samples**. Do not go looking for one. But the IDs are fully
**reconstructible**, and this was verified end-to-end — it is the method, not a plan.

There are only **three distinct `genes.tsv` files in the cohort** (checksum-verified,
byte-identical within group), each a positional dump of a public reference:

| rows | samples | reference | source GTF |
|---|---|---|---|
| 33538 | 37 | `refdata-cellranger-GRCh38-3.0.0` | Ensembl 93 |
| 33694 | 24 | `refdata-cellranger-GRCh38-1.2.0` | Ensembl 84 |
| 22185 | 1 | (a truncated 33694 write) | `56203_1`, **repaired and retained** — see the truncation section above |

**Reconstruction recipe (reproduces both files exactly):**
1. Fetch the Ensembl GTF (release 93 / 84) — ~41-44 MB each, not the 11 GB 10x tarball.
2. Take `feature == "gene"` rows **in GTF order**, keep first occurrence per `gene_id`,
   filtered to 10x's `mkgtf` biotype list (`protein_coding`, `lincRNA`, `antisense`,
   the 8 `IG_*` and 6 `TR_*` classes). This yields exactly 33538 / 33694 rows.
3. Apply the depositor's two transforms, in order: **`gsub("_", "-")`** then
   **R `make.unique`**. Both are Seurat artifacts — the deposited files were written
   through Seurat, which is what produced the `TBCE`/`TBCE.1` suffixes *and* the
   `RP11-442N24--B` (Ensembl writes `__`) spellings.
4. **Assert the result equals the deposited column position-for-position.** This is
   what makes the reconstruction self-certifying: a wrong biotype filter or a wrong
   release changes the row count or the order, and the assertion fails loudly.

**Result: 0 mismatches / 33538 rows and 0 mismatches / 33694 rows.** The mapping is
committed at `resources/gene_space/` (~1 MB gzipped, three TSVs: per-build
`row_index → deposited_symbol → ensembl_id`, plus the intersection table). It does not
need regenerating; the GTFs do not need re-downloading.

**The payoff is large and was not anticipated by either review:**

| join key | genes retained |
|---|---|
| raw symbols | 22,164 |
| symbols + 4-gene alias map | 22,168 |
| **Ensembl IDs** | **32,991** |

**+10,827 genes — a ~49% larger gene space**, because 11,140 intersected IDs carry a
*different symbol* in each build and were silently invisible to a symbol join. The
earlier note that the unique-symbol remainder was "dominated by annotation-version
noise" was too optimistic: most of it is real, recoverable genes.

**It also proves the mis-pairing risk was real, not theoretical:**

```
33538 build:  TBCE -> ENSG00000285053
33694 build:  TBCE -> ENSG00000116957     <- a DIFFERENT annotation entry
```

A symbol join silently merges those two rows. Nine canonical symbols remain duplicated
*within* the ID intersection (`COG8`, `CYB561D2`, `EMG1`, `LINC01238`, `LINC01505`,
`MATR3`, `PINX1`, `RGS5`, `TMSB15B`) — disambiguate as `SYMBOL__ENSGxxxxxxxxxxx`, never
by dropping or by `var_names_make_unique()`, which would re-introduce exactly the
positional mangling this section undoes.

**Implemented and tested 2026-08-21** in `src/mm_escape/gene_space.py` (with a partial
`config.py` carrying only the gene-space constants; the rest is filled in by the stages
that derive it). Public API: `attach_ensembl_ids` (per sample, pre-concat, verifies
position-for-position) -> `intersect_gene_space` (IDs only, refuses symbol-keyed input)
-> `to_canonical_symbols` (once, post-merge) -> `assert_required_genes`.
`rebuild_gene_map_from_gtf` regenerates the committed map and re-verifies it. Verified
against the real files: 33538+33694 -> 32,991 genes, 11,140 drifted symbols joined
correctly, all 64 required genes present, and four failure paths confirmed to raise —
reordered `var_names`, symbol-keyed intersection, a missing required gene (with the
legacy-symbol hint), and `var_names_make_unique()`-style mangling.

**Decision — `gene_space.py` joins on Ensembl ID; the alias map is demoted to an
assertion.** All 21 spot-checked required genes survive the ID intersection. The
four-gene alias dictionary (`WHSC1`→`NSD2` etc.) is kept **only** as a regression
assertion, not as the harmonization mechanism — it addressed 4 of 11,140 drifted
symbols. The earlier "drop the ~52 `make.unique`-ambiguous symbols" interim is
**superseded and must not be implemented**: those genes are now resolved correctly
rather than discarded.

**Index convention.** Use `ensembl_id` as `var_names` **through the merge only** — that
is where identity is load-bearing. Once the object is a single harmonized matrix the
mis-pairing risk is gone, so switch `var_names` to the canonical Ensembl-93 symbol
(with the 9 collisions suffixed as above) and retain `var["ensembl_id"]`,
`var["symbol_33538"]`, `var["symbol_33694"]`. Reason: every downstream consumer is
symbol-native — `score_genes` marker panels, dotplots, `celltypist`, `decoupler`,
`liana` — and keeping IDs as the index buys no further correctness while making every
figure unreadable and every gene-set call a translation step.

**AnnData, not MuData.** MuData exists for multiple modalities over a shared `obs`
axis (RNA + ADT + ATAC). This project has one cell-level modality. The stage-09 bulk
RNA-seq is *sample*-level, so it shares no `obs` axis with the single-cell object and
belongs as a plain DataFrame joined on `sample_id`, not as a MuData modality. Adding
MuData here buys nothing and costs compatibility across the scanpy stack.

**The 33538- and 33694-gene reference sets share only 22164 genes** *on raw symbols*
(22,168 after harmonization). A union merge
across retained samples would make ~11k genes structurally zero in whole sample
cohorts — indistinguishable downstream from a true biological zero, which is
exactly the quantity this project measures. **Decision: intersect gene sets across
retained samples, never union**, before any concatenation (`anndata.concat` with
`join="inner"`, not `join="outer"`). Assert that all markers needed downstream
survive the intersection — BCMA (`TNFRSF17`), GPRC5D, backup antigens (`SLAMF7`,
`FCRL5`), the full annotation marker panel, and the light-chain genes (`IGKC`,
`IGLC1-7`) — hard-fail with the specific missing gene name(s) if any assertion
fails, rather than silently proceeding on a partial marker set.

**Verified 2026-08-20 — every target gene survives the intersection.** Spot-checked
directly against `genes.tsv` in a 33538-gene sample (`MMRF_1695`), a 33694-gene
sample (`27522_1`), and a normal-BM control (`BM4`): `TNFRSF17`, `GPRC5D`,
`SLAMF7`, `FCRL5`, `SDC1`, `CD38`, `ITGB7`, `NCSTN`, `IGKC` are present in **both**
reference builds. The 33538 build additionally carries `GPRC5D-AS1` (an antisense
transcript, absent from 33694) — harmless, dropped by the intersection, and must
**not** be substituted for `GPRC5D` anywhere. The hard assertions in
`gene_space.py` are still required as a regression guard, but are expected to pass.

### GEO series metadata — added and parsed 2026-08-24

`raw/GSE22306{0,1}_family.soft.gz` (~8 KB each) carry per-sample facts that are **not
derivable from filenames** and were absent from this document until 2026-08-24. Parsed
into committed tables at `resources/sample_metadata/` (`raw/` is gitignored, same
pattern as `resources/gene_space/`); regenerate with
`io.rebuild_sample_metadata_from_soft`, which asserts sample counts and cohort
resolution so a revised deposit fails loudly.

**Cohort and prep, from `!Sample_extract_protocol_ch1`:**

| cohort | n | chemistry | dead-cell removal | reference build |
|---|---|---|---|---|
| WashU 1 | 23 | 10x 3′ **v2** | **no** | 20× 33694, 2× 33538, 1 truncated |
| WashU 2 | 13 | 10x 3′ v3.2 | yes | 13× 33538 |
| MMRF | 18 | 10x 3′ v3.3 | yes | 18× 33538 |
| Donors | 8 | 10x 3′ v3.2 | yes | 4× 33694 (`ND_*`), 4× 33538 (`BM*`) |

**This is a confounder for the headline metric, and it is measured rather than
assumed.** Sample-level medians of genes detected per pre-QC cell:

    MMRF   v3.3   1916 genes/cell
    WU2    v3.2   1210
    Donor  v3.2   1103
    WU1    v2     1023

    v2 vs all-v3: 1023 vs 1408 = 1.38x, Mann-Whitney p = 6.5e-05,
    but the sample distributions OVERLAP (v2 max 1602 > v3 min 793).

**Do not quote a "2-3x v2-vs-v3 chemistry effect" — this cohort does not show one.**
The axis that separates is **cohort** (MMRF ≈ 1.9× the others), of which chemistry
version is one component alongside site and protocol. It still must be modelled:
`frac_double_negative` is a fraction of zeros on a low-abundance transcript, so a 1.9×
depth spread that tracks cohort will move it and read as biology. Carry `cohort` and
`chemistry` as covariates in stage 08's depth regression and stage 10's null.

**`n_genes_ref` is NOT a proxy for this** — the build split cuts across cohorts (two
WU1 samples on 33538, the four `ND_*` donors on 33694). Stage 05's Harmony covariate
needs `cohort`/`chemistry` in addition to, not instead of, `n_genes_ref`.

**A free control:** the 8 donors span both references and are normal marrow, so
stage 07's negative control doubles as a build/chemistry control in a population with
no clone to confound it.

**Raw data exists, under controlled access — "no raw data" was too strong.** The
series design says raw data goes to dbGaP for patient privacy, with `ND_*` already
available at **`phs000159`** and MMRF bulk at **`phs000748`** (BioProjects
`PRJNA924769` / `PRJNA924778`). The practical conclusion is unchanged — no unfiltered
matrices without a DAC application, so SoupX/DecontX remains unavailable — but
"doesn't exist" and "exists behind a data-access committee" are different claims.

**The deposit's own processing metadata is unreliable.** `!Sample_data_processing`
claims Cell Ranger **v3.0.0 for all 62 samples**, which the files contradict for 24 of
them (Ensembl 84 / CR 1.2.0). The files win — that reconstruction is checksum- and
position-verified. Treat the rest of the SOFT text as evidence to check, not fact;
every claim taken from it above was verified against the data.

**Matched bulk RNA-seq (GSE223061) — already downloaded, previously unused.**
`raw/unpacked_bulk/` holds **29 usable bulk samples** (inventory corrected
2026-08-21 by direct count in `notebooks/01_download_data.py`, which now asserts
these numbers): **18** MMRF samples as `<GSM>_<sample>_tpm.tsv.gz` (gene × TPM
tables, GSM6939103-120) of which 2 are empty stubs, plus **13** WashU samples as
`<GSM>_<sample>.tar.gz` (GSM6939090-102, all 4.5-5.4 MB). The earlier "18 + 12 = 30
usable" was wrong twice over — the WashU count was 12 and the total did not subtract
the stubs. Correct arithmetic: (18 - 2) + 13 = 29.
**Overlap COMPUTED 2026-08-24 from the GEO titles, no longer inherited or
S1-gated: 26 bulk samples have an exact scRNA match** (not "~28"). 31 bulk GSMs, minus
2 empty stubs, minus 3 with no scRNA counterpart.

**The two bulk cohorts are not the same assay, and stage 09 must not pool them:**

| bulk cohort | n | prep | matching pseudobulk |
|---|---|---|---|
| MMRF | 18 (16 usable) | **CD138+ sorted** | malignant-cell pseudobulk |
| WashU 1 | 13 | **unsorted BMMC** | **whole-sample** pseudobulk |

Correlating malignant-PC pseudobulk against *unsorted* bulk measures tumour burden,
not antigen abundance — the dilution by non-plasma cells is proportional to burden,
which varies per patient and is itself correlated with the metric. That would have
silently corrupted 10 of the 26 matched comparisons. Split by cohort and use the
matching comparator.

Other gotchas, confirmed by direct inspection:
- **Two files are empty 114-byte gzip stubs and must be excluded**:
  `GSM6939104_MMRF_1505_tpm.tsv.gz`, `GSM6939120_MMRF_2259_tpm.tsv.gz`. Both are real
  GEO sample entries — a 114-byte read is a failed deposit, not "zero expression".
- **The three bulk/sc ID mismatches, partly resolved.** `59114` was overstated: bulk
  carries `59114_1` *and* `59114_2`, and `59114_1` matches scRNA exactly — only `_2`
  is orphaned. `47499` and `98433` have no scRNA counterpart; note both are WashU
  **cohort 1** bulk while `MMY98423` is cohort **2**, so the assumed
  `98433` ↔ `MMY98423` pairing is doubtful rather than a simple typo. Do not guess
  these; they are 3 of 29 and stage 09 works without them.

**Patient mapping — SOLVED 2026-08-24 by Supplementary Table S1.** S1 landed in
`raw/` and closes this. The naive rule (strip a trailing `_<digits>` only when the
stem is purely numeric, e.g. `27522_1` -> `27522`) gives **43 patients from the 54
myeloma samples** against the paper's 41 / 53. Both gaps close exactly, and the two
corrections are independent of each other:

- **`25183` is deposited but appears in NO supplementary table** — not the clinical
  summary, not the disease-stage sheet. It is what the 53-vs-54 gap is made of.
  `54 - 1 = 53`. It is **not dropped**: the data is real and stage 07 can use it. It
  carries `in_paper_cohort == False` and `clinical_source == "none"` so a per-patient
  aggregate excludes it deliberately rather than by accident.
- **`83942` (WashU 1) and `MMY83942` (WashU 2) are one patient.** S1 lists them as
  two, but with identical age (63), gender (Male), race (White), ISS stage (3) and
  treatment (Unknown). One patient sampled under both WashU protocols.
  `43 - 1 - 1 = 41`.

Implemented in `io.rebuild_clinical_metadata_from_s1` -> `resources/sample_metadata/`
(`patients_clinical.tsv`, 43 rows; `sample_disease_stage.tsv`, 22). `io.s1_patient_id`
is the mapping; `load_manifest` now emits `patient_id_source == "S1"` and keeps the
naive answer as `patient_id_naive` for comparison. `_assert_s1_reproduces_the_paper`
checks all three counts off the committed GEO table, so a revised S1 fails loudly
instead of quietly moving the denominator of `frac_double_negative`.

**What else S1 carries.** Per patient: age, sex, race, ISS stage, treatment regimen,
time-to-progression-or-death. Per sample: disease stage — but **WashU cohort 1 only**.
MMRF and WashU cohort 2 get `disease_stage = NA` and one is *not* imputed;
`newly_diagnosed` is a guess for those cohorts, not a datum. Two S1 patients
(`47499`, `98433`) are bulk-only and have no scRNA sample — they stay in the patient
table for stage 09.

**The `_N` suffixes are serial disease-course timepoints. SETTLED, not inferred.**
S1 sheet 2 reads:

    27522_1  Primary     27522_2  Remission-1   27522_3  Relapse-1
    27522_4  Relapse-2   27522_5  Remission-2   27522_6  Relapse-3
    47491_1  SMM         47491_2  Primary       58408_1  SMM  -> 58408_2 Primary

This confirms outright what the 2026-08-24 bulk/scRNA suffix argument had inferred, and
**the longitudinal arm (the S1-gated additions in the main project document) is now real
rather than speculative**. It also explains the lone non-`_1` samples: `37692_2` and
`57075_3` are later timepoints whose earlier draws were not deposited.

**An earlier figure of 47 patients / 57 disease samples is wrong and is superseded:
it counted the four `ND_*` samples as disease. SETTLED 2026-08-24 against the GEO
metadata**, which is now in the repo (see "GEO
series metadata" below). `!Sample_source_name_ch1` reads **`Donor BMMC, aspirate,
scRNAseq`** for all four `ND_*` and all four `BM*`, and those eight are exactly the
samples carrying no `diagnosis` characteristic; the other **54 read
`diagnosis: Multiple myeloma (MM)`**. `ND` = normal donor, and the suffixes are
collection dates.

| | samples | naive patients |
|---|---|---|
| `ND_*` counted as disease (the old figure) | 57 | 47 |
| **`ND_*` as donors, `56203_1` repaired and retained** | **54** | **43** |

The series summary states "53 bone marrow (BM) aspirates from 41 MM patients"
verbatim, matching the paper — and matching what the S1 mapping above now produces.
The denominator is settled; there is no remaining uncertainty in it.

**The other five supplementary tables are committed too** (~500 KB, `raw/`, see
`.gitignore` for which is which). Two are directly useful downstream and neither is
yet consumed:
- **Table S3** (file `s3`) — a 38x38 Pearson correlation matrix over the paper's
  candidate target genes, plus per-cohort co-expressed and mutually-exclusive gene
  pairs. It includes `GPRC5D` x `TNFRSF17`, and it is a caution as much as a
  comparator: the pooled r is **0.064**, but per cohort it is **MMRF +0.62,
  WU2 +0.54, WU1 -0.09** — the pair appears in the paper's co-expressed *and*
  mutually-exclusive lists depending on cohort, and the sign tracks the cohort depth
  ordering exactly. Whatever stage 08's co-escape enrichment finds must be read
  against that. (Note these are sample-level correlations, not the per-cell
  co-negativity stage 08 computes — related question, different unit.)
- **Table S5** (file `s6`) — the paper's own bulk-vs-scRNA Pearson r per gene for
  MMRF. A direct external comparator for stage 09.

Not useful here: Table S2 (target novelty tiers, mild framing value), Table S4
(NetMHC peptide-HLA binding), Table S6 (recurrent mutations, over different MMRF
samples). **File `s5` is strict-OOXML** — `openpyxl` reports zero sheets for it and it
is not corrupt; and files `s5`/`s6` are numbered off-by-one against their in-file
titles.

**Ambient RNA correction (SoupX/DecontX) is not possible for this dataset.** Both
require the *unfiltered* Cell Ranger matrix (including empty droplets) to estimate
the background contamination profile; GEO only hosts the *filtered* per-sample
matrices. This is a real, stated limitation, not an oversight — see "QC
methodology" below for the mitigation actually used instead.

**Datasets explicitly evaluated and rejected as data sources:**
- **He et al. 2022** (`ctm2.757`, paired NDMM/RRMM scVDJ-seq) — confirmed by
  directly reading the paper (PMC8926895): no data availability statement, no
  accession of any kind. Not usable as a data source under any framing. A
  fabricated-looking accession (`GSE124310`) was floated for this paper at one
  point and independently verified to actually belong to an unrelated
  glioblastoma study — do not resurrect that number.
- **GSE118900** (MGUS/SMM/NDMM/RRMM staging cohort) — real and open, but rejected
  as the Phase 2 validation dataset: zero healthy controls, only 597 cells across
  15 patients (as low as 7-24 per patient), too statistically thin for a per-cell
  classification metric. `GSE117156` (see the main project document's Phase 2
  section) is superior on every axis that matters here.
