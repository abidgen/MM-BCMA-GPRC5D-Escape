# RESUME HERE — MM Dual-Antigen pipeline (Python rebuild), session state

**Last updated:** 2026-08-21
**Branch:** `review3-stage08-corrections` — one unmerged commit (doc corrections).
`main` holds everything else. Working tree clean.

**To resume:** `git checkout main && git merge --ff-only review3-stage08-corrections`
(or review the diff first), then start at "Immediate next actions".

Read `CLAUDE.md` first for the settled decisions and data ground truth. This file
covers only *where execution stands* and *what to do next*.

---

## TL;DR

Clean start. Previously built substantially in R (data acquisition solved and
verified on all 62 samples, QC/doublet-removal run on the full 61-sample cohort,
integration not yet run) before switching to Python. **None of the R code is being
ported.** All dataset knowledge carries forward via `CLAUDE.md`.

**Stages 01-03 are complete: written, executed, asserted, committed.** All three
notebooks exist and run green (62/62 `triplet-ok`, manifest byte-identical via CLI and
notebook). The four conda envs are built and their kernels registered. `src/mm_escape/`
has `config.py` + `gene_space.py` scaffolded and tested against the real files.

**Stage 04 onward does not exist yet.** The next artifact to write is
`src/mm_escape/io.py`.

Three real defects have surfaced so far, all fixed: cross-reference HGNC symbol drift
(solved properly via Ensembl-ID reconstruction), a wrong bulk RNA-seq inventory in
`CLAUDE.md`, and a circular dropout-correction formula in the stage 08 plan. Details
below.

---

## 2026-08-20 — scope expansion + working-tree cleanup

Two things happened this session. A design review of the analysis plan ran and
**all proposed additions were adopted**, with all five `.md` files updated to match.
Then the R build was removed from the working tree. No Python was written. Full
reasoning lives in `CLAUDE.md` — the summary:

**The motivating problem.** `frac_double_negative` is a *fraction of zeros*, the
noisiest quantity scRNA-seq produces, and the plan bounded only one of its two error
directions. Ambient RNA (documented) deflates it; **dropout (undocumented) inflates
it, and is the larger effect here** because `GPRC5D` is a low-abundance transcript
and this cohort's median cell has ~2,044 detected genes.

**What changed:**
- **Stages 09 and 10 added, and the whole sequence renumbered** so number order
  equals execution order. 09 = escape robustness (matched bulk RNA-seq validation,
  normal-PC antigen baseline, permutation null). 10 = subclone + phenotype (is the
  DN population a real subclone or scattered noise — **reframed 2026-08-21, see the
  second design review below**; pseudobulk DE; pre-registered γ-secretase hypothesis). Cell-cell communication moved 09 → **11**; decision
  packet moved 10 → **12**.
- **Stage 08 gained its defense layer** — threshold sensitivity band, depth/dropout
  falsification test, expression-matched false-negative floor, bootstrap CIs with a
  declared minimum-cell rule, and a multi-antigen combinatorial coverage matrix.
- **Stage 07 hardened** — `infercnvpy` promoted from optional to required with a
  reported agreement rate; ratio-based (not presence-based) light-chain calls;
  normal-BM samples used as a negative control.
- **Stage 11 statistics corrected** — the tertile split was pseudoreplication;
  replaced with per-patient scores, continuous escape fraction, T/NK abundance
  controlled.
- **Stage 05 integration scoped** — Harmony for the immune compartment only;
  malignant subclustering per-patient and un-integrated.
- **Stage 06** now emits per-patient composition as a first-class output; `scCODA`
  for compositional comparisons.
- **Number order is execution order**, no exceptions: 04 → 05 → 06 → 07 → 08 → 09
  → 10 → 11 → 12.

**Two defaults applied this session (reversible, flag if you disagree):**
1. **Numbering** — the sequence was renumbered so stage numbers run in true
   execution order, rather than appending the new stages at 11/12 and leaving the
   decision packet numbered 10 but running last.
2. **Supp. Table S1** — build everything S1-independent first on the provisional
   mapping, labelling every S1-dependent number provisional *in the output itself*;
   attempt retrieval when stage 08 aggregation is reached. S1 no longer blocks the
   pipeline, only four specific items (see `CLAUDE.md`'s S1 policy).

## 2026-08-20 — architecture: the whole analysis runs in notebooks (01-12)

**Decision.** Every stage is a notebook, 01 through 12 — nothing left as a bare CLI
script you can't open and step through. `scripts/01-03` are **kept as a headless CLI
fallback** (fresh clone, remote box, CI) and are *wrapped* by notebooks 01-03, never
reimplemented; byte-identical output is verified both directions.

**`src/mm_escape/` is retained.** The earlier rule "notebooks are thin orchestration" is
deliberately relaxed: notebooks now carry the analysis — narrative, plots, intermediate
inspection. `src/mm_escape/` holds what earns being a library: **reusable** (used by
more than one stage), **testable** (worth asserting on independently), or **fiddly**
(the `read_mtx` loader, symbol harmonization, the noise-floor derivation). The test is
reuse and testability, **not line count** — a one-off plot belongs in the notebook; a
threshold used by three stages does not.

The two things this protects, unchanged from the original rationale: Codex reviews `.py`
diffs rather than notebook JSON, and logic with a single home cannot drift between
copies.

**To write:** `notebooks/01_download_data.ipynb` and `02_check_files.ipynb` (wrapping
the two bash scripts — likely thin, since those scripts are pure file handling).
`03_build_manifest` already exists and is the reference pattern for how 01/02 should
wrap their scripts.

---

## 2026-08-20 — stages 01-03 executed; symbol-drift defect found

**First actual execution this rebuild.** `scripts/01_download_data.sh` and
`02_check_files.sh` ran as the predicted no-op confirmation: 62/62 archives already
present, 0 re-downloaded, 0 re-extracted, all **62 `triplet-ok`**, no INCOMPLETE.
(`01`'s summary prints `Matrix files: 0` — its counter globs `*matrix*` while this
archive uses `counts.mtx`. Cosmetic, not missing data.)

**`notebooks/03_build_manifest.py` added** (percent-format; `.ipynb` is gitignored and
generated from it). At the time `CLAUDE.md` still ruled out notebook-ifying stages
01-03, so this began as an *additional* interactive view — **since superseded by the
full-notebook decision above**, which makes it the reference pattern for how notebooks
01-03 wrap their scripts. Either way it imports `build_manifest()` from the script and
writes the manifest on the script's exact schema. **Byte-identical
output verified in both directions.** Two traps documented in it, both real and both
hit during development:
- Under a Jupyter kernel `sys.argv[1]` is `-f`, so the script's module-level `RAW_DIR`
  becomes `Path('-f')`. Pass paths explicitly; never use `mf.RAW_DIR`.
- The manifest must hold repo-root-relative paths (it is committed and read elsewhere),
  so the notebook normalizes the absolute paths it has to pass in.

### The finding: cross-reference HGNC symbol drift

The notebook's required-gene assertions **failed on `NSD2`** — and the cause was not a
missing gene. The two Cell Ranger references use different HGNC symbol vintages, so
genes present in *both* builds under different names are dropped by a naive symbol
intersection:

| 33538 (newer) | 33694 (older) | consequence |
|---|---|---|
| `NSD2`    | `WHSC1`   | **t(4;14) uncallable** — highest-risk MM translocation, and stage 10's `MS` TC class |
| `TENT5C`  | `FAM46C`  | loses a recurrently-deleted MM tumour suppressor (1p12) |
| `NSD3`    | `WHSC1L1` | a *different* gene from NSD2 — do not conflate |
| `ATP5F1A` | `ATP5A1`  | OXPHOS program member |

**Fix: canonicalize symbols before intersecting** (22,164 → 22,168 genes). Belongs in
`src/mm_escape/gene_space.py`; the notebook carries a targeted version covering the
genes this project needs. A fuller HGNC reconciliation is worth doing there, since the
same drift certainly affects genes outside the required set.

**Rule going forward: a missing required gene means "check for a legacy symbol" before
concluding "biologically absent".** Two prior builds of this project missed this; the
assertions caught it in minutes, which is the argument for keeping them loud.

Also confirmed: the ~11.4k/11.5k symbols unique to each build are dominated by
annotation-version noise (`AC000032.1` vs `AC000032.2`), so **22,164 understates the
recoverable gene space** — fine for lncRNA/clone entries, not fine for named genes.

**Stopgap kernel.** `envs/` is still unscaffolded, so the notebook was validated
against the leftover `mm-dual-antigen` conda env (pandas 2.2.2, py3.11) and a kernel of
that name was registered. **Replace with `mm-qc`/`mm-core`/`mm-annotation`/
`mm-communication` once `envs/*.yml` are built** — the notebook's kernelspec will need
repointing.

---

## 2026-08-20 — stage 06 annotation methodology settled

`CLAUDE.md`'s stage 06 said "`celltypist` **and/or** marker-panel scoring." That
`and/or` was an unmade decision sitting mid-pipeline. It is now settled and documented
across all four docs. No code written — `src/` and `envs/` are still unscaffolded.

- **Three methods, compared**: manual marker panel, `celltypist`
  (`Immune_All_Low`/`High`, `majority_voting`, over-clustered on the stage-05 Leiden
  key), and `SingleR` against `celldex`'s `NovershternHematopoieticData`. SingleR was
  chosen over scArches specifically because `Immune_All_*` is immune-only and its
  predictable blind spot is erythroid/HSPC — a sorted-hematopoietic reference targets
  exactly that gap.
- **Decision is per class, not per method** (hybrid), with F1 thresholds declared
  before looking: PlasmaCell ≥ 0.95, T/NK/Myeloid ≥ 0.90, rest ≥ 0.85. Written to
  `results/06_annotation/annotation_decision.md`.
- **`obs["cell_type"]` is the only load-bearing output**; `cell_type_fine`,
  `annotation_source`, `annotation_conf` are provenance/convenience. Stages 07-12 read
  `cell_type` and nothing else, so the annotation choice stays reversible.
- **Identity vs. state separated**: cell-cycle / IFN / antigen-presentation / UPR /
  hypoxia are continuous `obs` floats, never categorical, never merged into
  `cell_type`.
- **New env `mm-annotation`** (`envs/env-annotation.yml`) — SingleR is R, and R stays
  isolated the way `env-qc` isolates scDblFinder. Stage 06's env changed from
  `mm-core` to `mm-annotation` in every doc. Four envs now, not three.
- **Stage 10 gained** malignant-cell program scoring (adds **MYC** and **OXPHOS** to
  the stage-06 programs) and **TC molecular subgroup** assignment per patient
  (`CCND1`/`CCND2`/`CCND3`/`NSD2`/`FGFR3`/`MAF`/`MAFB` + `CKS1B` for 1q21 gain, which
  also cross-checks `infercnvpy`).
- **UAMS 7-group explicitly rejected** — needs bulk-array signature sourcing and
  splits n≈41 into unpowerable bins; TC gives most of the interpretive value for far
  less. TC is used **descriptively**, not as a statistical stratifier at this n.
- **No custom `celltypist` model** for malignant states — a linear classifier bins
  continuous tumor substructure and hides intermediates.

**On external LLM suggestions (Gemini/ChatGPT).** A first round was answered against a
**melanoma** dataset and none of its gene-level content applies here (`MITF`, `MLANA`,
`AXL`, `NGFR`, neural-crest states, Tirosh/Jerby-Arnon/GSE115978 — melanocyte biology
with no plasma-cell counterpart). Its *structural* advice did transfer and was adopted:
continuous non-exclusive program scores, and no custom classifier. A second round,
correctly about myeloma, **converged on the architecture this project already had**
(automated labels for immune → light-chain + CNV for malignancy → program scoring for
substructure = stages 06 → 07 → 10). Genuinely additive from it: the MYC and OXPHOS
programs, and `CKS1B`/1q21. Treat these tools' disease context as unverified until
checked.

---

**New data facts confirmed this session (read-only inspection):**
- All target genes (`TNFRSF17`, `GPRC5D`, `SLAMF7`, `FCRL5`, `SDC1`, `CD38`, `ITGB7`,
  `NCSTN`, `IGKC`) present in **both** reference builds — the intersection costs no
  markers. 33538-gene build also has `GPRC5D-AS1`; do not substitute it for `GPRC5D`.
- `raw/unpacked_bulk/` holds **30 usable GSE223061 bulk samples** (18 MMRF TPM tables
  + 12 WashU archives), ~28 overlapping the sc cohort. Previously unused by the plan;
  now the stage 09 validation set.
- **Two bulk files are empty 114-byte stubs** — `MMRF_1505`, `MMRF_2259`. Exclude;
  do not read as zero expression.
- **Three bulk/sc ID mismatches** for S1 to settle: bulk `47499` vs sc `47491_1/2`;
  bulk `98433` vs sc `MMY98423`; bulk `59114_2` vs sc `59114_1`/`59114_4`. The suffix
  misalignment across assays is evidence the `_N` suffixes are **not** simple
  timepoint indices — don't assume they are.
- R-build QC numbers carried forward as cohort context: 61 samples, **181,336 post-QC
  cells**, median 2,555 cells/sample (min 480, max 7,937), median `nFeature` ≈ 2,044.
  Source `results/qc_summary_per_sample.csv` now lives only in the snapshot —
  `git show r-build-snapshot:results/qc_summary_per_sample.csv` to re-read it.

**Working-tree cleanup (done this session).** `env_creation/`, `scripts/04{a,b,c}_*.R`,
`scripts/lib/`, the two `scripts/__*` pre-fix backups and the entire `results/`
directory were removed. Everything is recoverable from the **`r-build-snapshot`** tag
(commit `8ee8624`), which captured all 642 MB including the 638 MB of per-sample QC
`.rds` checkpoints. `raw/` was not touched.

---

## Immediate next actions, in order

1. ~~**Re-run `scripts/01-03`**~~ — **DONE 2026-08-20.** 62/62 `triplet-ok`, manifest
   verified byte-identical via both the CLI and `notebooks/03_build_manifest.ipynb`.
2. ~~**Scaffold the repo**~~ — **DONE 2026-08-21.** Five env ymls written; four built.
   `src/mm_escape/` has `__init__.py`, `config.py` (partial — gene-space constants
   only) and `gene_space.py`. Notebooks 01 and 02 written, executed, committed.
3. ~~**Build `env-qc`**, register kernels~~ — **DONE 2026-08-21.** All four built and
   registered; every key package import-verified, both R bridges confirmed working.
4. **>> START HERE << Write `src/mm_escape/io.py`**: the loader replacing `scanpy.read_10x_mtx()`.
   Handles `counts.mtx` (not `matrix.mtx`), single-column `genes.tsv`, the extra
   nesting level per sample. Validate against 2-3 real samples before scaling to
   all 61 — most likely place for a silent format-assumption bug to hide.
5. **Write `src/mm_escape/qc.py`** and `notebooks/04_qc.ipynb` -> `results/04_qc/`:
   QC metrics, MAD-based outlier calling (start from 5 MAD, document the actual
   resulting thresholds for this cohort), `scDblFinder` via `rpy2`. Exclude
   `56203_1`. Checkpoint per-sample.
6. Continue through stages 05-12 per `CLAUDE.md`'s pipeline section — notebook
   number N always writes to `results/N_*/`, nothing else. Run them in numeric
   order; number order is execution order.
7. **Apply the review corrections** as each stage is written — they are documented in
   place in `CLAUDE.md`, not collected in one section. Summaries in "Second design
   review" and "Third design review" below. Stage 08 carries the most of them; read
   both summaries before writing it.

### Presentable stopping points

The stage ordering was chosen so the project is showable partway through rather than
all-or-nothing:
- **After 04-08:** escape fractions plus co-escape enrichment per patient.
- **After 09:** results that survive hostile questioning.
- **After 10:** the actual scientific finding (is the escape population structured).
- **After 11-12:** the full decision packet.
- **Then** S1 retrieval (un-flag provisional numbers, add longitudinal + cytogenetics),
  and only then Phase 2.

---

## Environments (built 2026-08-21)

`envs/env-{qc,core,annotation,communication}.yml` are written and built as conda envs
`mm-qc` / `mm-core` / `mm-annotation` / `mm-communication`. `envs/env-composition.yml`
(scCODA + TensorFlow) is **written but not built** — on demand only, per CLAUDE.md.

Register kernels once per env: `python -m ipykernel install --user --name mm-<env>`.

Kernels are registered for all four; each kernelspec was confirmed to point at its own
interpreter. **All four verified to import their key packages, and both R bridges work**
(mm-qc: R 4.3.3 + scDblFinder 1.16.0; mm-annotation: SingleR 2.4.0 + celldex 1.12.0 with
both Novershtern and HumanPrimaryCellAtlas references present). mm-core has infercnvpy
0.6.1, pydeseq2 0.5.4, decoupler 2.2.0, harmonypy 2.0.0; mm-communication has liana
1.8.1 with both `mt.cellchat` and `mt.rank_aggregate`.

Five deviations from CLAUDE.md's original env specs, all found by actually building and
importing rather than trusting the yml:
- `infercnvpy` is **pip-installed inside `env-core`** — not packaged for conda anywhere.
- `celldex` is really `bioconductor-celldex`.
- `decoupler` is **pip 2.2.0**, because bioconda's is stuck at 1.5.0, pins `numpy<2`, and
  fails to import under numba 0.67. **2.x is an API rewrite (`dc.mt.*`/`dc.op.*`) — do
  not write stage 10 against 1.x tutorials.**
- `liana` needed `pydantic` + `mudata<0.4` pinned; bioconda's recipe under-declares and
  `import liana` failed outright without them.
- `jupytext` added to all four (every notebook is paired).

**Do not `pip install` into these envs casually.** Trying `decoupler==1.8.0` during setup
silently downgraded numpy 2.5.2 → 1.26.4 and numba, breaking scanpy/scipy/pydeseq2/zarr
at once; the fix was to delete and recreate mm-core from the yml. Both legitimate pip
entries live in the yml so a rebuild reproduces them.

All pinned versions (`scanpy=1.11`, `r-base=4.3.3`, `rpy2=3.5.11`,
`bioconductor-scdblfinder=1.16.0`) exist and were confirmed pre-build.

---

## What's reused unchanged from the R build

`scripts/01_download_data.sh`, `02_check_files.sh`, `03_build_manifest.py` and the
`raw/` directory's contents — already solved and verified against the real
62-sample archive.

## What's being rebuilt from scratch

Everything from stage 04 onward: the loader, QC, integration, annotation,
malignant calling, antigen scoring, escape-fraction computation, and the
cell-cell-communication step (LIANA+ instead of CellChat directly, though LIANA+
natively reimplements CellChat's own algorithm).

---

## Known blockers / decisions already made (do not re-litigate — see CLAUDE.md for full reasoning)

- `56203_1` excluded (bad reference build, missing BCMA).
- Gene sets intersected across samples, never unioned.
- Ambient RNA correction (SoupX/DecontX) is impossible for this dataset — no
  unfiltered matrices exist. Mitigated via an empirical noise-floor threshold on
  antigen positivity instead of a naive `>0` call.
- **Dropout is the opposite-signed, larger bias and must be bounded too** — the
  headline metric is reported as an interval with a threshold sensitivity band.
  Ranking stability is the robustness *diagnostic*; **risk tiers** (robust-high /
  uncertain / robust-low) are the stage-12 deliverable, not an ordinal ranking.
- Patient-ID mapping is unresolved (47 vs. 41 patients) — fixed against
  Supplementary Table S1 before stage 08's per-patient aggregation is *final*, but
  per the 2026-08-20 policy this **no longer blocks execution**: proceed
  provisionally with provisional labelling in the output.
- Supplementary Table S1 itself is still not in the repo.
- `infercnvpy` is required (not optional) and normal-BM samples are controls
  (not filler).
- Stage 11 uses patient — not cell — as the unit of replication, and is **explicitly
  exploratory** (ninth in the scientific hierarchy, not co-equal with the antigen work).

### Second design review (2026-08-21) — five corrections + one added analysis

All adopted; documented in place in `CLAUDE.md`. Do not re-litigate.

1. **Transcriptional clustering ≠ clonality.** `clonality-of-escape` is retired for a
   three-level **DN-coherence hierarchy**; the word *subclone* requires CNV support
   (level 3), and level 3 reports **supported / not evaluable** — never "no subclone",
   because within-clone CNV resolution is often underpowered at ~2,044 genes/cell.
2. **Bulk RNA-seq validates antigen *abundance*, not the DN fraction** — bulk averages
   cells and destroys the joint distribution the metric is made of.
3. **A dropout-adjusted probabilistic DN estimate runs alongside the binary call.**
   Binary stays primary; **imputation/denoising stays forbidden** for positivity.
4. **The "label-permutation null" moved 09 → 08** — it was never a no-signal null; it
   tests independence between the two antigen-negativities.
5. **Normal-BM gives marrow expression context, not a safety axis** — GPRC5D's decisive
   off-tumor site (keratinized tissue) is invisible to a marrow dataset.

**Added: BCMA/GPRC5D co-negativity enrichment (stage 08)** — per-patient 2×2, observed
DN vs. `P(BCMA⁻)×P(GPRC5D⁻)`. Distinguishes two independent partial failures from a
coordinated antigen-low phenotype, which is what determines whether a second binder
helps at all.

**Two caveats added on top of the review, both load-bearing:**
- The independence null **must be depth-stratified**. Shallow cells read zero for both
  genes, so an unconditioned permutation manufactures co-escape enrichment from library
  size — biased toward the project's own hypothesis. Report the unconditioned ratio
  next to the conditioned one; the gap is the artifact.
- `infercnvpy` may not resolve **sub**clonal CNV within one patient's tumor. Underpowered
  ≠ negative. Same guard applies to stage 10's level-1 enrichment test, which is also
  depth-stratified.

**Other adjustments:** stage 07 emits malignant-confidence tiers (`high`/`probable`/
`uncertain`) with the stage-08 result re-run on `high` only as a sensitivity analysis;
stage 06's F1 numbers are **concordance**, not validation accuracy (marker-coverage test
is the biological evidence and can veto a class); stage 08's minimum-cell rule is
re-derived from needed resolution (expect 100-200, not 50) with hierarchical
patient→sample→cell bootstrapping; TC classes are labelled **TC-like expression
subtypes**, never translocation calls.
- Malignant subclustering is per-patient and un-integrated; Harmony is for the
  immune compartment only.
- GSE117156 is the confirmed Phase 2 validation dataset (own `phase2_NN_*`
  numbered notebooks, never merged with the primary cohort), sequenced strictly
  after Phase 1 (stages 01-12) completes; GSE118900 and He et al. 2022 were both
  evaluated and rejected as data sources.
- Notebooks and `results/` subdirectories are numbered `04`-`12`, one-to-one,
  continuing on from `scripts/01-03`. `src/mm_escape/` modules are named by
  function, not numbered.

## 2026-08-21 — envs built, gene space solved, notebooks 01-02, third review

Four things happened. Everything is committed.

### 1. Ensembl-ID reconstruction — the symbol-drift problem is properly solved

The 2026-08-20 fix was a 4-gene alias map. That was a patch. The real fix: the deposit
has no ID column (`genes.tsv` is symbol-only, **zero `ENSG` strings across all 62
samples** — do not go looking), but there are only **three distinct gene files in the
whole cohort**, byte-identical within group, each a positional dump of a public
reference. So the IDs are reconstructible and this was verified end to end:

    Ensembl 93 GTF -> 10x mkgtf biotype filter, GTF order -> Seurat's gsub("_","-")
      -> R make.unique -> assert == deposited column, position for position

**0 mismatches / 33538 rows and 0 mismatches / 33694 rows.** Self-certifying: a wrong
release or filter changes the count or order and fails loudly.

| join key | genes retained |
|---|---|
| raw symbols | 22,164 |
| symbols + 4-gene alias map | 22,168 |
| **Ensembl IDs** | **32,991**  (+10,827, ~49% more) |

11,140 intersected IDs carry a different symbol in each build. It also proved the
mis-pairing risk is real, not theoretical: `TBCE` is `ENSG00000285053` in the 33538
build and `ENSG00000116957` in the 33694 build — a symbol join silently merges two
different annotation entries.

Committed: `resources/gene_space/` (~1 MB gzipped, so the 41-44 MB GTFs never need
re-downloading) and `src/mm_escape/gene_space.py`. Pipeline is four calls in order:
`attach_ensembl_ids` -> `intersect_gene_space` -> `to_canonical_symbols` ->
`assert_required_genes`. Four failure paths confirmed to raise.

**`var_names_make_unique()` is banned on these objects** — it assigns bare-vs-suffixed
names by row position, exactly the mangling being undone. The 9 symbols still colliding
after the ID intersection get `SYMBOL__ENSG...`.

### 2. Environments built and verified — see the Environments section above

Two of the four were unusable as specified (`import liana` failed outright;
`decoupler` failed under numba 0.67). Both fixed in the ymls.

### 3. Notebooks 01 and 02 written; bulk inventory corrected

Both execute clean against the real data. Scripts 01/02 are bash, so these notebooks
`subprocess` out to them rather than importing a function as notebook 03 does.

Notebook 01's assertion immediately caught an error in `CLAUDE.md`'s ground truth:
recorded as "18 MMRF + 12 WashU = 30 usable" bulk samples; actually **13** WashU
`.tar.gz` and 18 MMRF TPM (2 of which are the known 114-byte stubs), so **(18-2)+13 =
29 usable**. The old figure was wrong twice — the WashU count, and the total never
subtracted the stubs. Both corrected.

The inherited "~28 samples overlap with the scRNA cohort" is now flagged **unverified**
— it depends on the three bulk/sc ID mismatches and the unresolved S1 mapping. Stage 09
recomputes it rather than quoting it.

Notebook 02 also asserts there are exactly **3 distinct `genes.tsv` checksums**, so a
new reference build fails loudly instead of merging silently.

Notebook 03's kernelspec was an R-build leftover (`mm-dual-antigen`) -> `mm-qc`.

### 4. Third design review — one real methodological error, fixed

Details in `CLAUDE.md` stage 08. The headline: the proposed "dropout-adjusted expected
DN", `sum_i P_i(BCMA-) * P_i(GPRC5D-)`, was **circular** — it assumes exactly the
independence the co-escape test exists to interrogate. Consequences:

- That quantity and the "expected DN under depth-conditioned independence" are the
  **same number**, specified twice. Reported once now, as a technical baseline.
- **No dropout-corrected DN point estimate is produced, and none is claimed.** Dropout
  is *bounded* (sensitivity band, false-negative floor, depth regression, downsampling),
  not corrected. The honest correction is a latent-class/EM model over the four true
  states — specified in `CLAUDE.md`, deliberately deferred, not on the critical path.
- **Bootstrap level was wrong**: a CI *for patient A* conditions on patient A, so
  patient is fixed. Per-patient CI = **sample -> cell**; cohort-level =
  **patient -> sample -> cell**. Single-sample patients get optimistic CIs; say so.
- **Co-escape was over-read**: it measures *eroded complementarity*, not futility.
  Adding GPRC5D to BCMA still moves 30% -> 15% under strong co-loss. So
  **incremental coverage gain** `P(A-) - P(A- and B-)` is now reported beside it and
  sits above co-escape in the hierarchy — it is what a target decision turns on.
- **No composite risk score anywhere**, including the coverage matrix.

---

## Status

Stages 01-03 complete and green. Envs built, kernels registered, gene space solved and
committed. Architecture stable — the third review round hit wording and interpretation
rather than method, which is the signal the plan is ready to implement against.

**Next artifact: `src/mm_escape/io.py`.** Validate on `MMRF_1695` (33538 build),
`27522_1` (33694 build) and `BM4` (normal-BM control) before scaling to all 61 — those
three cover the failure modes. It is also the first genuine integration test of
`gene_space.py`, which has never touched a real count matrix.

This file will keep growing the way the R build's did — exact numbers, bugs found and
fixed, open decisions — as each stage actually runs.
