# RESUME HERE — MM Dual-Antigen pipeline (Python rebuild), session state

**Last updated:** 2026-08-24 (S1 + stage 04)
**Branch:** `main`, everything merged, no feature branches open, working tree clean.

**To resume — three commands, then read on:**

```bash
cd /media/wrath/CART_mm_dual_antigen
git status --short                   # expect no output (clean tree)
git branch                          # expect only `main`
conda run -n mm-core pytest -q       # expect 155 passed, 2 skipped, ~27 s
ls results/04_qc/samples | wc -l     # expect 62 (stage 04 checkpoints)
ls -la results/05_integration/integrated.h5ad   # expect ~1.3 GB
```

If those are green the environment and data are intact and nothing needs rebuilding.
Then go to "Immediate next actions" — item 5 is the start point.

**Note:** `main` is ~10 commits ahead of `origin` and has not been pushed. Everything
is committed locally, so nothing is lost, but the work exists only on this machine
until `git push`.

Read `CLAUDE.md` first for the settled decisions and data ground truth. This file
covers only *where execution stands* and *what to do next*.

---

## TL;DR

Clean start. Previously built substantially in R (data acquisition solved and
verified on all 62 samples, QC/doublet-removal run on the full 61-sample cohort,
integration not yet run) before switching to Python. **None of the R code is being
ported.** All dataset knowledge carries forward via `CLAUDE.md`.

**Stages 01-03 are complete, and so is stage 04's loader.** All three acquisition
notebooks run green (62/62 `triplet-ok`, manifest byte-identical via CLI and notebook).
The **five** conda envs are built with kernels registered (`mm-integration` joined on
2026-08-24 for stage 05b). `src/mm_escape/` holds `config.py`, `gene_space.py`, `io.py`,
`qc.py`, `integration.py` and `benchmark.py`, covered by a **155-test** suite in
`tests/`. The whole cohort — **62 samples, 204,040 pre-QC cells** — loads in ~2 s and
harmonizes to 32,991 genes with all 65 required genes present.

**Stage 04 is now done too, and Supplementary Table S1 has landed.** 204,040 pre-QC
cells -> **172,940 kept (84.8%)**, per-cohort MAD thresholds derived and documented,
`scDblFinder` run on all 62 samples, one checkpoint per sample under
`results/04_qc/samples/` with every barcode retained and `obs["keep"]` set. S1 closed
the patient mapping: **41 patients over 53 in-cohort samples**, exactly the paper's
numbers, and confirmed the `_N` suffixes are serial disease-course timepoints.

**Stage 05 is done too.** 172,940 cells x **32,991 genes**, 30 Leiden clusters,
`results/05_integration/integrated.h5ad`. Harmony converged in 4 iterations. The gene
space came out exactly as predicted: 22,164 on symbols -> 32,991 on Ensembl IDs.

**Stage 05b (integration benchmark) is also done.** Seven arms scored with
`scib-metrics`; **no arm qualified and Harmony's stage-05 configuration stays**. The
arms that scored best on conventional scIB were the ones merging the cohort-censored
plasma populations (up to 20x), which is what the immune-only scoring rule existed to
catch. Stage 05's output is unchanged.

**>> The next artifact is `src/mm_escape/annotation.py` + `notebooks/06_annotation.ipynb`, in `mm-annotation`. <<**
Nothing from stage 06 onward exists.

Stage 04 turned up two things that change stage 08 and are written up in `CLAUDE.md`:
`pct_counts_in_top_20_genes` cannot be used as a filter here (it deletes
antigen-positive plasma cells), and the deposit is pre-filtered differently per
cohort, with a WashU 10,000-UMI ceiling that censors a band enriched 20-70x for
`GPRC5D`. Read the stage-04 entry in `CLAUDE.md` before writing stage 08.

Five real defects have surfaced so far, all fixed:
1. Cross-reference HGNC symbol drift — solved properly via Ensembl-ID reconstruction
   (32,991 genes vs. 22,164, and it prevents *mis-pairing*, not just loss).
2. A wrong bulk RNA-seq inventory in `CLAUDE.md` (30 usable -> 29, of which 26 match).
3. A circular dropout-correction formula in the stage 08 plan.
4. `56203_1` excluded for years on a misdiagnosis — it is a truncated deposit, now
   repaired on read and retained.
5. `ND_*` counted as disease rather than donors, which is where the "47 vs 41 patients"
   gap mostly came from.

Details below, newest session last.

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
- `raw/unpacked_bulk/` holds **29 usable GSE223061 bulk samples**, **26** of which
  have an exact scRNA match (18 MMRF TPM tables
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
4. ~~**Write `src/mm_escape/io.py`**~~ — **DONE 2026-08-24.** Loads all 62 samples in
   ~2 s; validated on the four failure-mode samples and end-to-end into `gene_space`.
   `pyproject.toml` added so `import mm_escape` works (`pip install -e . --no-deps`,
   done for `mm-qc` and `mm-core`; `mm-annotation`/`mm-communication` need it when
   stages 06 and 11 arrive). `tests/` added — 89 pass in `mm-core`.
5. ~~**Write `src/mm_escape/qc.py`** and `notebooks/04_qc.ipynb`~~ — **DONE
   2026-08-24.** Ran on all 62 samples: 204,040 -> 172,940 cells. Per-cohort MAD
   thresholds, one-sided `pct_counts_mt`, `scDblFinder` over the rpy2 bridge,
   per-sample checkpoints that keep every barcode. `tests/test_qc.py` adds 23 tests,
   most of them data-free. See "2026-08-24 — S1 lands, stage 04 runs" below for the
   two findings that came out of it.

6. ~~**Write `src/mm_escape/integration.py`** and
   `notebooks/05_integration_clustering.ipynb`~~ — **DONE 2026-08-24.** 172,940 x
   32,991, 30 clusters, ~190 s, **peak ~20 GB RAM** (the one stage that concatenates
   the matrices — it is the project's machine-size constraint). Donors and `25183`
   are kept in the embedding and excluded at the patient-level aggregation instead.
   See "2026-08-24 — stage 05" below for the compartment-specific integration result.

7. **>> START HERE << Write `src/mm_escape/annotation.py`** and
   `notebooks/06_annotation.ipynb` -> `results/06_annotation/`, in **`mm-annotation`**.

   - **First: `pip install -e . --no-deps` in `mm-annotation`.** It has only been done
     for `mm-qc` and `mm-core`, so `import mm_escape` will fail there.
   - **Three methods, compared, with the F1 thresholds declared before looking** —
     PlasmaCell 0.95, T/NK/myeloid 0.90, rest 0.85. Write
     `results/06_annotation/annotation_decision.md` with the per-class table.
   - **The marker-coverage test is the load-bearing evidence, not concordance.** A
     failed marker test vetoes a class regardless of how well the methods agree.
   - **Interface contract**: `cell_type` (the only load-bearing output),
     `cell_type_fine`, `annotation_source`, `annotation_conf`, and
     `config.ANNOTATION_DECISION`.
   - **State programs are continuous scores**, never `cell_type` labels.
   - Stage 05 split plasma cells into three cohort-specific clusters (see below).
     Harmless for annotation, but run the plasma-cell marker-coverage check **on
     myeloma marrows specifically**, not only the donors.

8. Continue through stages 06-12 per `CLAUDE.md`'s pipeline section — notebook
   number N always writes to `results/N_*/`, nothing else. Run them in numeric
   order; number order is execution order.
9. **Apply the review corrections** as each stage is written — they are documented in
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
*(superseded 2026-08-24: computed from the GEO metadata as exactly 26)*
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

## 2026-08-24 — `src/mm_escape/io.py` written and validated on the real data

Merged `review3-stage08-corrections` into `main` (fast-forward), then built the loader.

### What it does

`load_manifest()` -> `read_sample()` / `read_samples()`. The manifest's column schema
(from `scripts/03_build_manifest.py`) is the contract; the loader never auto-detects a
10x directory, because `read_10x_mtx()` hardcodes filenames this deposit does not use.
`load_manifest` also resolves the manifest's repo-root-relative paths to absolute,
checks they exist up front, and adds `gsm_id`, `sample_name`, `patient_id`,
`patient_id_source`, `sample_type`, `sample_type_certain`, `excluded`.

`read_sample` returns cells x genes, CSR float32, `var_names` = the **deposited symbols
in deposited order, untouched** (that is `attach_ensembl_ids`'s precondition),
`obs_names` = `<sample_name>_<barcode>` (raw barcodes repeat across samples), and
`n_genes_ref` in `.obs` so it survives concat and can be Harmony's covariate at stage 05.

### Validated against the real files, not just imported

Three samples covering the failure modes — `MMRF_1695` (33538), `27522_1` (33694),
`BM4` (normal-BM control):

- **The transpose.** The `.mtx` is genes x cells (`33538 1007 2604146`); AnnData is
  cells x genes. Checked by re-reading raw triplets and asserting
  `X[cell-1, gene-1] == count`, plus exact equality of total counts (13,357,462) and
  nnz (2,604,146). A transposed object still looks plausible, so this is asserted, not
  eyeballed.
- **The io -> gene_space handoff** (the point of the exercise — `gene_space.py` had
  never touched a real count matrix): `attach_ensembl_ids` accepted all three
  position-for-position, `intersect_gene_space` gave **32,991** genes as predicted,
  `to_canonical_symbols` + `assert_required_genes` passed with all 65 required genes
  and 11,140 drifted symbols correctly joined (`NSD2` in 33538 = `WHSC1` in 33694).
- **Failure paths raise**: excluded sample by name and by row, unknown sample,
  a non-GSM `sample_id`, an unknown name in `read_samples`, and — via `gene_space` —
  a reordered gene axis.
- **Scale**: all 61 retained samples load in ~4 s. **202,203 pre-QC cells**
  (154,053 disease / 48,150 normal-BM); 509 / 2,767 / 9,328 min/median/max cells per
  sample, an **18.3x** spread. (The R build's 480 / 2,555 / 7,937 ~15x figures in
  CLAUDE.md are *post*-QC — different quantity, both correct.)

### Real finding: the 47-vs-41 patient gap is mostly the `ND_*` samples

*(Superseded 2026-08-24 — GEO confirms `ND_*` are donors; the real figure is 43
patients from 54 myeloma samples. Kept as the session record.)*
The inherited "naive rule yields 47 patients from 57 disease samples vs. the paper's
41 / 53" counts the four `ND_*` samples as **disease**. CLAUDE.md itself treats them as
controls in two other places, including stage 07's negative control. Counting them
consistently as controls (and excluding `56203_1`):

| | samples | naive patients |
|---|---|---|
| `ND_*` as disease (inherited) | 57 | 47 |
| **`ND_*` as controls** | **53** | **43** |

53 is the paper's disease-sample count exactly, and 43 is **two** collapses short of 41
rather than six — `83942`/`MMY83942` being the obvious remaining pair. The four names
(`ND_083017`, `ND_090617`, `ND_170531`, `ND_170607`) carry collection **dates**, not
patient IDs, which is weak independent support.

**Not resolved — a hypothesis.** `ND` could read "newly diagnosed"; the deposit does not
say; and the 53 also leans on dropping `56203_1`, which the paper had no reason to drop
for our BCMA-reference reason. So `io.py` encodes the uncertainty in the data rather
than in a comment: `sample_type == "normal_bm"` with `sample_type_certain == False`.
S1 still settles it, and `patient_id_source == "naive"` still rides on every cell.

### `pyproject.toml` added — `import mm_escape` did not work in any env

The package was never installed anywhere, so a stage-04 notebook could not have
imported the loader without a `sys.path` hack. Added a minimal `pyproject.toml`
(setuptools, src layout) with **`dependencies = []` on purpose** — `envs/*.yml` is the
dependency manifest, and pip must never resolve anything into these envs (that is what
downgraded numpy 2.5.2 -> 1.26.4 during setup and broke four packages at once).

Installed into `mm-qc` and `mm-core` only, and only as:

    conda run -n <env> pip install -e . --no-deps

Both verified afterwards: `mm_escape.io.load_manifest()` returns the full cohort, and
`mm-core` still reports numpy 2.5.2 / scanpy 1.11.5 / scipy 1.18.0 / pydeseq2 0.5.4 /
decoupler 2.2.0 / numba 0.67.0 — untouched. `mm-annotation` and `mm-communication`
need the same one-liner when stages 06 and 11 arrive.

### Next artifact: `src/mm_escape/qc.py` + `notebooks/04_qc.ipynb`

MAD-based outlier calling (5 MADs on `log1p_total_counts`,
`log1p_n_genes_by_counts`, `pct_counts_in_top_20_genes`, plus a cohort-derived
`pct_counts_mt` cap — **re-derive, do not copy the tutorial's 8%**), then `scDblFinder`
over the `rpy2` bridge in `mm-qc`, checkpointing each sample's post-QC AnnData
individually. The loader underneath it is done and green.

---

## 2026-08-24 (same session) — GEO metadata arrived and rewrote three ground truths

Six GEO metadata files appeared in `raw/`. They are ~8 KB each and they overturned more
of the inherited ground truth than anything since the gene-space work. All committed.

### `56203_1` recovered — the exclusion was a misdiagnosis

    counts.mtx header:            33694 1837 2135520   <- normal 33694 matrix
    genes.tsv:                    22185 rows, ends 'KBTBD', no trailing newline
    33694 reference row 22185:    'KBTBD7'
    rows 1..22184 vs reference:   identical (strict prefix)

The gene-file write failed part-way. `TNFRSF17` (row 25539) and `IGLC1/2/3` (rows
32548-32552) were **past the cut, not absent from a reference**. The "22184 genes"
figure was a `wc -l` artifact of the missing trailing newline.

`io.read_sample` substitutes the canonical column from the committed gene map behind a
prefix assertion — every written row must match, and the final partial row must be a
prefix of the symbol it was cut from, or it raises. Both failure paths tested.
**`EXCLUDED_SAMPLES` is now empty. Cohort: 62 samples, 204,040 pre-QC cells.**

### `ND_*` settled — they are donors

GEO gives all four `ND_*` and all four `BM*` `source_name = "Donor BMMC"` and **no
`diagnosis` characteristic**; the other 54 read "Multiple myeloma (MM)". The morning's
hypothesis is now fact, and `sample_type_certain` is True everywhere. Naive mapping is
**54 samples / 43 patients** vs. the paper's 53 / 41 — two collapses short, not six.

### Cohort/chemistry — and a number I had to walk back

`resources/sample_metadata/{scrna,bulk}_samples.tsv` are committed (parsed by
`io.rebuild_sample_metadata_from_soft`); `load_manifest` joins them, so every cell
carries `cohort`, `chemistry`, `dead_cell_removal`, `diagnosis`.

**I first said v2-vs-v3 was a 2-3x sensitivity gap. Measured, it is 1.38x with
overlapping distributions.** Sample-level medians of genes per pre-QC cell:

    MMRF   v3.3   1916      WU2   v3.2   1210
    Donor  v3.2   1103      WU1   v2     1023

    v2 vs all-v3: 1023 vs 1408 = 1.38x, Mann-Whitney p = 6.5e-05
    v2 max 1602 > v3 min 793, so the distributions OVERLAP

The axis that separates is **cohort** (MMRF ~1.9x the rest), of which chemistry is one
component alongside site and protocol. Still must be modelled — a 1.9x depth spread
tracking cohort will move a fraction-of-zeros metric — but **do not quote a "2-3x
chemistry effect"**; CLAUDE.md now says so explicitly.

**`n_genes_ref` cannot stand in for it**: the build split cuts across cohorts (two WU1
samples on 33538, the four `ND_*` on 33694). Stage 05 needs both covariates.

### Three more corrections, all in CLAUDE.md

- **Bulk overlap is 26 exact matches**, computed, not the inherited "~28". And the two
  bulk cohorts are different assays: **MMRF is CD138+ sorted, WashU 1 is unsorted
  BMMC**. Pooling them would make 10 of 26 comparisons measure tumour burden instead
  of antigen abundance. Stage 09 pairs sorted bulk with malignant pseudobulk and
  unsorted bulk with whole-sample pseudobulk.
- **Raw data exists under controlled access** (dbGaP `phs000159` / `phs000748`), not
  "does not exist". Conclusion unchanged — still no unfiltered matrices — but the
  claim was too strong.
- **The `_N` suffix reading reversed.** Bulk suffixes are always a subset or overlap of
  the scRNA ones, never a different scheme, and `37692_2` / `57075_3` are *lone*
  samples with non-`_1` suffixes — which a fraction/sort/replicate label would not be.
  A serial per-patient index explains all of it. Not proof of *timepoint*, S1 still
  settles it, but the longitudinal arm is worth planning for now.

Caveat on all of the above: **the deposit's own metadata is not self-consistent** — it
claims Cell Ranger v3.0.0 for all 62 samples, which the files contradict for 24. Every
claim taken from the SOFT files was verified against the data before being written down.

### Test suite added — `tests/`, run in `mm-core`

Two-tier, because the loader's data lives on one machine and the invariants worth
protecting mostly do not need it:

    pytest                     89 pass,  1 skip   (~7 s, deposit present)
    pytest (no raw/)           51 pass, 39 skip   (~0.5 s, fresh clone)
    pytest -m "not slow"       skips the two full-cohort passes

Data-backed tests are gated on `conftest.requires_data` and **skip rather than fail**.
They run over `conftest.CANONICAL_SAMPLES` — `MMRF_1695` (33538), `27522_1` (33694,
spells NSD2 as WHSC1), `BM4` (donor), `56203_1` (truncated) — which cover the
deposit's failure modes. Add cases there rather than picking samples ad hoc.

Assertions pin the known invariants (204,040 pre-QC cells, 32,991 intersected genes,
11,140 drifted symbols, 26 matched bulk samples, the per-cohort depth ordering), so a
regression in the loader, the gene map or the metadata tables fails a test instead of
quietly changing a result.

Two fixes the suite surfaced: `scipy.io.mmread` was called without `spmatrix=`, which
scipy 1.20 flips to sparse *arrays* (now pinned to `spmatrix=True`, since the
scanpy/anndata stack is spmatrix-native); and `pytest` was only in `mm-core`
transitively, now declared in `envs/env-core.yml`.

### Next artifact: `src/mm_escape/qc.py` + `notebooks/04_qc.ipynb`

Unchanged, with one addition: **MAD thresholds are derived per cohort, not pooled.** A
pooled MAD across a 1.9x depth spread would flag much of WashU cohort 1 as low-quality
for a batch reason. The donors span both reference builds, so stage 07's negative
control doubles as a build control.

---

## 2026-08-24 — S1 lands, stage 04 runs

Supplementary Tables S1-S6 appeared in `raw/`. Two pieces of work came out of it:
S1 parsed and wired in, then stage 04 written and run end to end.

### S1 closed the patient mapping, exactly

The naive rule gave 43 patients over 54 deposited myeloma samples; the paper reports
41 / 53. Two independent corrections close both gaps:

- **`25183`** is deposited (scRNA *and* bulk) but appears in **no** supplementary
  table. That is the whole 53-vs-54 gap. It is **not dropped** — `in_paper_cohort ==
  False`, `clinical_source == "none"`, so an aggregate excludes it on purpose.
- **`83942` and `MMY83942` are one patient** — S1 lists them separately but with
  identical age/sex/race/ISS/treatment, sampled under both WashU protocols.

`io._assert_s1_reproduces_the_paper` checks all three counts off the committed GEO
table, so a revised S1 fails loudly rather than quietly moving the denominator of
`frac_double_negative`.

**The `_N` suffixes are serial disease-course timepoints.** S1 sheet 2:
`27522_1` Primary -> `_2` Remission-1 -> `_3` Relapse-1 -> `_4` Relapse-2 ->
`_5` Remission-2 -> `_6` Relapse-3. The longitudinal arm is real, not speculative.
Coverage is **WashU cohort 1 only** — MMRF and WU2 get `disease_stage = NA` and one is
not imputed. **S1 carries no cytogenetics at all**, so stage 10's TC proxy still has
nothing in this deposit to validate against.

New: `resources/sample_metadata/patients_clinical.tsv` (43 rows) and
`sample_disease_stage.tsv` (22). `load_manifest` emits age/sex/race/iss_stage/
treatment/ttpd_months/disease_stage/disease_phase/timepoint/clinical_source/
in_paper_cohort, and every cell carries them.

All six supplementary tables are committed (~500 KB). Table S3 and Table S5 (file
`s6`) are useful later — see `CLAUDE.md`. Note file `s5` is strict-OOXML so
`openpyxl` reports zero sheets for it; it is not corrupt.

### Stage 04 ran, and found two things that change stage 08

**204,040 -> 172,940 cells (84.8% kept)** over all 62 samples. Two passes: per sample
(metrics + `scDblFinder`, ~210 s) then per cohort (MAD thresholds, ~11 s, computed off
`obs` alone so the stage never concatenates the matrices).

**1. `pct_counts_in_top_20_genes` cannot be used as a filter in this tissue.** A
5-MAD band flags 17% of MMRF and 15% of WU1 against 3% of WU2. Those cells are two
populations — `IGKC` at ~25% of counts (plasma cells) and `HBB`/`HBA1/2` at ~32%
(erythroid debris) — and the plasma-cell half is the project's subject: `TNFRSF17`
detected in **21.8%** of the flagged decile vs **0.8%** elsewhere. An Ig-dominated
library is a plasma cell's normal state. Filtering on it deletes antigen-**positive**
malignant cells and inflates the escape fraction. The flag is computed and reported
(it is one of the few ambient-Ig handles available, since SoupX needs unfiltered
matrices this deposit lacks) but is not in `qc.DEFAULT_FILTERS`.

**2. The deposit is pre-filtered, differently per cohort — correcting an earlier
claim in this file.** The note that the depositors' 10,000-UMI cut "was not applied"
came from a cohort-wide average that pooled MMRF with WashU. Per cohort:

    WU1, WU2   UMI < 10,000   UMI >= 1,000   pct_mt < 20%   genes >= 200
    MMRF       uncensored     UMI >= 1,000   pct_mt < 10%   genes >= 200
    Donor      uncensored     uncensored     pct_mt < 20%   genes >= 200

The WashU ceiling is the consequential one. Malignant plasma cells are the
highest-RNA-content cells in marrow, so measured in the uncensored cohorts the band
above 10,000 UMIs is enriched **3-21x for `TNFRSF17`** and **20-70x for `GPRC5D`**.
**36 of 54 myeloma samples had the antigen-positive tail of their own tumours removed
before deposit**, inflating `frac_double_negative` for WU1/WU2 — biased toward the
project's hypothesis, and unfixable. Not corrected at stage 04 (that would mean
discarding 42% of MMRF's cells); **stage 08 owes a truncate-all-at-10k sensitivity
analysis**.

Thresholds produced (`results/04_qc/qc_thresholds.csv` has all of them):

| cohort | log1p_total_counts MAD | pct_mt cap | % removed |
|---|---|---|---|
| MMRF | 0.94 (wide — removes doublets only) | 12.6% | 4.7% |
| WU1 | 0.41 | 8.6% | 14.1% |
| WU2 | 0.45 | 12.9% | 18.0% |
| Donor | 0.37 | 11.7% | 22.6% |

Stable from 4 to 6 MADs (16.5% -> 14.6% overall), so nothing hangs on the exact count.

### Test suite

**117 passed, 2 skipped in `mm-core`** (was 89/1); **70 passed, 49 skipped on a fresh
clone with no `raw/`** (was 51/39). Two new gates in `conftest.py`: `requires_s1`
(the xlsx is a journal file, not part of the GEO deposit) and `requires_r`
(scDblFinder lives only in `mm-qc`, and `mm-core` carries no R). `pytest` was added to
`envs/env-qc.yml`, but the currently-built `mm-qc` predates that line — the bridge is
exercised by `notebooks/04_qc.ipynb` over all 62 samples, which is the stronger check
anyway.

---

## 2026-08-24 — stage 05: integration, and the censoring shows up again

**172,940 cells x 32,991 genes, 30 Leiden clusters, ~190 s, peak ~20 GB RAM.** That
peak is the project's machine-size constraint — stage 05 is the only stage that
concatenates the count matrices (stage 04 deliberately worked off `obs` alone).

The gene space came out exactly as `gene_space.py` predicted on four samples:
**22,164 genes on raw symbols -> 32,991 on Ensembl IDs (+10,827)**, 11,140 drifted
symbols joined correctly, `NSD2` resolving against `WHSC1`. Harmony converged in 4
iterations.

### Harmony fixed the immune compartment and did not fix the plasma cells

| | clusters | cells | median cohort-mixing entropy |
|---|---|---|---|
| plasma-cell-like (`MZB1` > 40%) | 11 | 39,893 | **0.105** |
| everything else | 19 | 133,047 | **0.751** |

Uncorrected PCA for reference: 0.341 over 54 clusters. Harmony overall: 0.621 over 30.
So the correction does real work — unevenly.

The three largest plasma-cell clusters are **one per cohort**, each spanning ~30
patients. That kills the benign reading: a patient-private clone would fragment into
~41 clusters, not three cohort-shaped ones. The likely cause is stage 04's finding —
WashU was cut at 10,000 UMIs and MMRF was not, and plasma cells are the
highest-RNA-content cells in marrow, so **WashU's plasma cells are a truncated subset
of the plasma-cell distribution**. No correction method restores cells that were never
deposited. It is compartment-specific for the same reason: T/NK/myeloid/B cells sit
well below 10,000 UMIs in every cohort, so the ceiling never touched them.

**Contained, not fatal**, and contained by decisions made before it was observed:
antigen calls are raw counts and never touch the embedding; stage 10 is per-patient
and un-integrated; stage 06 annotates at cluster level, where three plasma-cell
clusters all annotate as PlasmaCell. **What it forbids: reading any cross-cohort
comparison of malignant-cell state off this embedding.** And the
truncate-all-at-10,000 sensitivity analysis stage 08 owes is now owed twice — two
independent signs of one problem.

### Preprocessing diagnostics added, and they paid for themselves

Prompted by "don't we need normalization, feature selection and dim reduction?" — all
three had run (inside `normalize_and_hvg` and `run_pca_harmony`) but **none was
reported**: the notebook had one line printing an HVG count. Added a section covering
all three, read-only over the cached `integrated.h5ad` (verified byte-identical after
re-execution). Two findings came out of it:

**`GPRC5D` is not a highly variable gene** — mean 0.061 vs `TNFRSF17`'s 0.492 (8x),
HVG in 6 of ~50 patients. Affects nothing (embedding does not need it; stage 08 reads
raw counts) but it is the first number from this cohort behind the standing claim that
GPRC5D is dropout-prone. **GPRC5D-negative calls deserve more scepticism than
BCMA-negative ones.** The panel is deliberately not forced into the HVG set.

**The plasma-cell integration failure is the stage-04 censoring, measured:**

    compartment    MMRF     WU1     WU2    MMRF/WU1
    non-plasma     5,829   3,273   2,879     1.8x
    plasma-like   22,477   5,036   4,888     4.5x

MMRF's two biggest plasma clusters are 68% and 88% above 10,000 UMIs — cells WashU
cannot contain. WashU's press against the ceiling instead. Harmony is not failing: what
separates the compartments is a **non-recoverable sampling/censoring asymmetry**, not
established cohort biology — WashU's observed plasma distribution is missing its
high-RNA portion, so no one-to-one correspondence remains to recover. Two regression
tests pin this. (Earlier wording said "the populations differ", which implied a
biological claim the data does not support.)

### One defect found and fixed

`gene_space.to_canonical_symbols` named the `var` index `canonical_symbol` while also
keeping a `canonical_symbol` **column** holding the unsuffixed symbol — which differs
for the 9 collision-resolved genes. AnnData refuses to write such an index, so this
failed only at `write_h5ad`, after every in-memory test had passed. The index is now
named `symbol`, and `tests/test_integration.py` round-trips through `.h5ad` so a
serialization-only bug cannot hide again. **Lesson for later stages: anything that
only fails on write needs a test that writes.**

Test suite: **131 passed, 2 skipped** in `mm-core` (133 collected); **76 passed, 57
skipped** on a fresh clone with no `raw/` and no `results/`.

Note `tests/test_integration.py` already existed — it held 11 io->gene_space
end-to-end tests — and writing the stage-05 tests over it destroyed them. Recovered
from `db026e8` as **`tests/test_io_gene_space_e2e.py`**, which is the better name
anyway: "integration" there means the *software* sense, and in `test_integration.py`
it now means the *biological* sense (Harmony). The collision is what caused the
clobber. Every test file now maps 1:1 to a `src/mm_escape/` module, plus that one
cross-module file.

---

## 2026-08-24 — stage 05b: the integration benchmark, and what it caught

Stage 05 used Harmony because it was the default. `sc-best-practices`' integration
chapter says to run several methods and score them with scIB rather than assume, so
that gap is now closed. **The incumbent survived — and how it survived is the
interesting part.**

### Setup

`envs/env-integration.yml` -> **`mm-integration`**: every integration method under
comparison plus the scoring stack in one env (harmonypy, scvi-tools, scanorama, bbknn,
scib-metrics, celltypist). torch 2.13.0+cu130 on the RTX 5070 (sm_120) works; anndata
0.13.2 matches `mm-core` exactly. Two build traps: **`cxx-compiler` is load-bearing**
(Scanorama needs `annoy`, which has no cp312 wheel and builds from source; there is no
system g++ and conda-forge's `python-annoy` is py39-only), and **`bbknn` is installed
but not scored** (it yields a graph, not an embedding, so `Benchmarker` cannot place it
on the same footing).

### The result

Seven arms, scored on the **immune compartment** against provisional CellTypist labels,
all scored on `cohort` regardless of what they corrected on:

    arm                batch    bio    depth R2   plasma mixing   vs incumbent
    unintegrated       0.450   0.691     0.660        0.014       depth +0.291
    harmony_sample     0.615   0.718     0.509        0.515       depth +0.140, plasma 13.5x
    scvi_sample        0.570   0.701     0.541        0.452       depth +0.172, plasma 11.8x
    scanorama_sample   0.450   0.723     0.690        0.161       depth +0.321
    harmony_stage05    0.427   0.700     0.369        0.038       -- (incumbent)
    harmony_cohort     0.591   0.706     0.607        0.771       depth +0.238, plasma 20.2x
    scvi_cohort        0.492   0.690     0.576        0.017       depth +0.207

**No arm qualified. `harmony_stage05` stays.**

**The arms that win on conventional scIB are precisely the arms that merge the censored
plasma populations.** `harmony_sample` posts the best batch *and* bio scores while
mixing plasma **13.5x** harder than the incumbent and encoding more depth;
`harmony_cohort` reaches **20.2x**. The two arms that leave those populations apart
(`unintegrated`, `scvi_cohort`) are the ones with no real batch gain. **A standard
global scIB benchmark would have picked `harmony_sample`** — which buys its score by
fusing populations that cannot be fused. That is exactly the failure mode the scoring
design was built to catch, and it happened.

The incumbent is simultaneously the **worst batch corrector** (0.427, below
unintegrated's 0.450) and by a wide margin the **least depth-encoding** (R² 0.369 vs
0.51-0.69) and **least plasma-merging** (0.038).

**Stated honestly:** `depth_ok` did all the gating — every non-incumbent arm failed the
+0.05 tolerance, which was fixed before the 0.37-0.69 spread was known. Two things stop
that being a threshold artifact: `harmony_sample`, `scvi_sample` and `harmony_cohort`
**independently fail `overcorrection_ok`** and lose even with depth removed entirely;
and relaxing depth enough to admit anything admits only the two weakest-batch arms. The
tolerance is **not** re-tuned after the fact.

**Bonus finding: scVI encodes depth in plasma cells.** Plasma R²(depth ~ latent) is
0.793 / 0.850 for the scVI arms against the incumbent's 0.528 and `harmony_sample`'s
0.319 — its explicit library-size model appears to put depth *into* the latent space for
the one compartment where depth is the confound, the opposite of why it was the
principled candidate.

### Design points worth not re-deriving

- **Immune scored, plasma diagnosed.** Batch metrics cannot tell "correctly left apart"
  from "failed to merge"; plasma mixing never contributes positively.
- **`R²(depth ~ latent)`, fixed in advance**, because R² depends only on the column span
  and is rotation-invariant — latent axes are arbitrary across methods.
- **Labels are CellTypist with `majority_voting=False`** — voting uses an
  over-clustering, which would smuggle an embedding back into the labels. `ILC` was
  verified to be NK (NKG7 98.7%), and `Immune_All_High` *does* cover erythroid/HSPC
  here, so no hand-set marker thresholds enter the benchmark.
- **scANVI/scGen deferred**, not rejected: they need stage-06 labels, and stage 06
  consumes the embedding under selection.
- **The benchmark cannot undo the censoring.** Stage 08 still owes its
  truncate-all-cohorts-at-10,000 sensitivity analysis whatever wins here.

Test suite: **155 passed, 2 skipped** in `mm-core`; **100 passed, 57 skipped** on a
fresh clone.

---

## Status

Stages 01-05 plus 05b complete and green. Envs built, kernels registered, gene space solved and
committed, S1 parsed, QC run on the full cohort. Architecture stable.

**Next artifact: `src/mm_escape/annotation.py` + `notebooks/06_annotation.ipynb`,
in `mm-annotation`** (which needs `pip install -e . --no-deps` first). Three methods
compared per class against F1 thresholds declared in advance, with the marker-coverage
test as the load-bearing evidence. First presentable state is still stages 04-08:
escape fractions with co-escape enrichment.

This file will keep growing the way the R build's did — exact numbers, bugs found and
fixed, open decisions — as each stage actually runs.
