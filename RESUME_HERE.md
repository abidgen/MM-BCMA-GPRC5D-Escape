# RESUME HERE — MM Dual-Antigen pipeline (Python rebuild), session state

**Last updated:** 2026-08-20
**Branch:** `main` (working tree clean — R build removed, tagged `r-build-snapshot`)

Read `CLAUDE.md` first for the settled decisions and data ground truth. This file
covers only *where execution stands* and *what to do next*.

---

## TL;DR

Clean start. Previously built substantially in R (data acquisition solved and
verified on all 62 samples, QC/doublet-removal run on the full 61-sample cohort,
integration not yet run) before switching to Python. **None of the R code is being
ported.** All dataset knowledge carries forward via `CLAUDE.md`.

**Nothing has been run in Python yet.** The tree is ready for it: `raw/` intact at
62 samples, `scripts/01-03` in place, nothing else in the way.

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
  DN population a real subclone or scattered noise; pseudobulk DE; pre-registered
  γ-secretase hypothesis). Cell-cell communication moved 09 → **11**; decision
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

1. **Re-run `scripts/01-03`** to confirm `raw/sample_manifest.csv` still comes out
   clean (62 samples, no INCOMPLETE) — a no-op confirmation, not new debugging.
2. **Scaffold the repo**: `src/mm_escape/` package skeleton, `envs/env-qc.yml`,
   `envs/env-core.yml`, `envs/env-communication.yml`, `notebooks/`.
3. **Build `env-qc`**, register its Jupyter kernel (`mm-qc`).
4. **Write `src/mm_escape/io.py`**: the loader replacing `scanpy.read_10x_mtx()`.
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

### Presentable stopping points

The stage ordering was chosen so the project is showable partway through rather than
all-or-nothing:
- **After 04-08:** a working escape-fraction ranking.
- **After 09:** a ranking that survives hostile questioning.
- **After 10:** the actual scientific finding (subclone vs. noise).
- **After 11-12:** the full decision packet.
- **Then** S1 retrieval (un-flag provisional numbers, add longitudinal + cytogenetics),
  and only then Phase 2.

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
  headline metric is reported as an interval with a threshold sensitivity band, and
  the claim is ranking stability, not any single value.
- Patient-ID mapping is unresolved (47 vs. 41 patients) — fixed against
  Supplementary Table S1 before stage 08's per-patient aggregation is *final*, but
  per the 2026-08-20 policy this **no longer blocks execution**: proceed
  provisionally with provisional labelling in the output.
- Supplementary Table S1 itself is still not in the repo.
- `infercnvpy` is required (not optional) and normal-BM samples are controls
  (not filler).
- Stage 11 uses patient — not cell — as the unit of replication.
- Malignant subclustering is per-patient and un-integrated; Harmony is for the
  immune compartment only.
- GSE117156 is the confirmed Phase 2 validation dataset (own `phase2_NN_*`
  numbered notebooks, never merged with the primary cohort), sequenced strictly
  after Phase 1 (stages 01-12) completes; GSE118900 and He et al. 2022 were both
  evaluated and rejected as data sources.
- Notebooks and `results/` subdirectories are numbered `04`-`12`, one-to-one,
  continuing on from `scripts/01-03`. `src/mm_escape/` modules are named by
  function, not numbered.

## Status

Working tree clean, architecture decided, **scope expanded and all five `.md` files
updated (2026-08-20)**, no Python executed yet. This file will keep growing the way
the R build's did — exact numbers, bugs found and fixed, open decisions — as each
stage actually runs.
