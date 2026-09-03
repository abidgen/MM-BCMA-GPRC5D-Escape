# Environments — full specifications and build history

Detail split out of the main project document on 2026-08-24 to keep that file under its
size budget. That document keeps the env→stage table and the two hard rules; everything
below is the full record. **Five environments are built; six YAML specifications exist** —
`envs/env-composition.yml` (scCODA) is written but deliberately not built. The five built
bodies follow, what the channels actually
carried versus what the specs asked for, the verification results, and the build
traps that were hit for real.

**Five environments exist.** `mm-qc`, `mm-core`, `mm-annotation`,
`mm-communication` and `mm-integration` are built with kernels registered;
`envs/env-composition.yml` is written but deliberately **not** built (see the
scCODA note below).

---


The core four, split by actual dependency-conflict risk (not one per file);
`mm-integration` was added later for stage 05b, bringing the built total to five. Two of
them exist to quarantine R — `env-qc` for `scDblFinder`, `env-annotation` for
`SingleR` — so the pure-Python stack never carries an R dependency it doesn't need:

**`envs/env-qc.yml`** — stage 04 only (loading, QC, doublet detection). Isolates
the R/rpy2 bridge so no other environment needs to carry R at all.
```yaml
name: mm-qc
channels: [conda-forge, bioconda]
dependencies:
  - python=3.12
  - scanpy=1.11
  - anndata
  - pandas
  - numpy
  - scipy
  - seaborn
  - matplotlib
  - jupyterlab
  - ipykernel
  - r-base=4.3.3
  - rpy2=3.5.11
  - bioconductor-scdblfinder=1.16.0
  - bioconductor-singlecellexperiment
  - anndata2ri
```
(Adapted directly from `sc-best-practices`'s own published environment for this
exact chapter.)

**`envs/env-core.yml`** — stages 05-10 and 12 (integration, annotation, malignant
calling, antigen scoring/escape fraction, robustness, subclone/phenotype, decision
packet — everything except stage 04's QC and stage 11's LIANA+). Shared
dependencies, no conflicts between them.
```yaml
name: mm-core
channels: [conda-forge, bioconda]
dependencies:
  - python=3.12
  - scanpy=1.11
  - anndata
  - harmonypy
  - leidenalg
  - python-igraph
  - celltypist
  - infercnvpy
  - pydeseq2          # installed in mm-core; not used for the frozen Stage-10 DE result
  - decoupler         # stage 10 pathway/TF activity (Hallmark/PROGENy/CollecTRI)
  - pandas
  - numpy
  - scipy
  - scikit-learn
  - statsmodels       # stages 08/09 regression, bootstrap CIs, confounder models
  - seaborn
  - matplotlib
  - jupyterlab
  - ipykernel
```

**`scCODA` is deliberately NOT in `env-core`.** It pulls TensorFlow, which is a
heavyweight and version-brittle dependency with no relationship to the rest of the
core stack — exactly the conflict-risk criterion these env splits exist to respect.
If the compositional analysis (stage 06) is actually run, give it
`envs/env-composition.yml` with `sccoda` + `anndata` + `ipykernel` only. Do not
destabilize `mm-core` to save an environment.

**`envs/env-communication.yml`** — stage 11 only (LIANA+). Isolated because
`liana`/`omnipath`'s dependency tree is version-sensitive and unrelated to the
core scanpy stack.
```yaml
name: mm-communication
channels: [conda-forge, bioconda]
dependencies:
  - python=3.12
  - liana
  - omnipath
  - anndata
  - scanpy
  - jupyterlab
  - ipykernel
```

**`envs/env-annotation.yml`** — stage 06 only. Isolated for the same reason `env-qc`
is: `SingleR` is R, and R stays quarantined in the environments that actually need it
rather than being pulled into `mm-core`.
```yaml
name: mm-annotation
channels: [conda-forge, bioconda]
dependencies:
  - python=3.12
  - scanpy=1.11
  - anndata
  - celltypist
  - r-base=4.3.3
  - rpy2=3.5.11
  - bioconductor-singler
  - celldex
  - anndata2ri
  - scikit-learn        # ARI / F1 for the annotation comparison
  - pandas
  - seaborn
  - matplotlib
  - jupyterlab
  - ipykernel
```
`celltypist` also remains in `env-core` — pure Python, no conflict, and convenient if
labels ever need re-deriving outside stage 06.

Register a distinct Jupyter kernel per env: `python -m ipykernel install --user
--name mm-qc` (and `mm-core`, `mm-annotation`, `mm-communication`).

**Written and built 2026-08-21** as `envs/env-{qc,core,annotation,communication}.yml`,
plus `envs/env-composition.yml` which is **written but deliberately not built** — per
the scCODA note above, it is created on demand only if stage 06's compositional
comparison is actually run. Three deviations from the specs above, all forced by what
the channels actually carry (verified against `conda-forge`/`bioconda`, not assumed):

| spec said | reality | resolution |
|---|---|---|
| `infercnvpy` (conda) | **not packaged on any conda channel** | installed via `pip:` inside `env-core.yml` — still REQUIRED, not optional |
| `celldex` | no such conda package | the R package is `bioconductor-celldex` |
| `decoupler` (conda) | bioconda is stuck at **1.5.0 (2023)**, pins `numpy<2`, and fails to import under numba 0.67 | `pip: decoupler==2.2.0` in `env-core.yml` |
| `liana` (conda) | bioconda recipe **under-declares**: `import liana` fails outright | pin `pydantic` + `mudata<0.4` in `env-communication.yml` |
| (unstated) | notebooks are jupytext-paired in every env | `jupytext` added to all four |

Every pinned version in the specs does exist and was confirmed before building:
`scanpy=1.11`, `r-base=4.3.3`, `rpy2=3.5.11`, `bioconductor-scdblfinder=1.16.0`.

**`decoupler` 2.x is an API rewrite** (`dc.mt.*` / `dc.op.*`, not `dc.run_mlm`). Stage 10
must be written against 2.x — do **not** follow 1.x tutorials, including
`sc-best-practices`'s, without checking the call names. Downgrading is not an option:
1.x pins `numpy<2`, which collides with `pydeseq2`/`scipy`/`zarr` in the same env.

**Never `pip install` into these envs casually.** Installing `decoupler==1.8.0` during
setup silently downgraded `numpy` 2.5.2 → 1.26.4 *and* `numba`, breaking `scanpy`,
`scipy`, `pydeseq2` and `zarr` at once; the repair was to delete and recreate the env
from the yml. The two legitimate pip entries (`infercnvpy`, `decoupler`) live in the yml
so a rebuild reproduces them. If a pip install is unavoidable, rebuild the env
afterwards rather than patching it.

**Verified 2026-08-21 — all four envs import their key packages, both R bridges work:**

| env | verified |
|---|---|
| `mm-qc` | scanpy 1.11.5, rpy2 3.5.11, anndata2ri 2.0.1; **R 4.3.3 + scDblFinder 1.16.0 + SingleCellExperiment 1.24.0**, `scDblFinder()` callable |
| `mm-core` | scanpy 1.11.5, numpy 2.5.2, harmonypy 2.0.0, celltypist 1.7.1, infercnvpy 0.6.1, pydeseq2 0.5.4, decoupler 2.2.0; `cnv.tl.infercnv` + `run_harmony` callable |
| `mm-annotation` | celltypist 1.7.1, rpy2 3.5.11; **R 4.3.3 + SingleR 2.4.0 + celldex 1.12.0**, both `NovershternHematopoieticData` and `HumanPrimaryCellAtlasData` present |
| `mm-communication` | liana 1.8.1 with **both `liana.mt.cellchat` and `liana.mt.rank_aggregate`**, omnipath 1.0.12 |

Kernels registered for all four (`mm-qc`, `mm-core`, `mm-annotation`, `mm-communication`)
— plus `mm-integration` on 2026-08-24 (stage 05b),
each kernelspec confirmed to point at its own interpreter.

**Note `mm-communication` deliberately runs a different scanpy/anndata** (1.12.3 /
0.12.19) from the other three (1.11.5 / 0.13.2) — the spec leaves scanpy unpinned there,
and pinning it would over-constrain liana's already version-sensitive tree. Stage 11
reads `.h5ad` written by `mm-core`, so the two only meet on disk; if a forward/backward
`.h5ad` compatibility problem ever appears, that version gap is the first place to look.

**`envs/env-integration.yml`** — stage 05b only (the integration-method benchmark).
**Built 2026-08-24**, fulfilling the reservation this section previously carried ("if
scVI-based integration is ever considered as an alternative to Harmony, it gets its own
separate env (`env-scvi.yml`) — not created yet, only if actually needed"). It is named
for what it holds rather than for one member: **every** integration method under
comparison lives here, plus the scoring stack — `harmonypy` (the incumbent, which must
run beside its rivals), `scvi-tools`, `scanorama`, `bbknn`, `scib-metrics`, and
`celltypist` for the benchmark's provisional labels.

Verified: torch **2.13.0+cu130**, CUDA available, RTX 5070 at compute capability
**12.0 (Blackwell / sm_120)**; scvi-tools 1.5.0, scib-metrics 0.6.0, scanorama 1.7.4,
bbknn 1.6.0, celltypist 1.7.1, harmonypy 2.0.0; scanpy 1.11.5 and **anndata 0.13.2,
matching `mm-core` exactly** so the stage-05 `.h5ad` reads without a version gap.

Two build traps, both hit for real:
- **`cxx-compiler` is load-bearing.** Scanorama depends on `annoy`, which publishes no
  cp312 wheel and builds from source. Without a compiler *in the env* the pip stage
  dies with `error: [Errno 2] No such file or directory: 'g++'` and takes the whole env
  with it. There is no system `g++` on this machine and conda-forge's `python-annoy` is
  py39-only, so the compiler has to come from the env itself.
- **`bbknn` is installed but not scored.** `scib-metrics`' `Benchmarker` scores `obsm`
  embeddings; BBKNN yields a corrected neighbour *graph* and no embedding, so it cannot
  be placed on the same footing. It is installed anyway because it is small and pure
  Python, keeping a graph-only side diagnostic possible.

---

## As-built exports (pre-Stage-12 freeze, 2026-08-26)

The bodies above are the **build specifications**. The **as-built resolved state** of all
five environments at the Stage 06–11b freeze is exported under `provenance/environments/`
(`conda list --explicit`, `conda env export --no-builds`, `pip freeze`, and R package
tables for `mm-qc`/`mm-annotation`), with key versions summarised in
`provenance/environments/ENVIRONMENT_SUMMARY.md`. Where a spec and an export disagree,
**the export describes what actually ran**.

Note also the writable-cache requirement the Codex audit hit: under a restricted install,
Scanpy/Numba collection fails with "no locator available" unless `NUMBA_CACHE_DIR` points
somewhere writable. That is environmental, not a test failure.
