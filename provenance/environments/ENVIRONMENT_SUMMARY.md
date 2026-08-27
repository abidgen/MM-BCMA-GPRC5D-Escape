# Environment provenance — pre-Stage-12 freeze checkpoint

Captured 2026-08-26 during the pre-Stage-12 reproducibility/provenance repair pass.
**No environment was created, modified, upgraded or rebuilt to produce this record.**
These are exports of the environments as they stood when stages 04–11b were run.

Machine: Linux 7.0.0-30-generic · GPU NVIDIA GeForce RTX 5070 (driver 595.84).
The GPU is relevant only to stage 05b (`scVI` arm); every frozen Stage 07–10 result is CPU work.

## Files per environment

| file | produced by |
|---|---|
| `<env>.conda-explicit.txt` | `conda list -n <env> --explicit` |
| `<env>.conda-env-no-builds.yml` | `conda env export -n <env> --no-builds` |
| `<env>.pip-freeze.txt` | `conda run -n <env> python -m pip freeze` |
| `<env>.R-packages.tsv` | `Rscript -e 'installed.packages()[,c("Package","Version")]'` (R envs only) |

The `envs/*.yml` files in the repo root are the **build specifications**; the files here are
the **as-built resolved state**. When they disagree, the files here describe what actually ran.

## Key versions

| env | stages | python | scanpy | anndata | numpy | pandas | scipy |
|---|---|---|---|---|---|---|---|
| `mm-qc` | 04 | 3.12.13 | 1.11.5 | 0.13.2 | 2.5.2 | 3.0.5 | 1.18.0 |
| `mm-core` | 05, 07–10, 12, tests | 3.12.13 | 1.11.5 | 0.13.2 | 2.5.2 | 3.0.5 | 1.18.0 |
| `mm-annotation` | 06 | 3.12.13 | 1.11.5 | 0.13.2 | 2.5.2 | 3.0.5 | 1.18.0 |
| `mm-integration` | 05b | 3.12.14 | 1.11.5 | 0.13.2 | 2.5.2 | 3.0.5 | 1.18.0 |
| `mm-communication` | 11, 11b | 3.12.13 | **1.12.3** | **0.12.19** | **2.0.2** | **2.3.3** | 1.18.0 |

`mm-communication` deliberately runs a different scanpy/anndata generation from the others.
It meets `mm-core` only on disk, through `.h5ad`/CSV files.

### Analysis-critical packages

| package | env | version | used by |
|---|---|---|---|
| `scDblFinder` | `mm-qc` (R 4.3.3) | 1.16.0 | stage 04 doublet calls |
| `SingleCellExperiment` | `mm-qc` (R 4.3.3) | 1.24.0 | stage 04 R bridge |
| `rpy2` | `mm-qc`, `mm-annotation` | 3.5.11 | R bridge |
| `SingleR` | `mm-annotation` (R 4.3.3) | 2.4.0 | stage 06 |
| `celldex` | `mm-annotation` (R 4.3.3) | 1.12.0 | stage 06 references |
| `celltypist` | `mm-core`/`mm-annotation` | 1.7.1 | stage 06 |
| `harmonypy` | `mm-core` | 2.0.0 | stage 05 integration |
| `infercnvpy` | `mm-core` | 0.6.1 | stage 07 CNV (rejected as NOT_EVALUABLE) |
| `pydeseq2` | `mm-core` | 0.5.4 | stage 10 pseudobulk DE |
| `decoupler` | `mm-core` | **2.2.0** | stage 10 Hallmark/PROGENy/CollecTRI |
| `scib-metrics` | `mm-integration` | 0.6.0 | stage 05b benchmark |
| `scvi-tools` | `mm-integration` | 1.5.0 | stage 05b scVI arms |
| `liana` | `mm-communication` | **1.8.1** | stage 11b verification arm |
| `decoupler` | `mm-communication` | 2.1.6 | (not used for stage-10 results) |
| `statsmodels` | all | 0.14.6 | OLS / BH correction |

`decoupler` 2.x is an API rewrite (`dc.mt.*` / `dc.op.*`). Stage 10 was written against 2.2.0.
Note the two environments carry **different** decoupler versions; the frozen Stage-10
pathway results come from `mm-core`'s 2.2.0, never `mm-communication`'s 2.1.6.

## Known environment trap (recorded by the Codex audit)

Under a read-only or restricted install, Scanpy/Numba collection fails with
"no locator available" because Numba cannot write its default cache. Set a writable cache:

```bash
NUMBA_CACHE_DIR=/tmp/mm_numba_cache conda run -n mm-core pytest -q
```

This is environmental, not a test failure.
