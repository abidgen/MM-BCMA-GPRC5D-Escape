# `provenance/` — the freeze record for stages 04–11b

**Created 2026-08-26 during the pre-Stage-12 reproducibility/provenance repair pass.**
Nothing in `results/` was recomputed, altered or moved to build this. No scientific
result, threshold, tier, patient state or interpretation changed.

## What is here

| file | what it is |
|---|---|
| `frozen_artifacts_pre_stage12.tsv` | the freeze manifest — one row per frozen artifact (393 rows, stages 04–11b) |
| `frozen_artifacts_pre_stage12.sha256` | the same hashes in `sha256sum -c` format |
| `environments/` | as-built environment exports for all five conda environments, plus `ENVIRONMENT_SUMMARY.md` |
| `README.md` | this file |

## Why the large artifacts are not in Git

`results/` holds ~13 GB, including three `.h5ad` checkpoints of 0.8–1.3 GB each and 62
per-sample QC checkpoints. Commit `4a2b809` removed generated results from version control
and `.gitignore` excludes `results/` wholesale. That decision stands: ordinary Git is the
wrong store for multi-gigabyte binary matrices that change wholesale on every rerun, and
re-adding them would make the repository unusable to clone.

What the audit correctly objected to was not the exclusion — it was that **nothing
committed authenticated the excluded state**. Digest files existed (`frozen_upstream_digests.json`,
`frozen_state_digests.json`) but both the digests and their targets lived under ignored
`results/`, so a later local rewrite could replace artifact and digest together and Git
would see nothing. That is what this directory fixes.

## How the manifest authenticates the frozen state

Every row carries:

| column | meaning |
|---|---|
| `stage` | `04`, `05`, `05b`, `06`, `07`, `08`, `08c`, `09`, `09b`, `10`, `11`, `11b` |
| `path` | repo-relative path under `results/` |
| `bytes` | size at freeze time |
| `sha256` | content hash at freeze time |
| `artifact_role` | what the file is (result table, checkpoint, frozen design document, figure, …) |
| `producer_source_path` | the **committed** code that produced it (`production/…`, `notebooks/…`) |
| `producer_code_commit` | the commit that carries that producer, or `pending-commit(this-repair-pass)` for newly recovered drivers |
| `frozen_date` | the date that stage was accepted/frozen |
| `environment` | which conda environment ran it |
| `reproducibility` | see below |

### The four reproducibility classes

| class | count | meaning |
|---|---|---|
| `REPRODUCIBLE_FROM_COMMITTED_CODE` | 191 | producer is committed and every input is either a committed resource or another manifested frozen artifact |
| `REPRODUCIBLE_WITH_EXTERNAL_INPUT` | 180 | producer is committed but a rerun additionally needs something outside the repo — the ignored `raw/` GEO deposit (stages 04–05b), `celldex`/`SingleR` reference downloads (06), or LIANA's bundled `consensus` resource (11b) |
| `ARCHIVED_ONLY` | 22 | frozen design/summary documents authored in place; there is no analysis code to rerun, and the document *is* the artifact |
| `NOT_FULLY_REPRODUCIBLE` | 0 | reserved; nothing currently falls here |

`REPRODUCIBLE_*` is a statement about the **derivation being reviewable and re-runnable in
principle**, not a promise of bit-identical output. Harmony, Leiden, scVI and R packages
can move with library and hardware versions; that is why `environments/` exists and why
project-level seeds are fixed.

## Verifying the freeze locally

```bash
cd /media/wrath/CART_mm_dual_antigen
sha256sum -c provenance/frozen_artifacts_pre_stage12.sha256
```

or, as part of the suite (this is the slow path — it reads ~6 GB):

```bash
NUMBA_CACHE_DIR=/tmp/mm_numba_cache conda run -n mm-core \
  pytest tests/test_production_paths.py -q -m slow
```

**If a hash does not match, stop.** Do not regenerate the manifest to agree with the file.
A changed hash means a frozen artifact mutated, which is exactly the event this record
exists to detect. Report it before doing anything else.

## Which artifacts Stage 12 requires

Stage 12 consumes five evidence layers plus the bias table. The manifest rows below are
its required inputs; all are present and hashed:

| layer | key artifacts |
|---|---|
| measurement (08 / 09b) | `08_dual_antigen_escape/patient_antigen_states_{primary,sensitivity}.csv`, `patient_bootstrap_intervals.csv`, `patient_conegativity_enrichment.csv`, `patient_evidence_states.csv`, `risk_tier_provisional/risk_tiers_provisional.csv` |
| DN coherence (10) | `10_dn_coherence/stage10_evidence_levels.csv`, `dn_coherence_final_states.csv`, `level2_program_cohort_tests.csv`, `pseudobulk_de_results.csv`, `gamma_secretase_hypothesis.csv`, `tc_like_subtype.csv` |
| multi-antigen coverage (08c) | `08_dual_antigen_escape/multi_antigen_coverage/stage12_multi_antigen_interface.csv`, `target_measurement_qc.csv` |
| bulk context (09) | `09_bulk_validation/bulk_vs_sc_by_cohort.csv`, `normal_marrow_antigen_context.csv` |
| immune context (11 / 11b) | `11_immune_context/immune_vs_dn_measurement.csv`, `liana_verification/liana_vs_dn_associations.csv` |
| bias direction | `04_qc/umi_censoring_effect.csv`, `08_dual_antigen_escape/noise_floor_{ambient,technical_zero}.csv`, `truncate10k_sensitivity.csv` |

## Regenerated vs archived

- **Regenerated** (in principle, from committed code + the inputs named in the class above):
  everything marked `REPRODUCIBLE_*`. The producer is named per row; the recovered drivers
  live in `production/` and their run order is `production/README.md`.
- **Archived only**: the 22 frozen design and summary documents. These were authored in
  place — `stage08_predeclaration.md`, `stage10_design.md`, `multi_antigen_design.md`,
  `CNV_FROZEN_NOT_EVALUABLE.md`, `README_PROVISIONAL.md` and similar. Several are
  *predeclarations written before the corresponding results existed*, which is precisely
  why they must be preserved verbatim rather than regenerated.

## Where immutable external storage should live

Not yet provisioned, and deliberately not invented here. When it is, the requirements are:

1. **Write-once or versioned** object storage (S3/GCS with object versioning and a
   retention lock, an institutional archive, or a Zenodo deposit for the publishable subset).
2. Upload the **whole** `results/` tree for the manifested stages, preserving relative paths,
   so `frozen_artifacts_pre_stage12.sha256` verifies unchanged against the download.
3. Record the bucket/DOI and upload date **in this file**, plus the manifest's own SHA256,
   so the pointer is itself version-controlled.
4. The `.h5ad` checkpoints (`05_integration/integrated.h5ad`, `06_annotation/annotated.h5ad`,
   `08_dual_antigen_escape/antigen_states.h5ad`, `04_qc/samples/*.h5ad`) dominate the
   footprint at ~5 GB; everything else is a few hundred MB of CSV and can be archived
   independently if cost matters.

Until that exists, the honest statement is: **the frozen artifacts are local, and the
committed manifest is what makes their present state checkable.**

## What this directory is not

It is not a licence to commit results. `results/` stays ignored. The rule is:

> committed code + committed manifests + committed environment metadata,
> with large artifacts external/ignored and authenticated by hash.
