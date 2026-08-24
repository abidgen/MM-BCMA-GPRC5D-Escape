# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.5
#   kernelspec:
#     display_name: mm-qc
#     language: python
#     name: mm-qc
# ---

# %% [markdown]
# # 03 — Build the sample manifest (interactive)
#
# Interactive companion to `scripts/03_build_manifest.py`. Maps each of the 62 GEO
# sample directories to the exact paths of its `counts.mtx` / `barcodes.tsv` /
# `genes.tsv`, and writes `raw/sample_manifest.csv`.
#
# **The script remains the single source of truth.** This notebook *imports*
# `build_manifest()` from it rather than re-implementing it, so the two cannot drift.
# `CLAUDE.md` says stages 01-03 are reused as-is and not notebook-ified; this notebook
# is an additional interactive view, not a replacement. The CLI path still works:
#
# ```bash
# python scripts/03_build_manifest.py raw/samples
# ```
#
# Beyond what the script prints, this notebook adds three checks that are useful to
# eyeball before stage 04: the **Cell Ranger reference split**, the **required-gene
# presence assertions**, and a **diff against the existing manifest**.

# %%
import sys
from pathlib import Path
import importlib.util
import pathlib

import pandas as pd

pd.set_option("display.max_rows", 80)
pd.set_option("display.width", 200)

# Repo root, whether the kernel started in notebooks/ or the repo root
ROOT = Path.cwd()
if ROOT.name == "notebooks":
    ROOT = ROOT.parent

RAW_DIR = ROOT / "raw" / "samples"
OUT_CSV = ROOT / "raw" / "sample_manifest.csv"

print("repo root :", ROOT)
print("raw dir   :", RAW_DIR, "(exists:", RAW_DIR.is_dir(), ")")
print("output    :", OUT_CSV)

# %% [markdown]
# ## Import the builder from the script
#
# Loaded by file path so the notebook does not depend on `scripts/` being importable.
#
# **Do not use the script's module-level `RAW_DIR` / `OUT_CSV`.** They are derived from
# `sys.argv`, and under a Jupyter kernel `sys.argv[1]` is `-f` (the connection-file
# flag) — so `mf.RAW_DIR` evaluates to `Path('-f')`. Verified, not hypothetical. Always
# pass paths explicitly, as below.

# %%
_spec = importlib.util.spec_from_file_location(
    "build_manifest_mod", ROOT / "scripts" / "03_build_manifest.py"
)
mf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mf)

print("imported:", mf.__file__)
print("script's module-level RAW_DIR is", repr(str(mf.RAW_DIR)), "-> ignored on purpose")

# %% [markdown]
# ## Build

# %%
df = mf.build_manifest(RAW_DIR)

# `build_manifest` records whatever path it was handed. Passing an absolute RAW_DIR
# (which the notebook must, since the kernel's cwd is not guaranteed) therefore yields
# absolute paths — but the manifest is committed to git and read on other machines, so
# it must hold paths relative to the repo root, exactly as the CLI produces when run
# from there. Normalize rather than leaving the notebook and script to disagree.
PATH_COLS = ["h5_path", "matrix_path", "barcodes_path", "genefeat_path"]
for col in PATH_COLS:
    df[col] = df[col].map(
        lambda v: str(pathlib.PurePath(v).relative_to(ROOT)) if v else ""
    )

print(f"{len(df)} samples")
df.head(10)

# %% [markdown]
# ## Format breakdown
#
# Expect **62 `triplet-ok`**, zero `INCOMPLETE`. Anything else means the extraction in
# stage 01 is incomplete — go back to `scripts/01_download_data.sh` before continuing.

# %%
print(df["format"].value_counts().to_string())

incomplete = df[df["format"] == "INCOMPLETE"]
if len(incomplete):
    print(f"\n{len(incomplete)} INCOMPLETE sample(s):")
    display(incomplete[["sample_id", "n_mtx_found", "n_barcode_found", "n_genefeat_found"]])
else:
    print("\nOK: every sample has a complete triplet.")

# %% [markdown]
# ## Ambiguous matches
#
# More than one candidate file for a role means the filename matching was too loose for
# that sample; the manifest silently takes the first match, so these need a human look.

# %%
ambiguous = df[(df["n_mtx_found"] > 1) | (df["n_barcode_found"] > 1) | (df["n_genefeat_found"] > 1)]
if len(ambiguous):
    print(f"{len(ambiguous)} ambiguous sample(s) — manifest took the first match:")
    display(ambiguous[["sample_id", "n_mtx_found", "n_barcode_found", "n_genefeat_found"]])
else:
    print("OK: exactly one file matched per role for all samples.")

# %% [markdown]
# ## Cell Ranger reference split
#
# Not in the script — but this is the single most consequential property of this
# dataset, so it is worth seeing every time the manifest is rebuilt.
#
# The 62 samples were processed against **three different references**, identifiable by
# the row count of `genes.tsv`. Expected (confirmed ground truth, `CLAUDE.md`):
#
# | genes | samples | note |
# |---|---|---|
# | 33538 | 37 | |
# | 33694 | 24 | |
# | 22185 | 1  | `56203_1` only — a **truncated** 33694 file, repaired on read |
#
# `56203_1`'s `genes.tsv` write failed at row 22185 (ends `KBTBD`, where the 33694
# reference has `KBTBD7`); its matrix declares the full 33694 rows. `TNFRSF17` was
# past the cut, not absent. Repaired on read — see `config.TRUNCATED_GENE_FILES`.

# %%
def count_lines(path: str) -> int:
    # manifest paths are repo-root-relative; the kernel's cwd may be notebooks/
    with open(ROOT / path, "rb") as fh:
        return sum(1 for _ in fh)

df["n_genes_ref"] = df["genefeat_path"].map(count_lines)

print(df["n_genes_ref"].value_counts().sort_index().to_string())
print("\nSamples on a minority reference:")
display(
    df[df["n_genes_ref"] == df["n_genes_ref"].value_counts().idxmin()][
        ["sample_id", "n_genes_ref"]
    ]
)

# %% [markdown]
# ### Gene-space intersection preview — with symbol harmonization
#
# Stage 05 intersects gene sets across retained samples (never unions — a union would
# make ~11k genes structurally zero in whole cohorts, indistinguishable from the
# biological zeros this project measures).
#
# **But a naive set intersection is wrong here, and this notebook found out why.** The
# two references use *different HGNC symbol vintages*. Some genes are present in both
# builds under different names, so intersecting raw symbols drops them entirely:
#
# | 33538 (newer) | 33694 (older) | consequence if unharmonized |
# |---|---|---|
# | `NSD2`   | `WHSC1`   | **t(4;14) becomes uncallable** — the highest-risk MM translocation |
# | `TENT5C` | `FAM46C`  | loses a recurrently-deleted MM tumour suppressor |
# | `NSD3`   | `WHSC1L1` | (distinct gene from NSD2 — do not conflate) |
# | `ATP5F1A`| `ATP5A1`  | OXPHOS program member |
#
# So: **canonicalize symbols first, then intersect.** The map below covers the genes
# this project depends on. A full HGNC reconciliation belongs in
# `src/mm_escape/gene_space.py`, not here — this is the targeted version needed to
# validate the manifest.

# %%
# legacy symbol (older 33694 build) -> current symbol (newer 33538 build)
LEGACY_TO_CURRENT = {
    "WHSC1":   "NSD2",
    "WHSC1L1": "NSD3",
    "FAM46C":  "TENT5C",
    "ATP5A1":  "ATP5F1A",
}

def canonicalize(symbols: set) -> set:
    return {LEGACY_TO_CURRENT.get(s, s) for s in symbols}

retained = df
print(f"retained samples: {len(retained)} of {len(df)}  (56203_1 repaired, not excluded)")

raw_sets, canon_sets = {}, {}
for ref, sub in retained.groupby("n_genes_ref"):
    p = ROOT / sub.iloc[0]["genefeat_path"]
    raw = {ln.strip() for ln in open(p) if ln.strip()}
    raw_sets[ref], canon_sets[ref] = raw, canonicalize(raw)
    print(f"  ref {ref}: {len(raw):,} symbols   (example sample {sub.iloc[0]['sample_id']})")

naive = set.intersection(*raw_sets.values())
intersection = set.intersection(*canon_sets.values())

print(f"\nnaive intersection      : {len(naive):,} genes")
print(f"harmonized intersection : {len(intersection):,} genes  (+{len(intersection) - len(naive)})")
recovered = sorted(intersection - naive)
print(f"recovered by harmonization: {recovered}")

# %% [markdown]
# Note the ~11k symbols unique to each build are dominated by **annotation-version
# noise** — versioned clone identifiers such as `AC000032.1` vs `AC000032.2` — rather
# than genuinely absent genes. The 22,164-gene naive intersection therefore *understates*
# the recoverable gene space. That is tolerable (HVG selection never reaches those
# lncRNA/clone entries), but it is not tolerable for named genes the pipeline requires,
# which is what the harmonization above and the assertions below exist to guarantee.

# %% [markdown]
# ## Write the manifest
#
# Diffed against the existing file first, so a rebuild that silently changes paths or
# drops samples is visible rather than overwritten in place.

# %%
if OUT_CSV.exists():
    old = pd.read_csv(OUT_CSV)
    shared = [c for c in old.columns if c in df.columns and c != "n_genes_ref"]
    same = (
        len(old) == len(df)
        and old["sample_id"].tolist() == df["sample_id"].tolist()
        and old[shared].fillna("").astype(str).equals(df[shared].fillna("").astype(str))
    )
    print("identical to existing manifest" if same
          else f"DIFFERS from existing manifest (old {len(old)} rows / new {len(df)} rows) — inspect before writing")
else:
    print("no existing manifest; this will create it")

# %%
# `n_genes_ref` is computed above as a diagnostic and is deliberately NOT written:
# the manifest's schema must stay byte-identical to what
# `python scripts/03_build_manifest.py raw/samples` produces, or the notebook and the
# CLI would silently disagree about the pipeline's own input contract. Persisting it
# would mean adding the column to the script too — a change `CLAUDE.md` rules out for
# stages 01-03. Stage 05 recomputes it cheaply from `genefeat_path`.
SCRIPT_SCHEMA = [c for c in df.columns if c != "n_genes_ref"]
df[SCRIPT_SCHEMA].to_csv(OUT_CSV, index=False)

print(f"wrote {len(df)} rows x {len(SCRIPT_SCHEMA)} cols -> {OUT_CSV}")
print("(n_genes_ref shown above, not persisted — see comment)")
df[["sample_id", "format", "n_genes_ref", "matrix_path"]].head(10)

# %% [markdown]
# ## Next
#
# `notebooks/04_qc.ipynb` (env `mm-qc`) — per-sample load via `src/mm_escape/io.py`,
# QC metrics, MAD-based outlier filtering (per cohort), `scDblFinder`,
# per-sample checkpoints to `results/04_qc/`.
