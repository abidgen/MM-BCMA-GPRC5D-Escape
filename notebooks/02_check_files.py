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
#     display_name: Python (mm-qc)
#     language: python
#     name: mm-qc
# ---

# %% [markdown]
# # 02 — Inspect the extracted files (interactive)
#
# Interactive companion to `scripts/02_check_files.sh`. Classifies every extracted
# sample directory as `triplet-ok` / `h5` / `INCOMPLETE` **before** any loading code
# is written against it.
#
# **The script remains the single source of truth.** This notebook *invokes* it and
# then parses its output into a DataFrame for display and assertions — it does not
# re-implement the classification. Script 02 is bash, so there is no importable
# function as in notebook 03; running it and presenting its output is what "wrap,
# don't duplicate" means here. The CLI path still works:
#
# ```bash
# ./scripts/02_check_files.sh raw
# ```
#
# What this notebook adds: **assertions** on the classification, the **Cell Ranger
# reference split**, and the **checksum finding** that turns out to matter more than
# anything else in this stage.
#
# This stage is read-only. It writes nothing.

# %%
import hashlib
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

pd.set_option("display.max_rows", 80)
pd.set_option("display.width", 200)

ROOT = Path.cwd()
if ROOT.name == "notebooks":
    ROOT = ROOT.parent

SCRIPT = ROOT / "scripts" / "02_check_files.sh"
SAMPLES_DIR = ROOT / "raw" / "samples"
MANIFEST = ROOT / "raw" / "sample_manifest.csv"

print("repo root :", ROOT)
print("script    :", SCRIPT, "(exists:", SCRIPT.is_file(), ")")
print("samples   :", SAMPLES_DIR, "(exists:", SAMPLES_DIR.is_dir(), ")")

# %% [markdown]
# ## Run the script
#
# Output is captured rather than streamed — this one is fast and the interesting part
# is the table, which is parsed below. The raw text is printed first so the notebook
# shows exactly what the CLI shows.

# %%
proc = subprocess.run(
    ["bash", str(SCRIPT), "raw"],
    cwd=str(ROOT),
    capture_output=True,
    text=True,
)
print(proc.stdout)
if proc.stderr.strip():
    print("--- stderr ---", file=sys.stderr)
    print(proc.stderr, file=sys.stderr)

assert proc.returncode == 0, f"scripts/02_check_files.sh failed with exit code {proc.returncode}"

# %% [markdown]
# ## Parse the classification into a table
#
# The script prints a fixed-width `SAMPLE / FORMAT / N_FILES` block. Parsing that is
# presentation, not duplicated logic: the *decision* about what counts as
# `triplet-ok` stays in the script, and this notebook never re-derives it.

# %%
ROW = re.compile(r"^(\S+)\s+(triplet-ok|h5|INCOMPLETE)\s+(\d+)\s*$")

rows = [m.groups() for line in proc.stdout.splitlines() if (m := ROW.match(line))]
files_df = pd.DataFrame(rows, columns=["sample_id", "format", "n_files"])
files_df["n_files"] = files_df["n_files"].astype(int)

print(f"parsed {len(files_df)} sample rows\n")
print(files_df["format"].value_counts().to_string())
files_df.head(10)

# %% [markdown]
# ## Assertions
#
# Expected, per `CLAUDE.md`'s confirmed ground truth: **62 samples, all
# `triplet-ok`, zero `INCOMPLETE`.** An `INCOMPLETE` here means stage 01's extraction
# did not finish — go back to `01_download_data.ipynb` rather than continuing, since
# stage 03 would happily build a manifest with missing paths.

# %%
EXPECTED_SAMPLES = 62

incomplete = files_df[files_df["format"] == "INCOMPLETE"]
if len(incomplete):
    print(f"{len(incomplete)} INCOMPLETE sample(s):")
    display(incomplete)
else:
    print("OK: no INCOMPLETE samples.")

assert len(files_df) == EXPECTED_SAMPLES, (
    f"expected {EXPECTED_SAMPLES} samples, script classified {len(files_df)}"
)
assert incomplete.empty, f"{len(incomplete)} sample(s) missing part of the triplet"
assert (files_df["format"] == "triplet-ok").all(), (
    "expected every sample to be triplet-ok (no .h5 in this dataset)"
)
print(f"OK: {len(files_df)}/{EXPECTED_SAMPLES} triplet-ok.")

# %% [markdown]
# ## What the files are actually called
#
# The naming is the whole reason this stage exists. Old Cell Ranger v2 style,
# uncompressed, one directory level deeper than the GSM folder:
#
# ```
# raw/samples/<GSM>_<sample>/<sample>/barcodes.tsv    (no .gz)
# raw/samples/<GSM>_<sample>/<sample>/counts.mtx      (NOT matrix.mtx)
# raw/samples/<GSM>_<sample>/<sample>/genes.tsv       (NOT features.tsv)
# ```
#
# `scanpy.read_10x_mtx()` hardcodes `matrix.mtx`/`features.tsv` and will silently
# fail to find these. Stage 04 loads with `scanpy.read_mtx()` on the explicit paths
# from the stage 03 manifest.

# %%
name_counts = {}
for sample_dir in sorted(p for p in SAMPLES_DIR.iterdir() if p.is_dir()):
    for path in sample_dir.rglob("*"):
        if path.is_file() and path.name != ".extraction_complete":
            name_counts[path.name] = name_counts.get(path.name, 0) + 1

print("distinct filenames across all samples:")
for name, count in sorted(name_counts.items(), key=lambda kv: -kv[1]):
    print(f"  {count:>3} x  {name}")

# %% [markdown]
# ## Cell Ranger reference split
#
# The single most consequential property of this dataset. The 62 samples were
# processed against **two different references**, identifiable by the row count of
# `genes.tsv`:
#
# | rows | samples | reference | note |
# |---|---|---|---|
# | 33538 | 37 | GRCh38-3.0.0 (Ensembl 93) | |
# | 33694 | 24 | GRCh38-1.2.0 (Ensembl 84) | |
# | 22185 | 1  | GRCh38-1.2.0, **truncated** | `56203_1` only — repaired on read |
#
# **`56203_1` is not a third reference (corrected 2026-08-24).** Its `counts.mtx`
# header reads `33694 1837 2135520` — a normal 33694-build matrix — but its
# `genes.tsv` write stopped at row 22185, ending `KBTBD` where the reference has
# `KBTBD7`, with no trailing newline. `TNFRSF17` (canonical row 25539) was never
# absent from a reference; it was past the cut. `io.read_sample` substitutes the
# canonical column behind a prefix assertion. See `config.TRUNCATED_GENE_FILES`.
#
# Note the row count is **22185**, not the 22184 recorded previously — that figure
# was a `wc -l` artifact of the missing trailing newline.

# %%
gene_files = {}
for sample_dir in sorted(p for p in SAMPLES_DIR.iterdir() if p.is_dir()):
    hits = list(sample_dir.rglob("genes.tsv"))
    if hits:
        gene_files[sample_dir.name] = hits[0]

ref_df = pd.DataFrame(
    [
        {
            "sample_id": name,
            "n_genes_ref": sum(1 for _ in open(path, "rb")),
            "genes_md5": hashlib.md5(path.read_bytes()).hexdigest(),
        }
        for name, path in gene_files.items()
    ]
)

print(ref_df["n_genes_ref"].value_counts().sort_index().to_string())
print("\nsample(s) not on a full reference:")
# 22185, not 22184: n_genes_ref counts lines by iteration, while `wc -l` counts
# newlines and this file has none at the end. Matching on the wrong number here
# displayed an empty table rather than failing.
display(ref_df[~ref_df["n_genes_ref"].isin([33538, 33694])])
assert (ref_df["n_genes_ref"] == 22185).sum() == 1, (
    "expected exactly one truncated gene file (56203_1, 22185 rows) — see "
    "config.TRUNCATED_GENE_FILES"
)

# %% [markdown]
# ### Only three distinct gene files exist
#
# Not in the script, and the most useful thing this stage produces. Every sample on a
# given reference has a **byte-identical** `genes.tsv` — 62 files, 3 unique
# checksums.
#
# That is what makes the Ensembl-ID reconstruction in `src/mm_escape/gene_space.py`
# possible: each file is a *positional dump* of a public reference, so gene IDs can
# be restored by row-index join against the Ensembl GTF and verified
# position-for-position. Recovering them lifts the usable gene space from **22,164
# genes (symbol join) to 32,991 (ID join)**, and prevents a symbol join from pairing
# genuinely different genes that happen to share a symbol.
#
# If this ever shows more than three checksums, a new reference has appeared and
# `config.BUILDS` must be extended before the data is used.

# %%
checksums = ref_df.groupby(["genes_md5", "n_genes_ref"]).agg(
    n_samples=("sample_id", "size"),
    example=("sample_id", "min"),
).reset_index().sort_values("n_samples", ascending=False)

display(checksums)

assert len(checksums) == 3, (
    f"expected exactly 3 distinct genes.tsv files, found {len(checksums)} — "
    "a new Cell Ranger reference has appeared; extend config.BUILDS before proceeding"
)
print("OK: exactly 3 distinct gene files, as documented.")

# %% [markdown]
# ## Cross-check against the manifest
#
# Only if stage 03 has already run. Catches the case where the manifest is stale
# relative to what is actually on disk — samples added, removed, or re-extracted
# since it was last built.

# %%
if not MANIFEST.exists():
    print(f"{MANIFEST.name} not built yet — run 03_build_manifest.ipynb next.")
else:
    manifest = pd.read_csv(MANIFEST)
    on_disk = set(files_df["sample_id"])
    in_manifest = set(manifest["sample_id"])

    print(f"manifest rows: {len(manifest)}   on disk: {len(on_disk)}")
    only_disk = on_disk - in_manifest
    only_manifest = in_manifest - on_disk

    if only_disk:
        print(f"  on disk but NOT in manifest ({len(only_disk)}): {sorted(only_disk)[:5]}")
    if only_manifest:
        print(f"  in manifest but NOT on disk ({len(only_manifest)}): {sorted(only_manifest)[:5]}")
    if not only_disk and not only_manifest:
        print("  OK: manifest and disk agree exactly.")

    missing_paths = [
        p for col in ("matrix_path", "barcodes_path", "genefeat_path")
        for p in manifest[col].dropna()
        if p and not (ROOT / p).exists()
    ]
    print(f"\nmanifest paths that do not resolve: {len(missing_paths)}")
    for p in missing_paths[:5]:
        print("   ", p)

# %% [markdown]
# ## Next
#
# `notebooks/03_build_manifest.ipynb` — map each sample to the exact paths of its
# `counts.mtx` / `barcodes.tsv` / `genes.tsv` and write `raw/sample_manifest.csv`,
# the contract every later stage reads.
