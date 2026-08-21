# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#   kernelspec:
#     display_name: mm-qc
#     language: python
#     name: mm-qc
# ---

# %% [markdown]
# # 01 — Download and unpack GSE223060 / GSE223061 (interactive)
#
# Interactive companion to `scripts/01_download_data.sh`. Downloads the scRNA-seq
# archive (GSE223060, ~970 MB) and the matched bulk RNA-seq archive (GSE223061,
# ~77 MB), then extracts the 62 per-sample archives into `raw/samples/`.
#
# **The script remains the single source of truth.** This notebook *invokes* it via
# `subprocess` rather than re-implementing any of it. Script 01 is bash, so unlike
# notebook 03 — which imports `build_manifest()` as a Python function — there is
# nothing importable here; shelling out is what "wrap, don't duplicate" means for a
# shell script. The CLI path still works and produces the same result:
#
# ```bash
# ./scripts/01_download_data.sh raw
# DOWNLOAD_BULK=0 ./scripts/01_download_data.sh raw     # skip the bulk archive
# ```
#
# What this notebook adds on top of the script: a **pre-flight** that says whether
# this run will download or skip (so nobody starts a 1 GB transfer by accident), and
# **post-conditions** checked against the counts `CLAUDE.md` records as ground truth.
#
# ## Safe to re-run
#
# The script is idempotent by design and this is the normal case, not an edge case:
#
# | state | behaviour |
# |---|---|
# | archive present and `tar -tf` readable | skips the download |
# | archive present but truncated/corrupt | resumes via `wget -c` |
# | `raw/unpacked/.outer_extraction_complete` exists | skips outer extraction |
# | `raw/samples/<s>/.extraction_complete` exists | skips that sample |
#
# So on a machine that already has the data, running this notebook end to end is a
# no-op confirmation. On a fresh clone it does the full download.

# %%
import os
import shutil
import subprocess
import sys
from pathlib import Path

# Repo root, whether the kernel started in notebooks/ or the repo root
ROOT = Path.cwd()
if ROOT.name == "notebooks":
    ROOT = ROOT.parent

SCRIPT = ROOT / "scripts" / "01_download_data.sh"
RAW_DIR = ROOT / "raw"
SAMPLES_DIR = RAW_DIR / "samples"

# Set to "0" to skip the GSE223061 bulk archive. Keep it on: CLAUDE.md promotes the
# bulk data from "optional" to the stage 09 orthogonal validation set.
DOWNLOAD_BULK = "1"

print("repo root :", ROOT)
print("script    :", SCRIPT, "(exists:", SCRIPT.is_file(), ")")
print("output    :", RAW_DIR)


# %% [markdown]
# ## A streaming runner
#
# `subprocess.run(capture_output=True)` would hide everything until the download
# finished, which is useless for a ~1 GB transfer. Stream line by line instead so
# wget's progress and the per-sample `EXTRACT`/`SKIP` lines appear as they happen.

# %%
def run_bash(script: Path, *args: str, env_extra: dict | None = None) -> int:
    """Run a bash script, streaming stdout/stderr into the notebook as it goes."""
    env = {**os.environ, **(env_extra or {})}
    proc = subprocess.Popen(
        ["bash", str(script), *args],
        cwd=str(ROOT),                     # scripts assume repo-root-relative paths
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    for line in proc.stdout:
        print(line, end="")
        sys.stdout.flush()
    return proc.wait()


# %% [markdown]
# ## Pre-flight
#
# What *will* this run do — download, or skip? Checked before anything starts, since
# the answer is the difference between a few seconds and an hour.

# %%
ARCHIVES = {
    "GSE223060 scRNA-seq": RAW_DIR / "GSE223060_RAW.tar",
    "GSE223061 bulk": RAW_DIR / "GSE223061_RAW.tar",
}

for tool in ("wget", "tar", "find"):
    print(f"  {tool:<6} {'found' if shutil.which(tool) else 'MISSING — the script will abort'}")

free_gb = shutil.disk_usage(ROOT).free / 1e9
print(f"\nfree disk: {free_gb:.1f} GB  (need ~3 GB for archives + extraction)")
if free_gb < 3:
    print("  WARNING: low disk. The extraction step roughly doubles the archive footprint.")

print()
for label, path in ARCHIVES.items():
    if not path.exists():
        print(f"  {label:<22} absent            -> WILL DOWNLOAD")
        continue
    size_mb = path.stat().st_size / 1e6
    readable = subprocess.run(
        ["tar", "-tf", str(path)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    ).returncode == 0
    verdict = "will skip" if readable else "corrupt/truncated -> WILL RESUME"
    print(f"  {label:<22} {size_mb:>8.0f} MB  {verdict}")

n_existing = len(list(SAMPLES_DIR.glob("*/"))) if SAMPLES_DIR.is_dir() else 0
n_marked = len(list(SAMPLES_DIR.glob("*/.extraction_complete"))) if SAMPLES_DIR.is_dir() else 0
print(f"\nalready extracted: {n_existing} sample dirs, {n_marked} with completion markers")

# %% [markdown]
# ## Run
#
# Expect either a fast sequence of `SKIP` lines (data already present) or a long
# download followed by 62 `EXTRACT` lines.

# %%
rc = run_bash(SCRIPT, "raw", env_extra={"DOWNLOAD_BULK": DOWNLOAD_BULK})
print(f"\nexit code: {rc}")
assert rc == 0, f"scripts/01_download_data.sh failed with exit code {rc}"

# %% [markdown]
# ## Post-conditions
#
# The script reports its own counts; this checks them against what `CLAUDE.md`
# records as confirmed ground truth, so a partial extraction fails here rather than
# surfacing as a confusing error in stage 04.
#
# **Expected: 62 sample directories, each with a completion marker.** 62 rather than
# 53 because the GEO archive bundles normal/control BM beyond the paper's 53 disease
# samples — those are not junk, they are stage 07's negative control and stage 09's
# normal-plasma-cell baseline.

# %%
EXPECTED_SAMPLES = 62

sample_dirs = sorted(p for p in SAMPLES_DIR.iterdir() if p.is_dir()) if SAMPLES_DIR.is_dir() else []
marked = [p for p in sample_dirs if (p / ".extraction_complete").exists()]

print(f"sample directories : {len(sample_dirs)}  (expected {EXPECTED_SAMPLES})")
print(f"completion markers : {len(marked)}")

unmarked = [p.name for p in sample_dirs if p not in marked]
if unmarked:
    print(f"\n{len(unmarked)} directory(ies) without a completion marker — extraction did not finish:")
    for name in unmarked[:10]:
        print("   ", name)

assert len(sample_dirs) == EXPECTED_SAMPLES, (
    f"expected {EXPECTED_SAMPLES} sample directories, found {len(sample_dirs)}"
)
assert not unmarked, f"{len(unmarked)} sample(s) extracted incompletely"
print("\nOK: all 62 samples extracted with completion markers.")

# %% [markdown]
# ### The nesting level
#
# Each GSM directory contains **one extra subdirectory** named after the sample, and
# the files inside use old Cell Ranger v2 naming — `counts.mtx`, not `matrix.mtx`;
# `genes.tsv`, not `features.tsv`; uncompressed, no `.gz`.
#
# This is why `scanpy.read_10x_mtx()` cannot read this dataset: it hardcodes the
# filenames it expects and will simply not find these. Stage 04's loader uses
# `scanpy.read_mtx()` on the explicit paths the stage 03 manifest records.

# %%
first = sample_dirs[0]
print(first.name)
for path in sorted(first.rglob("*")):
    if path.name == ".extraction_complete":
        continue
    rel = path.relative_to(first)
    kind = "dir " if path.is_dir() else f"{path.stat().st_size / 1e6:>7.1f} MB"
    print(f"  {kind}  {rel}")

# %% [markdown]
# ## Matched bulk RNA-seq (GSE223061)
#
# Not optional any more — this is the **stage 09 orthogonal validation set**, and the
# only independent check available on the antigen quantification. Two gotchas,
# confirmed by direct inspection and re-checked here:
#
# - **Two files are empty 114-byte gzip stubs** and must be excluded:
#   `GSM6939104_MMRF_1505_tpm.tsv.gz`, `GSM6939120_MMRF_2259_tpm.tsv.gz`. A 114-byte
#   read is a **failed deposit**, not zero expression — treating it as the latter
#   would fabricate a data point.
# - **Three bulk/sc sample-ID mismatches** to reconcile against Supplementary Table
#   S1 (do not guess these pairings): bulk `47499` vs sc `47491_1`/`47491_2`; bulk
#   `98433` vs sc `MMY98423`; bulk `59114_2` vs sc `59114_1`/`59114_4`.
#
# **Inventory corrected 2026-08-21.** `CLAUDE.md` recorded "18 MMRF + 12 WashU = 30
# usable". The archive actually holds **13** WashU `.tar.gz` (GSM6939090-102, all
# 4.5-5.4 MB) and 18 MMRF TPM tables (GSM6939103-120), so after excluding the two
# stubs the usable count is **29**, not 30. The assertion below pins the corrected
# numbers so a future re-download that drops a file is caught here.

# %%
BULK_DIR = RAW_DIR / "unpacked_bulk"
KNOWN_STUBS = {"GSM6939104_MMRF_1505_tpm.tsv.gz", "GSM6939120_MMRF_2259_tpm.tsv.gz"}

EXPECTED_TPM, EXPECTED_WASHU = 18, 13

if not BULK_DIR.is_dir():
    print(f"{BULK_DIR} absent — bulk archive not downloaded (DOWNLOAD_BULK={DOWNLOAD_BULK}).")
    print("Stage 09 needs it; re-run with DOWNLOAD_BULK='1'.")
else:
    # Exclude dotfiles: the script drops a zero-byte `.outer_extraction_complete`
    # marker here, which is bookkeeping, not a failed deposit.
    files = sorted(p for p in BULK_DIR.iterdir() if p.is_file() and not p.name.startswith("."))
    tpm = [p for p in files if p.name.endswith("_tpm.tsv.gz")]
    tars = [p for p in files if p.name.endswith(".tar.gz")]
    stubs = [p for p in files if p.stat().st_size <= 200]

    print(f"bulk sample files: {len(files)}")
    print(f"  MMRF TPM tables: {len(tpm):>3}  (expected {EXPECTED_TPM})")
    print(f"  WashU archives : {len(tars):>3}  (expected {EXPECTED_WASHU})")
    print(f"\nempty/stub files (<=200 bytes): {len(stubs)}")
    for p in stubs:
        flag = "known stub" if p.name in KNOWN_STUBS else "NEW — investigate"
        print(f"   {p.stat().st_size:>4} B  {p.name}   [{flag}]")

    usable = len(files) - len(stubs)
    print(f"\nusable bulk samples: {usable}  = {len(tpm)} TPM - {len(stubs)} stubs + {len(tars)} WashU")

    unexpected = {p.name for p in stubs} - KNOWN_STUBS
    if unexpected:
        print(f"\nWARNING: undocumented empty file(s): {sorted(unexpected)}")
        print("Add them to the stage 09 exclusion list before correlating anything.")
    else:
        print("OK: exactly the two documented stubs, no new ones.")

    assert len(tpm) == EXPECTED_TPM, f"expected {EXPECTED_TPM} MMRF TPM tables, found {len(tpm)}"
    assert len(tars) == EXPECTED_WASHU, f"expected {EXPECTED_WASHU} WashU archives, found {len(tars)}"
    assert not unexpected, f"undocumented empty deposit(s): {sorted(unexpected)}"

# %% [markdown]
# ## Next
#
# `notebooks/02_check_files.ipynb` — classify every sample's file triplet and confirm
# **62/62 `triplet-ok`** before the manifest is built.
#
# Nothing here writes into `results/`; stages 01-03 produce `raw/` and
# `raw/sample_manifest.csv` only.
