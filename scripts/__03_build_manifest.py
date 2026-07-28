#!/usr/bin/env python3
"""
Build a sample -> file-path manifest from the unpacked GEO archive, handling
either archive layout (per-sample subdirectories, or flat GSM-prefixed files).

Run this AFTER 02_check_files.sh has shown you which pattern you have — this
script auto-detects it, but it's worth confirming visually first rather than
trusting auto-detection blindly on a dataset you haven't inspected yet.

Usage:
    python 03_build_manifest.py [raw_dir]

Default raw_dir: raw/unpacked
Output: raw/sample_manifest.csv
"""

import re
import sys
from pathlib import Path
import pandas as pd

RAW_DIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("raw/unpacked")
OUT_CSV = RAW_DIR.parent / "sample_manifest.csv"

GSM_PATTERN = re.compile(r"^(GSM\d+)_?(.*)")


def detect_pattern(raw_dir: Path) -> str:
    subdirs = [p for p in raw_dir.iterdir() if p.is_dir()]
    files = [p for p in raw_dir.iterdir() if p.is_file()]
    if subdirs and not files:
        return "subdirs"
    if files and not subdirs:
        return "flat"
    if subdirs and files:
        return "mixed"
    raise RuntimeError(f"No files or subdirectories found under {raw_dir}")


def build_manifest_subdirs(raw_dir: Path) -> pd.DataFrame:
    rows = []
    for sample_dir in sorted(p for p in raw_dir.iterdir() if p.is_dir()):
        h5_files = list(sample_dir.glob("*.h5"))
        rows.append({
            "sample_id": sample_dir.name,
            "path": str(sample_dir),
            "format": "h5" if h5_files else "triplet",
            "h5_file": str(h5_files[0]) if h5_files else "",
        })
    return pd.DataFrame(rows)


def build_manifest_flat(raw_dir: Path) -> pd.DataFrame:
    # Group files by GSM accession prefix
    groups: dict[str, list[Path]] = {}
    for f in sorted(raw_dir.iterdir()):
        if not f.is_file():
            continue
        m = GSM_PATTERN.match(f.name)
        if not m:
            print(f"WARNING: filename doesn't match expected GSM prefix pattern: {f.name}")
            continue
        gsm, rest = m.groups()
        groups.setdefault(gsm, []).append(f)

    rows = []
    for gsm, files in groups.items():
        h5_files = [f for f in files if f.suffix == ".h5"]
        has_barcodes = any("barcode" in f.name.lower() for f in files)
        has_features = any(("feature" in f.name.lower() or "gene" in f.name.lower()) for f in files)
        has_matrix = any("matrix" in f.name.lower() for f in files)

        # Try to recover a human-readable sample name from the filename remainder
        # after the GSM prefix (varies by submission — inspect manually if this
        # looks wrong for your specific series)
        sample_name_guess = files[0].name[len(gsm) + 1:]
        sample_name_guess = re.sub(
            r"_(barcodes|features|genes|matrix)\.(tsv|mtx)(\.gz)?$", "",
            sample_name_guess, flags=re.IGNORECASE
        )

        fmt = "h5" if h5_files else (
            "triplet" if (has_barcodes and has_features and has_matrix) else "incomplete"
        )

        rows.append({
            "sample_id": gsm,
            "sample_name_guess": sample_name_guess,
            "format": fmt,
            "n_files": len(files),
            "h5_file": str(h5_files[0]) if h5_files else "",
            "all_files": ";".join(str(f) for f in files),
        })
    return pd.DataFrame(rows)


def main():
    pattern = detect_pattern(RAW_DIR)
    print(f"Detected archive layout: {pattern}")

    if pattern == "subdirs":
        df = build_manifest_subdirs(RAW_DIR)
    elif pattern == "flat":
        df = build_manifest_flat(RAW_DIR)
    else:
        print("Mixed layout detected — both subdirectories and loose files present.")
        print("This needs manual inspection; auto-building a manifest here would guess wrong.")
        sys.exit(1)

    df.to_csv(OUT_CSV, index=False)
    print(f"\nWrote manifest for {len(df)} samples to {OUT_CSV}")
    print("\nFormat breakdown:")
    print(df["format"].value_counts().to_string())

    incomplete = df[df["format"] == "incomplete"]
    if len(incomplete) > 0:
        print(f"\nWARNING: {len(incomplete)} sample(s) have an incomplete triplet "
              f"(missing barcodes/features/matrix) and no .h5 fallback. Inspect these "
              f"manually before running the Seurat loading script:")
        print(incomplete["sample_id"].to_string(index=False))

    print("\nPreview:")
    print(df.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
