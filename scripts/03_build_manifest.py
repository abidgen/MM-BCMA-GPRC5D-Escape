#!/usr/bin/env python3
"""
Build a sample -> exact-file-path manifest from the extracted GEO archive.

Confirmed ground truth for GSE223060 (from direct inspection):
    raw/samples/<GSM_sample>/<inner_dir>/barcodes.tsv   (uncompressed, no .gz)
    raw/samples/<GSM_sample>/<inner_dir>/counts.mtx     (NOT "matrix.mtx")
    raw/samples/<GSM_sample>/<inner_dir>/genes.tsv      (NOT "features.tsv")

This is old Cell Ranger v2-style naming, uncompressed, with an extra nesting level
(each GSM folder contains one inner subdirectory named after the sample, e.g.
GSM6939028_MMRF_1695/MMRF_1695/). Because the matrix file isn't named matrix.mtx,
Seurat's Read10X() (which hardcodes expected filenames) will NOT find it even when
pointed at the right directory. Use ReadMtx() with the explicit paths this script
records instead.

This script recurses into each top-level GSM sample directory (regardless of nesting
depth) and matches files by role — barcodes / gene-features / matrix — using
extension and filename-substring checks broad enough to catch either the standard
(barcodes.tsv[.gz], features.tsv[.gz]/genes.tsv[.gz], matrix.mtx[.gz]) or .h5 formats.

Usage:
    python 03_build_manifest.py [raw_dir]

Default raw_dir: raw/samples
Output: raw/sample_manifest.csv
"""

import sys
from pathlib import Path
import pandas as pd

RAW_DIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("raw/samples")
OUT_CSV = RAW_DIR.parent / "sample_manifest.csv"

MTX_EXTS = (".mtx", ".mtx.gz")
BARCODE_MARKERS = ("barcode",)
GENEFEAT_MARKERS = ("gene", "feature")
TSV_EXTS = (".tsv", ".tsv.gz")


def is_mtx(f: Path) -> bool:
    return f.name.lower().endswith(MTX_EXTS)


def is_barcodes(f: Path) -> bool:
    name = f.name.lower()
    return any(m in name for m in BARCODE_MARKERS) and name.endswith(TSV_EXTS)


def is_genefeat(f: Path) -> bool:
    name = f.name.lower()
    return any(m in name for m in GENEFEAT_MARKERS) and name.endswith(TSV_EXTS)


def build_manifest(raw_dir: Path) -> pd.DataFrame:
    sample_dirs = sorted(p for p in raw_dir.iterdir() if p.is_dir())
    if not sample_dirs:
        raise RuntimeError(f"No sample subdirectories found under {raw_dir}")

    rows = []
    for sample_dir in sample_dirs:
        gsm_sample_id = sample_dir.name
        # Recurse fully — depth of the inner sample subdirectory can vary
        all_files = [f for f in sample_dir.rglob("*")
                     if f.is_file() and f.name != ".extraction_complete"]

        h5_files = [f for f in all_files if f.name.lower().endswith(".h5")]
        mtx_files = [f for f in all_files if is_mtx(f)]
        barcode_files = [f for f in all_files if is_barcodes(f)]
        genefeat_files = [f for f in all_files if is_genefeat(f)]

        if h5_files:
            fmt = "h5"
        elif mtx_files and barcode_files and genefeat_files:
            fmt = "triplet-ok"
        else:
            fmt = "INCOMPLETE"

        rows.append({
            "sample_id": gsm_sample_id,
            "format": fmt,
            "n_files_total": len(all_files),
            "h5_path": str(h5_files[0]) if h5_files else "",
            "matrix_path": str(mtx_files[0]) if mtx_files else "",
            "barcodes_path": str(barcode_files[0]) if barcode_files else "",
            "genefeat_path": str(genefeat_files[0]) if genefeat_files else "",
            "n_mtx_found": len(mtx_files),
            "n_barcode_found": len(barcode_files),
            "n_genefeat_found": len(genefeat_files),
        })

    return pd.DataFrame(rows)


def main():
    df = build_manifest(RAW_DIR)

    df.to_csv(OUT_CSV, index=False)
    print(f"Wrote manifest for {len(df)} samples to {OUT_CSV}")
    print("\nFormat breakdown:")
    print(df["format"].value_counts().to_string())

    # Flag anything with more than one candidate match per role — means the
    # extension/substring matching was too loose for that sample and needs a look
    ambiguous = df[(df["n_mtx_found"] > 1) | (df["n_barcode_found"] > 1) | (df["n_genefeat_found"] > 1)]
    if len(ambiguous) > 0:
        print(f"\nWARNING: {len(ambiguous)} sample(s) matched more than one candidate "
              f"file for a role — inspect manually, the manifest picked the first match:")
        print(ambiguous[["sample_id", "n_mtx_found", "n_barcode_found", "n_genefeat_found"]]
              .to_string(index=False))

    incomplete = df[df["format"] == "INCOMPLETE"]
    if len(incomplete) > 0:
        print(f"\nWARNING: {len(incomplete)} sample(s) are missing one or more required "
              f"files. These will fail to load — inspect before running the QC script:")
        print(incomplete["sample_id"].to_string(index=False))
    else:
        print("\nAll samples have a complete triplet or .h5 file. Good to proceed.")

    print("\nPreview:")
    print(df[["sample_id", "format", "matrix_path"]].head(10).to_string(index=False))


if __name__ == "__main__":
    main()
