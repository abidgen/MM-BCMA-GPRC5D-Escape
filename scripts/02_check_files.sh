#!/bin/bash
# Inspect the per-sample extracted directories before writing any Seurat loading code.
# Confirmed ground truth from this dataset (GSE223060):
#   raw/samples/<GSM_sample>/<inner_dir>/barcodes.tsv   (uncompressed, no .gz)
#   raw/samples/<GSM_sample>/<inner_dir>/counts.mtx     (NOT "matrix.mtx")
#   raw/samples/<GSM_sample>/<inner_dir>/genes.tsv      (NOT "features.tsv")
# i.e. one extra nesting level, old Cell Ranger v2-style naming, uncompressed.
# This means Read10X() won't work directly (it hardcodes matrix.mtx/features.tsv) —
# ReadMtx() with explicit paths is required. This script now searches recursively
# and matches by file role (barcodes / gene-features / matrix) rather than assuming
# any particular filename.
#
# Usage: ./02_check_files.sh [raw_dir]
# Default raw_dir: ./raw

set -euo pipefail

RAWDIR="${1:-raw}/samples"

if [ ! -d "$RAWDIR" ]; then
  echo "ERROR: $RAWDIR not found. Run 01_download_data.sh first (it populates raw/samples/)."
  exit 1
fi

N_SAMPLE_DIRS=$(find "$RAWDIR" -maxdepth 1 -mindepth 1 -type d | wc -l)
echo "== Sample directories found in $RAWDIR: $N_SAMPLE_DIRS =="
echo ""

if [ "$N_SAMPLE_DIRS" -eq 0 ]; then
  echo "ERROR: no per-sample subdirectories found. Did 01_download_data.sh's"
  echo "per-sample extraction step complete? Check for .extraction_complete markers."
  exit 1
fi

echo "== Full recursive listing of the first sample directory =="
FIRST_DIR=$(find "$RAWDIR" -maxdepth 1 -mindepth 1 -type d | sort | head -1)
echo "  $FIRST_DIR"
find "$FIRST_DIR" -type f
echo ""

echo "== Per-sample format classification (recursive, extension-based matching) =="
printf "%-40s %-12s %s\n" "SAMPLE" "FORMAT" "N_FILES"
for d in $(find "$RAWDIR" -maxdepth 1 -mindepth 1 -type d | sort); do
  sample=$(basename "$d")
  n_files=$(find "$d" -type f ! -name ".extraction_complete" | wc -l)

  # Match by extension/role, not by the literal word "matrix" — this dataset uses
  # counts.mtx / genes.tsv, which don't contain "matrix" or "features" as substrings.
  n_h5=$(find "$d" -iname "*.h5" | wc -l)
  n_mtx=$(find "$d" \( -iname "*.mtx" -o -iname "*.mtx.gz" \) | wc -l)
  n_barcodes=$(find "$d" -iname "*barcode*" | wc -l)
  n_genefeat=$(find "$d" \( -iname "*gene*" -o -iname "*feature*" \) -iname "*.tsv*" | wc -l)

  if [ "$n_h5" -gt 0 ]; then
    fmt="h5"
  elif [ "$n_mtx" -gt 0 ] && [ "$n_barcodes" -gt 0 ] && [ "$n_genefeat" -gt 0 ]; then
    fmt="triplet-ok"
  else
    fmt="INCOMPLETE"
  fi

  printf "%-40s %-12s %s\n" "$sample" "$fmt" "$n_files"
done

echo ""
echo "== Aggregate counts (extension/role-based, recursive) =="
echo "Total .h5 files:           $(find "$RAWDIR" -iname "*.h5" | wc -l)"
echo "Total .mtx / .mtx.gz files: $(find "$RAWDIR" \( -iname "*.mtx" -o -iname "*.mtx.gz" \) | wc -l)"
echo "Total barcode files:       $(find "$RAWDIR" -iname "*barcode*" | wc -l)"
echo "Total gene/feature files:  $(find "$RAWDIR" \( -iname "*gene*" -o -iname "*feature*" \) -iname "*.tsv*" | wc -l)"
echo ""
echo "Sample filenames actually used (first 3, from the first sample dir):"
find "$FIRST_DIR" -type f ! -name ".extraction_complete" -exec basename {} \;

echo ""
echo "Next: run 03_build_manifest.py raw/samples — it now recurses into the nested"
echo "sample subdirectory and records exact file paths (for Seurat's ReadMtx(),"
echo "since Read10X() won't recognize this naming convention)."
