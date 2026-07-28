#!/usr/bin/env bash
# Download and unpack GSE223060 (scRNA-seq) and, optionally, GSE223061
# (matched bulk RNA-seq). The GSE223060 outer archive contains one .tar.gz
# archive per sample; each sample archive is extracted into its own directory.
#
# Usage:
#   ./scripts/01_download_data.sh [output_dir]
#
# Examples:
#   ./scripts/01_download_data.sh
#   ./scripts/01_download_data.sh raw
#   DOWNLOAD_BULK=0 ./scripts/01_download_data.sh raw
#
# Default output_dir: ./raw
# Default DOWNLOAD_BULK: 1

set -euo pipefail

OUTDIR="${1:-raw}"
DOWNLOAD_BULK="${DOWNLOAD_BULK:-1}"

SCRNA_URL="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE223060&format=file"
BULK_URL="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE223061&format=file"

SCRNA_ARCHIVE="$OUTDIR/GSE223060_RAW.tar"
BULK_ARCHIVE="$OUTDIR/GSE223061_RAW.tar"
SCRNA_OUTER_DIR="$OUTDIR/unpacked"
SCRNA_SAMPLE_DIR="$OUTDIR/samples"
BULK_DIR="$OUTDIR/unpacked_bulk"

mkdir -p "$OUTDIR"

require_command() {
    if ! command -v "$1" >/dev/null 2>&1; then
        echo "ERROR: required command not found: $1" >&2
        exit 1
    fi
}

require_command wget
require_command tar
require_command find

archive_is_readable() {
    tar -tf "$1" >/dev/null 2>&1
}

download_if_missing() {
    local url="$1"
    local destination="$2"
    local label="$3"

    if [[ -s "$destination" ]] && archive_is_readable "$destination"; then
        echo "== $label already present; skipping download =="
        echo "   $destination"
        return
    fi

    if [[ -e "$destination" ]]; then
        echo "== Existing $label archive is empty or incomplete; resuming download =="
    else
        echo "== Downloading $label =="
    fi

    wget -c "$url" -O "$destination"

    if ! archive_is_readable "$destination"; then
        echo "ERROR: downloaded archive is not readable: $destination" >&2
        exit 1
    fi
}

extract_outer_archive() {
    local archive="$1"
    local destination="$2"
    local marker="$destination/.outer_extraction_complete"
    local label="$3"

    mkdir -p "$destination"

    if [[ -f "$marker" ]]; then
        echo "== $label outer archive already extracted; skipping =="
        return
    fi

    echo "== Extracting $label outer archive =="
    tar -xf "$archive" -C "$destination"
    touch "$marker"
}

extract_sample_archives() {
    local source_dir="$1"
    local destination_root="$2"
    local found=0
    local extracted=0
    local skipped=0
    local failed=0

    mkdir -p "$destination_root"

    echo "== Extracting per-sample scRNA-seq archives =="

    while IFS= read -r -d '' archive; do
        found=$((found + 1))

        local filename sample destination marker
        filename="$(basename "$archive")"
        sample="${filename%.tar.gz}"
        destination="$destination_root/$sample"
        marker="$destination/.extraction_complete"

        if [[ -f "$marker" ]]; then
            echo "SKIP  $sample"
            skipped=$((skipped + 1))
            continue
        fi

        mkdir -p "$destination"
        echo "EXTRACT  $sample"

        if tar -xzf "$archive" -C "$destination"; then
            touch "$marker"
            extracted=$((extracted + 1))
        else
            echo "ERROR: failed to extract $archive" >&2
            failed=$((failed + 1))
        fi
    done < <(find "$source_dir" -maxdepth 1 -type f -name '*.tar.gz' -print0 | sort -z)

    if [[ "$found" -eq 0 ]]; then
        echo "ERROR: no per-sample .tar.gz archives found in $source_dir" >&2
        exit 1
    fi

    echo ""
    echo "Per-sample extraction summary:"
    echo "  Archives found : $found"
    echo "  Newly extracted: $extracted"
    echo "  Already present: $skipped"
    echo "  Failed         : $failed"

    if [[ "$failed" -gt 0 ]]; then
        exit 1
    fi
}

# 1. Download only when an archive is absent, empty, corrupt, or incomplete.
download_if_missing "$SCRNA_URL" "$SCRNA_ARCHIVE" "GSE223060 scRNA-seq"

if [[ "$DOWNLOAD_BULK" == "1" ]]; then
    download_if_missing "$BULK_URL" "$BULK_ARCHIVE" "GSE223061 bulk RNA-seq"
else
    echo "== Skipping GSE223061 bulk RNA-seq download (DOWNLOAD_BULK=$DOWNLOAD_BULK) =="
fi

# 2. Extract the GEO outer archives only once.
extract_outer_archive "$SCRNA_ARCHIVE" "$SCRNA_OUTER_DIR" "GSE223060 scRNA-seq"

if [[ "$DOWNLOAD_BULK" == "1" ]]; then
    extract_outer_archive "$BULK_ARCHIVE" "$BULK_DIR" "GSE223061 bulk RNA-seq"
fi

# 3. Extract every GSE223060 sample archive into raw/samples/<sample>/.
extract_sample_archives "$SCRNA_OUTER_DIR" "$SCRNA_SAMPLE_DIR"

echo ""
echo "== Final scRNA-seq layout =="
echo "Outer sample archives: $SCRNA_OUTER_DIR"
echo "Extracted samples     : $SCRNA_SAMPLE_DIR"
echo "Sample directories    : $(find "$SCRNA_SAMPLE_DIR" -mindepth 1 -maxdepth 1 -type d | wc -l)"
echo "Barcode files         : $(find "$SCRNA_SAMPLE_DIR" -type f -iname '*barcodes*' | wc -l)"
echo "Feature/gene files    : $(find "$SCRNA_SAMPLE_DIR" -type f \( -iname '*features*' -o -iname '*genes*' \) | wc -l)"
echo "Matrix files          : $(find "$SCRNA_SAMPLE_DIR" -type f -iname '*matrix*' | wc -l)"
echo "H5 files              : $(find "$SCRNA_SAMPLE_DIR" -type f -iname '*.h5' | wc -l)"

echo ""
echo "Done. Build the manifest from: $SCRNA_SAMPLE_DIR"
