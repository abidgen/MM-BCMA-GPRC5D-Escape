#!/bin/bash
# Download and unpack GSE223060 (scRNA-seq) and GSE223061 (matched bulk RNA-seq, optional)
# from GEO's FTP. No SRA/raw FASTQ available for this series — this pulls the
# processed Cell Ranger output the authors deposited.
#
# Usage: ./01_download_data.sh [output_dir]
# Default output_dir: ./raw

set -euo pipefail

OUTDIR="${1:-raw}"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "== Downloading GSE223060 (scRNA-seq) =="
wget -c "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE223060&format=file" \
     -O GSE223060_RAW.tar

echo "== Downloading GSE223061 (matched bulk RNA-seq — optional, comment out if not needed) =="
wget -c "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE223061&format=file" \
     -O GSE223061_RAW.tar

echo "== Unpacking scRNA-seq archive =="
mkdir -p unpacked
tar -xf GSE223060_RAW.tar -C unpacked

echo "== Unpacking bulk RNA-seq archive =="
mkdir -p unpacked_bulk
tar -xf GSE223061_RAW.tar -C unpacked_bulk

echo "== Done. Contents of unpacked/ (first 20 entries): =="
ls unpacked/ | head -20

echo ""
echo "Total scRNA-seq files:"
find unpacked/ -type f | wc -l

echo ""
echo "Run the file-checking script next (02_check_files.sh) before writing any loading code —"
echo "GEO archives from this era can bundle either a 10x barcodes/features/matrix triplet"
echo "per sample or a single filtered_feature_bc_matrix.h5 per sample, and the two need"
echo "different Seurat loading calls."
