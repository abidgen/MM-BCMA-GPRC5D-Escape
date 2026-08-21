# ==============================================================================
# lib/load_sample.R -- read one sample's matrix triplet into a Seurat object
# ==============================================================================

# ------------------------------------------------------------------------------
# read_sample_counts() -- ReadMtx(), NOT Read10X().
#
# This archive uses CellRanger v2-era naming: counts.mtx / genes.tsv / barcodes.tsv
# rather than matrix.mtx / features.tsv / barcodes.tsv. Read10X() hardcodes the
# modern filenames and will not find these files even when pointed at the right
# directory. (Established in CLAUDE.md as non-negotiable.)
#
# feature.column = 1 because genes.tsv is a single symbol column here; ReadMtx()
# defaults to column 2, which does not exist in these files.
# ------------------------------------------------------------------------------
read_sample_counts <- function(row) {
  Seurat::ReadMtx(
    mtx            = row$matrix_path,
    cells          = row$barcodes_path,
    features       = row$genefeat_path,
    feature.column = 1,
    cell.column    = 1,
    feature.sep    = "\t",
    cell.sep       = "\t"
  )
}

# ------------------------------------------------------------------------------
# make_sample_object() -- counts -> QC-ready Seurat object.
#
# Restricts to the reconciled gene space (see lib/gene_space.R) and prefixes
# barcodes with the sample_id so they stay unique across the 61-sample merge.
# ------------------------------------------------------------------------------
make_sample_object <- function(counts, row, common_genes,
                               min_cells    = MIN_CELLS_GENE,
                               min_features = MIN_FEATURES) {

  missing <- setdiff(common_genes, rownames(counts))
  if (length(missing)) {
    stop("Sample ", row$sample_id, " lacks ", length(missing),
         " genes from the common gene space (e.g. ",
         paste(utils::head(missing, 5), collapse = ", "),
         ") -- gene space was built from a different manifest.")
  }
  counts <- counts[common_genes, , drop = FALSE]

  colnames(counts) <- paste(row$sample_id, colnames(counts), sep = "|")

  # min.cells defaults to 0 here on purpose -- gene filtering happens once after
  # the merge (04b / MERGED_MIN_CELLS), not per sample. See lib/00_config.R for
  # why per-sample gene filtering would reintroduce structural zeros.
  obj <- Seurat::CreateSeuratObject(
    counts       = counts,
    project      = row$sample_id,
    min.cells    = min_cells,
    min.features = min_features
  )

  obj$sample_id    <- row$sample_id
  obj$sample_name  <- row$sample_name
  obj$patient_id   <- row$patient_id
  obj$is_normal_bm <- row$is_normal_bm
  # reference version retained as a covariate so 04c/05 can check that no
  # cluster is driven by processing batch rather than biology
  obj$n_genes_ref  <- row$n_genes_ref

  obj
}

# ------------------------------------------------------------------------------
# load_sample() -- convenience wrapper: manifest row -> Seurat object.
# ------------------------------------------------------------------------------
load_sample <- function(row, common_genes) {
  counts <- read_sample_counts(row)
  n_raw  <- ncol(counts)
  obj    <- make_sample_object(counts, row, common_genes)
  rm(counts); invisible(gc())
  attr(obj, "n_cells_raw") <- n_raw
  obj
}
