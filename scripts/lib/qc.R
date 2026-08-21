# ==============================================================================
# lib/qc.R -- per-sample QC metrics, threshold filtering, doublet removal
#
# QC runs PER SAMPLE, before the merge, because:
#   - scDblFinder models doublet rate within a capture; pooling captures first
#     would mis-specify it, and
#   - filtering early keeps peak memory workable on a 32GB workstation
#     (~197k cells pre-QC across 61 samples).
# ==============================================================================

add_qc_metrics <- function(obj, mt_pattern = MT_PATTERN) {
  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = mt_pattern)
  obj
}

# ------------------------------------------------------------------------------
# filter_cells() -- apply the PROVISIONAL thresholds from lib/00_config.R.
# Callers should have written the pre-filter violin plot first (see lib/plots.R)
# so the thresholds can be reviewed against the data rather than assumed.
# ------------------------------------------------------------------------------
filter_cells <- function(obj,
                         min_features = MIN_FEATURES,
                         max_features = MAX_FEATURES,
                         max_mt       = MAX_PERCENT_MT) {
  keep <- obj$nFeature_RNA > min_features &
          obj$nFeature_RNA < max_features &
          obj$percent.mt   < max_mt
  subset(obj, cells = colnames(obj)[keep])
}

# ------------------------------------------------------------------------------
# remove_doublets() -- scDblFinder on a single capture.
#
# Returns list(obj, n_doublets, status). scDblFinder can fail on very small or
# low-complexity captures; rather than losing the sample we keep all cells,
# label them NA, and record the failure so it surfaces in the QC summary.
# ------------------------------------------------------------------------------
remove_doublets <- function(obj, sample_id = NULL) {
  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = SeuratObject::LayerData(obj, layer = "counts"))
  )

  res <- tryCatch(
    scDblFinder::scDblFinder(sce, verbose = FALSE),
    error = function(e) e
  )
  rm(sce); invisible(gc())

  if (inherits(res, "error")) {
    warning("scDblFinder failed for ", sample_id %||% "sample", ": ",
            conditionMessage(res), " -- keeping all cells, class set to NA.")
    obj$scDblFinder_class <- NA_character_
    return(list(obj = obj, n_doublets = NA_integer_, status = "dbl_failed"))
  }

  cls <- as.character(res$scDblFinder.class)
  names(cls) <- colnames(res)
  obj$scDblFinder_class <- unname(cls[colnames(obj)])
  rm(res); invisible(gc())

  n_doublets <- sum(obj$scDblFinder_class == "doublet", na.rm = TRUE)
  singlets   <- colnames(obj)[!is.na(obj$scDblFinder_class) &
                                obj$scDblFinder_class == "singlet"]

  list(obj        = subset(obj, cells = singlets),
       n_doublets = n_doublets,
       status     = "ok")
}

# ------------------------------------------------------------------------------
# qc_record() -- one row of the per-sample QC audit trail.
# ------------------------------------------------------------------------------
qc_record <- function(sample_id, patient_id,
                      n_cells_raw, n_cells_loaded, n_cells_thresholded,
                      n_doublets, n_cells_final,
                      median_nFeature, median_percent_mt, status) {
  data.frame(
    sample_id           = sample_id,
    patient_id          = patient_id,
    n_cells_raw         = n_cells_raw,
    n_cells_loaded      = n_cells_loaded,
    n_cells_thresholded = n_cells_thresholded,
    n_doublets          = n_doublets,
    n_cells_final       = n_cells_final,
    median_nFeature     = median_nFeature,
    median_percent_mt   = median_percent_mt,
    status              = status,
    stringsAsFactors    = FALSE
  )
}

# ------------------------------------------------------------------------------
# process_sample_qc() -- full per-sample path: load -> metrics -> plot ->
# threshold -> doublets. Returns list(obj, qc). obj is NULL if the sample was
# dropped for having too few surviving cells.
# ------------------------------------------------------------------------------
process_sample_qc <- function(row, common_genes, plot_dir = QC_DIR) {
  obj    <- load_sample(row, common_genes)
  n_raw  <- attr(obj, "n_cells_raw")
  obj    <- add_qc_metrics(obj)

  n_loaded <- ncol(obj)
  med_feat <- stats::median(obj$nFeature_RNA)
  med_mt   <- stats::median(obj$percent.mt)

  plot_qc_violin(obj, row$sample_id, plot_dir)

  obj <- filter_cells(obj)
  n_thresh <- ncol(obj)

  if (n_thresh < MIN_CELLS_KEEP) {
    log_msg("    DROPPED -- ", n_thresh, " cells survived thresholds (< ",
            MIN_CELLS_KEEP, ")")
    return(list(
      obj = NULL,
      qc  = qc_record(row$sample_id, row$patient_id, n_raw, n_loaded, n_thresh,
                      NA_integer_, 0L, med_feat, med_mt, "dropped_low_cells")
    ))
  }

  dbl <- remove_doublets(obj, row$sample_id)
  obj <- dbl$obj
  n_final <- ncol(obj)

  log_msg(sprintf("    raw %d -> loaded %d -> QC %d -> singlets %d  (%s doublets, mt median %.1f%%)",
                  n_raw, n_loaded, n_thresh, n_final,
                  ifelse(is.na(dbl$n_doublets), "NA", dbl$n_doublets), med_mt))

  list(
    obj = obj,
    qc  = qc_record(row$sample_id, row$patient_id, n_raw, n_loaded, n_thresh,
                    dbl$n_doublets, n_final, med_feat, med_mt,
                    if (dbl$status == "ok") "kept" else dbl$status)
  )
}
