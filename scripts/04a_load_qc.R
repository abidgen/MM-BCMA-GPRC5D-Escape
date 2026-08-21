#!/usr/bin/env Rscript
# ==============================================================================
# 04a_load_qc.R -- load every sample, QC it, remove doublets, checkpoint it.
#
# MM Dual-Antigen (BCMA/GPRC5D) Escape Risk Analysis -- stage 1 of 3.
#
#   04a_load_qc.R  ->  04b_merge_normalize.R  ->  04c_integrate_cluster.R
#
# This is the expensive, IO-bound stage (~61 samples, ~197k cells). Each sample
# is checkpointed individually, so an interrupted run resumes where it stopped
# (RESUME=1, the default) and a single bad sample can be re-done in isolation.
#
# Usage:
#   mamba activate mm-dual-antigen
#   Rscript scripts/04a_load_qc.R                 # full run
#   SMOKE_TEST=1 Rscript scripts/04a_load_qc.R    # 4 samples, ~2 min
#   RESUME=0 Rscript scripts/04a_load_qc.R        # force re-load everything
#   MAX_PERCENT_MT=10 Rscript scripts/04a_load_qc.R   # threshold sweep
#
# Inputs : raw/sample_manifest.csv
# Outputs: results/checkpoints/04a_samples/<sample_id>.rds  (+ .qc.csv)
#          results/checkpoints/04a_meta.rds
#          results/qc_summary_per_sample.csv
#          results/sample_patient_map_PROVISIONAL.csv
#          results/qc/qc_prefilter_<sample_id>.png
# ==============================================================================

source("scripts/lib/init.R")

suppressPackageStartupMessages({
  library(scDblFinder)
  library(SingleCellExperiment)
})

log_step("04a: load + QC")

# ------------------------------------------------------------------------------
# 1. Manifest -> cohort definition
# ------------------------------------------------------------------------------
mf <- read_manifest()
mf <- apply_exclusions(mf)
mf <- subset_for_smoke_test(mf)
mf <- summarise_cohort(mf)

# ------------------------------------------------------------------------------
# 2. Reconcile gene space across the differing CellRanger references.
#    (See lib/gene_space.R for why this is an intersection, not a union.)
# ------------------------------------------------------------------------------
gs <- build_gene_space(mf)
common_genes <- gs$common_genes
mf$n_genes_ref <- as.integer(gs$ref_sizes[mf$sample_id])

# Cached per-sample checkpoints are only reusable if they were built with the
# same gene space AND the same QC thresholds. Without this check a threshold
# sweep run with RESUME=1 would silently mix samples filtered at different
# cutoffs, which is worse than either cutoff alone. Refuse rather than resume.
load_fingerprint <- list(
  n_common_genes = length(common_genes),
  min_cells_gene = MIN_CELLS_GENE,
  min_features   = MIN_FEATURES,
  max_features   = MAX_FEATURES,
  max_percent_mt = MAX_PERCENT_MT,
  min_cells_keep = MIN_CELLS_KEEP
)

meta_path <- ckpt_path("04a_meta")
if (isTRUE(RESUME) && file.exists(meta_path)) {
  prev <- readRDS(meta_path)
  drift <- character(0)
  if (!identical(prev$common_genes, common_genes)) {
    drift <- c(drift, sprintf("gene space (%d genes cached vs %d now)",
                              length(prev$common_genes), length(common_genes)))
  }
  if (!identical(prev$load_fingerprint, load_fingerprint)) {
    changed <- names(load_fingerprint)[
      !mapply(identical, load_fingerprint, prev$load_fingerprint[names(load_fingerprint)])]
    drift <- c(drift, paste0("QC parameters: ", paste(changed, collapse = ", ")))
  }
  if (length(drift)) {
    stop("Cached checkpoints do not match the current configuration:\n  - ",
         paste(drift, collapse = "\n  - "),
         "\nDelete ", sample_ckpt_dir(), " and ", meta_path,
         " and re-run, or set RESUME=0.")
  }
}

# ------------------------------------------------------------------------------
# 3. Per-sample: load -> QC metrics -> threshold -> doublet removal -> checkpoint
# ------------------------------------------------------------------------------
ensure_dir(sample_ckpt_dir())
qc_rows <- list()

for (i in seq_len(nrow(mf))) {
  row <- mf[i, ]
  sid <- row$sample_id
  log_msg(sprintf("[%2d/%2d] %s", i, nrow(mf), sid))

  # --- resume: reuse a completed sample rather than re-loading it ------------
  if (isTRUE(RESUME) && file.exists(sample_ckpt_obj(sid)) &&
      file.exists(sample_ckpt_qc(sid))) {
    qc_rows[[sid]] <- utils::read.csv(sample_ckpt_qc(sid), stringsAsFactors = FALSE)
    log_msg("    cached -- skipping (RESUME=1)")
    next
  }

  res <- process_sample_qc(row, common_genes)
  qc_rows[[sid]] <- res$qc

  # QC row is written alongside the object so the cohort summary can be rebuilt
  # on resume without deserialising 61 Seurat objects.
  utils::write.csv(res$qc, sample_ckpt_qc(sid), row.names = FALSE)

  if (is.null(res$obj)) {
    # dropped for too few cells -- record the qc row, write no object
    if (file.exists(sample_ckpt_obj(sid))) unlink(sample_ckpt_obj(sid))
  } else {
    saveRDS(res$obj, sample_ckpt_obj(sid))
  }

  rm(res); invisible(gc())
}

# ------------------------------------------------------------------------------
# 4. Cohort QC summary
# ------------------------------------------------------------------------------
qc_df <- do.call(rbind, qc_rows)
rownames(qc_df) <- NULL

qc_out <- file.path(RESULTS_DIR,
                    paste0("qc_summary_per_sample", smoke_suffix(), ".csv"))
utils::write.csv(qc_df, qc_out, row.names = FALSE)
log_msg("Wrote ", basename(qc_out))

plot_qc_summary(qc_df)

kept <- qc_df[qc_df$n_cells_final > 0, ]
log_msg("Samples kept: ", nrow(kept), "/", nrow(qc_df),
        " | cells retained: ", sum(kept$n_cells_final),
        " (from ", sum(qc_df$n_cells_raw), " raw)")

dropped <- qc_df[qc_df$n_cells_final == 0, ]
if (nrow(dropped)) {
  log_msg("Dropped samples: ", paste(dropped$sample_id, collapse = ", "))
}
failed_dbl <- qc_df[qc_df$status == "dbl_failed", ]
if (nrow(failed_dbl)) {
  log_msg("scDblFinder failed (all cells kept, class NA): ",
          paste(failed_dbl$sample_id, collapse = ", "))
}

if (nrow(kept) < 2) stop("Fewer than 2 samples survived QC -- aborting.")

# ------------------------------------------------------------------------------
# 5. Stage metadata for 04b
# ------------------------------------------------------------------------------
save_ckpt(list(
  manifest         = mf,
  common_genes     = common_genes,
  load_fingerprint = load_fingerprint,
  gene_ref_tbl  = gs$n_by_ref,
  qc_summary    = qc_df,
  kept_samples  = kept$sample_id,
  provenance    = provenance(
    "scripts/04a_load_qc.R",
    n_samples_manifest = nrow(mf) + length(EXCLUDE_SAMPLES),
    n_samples_loaded   = nrow(mf),
    excluded_samples   = EXCLUDE_SAMPLES,
    gene_space         = "intersection across retained samples",
    n_genes_common     = length(common_genes),
    qc_thresholds      = list(min_features   = MIN_FEATURES,
                              max_features   = MAX_FEATURES,
                              max_percent_mt = MAX_PERCENT_MT,
                              min_cells_gene = MIN_CELLS_GENE,
                              min_cells_keep = MIN_CELLS_KEEP),
    patient_id_source  = "PROVISIONAL -- derived from sample names; confirm against Supplementary Table S1"
  )
), "04a_meta")

log_msg("04a done. Next: Rscript scripts/04b_merge_normalize.R")
