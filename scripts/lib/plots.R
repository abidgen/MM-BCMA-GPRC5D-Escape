# ==============================================================================
# lib/plots.R -- QC and diagnostic plot helpers
#
# All plotting is isolated here so the driver scripts stay analysis-only and a
# headless/plot-free run is a one-line change.
# ==============================================================================

PLOT_DPI <- 110

save_plot <- function(p, filename, dir = QC_DIR, width = 6, height = 5,
                      dpi = PLOT_DPI) {
  ensure_dir(dir)
  path <- file.path(dir, filename)
  suppressMessages(
    ggplot2::ggsave(path, p, width = width, height = height, dpi = dpi)
  )
  invisible(path)
}

# ------------------------------------------------------------------------------
# plot_qc_violin() -- PRE-FILTER distributions for one sample.
#
# Written before thresholds are applied so the provisional cutoffs in
# lib/00_config.R can be reviewed against real distributions (CLAUDE.md asks
# for per-sample inspection before locking them in). Threshold lines are drawn
# on the panels so over/under-cutting is visible at a glance.
# ------------------------------------------------------------------------------
plot_qc_violin <- function(obj, sample_id, dir = QC_DIR) {
  p <- Seurat::VlnPlot(obj,
                       features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                       ncol = 3, pt.size = 0)

  if (requireNamespace("patchwork", quietly = TRUE)) {
    p <- p + patchwork::plot_annotation(
      title = sprintf("%s  (pre-filter, n=%d cells)", sample_id, ncol(obj)),
      subtitle = sprintf("thresholds: nFeature %g-%g, percent.mt < %g",
                         MIN_FEATURES, MAX_FEATURES, MAX_PERCENT_MT))
  }

  save_plot(p, paste0("qc_prefilter_", sample_id, ".png"),
            dir = dir, width = 10, height = 4)
}

# ------------------------------------------------------------------------------
# plot_qc_summary() -- cohort-level view of the per-sample QC table, so an
# outlier sample is obvious without opening 61 violin plots.
# ------------------------------------------------------------------------------
plot_qc_summary <- function(qc_df, dir = QC_DIR) {
  d <- qc_df[order(qc_df$n_cells_final), ]
  d$sample_id <- factor(d$sample_id, levels = d$sample_id)

  p_cells <- ggplot2::ggplot(d, ggplot2::aes(x = sample_id, y = n_cells_final,
                                             fill = status)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Cells surviving QC per sample",
                  x = NULL, y = "cells after threshold + doublet removal") +
    ggplot2::theme_bw(base_size = 7)
  save_plot(p_cells, paste0("qc_cells_per_sample", smoke_suffix(), ".png"),
            dir = dir, width = 7, height = max(4, nrow(d) * 0.13))

  p_mt <- ggplot2::ggplot(d, ggplot2::aes(x = median_percent_mt,
                                          y = median_nFeature)) +
    ggplot2::geom_point(ggplot2::aes(size = n_cells_final), alpha = 0.6) +
    ggplot2::geom_vline(xintercept = MAX_PERCENT_MT, linetype = "dashed",
                        colour = "red") +
    ggplot2::labs(title = "Per-sample QC medians",
                  subtitle = "dashed line = percent.mt threshold",
                  x = "median percent.mt", y = "median nFeature_RNA") +
    ggplot2::theme_bw(base_size = 9)
  save_plot(p_mt, paste0("qc_medians_scatter", smoke_suffix(), ".png"),
            dir = dir, width = 6, height = 5)
}

plot_elbow <- function(obj, n_pcs = N_PCS, dir = QC_DIR) {
  save_plot(Seurat::ElbowPlot(obj, ndims = n_pcs),
            paste0("elbow_plot", smoke_suffix(), ".png"),
            dir = dir, width = 6, height = 4)
}

# ------------------------------------------------------------------------------
# plot_umap_by() -- UMAP coloured by a metadata column. Legend is suppressed for
# high-cardinality groupings (patient_id has ~45 levels).
# ------------------------------------------------------------------------------
plot_umap_by <- function(obj, group_by, reduction, filename, title,
                         dir = QC_DIR, label = FALSE) {
  n_levels <- length(unique(obj[[group_by, drop = TRUE]]))
  p <- Seurat::DimPlot(obj, reduction = reduction, group.by = group_by,
                       label = label) +
    ggplot2::ggtitle(title)
  if (n_levels > 20) p <- p + Seurat::NoLegend()
  save_plot(p, filename, dir = dir, width = 6.5, height = 5.5)
}
