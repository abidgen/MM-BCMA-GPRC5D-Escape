#!/usr/bin/env Rscript
# ==============================================================================
# 04b_merge_normalize.R -- merge QC'd samples, normalize, select HVGs, PCA.
#
# MM Dual-Antigen (BCMA/GPRC5D) Escape Risk Analysis -- stage 2 of 3.
#
#   04a_load_qc.R  ->  04b_merge_normalize.R  ->  04c_integrate_cluster.R
#
# Split from 04c so that re-running Harmony or re-clustering at a different
# resolution does not repeat the merge/normalize/PCA work.
#
# Usage:
#   mamba activate mm-dual-antigen
#   Rscript scripts/04b_merge_normalize.R
#   PLOT_PRE_HARMONY_UMAP=1 Rscript scripts/04b_merge_normalize.R  # + diagnostic
#   N_VAR_FEATURES=3000 Rscript scripts/04b_merge_normalize.R
#
# Inputs : results/checkpoints/04a_meta.rds
#          results/checkpoints/04a_samples/<sample_id>.rds
# Outputs: results/checkpoints/04b_merged_pca.rds
#          results/qc/elbow_plot.png
#          results/qc/umap_pre_harmony_by_patient.png   (if PLOT_PRE_HARMONY_UMAP=1)
# ==============================================================================

source("scripts/lib/init.R")

log_step("04b: merge + normalize + PCA")

meta <- load_ckpt("04a_meta", produced_by = "04a_load_qc.R")
log_msg("Loaded stage-04a metadata: ", length(meta$kept_samples), " samples kept, ",
        length(meta$common_genes), " genes in common space")

# ------------------------------------------------------------------------------
# 1. Read per-sample checkpoints
# ------------------------------------------------------------------------------
objs <- vector("list", length(meta$kept_samples))
names(objs) <- meta$kept_samples

for (i in seq_along(meta$kept_samples)) {
  sid <- meta$kept_samples[i]
  p <- sample_ckpt_obj(sid)
  require_file(p, hint = "Re-run scripts/04a_load_qc.R for this sample.")
  objs[[sid]] <- readRDS(p)
  if (i %% 10 == 0 || i == length(meta$kept_samples)) {
    log_msg(sprintf("  loaded %d/%d sample checkpoints", i, length(meta$kept_samples)))
  }
}

# ncol() on a Seurat object returns a double, not an integer
n_cells_in <- sum(vapply(objs, ncol, numeric(1)))
log_msg("  total cells to merge: ", n_cells_in)

# ------------------------------------------------------------------------------
# 2. Merge
# ------------------------------------------------------------------------------
merged <- timed("Merging objects",
                merge(objs[[1]], y = objs[-1], project = "MM_dual_antigen"))
rm(objs); invisible(gc())

# Seurat v5 keeps one counts layer per merged sample. We integrate with Harmony
# (a single corrected embedding), not Seurat's layer-wise anchor integration, so
# join the layers into one matrix and run the standard workflow on it.
merged <- timed("JoinLayers", JoinLayers(merged))
log_msg("  merged object: ", nrow(merged), " genes x ", ncol(merged), " cells")

stopifnot(ncol(merged) == n_cells_in)

# ------------------------------------------------------------------------------
# 2b. Cohort-wide gene detection filter.
#
# Applied ONCE here rather than per sample in 04a. Filtering per sample drops a
# different gene set in each sample, and the merge then zero-fills those genes
# for every cell of the samples that dropped them -- a cell would read as
# antigen-negative because its sample filtered the gene out, not because the
# antigen is absent. Doing it once keeps the gene space identical for every cell.
#
# REQUIRED_GENES are exempt: a gene that is genuinely near-zero cohort-wide is
# still a result we need to report, not one to silently discard.
# ------------------------------------------------------------------------------
cm <- SeuratObject::LayerData(merged, layer = "counts")
nnz_per_gene <- if (inherits(cm, "dgCMatrix")) {
  tabulate(cm@i + 1L, nbins = nrow(cm))          # cheap: nonzeros per row
} else {
  Matrix::rowSums(cm > 0)
}
names(nnz_per_gene) <- rownames(cm)
rm(cm); invisible(gc())

keep_genes <- names(nnz_per_gene)[nnz_per_gene >= MERGED_MIN_CELLS]
forced <- setdiff(intersect(REQUIRED_GENES, rownames(merged)), keep_genes)
if (length(forced)) {
  log_msg("  retaining ", length(forced), " required gene(s) below the ",
          MERGED_MIN_CELLS, "-cell threshold: ", paste(forced, collapse = ", "))
  keep_genes <- union(keep_genes, forced)
}

log_msg("  gene filter (detected in >= ", MERGED_MIN_CELLS, " cells cohort-wide): ",
        nrow(merged), " -> ", length(keep_genes), " genes")
merged <- subset(merged, features = keep_genes)

# the object 05-07 will score against must still contain every required gene
assert_required_genes(rownames(merged), REQUIRED_GENES)

# ------------------------------------------------------------------------------
# 3. Normalize -> HVG -> scale -> PCA
# ------------------------------------------------------------------------------
merged <- timed("NormalizeData",
                NormalizeData(merged, normalization.method = "LogNormalize",
                              scale.factor = 1e4, verbose = FALSE))

merged <- timed(paste0("FindVariableFeatures (n=", N_VAR_FEATURES, ")"),
                FindVariableFeatures(merged, selection.method = "vst",
                                     nfeatures = N_VAR_FEATURES, verbose = FALSE))

# ScaleData defaults to variable features only -- deliberate. Scaling all ~22k
# genes would produce a dense ~22k x 150k matrix and exhaust 32GB of RAM.
merged <- timed("ScaleData (variable features only)",
                ScaleData(merged, verbose = FALSE))

merged <- timed(paste0("RunPCA (npcs=", N_PCS, ")"),
                RunPCA(merged, npcs = N_PCS, verbose = FALSE))

plot_elbow(merged)

# ------------------------------------------------------------------------------
# 4. Optional pre-Harmony UMAP -- makes the patient batch effect (and therefore
#    Harmony's correction in 04c) visible. Costs several minutes at full size.
# ------------------------------------------------------------------------------
if (isTRUE(PLOT_PRE_HARMONY_UMAP)) {
  merged <- timed("Pre-Harmony UMAP (diagnostic)",
                  RunUMAP(merged, dims = 1:N_PCS, reduction = "pca",
                          reduction.name = "umap_prehm", verbose = FALSE))
  plot_umap_by(merged, "patient_id", "umap_prehm",
               paste0("umap_pre_harmony_by_patient", smoke_suffix(), ".png"),
               "Pre-Harmony (by patient)")
} else {
  log_msg("Skipping pre-Harmony UMAP (set PLOT_PRE_HARMONY_UMAP=1 to enable)")
}

# ------------------------------------------------------------------------------
# 5. Checkpoint
# ------------------------------------------------------------------------------
Misc(merged, "pipeline_04a") <- meta$provenance
Misc(merged, "pipeline_04b") <- provenance(
  "scripts/04b_merge_normalize.R",
  n_cells          = ncol(merged),
  n_genes          = nrow(merged),
  n_samples        = length(unique(merged$sample_id)),
  n_patients       = length(unique(merged$patient_id)),
  n_var_features   = N_VAR_FEATURES,
  n_pcs            = N_PCS,
  pre_harmony_umap = PLOT_PRE_HARMONY_UMAP
)

save_ckpt(merged, "04b_merged_pca")

log_msg("04b done -- ", ncol(merged), " cells x ", nrow(merged), " genes, ",
        N_PCS, " PCs.")
log_msg("Next: Rscript scripts/04c_integrate_cluster.R")
