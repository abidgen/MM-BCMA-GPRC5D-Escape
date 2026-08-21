#!/usr/bin/env Rscript
# ==============================================================================
# 04c_integrate_cluster.R -- Harmony integration, Leiden clustering, UMAP.
#
# MM Dual-Antigen (BCMA/GPRC5D) Escape Risk Analysis -- stage 3 of 3.
#
#   04a_load_qc.R  ->  04b_merge_normalize.R  ->  04c_integrate_cluster.R
#
# Cheapest stage to iterate on: re-run at a different CLUSTER_RES without
# repeating the load or the PCA.
#
# Usage:
#   mamba activate mm-dual-antigen
#   Rscript scripts/04c_integrate_cluster.R
#   CLUSTER_RES=0.8 Rscript scripts/04c_integrate_cluster.R
#
# Inputs : results/checkpoints/04b_merged_pca.rds
# Outputs: results/integrated_object.rds        <- consumed by scripts/05
#          results/qc/umap_post_harmony_by_{patient,cluster,sample_type,reference}.png
#          results/cluster_composition.csv
# ==============================================================================

source("scripts/lib/init.R")

suppressPackageStartupMessages(library(harmony))

log_step("04c: Harmony + clustering + UMAP")

merged <- load_ckpt("04b_merged_pca", produced_by = "04b_merge_normalize.R")
log_msg("Loaded: ", ncol(merged), " cells x ", nrow(merged), " genes")

# ------------------------------------------------------------------------------
# 1. Harmony over patient-of-origin
#
# Correcting on patient_id (not sample_id) is deliberate: the multi-sample
# patients (e.g. 27522_1..27522_6) are fractions/timepoints of one marrow, and
# collapsing them into one batch avoids over-correcting away real within-patient
# structure. Note this mapping is PROVISIONAL until Supplementary Table S1 is in
# the repo -- see lib/manifest.R.
# ------------------------------------------------------------------------------
merged <- timed("RunHarmony (group.by.vars = 'patient_id')",
                RunHarmony(merged,
                           group.by.vars  = "patient_id",
                           reduction.use  = "pca",
                           reduction.save = "harmony",
                           verbose        = TRUE))

# ------------------------------------------------------------------------------
# 2. Neighbours + clustering
# ------------------------------------------------------------------------------
merged <- timed("FindNeighbors (harmony)",
                FindNeighbors(merged, reduction = "harmony", dims = 1:N_PCS,
                              verbose = FALSE))

# CLAUDE.md specifies Leiden (algorithm = 4). Seurat routes that through
# reticulate to the Python leidenalg module, which is a common install gap.
# Check before the call rather than dying after Harmony has already run, and
# record which algorithm actually produced the clusters.
use_leiden <- requireNamespace("reticulate", quietly = TRUE) &&
  reticulate::py_module_available("leidenalg")

if (use_leiden) {
  cluster_algo <- "leiden"
  merged <- timed(paste0("FindClusters (Leiden, res=", CLUSTER_RES, ")"),
                  FindClusters(merged, resolution = CLUSTER_RES, algorithm = 4,
                               verbose = FALSE))
} else {
  cluster_algo <- "louvain_fallback"
  warning("Python 'leidenalg' not available -- falling back to Louvain ",
          "(algorithm = 1). To get the Leiden clustering specified in ",
          "CLAUDE.md:\n  mamba install -n mm-dual-antigen -c conda-forge ",
          "leidenalg python-igraph\nthen re-run this script.")
  log_msg("!! leidenalg unavailable -- using Louvain fallback")
  merged <- timed(paste0("FindClusters (Louvain FALLBACK, res=", CLUSTER_RES, ")"),
                  FindClusters(merged, resolution = CLUSTER_RES, algorithm = 1,
                               verbose = FALSE))
}
log_msg("  clusters found: ", nlevels(merged$seurat_clusters),
        " (", cluster_algo, ")")

# ------------------------------------------------------------------------------
# 3. UMAP on the corrected embedding
# ------------------------------------------------------------------------------
merged <- timed("RunUMAP (harmony)",
                RunUMAP(merged, reduction = "harmony", dims = 1:N_PCS,
                        reduction.name = "umap", verbose = FALSE))

# ------------------------------------------------------------------------------
# 4. Integration diagnostics
#
# The reference-version plot is the one that matters most here: if clusters
# track n_genes_ref rather than biology, the gene-space intersection failed to
# neutralise the processing batch and the antigen calls downstream are suspect.
# ------------------------------------------------------------------------------
sfx <- smoke_suffix()
plot_umap_by(merged, "patient_id", "umap",
             paste0("umap_post_harmony_by_patient", sfx, ".png"),
             "Post-Harmony (by patient)")
plot_umap_by(merged, "seurat_clusters", "umap",
             paste0("umap_post_harmony_by_cluster", sfx, ".png"),
             paste0("Post-Harmony clusters (", cluster_algo, ", res=", CLUSTER_RES, ")"),
             label = TRUE)
plot_umap_by(merged, "is_normal_bm", "umap",
             paste0("umap_post_harmony_by_sample_type", sfx, ".png"),
             "Normal BM control vs. MM sample")
plot_umap_by(merged, "n_genes_ref", "umap",
             paste0("umap_post_harmony_by_reference", sfx, ".png"),
             "CellRanger reference version (batch check)")

# Per-cluster composition: a cluster drawn from a single patient is usually a
# malignant clone (expected for plasma cells) rather than a failure of
# integration -- worth being able to tell those apart in step 05.
comp <- as.data.frame.matrix(table(merged$seurat_clusters, merged$patient_id))
comp$cluster    <- rownames(comp)
comp$n_cells    <- rowSums(comp[, setdiff(colnames(comp), "cluster"), drop = FALSE])
comp$n_patients <- rowSums(comp[, setdiff(colnames(comp), c("cluster", "n_cells")),
                                drop = FALSE] > 0)
comp_out <- file.path(RESULTS_DIR, paste0("cluster_composition", sfx, ".csv"))
utils::write.csv(comp[, c("cluster", "n_cells", "n_patients")], comp_out,
                 row.names = FALSE)
log_msg("Wrote ", basename(comp_out))

# ------------------------------------------------------------------------------
# 5. Provenance + final save
# ------------------------------------------------------------------------------
Misc(merged, "pipeline_04c") <- provenance(
  "scripts/04c_integrate_cluster.R",
  harmony_group_by  = "patient_id",
  n_pcs             = N_PCS,
  cluster_algorithm = cluster_algo,
  cluster_res       = CLUSTER_RES,
  n_clusters        = nlevels(merged$seurat_clusters),
  session           = utils::capture.output(utils::sessionInfo())
)

out_path <- file.path(RESULTS_DIR, paste0("integrated_object", sfx, ".rds"))
invisible(timed(paste0("Saving ", basename(out_path)), saveRDS(merged, out_path)))

log_step("04 pipeline complete")
log_msg("  cells      : ", ncol(merged))
log_msg("  genes      : ", nrow(merged), " (intersected gene space)")
log_msg("  samples    : ", length(unique(merged$sample_id)))
log_msg("  patients   : ", length(unique(merged$patient_id)), " (PROVISIONAL mapping)")
log_msg("  clusters   : ", nlevels(merged$seurat_clusters), " (", cluster_algo, ")")
log_msg("  written to : ", out_path)
log_msg("Next: scripts/05_annotate_celltypes.R")
