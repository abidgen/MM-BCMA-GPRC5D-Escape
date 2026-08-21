# ==============================================================================
# lib/init.R -- single entry point for every driver script.
#
#   source("scripts/lib/init.R")
#
# Loads packages, sources every lib module in dependency order, seeds the RNG,
# and creates the output directories. Driver scripts source this and nothing else.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
})

.lib_dir <- local({
  # Resolve lib/ whether invoked via Rscript, R CMD, or sourced interactively.
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    d <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    cand <- file.path(d, "lib")
    if (dir.exists(cand)) return(cand)     # driver in scripts/, lib in scripts/lib
    if (basename(d) == "lib") return(d)    # sourced from within lib/
  }
  for (cand in c("scripts/lib", "lib", ".")) {
    if (file.exists(file.path(cand, "utils.R"))) return(normalizePath(cand))
  }
  stop("Cannot locate scripts/lib -- run driver scripts from the project root.")
})

# utils.R first: 00_config.R depends on env_num()/env_flag().
for (f in c("utils.R", "00_config.R", "manifest.R", "gene_space.R",
            "load_sample.R", "qc.R", "plots.R")) {
  source(file.path(.lib_dir, f))
}

set.seed(RANDOM_SEED)
options(future.globals.maxSize = FUTURE_MAX_GLOBALS_GB * 1024^3)

ensure_dir(RESULTS_DIR, QC_DIR, CKPT_DIR)

require_file(
  MANIFEST_PATH,
  hint = paste0("Expected the manifest built by scripts/03_build_manifest.py. ",
                "Run driver scripts from the project root, or set MM_PROJECT_ROOT.")
)

log_msg("project root : ", PROJECT_ROOT)
if (isTRUE(SMOKE_TEST)) {
  log_msg("*** SMOKE_TEST mode -- ", SMOKE_N_SAMPLES,
          " samples, artefacts suffixed _SMOKETEST ***")
}
