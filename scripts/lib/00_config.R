# ==============================================================================
# lib/00_config.R -- paths, QC thresholds, gene sets, sample exclusions
#
# Single source of truth for every tunable in the 04* pipeline. Driver scripts
# read these; they do not define their own. Any constant here can be overridden
# at run time by an environment variable of the same name, e.g.
#
#   MAX_PERCENT_MT=10 SMOKE_TEST=1 Rscript scripts/04a_load_qc.R
#
# so threshold sweeps do not require editing tracked files.
# ==============================================================================

# --- project root -------------------------------------------------------------
# Driver scripts are expected to run from the project root. MM_PROJECT_ROOT
# overrides that (useful for cron / out-of-tree invocation).
PROJECT_ROOT <- Sys.getenv("MM_PROJECT_ROOT", unset = getwd())

MANIFEST_PATH <- file.path(PROJECT_ROOT, "raw", "sample_manifest.csv")
RESULTS_DIR   <- file.path(PROJECT_ROOT, "results")
QC_DIR        <- file.path(RESULTS_DIR, "qc")
CKPT_DIR      <- file.path(RESULTS_DIR, "checkpoints")
SAMPLE_CKPT_DIR <- file.path(CKPT_DIR, "04a_samples")

# --- QC thresholds ------------------------------------------------------------
# PROVISIONAL (per CLAUDE.md): 04a writes pre-filter violin plots and a
# per-sample count table before applying these, so they can be revisited.
# Bone marrow aspirates vary more than solid-tumour dissociations.
MIN_FEATURES   <- env_num("MIN_FEATURES",   200)   # nFeature_RNA lower bound
MAX_FEATURES   <- env_num("MAX_FEATURES",   6000)  # nFeature_RNA upper bound
MAX_PERCENT_MT <- env_num("MAX_PERCENT_MT", 15)    # percent.mt upper bound
MIN_CELLS_KEEP <- env_num("MIN_CELLS_KEEP", 50)    # drop sample below N cells

# Gene detection filter. Deliberately 0 PER SAMPLE and applied once AFTER the
# merge instead (MERGED_MIN_CELLS, used in 04b).
#
# Filtering genes per sample would drop a different gene set in each sample --
# the merge then zero-fills those genes for every cell of the samples that
# dropped them. That is the same structural-zero artifact that lib/gene_space.R
# exists to prevent: a cell would read as antigen-negative because its sample
# filtered the gene out, not because the antigen is absent. Applying the filter
# once on the merged object keeps the gene space identical for every cell.
# (Costs essentially no memory: dropping all-zero rows from a sparse matrix
# saves nothing, since they store no values.)
MIN_CELLS_GENE   <- env_num("MIN_CELLS_GENE",   0)  # per sample -- keep at 0
MERGED_MIN_CELLS <- env_num("MERGED_MIN_CELLS", 3)  # post-merge, cohort-wide

# --- dimensionality / clustering ----------------------------------------------
N_PCS          <- env_num("N_PCS",          30)
N_VAR_FEATURES <- env_num("N_VAR_FEATURES", 2000)
CLUSTER_RES    <- env_num("CLUSTER_RES",    0.5)
RANDOM_SEED    <- env_num("RANDOM_SEED",    1234)

# --- run modes ----------------------------------------------------------------
# SMOKE_TEST: run the whole path on a handful of samples first (~2 min) before
# committing to the full ~197k-cell run. Writes *_SMOKETEST artefacts so a smoke
# run can never overwrite a real one.
SMOKE_TEST      <- env_flag("SMOKE_TEST",      FALSE)
SMOKE_N_SAMPLES <- env_num("SMOKE_N_SAMPLES",  4)

# RESUME: in 04a, skip samples whose checkpoint already exists. Makes an
# interrupted load restartable without repeating completed samples.
RESUME <- env_flag("RESUME", TRUE)

# Pre-Harmony UMAP in 04b is a diagnostic only (shows the batch effect Harmony
# then corrects). It costs several minutes on the full object; off by default.
PLOT_PRE_HARMONY_UMAP <- env_flag("PLOT_PRE_HARMONY_UMAP", FALSE)

# Seurat/future trips its 500MB default global-export limit on an object this size.
FUTURE_MAX_GLOBALS_GB <- env_num("FUTURE_MAX_GLOBALS_GB", 8)

# ==============================================================================
# Sample exclusions
#
# GSM6939056_56203_1 -- processed against a 22,184-gene reference (every other
#   sample uses 33,538 or 33,694). That truncated reference contains NO
#   TNFRSF17 (BCMA) row and no IGLC1/2/3 rows.
#
#   This project's core metric is the fraction of malignant cells NEGATIVE for
#   BCMA, and step 06 calls malignancy by kappa/lambda light-chain restriction.
#   Including this sample would score 100% of its cells BCMA-negative and 100%
#   of its plasma cells kappa-restricted -- both purely technical artefacts that
#   would place it at the top of the escape-risk ranking for the wrong reason.
#
#   Patient 56203 remains represented by GSM6939057_56203_2 (33,694-gene
#   reference, complete), so no patient is lost from the cohort.
# ==============================================================================
EXCLUDE_SAMPLES <- c(
  "GSM6939056_56203_1" = "truncated 22184-gene reference; missing TNFRSF17 and IGLC1/2/3"
)

# ==============================================================================
# Genes that MUST survive the cross-sample gene intersection (lib/gene_space.R)
# for steps 05-07 to be interpretable. 04a stops rather than silently building
# an object that cannot answer the project's question.
# ==============================================================================
REQUIRED_GENES <- c(
  # step 07 -- antigen scoring (primary readout)
  "TNFRSF17", "GPRC5D", "SLAMF7", "FCRL5",
  # step 06 -- light-chain restriction / malignant calling
  "IGKC", "IGLC1", "IGLC2", "IGLC3",
  # step 05 -- plasma cell identity
  "SDC1", "CD38", "MZB1", "XBP1", "IRF4"
)

# Marker panel for step 05, kept here so 05 and the 04c sanity plots agree.
MARKER_PANEL <- list(
  PlasmaCell = c("SDC1", "CD38", "MZB1", "XBP1", "IRF4"),
  Bcell      = c("MS4A1", "CD79A", "CD19"),
  Tcell      = c("CD3D", "CD3E", "CD8A", "CD4"),
  NK         = c("NCAM1", "NKG7", "GNLY"),
  Myeloid    = c("CD14", "LYZ", "ITGAM"),
  Erythroid  = c("HBB", "GYPA"),
  HSPC       = c("CD34", "KIT")
)

MT_PATTERN <- "^MT-"
