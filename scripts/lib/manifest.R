# ==============================================================================
# lib/manifest.R -- read/validate raw/sample_manifest.csv, derive patient_id
# ==============================================================================

MANIFEST_REQUIRED_COLS <- c("sample_id", "format", "matrix_path",
                            "barcodes_path", "genefeat_path")
PATH_COLS <- c("matrix_path", "barcodes_path", "genefeat_path")

# ------------------------------------------------------------------------------
# patient_id derivation  [PROVISIONAL -- see CLAUDE.md open questions]
#
# Supplementary Table S1 is not yet in this repo, so patient mapping is derived
# from sample naming. Rule: strip a trailing _<digits> ONLY when the stem is
# purely numeric -- those are WashU cohort fraction/timepoint suffixes
# (27522_1 .. 27522_6 -> patient 27522). Prefixed IDs keep their trailing
# digits, because there the digits ARE the identifier, not a suffix:
#   MMRF_1695  -> MMRF_1695   (not MMRF)
#   ND_083017  -> ND_083017   (not ND)
#   MMY18273   -> MMY18273
#   25183      -> 25183
#
# KNOWN UNRESOLVED: samples "83942" and "MMY83942" share a numeric stem and may
# be one patient under two identifiers. They are deliberately NOT merged here.
# Confirm against Table S1 before the per-patient aggregation in step 08 -- if
# they are the same patient their malignant cells must be pooled, not ranked as
# two independent entries.
# ------------------------------------------------------------------------------
derive_patient_id <- function(sample_id) {
  sample_name <- strip_gsm(sample_id)
  sub("^([0-9]+)_[0-9]+$", "\\1", sample_name)
}

strip_gsm <- function(sample_id) sub("^GSM[0-9]+_", "", sample_id)

# Normal bone marrow controls (BM2/BM4/BM5/BM6) are not MM patients. Tagged so
# steps 06-08 can exclude them from the escape-risk ranking while step 09
# (CellChat) may still use them as a non-malignant microenvironment baseline.
is_normal_bm_sample <- function(sample_name) grepl("^BM[0-9]+$", sample_name)

# ------------------------------------------------------------------------------
# read_manifest() -- load, validate, resolve paths, annotate.
#
# Validates eagerly: a missing file or bad format should fail here, in seconds,
# not 40 minutes into the per-sample load loop.
# ------------------------------------------------------------------------------
read_manifest <- function(path = MANIFEST_PATH, project_root = PROJECT_ROOT) {
  log_msg("Reading manifest: ", path)
  mf <- utils::read.csv(path, stringsAsFactors = FALSE)
  log_msg("  manifest rows: ", nrow(mf))

  missing_cols <- setdiff(MANIFEST_REQUIRED_COLS, colnames(mf))
  if (length(missing_cols)) {
    stop("Manifest is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  bad_format <- mf$sample_id[mf$format != "triplet-ok"]
  if (length(bad_format)) {
    stop("Manifest rows not in 'triplet-ok' format: ",
         paste(bad_format, collapse = ", "))
  }

  if (anyDuplicated(mf$sample_id)) {
    stop("Duplicate sample_id in manifest: ",
         paste(unique(mf$sample_id[duplicated(mf$sample_id)]), collapse = ", "))
  }

  # manifest paths are project-root-relative
  for (col in PATH_COLS) mf[[col]] <- file.path(project_root, mf[[col]])

  missing_files <- unlist(lapply(PATH_COLS, function(col) mf[[col]][!file.exists(mf[[col]])]))
  if (length(missing_files)) {
    stop("Manifest references files that do not exist:\n  ",
         paste(utils::head(missing_files, 20), collapse = "\n  "))
  }
  log_msg("  verified ", nrow(mf) * length(PATH_COLS), " referenced files exist")

  annotate_manifest(mf)
}

annotate_manifest <- function(mf) {
  mf$sample_name  <- strip_gsm(mf$sample_id)
  mf$patient_id   <- derive_patient_id(mf$sample_id)
  mf$is_normal_bm <- is_normal_bm_sample(mf$sample_name)
  mf
}

# ------------------------------------------------------------------------------
# apply_exclusions() -- drop samples listed in EXCLUDE_SAMPLES, logging why.
# ------------------------------------------------------------------------------
apply_exclusions <- function(mf, exclusions = EXCLUDE_SAMPLES) {
  hit <- mf$sample_id %in% names(exclusions)
  for (sid in mf$sample_id[hit]) {
    log_msg("  EXCLUDING ", sid, " -- ", exclusions[[sid]])
  }
  unknown <- setdiff(names(exclusions), mf$sample_id)
  if (length(unknown)) {
    warning("EXCLUDE_SAMPLES names not present in manifest (typo?): ",
            paste(unknown, collapse = ", "))
  }
  mf[!hit, , drop = FALSE]
}

subset_for_smoke_test <- function(mf, n = SMOKE_N_SAMPLES) {
  if (!isTRUE(SMOKE_TEST)) return(mf)
  utils::head(mf, n)
}

# ------------------------------------------------------------------------------
# summarise_cohort() -- log patient structure and write the provisional map for
# diffing against Supplementary Table S1 when it lands.
# ------------------------------------------------------------------------------
summarise_cohort <- function(mf, write_map = TRUE) {
  log_msg("  samples: ", nrow(mf),
          " | patients (provisional): ", length(unique(mf$patient_id)),
          " | normal BM controls: ", sum(mf$is_normal_bm))

  multi <- table(mf$patient_id)
  multi <- multi[multi > 1]
  if (length(multi)) {
    log_msg("  multi-sample patients: ",
            paste(sprintf("%s(n=%d)", names(multi), as.integer(multi)),
                  collapse = ", "))
  }

  if (write_map) {
    out <- file.path(RESULTS_DIR,
                     paste0("sample_patient_map_PROVISIONAL", smoke_suffix(), ".csv"))
    utils::write.csv(
      mf[, c("sample_id", "sample_name", "patient_id", "is_normal_bm")],
      out, row.names = FALSE)
    log_msg("  wrote ", basename(out), " -- diff against Table S1 before step 08")
  }
  invisible(mf)
}
