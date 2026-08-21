# ==============================================================================
# lib/utils.R -- logging, env-var config overrides, checkpoint I/O
#
# Sourced FIRST by lib/init.R: 00_config.R depends on env_num()/env_flag().
# ==============================================================================

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# --- config overrides from the environment ------------------------------------
env_num <- function(name, default) {
  v <- Sys.getenv(name, unset = NA_character_)
  if (is.na(v) || !nzchar(v)) return(default)
  n <- suppressWarnings(as.numeric(v))
  if (is.na(n)) stop(sprintf("Env var %s='%s' is not numeric", name, v))
  n
}

env_flag <- function(name, default) {
  v <- Sys.getenv(name, unset = NA_character_)
  if (is.na(v) || !nzchar(v)) return(default)
  tolower(v) %in% c("1", "true", "yes", "y", "t")
}

env_chr <- function(name, default) {
  v <- Sys.getenv(name, unset = NA_character_)
  if (is.na(v) || !nzchar(v)) default else v
}

# --- logging ------------------------------------------------------------------
log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
  flush.console()
}

log_step <- function(...) {
  msg <- paste0(...)
  cat(sprintf("\n[%s] === %s ===\n", format(Sys.time(), "%H:%M:%S"), msg))
  flush.console()
}

# Wrap a long-running block so the log records how long it actually took --
# useful when deciding whether a step is worth checkpointing separately.
timed <- function(label, expr) {
  t0 <- Sys.time()
  log_msg(label, " ...")
  res <- force(expr)
  log_msg(label, sprintf(" done (%.2f min)",
                         as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  res
}

# --- filesystem ---------------------------------------------------------------
ensure_dir <- function(...) {
  for (d in c(...)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
  invisible(NULL)
}

require_file <- function(path, hint = NULL) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path,
         if (!is.null(hint)) paste0("\n  ", hint) else "")
  }
  invisible(path)
}

# --- checkpoints --------------------------------------------------------------
# SMOKE_TEST artefacts get a distinct suffix so a smoke run can never overwrite
# a real one, and 04b/04c never silently consume a 4-sample checkpoint.
smoke_suffix <- function() if (isTRUE(SMOKE_TEST)) "_SMOKETEST" else ""

ckpt_path <- function(name) {
  file.path(CKPT_DIR, paste0(name, smoke_suffix(), ".rds"))
}

# Per-sample checkpoint directory (04a writes, 04b reads). Kept separate from
# the monolithic checkpoints so an interrupted load is restartable per sample.
sample_ckpt_dir <- function() paste0(SAMPLE_CKPT_DIR, smoke_suffix())

sample_ckpt_obj <- function(sid) file.path(sample_ckpt_dir(), paste0(sid, ".rds"))
sample_ckpt_qc  <- function(sid) file.path(sample_ckpt_dir(), paste0(sid, ".qc.csv"))

save_ckpt <- function(obj, name) {
  ensure_dir(CKPT_DIR)
  p <- ckpt_path(name)
  timed(paste0("Writing checkpoint ", basename(p)), saveRDS(obj, p))
  invisible(p)
}

load_ckpt <- function(name, produced_by) {
  p <- ckpt_path(name)
  require_file(p, hint = paste0("Run `Rscript scripts/", produced_by, "` first."))
  timed(paste0("Reading checkpoint ", basename(p)), readRDS(p))
}

# --- provenance ---------------------------------------------------------------
# Stamped into the Seurat object's Misc slot at each stage so the final RDS
# carries a full record of how it was built.
provenance <- function(script, ...) {
  c(list(script = script,
         run_at  = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
         host    = Sys.info()[["nodename"]],
         smoke_test = SMOKE_TEST),
    list(...))
}
