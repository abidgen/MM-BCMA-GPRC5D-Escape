# ==============================================================================
# lib/gene_space.R -- reconcile differing CellRanger references across samples
#
# WHY THIS MODULE EXISTS
#
# GSE223060 samples were not all processed against the same reference. Observed
# in this archive:
#     33,538 genes -- 37 samples
#     33,694 genes -- 24 samples
#     22,184 genes --  1 sample (GSM6939056_56203_1, excluded in 00_config.R)
# The two retained references share only ~22,164 symbols.
#
# Merging on the UNION of gene sets fills absent genes with zeros. For ~11k
# genes that means a STRUCTURAL zero across an entire cohort of samples, which:
#
#   (a) makes those genes spuriously "highly variable" along reference version
#       rather than biology -- and Harmony cannot correct a structural zero,
#       because there is no signal in the batch to align; and
#   (b) is indistinguishable downstream from a true biological zero, which is
#       precisely the quantity this project measures (antigen-negative cells).
#
# We therefore restrict to the INTERSECTION across retained samples. Every gene
# in the final object was measurable in every cell, so a zero means "not
# detected" rather than "not in this sample's reference".
# ==============================================================================

# genes.tsv in this archive is a SINGLE column of gene symbols -- not the usual
# 2-3 column Ensembl+symbol CellRanger features file (verified across all
# samples). Callers must not assume column 2 exists.
read_gene_symbols <- function(path) {
  utils::read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                    quote = "", comment.char = "")[[1]]
}

# ------------------------------------------------------------------------------
# build_gene_space() -- intersect gene sets across samples and assert that
# everything steps 05-07 depend on survived.
#
# Returns list(common_genes, ref_sizes, n_by_ref).
# ------------------------------------------------------------------------------
build_gene_space <- function(mf, required = REQUIRED_GENES) {
  log_msg("Reconciling gene sets across ", nrow(mf), " samples ...")

  gene_lists <- lapply(mf$genefeat_path, read_gene_symbols)
  names(gene_lists) <- mf$sample_id

  ref_sizes <- vapply(gene_lists, length, integer(1))
  n_by_ref  <- table(ref_sizes)
  log_msg("  reference sizes present: ",
          paste(sprintf("%s genes (n=%d)", names(n_by_ref), as.integer(n_by_ref)),
                collapse = "; "))

  common_genes <- Reduce(intersect, gene_lists)
  log_msg("  intersection: ", length(common_genes), " genes",
          " (largest reference: ", max(ref_sizes), ")")

  if (length(common_genes) == 0) {
    stop("Gene-set intersection is empty -- check that genes.tsv files parse ",
         "as a single symbol column.")
  }

  assert_required_genes(common_genes, required, gene_lists)

  rm(gene_lists)
  invisible(gc())

  list(common_genes = common_genes,
       ref_sizes    = ref_sizes,
       n_by_ref     = n_by_ref)
}

# ------------------------------------------------------------------------------
# assert_required_genes() -- hard stop if a gene the downstream analysis depends
# on is absent, naming the culprit sample(s) so the fix is obvious.
# ------------------------------------------------------------------------------
assert_required_genes <- function(common_genes, required, gene_lists = NULL) {
  missing <- setdiff(required, common_genes)
  if (!length(missing)) {
    log_msg("  all ", length(required), " required marker/antigen genes present")
    return(invisible(TRUE))
  }

  detail <- ""
  if (!is.null(gene_lists)) {
    culprits <- vapply(missing, function(g) {
      bad <- names(gene_lists)[!vapply(gene_lists, function(gl) g %in% gl, logical(1))]
      paste0(g, " absent from: ", paste(utils::head(bad, 5), collapse = ", "),
             if (length(bad) > 5) sprintf(" (+%d more)", length(bad) - 5) else "")
    }, character(1))
    detail <- paste0("\n  ", paste(culprits, collapse = "\n  "))
  }

  stop("Genes required for steps 05-07 are absent from the cross-sample ",
       "intersection: ", paste(missing, collapse = ", "), detail,
       "\nThe integrated object could not answer this project's question. ",
       "Add the offending sample(s) to EXCLUDE_SAMPLES in lib/00_config.R, ",
       "with a recorded reason.")
}
