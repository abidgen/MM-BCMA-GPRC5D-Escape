options(repos = c(CRAN = "https://cloud.r-project.org"))

required_nmf <- package_version("0.23.0")
installed_nmf <- tryCatch(packageVersion("NMF"), error = function(e) package_version("0.0.0"))

if (installed_nmf < required_nmf) {
  message("Installing/upgrading NMF from CRAN...")
  install.packages("NMF")
}

stopifnot(packageVersion("NMF") >= required_nmf)
message("NMF version: ", packageVersion("NMF"))

packages <- c(
  "Seurat", "harmony", "scDblFinder", "SingleR", "celldex",
  "infercnv", "glmGamPoi", "SingleCellExperiment", "scran", "CellChat"
)

failed <- character()
for (pkg in packages) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  if (ok) {
    message(sprintf("%-25s %s", pkg, as.character(packageVersion(pkg))))
  } else {
    message(sprintf("%-25s FAILED", pkg))
    failed <- c(failed, pkg)
  }
}

if (length(failed)) {
  stop("Packages that failed to load: ", paste(failed, collapse = ", "))
}
