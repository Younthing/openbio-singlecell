options(warn = 1)
suppressPackageStartupMessages(library(monocle))
args <- commandArgs(trailingOnly = TRUE)
cds <- readRDS(args[[1]])
versions <- vapply(c("monocle", "DDRTree", "igraph", "dplyr", "Matrix"), function(pkg) as.character(packageVersion(pkg)), character(1))
print(versions)
probe <- function(label, expression) {
  cat("\n---", label, "---\n")
  tryCatch({
    force(expression)
    cat("SUCCEEDED\n")
  }, error = function(error) {
    cat("FAILED:", conditionMessage(error), "\n")
  })
}
probe("estimateDispersions", estimateDispersions(cds))
probe("orderCells", orderCells(cds))
probe("BEAM", BEAM(cds[1:12, ], branch_point = 1, progenitor_method = "duplicate", cores = 1))
