# Optional native check: Rscript --vanilla tests/test_monocle2_compat.R
suppressPackageStartupMessages(library(monocle))
script_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
repo <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg))))
source(file.path(repo, "openbio_singlecell", "r", "monocle2_compat.R"))

# The same branching fixture as the R driver exercises a real dispersion fit.
set.seed(2026)
n <- 180L
g <- 400L
progression <- rep(seq(0, 1, length.out = n / 3), 3)
branch <- rep(c("progenitor", "branch_a", "branch_b"), each = n / 3)
latent <- cbind(progression = progression + rep(c(0, 1, 1), each = n / 3),
  branch_a = as.numeric(branch == "branch_a") * progression,
  branch_b = as.numeric(branch == "branch_b") * progression)
loadings <- matrix(rnorm(g * 3, sd = 0.8), ncol = 3)
means <- exp(1.5 + loadings %*% t(latent))
counts <- matrix(rnbinom(g * n, mu = as.vector(means), size = 5), nrow = g,
                 dimnames = list(paste0("gene_", seq_len(g)), paste0("cell_", seq_len(n))))
cds <- monocle::newCellDataSet(methods::as(counts, "sparseMatrix"),
  phenoData = methods::new("AnnotatedDataFrame", data = data.frame(
    branch = branch, row.names = colnames(counts))),
  featureData = methods::new("AnnotatedDataFrame", data = data.frame(
    gene_short_name = rownames(counts), row.names = rownames(counts))))
cds <- BiocGenerics::estimateSizeFactors(cds)
native_dispersion <- monocle:::estimateDispersionsForCellDataSet
backend <- monocle2_backend()
prepared <- withCallingHandlers(backend$estimateDispersions(cds), warning = function(w) {
  if (grepl("group_by_|select_", conditionMessage(w))) stop(w)
})
stopifnot(nrow(monocle::dispersionTable(prepared)) == g)
stopifnot(identical(monocle:::estimateDispersionsForCellDataSet, native_dispersion))
restored <- unserialize(serialize(prepared, NULL))
stopifnot(identical(monocle::dispersionTable(restored), monocle::dispersionTable(prepared)))
stopifnot(length(serialize(restored@dispFitInfo$blind$disp_func, NULL)) < 4096)
stopifnot(length(serialize(restored@dispFitInfo, NULL)) < 100000)
message("PASS dispersion: current dplyr API, native namespace isolation, compact saved fit")

# Native formula parsing accepts backticks around user annotation names.
Biobase::pData(cds)[["branch label"]] <- branch
grouped <- backend$estimateDispersions(cds, modelFormulaStr = "~`branch label`",
                                       remove_outliers = FALSE)
stopifnot(nrow(grouped@dispFitInfo$blind$disp_table) == g * 3)
message("PASS grouped dispersion: user formula syntax preserves quoted annotation names")

native_order <- monocle::orderCells
reduced <- monocle::reduceDimension(monocle::setOrderingFilter(prepared, rownames(prepared)),
                                    max_components = 2, reduction_method = "DDRTree")
ordered <- backend$orderCells(reduced)
stopifnot(all(is.finite(Biobase::pData(ordered)$Pseudotime)))
stopifnot(length(unique(Biobase::pData(ordered)$State)) > 1)
early_states <- table(Biobase::pData(ordered)$State[branch == "progenitor"])
root_state <- names(early_states)[which.max(early_states)]
rooted <- backend$orderCells(ordered, root_state = root_state)
stopifnot(as.character(Biobase::pData(rooted)$State[which.min(Biobase::pData(rooted)$Pseudotime)]) == root_state)
stopifnot(identical(monocle::orderCells, native_order))
stopifnot(length(serialize(rooted@auxOrderingData, NULL)) < 1000000)
message("PASS ordering: finite pseudotime and explicit root State using private igraph adapters")

native_beam <- monocle::BEAM
tested_genes <- rownames(rooted)[1:8]
beam <- backend$BEAM(rooted[tested_genes, ], branch_point = 1,
                     progenitor_method = "duplicate", cores = 1)
stopifnot(identical(rownames(beam), tested_genes), any(beam$status == "OK"))
stopifnot(identical(monocle::BEAM, native_beam))
message("PASS BEAM: branch-dependent fits retain selected genes and native fit status")
