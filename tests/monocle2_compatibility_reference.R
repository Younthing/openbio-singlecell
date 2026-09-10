# Independent native Monocle 2 reference and read-only RDS inspection for parity tests.
suppressPackageStartupMessages(library(monocle))
args <- commandArgs(trailingOnly = TRUE)
mode <- args[[1]]

if (mode == "reference") {
  input <- args[[2]]
  output <- args[[3]]
  expression <- methods::as(Matrix::readMM(file.path(input, "expression.mtx")), "CsparseMatrix")
  obs <- read.csv(file.path(input, "obs.csv"), row.names = 1, check.names = FALSE)
  var <- read.csv(file.path(input, "var.csv"), row.names = 1, check.names = FALSE)
  rownames(expression) <- rownames(var)
  colnames(expression) <- rownames(obs)
  cds <- monocle::newCellDataSet(expression,
    phenoData = methods::new("AnnotatedDataFrame", data = obs),
    featureData = methods::new("AnnotatedDataFrame", data = var),
    expressionFamily = VGAM::negbinomial.size(), lowerDetectionLimit = 0.1)
  cds <- BiocGenerics::estimateSizeFactors(cds)
  cds <- BiocGenerics::estimateDispersions(cds, modelFormulaStr = "~branch", remove_outliers = FALSE)
  cds <- monocle::detectGenes(cds)
  saveRDS(cds, file.path(output, "create.rds"))
  cds <- monocle::setOrderingFilter(cds, rownames(cds))
  cds <- monocle::reduceDimension(cds, reduction_method = "DDRTree", max_components = 2)
  saveRDS(cds, file.path(output, "ddrtree.rds"))
  cds <- monocle::orderCells(cds, reverse = FALSE)
  saveRDS(cds, file.path(output, "order_cells.rds"))
  states <- table(Biobase::pData(cds)$State[obs$branch == "progenitor"])
  root_state <- names(states)[which.max(states)]
  cds <- monocle::orderCells(cds, root_state = root_state, reverse = FALSE)
  saveRDS(cds, file.path(output, "reroot.rds"))
  writeLines(root_state, file.path(output, "root_state.txt"))
  selected <- cds[rownames(cds)[1:8], ]
  de <- monocle::differentialGeneTest(selected,
    fullModelFormulaStr = "~sm.ns(Pseudotime, df=3) + branch",
    reducedModelFormulaStr = "~branch", relative_expr = FALSE, cores = 1)
  write.csv(de, file.path(output, "differential_test.csv"))
  beam <- monocle::BEAM(selected, branch_point = 1, progenitor_method = "duplicate", cores = 1)
  write.csv(beam, file.path(output, "beam.csv"))
} else if (mode == "inspect") {
  # Loading each state in a new process also checks that the adapter has not
  # serialized closures that depend on its worker's temporary environment.
  cds <- readRDS(args[[2]])
  graph_snapshot <- function(graph) {
    edges <- igraph::as_data_frame(graph, what = "edges")
    edges <- edges[order(edges$from, edges$to), , drop = FALSE]
    list(vertices = igraph::V(graph)$name, edges = unname(as.matrix(edges[, c("from", "to")])),
      weights = edges$weight)
  }
  snapshot <- list(
    cells = colnames(cds), genes = rownames(cds),
    size_factors = unname(BiocGenerics::sizeFactors(cds)),
    dispersion = monocle::dispersionTable(cds),
    dispersion_coefficients = unname(attr(cds@dispFitInfo[["blind"]]$disp_func, "coefficients")))
  if (length(monocle::reducedDimS(cds))) {
    snapshot$coordinates <- unname(monocle::reducedDimS(cds))
    snapshot$principal_points <- unname(monocle::reducedDimK(cds))
    snapshot$graph <- graph_snapshot(monocle::minSpanningTree(cds))
    auxiliary <- cds@auxOrderingData[["DDRTree"]]
    if (length(auxiliary$root_cell)) {
      snapshot$root <- auxiliary$root_cell
      snapshot$state <- as.character(Biobase::pData(cds)$State)
      snapshot$pseudotime <- Biobase::pData(cds)$Pseudotime
      snapshot$branch_points <- auxiliary$branch_points
      snapshot$projection <- unname(auxiliary$pr_graph_cell_proj_dist)
      snapshot$closest_vertex <- unname(auxiliary$pr_graph_cell_proj_closest_vertex)
      graph <- auxiliary$pr_graph_cell_proj_tree
      snapshot$projected_graph <- graph_snapshot(graph)
      # DFS parents are observable topology, independent of the Monocle adapter.
      traversal <- igraph::dfs(graph, root = auxiliary$root_cell, mode = "all",
        unreachable = FALSE, father = TRUE)
      parent <- if ("parent" %in% names(traversal)) traversal$parent else traversal$father
      snapshot$parents <- unname(as.integer(parent))
    }
  }
  jsonlite::write_json(snapshot, args[[3]], auto_unbox = TRUE, digits = NA, na = "null", null = "null")
} else {
  stop(paste("Unknown reference mode:", mode))
}
