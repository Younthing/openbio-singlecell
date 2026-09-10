# Local API translations for Monocle 2; installed package namespaces stay intact.
monocle2_backend <- function() {
  backend <- new.env(parent = asNamespace("monocle"))
  backend$.select_columns <- function(.data, ..., .dots = character()) {
    # Legacy Monocle passes formula terms as strings. Retain their expression
    # parsing (including backticked names), then use the current dplyr API.
    caller <- parent.frame()
    terms <- lapply(c(as.list(.dots), list(...)), function(term) {
      if (is.character(term)) rlang::parse_quo(term, env = caller) else
        rlang::new_quosure(term, env = caller)
    })
    dplyr::select(.data, !!!terms)
  }
  backend$.group_columns <- function(.data, .dots = character()) {
    terms <- lapply(.dots, rlang::parse_quo, env = parent.frame())
    dplyr::group_by(.data, !!!terms)
  }
  backend$.dfs_parent <- if ("parent" %in% names(formals(igraph::dfs))) "parent" else "father"
  backend$.graph_dfs <- function(graph, root, neimode, unreachable, father) {
    arguments <- list(graph = graph, root = root, mode = neimode, unreachable = unreachable)
    arguments[[.dfs_parent]] <- father
    traversal <- do.call(igraph::dfs, arguments)
    traversal$father <- traversal[[.dfs_parent]]
    traversal
  }
  environment(backend$.graph_dfs) <- backend
  dot_nei <- exists(".nei", envir = asNamespace("igraph"), inherits = FALSE)
  translate_api <- function(expr) {
    if (missing(expr)) return(quote(expr = ))
    if (identical(expr, quote(dplyr::select_)) || identical(expr, quote(select_))) {
      return(quote(.select_columns))
    }
    if (identical(expr, quote(dplyr::group_by_))) return(quote(.group_columns))
    if (identical(expr, quote(graph.dfs))) return(quote(.graph_dfs))
    # .nei is evaluated by igraph's vertex-sequence data mask, not a normal
    # function lookup. A wrapper around nei() cannot cross that boundary.
    if (dot_nei && identical(expr, quote(nei))) return(quote(.nei))
    # orderCells and the dispersion method create state containers. Their
    # parent frames are not data and must not capture this backend in an RDS.
    if (is.call(expr) && identical(expr[[1]], quote(new.env))) expr$parent <- quote(emptyenv())
    if (is.call(expr)) return(as.call(lapply(expr, translate_api)))
    expr
  }
  for (name in c("orderCells", "project2MST", "extract_ddrtree_ordering",
                 "BEAM", "branchTest", "buildBranchCellDataSet", "plot_cell_trajectory")) {
    adapted <- get(name, asNamespace("monocle"))
    body(adapted) <- translate_api(body(adapted))
    environment(adapted) <- backend
    backend[[name]] <- adapted
  }
  fit <- get("estimateDispersionsForCellDataSet", asNamespace("monocle"))
  body(fit) <- translate_api(body(fit))
  environment(fit) <- backend
  backend$estimateDispersionsForCellDataSet <- fit
  # Monocle's selected S4 method embeds its implementation as a .local closure.
  # Calling the generic would dispatch back into the original namespace.
  method <- body(methods::getMethod("estimateDispersions", "CellDataSet"))[[2]][[3]]
  body(method) <- translate_api(body(method))
  environment(method) <- backend
  backend$.estimate_dispersions <- method
  backend$estimateDispersions <- function(object, ...) {
    object <- .estimate_dispersions(object, ...)
    # The native fit closure captures its complete fitting frame. Persist only
    # its coefficients so RDS never captures this adapter or the input matrix.
    fit <- object@dispFitInfo$blind$disp_func
    environment(fit) <- list2env(list(coefs = attr(fit, "coefficients")), parent = baseenv())
    object@dispFitInfo$blind$disp_func <- fit
    object
  }
  environment(backend$estimateDispersions) <- backend
  backend$compatibility <- list(adapter = "monocle2-api-v1",
    adaptations = as.list(c("dplyr: select_/group_by_ -> select/group_by with parsed formula terms",
      paste0("igraph: graph.dfs(neimode, father) -> dfs(mode, ", backend$.dfs_parent, ")"),
      if (dot_nei) "igraph: nei -> .nei")))
  backend
}
