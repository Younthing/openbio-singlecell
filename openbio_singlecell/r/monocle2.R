# One native Monocle2 operation per process. Python owns the complete AnnData.
suppressPackageStartupMessages(library(monocle))
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[[1]])
source(file.path(dirname(normalizePath(script_path)), "monocle2_compat.R"), local = TRUE)
backend <- monocle2_backend()

read_frame <- function(path) {
  metadata <- jsonlite::read_json(path, simplifyVector = FALSE)
  result <- data.frame(row.names = unlist(metadata$index, use.names = FALSE))
  for (name in names(metadata$columns)) {
    column <- metadata$columns[[name]]
    values <- unlist(lapply(column$values, function(value) {
      if (is.null(value)) NA else value
    }), use.names = FALSE)
    result[[name]] <- switch(column$type,
      category = if (is.null(column$levels)) factor(values, ordered = isTRUE(column$ordered)) else
        factor(values, levels = unlist(column$levels, use.names = FALSE),
               ordered = isTRUE(column$ordered)),
      integer = as.numeric(values),
      number = as.numeric(values),
      boolean = as.logical(values),
      string = as.character(values),
      stop(paste("Unknown metadata column type:", column$type)))
  }
  result
}

request <- jsonlite::read_json(commandArgs(trailingOnly = TRUE)[[1]], simplifyVector = TRUE)
parameters <- request$parameters
out <- function(name) file.path(request$output_dir, name)
write_frame <- function(frame, csv_name) {
  original_names <- colnames(frame)
  colnames(frame) <- make.unique(original_names)
  columns <- lapply(seq_along(frame), function(i) {
    column <- frame[[i]]
    type <- if (is.factor(column)) "category" else if (is.logical(column)) "boolean" else
      if (is.integer(column)) "integer" else if (is.numeric(column)) "number" else "string"
    values <- if (type %in% c("category", "string")) as.character(column) else column
    values <- unname(lapply(values, function(value) {
      if (is.na(value)) NULL else if (is.numeric(value) && is.infinite(value)) {
        if (value > 0) "Inf" else "-Inf"
      } else value
    }))
    record <- list(type = type, values = values)
    if (type == "category") {
      record$levels <- as.list(levels(column))
      record$ordered <- is.ordered(column)
    }
    if (colnames(frame)[i] != original_names[i]) record$original_name <- original_names[i]
    record
  })
  names(columns) <- colnames(frame)
  jsonlite::write_json(list(index = as.list(rownames(frame)), columns = columns),
    out(sub("\\.csv$", ".json", csv_name)), auto_unbox = TRUE, null = "null", na = "null", digits = NA)
  utils::write.csv(frame, out(csv_name), row.names = FALSE, na = "")
}
warnings <- character()
effective_parameters <- parameters
native_parameters <- function(exclude = character()) {
  arguments <- parameters[setdiff(names(parameters), c("extra_parameters", exclude))]
  if (!is.null(parameters$extra_parameters)) {
    arguments <- utils::modifyList(arguments, parameters$extra_parameters, keep.null = TRUE)
  }
  if (any(names(arguments) %in% c("cds", "cds_subset", "object"))) {
    stop("The input CellDataSet is supplied by the artifact, not by parameters.")
  }
  effective_parameters <<- arguments
  arguments
}
save_plot <- function(plot) {
  width <- if (is.null(parameters$width)) 8 else parameters$width
  height <- if (is.null(parameters$height)) 6 else parameters$height
  dpi <- if (is.null(parameters$dpi)) 150 else parameters$dpi
  grDevices::png(out("plot.png"), width = width, height = height, units = "in", res = dpi)
  tryCatch(print(plot), finally = grDevices::dev.off())
}

withCallingHandlers({
  cds <- switch(request$operation,
    create = {
      expression <- methods::as(Matrix::readMM(request$matrix_path), "CsparseMatrix")
      obs <- read_frame(request$obs_path)
      var <- read_frame(request$var_path)
      rownames(expression) <- rownames(var)
      colnames(expression) <- rownames(obs)
      family_name <- parameters$expression_family
      if (is.null(family_name)) family_name <- "negbinomial.size"
      families <- c("negbinomial.size", "negbinomial", "tobit", "gaussianff")
      if (!family_name %in% families) stop(paste("Unknown expression family:", family_name))
      limit <- parameters$lowerDetectionLimit
      if (is.null(limit)) limit <- 0.1
      cds <- monocle::newCellDataSet(expression,
        phenoData = methods::new("AnnotatedDataFrame", data = obs),
        featureData = methods::new("AnnotatedDataFrame", data = var),
        expressionFamily = getExportedValue("VGAM", family_name)(),
        lowerDetectionLimit = limit)
      if (isTRUE(parameters$estimate_size_factors)) {
        cds <- do.call(BiocGenerics::estimateSizeFactors,
                       c(list(object = cds), parameters$size_factor_parameters))
      }
      if (isTRUE(parameters$estimate_dispersions)) {
        cds <- do.call(backend$estimateDispersions,
                       c(list(object = cds), parameters$dispersion_parameters))
      }
      if (isTRUE(parameters$detect_genes)) {
        cds <- do.call(monocle::detectGenes, c(list(cds = cds), parameters$detect_parameters))
      }
      cds
    },
    ordering_genes = {
      cds <- readRDS(request$input_rds)
      genes <- switch(parameters$method,
        explicit = parameters$genes,
        var_column = {
          if (!parameters$var_column %in% colnames(Biobase::fData(cds))) {
            stop(paste("Feature annotation column not found:", parameters$var_column))
          }
          rownames(cds)[Biobase::fData(cds)[[parameters$var_column]] %in% TRUE]
        },
        dispersion = {
          dispersions <- monocle::dispersionTable(cds)
          mean_expression <- parameters$mean_expression
          if (is.null(mean_expression)) mean_expression <- 0.5
          dispersion_fold <- parameters$dispersion_fold
          if (is.null(dispersion_fold)) dispersion_fold <- 1
          dispersions$gene_id[dispersions$mean_expression >= mean_expression &
            dispersions$dispersion_empirical >= dispersion_fold * dispersions$dispersion_fit]
        },
        stop(paste("Unknown ordering gene method:", parameters$method)))
      monocle::setOrderingFilter(cds, genes)
    },
    ddrtree = {
      cds <- readRDS(request$input_rds)
      arguments <- native_parameters()
      if (!is.null(arguments$reduction_method) && !identical(arguments$reduction_method, "DDRTree")) {
        stop("Advanced parameters cannot change the DDRTree algorithm: reduction_method must be 'DDRTree'.")
      }
      had_ordering <- length(cds@dim_reduce_type) &&
        length(cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell) > 0
      # Native reduceDimension rebuilds graph projections but leaves old ordering auxiliaries behind.
      cds@auxOrderingData[["DDRTree"]] <- NULL
      cds <- do.call(monocle::reduceDimension, c(list(cds = cds), arguments))
      if (had_ordering) {
        Biobase::pData(cds)$State <- NULL
        Biobase::pData(cds)$Pseudotime <- NULL
        warning("The new DDRTree graph invalidated the previous ordering; State and Pseudotime were removed from this output. Run Order Cells again.")
      }
      cds
    },
    order_cells = {
      cds <- readRDS(request$input_rds)
      do.call(backend$orderCells, c(list(cds = cds), native_parameters()))
    },
    trajectory_plot = {
      cds <- readRDS(request$input_rds)
      plot <- do.call(backend$plot_cell_trajectory,
        c(list(cds = cds), native_parameters(c("width", "height", "dpi"))))
      save_plot(plot)
      cds
    },
    differential_test = {
      cds <- readRDS(request$input_rds)
      test_cds <- if (length(parameters$genes)) cds[parameters$genes, ] else cds
      result <- do.call(monocle::differentialGeneTest,
        c(list(cds = test_cds), native_parameters("genes")))
      result <- data.frame(gene_id = rownames(result), result, check.names = FALSE)
      write_frame(result, "table.csv")
      cds
    },
    gene_trends = {
      cds <- readRDS(request$input_rds)
      plot_cds <- if (length(parameters$genes)) cds[parameters$genes, ] else cds
      plot <- do.call(monocle::plot_genes_in_pseudotime,
        c(list(cds_subset = plot_cds), native_parameters(c("genes", "width", "height", "dpi"))))
      save_plot(plot)
      cds
    },
    beam = {
      cds <- readRDS(request$input_rds)
      test_cds <- if (length(parameters$genes)) cds[parameters$genes, ] else cds
      result <- do.call(backend$BEAM,
        c(list(cds = test_cds), native_parameters("genes")))
      result <- data.frame(gene_id = rownames(result), result, check.names = FALSE)
      write_frame(result, "table.csv")
      cds
    },
    export = {
      cds <- readRDS(request$input_rds)
      if (!length(monocle::reducedDimS(cds)) || !length(cds@dim_reduce_type) ||
          !length(cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell)) {
        stop("Run Monocle 2 Order Cells before export; native ordering results are not available in this CellDataSet.")
      }
      cds
    },
    stop(paste("Unknown Monocle2 operation:", request$operation)))

  if (request$operation %in% c("create", "ordering_genes", "ddrtree", "order_cells")) {
    saveRDS(cds, out("state.rds"))
  }
  cell_data <- data.frame(cell_id = colnames(cds), row.names = colnames(cds), check.names = FALSE)
  for (name in intersect(c("Size_Factor", "num_genes_expressed", "Pseudotime", "State"),
                         colnames(Biobase::pData(cds)))) {
    cell_data[[name]] <- Biobase::pData(cds)[[name]]
  }
  write_frame(cell_data, "cell_metadata.csv")
  gene_data <- data.frame(gene_id = rownames(cds), row.names = rownames(cds), check.names = FALSE)
  for (name in intersect(c("gene_short_name", "use_for_ordering", "num_cells_expressed"),
                         colnames(Biobase::fData(cds)))) {
    gene_data[[name]] <- Biobase::fData(cds)[[name]]
  }
  write_frame(gene_data, "ordering_genes.csv")
  if (length(monocle::reducedDimS(cds))) {
    embedding <- as.data.frame(t(monocle::reducedDimS(cds)))
    colnames(embedding) <- paste0("DDRTree", seq_len(ncol(embedding)))
    embedding <- data.frame(cell_id = rownames(embedding), embedding, check.names = FALSE)
    write_frame(embedding, "embedding.csv")
  }
  summary <- list(operation = request$operation, n_obs = ncol(cds), n_vars = nrow(cds),
                  versions = list(R = as.character(getRversion()),
                    monocle = as.character(utils::packageVersion("monocle")),
                    igraph = as.character(utils::packageVersion("igraph")),
                    dplyr = as.character(utils::packageVersion("dplyr")),
                    DDRTree = as.character(utils::packageVersion("DDRTree"))),
                  compatibility = backend$compatibility,
                  effective_parameters = effective_parameters,
                  ordering_gene_count = sum(Biobase::fData(cds)$use_for_ordering),
                  warnings = as.list(unique(warnings)))
  if ("State" %in% colnames(Biobase::pData(cds))) {
    summary$state_counts <- as.list(table(Biobase::pData(cds)$State))
  }
  if (length(cds@dim_reduce_type)) {
    summary$root_cell <- cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell
  }
  if (request$operation == "order_cells") {
    summary$root_method <- if (!is.null(effective_parameters$root_state)) "state" else
      if (isTRUE(effective_parameters$reverse)) "reverse" else "automatic"
  }
  jsonlite::write_json(summary, out("summary.json"), auto_unbox = TRUE,
                       null = "null", na = "null", digits = NA)
}, warning = function(warning) {
  warnings <<- c(warnings, conditionMessage(warning))
})
