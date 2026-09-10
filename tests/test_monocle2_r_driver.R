# Optional native integration test: Rscript --vanilla tests/test_monocle2_r_driver.R
suppressPackageStartupMessages(library(monocle))
script_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
repo <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg))))
driver <- file.path(repo, "openbio_singlecell", "r", "monocle2.R")
work <- tempfile("monocle2-driver-")
dir.create(work)
args <- commandArgs(trailingOnly = TRUE)
phase <- if (length(args)) args[[1]] else "all"
write_json <- function(value, path) {
  jsonlite::write_json(value, path, auto_unbox = TRUE, null = "null", na = "null")
}
run_driver <- function(operation, parameters = list(), input_rds = NULL, ...) {
  out <- tempfile(paste0(operation, "-"), tmpdir = work)
  dir.create(out)
  request <- c(list(operation = operation, parameters = parameters,
                    output_dir = out, input_rds = input_rds), list(...))
  request_path <- paste0(out, ".json")
  write_json(request, request_path)
  result <- system2(file.path(R.home("bin"), "Rscript"),
                    c("--vanilla", shQuote(driver), shQuote(request_path)),
                    stdout = TRUE, stderr = TRUE)
  status <- attr(result, "status")
  if (!is.null(status) && status != 0) stop(paste(result, collapse = "\n"))
  out
}

# Creation preserves the user-selected floating expression values and metadata.
mat <- matrix(seq(0.25, 12, by = 0.25), nrow = 8)
genes <- paste0("gene", seq_len(nrow(mat)))
cells <- paste0("cell", seq_len(ncol(mat)))
rownames(mat) <- genes
colnames(mat) <- cells
matrix_path <- file.path(work, "expression.mtx")
invisible(Matrix::writeMM(methods::as(mat, "sparseMatrix"), matrix_path))
obs_path <- file.path(work, "obs.json")
var_path <- file.path(work, "var.json")
write_json(list(index = cells, columns = list(
  subtype = list(type = "category", values = list("early", "late", NULL, "early", "late", "early"),
                 levels = c("early", "late", "unused"), ordered = TRUE),
  cohort = list(type = "category", values = c("A", "B", "A", "B", "A", "B")),
  State = list(type = "string", values = rep("input_annotation", 6)),
  note = list(type = "string", values = c("α", "b", "c", "d", "e", "f")))), obs_path)
write_json(list(index = genes, columns = list(
  gene_short_name = list(type = "string", values = list("001", "", "NA", NULL, "gene5", "gene6", "gene7", "gene8")),
  highly_variable = list(type = "boolean", values = rep(c(TRUE, FALSE), 4)))), var_path)
created <- run_driver("create", list(estimate_size_factors = FALSE,
  estimate_dispersions = FALSE, detect_genes = FALSE), matrix_path = matrix_path,
  obs_path = obs_path, var_path = var_path)
cds <- readRDS(file.path(created, "state.rds"))
stopifnot(identical(as.matrix(Biobase::exprs(cds)), mat))
stopifnot(identical(levels(Biobase::pData(cds)$subtype), c("early", "late", "unused")))
stopifnot(is.ordered(Biobase::pData(cds)$subtype), is.na(Biobase::pData(cds)$subtype[3]))
stopifnot(Biobase::pData(cds)$note[1] == "α")
stopifnot(identical(Biobase::pData(cds)$State, rep("input_annotation", 6)))
stopifnot(identical(as.character(Biobase::pData(cds)$cohort), c("A", "B", "A", "B", "A", "B")))
summary <- jsonlite::read_json(file.path(created, "summary.json"), simplifyVector = TRUE)
stopifnot(summary$n_obs == ncol(mat), summary$n_vars == nrow(mat))
typed_genes <- jsonlite::read_json(file.path(created, "ordering_genes.json"), simplifyVector = FALSE)
stopifnot(identical(unlist(typed_genes$index), genes))
stopifnot(typed_genes$columns$gene_short_name$type == "string")
stopifnot(identical(typed_genes$columns$gene_short_name$values[1:4], list("001", "", "NA", NULL)))
typed_cells <- jsonlite::read_json(file.path(created, "cell_metadata.json"), simplifyVector = FALSE)
stopifnot(identical(unlist(typed_cells$index), cells))
stopifnot(typed_cells$columns$Size_Factor$type == "number")
stopifnot(length(typed_cells$columns$Size_Factor$values) == length(cells),
          all(vapply(typed_cells$columns$Size_Factor$values, is.null, logical(1))))
message("PASS create: floating expression and categorical metadata survive native RDS")
if (phase == "create") quit(status = 0)

# An explicit ordering gene list keeps native unknown-ID behavior and all genes.
ordered_genes <- run_driver("ordering_genes", list(method = "explicit",
  genes = c("gene2", "gene5", "not_in_data")), file.path(created, "state.rds"))
selected <- readRDS(file.path(ordered_genes, "state.rds"))
stopifnot(identical(rownames(Biobase::fData(selected))[Biobase::fData(selected)$use_for_ordering],
                    c("gene2", "gene5")))
stopifnot(identical(as.matrix(Biobase::exprs(selected)), mat))
stopifnot(file.exists(file.path(ordered_genes, "ordering_genes.csv")))
typed_selection <- jsonlite::read_json(file.path(ordered_genes, "ordering_genes.json"), simplifyVector = FALSE)
stopifnot(typed_selection$columns$use_for_ordering$type == "boolean",
          identical(unlist(typed_selection$columns$use_for_ordering$values), c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE)))
message("PASS ordering genes: explicit selections preserve the complete expression matrix")
if (phase == "ordering_genes") quit(status = 0)

# Existing highly-variable annotations are a user-owned ordering selection.
var_genes <- run_driver("ordering_genes", list(method = "var_column",
  var_column = "highly_variable"), file.path(created, "state.rds"))
var_selected <- readRDS(file.path(var_genes, "state.rds"))
stopifnot(identical(which(Biobase::fData(var_selected)$use_for_ordering), c(1L, 3L, 5L, 7L)))
missing_column <- tryCatch({
  suppressWarnings(run_driver("ordering_genes", list(method = "var_column", var_column = "missing"),
                              file.path(created, "state.rds")))
  NULL
}, error = conditionMessage)
stopifnot(is.character(missing_column), grepl("Feature annotation column not found", missing_column))
message("PASS ordering var column: existing feature annotations select ordering genes")
if (phase == "var_column") quit(status = 0)

# A reproducible branching population supplies real native fits for later stages.
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
counts <- matrix(rnbinom(g * n, mu = as.vector(means), size = 5), nrow = g)
rownames(counts) <- paste0("gene_", seq_len(g))
colnames(counts) <- paste0("cell_", seq_len(n))
invisible(Matrix::writeMM(methods::as(counts, "sparseMatrix"), matrix_path))
write_json(list(index = colnames(counts), columns = list(
  branch = list(type = "string", values = branch),
  State = list(type = "string", values = rep("old_annotation", n)),
  Pseudotime = list(type = "number", values = rep(42, n)),
  true_time = list(type = "number", values = latent[, "progression"]))), obs_path)
write_json(list(index = rownames(counts), columns = list(
  gene_short_name = list(type = "string", values = rownames(counts)))), var_path)
prepared <- run_driver("create", list(estimate_size_factors = TRUE,
  estimate_dispersions = TRUE, detect_genes = TRUE), matrix_path = matrix_path,
  obs_path = obs_path, var_path = var_path)
prepared_cds <- readRDS(file.path(prepared, "state.rds"))
stopifnot(all(is.finite(BiocGenerics::sizeFactors(prepared_cds))))
stopifnot(nrow(monocle::dispersionTable(prepared_cds)) > 0)

# Dispersion thresholds remain user-selectable, including permissive thresholds.
dispersion_genes <- run_driver("ordering_genes", list(method = "dispersion",
  mean_expression = 0, dispersion_fold = 0), file.path(prepared, "state.rds"))
dispersion_cds <- readRDS(file.path(dispersion_genes, "state.rds"))
stopifnot(sum(Biobase::fData(dispersion_cds)$use_for_ordering) == g)
message("PASS dispersion: native fit survives RDS and user thresholds select all genes")
if (phase == "dispersion") quit(status = 0)

# DDRTree executes in a fresh R process and returns coordinates plus native state.
changed_algorithm <- tryCatch({
  suppressWarnings(run_driver("ddrtree", list(extra_parameters = list(reduction_method = "ICA")),
                               file.path(dispersion_genes, "state.rds")))
  NULL
}, error = conditionMessage)
stopifnot(is.character(changed_algorithm), grepl("cannot change the DDRTree algorithm", changed_algorithm))
reduced <- run_driver("ddrtree", list(max_components = 3,
  extra_parameters = list(max_components = 2, reduction_method = "DDRTree")),
  file.path(dispersion_genes, "state.rds"))
reduced_cds <- readRDS(file.path(reduced, "state.rds"))
embedding <- read.csv(file.path(reduced, "embedding.csv"), check.names = FALSE)
stopifnot(identical(dim(monocle::reducedDimS(reduced_cds)), c(2L, n)))
stopifnot(identical(embedding$cell_id, colnames(counts)), ncol(embedding) == 3)
typed_embedding <- jsonlite::read_json(file.path(reduced, "embedding.json"), simplifyVector = FALSE)
stopifnot(identical(unlist(typed_embedding$index), colnames(counts)))
stopifnot(typed_embedding$columns$DDRTree1$type == "number", is.list(typed_embedding$columns$DDRTree1$values))
stopifnot(igraph::vcount(monocle::minSpanningTree(reduced_cds)) > 0)
summary <- jsonlite::read_json(file.path(reduced, "summary.json"), simplifyVector = TRUE)
stopifnot(summary$effective_parameters$max_components == 2)
message("PASS DDRTree: native tree and cell coordinates survive subprocess RDS")
if (phase == "ddrtree") quit(status = 0)

# Input columns with familiar names are not newly computed Monocle pseudotime.
unordered_export <- tryCatch({
  suppressWarnings(run_driver("export", input_rds = file.path(reduced, "state.rds")))
  NULL
}, error = conditionMessage)
stopifnot(is.character(unordered_export), grepl("Run Monocle 2 Order Cells before export", unordered_export))
message("PASS export precondition: preexisting State/Pseudotime do not masquerade as an ordering")
if (phase == "export_precondition") quit(status = 0)

# Initial native ordering is exploratory; a later explicit State chooses direction.
initial_order <- run_driver("order_cells", list(), file.path(reduced, "state.rds"))
initial_cds <- readRDS(file.path(initial_order, "state.rds"))
stopifnot(all(is.finite(Biobase::pData(initial_cds)$Pseudotime)))
stopifnot(length(unique(Biobase::pData(initial_cds)$State)) > 1)
initial_export <- run_driver("export", input_rds = file.path(initial_order, "state.rds"))
stopifnot(file.exists(file.path(initial_export, "cell_metadata.csv")))
early_states <- table(Biobase::pData(initial_cds)$State[branch == "progenitor"])
root_state <- names(early_states)[which.max(early_states)]
rooted <- run_driver("order_cells", list(root_state = root_state),
                     file.path(initial_order, "state.rds"))
rooted_cds <- readRDS(file.path(rooted, "state.rds"))
metadata <- read.csv(file.path(rooted, "cell_metadata.csv"), check.names = FALSE)
stopifnot(identical(metadata$cell_id, colnames(counts)))
stopifnot(as.character(metadata$State[which.min(metadata$Pseudotime)]) == root_state)
summary <- jsonlite::read_json(file.path(rooted, "summary.json"), simplifyVector = TRUE)
stopifnot(sum(unlist(summary$state_counts)) == n, summary$root_method == "state")
message("PASS order cells: default exploration and explicit root State survive subprocesses")
if (phase == "order_cells") quit(status = 0)

# Recomputing a graph invalidates ordering from the old graph in the new output.
rooted_path <- file.path(rooted, "state.rds")
rooted_hash <- unname(tools::md5sum(rooted_path))
recomputed <- run_driver("ddrtree", list(max_components = 2), rooted_path)
recomputed_cds <- readRDS(file.path(recomputed, "state.rds"))
recomputed_export <- tryCatch({
  suppressWarnings(run_driver("export", input_rds = file.path(recomputed, "state.rds")))
  NULL
}, error = conditionMessage)
stopifnot(is.character(recomputed_export), grepl("Run Monocle 2 Order Cells before export", recomputed_export))
stopifnot(!any(c("State", "Pseudotime") %in% colnames(Biobase::pData(recomputed_cds))))
stopifnot(!length(recomputed_cds@auxOrderingData[["DDRTree"]]$root_cell))
stopifnot(unname(tools::md5sum(rooted_path)) == rooted_hash)
recomputed_summary <- jsonlite::read_json(file.path(recomputed, "summary.json"), simplifyVector = TRUE)
stopifnot(any(grepl("previous ordering", recomputed_summary$warnings)))
reordered <- run_driver("order_cells", list(), file.path(recomputed, "state.rds"))
reordered_export <- run_driver("export", input_rds = file.path(reordered, "state.rds"))
stopifnot(file.exists(file.path(reordered_export, "cell_metadata.csv")))
message("PASS graph recomputation: stale ordering is removed, new native ordering restores export")
if (phase == "recompute") quit(status = 0)

# Native trajectory plotting runs after loading the RDS in another process.
trajectory <- run_driver("trajectory_plot", list(color_by = "branch",
  show_state_number = TRUE, width = 5, height = 4, dpi = 100),
  file.path(rooted, "state.rds"))
stopifnot(file.info(file.path(trajectory, "plot.png"))$size > 1000)
stopifnot(identical(readBin(file.path(trajectory, "plot.png"), "raw", 8),
                    as.raw(c(137, 80, 78, 71, 13, 10, 26, 10))))
stopifnot(!file.exists(file.path(trajectory, "state.rds")))
message("PASS trajectory plot: native Monocle figure renders from the saved state")
if (phase == "trajectory_plot") quit(status = 0)

# A selected set of genes is tested without silently dropping failed fits.
tested_genes <- rownames(counts)[1:12]
typed_cds <- readRDS(file.path(rooted, "state.rds"))
Biobase::fData(typed_cds)$numeric_annotation <- c(0.125, NA_real_, Inf, -Inf, rep(1.5, g - 4))
Biobase::fData(typed_cds)$integer_annotation <- c(1L, NA_integer_, rep(2L, g - 2))
Biobase::fData(typed_cds)$label <- c("001", "", "NA", NA_character_, rep("label", g - 4))
Biobase::fData(typed_cds)$category <- factor(c("001", NA_character_, rep("later", g - 2)),
  levels = c("001", "later", "unused"), ordered = TRUE)
Biobase::fData(typed_cds)$flag <- c(TRUE, NA, FALSE, rep(TRUE, g - 3))
Biobase::fData(typed_cds)$gene_id <- paste0("external_", seq_len(g))
typed_state <- file.path(work, "typed-cds.rds")
saveRDS(typed_cds, typed_state)
de <- run_driver("differential_test", list(genes = tested_genes,
  fullModelFormulaStr = "~sm.ns(Pseudotime, df=3)", cores = 1),
  typed_state)
de_table <- read.csv(file.path(de, "table.csv"), check.names = FALSE)
stopifnot(setequal(de_table$gene_id, tested_genes), nrow(de_table) == 12)
stopifnot(all(c("status", "pval", "qval") %in% names(de_table)))
stopifnot(any(de_table$status == "OK"))
stopifnot(!file.exists(file.path(de, "state.rds")))
typed_table <- jsonlite::read_json(file.path(de, "table.json"), simplifyVector = FALSE)
stopifnot(identical(unlist(typed_table$index), tested_genes))
stopifnot(typed_table$columns$numeric_annotation$type == "number",
  identical(typed_table$columns$numeric_annotation$values[1:4], list(0.125, NULL, "Inf", "-Inf")))
stopifnot(typed_table$columns$integer_annotation$type == "integer",
  identical(typed_table$columns$integer_annotation$values[1:2], list(1L, NULL)))
stopifnot(typed_table$columns$label$type == "string",
  identical(typed_table$columns$label$values[1:4], list("001", "", "NA", NULL)))
stopifnot(typed_table$columns$category$type == "category", isTRUE(typed_table$columns$category$ordered),
  identical(typed_table$columns$category$levels, list("001", "later", "unused")))
stopifnot(typed_table$columns$flag$type == "boolean",
  identical(typed_table$columns$flag$values[1:3], list(TRUE, NULL, FALSE)))
stopifnot(typed_table$columns$gene_id.1$original_name == "gene_id",
  typed_table$columns$gene_id.1$values[[1]] == "external_1", typed_table$columns$gene_id$values[[1]] == "gene_1")
single_gene <- run_driver("differential_test", list(genes = tested_genes[1], cores = 1), typed_state)
typed_single <- jsonlite::read_json(file.path(single_gene, "table.json"), simplifyVector = FALSE)
stopifnot(is.list(typed_single$index), length(typed_single$index) == 1,
  is.list(typed_single$columns$gene_id$values), length(typed_single$columns$gene_id$values) == 1)
message("PASS differential test: complete native result rows preserve gene IDs and fit status")
if (phase == "differential_test") quit(status = 0)

# Marker trends use the saved size-factor and dispersion fit for native smoothing.
trends <- run_driver("gene_trends", list(genes = tested_genes[1:3],
  color_by = "branch", ncol = 3, width = 8, height = 3), file.path(rooted, "state.rds"))
stopifnot(file.info(file.path(trends, "plot.png"))$size > 1000)
stopifnot(!file.exists(file.path(trends, "state.rds")))
message("PASS gene trends: native pseudotime smoothing renders from saved dispersion fits")
if (phase == "gene_trends") quit(status = 0)

# BEAM uses the native branch graph and leaves statistical failures in its table.
beam <- run_driver("beam", list(genes = tested_genes,
  branch_point = 1, progenitor_method = "duplicate", cores = 1),
  file.path(rooted, "state.rds"))
beam_table <- read.csv(file.path(beam, "table.csv"), check.names = FALSE)
stopifnot(setequal(beam_table$gene_id, tested_genes), nrow(beam_table) == 12)
stopifnot(all(c("status", "pval", "qval") %in% names(beam_table)))
stopifnot(any(beam_table$status == "OK"))
message("PASS BEAM: native branch-dependent test executes from the saved graph")
if (phase == "beam") quit(status = 0)

# Export exposes aligned additions for Python while leaving the input immutable.
source_state <- file.path(rooted, "state.rds")
source_digest <- unname(tools::md5sum(source_state))
exported <- run_driver("export", input_rds = source_state)
exported_cells <- read.csv(file.path(exported, "cell_metadata.csv"), check.names = FALSE)
exported_embedding <- read.csv(file.path(exported, "embedding.csv"), check.names = FALSE)
stopifnot(identical(exported_cells$cell_id, colnames(counts)))
stopifnot(identical(exported_embedding$cell_id, colnames(counts)))
stopifnot(isTRUE(all.equal(exported_cells$Pseudotime, Biobase::pData(rooted_cds)$Pseudotime)))
stopifnot(unname(tools::md5sum(source_state)) == source_digest)
stopifnot(!file.exists(file.path(exported, "state.rds")))
message("PASS export: cell additions align by ID and the source artifact stays immutable")
