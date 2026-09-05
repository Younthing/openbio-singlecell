# THROWAWAY: consume a cells-by-genes matrix; return only this operation's results.
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(jsonlite))

manifest <- fromJSON(commandArgs(trailingOnly = TRUE)[1])
expression <- readMM(manifest$matrix)
stopifnot(identical(dim(expression), c(length(manifest$obs_names), length(manifest$var_names))))
rownames(expression) <- manifest$obs_names
colnames(expression) <- manifest$var_names

# Densification is intentional for this tiny prcomp probe, not a production strategy.
fit <- prcomp(as.matrix(expression), center = TRUE, scale. = FALSE, rank. = 2L)
result <- list(
  obs_names = unname(as.list(rownames(expression))),
  totals = unname(as.list(as.numeric(rowSums(expression)))),
  scores = unname(lapply(seq_len(nrow(fit$x)), function(i) as.list(unname(fit$x[i, ])))),
  sdev = unname(as.list(fit$sdev)),
  r_pid = Sys.getpid(),
  versions = list(R = as.character(getRversion()), Matrix = as.character(packageVersion("Matrix")),
                  jsonlite = as.character(packageVersion("jsonlite")))
)
write_json(result, manifest$result, auto_unbox = TRUE, digits = NA, pretty = TRUE)
saveRDS(fit, manifest$model)
restored <- readRDS(manifest$model)
stopifnot(identical(restored$x, fit$x))
png(manifest$plot, width = 720, height = 480)
plot(fit$x[, 1], fit$x[, 2], xlab = "R PC1", ylab = "R PC2", pch = 19,
     col = "#246B8E", main = "Real R PCA through the existing file-artifact workflow")
dev.off()
cat("R PCA finished for", nrow(expression), "cells; PID", Sys.getpid(), "\n")
