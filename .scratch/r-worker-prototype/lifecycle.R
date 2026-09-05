args <- commandArgs(trailingOnly = TRUE)
mode <- args[[1]]
output_dir <- args[[2]]
evidence_dir <- args[[3]]
delay_seconds <- as.numeric(args[[4]])

writeLines(as.character(Sys.getpid()), file.path(evidence_dir, "r.pid"))
writeLines(R.version.string, file.path(evidence_dir, "r-version.txt"))
writeLines("partial R output", file.path(output_dir, "result.txt"))
writeLines("R has started", file.path(evidence_dir, "started.txt"))

if (mode == "fail") {
  stop("Intentional R lifecycle failure after writing a partial output.")
}
if (mode == "sleep") {
  Sys.sleep(delay_seconds)
  writeLines("R survived cancellation", file.path(evidence_dir, "late.txt"))
}
writeLines("complete R output", file.path(output_dir, "result.txt"))
