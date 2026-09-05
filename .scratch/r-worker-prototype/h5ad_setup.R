# THROWAWAY PROTOTYPE: install only binary runtime dependencies in disposable R.
options(timeout = 180, repos = c(
  BioC = "https://bioconductor.org/packages/3.22/bioc",
  CRAN = "https://cloud.r-project.org"
))
stopifnot(grepl("openbio-r-prototype", .libPaths()[1], fixed = TRUE))
install.packages(c("anndataR", "rhdf5"), type = "win.binary", dependencies = NA)
stopifnot(requireNamespace("anndataR"), requireNamespace("rhdf5"))
print(installed.packages()[, c("Package", "Version")])
