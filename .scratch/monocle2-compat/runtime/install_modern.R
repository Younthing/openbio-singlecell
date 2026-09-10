install.packages(
  c("https://cran.r-project.org/src/contrib/igraph_2.3.3.tar.gz",
    "https://cran.r-project.org/src/contrib/dplyr_1.2.1.tar.gz"),
  repos = NULL,
  type = "source",
  lib = .Library,
  Ncpus = 4L
)
stopifnot(packageVersion("igraph") == "2.3.3", packageVersion("dplyr") == "1.2.1")
