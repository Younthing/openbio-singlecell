# THROWAWAY: isolate native package-load failure from H5AD behavior.
p <- commandArgs(trailingOnly = TRUE)[[1]]
library(p, character.only = TRUE)
cat(p, as.character(packageVersion(p)), "OK\n")
