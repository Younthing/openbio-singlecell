request <- jsonlite::fromJSON(commandArgs(trailingOnly = TRUE)[[1]], simplifyVector = FALSE)
for (package in unlist(request$packages)) {
    loadNamespace(package)
}
packages <- lapply(sort(loadedNamespaces()), function(package) {
    path <- normalizePath(find.package(package), winslash = "/", mustWork = TRUE)
    list(version = as.character(utils::packageVersion(package)), path = path)
})
names(packages) <- sort(loadedNamespaces())
jsonlite::write_json(
    list(r_version = R.version.string, packages = packages),
    file.path(request$output_dir, "summary.json"), auto_unbox = TRUE, pretty = TRUE
)
