project_library <- Sys.getenv("NEFF_RLIB", "Rlib")
dependency_library <- Sys.getenv("NEFF_DEPENDENCY_RLIB", "")
if (nzchar(dependency_library)) {
  dependency_paths <- strsplit(dependency_library, .Platform$path.sep, fixed = TRUE)[[1]]
  .libPaths(c(dependency_paths, .libPaths()))
}
dir.create(project_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(project_library, .libPaths()))
install.packages("provenance/ontophylo", repos = NULL, type = "source", lib = project_library)
source("scripts/load-packages.R")
capture.output(sessionInfo(), file = "provenance/installed-session.txt")
packages <- c("ontophylo", "ape", "diversitree", "corHMM", "phytools", "digest", "R.utils")
versions <- data.frame(package = packages, version = NA_character_)
for (i in seq_along(packages)) versions$version[i] <- as.character(packageVersion(packages[i]))
write.csv(versions, "provenance/package-versions.csv", row.names = FALSE)
print(versions)
