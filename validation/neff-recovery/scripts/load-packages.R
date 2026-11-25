# Install the frozen source into Rlib before running this workflow.
project_library <- Sys.getenv("NEFF_RLIB", "Rlib")
dependency_library <- Sys.getenv("NEFF_DEPENDENCY_RLIB", "")
if (nzchar(dependency_library)) {
  dependency_paths <- strsplit(dependency_library, .Platform$path.sep, fixed = TRUE)[[1]]
  .libPaths(c(dependency_paths, .libPaths()))
}
.libPaths(c(project_library, .libPaths()))
library(ape)
library(diversitree)
library(ontophylo)
stopifnot(all(c("fit_models", "paramo_stm", "branch_rate_pipeline", "phenomic_metrics") %in%
              getNamespaceExports("ontophylo")))
stopifnot(normalizePath(find.package("ontophylo")) ==
          normalizePath(file.path(project_library, "ontophylo")))
