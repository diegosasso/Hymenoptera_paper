# Use the project library, with an optional shared dependency library.
project_library <- Sys.getenv("PHYLO_SAMPLING_RLIB", "Rlib")
dependency_library <- Sys.getenv("PHYLO_SAMPLING_DEPENDENCY_RLIB", "")
if (nzchar(dependency_library)) .libPaths(c(dependency_library, .libPaths()))
.libPaths(c(project_library, .libPaths()))
library(ape)
library(ontophylo)
stopifnot(all(c("fit_models", "paramo_stm", "read_stm", "branch_rate_pipeline",
               "branch_rate_thru_time") %in% getNamespaceExports("ontophylo")))

