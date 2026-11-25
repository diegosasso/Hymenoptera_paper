# Run one prespecified HAIKA BM fit from the comparison directory.
# Usage: Rscript --vanilla scripts/04-fit-bm.R TASK_ID

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1) {
  stop("Provide exactly one task ID from data/fit-manifest.csv.")
}
task_id <- suppressWarnings(as.integer(arguments[1]))
if (is.na(task_id) || task_id < 1 || as.character(task_id) != arguments[1]) {
  stop("The task ID must be a positive integer.")
}

project_library <- Sys.getenv("HAIKA_RLIB",
                              "/projappl/project_2004996/haika-wings/Rlib")
.libPaths(c(project_library,
            "/projappl/project_2004996/phylo-sampling/Rlib",
            "/projappl/project_2004996/hymenoptera-node-pattern-bias/Rlib",
            .libPaths()))

dir.create("results/fits", recursive = TRUE, showWarnings = FALSE)
result_file <- sprintf("results/fits/task-%04d.rds", task_id)
diagnostic_file <- sprintf("results/fits/task-%04d-diagnostics.csv", task_id)
summary_file <- sprintf("results/fits/task-%04d-summary.csv", task_id)
time_started <- Sys.time()
elapsed_started <- proc.time()[["elapsed"]]
manifest_row <- NULL
dataset_metadata <- NULL
fitted_model <- NULL
error_message <- ""

tryCatch({
  source("scripts/bm-model.R")
  manifest <- read.csv("data/fit-manifest.csv", stringsAsFactors = FALSE,
                       check.names = FALSE)
  required_fields <- c("task_id", "dataset", "axis", "lambda", "seed")
  if (!all(required_fields %in% names(manifest))) {
    stop("The fit manifest does not contain all required columns.")
  }
  selected_rows <- which(manifest$task_id == task_id)
  if (length(selected_rows) != 1) {
    stop("The requested task ID must match exactly one manifest row.")
  }
  manifest_row <- manifest[selected_rows, , drop = FALSE]
  prepared_datasets <- readRDS("data/prepared-datasets.rds")
  dataset_name <- manifest_row$dataset[1]
  axis_name <- manifest_row$axis[1]
  if (!dataset_name %in% names(prepared_datasets)) {
    stop("The manifest dataset is absent from prepared-datasets.rds.")
  }
  dataset <- prepared_datasets[[dataset_name]]
  if (is.null(dataset$tree) || is.null(dataset$traits)) {
    stop("The prepared dataset must contain tree and traits.")
  }
  if (!axis_name %in% colnames(dataset$traits)) {
    stop("The manifest axis is absent from the prepared trait matrix.")
  }
  if (is.null(rownames(dataset$traits))) {
    stop("Prepared trait matrices must have taxon row names.")
  }
  trait_values <- dataset$traits[, axis_name]
  trait_values <- as.numeric(trait_values)
  names(trait_values) <- rownames(dataset$traits)
  dataset_metadata <- dataset$metadata

  cat("Task:", task_id, "Dataset:", dataset_name, "Axis:", axis_name,
      "Lambda:", manifest_row$lambda[1], "Tips:", length(trait_values), "\n")
  cat("Started:", format(time_started, tz = "UTC", usetz = TRUE), "\n")
  flush.console()
  fitted_model <- fit_bm_bounded(dataset$tree, trait_values,
                                 lambda = manifest_row$lambda[1],
                                 n_starts = 3L,
                                 seed = manifest_row$seed[1],
                                 maxit = 2000L)
  if (!isTRUE(fitted_model$success)) {
    if (!is.null(fitted_model$error)) {
      error_message <- fitted_model$error
    } else {
      error_message <- paste("Fit failed convergence or likelihood validation;",
                             "convergence code", fitted_model$convergence)
    }
  }
}, error = function(error) {
  error_message <<- conditionMessage(error)
})

elapsed_seconds <- proc.time()[["elapsed"]] - elapsed_started
success <- !is.null(fitted_model) && isTRUE(fitted_model$success) &&
  identical(error_message, "")
task_summary <- data.frame(task_id = task_id, dataset = NA_character_,
                           axis = NA_character_, lambda = NA_real_, seed = NA_integer_,
                           success = success, elapsed_seconds = elapsed_seconds,
                           started_utc = format(time_started, tz = "UTC", usetz = TRUE),
                           ended_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
                           error = error_message, stringsAsFactors = FALSE)
if (!is.null(manifest_row)) {
  task_summary$dataset <- manifest_row$dataset[1]
  task_summary$axis <- manifest_row$axis[1]
  task_summary$lambda <- manifest_row$lambda[1]
  task_summary$seed <- manifest_row$seed[1]
}

if (!is.null(fitted_model$diagnostics)) {
  diagnostics <- fitted_model$diagnostics
} else {
  diagnostics <- data.frame(start = NA_integer_, objective = NA_real_,
                             convergence = NA_integer_,
                             message = error_message, stringsAsFactors = FALSE)
}
diagnostics$task_id <- task_id
diagnostics$dataset <- task_summary$dataset
diagnostics$axis <- task_summary$axis
diagnostics$lambda <- task_summary$lambda
diagnostics$seed <- task_summary$seed
diagnostics$task_success <- success
diagnostics$total_elapsed_seconds <- elapsed_seconds

if (success) {
  result <- fitted_model
} else {
  # Failed candidates are retained only as diagnostics, never as usable rates.
  result <- list(success = FALSE, error = error_message,
                 convergence = fitted_model$convergence,
                 diagnostics = diagnostics,
                 likelihood_check = fitted_model$likelihood_check,
                 objective_check = fitted_model$objective_check,
                 phytools_version = fitted_model$phytools_version)
}
result$task_id <- task_id
result$manifest <- manifest_row
result$dataset_metadata <- dataset_metadata
result$task_summary <- task_summary
result$session_info <- utils::sessionInfo()
result$input_md5 <- tools::md5sum(c("data/fit-manifest.csv",
                                   "data/prepared-datasets.rds",
                                   "scripts/bm-model.R", "scripts/bm-cached-objective.R",
                                   "scripts/04-fit-bm.R"))
saveRDS(result, result_file)
write.csv(diagnostics, diagnostic_file, row.names = FALSE)
write.csv(task_summary, summary_file, row.names = FALSE)
print(task_summary)

if (!success) {
  quit(save = "no", status = 1L)
}
