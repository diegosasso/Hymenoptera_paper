# Run on Roihu: Rscript --vanilla scripts/03-run-trial.R TASK_ID
source("scripts/load-packages.R")
source("scripts/rate-summaries.R")
arguments <- commandArgs(trailingOnly = TRUE)
task <- as.integer(arguments[1])
workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
stopifnot(is.finite(task), task >= 1L, task <= 1500L)
design <- readRDS("design/design.rds")
dataset <- readRDS(sprintf("simulated/task-%04d.rds", task))
selected <- dataset$selected_nodes
ids <- selected$character
retained_ids <- selected$character[selected$variable_study]
out_dir <- sprintf("results/task-%04d", task)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
settings <- list(models = "ER", root_p = "maddfitz", nstarts = 5L, n_maps = 100L,
  normalize = "none", mapping_seed = 420261120L + task * 10L,
  simulation = tools::md5sum(sprintf("simulated/task-%04d.rds", task)),
  full_fit_bank = tools::md5sum("design/node-fits-full.rds"),
  study_fit_bank = tools::md5sum("design/node-fits-study.rds"),
  source = tools::md5sum("provenance/ontophylo-source-sha256.json"))
settings_file <- file.path(out_dir, "settings.rds")
if (file.exists(settings_file)) stopifnot(identical(readRDS(settings_file), settings))
if (!file.exists(settings_file)) saveRDS(settings, settings_file)
result_file <- file.path(out_dir, "result.rds")
if (file.exists(result_file)) {
  stopifnot(identical(readRDS(result_file)$settings, settings))
  cat("Keeping completed task", task, "\n")
  quit(save = "no", status = 0)
}
start_time <- Sys.time()
groups <- list()
diagnostics <- list()
capture.output(sessionInfo(), file = file.path(out_dir, "session-info.txt"))
for (tree_name in c("full", "study")) {
  tree <- design[[paste0(tree_name, "_tree")]]
  fits <- readRDS(sprintf("design/node-fits-%s.rds", tree_name))
  fits$characters <- fits$characters[ids]
  fits$status <- fits$status[match(ids, fits$status$character), , drop = FALSE]
  expected <- ids
  if (tree_name == "study") expected <- retained_ids
  stopifnot(fits$settings$root_p == "maddfitz",
            !any(fits$status$status == "failed"),
            setequal(fits$status$character[fits$status$status == "ok"], expected))
  mapping_seed <- settings$mapping_seed + match(tree_name, c("full", "study"))
  map_file <- file.path(out_dir, paste0("maps-", tree_name, ".rds"))
  cat("Task", task, tree_name, "mapping", length(expected), "characters",
      format(Sys.time()), "\n")
  map_start <- Sys.time()
  if (file.exists(map_file)) {
    mapped <- readRDS(map_file)
  } else {
    # Concurrent cluster startup can occupy a port after it has been probed.
    # Retry only that startup error, always with the same explicit mapping seed.
    for (attempt in 1:3) {
      if (workers > 1L && nzchar(Sys.getenv("SLURM_JOB_ID"))) {
        worker_port <- system2("python3", "scripts/select-worker-port.py", stdout = TRUE)
        stopifnot(length(worker_port) == 1L, !is.na(as.integer(worker_port)))
        # ontophylo uses parallel's default port rather than exposing a port argument.
        parallel:::setDefaultClusterOptions(port = as.integer(worker_port))
        cat("Worker port", worker_port, "startup attempt", attempt, "\n")
      }
      mapped <- tryCatch(paramo_stm(fits, n_stm = settings$n_maps, workers = workers,
        seed = mapping_seed, batch_size = 25L, store = "branch_counts",
        on_error = "continue"), error = identity)
      if (!inherits(mapped, "error")) break
      port_error <- grepl("creation of server socket failed: port", conditionMessage(mapped),
                          fixed = TRUE)
      if (!port_error || attempt == 3L) stop(mapped)
      cat("Worker port occupied at startup; retrying with the same mapping seed.\n")
    }
    saveRDS(mapped, paste0(map_file, ".tmp"))
    stopifnot(file.rename(paste0(map_file, ".tmp"), map_file))
  }
  stopifnot(identical(names(mapped$characters), ids),
            mapped$settings$seed == mapping_seed,
            mapped$settings$n_stm == settings$n_maps,
            !any(mapped$status$status == "failed"),
            setequal(mapped$status$character[mapped$status$status == "ok"], expected),
            all(mapped$status$n_maps[mapped$status$status == "ok"] == settings$n_maps))
  write.csv(fits$status, file.path(out_dir, paste0("fits-", tree_name, ".csv")), row.names = FALSE)
  write.csv(mapped$status, file.path(out_dir, paste0("maps-", tree_name, ".csv")), row.names = FALSE)
  diagnostics[[tree_name]] <- list(fit_status = fits$status, mapping_status = mapped$status,
    map_seconds = as.numeric(difftime(Sys.time(), map_start, units = "secs")))

  ages <- design[[paste0(tree_name, "_age")]]
  exposure <- design[[paste0(tree_name, "_exposure")]]
  apo_edges <- design[[paste0(tree_name, "_apo_edges")]]
  origins <- selected$edge
  if (tree_name == "study") origins <- selected$study_edge[selected$variable_study]
  groups[[paste0(tree_name, "_all")]] <- summarize_group(tree, mapped, expected,
    origins, ages, exposure, apo_edges, design$age_grid, dataset$epoch)
  if (tree_name == "full" && length(retained_ids) > 0L) {
    groups$full_matched <- summarize_group(tree, mapped, retained_ids,
      selected$edge[selected$variable_study], ages, exposure, apo_edges,
      design$age_grid, dataset$epoch)
  }
  rm(mapped, fits)
  gc()
}
result <- list(task = task, condition = dataset$condition, epoch = dataset$epoch,
  proportion = dataset$proportion, trial = dataset$trial, settings = settings,
  n_characters = length(ids), n_retained = length(retained_ids),
  selected_nodes = selected, groups = groups, diagnostics = diagnostics,
  elapsed_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs")))
saveRDS(result, paste0(result_file, ".tmp"))
stopifnot(file.rename(paste0(result_file, ".tmp"), result_file))
cat("Completed task", task, "in", result$elapsed_seconds, "seconds\n")

