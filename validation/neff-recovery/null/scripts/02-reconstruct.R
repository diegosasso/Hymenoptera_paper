# Run from this project: Rscript --vanilla scripts/02-reconstruct.R TRIAL
# Only tree and terminal observations are supplied to inference.
source("scripts/load-packages.R")
arguments <- commandArgs(trailingOnly = TRUE)
trial <- as.integer(arguments[1])
stopifnot(trial >= 1L, trial <= 100L)
sim_dir <- sprintf("simulated/trial-%03d", trial)
out_dir <- sprintf("results/trial-%03d", trial)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
simulation <- readRDS(file.path(sim_dir, "simulation.rds"))
simulation_settings <- readRDS(file.path(sim_dir, "settings.rds"))
observations <- simulation$scenarios$root_er$observations
workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
constant_policy <- Sys.getenv("NEFF_CONSTANT_POLICY", "omit")
stopifnot(constant_policy %in% c("omit", "unavailable"))
settings <- list(trial = trial, models = "ER", root_p = "maddfitz", nstarts = 3L,
  n_maps = 100L, fit_seed = simulation_settings$character_seed + 1L,
  mapping_seed = simulation_settings$character_seed + 2L,
  constant_policy = constant_policy,
  simulation_md5 = unname(tools::md5sum(file.path(sim_dir, "simulation.rds"))),
  package_md5 = unname(tools::md5sum("provenance/ontophylo-source-sha256.json")))
settings_file <- file.path(out_dir, "inference-settings.rds")
if (file.exists(settings_file)) stopifnot(identical(settings, readRDS(settings_file)))
if (!file.exists(settings_file)) saveRDS(settings, settings_file)
if (file.exists(file.path(out_dir, "inference-status.csv"))) {
  cat("Keeping completed inference attempt for trial", trial, "\n")
  quit(save = "no", status = 0)
}
started <- Sys.time()
capture.output(sessionInfo(), file = file.path(out_dir, "session-info.txt"))
checks <- simulation$scenarios$root_er$checks
variable <- checks$character[checks$tips_0 > 0 & checks$tips_1 > 0]
retention <- data.frame(character = checks$character, unit = checks$unit,
  variable = checks$character %in% variable, tips_0 = checks$tips_0, tips_1 = checks$tips_1)
write.csv(retention, file.path(out_dir, "character-retention.csv"), row.names = FALSE)
# Constant observations do not imply zero historical changes. Their simulated
# histories are retained, but no fitted Q or zero-change map is invented.
if (constant_policy == "unavailable" && length(variable) < ncol(observations)) {
  status <- data.frame(trial, status = "constant_characters", n_variable = length(variable),
    n_fit_failed = 0L, n_map_failed = 0L, n_maps = 0L, elapsed_seconds = 0)
  write.csv(status, file.path(out_dir, "inference-status.csv"), row.names = FALSE)
  quit(save = "no", status = 0)
}
# Retry only a worker-port collision, with unchanged scientific seeds.
for (attempt in 1:3) {
  if (workers > 1L) {
    port <- as.integer(system2("python3", "scripts/select-worker-port.py", stdout = TRUE))
    parallel:::setDefaultClusterOptions(port = port)
  }
  fits <- tryCatch(ontophylo::fit_models(simulation$tree, observations,
    models = "ER", root_p = "maddfitz", nstarts = settings$nstarts,
    workers = workers, seed = settings$fit_seed, outdir = file.path(out_dir, "fit-checkpoints"),
    resume = TRUE, on_error = "continue"), error = identity)
  if (!inherits(fits, "error")) break
  collision <- grepl("creation of server socket failed: port", conditionMessage(fits), fixed = TRUE)
  if (!collision || attempt == 3L) stop(fits)
}
stopifnot(setequal(fits$status$character[fits$status$status == "skipped"],
                  setdiff(colnames(observations), variable)))
saveRDS(fits, file.path(out_dir, "fits.rds.tmp"))
stopifnot(file.rename(file.path(out_dir, "fits.rds.tmp"), file.path(out_dir, "fits.rds")))
write.csv(fits$status, file.path(out_dir, "fit-status.csv"), row.names = FALSE)
fit_diagnostics <- data.frame()
for (id in variable) {
  current <- fits$characters[[id]]
  if (current$status != "ok") next
  model <- current$model_table
  fit_diagnostics <- rbind(fit_diagnostics, data.frame(character = id,
    fitted_q = current$Q_matrix[1, 2], model))
}
write.csv(fit_diagnostics, file.path(out_dir, "fit-diagnostics.csv"), row.names = FALSE)
for (attempt in 1:3) {
  if (workers > 1L) {
    port <- as.integer(system2("python3", "scripts/select-worker-port.py", stdout = TRUE))
    parallel:::setDefaultClusterOptions(port = port)
  }
  mapped <- tryCatch(ontophylo::paramo_stm(fits, n_stm = settings$n_maps,
    batch_size = 25L, store = "branch_counts", workers = workers,
    seed = settings$mapping_seed, outdir = file.path(out_dir, "map-checkpoints"),
    resume = TRUE, on_error = "continue"), error = identity)
  if (!inherits(mapped, "error")) break
  collision <- grepl("creation of server socket failed: port", conditionMessage(mapped), fixed = TRUE)
  if (!collision || attempt == 3L) stop(mapped)
}
write.csv(mapped$status, file.path(out_dir, "map-status.csv"), row.names = FALSE)
# Materialize counts so the final object is portable; checkpoints remain on disk.
successful <- mapped$status$character[mapped$status$status == "ok"]
counts <- ontophylo::read_stm(mapped, successful, type = "branch_counts")
for (id in successful) {
  mapped$characters[[id]]$branch_counts <- counts[[id]]
  mapped$characters[[id]]$files <- NULL
}
saveRDS(mapped, file.path(out_dir, "maps.rds.tmp"))
stopifnot(file.rename(file.path(out_dir, "maps.rds.tmp"), file.path(out_dir, "maps.rds")))
status_label <- "ok"
if (length(variable) < ncol(observations)) status_label <- "ok_constants_omitted"
if (!setequal(successful, variable)) status_label <- "inference_failed"
status <- data.frame(trial, status = status_label, n_variable = length(variable),
  n_fit_failed = sum(fits$status$status == "failed"),
  n_map_failed = sum(mapped$status$status == "failed"),
  n_maps = sum(mapped$status$n_maps),
  elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")))
write.csv(status, file.path(out_dir, "inference-status.csv"), row.names = FALSE)
print(status)
