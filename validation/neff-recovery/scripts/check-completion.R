# Technical completion gate only; correlations never determine acceptance.
arguments <- commandArgs(trailingOnly = TRUE)
trial_ids <- as.integer(arguments)
if (!length(trial_ids)) trial_ids <- 1:100
manifest <- read.csv("design/task-manifest.csv")
for (trial in trial_ids) {
  out_dir <- sprintf("results/trial-%03d", trial)
  settings <- readRDS(file.path(out_dir, "inference-settings.rds"))
  status <- read.csv(file.path(out_dir, "inference-status.csv"))
  retained <- read.csv(file.path(out_dir, "character-retention.csv"))
  fits <- read.csv(file.path(out_dir, "fit-status.csv"))
  maps <- read.csv(file.path(out_dir, "map-status.csv"))
  seed_row <- manifest[manifest$trial == trial, ]
  stopifnot(nrow(retained) == 200L, nrow(fits) == 200L, nrow(maps) == 200L,
    status$status %in% c("ok", "ok_constants_omitted"),
    settings$n_maps == 100L, settings$nstarts == 3L,
    settings$fit_seed == seed_row$fit_seed, settings$mapping_seed == seed_row$mapping_seed,
    all(maps$n_maps[maps$status == "ok"] == 100L),
    !any(maps$status == "failed"), !any(fits$status == "failed"),
    setequal(retained$character[retained$variable], fits$character[fits$status == "ok"]),
    setequal(retained$character[retained$variable], maps$character[maps$status == "ok"]),
    file.exists(file.path(out_dir, "summary-complete.txt")))
  cat("Trial", trial, ":", sum(retained$variable), "variable characters,",
      sum(maps$n_maps), "maps; technical checks passed.\n")
}
