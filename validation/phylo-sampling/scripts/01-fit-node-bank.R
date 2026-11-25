# Fit each unique node pattern once on each tree; repeated trials reuse its MLE.
# Slurm array 1-19, 40 patterns per chunk; character seeds do not depend on chunks.
source("scripts/load-packages.R")
arguments <- commandArgs(trailingOnly = TRUE)
chunk <- as.integer(arguments[1])
workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
design <- readRDS("design/design.rds")
characters <- readRDS("design/node-characters.rds")
first <- (chunk - 1L) * 40L + 1L
last <- min(chunk * 40L, nrow(design$nodes))
stopifnot(first <= last, first >= 1L)
ids <- design$nodes$character[first:last]
out_dir <- sprintf("fits/chunk-%02d", chunk)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
for (tree_name in c("full", "study")) {
  tree <- design[[paste0(tree_name, "_tree")]]
  states <- characters[[tree_name]][, ids, drop = FALSE]
  checkpoint_dir <- file.path(out_dir, tree_name)
  fitting_seed <- 320261120L + match(tree_name, c("full", "study"))
  fits <- fit_models(tree, states, models = "ER", root_p = "maddfitz",
    nstarts = 5L, workers = workers, seed = fitting_seed,
    outdir = checkpoint_dir, resume = TRUE, on_error = "continue")
  successful <- fits$status$character[fits$status$status == "ok"]
  variable <- colSums(states == "1") > 0L & colSums(states == "1") < Ntip(tree)
  stopifnot(!any(fits$status$status == "failed"), setequal(successful, ids[variable]))
  saveRDS(fits, file.path(out_dir, paste0(tree_name, ".rds")))
  write.csv(fits$status, file.path(out_dir, paste0(tree_name, "-status.csv")), row.names = FALSE)
}
cat("Completed fit chunk", chunk, "\n")

