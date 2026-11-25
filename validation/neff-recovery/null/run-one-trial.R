# Open null-recovery.Rproj, then Source this file in RStudio.
# The corresponding associated tree must exist in ../simulated/trial-001/.
# Install the parent's frozen OntoPhylo first, or set NEFF_RLIB to that library.
trial <- 1L
rscript <- file.path(R.home("bin"), "Rscript")
scripts <- c("scripts/01-simulate.R", "scripts/02-reconstruct.R",
  "scripts/03-summarize.R", "scripts/check-completion.R", "scripts/04-audit.R")
for (script in scripts) {
  status <- system2(rscript, c("--vanilla", script, trial))
  if (status != 0L) stop("Execution failed: ", script)
}
