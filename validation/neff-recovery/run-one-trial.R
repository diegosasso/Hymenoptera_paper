# Open neff-recovery.Rproj in RStudio, then Source this file.
# The delivered results already exist. They will be reused, not redrawn.
# For a fresh experiment, copy scripts/, design/ and provenance/ to a new folder.

# 1. Choose one of the prespecified trials.
trial <- 1L
stopifnot(trial >= 1L, trial <= 100L)

# 2. Install the frozen package once if needed (dependencies must be available).
# source("scripts/install-package.R")

# 3. Run the three readable analysis stages with this trial ID.
rscript <- file.path(R.home("bin"), "Rscript")
status <- system2(rscript, c("--vanilla", "scripts/01-simulate.R", as.character(trial)))
stopifnot(status == 0L)
status <- system2(rscript, c("--vanilla", "scripts/02-reconstruct.R", as.character(trial)))
stopifnot(status == 0L)
status <- system2(rscript, c("--vanilla", "scripts/03-summarize.R", as.character(trial)))
stopifnot(status == 0L)

# 4. Inspect this trial's correlations.
correlations <- read.csv(sprintf("results/trial-%03d/correlations.csv", trial))
print(correlations)
