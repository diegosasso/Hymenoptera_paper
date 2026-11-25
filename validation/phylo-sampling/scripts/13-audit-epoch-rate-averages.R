# Compare lineage-time averages with time averages of the marginal curves.
# Run from validation/phylo-sampling. Original results and figures are unchanged.
design <- readRDS("design/design.rds")
trials <- read.csv("summary/paired-trials.csv")
conditions <- read.csv("summary/condition-summary.csv")

# Within each interval between node ages, the marginal rate is constant.
# Integrate these intervals exactly, including the boundaries of each epoch.
weights <- list()
lineage_times <- list()
grid_weights <- list()
maximum_tip_roundoff <- 0
for (tree_name in c("full", "study")) {
  tree <- design[[paste0(tree_name, "_tree")]]
  ages <- design[[paste0(tree_name, "_age")]]
  # Treat extant tips as age zero; stored ages differ by at most 1.2e-7 Myr.
  tip_roundoff <- max(abs(ages[seq_along(tree$tip.label)]))
  stopifnot(tip_roundoff < 1e-6)
  maximum_tip_roundoff <- max(maximum_tip_roundoff, tip_roundoff)
  ages[seq_along(tree$tip.label)] <- 0
  parent_age <- ages[tree$edge[, 1]]
  child_age <- ages[tree$edge[, 2]]
  breaks <- sort(unique(c(0, 80, 160, 280, ages[ages > 0 & ages < 280])))
  time_weights <- matrix(0, nrow(tree$edge), 3)
  exposure <- matrix(0, nrow(tree$edge), 3)

  for (j in seq_len(length(breaks) - 1L)) {
    duration <- breaks[j + 1L] - breaks[j]
    midpoint <- (breaks[j + 1L] + breaks[j]) / 2
    epoch <- findInterval(midpoint, design$epochs$lower_ma)
    epoch_duration <- design$epochs$upper_ma[epoch] - design$epochs$lower_ma[epoch]
    active <- which(parent_age > midpoint & child_age < midpoint)
    stopifnot(length(active) > 0L)
    time_weights[active, epoch] <- time_weights[active, epoch] +
      duration / length(active) / epoch_duration
    exposure[active, epoch] <- exposure[active, epoch] + duration
  }

  stopifnot(max(abs(colSums(time_weights) - 1)) < 1e-10)
  original_exposure <- design[[paste0(tree_name, "_exposure")]]
  stopifnot(max(abs(exposure - original_exposure)) <= tip_roundoff + 1e-8)
  # Retain the original exposure convention when reproducing the old statistic.
  lineage_weights <- sweep(original_exposure, 2, colSums(original_exposure), "/")
  weights[[tree_name]] <- list(time = time_weights, lineage = lineage_weights)
  lineage_times[[tree_name]] <- colSums(exposure)

  # Independently reproduce the saved 1-Myr marginal values for audit trials.
  grid <- matrix(0, nrow(tree$edge), length(design$age_grid))
  for (j in seq_along(design$age_grid)) {
    age <- design$age_grid[j]
    active <- which(parent_age >= age & child_age < age)
    if (age == 0) active <- which(tree$edge[, 2] <= length(tree$tip.label))
    stopifnot(length(active) > 0L)
    grid[active, j] <- 1 / length(active)
  }
  grid_weights[[tree_name]] <- grid
}

audit <- trials[, c("task", "condition", "epoch", "proportion", "trial")]
audit$full_lineage_mean <- audit$study_lineage_mean <- NA_real_
audit$full_time_mean <- audit$study_time_mean <- NA_real_
maximum_branch_error <- maximum_grid_error <- 0
for (i in seq_len(nrow(audit))) {
  result <- readRDS(sprintf("results/task-%04d/result.rds", audit$task[i]))
  epoch <- audit$epoch[i]
  for (tree_name in c("full", "study")) {
    group <- result$groups[[paste0(tree_name, "_all")]]
    branch_mean <- rowMeans(group$rate_samples)
    maximum_branch_error <- max(maximum_branch_error,
      abs(branch_mean - group$branch_summary$mean))
    audit[i, paste0(tree_name, "_lineage_mean")] <-
      sum(branch_mean * weights[[tree_name]]$lineage[, epoch])
    audit[i, paste0(tree_name, "_time_mean")] <-
      sum(branch_mean * weights[[tree_name]]$time[, epoch])
    if (result$trial == 1L || result$condition == 5L) {
      grid_mean <- as.vector(crossprod(grid_weights[[tree_name]], branch_mean))
      maximum_grid_error <- max(maximum_grid_error,
        abs(grid_mean - group$marginal_summary$mean))
    }
  }
  if (i %% 250L == 0L) cat("Checked", i, "of", nrow(audit), "paired trials\n")
}
audit$lineage_ratio <- audit$study_lineage_mean / audit$full_lineage_mean
audit$time_ratio <- audit$study_time_mean / audit$full_time_mean
maximum_existing_error <- max(abs(audit$full_lineage_mean - trials$full_epoch_rate),
  abs(audit$study_lineage_mean - trials$study_epoch_rate),
  abs(audit$lineage_ratio - trials$study_full_rate_ratio))
stopifnot(nrow(audit) == 1500L, maximum_existing_error < 1e-9,
          maximum_branch_error < 1e-9, maximum_grid_error < 1e-9)

summary <- conditions[, c("condition", "epoch", "proportion", "n_trials")]
summary$full_lineage_mean <- summary$study_lineage_mean <- NA_real_
summary$lineage_ratio_mean <- summary$full_time_mean <- summary$study_time_mean <- NA_real_
summary$time_ratio_mean <- summary$time_ratio_lower <- summary$time_ratio_upper <- NA_real_
summary$ratio_of_mean_curves <- NA_real_
for (i in seq_len(nrow(summary))) {
  one <- audit[audit$condition == summary$condition[i], ]
  summary$full_lineage_mean[i] <- mean(one$full_lineage_mean)
  summary$study_lineage_mean[i] <- mean(one$study_lineage_mean)
  summary$lineage_ratio_mean[i] <- mean(one$lineage_ratio)
  summary$full_time_mean[i] <- mean(one$full_time_mean)
  summary$study_time_mean[i] <- mean(one$study_time_mean)
  summary$time_ratio_mean[i] <- mean(one$time_ratio)
  summary$time_ratio_lower[i] <- quantile(one$time_ratio, .025, names = FALSE)
  summary$time_ratio_upper[i] <- quantile(one$time_ratio, .975, names = FALSE)
  summary$ratio_of_mean_curves[i] <- mean(one$study_time_mean) / mean(one$full_time_mean)
}
stopifnot(max(abs(summary$lineage_ratio_mean - conditions$rate_ratio_mean)) < 1e-9)

write.csv(audit, "summary/epoch-average-audit-trials.csv", row.names = FALSE)
write.csv(summary, "summary/epoch-average-audit-summary.csv", row.names = FALSE)
report <- c("Audit: different weights explain the epoch-rate discrepancy.",
  "Lineage-time average = integral(N(t) * marginal(t)) / integral(N(t)).",
  "Time average of marginal curve = integral(marginal(t)) / epoch duration.",
  "Integrals are exact over intervals between node ages, not 1-Myr approximations.",
  paste("Paired trials checked:", nrow(audit)),
  paste("Maximum difference from existing epoch rates and ratios:", maximum_existing_error),
  paste("Maximum branch-mean discrepancy:", maximum_branch_error),
  paste("Maximum discrepancy from saved marginal curves:", maximum_grid_error),
  paste("Maximum extant-tip age roundoff set to zero for integration (Myr):", maximum_tip_roundoff),
  "", capture.output(print(summary, row.names = FALSE)),
  "", "Original simulation results, summaries and figures were not overwritten.")
writeLines(report, "provenance/epoch-average-audit.txt")
cat(report, sep = "\n")
