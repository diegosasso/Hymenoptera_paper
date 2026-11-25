# Collect every prespecified trial. No filtering or replacement by significance.
manifest <- read.csv("design/task-manifest.csv")
stopifnot(nrow(manifest) == 100L, identical(manifest$trial, 1:100))
trials <- correlations <- recovery <- curves <- bins <- marginal <- inference <- data.frame()
missing <- integer()
for (trial in manifest$trial) {
  out_dir <- sprintf("results/trial-%03d", trial)
  if (!file.exists(file.path(out_dir, "summary-complete.txt"))) {
    missing <- c(missing, trial)
    next
  }
  trials <- rbind(trials, read.csv(file.path(out_dir, "trial-summary.csv")))
  correlations <- rbind(correlations, read.csv(file.path(out_dir, "correlations.csv")))
  recovered <- read.csv(file.path(out_dir, "recovery.csv"))
  if (nrow(recovered)) recovery <- rbind(recovery, recovered)
  curves <- rbind(curves, read.csv(file.path(out_dir, "neff-curves.csv")))
  bins <- rbind(bins, read.csv(file.path(out_dir, "neff-20Ma-bins.csv")))
  marginal <- rbind(marginal, read.csv(file.path(out_dir, "marginal-rates.csv")))
  inference <- rbind(inference, read.csv(file.path(out_dir, "inference-status.csv")))
}
write.csv(data.frame(trial = missing), "summary/missing-trials.csv", row.names = FALSE)
if (length(missing)) stop("Incomplete batch: see summary/missing-trials.csv")
stopifnot(nrow(trials) == 100L, !anyDuplicated(trials$trial))
write.csv(trials, "summary/trials.csv", row.names = FALSE)
write.csv(correlations, "summary/correlations.csv", row.names = FALSE)
write.csv(recovery, "summary/recovery.csv", row.names = FALSE)
write.csv(curves, "summary/neff-curves.csv", row.names = FALSE)
write.csv(bins, "summary/neff-20Ma-bins.csv", row.names = FALSE)
write.csv(inference, "summary/inference-status.csv", row.names = FALSE)
saveRDS(marginal, "summary/marginal-rates.rds")
# Two-sided thresholds and positive-direction detections are both reported.
percentages <- data.frame()
for (source in c("truth", "truth_retained", "reconstructed")) {
  for (cutoff in c("strict", "inclusive")) {
    d <- correlations[correlations$source == source & correlations$cutoff == cutoff, ]
    finite <- is.finite(d$nominal_p) & is.finite(d$Pearson_r)
    n_estimable <- sum(finite)
    for (alpha in c(.05, .01, .001)) {
      for (direction in c("two_sided_any_sign", "positive_r")) {
        detected <- finite & d$nominal_p < alpha
        if (direction == "positive_r") detected <- detected & d$Pearson_r > 0
        n_detected <- sum(detected)
        interval <- binom.test(n_detected, 100L)$conf.int * 100
        percent_estimable <- lower_estimable <- upper_estimable <- NA_real_
        if (n_estimable > 0L) {
          percent_estimable <- 100 * n_detected / n_estimable
          conditional_interval <- binom.test(n_detected, n_estimable)$conf.int * 100
          lower_estimable <- conditional_interval[1]
          upper_estimable <- conditional_interval[2]
        }
        percentages <- rbind(percentages, data.frame(source, cutoff, alpha, direction,
          n_planned = 100L, n_estimable, n_unavailable = 100L - n_estimable,
          n_detected, percent_all_trials = n_detected,
          lower_all_trials = interval[1], upper_all_trials = interval[2],
          percent_estimable, lower_estimable, upper_estimable))
      }
    }
  }
}
write.csv(percentages, "summary/p-value-percentages.csv", row.names = FALSE)
# Across-trial bands summarize variation in trial means, not map uncertainty or
# confidence intervals for the grand mean. Their source is stated in captions.
trajectory_summary <- data.frame()
for (source in unique(curves$source)) {
  for (cutoff in unique(curves$cutoff)) {
    for (age in sort(unique(curves$age), decreasing = TRUE)) {
      x <- curves$neff[curves$source == source & curves$cutoff == cutoff & curves$age == age]
      bounds <- quantile(x, c(.025, .5, .975), names = FALSE)
      trajectory_summary <- rbind(trajectory_summary, data.frame(source, cutoff, age,
        n_trials = length(x), mean = mean(x), median = bounds[2], lower = bounds[1], upper = bounds[3]))
    }
  }
}
write.csv(trajectory_summary, "summary/neff-trajectory-summary.csv", row.names = FALSE)
# Save a machine-readable compact overview for manuscripts and plot checks.
overview <- data.frame(metric = c("n_trials", "n_reconstructed", "n_constant_characters",
  "trials_with_constants", "median_recovery_r", "median_RMSE", "median_truth_r_div",
  "median_reconstructed_r_div", "n_reconstructed_positive_p05"), value = NA_real_)
primary <- correlations[correlations$source == "reconstructed" & correlations$cutoff == "strict", ]
true <- correlations[correlations$source == "truth" & correlations$cutoff == "strict", ]
r <- recovery[recovery$reference == "truth", ]
overview$value <- c(nrow(trials), nrow(primary), sum(trials$n_constant), sum(trials$n_constant > 0),
  median(r$trajectory_r), median(r$RMSE), median(true$Pearson_r), median(primary$Pearson_r),
  sum(primary$Pearson_r > 0 & primary$nominal_p < .05, na.rm = TRUE))
write.csv(overview, "summary/overview.csv", row.names = FALSE)
print(percentages[percentages$cutoff == "strict" & percentages$direction == "positive_r", ])
print(overview)
