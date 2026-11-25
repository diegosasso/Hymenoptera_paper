# Descriptive recovery of a recent increase in the saved simulations.
# This is a post hoc summary, not a calibrated rate-shift significance test.
# Use the existing lineage-time-weighted epoch rates, comparing adjacent
# equal-duration epochs: 0-80 Ma versus 80-160 Ma.
manifest <- read.csv("design/task-manifest.csv")
recent <- manifest[manifest$epoch == 1L, ]
rows <- list()
for (i in seq_len(nrow(recent))) {
  task <- recent$task[i]
  result <- readRDS(sprintf("results/task-%04d/result.rds", task))
  for (tree in c("full", "study")) {
    group <- result$groups[[paste0(tree, "_all")]]
    rates <- group$epoch_rate_samples
    recent_mean <- mean(rates[1, ])
    previous_mean <- mean(rates[2, ])
    difference <- rates[1, ] - rates[2, ]
    limits <- quantile(difference, c(.025, .975), names = FALSE)
    rows[[length(rows) + 1L]] <- data.frame(task = task,
      condition = result$condition, proportion = result$proportion,
      tree = tree, n_characters = result$n_characters,
      n_retained = result$n_retained,
      recent_rate = recent_mean, previous_rate = previous_mean,
      rate_ratio = recent_mean / previous_mean,
      difference_mean = mean(difference), difference_lower = limits[1],
      difference_upper = limits[2],
      recent_higher = mean(difference) > 1e-12,
      mapping_interval_positive = limits[1] > 1e-12)
  }
}
trials <- do.call(rbind, rows)
write.csv(trials, "summary/recent-signal-trials.csv", row.names = FALSE)
rows <- list()
for (condition in 1:5) {
  for (tree in c("full", "study")) {
    one <- trials[trials$condition == condition & trials$tree == tree, ]
    stopifnot(nrow(one) == 100L)
    rows[[length(rows) + 1L]] <- data.frame(condition = condition,
      proportion = one$proportion[1], tree = tree,
      n_characters = one$n_characters[1], n_trials = nrow(one),
      retained_mean = mean(one$n_retained), retained_min = min(one$n_retained),
      retained_max = max(one$n_retained),
      recent_higher = sum(one$recent_higher),
      mapping_interval_positive = sum(one$mapping_interval_positive),
      ratio_median = median(one$rate_ratio),
      ratio_lower = quantile(one$rate_ratio, .025, names = FALSE),
      ratio_upper = quantile(one$rate_ratio, .975, names = FALSE))
  }
}
summary <- do.call(rbind, rows)
write.csv(summary, "summary/recent-signal-summary.csv", row.names = FALSE)
print(summary, row.names = FALSE)
cat("\nIntervals for each trial condition on fitted Q and the tree.\n",
    "No constant-rate null simulations or shift-model test were performed.\n",
    "These counts describe recovery under the imposed clade-origin design.\n", sep = "")
