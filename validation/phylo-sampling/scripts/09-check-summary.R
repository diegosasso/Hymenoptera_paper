# Independently verify the numerical data used in the final publication figures.
trials <- read.csv("summary/paired-trials.csv")
metrics <- read.csv("summary/metric-trials.csv")
curves <- readRDS("summary/marginal-trials.rds")
curve_summary <- read.csv("summary/marginal-summary.csv")
diagnostics <- read.csv("summary/figure-diagnostics.csv")
conditions <- readRDS("design/design.rds")$conditions
stopifnot(nrow(trials) == 1500L, length(unique(trials$task)) == 1500L,
          all(table(trials$condition) == 100L), nrow(curve_summary) == 16860L,
          nrow(diagnostics) == 60L)

# Split once, then check every age, tree, comparison and condition independently.
curve_key <- paste(curves$condition, curves$group, curves$time, sep = ":")
curve_rows <- split(seq_len(nrow(curves)), curve_key)
maximum_curve_error <- 0
for (i in seq_len(nrow(curve_summary))) {
  expected <- curve_summary[i, ]
  key <- paste(expected$condition, expected$group, expected$age, sep = ":")
  one <- curves[curve_rows[[key]], ]
  stopifnot(nrow(one) == 100L, length(unique(one$trial)) == 100L,
            expected$n_trials == 100L, all(is.finite(one$mean)))
  actual <- c(mean(one$mean), quantile(one$mean, c(.025, .975)),
              mean(one$origin_reference))
  recorded <- c(expected$mean, expected$lower, expected$upper,
                expected$origin_reference)
  maximum_curve_error <- max(maximum_curve_error, abs(actual - recorded))
}

# The matched study curves must be identical to the all-character study curves.
all_study <- curves[curves$group == "study_all", ]
matched_study <- curves[curves$group == "study_matched", ]
all_study <- all_study[order(all_study$task, all_study$time), ]
matched_study <- matched_study[order(matched_study$task, matched_study$time), ]
stopifnot(identical(all_study$task, matched_study$task),
          identical(all_study$time, matched_study$time),
          identical(all_study$mean, matched_study$mean),
          identical(all_study$lower, matched_study$lower),
          identical(all_study$upper, matched_study$upper))

# Reconstruct every recovery-panel point and interval from the trial tables.
maximum_diagnostic_error <- 0
for (i in seq_len(nrow(diagnostics))) {
  expected <- diagnostics[i, ]
  one <- trials[trials$condition == expected$condition, ]
  one_metric <- metrics[metrics$condition == expected$condition &
                          metrics$group == "study_all", ]
  panel <- substr(expected$panel, 1, 1)
  if (panel == "a") values <- 100 * one$n_retained / one$n_characters
  if (panel == "b") values <- one$study_epoch_rate / one$full_epoch_rate
  if (panel == "c") values <- 100 * one_metric$mean[
    one_metric$metric == "target_epoch_fraction"]
  if (panel == "d") values <- 100 * one_metric$mean[
    one_metric$metric == "correct_origin_fraction"]
  stopifnot(length(values) == 100L, all(is.finite(values)),
            expected$n_trials == 100L,
            expected$percent == 100 * conditions$proportion[expected$condition])
  actual <- c(mean(values), quantile(values, c(.025, .975)))
  recorded <- c(expected$mean, expected$lower, expected$upper)
  maximum_diagnostic_error <- max(maximum_diagnostic_error, abs(actual - recorded))
}
stopifnot(maximum_curve_error < 1e-10, maximum_diagnostic_error < 1e-10)
lines <- c("PASS: all final marginal and recovery figure summaries verified.",
  "16,860 marginal points and 60 recovery points; 100 trial estimates per point.",
  paste("Maximum marginal-summary disagreement:", maximum_curve_error),
  paste("Maximum recovery-summary disagreement:", maximum_diagnostic_error),
  "Matched and all-character study curves agree exactly.",
  "Ribbons and bars are empirical central 95% ranges across trial means.")
writeLines(lines, "provenance/figure-data-verification.txt")
cat(lines, sep = "\n")
