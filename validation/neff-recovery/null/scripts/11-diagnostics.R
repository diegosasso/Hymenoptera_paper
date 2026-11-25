# Descriptive diagnostics; none is used to accept, reject or rerun a trial.
trials <- read.csv("summary/trials.csv")
correlations <- read.csv("summary/correlations.csv")
recovery <- read.csv("summary/recovery.csv")
fits <- data.frame()
for (trial in 1:100) {
  file <- sprintf("results/trial-%03d/fit-diagnostics.csv", trial)
  d <- read.csv(file)
  d$trial <- trial
  fits <- rbind(fits, d)
}
write.csv(fits, "summary/fit-diagnostics.csv", row.names = FALSE)
distributions <- data.frame()
quantities <- list(fitted_q = fits$fitted_q,
  recovery_r = recovery$trajectory_r[recovery$reference == "truth"],
  recovery_RMSE = recovery$RMSE[recovery$reference == "truth"],
  recovery_bias = recovery$mean_bias[recovery$reference == "truth"],
  true_diversification_r = correlations$Pearson_r[correlations$source == "truth" & correlations$cutoff == "strict"],
  reconstructed_diversification_r = correlations$Pearson_r[correlations$source == "reconstructed" & correlations$cutoff == "strict"])
for (quantity in names(quantities)) {
  values <- quantities[[quantity]]
  values <- values[is.finite(values)]
  bounds <- quantile(values, c(0, .025, .5, .975, 1), names = FALSE)
  distributions <- rbind(distributions, data.frame(quantity, n = length(values),
    minimum = bounds[1], lower = bounds[2], median = bounds[3], upper = bounds[4],
    maximum = bounds[5], mean = mean(values)))
}
write.csv(distributions, "summary/diagnostic-distributions.csv", row.names = FALSE)
simulation_check <- data.frame(n_characters = sum(trials$n_characters),
  realized_transitions = sum(trials$true_changes),
  expected_transitions = sum(trials$expected_changes),
  observed_over_expected = sum(trials$true_changes) / sum(trials$expected_changes))
write.csv(simulation_check, "summary/transition-count-check.csv", row.names = FALSE)
print(distributions)
print(simulation_check)
