# Descriptive concordance and retained-character sensitivity; no new tests.
c <- read.csv("summary/correlations.csv")
r <- read.csv("summary/recovery.csv")
t <- read.csv("summary/trials.csv")
truth <- c[c$source == "truth" & c$cutoff == "strict", ]
reconstructed <- c[c$source == "reconstructed" & c$cutoff == "strict", ]
paired <- merge(truth, reconstructed, by = "trial", suffixes = c("_truth", "_reconstructed"))
concordance <- data.frame()
for (alpha in c(.05, .01, .001)) {
  true_detected <- paired$Pearson_r_truth > 0 & paired$nominal_p_truth < alpha
  reconstructed_detected <- paired$Pearson_r_reconstructed > 0 & paired$nominal_p_reconstructed < alpha
  concordance <- rbind(concordance, data.frame(alpha,
    both_detected = sum(true_detected & reconstructed_detected),
    truth_only = sum(true_detected & !reconstructed_detected),
    reconstructed_only = sum(!true_detected & reconstructed_detected),
    neither_detected = sum(!true_detected & !reconstructed_detected)))
}
write.csv(concordance, "summary/detection-concordance.csv", row.names = FALSE)
recovery_summary <- data.frame()
for (reference in c("truth", "truth_retained")) {
  d <- r[r$reference == reference, ]
  for (metric in c("trajectory_r", "RMSE", "mean_bias")) {
    values <- d[[metric]]
    bounds <- quantile(values, c(.025, .5, .975), names = FALSE)
    recovery_summary <- rbind(recovery_summary, data.frame(reference, metric,
      n = length(values), mean = mean(values), median = bounds[2],
      central_95_lower = bounds[1], central_95_upper = bounds[3],
      minimum = min(values), maximum = max(values)))
  }
}
write.csv(recovery_summary, "summary/recovery-summary.csv", row.names = FALSE)
constant_summary <- data.frame()
for (source in c("truth", "truth_retained", "reconstructed")) {
  d <- c[c$source == source & c$cutoff == "inclusive", ]
  constant_summary <- rbind(constant_summary, data.frame(source, n_tests = nrow(d),
    n_undefined = sum(!is.finite(d$nominal_p))))
}
write.csv(constant_summary, "summary/inclusive-cutoff-diagnostics.csv", row.names = FALSE)
print(concordance)
print(recovery_summary)
print(constant_summary)
