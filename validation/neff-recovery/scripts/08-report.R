# Write numerical manuscript inserts directly from the complete result tables.
trials <- read.csv("summary/trials.csv")
p <- read.csv("summary/p-value-percentages.csv")
r <- read.csv("summary/recovery.csv")
correlations <- read.csv("summary/correlations.csv")
inference <- read.csv("summary/inference-status.csv")
primary <- p[p$source == "reconstructed" & p$cutoff == "strict" &
               p$direction == "positive_r" & p$alpha == .05, ]
truth <- p[p$source == "truth" & p$cutoff == "strict" &
             p$direction == "positive_r" & p$alpha == .05, ]
recovery <- r[r$reference == "truth", ]
true_cor <- correlations[correlations$source == "truth" & correlations$cutoff == "strict", ]
estimated_cor <- correlations[correlations$source == "reconstructed" & correlations$cutoff == "strict", ]
paired <- merge(true_cor, estimated_cor, by = "trial", suffixes = c("_truth", "_reconstructed"))
true_detected <- paired$Pearson_r_truth > 0 & paired$nominal_p_truth < .05
estimated_detected <- paired$Pearson_r_reconstructed > 0 & paired$nominal_p_reconstructed < .05
stopifnot(nrow(trials) == 100, nrow(primary) == 1, nrow(truth) == 1)
macros <- c(
  sprintf("\\newcommand{\\NeffBothDetected}{%d}", sum(true_detected & estimated_detected)),
  sprintf("\\newcommand{\\NeffAdditionalDetected}{%d}", sum(!true_detected & estimated_detected)),
  sprintf("\\newcommand{\\NeffRecoveryR}{%.3f}", median(recovery$trajectory_r)),
  sprintf("\\newcommand{\\NeffRecoveryRMSE}{%.3f}", median(recovery$RMSE)),
  sprintf("\\newcommand{\\NeffDetected}{%d}", primary$n_detected),
  sprintf("\\newcommand{\\NeffDetectedPercent}{%.0f}", primary$percent_all_trials),
  sprintf("\\newcommand{\\NeffDetectedLower}{%.1f}", primary$lower_all_trials),
  sprintf("\\newcommand{\\NeffDetectedUpper}{%.1f}", primary$upper_all_trials),
  sprintf("\\newcommand{\\NeffTruthDetected}{%d}", truth$n_detected),
  sprintf("\\newcommand{\\NeffConstantCharacters}{%d}", sum(trials$n_constant)),
  sprintf("\\newcommand{\\NeffConstantTrials}{%d}", sum(trials$n_constant > 0)),
  sprintf("\\newcommand{\\NeffReconstructedTrials}{%d}", primary$n_estimable))
writeLines(macros, "manuscript/results-macros.tex")
# Standalone snippets contain literal numbers for straightforward Overleaf use.
for (name in c("main-text-methods", "supplementary-methods")) {
  text <- readLines(paste0("manuscript/", name, ".tex"))
  for (macro in macros) {
    key <- sub(".*\\\\newcommand\\{(\\\\[^}]+)\\}.*", "\\1", macro)
    value <- sub(".*\\}\\{([^}]+)\\}$", "\\1", macro)
    text <- gsub(paste0(key, "{}"), value, text, fixed = TRUE)
  }
  writeLines(text, paste0("manuscript/", name, "-standalone.tex"))
}
table <- c("\\begin{table}[ht]", "\\centering",
  "\\caption{Positive correlations meeting nominal two-sided thresholds. All denominators are 100 trials.}",
  "\\label{tab:neff-detection}", "\\begin{tabular}{llrr}", "\\toprule",
  "Source & Threshold & Trials (\\%) & Exact 95\\% interval (\\%) \\\\", "\\midrule")
for (source in c("truth", "reconstructed")) {
  for (alpha in c(.05, .01, .001)) {
    row <- p[p$source == source & p$cutoff == "strict" & p$direction == "positive_r" & p$alpha == alpha, ]
    label <- if (source == "truth") "Simulation truth" else "ER reconstruction"
    table <- c(table, sprintf("%s & $p<%s$ & %d (%.0f) & %.1f--%.1f \\\\",
      label, alpha, row$n_detected, row$percent_all_trials, row$lower_all_trials, row$upper_all_trials))
  }
}
table <- c(table, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(table, "manuscript/results-table.tex")
summary <- c("NEFF RECOVERY: 100 PRESPECIFIED INDEPENDENT TRIALS", "",
  sprintf("Reconstructed trials: %d / 100", primary$n_estimable),
  sprintf("Characters: 20000 simulated; %d constant observations omitted; %d variable fitted", sum(trials$n_constant), sum(trials$n_variable)),
  sprintf("Mapping histories completed: %d", sum(inference$n_maps)),
  sprintf("Median Neff trajectory recovery r: %.6f", median(recovery$trajectory_r)),
  sprintf("Median Neff trajectory RMSE: %.6f", median(recovery$RMSE)),
  sprintf("Positive reconstructed Neff/diversification correlations with nominal p<.05: %d/100 (%.1f%%; exact 95%% CI %.1f-%.1f%%)",
    primary$n_detected, primary$percent_all_trials, primary$lower_all_trials, primary$upper_all_trials),
  sprintf("Corresponding true-history detections: %d/100", truth$n_detected), "",
  "All p-value thresholds, directions and denominators: summary/p-value-percentages.csv",
  "Nominal Pearson tests use 15 matched 20-Myr bins and do not correct temporal autocorrelation.",
  "This imposed-association experiment does not estimate null false-positive rates or causality.", "",
  "Truth uses all 200 characters/trial. Matched truth retains only characters used for ASR.",
  "No character or completed tree was redrawn because of its p-value or reconstruction outcome.",
  "The inclusive cutoff is a secondary sensitivity analysis; all 100 true and reconstructed trajectories were constant and their correlations undefined.")
writeLines(summary, "RESULTS.txt")
cat(paste(summary, collapse = "\n"), "\n")
