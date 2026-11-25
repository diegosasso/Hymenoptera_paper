# Numerical inserts come directly from the complete, independently checked tables.
p <- read.csv("summary/comparison-percentages.csv")
t <- read.csv("summary/trials.csv")
i <- read.csv("summary/inference-status.csv")
r <- read.csv("summary/recovery.csv")
c <- read.csv("summary/correlations.csv")
null <- p[p$scenario == "Root ER control" & p$cutoff == "strict" & p$alpha == .05, ]
positive <- null[null$source == "reconstructed" & null$direction == "positive_r", ]
any_sign <- null[null$source == "reconstructed" & null$direction == "two_sided_any_sign", ]
truth_positive <- null[null$source == "truth" & null$direction == "positive_r", ]
truth_any <- null[null$source == "truth" & null$direction == "two_sided_any_sign", ]
recovery <- r[r$reference == "truth", ]
inclusive <- c[c$cutoff == "inclusive", ]
stopifnot(nrow(t) == 100L, nrow(positive) == 1L, nrow(any_sign) == 1L)

values <- c(NullPositive = positive$n_detected, NullAny = any_sign$n_detected,
  NullNegative = any_sign$n_detected - positive$n_detected,
  NullTruthPositive = truth_positive$n_detected, NullTruthAny = truth_any$n_detected,
  NullPositiveLower = sprintf("%.1f", positive$lower),
  NullPositiveUpper = sprintf("%.1f", positive$upper),
  NullAnyLower = sprintf("%.1f", any_sign$lower), NullAnyUpper = sprintf("%.1f", any_sign$upper),
  NullRecoveryR = sprintf("%.3f", median(recovery$trajectory_r)),
  NullRecoveryRMSE = sprintf("%.3f", median(recovery$RMSE)),
  NullConstants = sum(t$n_constant), NullConstantTrials = sum(t$n_constant > 0),
  NullReconstructed = sum(i$status %in% c("ok", "ok_constants_omitted")),
  NullInclusiveTruth = sum(is.finite(inclusive$nominal_p[inclusive$source == "truth"])),
  NullInclusiveReconstructed = sum(is.finite(inclusive$nominal_p[inclusive$source == "reconstructed"])))
macros <- character()
for (name in names(values)) {
  macros <- c(macros, paste0("\\newcommand{\\", name, "}{", values[name], "}"))
}
writeLines(macros, "manuscript/results-macros.tex")
for (name in c("main-text-methods", "main-text-results", "supplementary-methods")) {
  lines <- readLines(paste0("manuscript/", name, ".tex"))
  for (key in names(values)) lines <- gsub(paste0("\\", key, "{}"), values[key], lines, fixed = TRUE)
  writeLines(lines, paste0("manuscript/", name, "-standalone.tex"))
}

row_end <- intToUtf8(c(92L, 92L))
table <- c("\\begin{table}[ht]", "\\centering",
  "\\caption{Trial counts meeting nominal two-sided thresholds. Each entry is also a percentage because each scenario has 100 planned trials. T: exact simulated histories; R: ER reconstruction. Strict acceleration cutoff.}",
  "\\label{tab:paired-neff}", "\\begin{tabular}{llrrrr}", "\\toprule",
  paste0(" & & \\multicolumn{2}{c}{Associated} & \\multicolumn{2}{c}{Root ER control} ", row_end),
  paste0("Direction & Threshold & T & R & T & R ", row_end), "\\midrule")
for (direction in c("two_sided_any_sign", "positive_r")) {
  for (alpha in c(.05, .01, .001)) {
    counts <- integer()
    for (scenario in c("Associated", "Root ER control")) {
      for (source in c("truth", "reconstructed")) {
        row <- p[p$scenario == scenario & p$source == source & p$cutoff == "strict" &
                   p$direction == direction & p$alpha == alpha, ]
        counts <- c(counts, row$n_detected)
      }
    }
    label <- if (direction == "positive_r") "Positive $r$" else "Either sign"
    table <- c(table, paste0(sprintf("%s & $p<%s$ & %d & %d & %d & %d ",
      label, alpha, counts[1], counts[2], counts[3], counts[4]), row_end))
  }
}
table <- c(table, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(table, "manuscript/results-table.tex")

summary <- c("HOMOGENEOUS ROOT ER CONTROL: 100 PAIRED TRIALS", "",
  "Q fixed before this batch: q01=q10=0.003 per Myr; root state 0; no redraws.",
  sprintf("Reconstructed trials: %s / 100", values["NullReconstructed"]),
  sprintf("Characters: 20000 simulated; %d constant observations omitted; %d variable fitted", sum(t$n_constant), sum(t$n_variable)),
  sprintf("Mapping histories completed: %d", sum(i$n_maps)),
  sprintf("Positive reconstructed correlations with nominal p<.05: %d/100; exact 95%% CI %.1f-%.1f%%", positive$n_detected, positive$lower, positive$upper),
  sprintf("Reconstructed correlations of either sign with p<.05: %d/100; exact 95%% CI %.1f-%.1f%%", any_sign$n_detected, any_sign$lower, any_sign$upper),
  sprintf("Significant negative reconstructed correlations: %d/100", any_sign$n_detected - positive$n_detected),
  sprintf("True-history controls at p<.05: %d positive; %d either sign", truth_positive$n_detected, truth_any$n_detected),
  "Paired associated benchmark: 78/100 reconstructed, 71/100 true-history positive detections.",
  sprintf("Median true-reconstructed trajectory r: %.6f; RMSE: %.6f", median(recovery$trajectory_r), median(recovery$RMSE)), "",
  "All thresholds, directions and intervals: summary/comparison-percentages.csv",
  "Every correlation: summary/comparison-correlations.csv",
  "Strict acceleration is primary. Inclusive-cutoff results are retained separately.",
  "P-values are nominal two-sided Pearson tests of 15 matched 20-Myr bins.",
  "Null frequencies are specific to homogeneous root ER on this conditioned tree ensemble.",
  "Different character generators/rates mean the scenario contrast does not isolate timing alone.",
  "No result-dependent character redraws, tree replacement or significance-based gate was used.")
writeLines(summary, "RESULTS.txt")
cat(paste(summary, collapse = "\n"), "\n")
