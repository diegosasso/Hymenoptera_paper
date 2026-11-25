# Branch-rate and Neff calculations from saved truth and reconstructed counts.
source("scripts/load-packages.R")
arguments <- commandArgs(trailingOnly = TRUE)
trial <- as.integer(arguments[1])
sim_dir <- sprintf("simulated/trial-%03d", trial)
out_dir <- sprintf("results/trial-%03d", trial)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
simulation <- readRDS(file.path(sim_dir, "simulation.rds"))
settings <- readRDS(file.path(sim_dir, "settings.rds"))
ages <- simulation$diversification$age
checks <- simulation$scenarios$root_er$checks
variable <- checks$character[checks$tips_0 > 0 & checks$tips_1 > 0]
groups <- simulation$groups
retained_groups <- groups
for (unit in names(groups)) retained_groups[[unit]] <- intersect(groups[[unit]], variable)
full_units_retained <- all(lengths(retained_groups) > 0L)
truth <- ontophylo::branch_rate_pipeline(simulation$scenarios$root_er$true_maps,
  groups = groups, normalize = "none", time_points = ages, time_scale = "age")
rate_results <- list(truth = truth)
if (full_units_retained) {
  rate_results$truth_retained <- ontophylo::branch_rate_pipeline(
    simulation$scenarios$root_er$true_maps, groups = retained_groups,
    normalize = "none", time_points = ages, time_scale = "age")
}
status_file <- file.path(out_dir, "inference-status.csv")
status <- "inference_not_run"
if (file.exists(status_file)) {
  status_table <- read.csv(status_file)
  status <- status_table$status
  if (status %in% c("ok", "ok_constants_omitted") && full_units_retained) {
    mapped <- readRDS(file.path(out_dir, "maps.rds"))
    stopifnot(all(mapped$status$n_maps[mapped$status$status == "ok"] == 100L))
    rate_results$reconstructed <- ontophylo::branch_rate_pipeline(mapped,
      groups = retained_groups, normalize = "none", time_points = ages, time_scale = "age")
  }
}
metrics <- curves <- bins <- correlations <- acceleration <- marginal <- data.frame()
strict_results <- list()
for (source_name in names(rate_results)) {
  rate_result <- rate_results[[source_name]]
  marginal <- rbind(marginal, data.frame(trial, source = source_name,
                                        rate_result$marginal_summary))
  for (cutoff in c("strict", "inclusive")) {
    ties <- if (cutoff == "strict") "exclude" else "include"
    constant <- if (cutoff == "strict") "unclassified" else "include"
    metric <- ontophylo::phenomic_metrics(rate_result, lower_q = .05, upper_q = .95,
      ties = ties, constant = constant, time_points = ages, time_scale = "age")
    if (cutoff == "strict") {
      strict_results[[source_name]] <- metric
      acceleration <- rbind(acceleration, data.frame(trial, source = source_name,
        metric$regimes$branch_summary))
    }
    d <- metric$summary[metric$summary$metric == "Neff_plus", ]
    curve <- data.frame(trial, source = source_name, cutoff, age = ages,
      time_bin = floor(ages / 20), neff = d$mean, lower = d$lower, upper = d$upper,
      diversification = simulation$diversification$rate)
    curves <- rbind(curves, curve)
    binned <- aggregate(cbind(neff, diversification) ~ time_bin, data = curve, FUN = mean)
    binned <- binned[order(binned$time_bin, decreasing = TRUE), ]
    binned$age <- binned$time_bin * 20 + 10
    bins <- rbind(bins, data.frame(trial, source = source_name, cutoff, binned))
    r <- p <- NA_real_
    reason <- "constant_neff"
    if (sd(binned$neff) > 1e-12 && sd(binned$diversification) > 1e-12) {
      test <- cor.test(binned$neff, binned$diversification, method = "pearson")
      r <- unname(test$estimate)
      p <- test$p.value
      reason <- "estimable"
    }
    correlations <- rbind(correlations, data.frame(trial, source = source_name,
      cutoff, n_bins = nrow(binned), Pearson_r = r, nominal_p = p, reason))
  }
}
recovery <- data.frame()
if ("reconstructed" %in% names(strict_results)) {
  estimated <- rowMeans(strict_results$reconstructed$samples$Neff_plus)
  for (reference in c("truth", "truth_retained")) {
    true <- rowMeans(strict_results[[reference]]$samples$Neff_plus)
    recovery <- rbind(recovery, data.frame(trial, reference,
      trajectory_r = cor(true, estimated), RMSE = sqrt(mean((estimated - true)^2)),
      mean_bias = mean(estimated - true), true_peak_age = ages[which.max(true)],
      reconstructed_peak_age = ages[which.max(estimated)]))
  }
}
trial_summary <- data.frame(trial, status, n_tips = ape::Ntip(simulation$tree),
  tree_proposals = nrow(read.csv(file.path(sim_dir, "tree-attempts.csv"))),
  homogeneous_q = settings$homogeneous_q, n_characters = nrow(checks),
  n_variable = length(variable), n_constant = nrow(checks) - length(variable),
  true_changes = sum(checks$total_changes),
  expected_changes = 200 * settings$homogeneous_q * sum(simulation$tree$edge.length))
write.csv(trial_summary, file.path(out_dir, "trial-summary.csv"), row.names = FALSE)
write.csv(curves, file.path(out_dir, "neff-curves.csv"), row.names = FALSE)
write.csv(bins, file.path(out_dir, "neff-20Ma-bins.csv"), row.names = FALSE)
write.csv(correlations, file.path(out_dir, "correlations.csv"), row.names = FALSE)
write.csv(recovery, file.path(out_dir, "recovery.csv"), row.names = FALSE)
write.csv(acceleration, file.path(out_dir, "branch-acceleration.csv"), row.names = FALSE)
write.csv(marginal, file.path(out_dir, "marginal-rates.csv"), row.names = FALSE)
saveRDS(list(trial = trial, rate_results = rate_results, strict = strict_results),
        file.path(out_dir, "metrics.rds.tmp"))
stopifnot(file.rename(file.path(out_dir, "metrics.rds.tmp"), file.path(out_dir, "metrics.rds")))
writeLines("complete", file.path(out_dir, "summary-complete.txt"))
print(correlations)
