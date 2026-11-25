# Empirical Neff+ branch-label randomization.
# In RStudio, open this script and set the working directory to its folder.
# From the repository root:
# Rscript validation/neff-recovery/empirical-branch-shuffle/run-shuffles.R
# Requires ape and Rcpp, plus a working C++ compiler for the small Neff helper.
# No tree, character fitting, stochastic mapping or PyRate fitting is rerun.

library(ape)
library(Rcpp)

# Settings ----
n_shuffles <- as.integer(Sys.getenv("NEFF_SHUFFLES", "9999"))
seed <- 20261123
stopifnot(n_shuffles >= 1)
set.seed(seed)

analysis_dir <- "."
if (file.exists("validation/neff-recovery/empirical-branch-shuffle/run-shuffles.R")) {
  analysis_dir <- "validation/neff-recovery/empirical-branch-shuffle"
}
analysis_dir <- normalizePath(analysis_dir)
repo_dir <- normalizePath(file.path(analysis_dir, "../../.."))
output_dir <- Sys.getenv("NEFF_OUTPUT_DIR", file.path(analysis_dir, "results"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
sourceCpp(file.path(analysis_dir, "neff-kernel.cpp"))

# Existing empirical inputs ----
tree <- readRDS(file.path(repo_dir, "data/hym_tree.RDS"))
branch_rates <- readRDS(file.path(repo_dir, "data_out/br_rates_all.RDS"))
neff_saved <- readRDS(file.path(repo_dir, "data_out/neff-per-replicate.RDS"))
pyrate <- readRDS(file.path(repo_dir, "data_out/pyrate-dive.RDS"))
units <- names(branch_rates)[names(branch_rates) != "larval"]
n_units <- length(units)
n_maps <- ncol(branch_rates[[1]])
n_edges <- nrow(tree$edge)
stopifnot(length(tree$tip.label) == 77, n_edges == 152,
          n_units == 15, n_maps == 1000)

# Accelerated branches: reproduce the original inclusive 95th percentile.
# Keep every map, including ties and any all-zero rate samples.
accelerated <- vector("list", n_units)
names(accelerated) <- units
counts <- matrix(0L, nrow = n_maps, ncol = n_units,
                 dimnames = list(NULL, units))
all_zero_rate_maps <- integer(n_units)
zero_threshold_maps <- integer(n_units)
observed_patterns <- matrix(0L, nrow = n_edges, ncol = n_maps)
for (unit in 1:n_units) {
  accelerated[[unit]] <- matrix(0L, nrow = n_edges, ncol = n_maps)
  for (map in 1:n_maps) {
    rates <- branch_rates[[units[unit]]][, map]
    stopifnot(all(is.finite(rates)))
    threshold <- quantile(rates, 0.95, names = FALSE)
    zero_threshold_maps[unit] <- zero_threshold_maps[unit] + (threshold == 0)
    accelerated[[unit]][, map] <- as.integer(rates >= threshold)
    counts[map, unit] <- sum(accelerated[[unit]][, map])
    all_zero_rate_maps[unit] <- all_zero_rate_maps[unit] + all(rates == 0)
  }
  # A 15-bit integer encodes each complete anatomical acceleration pattern.
  observed_patterns <- observed_patterns + accelerated[[unit]] * 2^(unit - 1)
}
storage.mode(observed_patterns) <- "integer"
count_summary <- data.frame(unit = units, minimum = NA_integer_,
                            median = NA_real_, maximum = NA_integer_,
                            all_zero_rate_maps = all_zero_rate_maps,
                            zero_threshold_maps = zero_threshold_maps,
                            all_branches_accelerated_maps = NA_integer_)
for (unit in 1:n_units) {
  count_summary$minimum[unit] <- min(counts[, unit])
  count_summary$median[unit] <- median(counts[, unit])
  count_summary$maximum[unit] <- max(counts[, unit])
  count_summary$all_branches_accelerated_maps[unit] <- sum(counts[, unit] == n_edges)
}
write.csv(count_summary, file.path(output_dir, "accelerated-branch-counts.csv"),
          row.names = FALSE)

# Exactly the original 1,000 slices, active-branch rule and age assignment ----
depth <- node.depth.edgelength(tree)
edge_start <- depth[tree$edge[, 1]]
edge_end <- depth[tree$edge[, 2]]
slice_time <- seq(0, max(depth), length.out = 1000)
slice_time[1000] <- slice_time[1000] - 1e-4
slice_age <- rev(slice_time)
slice_edges <- vector("list", 1000)
slice_key <- character(1000)
for (slice in 1:1000) {
  slice_edges[[slice]] <- which(slice_time[slice] >= edge_start &
                                slice_time[slice] <= edge_end)
  slice_key[slice] <- paste(slice_edges[[slice]], collapse = ",")
}
# Repeated slices with the same branches have identical Neff, so compute once.
group <- match(slice_key, unique(slice_key))
first_slice <- match(unique(slice_key), slice_key)
active_edges <- slice_edges[first_slice]
n_groups <- length(active_edges)

# Thirteen complete 20-Myr bins (260-0 Ma), as in the manuscript correlation.
bins <- data.frame(bin = 1:13, age_old = seq(260, 20, -20),
                   age_young = seq(240, 0, -20))
bins$age <- (bins$age_old + bins$age_young) / 2
bins$pyrate <- pyrate$mean[match(bins$bin, pyrate$bin)]
bin_weights <- matrix(0, nrow = 13, ncol = n_groups)
for (bin in 1:13) {
  keep <- slice_age >= bins$age_young[bin] & slice_age < bins$age_old[bin]
  bin_weights[bin, ] <- tabulate(group[keep], nbins = n_groups) / sum(keep)
}

# Verify the C++ helper against the readable definition on actual samples ----
reference <- numeric(n_groups)
for (map in 1:10) {
  for (g in 1:n_groups) {
    patterns <- observed_patterns[active_edges[[g]], map]
    positive <- patterns[patterns != 0L]
    if (length(positive) > 0) {
      probabilities <- as.numeric(table(positive)) / length(positive)
      H_plus <- -sum(probabilities * log(probabilities))
      f_plus <- length(positive) / length(patterns)
      reference[g] <- reference[g] + f_plus * exp(H_plus) / 10
    }
  }
}
fast_check <- mean_neff_by_group(observed_patterns[, 1:10], active_edges)
kernel_error <- max(abs(reference - fast_check))
stopifnot(kernel_error < 1e-12)

# Reproduce the saved mean curve before doing any randomization.
observed_group <- mean_neff_by_group(observed_patterns, active_edges)
saved_order <- match(slice_age, neff_saved$time)
stopifnot(!anyNA(saved_order))
curve_error <- max(abs(observed_group[group] - neff_saved$mean_entropy[saved_order]))
stopifnot(curve_error < 1e-10)
bins$neff <- as.numeric(bin_weights %*% observed_group)
bins$lineages <- NA_real_
for (bin in 1:13) {
  keep <- slice_age >= bins$age_young[bin] & slice_age < bins$age_old[bin]
  n_lineages <- integer(sum(keep))
  selected_slices <- which(keep)
  for (i in seq_along(selected_slices)) {
    n_lineages[i] <- length(slice_edges[[selected_slices[i]]])
  }
  bins$lineages[bin] <- mean(n_lineages)
}
write.csv(bins, file.path(output_dir, "observed-bins.csv"), row.names = FALSE)
cat("Reproduced saved Neff: maximum difference", curve_error, "\n")

# Primary: unrestricted shuffling, as requested.
# Sensitivity: shuffle only among branches in the same length quartile.
# In each unit, move its complete 1,000-map branch profile together. Thus each
# map retains its accelerated-branch count and the unit retains its uncertainty.
# Draw a different permutation for each unit; keep it fixed across all slices.
length_breaks <- quantile(tree$edge.length, seq(0, 1, 0.25), names = FALSE)
length_class <- cut(tree$edge.length, length_breaks, include.lowest = TRUE,
                    labels = FALSE)
length_audit <- data.frame(quartile = 1:4, lower_Myr = head(length_breaks, -1),
                           upper_Myr = tail(length_breaks, -1),
                           n_branches = NA_integer_, accelerated_fraction = NA_real_)
for (quartile in 1:4) {
  edges <- which(length_class == quartile)
  total <- 0
  for (unit in 1:n_units) total <- total + sum(accelerated[[unit]][edges, ])
  length_audit$n_branches[quartile] <- length(edges)
  length_audit$accelerated_fraction[quartile] <- total / (length(edges) * n_units * n_maps)
}
write.csv(length_audit, file.path(output_dir, "branch-length-audit.csv"), row.names = FALSE)
null_types <- c("unrestricted", "length_quartile")
null_groups <- vector("list", length(null_types))
names(null_groups) <- null_types
null_bins <- vector("list", length(null_types))
names(null_bins) <- null_types

for (null_type in null_types) {
  group_values <- matrix(NA_real_, nrow = n_groups, ncol = n_shuffles)
  bin_values <- matrix(NA_real_, nrow = 13, ncol = n_shuffles)
  started <- proc.time()[3]
  for (shuffle in 1:n_shuffles) {
    shuffled_patterns <- matrix(0L, nrow = n_edges, ncol = n_maps)
    for (unit in 1:n_units) {
      permutation <- seq_len(n_edges)
      if (null_type == "unrestricted") {
        permutation <- sample.int(n_edges)
      } else {
        for (quartile in 1:4) {
          edges <- which(length_class == quartile)
          permutation[edges] <- sample(edges)
        }
      }
      shuffled_unit <- accelerated[[unit]][permutation, , drop = FALSE]
      if (shuffle == 1) {
        stopifnot(identical(colSums(shuffled_unit), colSums(accelerated[[unit]])))
        if (null_type == "length_quartile") {
          for (quartile in 1:4) {
            edges <- which(length_class == quartile)
            stopifnot(identical(colSums(shuffled_unit[edges, , drop = FALSE]),
                               colSums(accelerated[[unit]][edges, , drop = FALSE])))
          }
        }
      }
      shuffled_patterns <- shuffled_patterns + shuffled_unit * 2^(unit - 1)
    }
    storage.mode(shuffled_patterns) <- "integer"
    group_values[, shuffle] <- mean_neff_by_group(shuffled_patterns, active_edges)
    bin_values[, shuffle] <- bin_weights %*% group_values[, shuffle]
    if (shuffle == 1 || shuffle %% 250 == 0 || shuffle == n_shuffles) {
      cat(null_type, shuffle, "/", n_shuffles, ":",
          round(proc.time()[3] - started, 1), "seconds\n")
      flush.console()
    }
  }
  null_groups[[null_type]] <- group_values
  null_bins[[null_type]] <- bin_values
  saveRDS(list(null_type = null_type, group_neff = group_values,
               binned_neff = bin_values, n_shuffles = n_shuffles, seed = seed),
          file.path(output_dir, paste0(null_type, "-shuffles.rds")))
}

# Test whole trajectories. No permutation of time bins or PyRate rates.
# Raw Pearson r is primary; detrending and adjacent changes are secondary.
statistic_names <- c("raw", "linear_detrended", "first_difference")
observed_tests <- data.frame(statistic = statistic_names, n = c(13, 13, 12),
                             r = NA_real_, nominal_p = NA_real_)
null_statistics <- data.frame(null_type = rep(null_types, each = n_shuffles * 3),
                               shuffle = rep(rep(1:n_shuffles, each = 3), 2),
                               statistic = rep(statistic_names, n_shuffles * 2),
                               r = NA_real_, nominal_p = NA_real_)
for (statistic in statistic_names) {
  x <- bins$neff
  y <- bins$pyrate
  if (statistic == "linear_detrended") {
    x <- resid(lm(x ~ bins$age))
    y <- resid(lm(y ~ bins$age))
  }
  if (statistic == "first_difference") {
    x <- diff(x)
    y <- diff(y)
  }
  test <- cor.test(x, y, method = "pearson")
  row <- match(statistic, statistic_names)
  observed_tests$r[row] <- unname(test$estimate)
  observed_tests$nominal_p[row] <- test$p.value
}
row <- 0
for (null_type in null_types) {
  for (shuffle in 1:n_shuffles) {
    for (statistic in statistic_names) {
      x <- null_bins[[null_type]][, shuffle]
      y <- bins$pyrate
      if (statistic == "linear_detrended") {
        x <- resid(lm(x ~ bins$age))
        y <- resid(lm(y ~ bins$age))
      }
      if (statistic == "first_difference") {
        x <- diff(x)
        y <- diff(y)
      }
      row <- row + 1
      r <- cor(x, y)
      df <- length(x) - 2
      t_value <- abs(r) * sqrt(df / (1 - r^2))
      null_statistics$r[row] <- r
      null_statistics$nominal_p[row] <- 2 * pt(t_value, df = df, lower.tail = FALSE)
    }
  }
}
stopifnot(all(is.finite(null_statistics$r)))
summary <- data.frame()
for (null_type in null_types) {
  for (statistic in statistic_names) {
    observed <- observed_tests[observed_tests$statistic == statistic, ]
    selected <- null_statistics$null_type == null_type &
                null_statistics$statistic == statistic
    null <- null_statistics[selected, ]
    n_at_least_observed <- sum(null$r >= observed$r)
    n_absolute_at_least_observed <- sum(abs(null$r) >= abs(observed$r))
    result <- data.frame(null_type = null_type, statistic = statistic,
                          n_shuffles = n_shuffles, observed_r = observed$r,
                          observed_nominal_p = observed$nominal_p,
                          null_r_median = median(null$r),
                          null_r_q025 = quantile(null$r, 0.025, names = FALSE),
                          null_r_q975 = quantile(null$r, 0.975, names = FALSE),
                          n_at_least_observed = n_at_least_observed,
                          permutation_p_positive = (n_at_least_observed + 1) / (n_shuffles + 1),
                          permutation_p_absolute = (n_absolute_at_least_observed + 1) / (n_shuffles + 1),
                          nominal_p_below_005_percent = 100 * mean(null$nominal_p < 0.05),
                          positive_nominal_p_below_005_percent =
                            100 * mean(null$nominal_p < 0.05 & null$r > 0))
    summary <- rbind(summary, result)
  }
}
write.csv(observed_tests, file.path(output_dir, "observed-tests.csv"), row.names = FALSE)
write.csv(null_statistics, file.path(output_dir, "null-correlations.csv"), row.names = FALSE)
write.csv(summary, file.path(output_dir, "permutation-summary.csv"), row.names = FALSE)
saveRDS(list(tree = tree, units = units, counts = counts,
             slice_age = slice_age, group = group, active_edges = active_edges,
             observed_group = observed_group, bins = bins, bin_weights = bin_weights,
             null_groups = null_groups, null_bins = null_bins,
             observed_tests = observed_tests, summary = summary,
             length_class = length_class, length_breaks = length_breaks,
             kernel_error = kernel_error, curve_error = curve_error,
             seed = seed, n_shuffles = n_shuffles),
        file.path(output_dir, "analysis.rds"))
inputs <- c("data/hym_tree.RDS", "data_out/br_rates_all.RDS",
            "data_out/neff-per-replicate.RDS", "data_out/pyrate-dive.RDS")
checksums <- tools::md5sum(file.path(repo_dir, inputs))
write.csv(data.frame(path = inputs, md5 = unname(checksums)),
          file.path(output_dir, "input-checksums.csv"), row.names = FALSE)
write.csv(data.frame(check = c("kernel_vs_plain_R", "saved_mean_curve"),
                      maximum_absolute_error = c(kernel_error, curve_error)),
          file.path(output_dir, "validation.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output_dir, "session-info.txt"))
print(summary, digits = 5, row.names = FALSE)
