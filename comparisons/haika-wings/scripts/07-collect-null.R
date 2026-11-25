# Compare observed associations with the conditional BM/proxy null.
# Run after 05-compare-rates.R; --final requires all 199 replicates.
arguments <- commandArgs(trailingOnly = TRUE)
if (length(setdiff(arguments, "--final")) > 0L) stop("Only --final is supported.")
final <- "--final" %in% arguments
objects <- readRDS("results/comparison-objects.rds")
datasets <- readRDS("data/prepared-datasets.rds")
crosswalk <- read.csv("data/haika-family-proxy-crosswalk.csv", stringsAsFactors = FALSE)
partitions <- c("fore_wing", "fore_wing_sclerites", "all_wings", "hind_wing")
dataset_ids <- c("seed3407_full77", "seed3407_family48_first")
input_md5 <- tools::md5sum(c("data/prepared-datasets.rds",
                             "data/haika-family-proxy-crosswalk.csv"))
spearman_statistic <- function(first, second) {
  if (length(first) != length(second) || length(first) < 3L ||
      any(!is.finite(first)) || any(!is.finite(second)) ||
      length(unique(first)) < 2L || length(unique(second)) < 2L) return(NA_real_)
  cor(first, second, method = "spearman")
}
reference <- list()
observed_rows <- list()
row_index <- 0L
for (dataset_id in dataset_ids) {
  projection_id <- objects$dataset_projection_ids[[dataset_id]]
  projected <- objects$ontophylo_projections[[projection_id]]
  cache <- objects$ontophylo_summaries[[projection_id]]
  tree <- projected$tree
  branch_families <- integer(nrow(tree$edge))
  for (edge_index in seq_len(nrow(tree$edge))) {
    descendants <- strsplit(projected$branches$descendant_tips[edge_index], "|", fixed = TRUE)[[1]]
    families <- crosswalk$selected_family[match(descendants, crosswalk$analysis_label)]
    stopifnot(!anyNA(families))
    branch_families[edge_index] <- length(unique(families))
  }
  scopes <- list(all_edges = rep(TRUE, nrow(tree$edge)),
                 exclude_root_edges = tree$edge[, 1] != length(tree$tip.label) + 1L,
                 interfamily_edges = branch_families >= 2L)
  for (scope in names(scopes)) {
    stopifnot(identical(which(scopes[[scope]]), cache$edge_subsets[[scope]]))
  }
  exposure <- projected$marginal$contracted_branch_exposure
  total_exposure <- rowSums(exposure)
  time_keep <- total_exposure > 0 & diff(projected$marginal$breaks_Ma) >= 5
  stopifnot(identical(which(time_keep), cache$temporal_subsets$full_width_bins))
  onto_reference <- list()
  for (partition in partitions) {
    summary <- cache$partitions[[partition]]
    # Reuse the exact means from 05, including its handling of invalid maps.
    # Reordering means and matrix products can perturb ties in Spearman ranks.
    branch_mean <- summary$relative_summary$mean
    temporal_mean <- summary$contracted_relative_summary$mean
    stopifnot(length(branch_mean) == nrow(tree$edge),
              length(temporal_mean) == nrow(exposure))
    onto_reference[[partition]] <- list(branch = branch_mean, temporal = temporal_mean)
  }
  reference[[dataset_id]] <- list(tree = tree, scopes = scopes, exposure = exposure,
                                  total_exposure = total_exposure,
                                  time_keep = time_keep, onto = onto_reference)
  axes <- c(colnames(datasets[[dataset_id]]$traits), "TOTAL_VARIANCE")
  for (axis in axes) {
    estimate_key <- paste(dataset_id, axis, 1, sep = "::")
    estimate <- objects$haika_estimates[[estimate_key]]
    if (is.null(estimate)) stop("Missing observed lambda=1 fit: ", estimate_key)
    for (partition in partitions) {
      for (scope in names(scopes)) {
        keep <- scopes[[scope]]
        row_index <- row_index + 1L
        observed_rows[[row_index]] <- data.frame(dataset = dataset_id, axis = axis,
          partition = partition, scope = scope, n = sum(keep),
          statistic = spearman_statistic(estimate$relative[keep],
                                         onto_reference[[partition]]$branch[keep]))
      }
      row_index <- row_index + 1L
      observed_rows[[row_index]] <- data.frame(dataset = dataset_id, axis = axis,
        partition = partition, scope = "time_full_width_bins", n = sum(time_keep),
        statistic = spearman_statistic(estimate$marginal_relative[time_keep],
                                        onto_reference[[partition]]$temporal[time_keep]))
    }
  }
}
observed <- do.call(rbind, observed_rows)
observed$key <- paste(observed$dataset, observed$axis, observed$partition, observed$scope, sep = "::")
stopifnot(!anyDuplicated(observed$key))
observed$statistic_defined <- is.finite(observed$statistic)
null_statistics <- matrix(NA_real_, nrow(observed), 199L,
                          dimnames = list(observed$key, paste0("rep", seq_len(199L))))
status <- data.frame(replicate = seq_len(199L), exists = FALSE, success = FALSE,
                      message = "Result file is absent.", stringsAsFactors = FALSE)
for (replicate_id in seq_len(199L)) {
  path <- sprintf("results/null/rep-%04d.rds", replicate_id)
  if (!file.exists(path)) next
  status$exists[replicate_id] <- TRUE
  replicate_error <- ""
  replicate_statistics <- setNames(rep(NA_real_, nrow(observed)), observed$key)
  tryCatch({
    simulated <- readRDS(path)
    if (!isTRUE(simulated$success)) stop("One or more null BM fits failed.")
    stopifnot(identical(as.integer(simulated$replicate), replicate_id),
              identical(simulated$source, "seed3407"), simulated$lambda == 1)
    for (input_file in names(input_md5)) {
      if (is.null(simulated$input_md5) || !input_file %in% names(simulated$input_md5) ||
          !identical(unname(simulated$input_md5[input_file]), unname(input_md5[input_file]))) {
        stop("Null input checksum mismatch: ", input_file)
      }
    }
    for (dataset_id in dataset_ids) {
      ref <- reference[[dataset_id]]
      current <- simulated$results[[dataset_id]]
      stopifnot(identical(current$tree$edge, ref$tree$edge),
                identical(current$tree$tip.label, ref$tree$tip.label),
                length(current$tree$edge.length) == length(ref$tree$edge.length),
                max(abs(current$tree$edge.length - ref$tree$edge.length)) <= 1e-10)
      rates <- current$axis_rates
      expected_axes <- colnames(datasets[[dataset_id]]$traits)
      stopifnot(is.matrix(rates), nrow(rates) == nrow(ref$tree$edge),
                !anyDuplicated(colnames(rates)), setequal(colnames(rates), expected_axes),
                all(is.finite(rates)), all(rates > 0))
      rates <- rates[, expected_axes, drop = FALSE]
      total_variance <- rep(0, nrow(rates))
      for (axis in expected_axes) total_variance <- total_variance + rates[, axis]
      rates <- cbind(rates, TOTAL_VARIANCE = total_variance)
      relative <- matrix(NA_real_, nrow(rates), ncol(rates), dimnames = dimnames(rates))
      curves <- matrix(NA_real_, nrow(ref$exposure), ncol(rates),
                        dimnames = list(NULL, colnames(rates)))
      # Preserve the observed pipeline's arithmetic ordering exactly: integrate
      # raw rates, divide by exposure, then divide by the treewide mean.
      for (axis in colnames(rates)) {
        mean_rate <- sum(rates[, axis] * ref$tree$edge.length) / sum(ref$tree$edge.length)
        relative[, axis] <- rates[, axis] / mean_rate
        raw_curve <- as.numeric(ref$exposure %*% rates[, axis]) / ref$total_exposure
        raw_curve[ref$total_exposure == 0] <- NA_real_
        curves[, axis] <- raw_curve / mean_rate
      }
      for (axis in colnames(relative)) {
        for (partition in partitions) {
          for (scope in names(ref$scopes)) {
            keep <- ref$scopes[[scope]]
            key <- paste(dataset_id, axis, partition, scope, sep = "::")
            replicate_statistics[key] <- spearman_statistic(relative[keep, axis],
              ref$onto[[partition]]$branch[keep])
          }
          key <- paste(dataset_id, axis, partition, "time_full_width_bins", sep = "::")
          replicate_statistics[key] <- spearman_statistic(curves[ref$time_keep, axis],
            ref$onto[[partition]]$temporal[ref$time_keep])
        }
      }
    }
    if (!all(is.finite(replicate_statistics[observed$statistic_defined]))) {
      stop("At least one required null association is undefined or nonfinite.")
    }
  }, error = function(error) {
    replicate_error <<- conditionMessage(error)
  })
  if (!identical(replicate_error, "")) {
    status$message[replicate_id] <- replicate_error
    next
  }
  status$success[replicate_id] <- TRUE
  status$message[replicate_id] <- "Validated inputs, trees, axes, rates and statistics."
  null_statistics[, replicate_id] <- replicate_statistics
}
write.csv(status, "results/null-status.csv", row.names = FALSE)
if (final && !all(status$success)) stop("Null collection incomplete: ", sum(status$success), "/199 valid replicates")
valid <- status$success
stopifnot(sum(valid) > 0)
summary <- observed
summary$n_null <- rowSums(is.finite(null_statistics[, valid, drop = FALSE]))
summary$null_complete <- all(status$success)
summary$null_median <- NA_real_
summary$null_lower <- NA_real_
summary$null_upper <- NA_real_
summary$p_positive_concordance <- NA_real_
summary$p_two_sided_tail <- NA_real_
for (statistic_index in seq_len(nrow(summary))) {
  if (!summary$statistic_defined[statistic_index]) next
  values <- null_statistics[statistic_index, valid]
  stopifnot(all(is.finite(values)))
  intervals <- quantile(values, c(0.025, 0.5, 0.975))
  summary$null_lower[statistic_index] <- intervals[1]
  summary$null_median[statistic_index] <- intervals[2]
  summary$null_upper[statistic_index] <- intervals[3]
  above <- (1 + sum(values >= summary$statistic[statistic_index])) / (length(values) + 1)
  below <- (1 + sum(values <= summary$statistic[statistic_index])) / (length(values) + 1)
  summary$p_positive_concordance[statistic_index] <- above
  summary$p_two_sided_tail[statistic_index] <- min(1, 2 * min(above, below))
}
summary$interpretation <- "Plug-in homogeneous multivariate BM null with fixed PCA/tree and identical family-proxy duplication; unadjusted exploratory probabilities"
write.csv(summary, "results/null-association-summary.csv", row.names = FALSE)
saveRDS(list(summary = summary, null_statistics = null_statistics, status = status),
        "results/null-comparison.rds", compress = "xz")
cat("Collected", sum(valid), "of 199 valid conditional null replicates.\n")
