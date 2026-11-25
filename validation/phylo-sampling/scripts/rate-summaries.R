# Helpers for summaries; fitting and simulation remain in the numbered scripts.
sample_summary <- function(samples) {
  answer <- data.frame(mean = rowMeans(samples), median = NA_real_, sd = NA_real_,
                       lower = NA_real_, upper = NA_real_)
  for (i in seq_len(nrow(samples))) {
    answer$median[i] <- median(samples[i, ])
    answer$sd[i] <- sd(samples[i, ])
    answer$lower[i] <- quantile(samples[i, ], 0.025, names = FALSE)
    answer$upper[i] <- quantile(samples[i, ], 0.975, names = FALSE)
  }
  answer
}

summarize_group <- function(tree, mapped, ids, origin_edges, ages, exposure,
                            apo_edges, age_grid, epoch) {
  n_maps <- mapped$settings$n_stm
  n_edges <- nrow(tree$edge)
  counts <- matrix(0, n_edges, n_maps)
  correct_origin <- numeric(n_maps)
  origin_counts <- tabulate(origin_edges, nbins = n_edges)
  status <- "ok"
  if (length(ids) > 0L) {
    pipeline <- branch_rate_pipeline(mapped, groups = list(combined = ids),
      normalize = "none", time_points = age_grid, time_scale = "age")
    result <- pipeline$groups$combined
    stopifnot(result$normalization == "none", result$divisor == 1,
              identical(result$tree$edge, tree$edge))
    for (i in seq_along(ids)) {
      one <- read_stm(mapped, characters = ids[i], type = "branch_counts")[[1]]
      stopifnot(all(one >= 0), all(one == floor(one)),
                nrow(one) == n_edges, ncol(one) == n_maps)
      counts <- counts + one
      correct_origin <- correct_origin + one[origin_edges[i], ]
    }
    stopifnot(max(abs(result$rate_samples * tree$edge.length - counts)) < 1e-7)
  } else {
    # An empty sum is zero recoverable signal, NOT a fitted zero-rate history.
    status <- "no_variable_characters"
    rates <- matrix(0, n_edges, n_maps)
    zero_curve <- branch_rate_thru_time(tree, rep(0, n_edges),
      time_points = age_grid, time_scale = "age")
    marginal_samples <- matrix(0, length(age_grid), n_maps)
    branch <- cbind(data.frame(edge = seq_len(n_edges), parent = tree$edge[, 1],
      child = tree$edge[, 2], branch_length = tree$edge.length), sample_summary(rates))
    marginal <- cbind(zero_curve[, c("time", "n_branches")], sample_summary(marginal_samples))
    result <- list(tree = tree, branch_summary = branch, rate_samples = rates,
      marginal_summary = marginal, marginal_samples = marginal_samples,
      normalization = "none", divisor = 1)
  }
  reference_rate <- origin_counts / tree$edge.length
  reference_curve <- branch_rate_thru_time(tree, reference_rate,
    time_points = age_grid, time_scale = "age")
  result$branch_summary$origin_count <- origin_counts
  result$branch_summary$origin_reference <- reference_rate
  result$branch_summary$parent_age <- ages[tree$edge[, 1]]
  result$branch_summary$child_age <- ages[tree$edge[, 2]]
  result$marginal_summary$origin_reference <- reference_curve$rate

  changes <- colSums(counts)
  epoch_changes <- crossprod(exposure, result$rate_samples)
  epoch_rate <- epoch_changes / colSums(exposure)
  midpoint <- (ages[tree$edge[, 1]] + ages[tree$edge[, 2]]) / 2
  safe_total <- changes
  safe_total[safe_total == 0] <- NA_real_
  metrics <- data.frame(map = seq_len(n_maps), changes = changes,
    target_epoch_rate = epoch_rate[epoch, ],
    target_epoch_fraction = epoch_changes[epoch, ] / safe_total,
    mean_change_age = as.vector(crossprod(midpoint, counts)) / safe_total,
    correct_origin_fraction = correct_origin / safe_total,
    outside_apocrita_fraction = colSums(counts[!apo_edges, , drop = FALSE]) / safe_total)
  n_origins <- sum(origin_counts)
  metrics$branch_count_error <- NA_real_
  if (n_origins > 0) {
    metrics$branch_count_error <- colSums(abs(counts - origin_counts)) / n_origins
  }
  reference_epoch_changes <- as.vector(crossprod(exposure, reference_rate))
  reference <- data.frame(changes = n_origins,
    target_epoch_rate = reference_epoch_changes[epoch] / sum(exposure[, epoch]),
    target_epoch_fraction = NA_real_, mean_change_age = NA_real_)
  if (n_origins > 0) {
    reference$target_epoch_fraction <- reference_epoch_changes[epoch] / n_origins
    reference$mean_change_age <- sum(midpoint * origin_counts) / n_origins
  }
  result$tree <- NULL
  result$status <- status
  result$n_characters <- length(ids)
  result$character_ids <- ids
  result$metrics <- metrics
  result$reference <- reference
  result$epoch_rate_samples <- epoch_rate
  result
}

