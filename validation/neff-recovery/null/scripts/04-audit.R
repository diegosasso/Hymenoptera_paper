# Independent checks: do not call OntoPhylo metric functions here.
library(ape)
arguments <- commandArgs(trailingOnly = TRUE)
trials <- as.integer(arguments)
if (!length(trials)) trials <- 1:100
checks <- data.frame()
for (trial in trials) {
  sim_dir <- sprintf("simulated/trial-%03d", trial)
  out_dir <- sprintf("results/trial-%03d", trial)
  s <- readRDS(file.path(sim_dir, "simulation.rds"))
  settings <- readRDS(file.path(sim_dir, "settings.rds"))
  current <- readRDS(file.path(out_dir, "metrics.rds"))
  tree <- s$tree
  counts <- s$scenarios$root_er$true_counts
  chars <- s$scenarios$root_er$checks
  events <- s$scenarios$root_er$events
  n_edges <- nrow(tree$edge)
  associated <- readRDS(sprintf("../simulated/trial-%03d/simulation.rds", trial))
  stopifnot(identical(tree, associated$tree), identical(s$groups, associated$groups),
    identical(s$diversification, associated$diversification),
    Ntip(tree) >= 120L, Ntip(tree) <= 250L, is.binary(tree), is.ultrametric(tree),
    min(tree$edge.length) > 1e-5, ncol(counts) == 200L,
    all(events$type == "ER_transition"), sum(counts) == nrow(events),
    all(chars$proposals == 1L), identical(settings$lambda, c(.010, .028, .010)),
    settings$homogeneous_q == .003, settings$root_state == 0L,
    settings$n_maps == 100L, settings$mu == 0,
    length(s$focal_edges$root_er) == 0L,
    s$scenarios$root_er$Q[1, 2] == .003, s$scenarios$root_er$Q[2, 1] == .003)
  # Check full histories against event counts, inherited node states and tips.
  for (id in colnames(counts)) {
    history <- s$scenarios$root_er$true_maps[[id]][[1]]
    states <- rep("0", Ntip(tree) + tree$Nnode)
    for (edge in seq_len(n_edges)) {
      dwell <- history$maps[[edge]]
      stopifnot(abs(sum(dwell) - tree$edge.length[edge]) < 1e-8,
        names(dwell)[1] == states[tree$edge[edge, 1]],
        length(dwell) - 1L == counts[edge, id])
      states[tree$edge[edge, 2]] <- tail(names(dwell), 1)
    }
    stopifnot(identical(unname(s$scenarios$root_er$observations[, id]), states[seq_len(Ntip(tree))]))
  }
  for (id in colnames(counts)) {
    observed <- tabulate(events$edge[events$character == id], nbins = n_edges)
    stopifnot(identical(as.integer(observed), as.integer(counts[, id])))
  }
  mapped <- NULL
  if ("reconstructed" %in% names(current$strict)) mapped <- readRDS(file.path(out_dir, "maps.rds"))
  maximum_error <- 0
  evaluated <- 0L
  for (source_name in names(current$strict)) {
    rates <- current$rate_results[[source_name]]
    metric <- current$strict[[source_name]]
    nmaps <- ncol(rates$groups[[1]]$rate_samples)
    chosen_maps <- unique(c(1L, ceiling(nmaps / 2), nmaps))
    for (m in chosen_maps) {
      mask <- matrix(0L, n_edges, length(rates$groups))
      for (u in seq_along(rates$groups)) {
        rate <- rates$groups[[u]]$rate_samples[, m]
        threshold <- quantile(rate, .95, names = FALSE, type = 7)
        mask[, u] <- as.integer(rate > threshold)
        if (min(rate) == max(rate)) mask[, u] <- 0L
        if (source_name == "reconstructed") {
          ids <- rates$memberships[[u]]
          independent_counts <- numeric(n_edges)
          for (id in ids) {
            independent_counts <- independent_counts + mapped$characters[[id]]$branch_counts[, m]
          }
          stopifnot(max(abs(independent_counts / tree$edge.length - rate)) < 1e-10)
        }
        if (source_name %in% c("truth", "truth_retained")) {
          ids <- rates$memberships[[u]]
          independent_rate <- rowSums(counts[, ids, drop = FALSE]) / tree$edge.length
          stopifnot(max(abs(independent_rate - rate)) < 1e-10)
        }
      }
      for (i in seq_along(s$diversification$age)) {
        age <- s$diversification$age[i]
        active <- which(s$branch_design$age_old >= age & s$branch_design$age_young < age)
        patterns <- mask[active, , drop = FALSE]
        nonzero <- rowSums(patterns) > 0
        expected <- 0
        if (any(nonzero)) {
          codes <- character(sum(nonzero))
          selected <- patterns[nonzero, , drop = FALSE]
          for (j in seq_len(nrow(selected))) codes[j] <- paste(selected[j, ], collapse = "")
          frequency <- as.numeric(table(codes)) / length(codes)
          expected <- mean(nonzero) * exp(-sum(frequency * log(frequency)))
        }
        difference <- abs(expected - metric$samples$Neff_plus[i, m])
        maximum_error <- max(maximum_error, difference)
        evaluated <- evaluated + 1L
      }
    }
  }
  stopifnot(maximum_error < 1e-10)
  corr <- read.csv(file.path(out_dir, "correlations.csv"))
  bins <- read.csv(file.path(out_dir, "neff-20Ma-bins.csv"))
  for (i in seq_len(nrow(corr))) {
    b <- bins[bins$source == corr$source[i] & bins$cutoff == corr$cutoff[i], ]
    stopifnot(nrow(b) == 15L)
    if (is.finite(corr$nominal_p[i])) {
      check <- cor.test(b$neff, b$diversification, method = "pearson")
      stopifnot(abs(check$p.value - corr$nominal_p[i]) < 1e-10,
                abs(unname(check$estimate) - corr$Pearson_r[i]) < 1e-10)
    } else stopifnot(sd(b$neff) <= 1e-12)
  }
  checks <- rbind(checks, data.frame(trial, independent_neff_values = evaluated,
    maximum_absolute_error = maximum_error, status = "passed"))
}
write.csv(checks, "provenance/metric-audit.csv", row.names = FALSE)
print(checks)
