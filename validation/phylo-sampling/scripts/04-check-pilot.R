# Validate the preselected screening trials before the full batch.
source("scripts/load-packages.R")
design <- readRDS("design/design.rds")
characters <- readRDS("design/node-characters.rds")
pilot_tasks <- c(1L, 401L, 501L, 901L, 1001L, 1401L)
dir.create("summary", showWarnings = FALSE)

# Independent binary ER pruning likelihood with the maddfitz root weighting.
er_log_likelihood <- function(rate, tree, states) {
  n_tips <- Ntip(tree)
  conditional <- matrix(1, n_tips + tree$Nnode, 2L)
  conditional[seq_len(n_tips), 1] <- states == "0"
  conditional[seq_len(n_tips), 2] <- states == "1"
  probability_change <- -0.5 * expm1(-2 * rate * tree$edge.length)
  log_scale <- 0
  for (edge in rev(seq_len(nrow(tree$edge)))) {
    parent <- tree$edge[edge, 1]
    child <- tree$edge[edge, 2]
    p <- probability_change[edge]
    message <- c((1-p)*conditional[child,1] + p*conditional[child,2],
                  p*conditional[child,1] + (1-p)*conditional[child,2])
    conditional[parent, ] <- conditional[parent, ] * message
    scale <- sum(conditional[parent, ])
    if (scale <= 0) return(-Inf)
    conditional[parent, ] <- conditional[parent, ] / scale
    log_scale <- log_scale + log(scale)
  }
  root <- conditional[n_tips + 1L, ]
  log_scale + log(sum(root^2) / sum(root))
}
screen_rows <- which(design$nodes$epoch == 3L)
for (epoch in 1:2) {
  rows <- which(design$nodes$epoch == epoch)
  rows <- rows[order(design$nodes$age_ma[rows])]
  screen_rows <- c(screen_rows, rows[round(seq(1, length(rows), length.out = 6))])
}
screen_ids <- design$nodes$character[unique(screen_rows)]
bank_hashes <- tools::md5sum(c("design/node-fits-full.rds", "design/node-fits-study.rds"))
check_file <- "summary/pilot-likelihood-check.rds"
if (file.exists(check_file)) {
  saved <- readRDS(check_file)
  stopifnot(identical(saved$bank_hashes, bank_hashes))
  checks <- saved$checks
} else {
  checks <- data.frame()
  for (tree_name in c("full", "study")) {
    tree <- design[[paste0(tree_name, "_tree")]]
    fits <- readRDS(sprintf("design/node-fits-%s.rds", tree_name))
    stopifnot(!any(fits$status$status == "failed"), !any(fits$status$boundary))
    for (id in screen_ids) {
      fit <- fits$characters[[id]]
      if (fit$status != "ok") next
      states <- characters[[tree_name]][, id]
      rate <- fit$Q_matrix["0", "1"]
      ll <- er_log_likelihood(rate, tree, states)
      backend <- fit$model_table$logLik[fit$model_table$model == "ER"]
      grid <- exp(seq(log(1e-9), log(100), length.out = 80))
      grid_ll <- numeric(length(grid))
      for (i in seq_along(grid)) grid_ll[i] <- er_log_likelihood(grid[i], tree, states)
      best <- which.max(grid_ll)
      bounds <- grid[c(max(1L,best-1L), min(length(grid),best+1L))]
      objective <- function(log_rate) er_log_likelihood(exp(log_rate), tree, states)
      optimum <- optimize(objective, interval = log(bounds), maximum = TRUE, tol = 1e-9)
      checks <- rbind(checks, data.frame(tree = tree_name, character = id,
        fitted_q = rate, independent_q = exp(optimum$maximum),
        likelihood_difference = ll-backend, optimum_improvement = optimum$objective-ll))
    }
  }
  stopifnot(all(abs(checks$likelihood_difference) < 1e-5),
            all(checks$optimum_improvement < 1e-4))
  saveRDS(list(bank_hashes = bank_hashes, checks = checks), check_file)
}
write.csv(checks, "summary/pilot-likelihood-check.csv", row.names = FALSE)
stopifnot(all(abs(checks$likelihood_difference) < 1e-5),
          all(checks$optimum_improvement < 1e-4))
if ("--likelihood-only" %in% commandArgs(trailingOnly = TRUE)) {
  cat("PASS:", nrow(checks), "independent likelihood checks.\n")
  quit(save = "no", status = 0)
}

# Check each simulated pattern and independently marginalize its summed maps.
audit <- data.frame()
for (task in pilot_tasks) {
  result <- readRDS(sprintf("results/task-%04d/result.rds", task))
  dataset <- readRDS(sprintf("simulated/task-%04d.rds", task))
  selected <- dataset$selected_nodes
  stopifnot(!anyDuplicated(selected$node),
            all(selected$epoch == dataset$epoch),
            nrow(selected) == design$conditions$n_characters[dataset$condition],
            result$settings$root_p == "maddfitz", result$settings$n_maps == 100L)
  for (id in selected$character) {
    node <- selected$node[match(id, selected$character)]
    descendants <- extract.clade(design$full_tree, node)$tip.label
    ones <- rownames(characters$full)[characters$full[, id] == "1"]
    stopifnot(setequal(descendants, ones))
  }
  max_error <- 0
  for (group_name in names(result$groups)) {
    group <- result$groups[[group_name]]
    tree_name <- if (group_name == "study_all") "study" else "full"
    tree <- design[[paste0(tree_name, "_tree")]]
    ages <- design[[paste0(tree_name, "_age")]]
    mapped <- readRDS(sprintf("results/task-%04d/maps-%s.rds", task, tree_name))
    stopifnot(!any(mapped$status$status == "failed"))
    total <- matrix(0, nrow(tree$edge), 100L)
    for (id in group$character_ids) {
      one <- read_stm(mapped, characters = id, type = "branch_counts")[[1]]
      total <- total + one
    }
    expected_rates <- total / tree$edge.length
    stopifnot(max(abs(expected_rates-group$rate_samples)) < 1e-9)
    for (i in seq_along(design$age_grid)) {
      age <- design$age_grid[i]
      active <- ages[tree$edge[, 1]] >= age & ages[tree$edge[, 2]] < age
      # At the present use the limiting terminal-branch rates.
      if (age == 0) active <- tree$edge[, 2] <= Ntip(tree)
      expected <- colMeans(expected_rates[active, , drop = FALSE])
      max_error <- max(max_error, abs(expected-group$marginal_samples[i, ]))
    }
    stopifnot(max_error < 1e-9,
              all(group$branch_summary$lower <= group$branch_summary$upper),
              all(group$marginal_summary$lower <= group$marginal_summary$upper))
  }
  audit <- rbind(audit, data.frame(task = task, condition = result$condition,
    characters = result$n_characters, retained = result$n_retained,
    marginal_max_error = max_error, elapsed_seconds = result$elapsed_seconds))
}
write.csv(audit, "summary/pilot-audit.csv", row.names = FALSE)
print(audit)
writeLines(c("PASS: six screening trials validated.",
  paste("Independent likelihood checks:",nrow(checks)),
  paste("Maximum likelihood disagreement:",max(abs(checks$likelihood_difference))),
  paste("Maximum independent likelihood improvement:",max(checks$optimum_improvement)),
  paste("Maximum marginal disagreement:",max(audit$marginal_max_error))),
  "summary/pilot-check.txt")

