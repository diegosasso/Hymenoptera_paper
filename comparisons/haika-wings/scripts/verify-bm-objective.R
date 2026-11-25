# Reproducible, small synthetic checks of the cached multirateBM objective.
# Run from comparisons/haika-wings. No HAIKA or ontophylo data are fitted.
project_library <- Sys.getenv("HAIKA_RLIB",
                              "/projappl/project_2004996/haika-wings/Rlib")
.libPaths(c(project_library,
            "/projappl/project_2004996/phylo-sampling/Rlib",
            "/projappl/project_2004996/hymenoptera-node-pattern-bias/Rlib",
            .libPaths()))
library(ape)
library(phytools)
source("scripts/bm-cached-objective.R")

set.seed(6161)
verification_rows <- list()
row_index <- 0L
for (n_tips in c(6L, 8L, 10L)) {
  example_tree <- ape::rtree(n_tips)
  example_tree$edge.length <- example_tree$edge.length + 0.1
  example_x <- phytools::fastBM(example_tree, sig2 = 0.8)
  n_nodes <- n_tips + example_tree$Nnode
  parameter_sets <- list(
    random = rnorm(n_nodes, -0.5, 0.7),
    equal = rep(log(0.7), n_nodes),
    nearly_equal = rep(log(0.7), n_nodes) + rnorm(n_nodes, 0, 1e-7)
  )
  for (lambda in c(0, 0.1, 1, 10)) {
    objective <- make_bm_objective(example_tree, example_x, lambda)
    for (parameter_case in names(parameter_sets)) {
      parameters <- parameter_sets[[parameter_case]]
      actual_value <- objective$fn(parameters)
      reference_value <- phytools:::log_lik(parameters, example_tree,
                                             example_x, lambda)
      # The original logarithmic mean loses some precision near equal rates.
      objective_error <- abs(actual_value - reference_value)
      analytic_gradient <- objective$gr(parameters)
      numerical_gradient <- numeric(n_nodes)
      step <- 1e-6
      for (node_index in seq_len(n_nodes)) {
        upper <- parameters
        lower <- parameters
        upper[node_index] <- upper[node_index] + step
        lower[node_index] <- lower[node_index] - step
        numerical_gradient[node_index] <-
          (objective$fn(upper) - objective$fn(lower)) / (2 * step)
      }
      gradient_error <- max(abs(analytic_gradient - numerical_gradient))
      objective_tolerance <- 1e-7 * (1 + abs(reference_value))
      gradient_tolerance <- 1e-6 * (1 + max(abs(analytic_gradient)))
      row_index <- row_index + 1L
      verification_rows[[row_index]] <- data.frame(
        check = "objective_and_gradient", n_tips = n_tips,
        lambda = lambda, parameter_case = parameter_case,
        objective_cached = actual_value, objective_reference = reference_value,
        objective_absolute_error = objective_error,
        objective_tolerance = objective_tolerance,
        gradient_maximum_absolute_error = gradient_error,
        gradient_tolerance = gradient_tolerance,
        optimum_parameter_maximum_absolute_error = NA_real_,
        cached_convergence = NA_integer_, reference_convergence = NA_integer_,
        cached_elapsed_seconds = NA_real_, reference_elapsed_seconds = NA_real_,
        passed = objective_error < objective_tolerance &&
          gradient_error < gradient_tolerance,
        stringsAsFactors = FALSE
      )
    }
  }
}

# One small optimization comparison verifies the gradient's practical use.
set.seed(6162)
example_tree <- ape::rtree(8)
example_tree$edge.length <- example_tree$edge.length + 0.2
example_x <- phytools::fastBM(example_tree, sig2 = 0.8)
n_tips <- length(example_tree$tip.label)
n_nodes <- n_tips + example_tree$Nnode
initial_rate <- mean(ape::pic(example_x, example_tree)^2) * (n_tips - 1) / n_tips
initial_values <- rep(log(initial_rate), n_nodes)
bound_center <- log(var(example_x) / max(phytools::nodeHeights(example_tree)))
lower_bound <- rep(bound_center - 10, n_nodes)
upper_bound <- rep(bound_center + 10, n_nodes)
objective <- make_bm_objective(example_tree, example_x, lambda = 1)
cached_time <- system.time(
  cached_fit <- optim(initial_values, objective$fn, gr = objective$gr,
                      method = "L-BFGS-B", lower = lower_bound, upper = upper_bound,
                      control = list(maxit = 1000))
)
reference_time <- system.time(
  reference_fit <- optim(initial_values, phytools:::log_lik,
                         tree = example_tree, x = example_x, lambda = 1,
                         method = "L-BFGS-B", lower = lower_bound, upper = upper_bound,
                         control = list(maxit = 1000))
)
optimum_objective_error <- abs(cached_fit$value - reference_fit$value)
optimum_parameter_error <- max(abs(cached_fit$par - reference_fit$par))
row_index <- row_index + 1L
verification_rows[[row_index]] <- data.frame(
  check = "optimizer_comparison", n_tips = n_tips,
  lambda = 1, parameter_case = "constant_PIC_start",
  objective_cached = cached_fit$value, objective_reference = reference_fit$value,
  objective_absolute_error = optimum_objective_error,
  objective_tolerance = 1e-6,
  gradient_maximum_absolute_error = NA_real_, gradient_tolerance = NA_real_,
  optimum_parameter_maximum_absolute_error = optimum_parameter_error,
  cached_convergence = cached_fit$convergence,
  reference_convergence = reference_fit$convergence,
  cached_elapsed_seconds = cached_time[["elapsed"]],
  reference_elapsed_seconds = reference_time[["elapsed"]],
  passed = cached_fit$convergence == 0 && reference_fit$convergence == 0 &&
    optimum_objective_error < 1e-6 && optimum_parameter_error < 1e-4,
  stringsAsFactors = FALSE
)

verification <- do.call(rbind, verification_rows)
verification$phytools_version <- as.character(utils::packageVersion("phytools"))
verification$ape_version <- as.character(utils::packageVersion("ape"))
verification$R_version <- R.version.string
verification$verified_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
verification$objective_md5 <- unname(tools::md5sum("scripts/bm-cached-objective.R"))
verification$verification_script_md5 <-
  unname(tools::md5sum("scripts/verify-bm-objective.R"))
dir.create("provenance", showWarnings = FALSE)
write.csv(verification, "provenance/bm-objective-verification.csv", row.names = FALSE)
cat("Checks passed:", sum(verification$passed), "of", nrow(verification), "\n")
cat("Maximum objective error:", max(verification$objective_absolute_error), "\n")
cat("Maximum gradient error:",
    max(verification$gradient_maximum_absolute_error, na.rm = TRUE), "\n")
cat("Optimizer objective error:", optimum_objective_error, "\n")
cat("Cached/reference optimization seconds:", cached_time[["elapsed"]],
    reference_time[["elapsed"]], "\n")
stopifnot(all(verification$passed))
