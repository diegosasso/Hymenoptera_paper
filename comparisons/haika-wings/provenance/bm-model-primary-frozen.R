# Bounded optimization of the exact phytools multirateBM likelihood.
# Verified against phytools 2.5.2 on 2026-11-20.
# Node rates are converted to branch averages using the logarithmic mean.
source("scripts/bm-cached-objective.R")

bm_logarithmic_mean <- function(parent_rate, child_rate) {
  if (length(parent_rate) != length(child_rate)) {
    stop("Parent and child rate vectors must have the same length.")
  }
  if (any(!is.finite(parent_rate)) || any(!is.finite(child_rate)) ||
      any(parent_rate <= 0) || any(child_rate <= 0)) {
    stop("Node rates must be finite and positive.")
  }

  log_difference <- log(child_rate) - log(parent_rate)
  answer <- parent_rate
  different <- abs(log_difference) > 1e-8
  answer[different] <- parent_rate[different] *
    expm1(log_difference[different]) / log_difference[different]
  close <- !different
  answer[close] <- parent_rate[close] *
    (1 + log_difference[close] / 2 + log_difference[close]^2 / 6)
  answer
}

bm_profile_loglik <- function(tree, x) {
  covariance <- ape::vcv.phylo(tree)
  covariance <- covariance[tree$tip.label, tree$tip.label, drop = FALSE]
  covariance_chol <- chol(covariance)
  inverse_one <- backsolve(covariance_chol,
                          forwardsolve(t(covariance_chol), rep(1, length(x))))
  ancestral_mean <- sum(inverse_one * x) / sum(inverse_one)
  centered_values <- x - ancestral_mean
  whitened_values <- forwardsolve(t(covariance_chol), centered_values)
  log_determinant <- 2 * sum(log(diag(covariance_chol)))
  log_likelihood <- -0.5 * (length(x) * log(2 * pi) + log_determinant +
                            sum(whitened_values^2))
  list(logLik = log_likelihood, ancestral_mean = ancestral_mean,
       quadratic = sum(whitened_values^2),
       log_determinant = log_determinant)
}

fit_bm_constant <- function(tree, x) {
  x <- x[tree$tip.label]
  unit_fit <- bm_profile_loglik(tree, x)
  sigma2 <- unit_fit$quadratic / length(x)
  if (!is.finite(sigma2) || sigma2 <= 0) {
    stop("Constant BM has a nonpositive fitted variance rate.")
  }
  log_likelihood <- -0.5 * (length(x) * log(2 * pi) +
                            length(x) * log(sigma2) +
                            unit_fit$log_determinant + length(x))
  list(sigma2 = sigma2, ancestral_mean = unit_fit$ancestral_mean,
       logLik = log_likelihood, method = "ML", n_parameters = 2L)
}

fit_bm_bounded <- function(tree, x, lambda = 1, n_starts = 3L,
                           seed = 1L, maxit = 2000L) {
  if (!requireNamespace("phytools", quietly = TRUE) ||
      !requireNamespace("ape", quietly = TRUE)) {
    stop("Packages phytools and ape are required.")
  }
  if (!inherits(tree, "phylo")) {
    stop("tree must be a phylo object.")
  }
  if (is.null(names(x)) || anyDuplicated(names(x)) ||
      anyDuplicated(tree$tip.label)) {
    stop("Trait names and tree tip labels must be unique.")
  }
  if (!setequal(names(x), tree$tip.label)) {
    stop("The named trait vector must match the tree tips exactly.")
  }
  x <- x[tree$tip.label]
  if (!is.numeric(x) || any(!is.finite(x)) || var(x) <= 0) {
    stop("Trait values must be finite numeric values with positive variance.")
  }
  if (length(tree$edge.length) != nrow(tree$edge) ||
      any(!is.finite(tree$edge.length)) || any(tree$edge.length <= 0)) {
    stop("This likelihood requires finite, strictly positive edge lengths.")
  }
  if (length(lambda) != 1 || !is.finite(lambda) || lambda <= 0) {
    stop("lambda must be one positive finite number.")
  }
  if (length(n_starts) != 1 || !is.finite(n_starts) ||
      n_starts < 1 || n_starts != as.integer(n_starts)) {
    stop("n_starts must be a positive integer.")
  }
  if (length(maxit) != 1 || !is.finite(maxit) ||
      maxit < 1 || maxit != as.integer(maxit)) {
    stop("maxit must be a positive integer.")
  }

  # Preserve the input topology and order; these are the original phytools bounds.
  n_tips <- length(tree$tip.label)
  n_nodes <- n_tips + tree$Nnode
  root_node <- n_tips + 1L
  if (!ape::is.rooted(tree) ||
      !identical(sort(unique(as.integer(tree$edge))), seq_len(n_nodes)) ||
      root_node %in% tree$edge[, 2]) {
    stop("The rooted tree must use standard ape node numbering.")
  }
  node_names <- c(tree$tip.label, as.character(seq_len(tree$Nnode) + n_tips))
  heights <- phytools::nodeHeights(tree)
  bound_center <- log(var(x) / max(heights))
  lower_bound <- rep(bound_center - 10, n_nodes)
  upper_bound <- rep(bound_center + 10, n_nodes)

  initial_tree <- ape::collapse.singles(tree)
  initial_tree <- ape::multi2di(initial_tree, random = FALSE)
  contrasts <- ape::pic(x, initial_tree)
  initial_rate <- mean(contrasts^2) * (n_tips - 1) / n_tips
  initial_log_rate <- log(initial_rate)
  exact_objective <- getFromNamespace("log_lik", "phytools")
  cached_objective <- make_bm_objective(tree, x, lambda)

  # Fixed numbers of independent starts avoid phytools' unbounded retry loop.
  set.seed(seed)
  fits <- vector("list", n_starts)
  diagnostics <- data.frame(start = seq_len(n_starts), objective = NA_real_,
                            convergence = NA_integer_, function_evaluations = NA_integer_,
                            gradient_evaluations = NA_integer_, elapsed_seconds = NA_real_,
                            n_near_lower_bound = NA_integer_,
                            n_near_upper_bound = NA_integer_,
                            message = rep("", n_starts), stringsAsFactors = FALSE)
  for (start_index in seq_len(n_starts)) {
    initial_values <- rep(initial_log_rate, n_nodes)
    if (start_index > 1) {
      initial_values <- initial_values + rnorm(n_nodes, mean = 0, sd = 0.5)
    }
    initial_values <- pmax(lower_bound, pmin(upper_bound, initial_values))
    time_started <- proc.time()[["elapsed"]]
    result <- tryCatch(
      optim(initial_values, cached_objective$fn, gr = cached_objective$gr,
            method = "L-BFGS-B",
            lower = lower_bound, upper = upper_bound,
            control = list(maxit = maxit)),
      error = function(error) error
    )
    diagnostics$elapsed_seconds[start_index] <-
      proc.time()[["elapsed"]] - time_started
    if (inherits(result, "error")) {
      diagnostics$message[start_index] <- conditionMessage(result)
    } else {
      fits[[start_index]] <- result
      diagnostics$objective[start_index] <- result$value
      diagnostics$convergence[start_index] <- result$convergence
      diagnostics$function_evaluations[start_index] <- result$counts[["function"]]
      diagnostics$gradient_evaluations[start_index] <- result$counts[["gradient"]]
      diagnostics$n_near_lower_bound[start_index] <-
        sum(result$par - lower_bound < 1e-4)
      diagnostics$n_near_upper_bound[start_index] <-
        sum(upper_bound - result$par < 1e-4)
      if (!is.null(result$message)) {
        diagnostics$message[start_index] <- result$message
      }
    }
  }

  finite_starts <- which(is.finite(diagnostics$objective))
  if (length(finite_starts) == 0) {
    return(list(success = FALSE, convergence = NA_integer_,
                diagnostics = diagnostics, lambda = lambda, seed = seed,
                n_starts = n_starts, maxit = maxit,
                phytools_version = as.character(utils::packageVersion("phytools")),
                error = "All optimization starts failed."))
  }
  best_start <- finite_starts[which.min(diagnostics$objective[finite_starts])]
  best_fit <- fits[[best_start]]
  node_sigma2 <- exp(best_fit$par)
  names(node_sigma2) <- node_names
  parent_rate <- node_sigma2[tree$edge[, 1]]
  child_rate <- node_sigma2[tree$edge[, 2]]
  edge_sigma2 <- bm_logarithmic_mean(parent_rate, child_rate)
  names(edge_sigma2) <- paste(tree$edge[, 1], tree$edge[, 2], sep = "->")

  # Check the exported branch averages against the likelihood used to optimize.
  log_likelihood <- -exact_objective(best_fit$par, tree, x, lambda = 0)
  penalized_objective <- exact_objective(best_fit$par, tree, x, lambda = lambda)
  scaled_tree <- tree
  scaled_tree$edge.length <- tree$edge.length * edge_sigma2
  stable_check <- bm_profile_loglik(scaled_tree, x)
  likelihood_difference <- stable_check$logLik - log_likelihood
  likelihood_tolerance <- 1e-6 * (1 + abs(log_likelihood))
  likelihood_check_passed <- abs(likelihood_difference) <= likelihood_tolerance
  objective_difference <- penalized_objective - best_fit$value
  objective_check_passed <- abs(objective_difference) <=
    1e-8 * (1 + abs(penalized_objective))
  constant_fit <- fit_bm_constant(tree, x)

  list(success = best_fit$convergence == 0 && likelihood_check_passed &&
         objective_check_passed,
       convergence = best_fit$convergence,
       node_sigma2 = node_sigma2, edge_sigma2 = edge_sigma2,
       logLik = log_likelihood, penalized_objective = penalized_objective,
       ancestral_mean = stable_check$ancestral_mean,
       lambda = lambda, method = "ML", best_start = best_start,
       diagnostics = diagnostics,
       likelihood_check = data.frame(internal_logLik = log_likelihood,
                                     stable_edge_logLik = stable_check$logLik,
                                     difference = likelihood_difference,
                                     tolerance = likelihood_tolerance,
                                     passed = likelihood_check_passed),
       objective_check = data.frame(optimizer_value = best_fit$value,
                                    recomputed_value = penalized_objective,
                                    difference = objective_difference,
                                    passed = objective_check_passed),
       constant_BM = constant_fit,
       lower_log_bound = lower_bound[1], upper_log_bound = upper_bound[1],
       seed = seed, n_starts = n_starts, maxit = maxit,
       phytools_version = as.character(utils::packageVersion("phytools")),
       objective_implementation = "Exact cached phytools likelihood with analytic gradient",
       tree = tree, x = x)
}
