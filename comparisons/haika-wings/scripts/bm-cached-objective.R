# Exact cached objective and analytic gradient for phytools::multirateBM.
# The model is unchanged from phytools 2.5.2; only its numerical evaluation differs.
#
# Write z_v = log(sigma_v^2) for each tree node. On edge e=(p,c), d=z_c-z_p,
# L_e=exp(z_p)*expm1(d)/d is the logarithmic mean rate, and v_e=l_e*L_e.
# Let A_ie=1 when edge e is ancestral to tip i. Then C=A diag(v) A'.
# The trait likelihood is MVN(x; alpha*1,C), with alpha profiled by GLS.
# For u=C^-1(x-alpha*1), the derivative of its negative log likelihood is
# dNLL/dv_e = 0.5*(a_e' C^-1 a_e - (a_e' u)^2).
# Profiling alpha adds no gradient term because dNLL/dalpha=0 at alpha_hat.
#
# The log-rate penalty in phytools is MVN(z_nonroot; z_root*1, VCV(tree)).
# Because this vector contains EVERY nonroot node, this density factors exactly
# into Brownian edge increments. Its negative log density is
# 0.5*sum_e(log(2*pi*l_e) + (z_c-z_p)^2/l_e).
# This is an exact factorization, not a different prior or an approximation.
#
# dL/dz_p=exp(z_p)*(expm1(d)-d)/d^2,
# dL/dz_c=exp(z_p)*(d+(d-1)*expm1(d))/d^2.
# Taylor evaluation near zero removes cancellation without changing the model.
# Fixed descendants and penalty constants are cached. fn/gr calls at the same
# parameter vector share one Cholesky factorization and covariance calculation.

make_bm_objective <- function(tree, x, lambda) {
  if (!inherits(tree, "phylo")) {
    stop("tree must be a phylo object.")
  }
  if (is.null(names(x)) || anyDuplicated(names(x)) ||
      anyDuplicated(tree$tip.label) || !setequal(names(x), tree$tip.label)) {
    stop("The named trait vector must match unique tree tip labels exactly.")
  }
  x <- x[tree$tip.label]
  if (!is.numeric(x) || any(!is.finite(x))) {
    stop("Trait values must be finite and numeric.")
  }
  if (length(lambda) != 1 || !is.finite(lambda) || lambda < 0) {
    stop("lambda must be one nonnegative finite number.")
  }
  edge_length <- tree$edge.length
  if (length(edge_length) != nrow(tree$edge) ||
      any(!is.finite(edge_length)) || any(edge_length <= 0)) {
    stop("The tree must have finite, strictly positive edge lengths.")
  }

  n_tips <- length(tree$tip.label)
  n_nodes <- n_tips + tree$Nnode
  root_node <- n_tips + 1L
  edge_parent <- tree$edge[, 1]
  edge_child <- tree$edge[, 2]
  n_edges <- nrow(tree$edge)
  node_numbers <- sort(unique(as.integer(tree$edge)))
  if (!identical(node_numbers, seq_len(n_nodes)) ||
      root_node %in% edge_child || anyDuplicated(edge_child) ||
      n_edges != n_nodes - 1L) {
    stop("The rooted tree must use standard ape node numbering.")
  }

  edge_to_child <- integer(n_nodes)
  edge_to_child[edge_child] <- seq_len(n_edges)
  descendants <- matrix(0, nrow = n_tips, ncol = n_edges)
  for (tip_index in seq_len(n_tips)) {
    current_node <- tip_index
    traversed <- 0L
    while (current_node != root_node) {
      edge_index <- edge_to_child[current_node]
      if (edge_index == 0L || traversed >= n_edges) {
        stop("The tree does not form connected paths to its root.")
      }
      descendants[tip_index, edge_index] <- 1
      current_node <- edge_parent[edge_index]
      traversed <- traversed + 1L
    }
  }
  inverse_edge_length <- 1 / edge_length
  penalty_constant <- sum(log(2 * pi * edge_length))
  trait_constant <- n_tips * log(2 * pi)
  cache <- new.env(parent = emptyenv())
  cache$parameter <- NULL
  cache$result <- NULL

  evaluate <- function(parameter) {
    parameter <- as.numeric(parameter)
    if (length(parameter) != n_nodes || any(!is.finite(parameter))) {
      stop("Supply one finite log variance rate for every tree node.")
    }
    if (identical(parameter, cache$parameter)) {
      return(cache$result)
    }

    parent_rate <- exp(parameter[edge_parent])
    difference <- parameter[edge_child] - parameter[edge_parent]
    close <- abs(difference) < 1e-4
    separated <- !close
    mean_factor <- numeric(n_edges)
    parent_factor <- numeric(n_edges)
    child_factor <- numeric(n_edges)

    # Integral of exp(s*d), with s on [0,1], and its endpoint derivatives.
    d <- difference[separated]
    exponential_minus_one <- expm1(d)
    mean_factor[separated] <- exponential_minus_one / d
    parent_factor[separated] <- (exponential_minus_one - d) / d^2
    child_factor[separated] <- (d + (d - 1) * exponential_minus_one) / d^2
    d <- difference[close]
    mean_factor[close] <- 1 + d / 2 + d^2 / 6 + d^3 / 24 +
      d^4 / 120 + d^5 / 720
    parent_factor[close] <- 1 / 2 + d / 6 + d^2 / 24 + d^3 / 120 +
      d^4 / 720 + d^5 / 5040
    child_factor[close] <- 1 / 2 + d / 3 + d^2 / 8 + d^3 / 30 +
      d^4 / 144 + d^5 / 840
    edge_rate <- parent_rate * mean_factor
    edge_variance <- edge_length * edge_rate
    if (any(!is.finite(edge_variance)) || any(edge_variance <= 0)) {
      stop("The proposed rates do not give finite positive branch variances.")
    }

    weighted_descendants <- sweep(descendants, 2, edge_variance, "*")
    covariance <- tcrossprod(weighted_descendants, descendants)
    covariance_chol <- chol(covariance)
    inverse_covariance <- chol2inv(covariance_chol)
    inverse_one <- rowSums(inverse_covariance)
    ancestral_mean <- sum(inverse_one * x) / sum(inverse_one)
    residual <- x - ancestral_mean
    inverse_residual <- as.vector(inverse_covariance %*% residual)
    quadratic <- sum(residual * inverse_residual)
    log_determinant <- 2 * sum(log(diag(covariance_chol)))
    trait_nll <- 0.5 * (trait_constant + log_determinant + quadratic)
    penalty_nll <- 0.5 * (penalty_constant +
                           sum(difference^2 * inverse_edge_length))
    objective <- trait_nll + lambda * penalty_nll

    inverse_descendants <- inverse_covariance %*% descendants
    edge_quadratic <- colSums(descendants * inverse_descendants)
    edge_residual <- as.vector(crossprod(descendants, inverse_residual))
    variance_gradient <- 0.5 * (edge_quadratic - edge_residual^2)
    parent_gradient <- variance_gradient * edge_length * parent_rate * parent_factor
    child_gradient <- variance_gradient * edge_length * parent_rate * child_factor
    increment_gradient <- lambda * difference * inverse_edge_length
    gradient <- numeric(n_nodes)
    for (edge_index in seq_len(n_edges)) {
      parent <- edge_parent[edge_index]
      child <- edge_child[edge_index]
      gradient[parent] <- gradient[parent] + parent_gradient[edge_index] -
        increment_gradient[edge_index]
      gradient[child] <- gradient[child] + child_gradient[edge_index] +
        increment_gradient[edge_index]
    }

    result <- list(value = objective, gradient = gradient,
                   logLik = -trait_nll, penalty_nll = penalty_nll,
                   ancestral_mean = ancestral_mean,
                   edge_sigma2 = edge_rate)
    cache$parameter <- parameter
    cache$result <- result
    result
  }

  objective_function <- function(parameter) {
    evaluate(parameter)$value
  }
  gradient_function <- function(parameter) {
    evaluate(parameter)$gradient
  }
  list(fn = objective_function, gr = gradient_function, evaluate = evaluate,
       n_nodes = n_nodes, n_edges = n_edges, tip_order = tree$tip.label,
       lambda = lambda)
}
