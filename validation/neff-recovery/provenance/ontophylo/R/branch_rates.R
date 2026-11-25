# Branch-specific reconstruction adapted from Sergei Tarasov's functions in
# Hymenoptera_paper/R/utils-branch-ama.R and R/utils-ST.R.

#' Count mapped changes on each branch
#'
#' Count state changes in an original, single-character stochastic map.
#' Adjacent segments with the same state are not counted as changes.
#'
#' @param simmap A stochastic map produced by [phytools::make.simmap()].
#' @return A numeric vector in the order of `simmap$edge`.
#' @details Use the original maps, before discretization or PARAMO stacking.
#' Discretized maps can omit transitions; counting stacked-state transitions
#' does not in general count all component-character changes.
#' @author Sergei Tarasov
#' @examples
#' data("hym_stm")
#' head(count_branch_changes(hym_stm[[1]][[1]]))
#' @md
#' @export
count_branch_changes <- function(simmap) {
  maps <- .check_branch_maps(simmap)
  if (length(maps) != 1L) {
    stop("Supply one stochastic map.", call. = FALSE)
  }
  .count_branch_changes(maps[[1]])
}

#' Extract branch-rate samples for one character
#'
#' Divide mapped change counts by branch duration separately for every map.
#'
#' @param smaps A list of original stochastic maps for one character, or a
#'   single `simmap`. All maps must have the same tree and edge order.
#' @return A list of numeric vectors, one per map replicate, in tree edge order.
#' @seealso [count_branch_changes()], [ama_branches_multi_char()]
#' @author Sergei Tarasov
#' @examples
#' data("hym_stm")
#' maps <- hym_stm[[1]][1:3]
#' samples <- get_branch_rate_sample_smaps(maps)
#' head(samples[[1]])
#' @md
#' @export
get_branch_rate_sample_smaps <- function(smaps) {
  smaps <- .check_branch_maps(smaps)
  branch.length <- smaps[[1]]$edge.length
  rate.samples <- vector("list", length(smaps))
  for (i in seq_along(smaps)) {
    changes <- .count_branch_changes(smaps[[i]])
    rate.samples[[i]] <- changes / branch.length
  }
  rate.samples
}

#' Amalgamate branch-specific rates across characters
#'
#' Sum character-specific mapped rates within each branch and map replicate.
#' This implements branch-specific amalgamation without discretizing maps or
#' constructing a stacked-state stochastic map.
#'
#' @param list.of.chars A nonempty list of characters. Each character contains
#'   an equally sized list of original stochastic maps (or one `simmap`).
#'   Trees, edge order and branch lengths must agree across all maps.
#' @return A list with one numeric branch-rate vector per map replicate.
#' @details For branch b and replicate m, the combined rate is the sum of
#'   mapped changes across characters divided by branch duration. Its units
#'   are total character changes per tree time unit. Divide by the number of
#'   characters for a mean rate per character. Replicate m is paired across
#'   characters: retain joint draw identities when histories are dependent.
#'   For independent character samples, this pairing samples their product
#'   distribution. Do not include the same character twice.
#'
#'   These are mapped change rates conditional on the supplied histories,
#'   not separately fitted branch-specific Q matrices.
#' @author Sergei Tarasov
#' @examples
#' data("hym_stm")
#' maps <- hym_stm[1:2]
#' samples <- ama_branches_multi_char(maps)
#' head(samples[[1]])
#' @md
#' @export
ama_branches_multi_char <- function(list.of.chars) {
  .check_character_list(list.of.chars)
  first.maps <- .check_branch_maps(list.of.chars[[1]])
  tree <- first.maps[[1]]
  n_maps <- length(first.maps)
  summed_across_chars <- vector("list", n_maps)
  for (i in seq_len(n_maps)) {
    summed_across_chars[[i]] <- numeric(nrow(tree$edge))
  }

  for (j in seq_along(list.of.chars)) {
    char.maps <- .check_branch_maps(list.of.chars[[j]], reference = tree)
    if (length(char.maps) != n_maps) {
      stop("All characters must have the same number of map replicates.",
           call. = FALSE)
    }
    for (i in seq_len(n_maps)) {
      changes <- .count_branch_changes(char.maps[[i]])
      character.rate <- changes / tree$edge.length
      summed_across_chars[[i]] <- summed_across_chars[[i]] + character.rate
    }
  }
  summed_across_chars
}

#' Summarize branch rates across stochastic maps
#'
#' Compute branch-wise statistics across map replicates for either one
#' character or an amalgamated set of characters.
#'
#' @param multi_simmap A list of original stochastic maps for one character,
#'   or a single `simmap`.
#' @param stat Statistic: `"mean"`, `"median"`, `"sd"`, or `"quantile"`.
#' @param probs Quantile probabilities, used when `stat = "quantile"`.
#' @return A vector in edge order, or a matrix with quantiles in rows and
#'   branches in columns. SD is `NA` when there is only one map.
#' @details Quantiles describe the distribution of reconstructed mapped rates,
#'   conditional on the tree, model and observations used to sample histories.
#'   With fixed Q and tree, they do not include parameter or tree uncertainty.
#' @author Sergei Tarasov
#' @examples
#' data("hym_stm")
#' head(summarize_branch_rates(hym_stm[[1]], stat = "mean"))
#' head(summarize_branch_rates_multi(hym_stm[1:2], stat = "mean"))
#' @md
#' @export
summarize_branch_rates <- function(multi_simmap,
                                   stat = c("mean", "median", "sd", "quantile"),
                                   probs = c(0.025, 0.5, 0.975)) {
  stat <- match.arg(stat)
  samples <- get_branch_rate_sample_smaps(multi_simmap)
  mat <- .branch_sample_matrix(samples)
  .branch_statistic(mat, stat, probs)
}

#' @rdname summarize_branch_rates
#' @inheritParams ama_branches_multi_char
#' @md
#' @export
summarize_branch_rates_multi <- function(list.of.chars,
                                         stat = c("mean", "median", "sd", "quantile"),
                                         probs = c(0.025, 0.5, 0.975)) {
  stat <- match.arg(stat)
  samples <- ama_branches_multi_char(list.of.chars)
  mat <- .branch_sample_matrix(samples)
  .branch_statistic(mat, stat, probs)
}

#' Marginalize branch rates through time
#'
#' Average branch-constant rates across sampled lineages present at each time.
#'
#' @param tree A rooted ultrametric `phylo` with positive branch lengths.
#' @param branch_rate A finite, nonnegative vector in `tree$edge` order.
#' @param n_points Number of equally spaced times when `time_points` is NULL.
#' @param time_points Optional times in the units of the tree. The supplied
#'   order is retained. Their interpretation is controlled by `time_scale`.
#' @param time_scale `"root"` for elapsed time from the root (the original
#'   interface), or `"age"` for time before the present. Returned `time` is
#'   always age before the present.
#' @return A data frame with `time`, `rate`, and `n_branches`.
#' @details A branch contributes from its parent time inclusive to its child
#'   time exclusive. At a node, descendant branches contribute instead of the
#'   ancestral branch. At the present, terminal branches contribute their
#'   limiting rates. The root stem, if supplied, is not included.
#'   This is a mean over sampled lineages, not the sum of their rates.
#' @author Sergei Tarasov
#' @examples
#' tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
#' branch_rate_thru_time(tree, c(1, 2, 3, 1),
#'                        time_points = c(2, 1, 0), time_scale = "age")
#' @md
#' @export
branch_rate_thru_time <- function(tree, branch_rate, n_points = 100,
                                  time_points = NULL,
                                  time_scale = c("root", "age")) {
  slices <- .branch_time_slices(tree, n_points, time_points, time_scale)
  .check_branch_vector(branch_rate, nrow(tree$edge))
  rate_vals <- numeric(length(slices$edges))
  for (i in seq_along(slices$edges)) {
    rate_vals[i] <- mean(branch_rate[slices$edges[[i]]])
  }
  data.frame(time = slices$time, rate = rate_vals,
             n_branches = lengths(slices$edges))
}

#' Summarize marginal rates across map replicates
#'
#' Marginalize each map replicate before calculating statistics across maps.
#'
#' @inheritParams branch_rate_thru_time
#' @param replicate_list A list of branch-rate vectors, as returned by
#'   [ama_branches_multi_char()]. For compatibility with the original workflow,
#'   each element can instead be a branches-by-groups matrix; its columns are
#'   summed first. Those groups must not contain overlapping characters.
#' @param return_samples If TRUE, return both the time-by-replicate `samples`
#'   matrix and the `summary` data frame. Otherwise return the summary alone.
#' @return A summary with age (`time`), number of contributing branches,
#'   `mean`, `median`, `sd`, `lower` and `upper`, or a list with `samples` and
#'   `summary`. Bounds are pointwise 2.5% and 97.5% map quantiles.
#' @seealso [branch_rates()]
#' @author Sergei Tarasov
#' @examples
#' tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
#' samples <- list(c(1, 2, 3, 1), c(1, 3, 2, 1))
#' curves <- branch_rate_thru_time_Multi(tree, samples, n_points = 5,
#'                                       return_samples = TRUE)
#' curves$summary
#' @md
#' @export
branch_rate_thru_time_Multi <- function(tree, replicate_list, n_points = 100,
                                        time_points = NULL,
                                        time_scale = c("root", "age"),
                                        return_samples = FALSE) {
  time_scale <- match.arg(time_scale)
  if (!is.list(replicate_list) || length(replicate_list) == 0L) {
    stop("replicate_list must be a nonempty list of rate samples.", call. = FALSE)
  }
  if (!is.logical(return_samples) || length(return_samples) != 1L ||
      is.na(return_samples)) {
    stop("return_samples must be TRUE or FALSE.", call. = FALSE)
  }
  n_maps <- length(replicate_list)
  for (i in seq_len(n_maps)) {
    branch_rate <- replicate_list[[i]]
    if (is.matrix(branch_rate)) {
      .check_rate_matrix(branch_rate)
      branch_rate <- rowSums(branch_rate)
    }
    one_curve <- branch_rate_thru_time(tree, branch_rate, n_points,
                                       time_points, time_scale)
    if (i == 1L) {
      samples <- matrix(NA_real_, nrow(one_curve), n_maps)
      times <- one_curve$time
      n_branches <- one_curve$n_branches
    }
    samples[, i] <- one_curve$rate
  }
  colnames(samples) <- paste0("map_", seq_len(n_maps))
  summary <- data.frame(time = times, n_branches = n_branches,
                         .branch_summary(samples))
  if (return_samples) {
    return(list(samples = samples, summary = summary))
  }
  summary
}

#' Reconstruct combined-character branch and marginal rates
#'
#' A wrapper around branch-specific amalgamation and its summaries. Original
#' character histories are combined by summing mapped rates on each branch
#' within each replicate. No map discretization or NHPP fitting is required.
#'
#' @inheritParams ama_branches_multi_char
#' @inheritParams branch_rate_thru_time
#' @param list.of.chars A named list of original character maps, or an
#'   `ontophylo_stm` result from [paramo_stm()]. The latter reads one character's
#'   counts at a time and also supports runs storing only branch counts.
#' @param marginal Also calculate marginal rates through time. Set FALSE for
#'   branch-only reconstruction, including on non-ultrametric trees.
#' @return A list containing `tree`, `characters`, `n_characters`, `n_maps`,
#'   `rate_samples` (branches by maps), `branch_summary` (edge IDs, parent and
#'   child nodes, branch lengths, mean, median, SD and 95% map quantiles),
#'   `marginal_samples` (times by maps) and `marginal_summary`. The last two
#'   elements are NULL when `marginal = FALSE`. Branch order matches the input.
#' @details This wrapper preserves map replicates through both calculations.
#'   Intervals represent conditional mapping uncertainty, not automatically
#'   uncertainty in fitted Q matrices, the phylogeny or taxon sampling.
#'   Select one unique set of characters to define the combined character.
#'   For ontology groups, subset the input using a result from [RAC_query()].
#' @examples
#' tree <- ape::read.tree(text = "(a:1,b:1);")
#' Q <- matrix(c(-1, 1, 1, -1), 2, 2)
#' dimnames(Q) <- list(c("0", "1"), c("0", "1"))
#' x <- c(a = "0", b = "1")
#' set.seed(1)
#' maps1 <- phytools::make.simmap(tree, x, Q = Q, pi = c(0.5, 0.5),
#'                               nsim = 5, message = FALSE)
#' maps2 <- phytools::make.simmap(tree, x, Q = Q, pi = c(0.5, 0.5),
#'                               nsim = 5, message = FALSE)
#' characters <- list(char1 = maps1, char2 = maps2)
#' result <- branch_rates(characters, time_points = c(1, 0.5, 0),
#'                         time_scale = "age")
#' result$branch_summary
#' result$marginal_summary
#' @md
#' @export
branch_rates <- function(list.of.chars, n_points = 100, time_points = NULL,
                          time_scale = c("root", "age"), marginal = TRUE) {
  time_scale <- match.arg(time_scale)
  if (!is.logical(marginal) || length(marginal) != 1L || is.na(marginal)) {
    stop("marginal must be TRUE or FALSE.", call. = FALSE)
  }
  if (inherits(list.of.chars, "ontophylo_stm")) {
    first.tree <- list.of.chars$tree
    characters <- names(list.of.chars$characters)
    total <- matrix(0, nrow(first.tree$edge), list.of.chars$settings$n_stm)
    for (id in characters) {
      one <- read_stm(list.of.chars, characters = id, type = "branch_counts")[[1]]
      if (!identical(dim(one), dim(total)) || any(!is.finite(one)) ||
          any(one < 0) || any(one != floor(one))) {
        stop("Invalid branch-count samples.", call. = FALSE)
      }
      total <- total + one
    }
    samples <- vector("list", ncol(total))
    for (i in seq_along(samples)) samples[[i]] <- total[, i] / first.tree$edge.length
  } else {
    samples <- ama_branches_multi_char(list.of.chars)
    first.maps <- .check_branch_maps(list.of.chars[[1]])
    first.tree <- first.maps[[1]]
    characters <- names(list.of.chars)
  }
  tree <- list(edge = first.tree$edge, edge.length = first.tree$edge.length,
               tip.label = first.tree$tip.label, Nnode = first.tree$Nnode)
  class(tree) <- "phylo"
  if (is.null(characters)) {
    characters <- paste0("char_", seq_along(list.of.chars))
  }
  rate_samples <- .branch_sample_matrix(samples)
  rownames(rate_samples) <- paste0("branch_", seq_len(nrow(tree$edge)))
  colnames(rate_samples) <- paste0("map_", seq_along(samples))
  branch_summary <- data.frame(edge = seq_len(nrow(tree$edge)),
                                parent = tree$edge[, 1], child = tree$edge[, 2],
                                branch_length = tree$edge.length,
                                .branch_summary(rate_samples), row.names = NULL)
  marginal_samples <- NULL
  marginal_summary <- NULL
  if (marginal) {
    marginal_result <- branch_rate_thru_time_Multi(
      tree, samples, n_points, time_points, time_scale, return_samples = TRUE)
    marginal_samples <- marginal_result$samples
    marginal_summary <- marginal_result$summary
  }
  list(tree = tree, characters = characters,
       n_characters = length(characters), n_maps = length(samples),
       rate_samples = rate_samples, branch_summary = branch_summary,
       marginal_samples = marginal_samples, marginal_summary = marginal_summary)
}

# Internal validation and shared summaries.
.check_branch_tree <- function(tree) {
  if (!inherits(tree, "phylo") || !is.matrix(tree$edge) ||
      ncol(tree$edge) != 2L || nrow(tree$edge) == 0L ||
      !is.numeric(tree$edge) || any(!is.finite(tree$edge))) {
    stop("Supply a phylo tree with a valid edge matrix.", call. = FALSE)
  }
  if (!is.numeric(tree$edge.length) ||
      length(tree$edge.length) != nrow(tree$edge) ||
      any(!is.finite(tree$edge.length)) || any(tree$edge.length <= 0)) {
    stop("All branches must have finite, positive lengths.", call. = FALSE)
  }
}

.check_character_list <- function(list.of.chars) {
  if (!is.list(list.of.chars) || length(list.of.chars) == 0L ||
      inherits(list.of.chars, "phylo")) {
    stop("list.of.chars must be a nonempty list of characters.", call. = FALSE)
  }
  ids <- names(list.of.chars)
  if (!is.null(ids) && (anyNA(ids) || any(ids == "") || anyDuplicated(ids))) {
    stop("Character names must be nonempty and unique.", call. = FALSE)
  }
}

.check_branch_maps <- function(smaps, reference = NULL) {
  if (inherits(smaps, "simmap")) {
    smaps <- list(smaps)
  }
  if (!is.list(smaps) || length(smaps) == 0L) {
    stop("Supply a nonempty list of original stochastic maps.", call. = FALSE)
  }
  if (is.null(reference)) {
    reference <- smaps[[1]]
  }
  .check_branch_tree(reference)
  for (i in seq_along(smaps)) {
    tree <- smaps[[i]]
    .check_branch_tree(tree)
    same.edges <- identical(dim(tree$edge), dim(reference$edge)) &&
      all(tree$edge == reference$edge)
    same.lengths <- isTRUE(all.equal(tree$edge.length, reference$edge.length,
                                     tolerance = 1e-8, check.attributes = FALSE))
    if (!same.edges || !identical(tree$tip.label, reference$tip.label) || !same.lengths) {
      stop("All maps must share the same tree, tip labels, edge order and branch lengths.",
           call. = FALSE)
    }
    if (!inherits(tree, "simmap") || !is.list(tree$maps) ||
        length(tree$maps) != nrow(tree$edge)) {
      stop("Each stochastic map must have one mapped history per branch.", call. = FALSE)
    }
    for (b in seq_along(tree$maps)) {
      segments <- tree$maps[[b]]
      states <- names(segments)
      if (!is.numeric(segments) || length(segments) == 0L ||
          any(!is.finite(segments)) || any(segments <= 0) ||
          is.null(states) || anyNA(states) || any(states == "")) {
        stop("Mapped segments must have state labels and finite, positive durations.",
             call. = FALSE)
      }
      if (abs(sum(segments) - tree$edge.length[b]) >
          1e-8 * max(1, tree$edge.length[b])) {
        stop("Mapped segment durations must sum to branch lengths.", call. = FALSE)
      }
    }
  }
  smaps
}

.count_branch_changes <- function(simmap) {
  changes <- numeric(length(simmap$maps))
  for (b in seq_along(simmap$maps)) {
    states <- names(simmap$maps[[b]])
    n_segments <- length(states)
    if (n_segments > 1L) {
      changes[b] <- sum(states[-1L] != states[-n_segments])
    }
  }
  changes
}

.check_branch_vector <- function(x, n_branches) {
  if (!is.numeric(x) || !is.null(dim(x)) || length(x) != n_branches ||
      any(!is.finite(x)) || any(x < 0)) {
    stop("Supply one finite, nonnegative rate per branch, in tree edge order.",
         call. = FALSE)
  }
}

.check_rate_matrix <- function(mat) {
  if (!is.matrix(mat) || !is.numeric(mat) || any(dim(mat) == 0L) ||
      any(!is.finite(mat)) || any(mat < 0)) {
    stop("Rate samples must be a nonempty, finite, nonnegative numeric matrix.",
         call. = FALSE)
  }
}

.branch_sample_matrix <- function(samples) {
  n_branches <- length(samples[[1]])
  mat <- matrix(NA_real_, n_branches, length(samples))
  for (i in seq_along(samples)) {
    .check_branch_vector(samples[[i]], n_branches)
    mat[, i] <- samples[[i]]
  }
  mat
}

.branch_statistic <- function(mat, stat, probs) {
  if (stat == "quantile") {
    if (!is.numeric(probs) || length(probs) == 0L || any(!is.finite(probs)) ||
        any(probs < 0 | probs > 1)) {
      stop("probs must contain probabilities between zero and one.", call. = FALSE)
    }
    out <- matrix(NA_real_, length(probs), nrow(mat))
    for (b in seq_len(nrow(mat))) {
      out[, b] <- stats::quantile(mat[b, ], probs = probs, names = FALSE)
    }
    rownames(out) <- names(stats::quantile(mat[1, ], probs = probs))
    return(out)
  }
  out <- numeric(nrow(mat))
  for (b in seq_len(nrow(mat))) {
    if (stat == "mean") out[b] <- mean(mat[b, ])
    if (stat == "median") out[b] <- stats::median(mat[b, ])
    if (stat == "sd") out[b] <- stats::sd(mat[b, ])
  }
  out
}

.branch_summary <- function(mat) {
  bounds <- .branch_statistic(mat, "quantile", c(0.025, 0.975))
  data.frame(mean = .branch_statistic(mat, "mean", NULL),
             median = .branch_statistic(mat, "median", NULL),
             sd = .branch_statistic(mat, "sd", NULL),
             lower = bounds[1, ], upper = bounds[2, ])
}

# Shared node-boundary convention for marginal rates and phenomic metrics.
.branch_time_slices <- function(tree, n_points, time_points,
                                 time_scale = c("root", "age")) {
  time_scale <- match.arg(time_scale)
  .check_branch_tree(tree)
  if (!ape::is.rooted(tree) || !ape::is.ultrametric(tree)) {
    stop("Marginal rates require a rooted ultrametric tree.", call. = FALSE)
  }
  heights <- phytools::nodeHeights(tree)
  tree_height <- max(heights)

  if (is.null(time_points)) {
    if (!is.numeric(n_points) || length(n_points) != 1L ||
        !is.finite(n_points) || n_points < 2 || n_points != as.integer(n_points)) {
      stop("n_points must be an integer of at least two.", call. = FALSE)
    }
    times <- seq(0, tree_height, length.out = n_points)
  } else {
    if (!is.numeric(time_points) || length(time_points) == 0L ||
        any(!is.finite(time_points)) ||
        any(time_points < 0 | time_points > tree_height)) {
      stop("time_points must lie between zero and the tree height.", call. = FALSE)
    }
    times <- time_points
    if (time_scale == "age") {
      times <- tree_height - time_points
    }
  }

  edge_sets <- vector("list", length(times))
  for (i in seq_along(times)) {
    t <- times[i]
    edge_ids <- which(t >= heights[, 1] & t < heights[, 2])
    if (t == tree_height) {
      # All tips represent the present, allowing numerical ultrametric error.
      edge_ids <- which(tree$edge[, 2] <= length(tree$tip.label))
    }
    edge_sets[[i]] <- edge_ids
  }
  list(time = tree_height - times, edges = edge_sets)
}
