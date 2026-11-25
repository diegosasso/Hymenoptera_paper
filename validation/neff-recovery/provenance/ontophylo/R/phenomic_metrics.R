#' Identify unusually high or low branch rates within each region
#'
#' Port of the Hymenoptera `get_enrichment_mask()` calculation. Quantiles are
#' calculated across branches separately for each column of one map replicate.
#'
#' @param rate_matrix Nonempty numeric matrix: branches by anatomical regions.
#'   Values must be finite and nonnegative; missing values are not discarded.
#' @param lower_q,upper_q Lower and upper quantile probabilities, with
#'   `0 <= lower_q < upper_q <= 1`. R's type-7 quantiles are used.
#' @param direction `"over"` for acceleration, `"under"` for deceleration,
#'   or `"both"` for their union.
#' @param return_type `"binary"` for a 0/1 mask, or `"masked"` to retain
#'   selected rates and set other entries to zero. A masked zero cannot
#'   distinguish a selected zero rate from an unselected branch: retain the mask.
#' @param ties `"include"` includes rates equal to the cutoff, matching the
#'   original directional analyses. `"exclude"` uses strict inequalities.
#' @param constant `"unclassified"` assigns no regime to a column with identical
#'   rates. `"include"` applies the cutoff even to constant columns; together
#'   with `ties = "include"` this reproduces the old directional behaviour.
#' @return A matrix with the input dimensions and dimension names.
#' @details These are relative tail classifications, not estimates of a temporal
#'   derivative or tests against a fitted homogeneous-rate null. Ties can select
#'   substantially more branches than the nominal tail probability; zero-rate
#'   ties can be common for deceleration. The two masks can overlap when the
#'   lower and upper cutoffs coincide. `"both"` is always their union (the old
#'   helper inconsistently used strict inequalities only for `"both"`).
#' @examples
#' rates <- cbind(head = c(0, 1, 2, 10), thorax = c(0, 1, 10, 2))
#' get_enrichment_mask(rates, lower_q = 0.25, upper_q = 0.75,
#'                     direction = "over")
#' get_enrichment_mask(rates, direction = "under", return_type = "masked")
#' @md
#' @export
get_enrichment_mask <- function(rate_matrix, lower_q = 0.025, upper_q = 0.975,
                                direction = c("both", "over", "under"),
                                return_type = c("binary", "masked"),
                                ties = c("include", "exclude"),
                                constant = c("unclassified", "include")) {
  .check_rate_matrix(rate_matrix)
  .check_regime_quantiles(lower_q, upper_q)
  direction <- match.arg(direction)
  return_type <- match.arg(return_type)
  ties <- match.arg(ties)
  constant <- match.arg(constant)
  mask <- matrix(0L, nrow(rate_matrix), ncol(rate_matrix),
                 dimnames = dimnames(rate_matrix))
  for (j in seq_len(ncol(rate_matrix))) {
    rates <- rate_matrix[, j]
    if (constant == "unclassified" && min(rates) == max(rates)) next
    cutoffs <- stats::quantile(rates, c(lower_q, upper_q), names = FALSE)
    high <- rates > cutoffs[2]
    low <- rates < cutoffs[1]
    if (ties == "include") {
      high <- rates >= cutoffs[2]
      low <- rates <= cutoffs[1]
    }
    selected <- high | low
    if (direction == "over") selected <- high
    if (direction == "under") selected <- low
    mask[, j] <- as.integer(selected)
  }
  if (return_type == "masked") return(rate_matrix * mask)
  mask
}

#' Classify branch rates for each stochastic-map replicate
#'
#' Obtain accelerated and decelerated anatomical-region patterns from a
#' [branch_rate_pipeline()] result without averaging away map uncertainty.
#'
#' @param result A result from [branch_rate_pipeline()].
#' @param groups Optional character vector selecting anatomical groups from
#'   `result$groups`, in the desired order. NULL selects every group. Select
#'   elementary regions, excluding a whole-phenome aggregate of those regions.
#' @inheritParams get_enrichment_mask
#' @return An `ontophylo_regimes` list with `accelerated` and `decelerated`
#'   integer arrays (branches x regions x maps), `thresholds` for each region
#'   and map, and a `branch_summary` giving classification frequencies across
#'   maps. `tree`, `groups`, `memberships`, `shared_characters`, and `settings`
#'   preserve provenance. Frequencies are fractions of conditional maps, not
#'   p-values or probabilities that an underlying Q changed.
#' @details Thresholds use all branches with equal weight, separately within
#'   each region and replicate, as in the Hymenoptera scripts. No pooling across
#'   regions, maps, or time slices occurs. Constant regions remain in the region
#'   denominator. Shared characters among selected regions trigger a warning and
#'   are recorded: overlapping regions are not independent anatomical units.
#'   Map columns are paired by position; if named, their names must also agree.
#' @examples
#' data("hym_stm")
#' maps <- hym_stm[1:2]
#' groups <- list(region_A = names(maps)[1], region_B = names(maps)[2])
#' rates <- branch_rate_pipeline(maps, groups, marginal = FALSE)
#' regimes <- classify_branch_rates(rates)
#' head(regimes$thresholds)
#' head(regimes$branch_summary)
#' @md
#' @export
classify_branch_rates <- function(result, groups = NULL,
                                   lower_q = 0.05, upper_q = 0.95,
                                   ties = c("include", "exclude"),
                                   constant = c("unclassified", "include")) {
  ties <- match.arg(ties)
  constant <- match.arg(constant)
  .check_regime_quantiles(lower_q, upper_q)
  if (!is.list(result) || is.null(result$groups)) {
    stop("Supply a branch_rate_pipeline result.", call. = FALSE)
  }
  .check_named_groups(result$groups, "result$groups")
  if (is.null(groups)) groups <- names(result$groups)
  if (!is.character(groups) || length(groups) == 0L || anyNA(groups) ||
      anyDuplicated(groups) || any(!groups %in% names(result$groups))) {
    stop("groups must select unique, known anatomical groups.", call. = FALSE)
  }
  first <- result$groups[[groups[1]]]
  .check_rate_result(first)
  tree <- first$tree
  n_edges <- nrow(tree$edge)
  n_maps <- ncol(first$rate_samples)
  for (group in groups) {
    current <- result$groups[[group]]
    .check_rate_result(current)
    same_tree <- identical(current$tree$tip.label, tree$tip.label) &&
      identical(dim(current$tree$edge), dim(tree$edge)) &&
      all(current$tree$edge == tree$edge) &&
      isTRUE(all.equal(current$tree$edge.length, tree$edge.length,
                       tolerance = 1e-10, check.attributes = FALSE))
    if (!same_tree || ncol(current$rate_samples) != n_maps ||
        !identical(colnames(current$rate_samples), colnames(first$rate_samples))) {
      stop("Groups must share the same tree, edge order and map columns.", call. = FALSE)
    }
  }
  memberships <- result$memberships[groups]
  shared <- character()
  if (!is.null(result$memberships)) {
    ids <- unlist(memberships, use.names = FALSE)
    shared <- unique(ids[duplicated(ids)])
    if (length(shared) > 0L) {
      warning("Selected groups share characters; phenomic metrics describe overlapping regions. ",
              "See shared_characters, and exclude whole-phenome aggregates.", call. = FALSE)
    }
  }
  map_names <- colnames(first$rate_samples)
  if (is.null(map_names)) map_names <- paste0("map_", seq_len(n_maps))
  dimensions <- list(paste0("branch_", seq_len(n_edges)), groups, map_names)
  accelerated <- array(0L, c(n_edges, length(groups), n_maps), dimnames = dimensions)
  decelerated <- accelerated
  thresholds <- vector("list", n_maps)
  branch_tables <- vector("list", length(groups))
  for (i in seq_len(n_maps)) {
    rates <- matrix(0, n_edges, length(groups), dimnames = dimensions[1:2])
    lower <- upper <- numeric(length(groups))
    is_constant <- logical(length(groups))
    for (j in seq_along(groups)) {
      rates[, j] <- result$groups[[groups[j]]]$rate_samples[, i]
      cutoffs <- stats::quantile(rates[, j], c(lower_q, upper_q), names = FALSE)
      lower[j] <- cutoffs[1]
      upper[j] <- cutoffs[2]
      is_constant[j] <- min(rates[, j]) == max(rates[, j])
    }
    high <- get_enrichment_mask(rates, lower_q, upper_q, "over", "binary", ties, constant)
    low <- get_enrichment_mask(rates, lower_q, upper_q, "under", "binary", ties, constant)
    accelerated[, , i] <- high
    decelerated[, , i] <- low
    thresholds[[i]] <- data.frame(map = i, group = groups, lower = lower,
                                  upper = upper, constant = is_constant,
                                  n_accelerated = colSums(high),
                                  n_decelerated = colSums(low),
                                  n_both = colSums(high * low), row.names = NULL)
  }
  for (j in seq_along(groups)) {
    high <- matrix(accelerated[, j, ], nrow = n_edges, ncol = n_maps)
    low <- matrix(decelerated[, j, ], nrow = n_edges, ncol = n_maps)
    branch_tables[[j]] <- data.frame(group = groups[j], edge = seq_len(n_edges),
      parent = tree$edge[, 1], child = tree$edge[, 2],
      accelerated_frequency = rowMeans(high), decelerated_frequency = rowMeans(low),
      both_frequency = rowMeans(high * low), row.names = NULL)
  }
  out <- list(tree = tree, groups = groups, memberships = memberships,
              shared_characters = shared, accelerated = accelerated,
              decelerated = decelerated, thresholds = do.call(rbind, thresholds),
              branch_summary = do.call(rbind, branch_tables),
              settings = list(lower_q = lower_q, upper_q = upper_q,
                              ties = ties, constant = constant))
  class(out) <- "ontophylo_regimes"
  out
}

#' Calculate effective diversity of anatomical regime patterns
#'
#' The manuscript's Neff calculation, previously implemented by
#' `diversity_metric(type = "Neff_scaled")` in the Hymenoptera scripts.
#'
#' @param enrichment A nonempty 0/1 matrix: lineages by anatomical regions.
#'   Use an accelerated or decelerated mask; keep all-zero rows.
#' @return A one-row data frame containing `Neff`, `fraction` of lineages with
#'   at least one selected region, `entropy` of their whole binary patterns,
#'   `effective_patterns = exp(entropy)`, and counts `n_lineages`, `n_enriched`,
#'   `n_patterns`. With no enriched lineages all diversity outputs are zero
#'   by convention, including `effective_patterns` (the entropy is then undefined
#'   mathematically, but recorded as zero).
#' @details Neff = fraction * exp(entropy), where natural-log Shannon entropy
#'   is calculated over distinct nonzero binary patterns across lineages. It is
#'   not entropy over individual region frequencies. The same formula gives
#'   Neff+ for accelerated masks and Neff- for decelerated masks.
#' @examples
#' patterns <- rbind(c(1, 0), c(0, 1), c(0, 0), c(0, 0))
#' phenomic_diversity(patterns) # fraction = 1/2, exp(H) = 2, Neff = 1
#' @md
#' @export
phenomic_diversity <- function(enrichment) {
  .check_enrichment_matrix(enrichment)
  .phenomic_diversity(enrichment)
}

#' Calculate phenomic coordination Pc
#'
#' The fraction of regions accelerated on each lineage, averaged over all
#' lineages, including those with no accelerated regions. This was previously
#' called enrichment strength (Se); it matches the old scaled
#' `relative_null_metric()` and the calculation in `make-Pc.Rmd`.
#'
#' @inheritParams phenomic_diversity
#' @return A list with `Pc` (one mean), `branch` (one proportion per input row),
#'   `n_lineages`, and `n_regions`.
#' @details Pc lies between zero and one. Pc uses the accelerated mask; applying
#'   the same calculation to a decelerated mask gives the mean fraction of
#'   decelerated regions, not Pc. This is a descriptive extent of co-occurrence,
#'   not a test of anatomical dependence relative to an independence null.
#' @examples
#' patterns <- rbind(c(1, 1), c(0, 0))
#' phenomic_coordination(patterns) # branch proportions 1 and 0; Pc = 1/2
#' @md
#' @export
phenomic_coordination <- function(enrichment) {
  .check_enrichment_matrix(enrichment)
  branch <- rowMeans(enrichment)
  list(Pc = mean(branch), branch = branch, n_lineages = nrow(enrichment),
       n_regions = ncol(enrichment))
}

#' Calculate Neff+, Neff- and Pc through time
#'
#' Integrate regime classification, phenomic diversity and coordination while
#' preserving map replicates until the final summary.
#'
#' @inheritParams classify_branch_rates
#' @inheritParams branch_rate_thru_time
#' @return A list with `regimes` from [classify_branch_rates()], `samples`
#'   (time x map matrices named `Neff_plus`, `Neff_minus`, `Pc`, `fraction_plus`,
#'   `fraction_minus`, `entropy_plus`, `entropy_minus`, `effective_patterns_plus`,
#'   `effective_patterns_minus`, `decelerated_fraction`), a long `summary`
#'   (metric, time, lineage count, mean, median, SD and 95% map quantiles), and
#'   `branch_coordination` (branch x map proportions), its summary and settings.
#' @details Calls [classify_branch_rates()] within each map, then calculates
#'   metrics across the branches intersecting each requested time. Descendants
#'   replace ancestral branches at nodes, as for marginal rates. Returned `time`
#'   is always age before present. All regions and all active branches remain in
#'   the denominator. Conditional map quantiles do not include uncertainty in Q,
#'   tree, taxon sampling or region definition and are not hypothesis tests.
#' @examples
#' data("hym_stm")
#' maps <- hym_stm[1:2]
#' groups <- list(region_A = names(maps)[1], region_B = names(maps)[2])
#' rates <- branch_rate_pipeline(maps, groups, marginal = FALSE)
#' metrics <- phenomic_metrics(rates, n_points = 10)
#' head(metrics$summary)
#' plot_marginal_rates(metrics$summary[metrics$summary$metric == "Pc", ],
#'                     ylab = "Phenomic coordination (Pc)")
#' @md
#' @export
phenomic_metrics <- function(result, groups = NULL, lower_q = 0.05, upper_q = 0.95,
                              ties = c("include", "exclude"),
                              constant = c("unclassified", "include"),
                              n_points = 100, time_points = NULL,
                              time_scale = c("root", "age")) {
  ties <- match.arg(ties)
  constant <- match.arg(constant)
  time_scale <- match.arg(time_scale)
  regimes <- classify_branch_rates(result, groups, lower_q, upper_q, ties, constant)
  slices <- .branch_time_slices(regimes$tree, n_points, time_points, time_scale)
  sizes <- dim(regimes$accelerated)
  n_times <- length(slices$edges)
  # Many requested times intersect the same branches; evaluate each slice once.
  slice_keys <- character(n_times)
  for (t in seq_len(n_times)) slice_keys[t] <- paste(slices$edges[[t]], collapse = ",")
  unique_times <- which(!duplicated(slice_keys))
  source_rows <- unique_times[match(slice_keys, slice_keys[unique_times])]
  n_maps <- sizes[3]
  metric_names <- c("Neff_plus", "Neff_minus", "Pc", "fraction_plus",
                     "fraction_minus", "entropy_plus", "entropy_minus",
                     "effective_patterns_plus", "effective_patterns_minus",
                     "decelerated_fraction")
  samples <- vector("list", length(metric_names))
  names(samples) <- metric_names
  sample_names <- list(as.character(slices$time), dimnames(regimes$accelerated)[[3]])
  for (metric in metric_names) {
    samples[[metric]] <- matrix(0, n_times, n_maps, dimnames = sample_names)
  }
  branch_coordination <- matrix(0, sizes[1], n_maps,
                                dimnames = dimnames(regimes$accelerated)[c(1, 3)])
  for (i in seq_len(n_maps)) {
    high <- matrix(regimes$accelerated[, , i], sizes[1], sizes[2])
    low <- matrix(regimes$decelerated[, , i], sizes[1], sizes[2])
    branch_coordination[, i] <- rowMeans(high)
    for (t in unique_times) {
      ids <- slices$edges[[t]]
      high_now <- high[ids, , drop = FALSE]
      low_now <- low[ids, , drop = FALSE]
      plus <- .phenomic_diversity(high_now)
      minus <- .phenomic_diversity(low_now)
      samples$Neff_plus[t, i] <- plus$Neff
      samples$Neff_minus[t, i] <- minus$Neff
      samples$Pc[t, i] <- mean(branch_coordination[ids, i])
      samples$fraction_plus[t, i] <- plus$fraction
      samples$fraction_minus[t, i] <- minus$fraction
      samples$entropy_plus[t, i] <- plus$entropy
      samples$entropy_minus[t, i] <- minus$entropy
      samples$effective_patterns_plus[t, i] <- plus$effective_patterns
      samples$effective_patterns_minus[t, i] <- minus$effective_patterns
      samples$decelerated_fraction[t, i] <- mean(low_now)
    }
  }
  for (metric in metric_names) {
    samples[[metric]] <- samples[[metric]][source_rows, , drop = FALSE]
    rownames(samples[[metric]]) <- sample_names[[1]]
  }
  summaries <- vector("list", length(metric_names))
  for (j in seq_along(metric_names)) {
    metric <- metric_names[j]
    summaries[[j]] <- data.frame(metric = metric, time = slices$time,
      n_branches = lengths(slices$edges), .branch_summary(samples[[metric]]),
      row.names = NULL)
  }
  branch_summary <- data.frame(edge = seq_len(sizes[1]),
    parent = regimes$tree$edge[, 1], child = regimes$tree$edge[, 2],
    .branch_summary(branch_coordination), row.names = NULL)
  list(regimes = regimes, samples = samples, summary = do.call(rbind, summaries),
       branch_coordination = branch_coordination,
       branch_coordination_summary = branch_summary,
       settings = c(regimes$settings, list(time_scale = time_scale,
                                           time_points = time_points, n_points = n_points)))
}

.check_regime_quantiles <- function(lower_q, upper_q) {
  if (!is.numeric(lower_q) || !is.numeric(upper_q) || length(lower_q) != 1L ||
      length(upper_q) != 1L || !is.finite(lower_q) || !is.finite(upper_q) ||
      lower_q < 0 || upper_q > 1 || lower_q >= upper_q) {
    stop("Require 0 <= lower_q < upper_q <= 1.", call. = FALSE)
  }
}

.check_enrichment_matrix <- function(enrichment) {
  if (!is.matrix(enrichment) || !(is.numeric(enrichment) || is.logical(enrichment)) ||
      any(dim(enrichment) == 0L) || anyNA(enrichment) ||
      any(!enrichment %in% c(0, 1))) {
    stop("enrichment must be a nonempty 0/1 matrix without missing values.", call. = FALSE)
  }
}

.phenomic_diversity <- function(enrichment) {
  keep <- rowSums(enrichment) > 0
  n_enriched <- sum(keep)
  fraction <- n_enriched / nrow(enrichment)
  entropy <- effective <- Neff <- n_patterns <- 0
  if (n_enriched > 0L) {
    selected <- enrichment[keep, , drop = FALSE]
    patterns <- character(n_enriched)
    for (i in seq_len(n_enriched)) {
      patterns[i] <- paste0(as.integer(selected[i, ]), collapse = "")
    }
    counts <- table(patterns)
    frequencies <- as.numeric(counts) / n_enriched
    n_patterns <- length(counts)
    entropy <- -sum(frequencies * log(frequencies))
    effective <- exp(entropy)
    Neff <- fraction * effective
  }
  data.frame(Neff = Neff, fraction = fraction, entropy = entropy,
              effective_patterns = effective, n_lineages = nrow(enrichment),
              n_enriched = n_enriched, n_patterns = n_patterns)
}
