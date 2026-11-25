#' Run the branch-rate pipeline for one or more character groups
#'
#' Combine original stochastic maps into branch-rate samples, summarize their
#' uncertainty, and calculate marginal rates for each selected character group.
#'
#' @param smaps A named list of characters, each containing equally many
#'   original stochastic maps, or an `ontophylo_stm` result from [paramo_stm()].
#'   Mapping results also support disk-backed and branch-count-only storage.
#' @param groups Optional named list of unique character IDs for each anatomical
#'   group, for example the output of [RAC_query()]. NULL combines all characters
#'   into a group named `combined`. Different groups may overlap; they are
#'   analyzed separately and are never added together.
#' @param normalize `"none"` for total character changes per time unit, or
#'   `"characters"` to divide rates by the number of characters in each group.
#' @param contrasts Optional named list of predefined branch comparisons. Each
#'   element is a list with `focal` and `reference` edge indices, passed to
#'   [branch_rate_contrast()]. Comparisons use branch-duration weighting.
#' @param phenomic NULL (default) skips phenomic metrics. Otherwise supply a
#'   named list of options for [phenomic_metrics()], such as
#'   `list(groups = c("region_A", "region_B"), upper_q = 0.95)`. An empty list
#'   uses all groups. Options are `groups`, `lower_q`, `upper_q`, `ties`, and
#'   `constant`. The pipeline time grid is reused. Select elementary regions,
#'   excluding a whole-phenome aggregate, to avoid counting it as another region.
#' @inheritParams branch_rates
#' @return A list with `groups` (one [branch_rates()] result per character group),
#'   combined `branch_summary` and `marginal_summary` tables, a `group_summary`
#'   recording membership and normalization, and `contrasts` containing their
#'   samples and summaries. Every combined summary row includes a `group` ID.
#'   Raw rate samples remain available in `groups[[name]]$rate_samples`.
#'   `input_status` retains the original status table for mapping-result input.
#'   When requested, `phenomic` contains Neff+, Neff-, Pc and regime results.
#' @details Group memberships are checked before analysis. All supplied maps
#'   must describe the same tree in the same edge order, with equally many
#'   replicates. Failed/skipped mapping results require explicit exclusion via
#'   `groups`; their status remains in `input_status`. Maps are paired by replicate
#'   position across characters; retain joint draw identities when necessary.
#'   Normalization scales all samples and summaries consistently. Conditional
#'   mapping intervals are not hypothesis tests of branch-specific Q matrices.
#' @examples
#' data("hym_stm")
#' maps <- hym_stm[1:2]
#' groups <- list(combined = names(maps), first = names(maps)[1])
#' result <- branch_rate_pipeline(maps, groups, normalize = "characters",
#'                                 n_points = 10)
#' head(result$branch_summary)
#' head(result$marginal_summary)
#' @md
#' @export
branch_rate_pipeline <- function(smaps, groups = NULL,
                                  normalize = c("none", "characters"),
                                  contrasts = NULL, n_points = 100,
                                  time_points = NULL,
                                  time_scale = c("root", "age"), marginal = TRUE,
                                  phenomic = NULL) {
  normalize <- match.arg(normalize)
  time_scale <- match.arg(time_scale)
  if (!is.null(phenomic)) {
    allowed <- c("groups", "lower_q", "upper_q", "ties", "constant")
    if (!is.list(phenomic) || (length(phenomic) > 0L &&
        (is.null(names(phenomic)) || anyNA(names(phenomic)) ||
         anyDuplicated(names(phenomic)) || any(!names(phenomic) %in% allowed)))) {
      stop("phenomic must be a named options list for phenomic_metrics().", call. = FALSE)
    }
  }
  mapping_result <- inherits(smaps, "ontophylo_stm")
  input_status <- NULL
  if (mapping_result) {
    character_names <- names(smaps$characters)
    .check_character_list(smaps$characters)
    reference_tree <- smaps$tree
    .check_branch_tree(reference_tree)
    n_maps <- smaps$settings$n_stm
    input_status <- smaps$status
  } else {
    .check_character_list(smaps)
    character_names <- names(smaps)
  }
  if (is.null(character_names)) {
    stop("smaps must have character names.", call. = FALSE)
  }
  if (is.null(groups)) groups <- list(combined = character_names)
  .check_named_groups(groups, "groups")
  for (i in seq_along(groups)) {
    ids <- groups[[i]]
    if (!is.character(ids) || length(ids) == 0L || anyNA(ids) ||
        anyDuplicated(ids) || any(!ids %in% character_names)) {
      stop("Each group must contain unique, known character names.", call. = FALSE)
    }
  }
  # Validate the full input, including characters not used by a group.
  if (!mapping_result) {
    reference_maps <- .check_branch_maps(smaps[[1]])
    reference_tree <- reference_maps[[1]]
    n_maps <- length(reference_maps)
    for (i in seq_along(smaps)) {
      maps <- .check_branch_maps(smaps[[i]], reference = reference_tree)
      if (length(maps) != n_maps) {
        stop("All characters must have the same number of map replicates.", call. = FALSE)
      }
    }
  }
  if (!is.null(contrasts)) {
    .check_named_groups(contrasts, "contrasts")
    for (i in seq_along(contrasts)) {
      specification <- contrasts[[i]]
      if (!is.list(specification) || length(specification) != 2L ||
          !setequal(names(specification), c("focal", "reference"))) {
        stop("Each contrast must have focal and reference edge indices.", call. = FALSE)
      }
      .check_contrast_edges(specification$focal, specification$reference,
                            nrow(reference_tree$edge))
    }
  }

  results <- vector("list", length(groups))
  names(results) <- names(groups)
  branch_tables <- vector("list", length(groups))
  marginal_tables <- vector("list", length(groups))
  comparison_results <- vector("list", length(groups))
  names(comparison_results) <- names(groups)
  group_summary <- data.frame(group = names(groups), n_characters = lengths(groups),
                              n_maps = n_maps, divisor = 1,
                              normalization = normalize, row.names = NULL)

  for (i in seq_along(groups)) {
    ids <- groups[[i]]
    selected <- smaps[ids]
    if (mapping_result) {
      selected <- smaps
      selected$characters <- smaps$characters[ids]
      selected$status <- smaps$status[match(ids, smaps$status$character), , drop = FALSE]
    }
    result <- branch_rates(selected, n_points, time_points, time_scale, marginal)
    divisor <- 1
    if (normalize == "characters") divisor <- length(ids)
    group_summary$divisor[i] <- divisor
    result$rate_samples <- result$rate_samples / divisor
    result$branch_summary[, c("mean", "median", "sd", "lower", "upper")] <-
      result$branch_summary[, c("mean", "median", "sd", "lower", "upper")] / divisor
    if (!is.null(result$marginal_samples)) {
      result$marginal_samples <- result$marginal_samples / divisor
      result$marginal_summary[, c("mean", "median", "sd", "lower", "upper")] <-
        result$marginal_summary[, c("mean", "median", "sd", "lower", "upper")] / divisor
      marginal_tables[[i]] <- data.frame(group = names(groups)[i],
                                          result$marginal_summary, row.names = NULL)
    }
    result$normalization <- normalize
    result$divisor <- divisor
    results[[i]] <- result
    branch_tables[[i]] <- data.frame(group = names(groups)[i], result$branch_summary,
                                      row.names = NULL)
    if (!is.null(contrasts)) {
      group_contrasts <- vector("list", length(contrasts))
      names(group_contrasts) <- names(contrasts)
      for (j in seq_along(contrasts)) {
        specification <- contrasts[[j]]
        group_contrasts[[j]] <- branch_rate_contrast(
          result, specification$focal, specification$reference)
      }
      comparison_results[[i]] <- group_contrasts
    }
  }
  marginal_summary <- NULL
  if (marginal) marginal_summary <- do.call(rbind, marginal_tables)
  out <- list(groups = results, group_summary = group_summary, memberships = groups,
       branch_summary = do.call(rbind, branch_tables),
       marginal_summary = marginal_summary,
       contrasts = if (is.null(contrasts)) NULL else comparison_results,
       input_status = input_status)
  if (!is.null(phenomic)) {
    arguments <- c(list(result = out, n_points = n_points,
                         time_points = time_points, time_scale = time_scale), phenomic)
    out$phenomic <- do.call(phenomic_metrics, arguments)
  }
  out
}

#' Compare mapped rates between two predefined branch sets
#'
#' Adapt the lineage-time-weighted clade comparison used in the Hymenoptera
#' simulation analysis. Average rates within each set for every map replicate,
#' then summarize their difference and ratio across replicates.
#'
#' @param result A result from [branch_rates()], or one element of the `groups`
#'   returned by [branch_rate_pipeline()].
#' @param focal,reference Nonempty, disjoint vectors of edge-row indices.
#' @param weighting `"duration"` weights by branch length, so the average is
#'   total mapped changes divided by total lineage time. `"equal"` gives each
#'   branch equal weight.
#' @return A list with `samples` (focal and reference rates, difference and ratio
#'   for each map), `summary` (means, medians, SDs, 95% map quantiles), branch
#'   indices, total lineage times, weighting, and `fraction_focal_greater`.
#' @details These comparisons describe reconstructed histories conditional on
#'   the mapping model. `fraction_focal_greater` is not a p-value. Mapping
#'   intervals are not calibrated frequentist confidence intervals for a true
#'   generating rate ratio. Define branch sets before looking at inferred rates.
#'
#'   When the reference rate is zero, the corresponding ratio is undefined and
#'   stored as NA. If any ratios are undefined, all ratio summary statistics
#'   are NA; `n_defined` and `n_undefined` make this explicit. Differences remain
#'   defined. There is no pseudocount or silent removal of zero denominators.
#' @examples
#' data("hym_stm")
#' result <- branch_rates(hym_stm[1:2], marginal = FALSE)
#' comparison <- branch_rate_contrast(result, focal = 1:5, reference = 6:10)
#' comparison$summary
#' @md
#' @export
branch_rate_contrast <- function(result, focal, reference,
                                  weighting = c("duration", "equal")) {
  weighting <- match.arg(weighting)
  .check_rate_result(result)
  .check_contrast_edges(focal, reference, nrow(result$tree$edge))
  focal_weights <- rep(1, length(focal))
  reference_weights <- rep(1, length(reference))
  if (weighting == "duration") {
    focal_weights <- result$tree$edge.length[focal]
    reference_weights <- result$tree$edge.length[reference]
  }
  focal_weights <- focal_weights / sum(focal_weights)
  reference_weights <- reference_weights / sum(reference_weights)
  focal_rates <- as.numeric(crossprod(focal_weights,
                                     result$rate_samples[focal, , drop = FALSE]))
  reference_rates <- as.numeric(crossprod(reference_weights,
                                         result$rate_samples[reference, , drop = FALSE]))
  ratio <- rep(NA_real_, length(focal_rates))
  defined <- reference_rates > 0
  ratio[defined] <- focal_rates[defined] / reference_rates[defined]
  samples <- data.frame(map = seq_along(focal_rates), focal = focal_rates,
                        reference = reference_rates,
                        difference = focal_rates - reference_rates, ratio = ratio)
  values <- t(as.matrix(samples[, c("focal", "reference", "difference", "ratio")]))
  summary <- data.frame(metric = rownames(values), n_defined = rowSums(is.finite(values)),
                        n_undefined = rowSums(!is.finite(values)), mean = NA_real_,
                        median = NA_real_, sd = NA_real_, lower = NA_real_, upper = NA_real_,
                        row.names = NULL)
  for (i in seq_len(nrow(values))) {
    if (all(is.finite(values[i, ]))) {
      statistics <- .branch_summary(values[i, , drop = FALSE])
      summary[i, c("mean", "median", "sd", "lower", "upper")] <- statistics
    }
  }
  list(samples = samples, summary = summary, focal = focal, reference = reference,
       focal_lineage_time = sum(result$tree$edge.length[focal]),
       reference_lineage_time = sum(result$tree$edge.length[reference]),
       weighting = weighting,
       fraction_focal_greater = mean(focal_rates > reference_rates))
}

#' Plot a phylogeny coloured by branch values
#'
#' Plot branch-specific rates or other finite values in tree edge order.
#'
#' @param tree A `phylo` tree.
#' @param branch_values One finite numeric value per edge, in tree edge order.
#' @param palette A vector of colours from low to high.
#' @param title Optional plot title. NULL leaves the plot untitled.
#' @param legend_title Title for the colour legend.
#' @param legend Show the legend.
#' @param ... Additional arguments passed to [ape::plot.phylo()].
#' @return Invisibly, a data frame containing edge IDs, values and colours.
#' @details Adapted from the Hymenoptera plotting helper. Constant values receive
#'   the middle palette colour. Plotting does not alter the input tree or rates.
#'   For the default right-facing rectangular tree, extra horizontal space is
#'   reserved for the legend unless `x.lim` is supplied explicitly.
#' @examples
#' tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
#' plot_branch_colored_tree(tree, c(0.1, 0.5, 0.4, 0.1), edge.width = 3)
#' @author Sergei Tarasov
#' @md
#' @export
plot_branch_colored_tree <- function(tree, branch_values,
                                      palette = grDevices::colorRampPalette(
                                        c("#440154", "#21908C", "#FDE725"))(100),
                                      title = NULL, legend_title = "Rate", legend = TRUE,
                                      ...) {
  .check_branch_tree(tree)
  if (!is.numeric(branch_values) || length(branch_values) != nrow(tree$edge) ||
      any(!is.finite(branch_values))) {
    stop("Supply one finite value per tree edge.", call. = FALSE)
  }
  if (!is.character(palette) || length(palette) < 2L || anyNA(palette)) {
    stop("palette must contain at least two colours.", call. = FALSE)
  }
  grDevices::col2rgb(palette)
  if (!is.logical(legend) || length(legend) != 1L || is.na(legend)) {
    stop("legend must be TRUE or FALSE.", call. = FALSE)
  }
  limits <- range(branch_values)
  scaled <- rep(0.5, length(branch_values))
  if (diff(limits) > 0) scaled <- (branch_values - limits[1]) / diff(limits)
  indices <- 1L + floor(scaled * (length(palette) - 1L))
  colours <- palette[indices]
  plot_arguments <- list(...)
  rectangular <- is.null(plot_arguments$type) ||
    plot_arguments$type %in% c("phylogram", "cladogram")
  rightwards <- is.null(plot_arguments$direction) ||
    plot_arguments$direction == "rightwards"
  if (legend && rectangular && rightwards && is.null(plot_arguments$x.lim)) {
    height <- max(ape::node.depth.edgelength(tree))
    plot_arguments$x.lim <- c(0, 1.5 * height)
  }
  do.call(ape::plot.phylo, c(list(x = tree, edge.color = colours), plot_arguments))
  if (!is.null(title)) graphics::title(main = title)
  if (legend) {
    legend_values <- unique(c(limits[1], mean(limits), limits[2]))
    legend_indices <- round(seq(1, length(palette), length.out = 3))
    if (length(legend_values) == 1L) legend_indices <- indices[1]
    graphics::legend("topright", legend = format(signif(legend_values, 3)),
                      fill = palette[legend_indices], title = legend_title,
                      border = NA, bty = "n", cex = 0.8)
  }
  invisible(data.frame(edge = seq_along(branch_values), value = branch_values,
                        colour = colours))
}

#' Plot a marginal-rate curve with its mapping interval
#'
#' @param summary A marginal summary from [branch_rates()] or
#'   [branch_rate_thru_time_Multi()]. Supply one character group at a time.
#' @param xlab,ylab Axis labels. Time is age before the present in tree units.
#' @param colour Line colour.
#' @param fill Interval colour.
#' @return A ggplot object; the band is the pointwise 95% map interval. Rates
#'   and bounds are drawn as steps on the supplied evaluation grid.
#' @examples
#' tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
#' rates <- list(c(1, 2, 2, 1), c(1, 3, 1, 1))
#' summary <- branch_rate_thru_time_Multi(tree, rates, n_points = 20)
#' plot_marginal_rates(summary)
#' @md
#' @export
plot_marginal_rates <- function(summary, xlab = "Age (tree time units)",
                                 ylab = "Mapped changes per time unit",
                                 colour = "#0072B2", fill = "#B3D9EB") {
  required <- c("time", "mean", "lower", "upper")
  if (!is.data.frame(summary) || !all(required %in% names(summary)) ||
      nrow(summary) == 0L || anyDuplicated(summary$time)) {
    stop("Supply a marginal summary for one group with unique times.", call. = FALSE)
  }
  for (name in required) {
    if (!is.numeric(summary[[name]]) || any(!is.finite(summary[[name]]))) {
      stop("Marginal times, means and bounds must be finite.", call. = FALSE)
    }
  }
  if (any(summary$lower > summary$upper)) {
    stop("Lower bounds must not exceed upper bounds.", call. = FALSE)
  }
  # Declare the column names used in ggplot's data mask.
  time <- mean <- lower <- upper <- NULL
  ordered <- summary[order(summary$time, decreasing = TRUE), , drop = FALSE]
  band <- ordered
  if (nrow(ordered) > 1L) {
    rows <- c(rep(seq_len(nrow(ordered) - 1L), each = 2), nrow(ordered))
    band <- ordered[rows, , drop = FALSE]
    band$time <- c(ordered$time[1], rep(ordered$time[-1], each = 2))
  }
  ggplot2::ggplot(summary, ggplot2::aes(x = time, y = mean)) +
    ggplot2::geom_ribbon(data = band, ggplot2::aes(ymin = lower, ymax = upper),
                         fill = fill, alpha = 0.6) +
    ggplot2::geom_step(colour = colour, direction = "hv") +
    ggplot2::scale_x_reverse() +
    ggplot2::labs(x = xlab, y = ylab) + ggplot2::theme_classic()
}

.check_named_groups <- function(groups, name) {
  ids <- names(groups)
  if (!is.list(groups) || length(groups) == 0L || is.null(ids) ||
      anyNA(ids) || any(ids == "") || anyDuplicated(ids)) {
    stop(name, " must be a nonempty list with unique names.", call. = FALSE)
  }
}

.check_contrast_edges <- function(focal, reference, n_edges) {
  selections <- list(focal, reference)
  for (ids in selections) {
    if (!is.numeric(ids) || length(ids) == 0L || any(!is.finite(ids)) ||
        any(ids != floor(ids)) || any(ids < 1 | ids > n_edges) || anyDuplicated(ids)) {
      stop("Branch sets must contain unique, valid edge-row indices.", call. = FALSE)
    }
  }
  if (length(intersect(focal, reference)) > 0L) {
    stop("Focal and reference branch sets must be disjoint.", call. = FALSE)
  }
}

.check_rate_result <- function(result) {
  if (!is.list(result) || is.null(result$tree) || is.null(result$rate_samples)) {
    stop("Supply a branch_rates result.", call. = FALSE)
  }
  .check_branch_tree(result$tree)
  .check_rate_matrix(result$rate_samples)
  if (nrow(result$rate_samples) != nrow(result$tree$edge)) {
    stop("Rate rows must match the tree edge order.", call. = FALSE)
  }
}
