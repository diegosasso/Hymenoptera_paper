# Project existing ontophylo rate samples onto a shared set of taxa.
# No models are refitted. The returned induced-tree rates conserve expected
# transition counts along contracted paths. Temporal outputs distinguish
# original segment rates from rates made constant across contracted branches.
project_ontophylo <- function(full, keep, bin_width_Myr = 10, tolerance = 1e-6) {
  tree <- full$tree
  if (!inherits(tree, "phylo")) stop("full$tree must be a phylo object.")
  if (!is.character(keep) || anyDuplicated(keep) || length(keep) < 2L) {
    stop("keep must contain at least two unique analysis labels.")
  }
  if (!all(keep %in% tree$tip.label)) stop("Some keep labels are absent from full$tree.")
  if (!is.numeric(bin_width_Myr) || length(bin_width_Myr) != 1L || bin_width_Myr <= 0) {
    stop("bin_width_Myr must be a positive number.")
  }
  n_edges <- nrow(tree$edge)
  n_tips <- length(tree$tip.label)
  if (!is.list(full$samples) || length(full$samples) == 0L) stop("full$samples must be a named list of matrices.")
  for (partition in names(full$samples)) {
    rates <- full$samples[[partition]]
    if (!is.matrix(rates) || nrow(rates) != n_edges || any(!is.finite(rates)) || any(rates < 0)) {
      stop("Each sample matrix must have finite nonnegative values in full-tree edge order.")
    }
  }

  # Build original descendant sets independently of edge order.
  full_descendants <- vector("list", n_edges)
  full_keys <- character(n_edges)
  for (edge_index in seq_len(n_edges)) {
    frontier <- tree$edge[edge_index, 2]
    descendants <- frontier
    while (length(frontier) > 0L) {
      frontier <- tree$edge[tree$edge[, 1] %in% frontier, 2]
      descendants <- c(descendants, frontier)
    }
    labels <- tree$tip.label[descendants[descendants <= n_tips]]
    full_descendants[[edge_index]] <- sort(labels)
    retained_labels <- sort(intersect(labels, keep))
    full_keys[edge_index] <- paste(retained_labels, collapse = "|")
  }

  # Keep the input phylogeny's branch lengths. ape removes the stem above the
  # retained MRCA; we also identify it explicitly by the all-retained key.
  pruned_tree <- ape::keep.tip(tree, keep)
  n_pruned_edges <- nrow(pruned_tree$edge)
  n_pruned_tips <- length(pruned_tree$tip.label)
  mrca_node <- ape::getMRCA(tree, keep)
  full_depth <- ape::node.depth.edgelength(tree)
  full_height <- max(full_depth[seq_len(n_tips)])
  retained_root_age <- full_height - full_depth[mrca_node]
  original_older <- full_height - full_depth[tree$edge[, 1]]
  original_younger <- full_height - full_depth[tree$edge[, 2]]
  pruned_keys <- character(n_pruned_edges)
  pruned_descendants <- vector("list", n_pruned_edges)
  for (edge_index in seq_len(n_pruned_edges)) {
    frontier <- pruned_tree$edge[edge_index, 2]
    descendants <- frontier
    while (length(frontier) > 0L) {
      frontier <- pruned_tree$edge[pruned_tree$edge[, 1] %in% frontier, 2]
      descendants <- c(descendants, frontier)
    }
    labels <- sort(pruned_tree$tip.label[descendants[descendants <= n_pruned_tips]])
    pruned_descendants[[edge_index]] <- labels
    pruned_keys[edge_index] <- paste(labels, collapse = "|")
  }
  if (anyDuplicated(pruned_keys)) stop("Induced tree has ambiguous clade keys.")
  original_to_pruned <- match(full_keys, pruned_keys)
  retained_edges <- which(!is.na(original_to_pruned))
  all_retained_key <- paste(sort(keep), collapse = "|")
  excluded_stem_edges <- which(full_keys == all_retained_key)
  edge_paths <- vector("list", n_pruned_edges)
  names(edge_paths) <- paste0("PRUNED_BR_", seq_len(n_pruned_edges))
  branches <- data.frame(
    edge_id = seq_len(n_pruned_edges), branch_id = names(edge_paths),
    parent = pruned_tree$edge[, 1], child = pruned_tree$edge[, 2],
    length_Myr = pruned_tree$edge.length,
    older_age_Ma = numeric(n_pruned_edges), younger_age_Ma = numeric(n_pruned_edges),
    descendant_tips = pruned_keys, n_original_edges = integer(n_pruned_edges),
    original_edge_ids = character(n_pruned_edges)
  )
  for (edge_index in seq_len(n_pruned_edges)) {
    path <- which(original_to_pruned == edge_index)
    path <- path[order(original_older[path], decreasing = TRUE)]
    if (length(path) == 0L) stop("Pruned edge has no matching full-tree path.")
    if (abs(sum(tree$edge.length[path]) - pruned_tree$edge.length[edge_index]) > tolerance) {
      stop("Contracted path length differs from pruned edge length.")
    }
    edge_paths[[edge_index]] <- path
    branches$older_age_Ma[edge_index] <- max(original_older[path])
    branches$younger_age_Ma[edge_index] <- min(original_younger[path])
    branches$n_original_edges[edge_index] <- length(path)
    branches$original_edge_ids[edge_index] <- paste(path, collapse = "|")
  }

  samples <- vector("list", length(full$samples))
  counts <- vector("list", length(full$samples))
  names(samples) <- names(counts) <- names(full$samples)
  count_errors <- numeric(length(full$samples))
  names(count_errors) <- names(full$samples)
  for (partition in names(full$samples)) {
    original_counts <- full$samples[[partition]] * tree$edge.length
    projected_counts <- matrix(0, n_pruned_edges, ncol(original_counts),
      dimnames = list(branches$branch_id, colnames(original_counts)))
    for (edge_index in seq_len(n_pruned_edges)) {
      projected_counts[edge_index, ] <- colSums(original_counts[edge_paths[[edge_index]], , drop = FALSE])
    }
    projected_rates <- projected_counts / pruned_tree$edge.length
    counts[[partition]] <- projected_counts
    samples[[partition]] <- projected_rates
    count_errors[partition] <- max(abs(colSums(projected_rates * pruned_tree$edge.length) -
      colSums(original_counts[retained_edges, , drop = FALSE])))
    if (count_errors[partition] > tolerance) stop("Projected transition counts were not conserved.")
  }

  # The grid retains the original 77-tip root age. Bins older than the
  # retained MRCA have zero exposure and NA rates, not zero rates.
  breaks <- seq(0, floor(full_height / bin_width_Myr) * bin_width_Myr, by = bin_width_Myr)
  if (tail(breaks, 1) < full_height) breaks <- c(breaks, full_height)
  n_bins <- length(breaks) - 1L
  segment_exposure <- matrix(0, n_bins, n_edges)
  contracted_exposure <- matrix(0, n_bins, n_pruned_edges)
  for (bin_index in seq_len(n_bins)) {
    segment_exposure[bin_index, retained_edges] <- pmax(0,
      pmin(original_older[retained_edges], breaks[bin_index + 1L]) -
      pmax(original_younger[retained_edges], breaks[bin_index]))
    contracted_exposure[bin_index, ] <- pmax(0,
      pmin(branches$older_age_Ma, breaks[bin_index + 1L]) -
      pmax(branches$younger_age_Ma, breaks[bin_index]))
  }
  exposure_total <- rowSums(segment_exposure)
  exposure_difference <- max(abs(exposure_total - rowSums(contracted_exposure)))
  if (exposure_difference > tolerance) stop("Original and contracted marginal exposure differs.")
  segment_weights <- segment_exposure
  contracted_weights <- contracted_exposure
  occupied <- exposure_total > 0
  segment_weights[occupied, ] <- segment_exposure[occupied, , drop = FALSE] / exposure_total[occupied]
  contracted_weights[occupied, ] <- contracted_exposure[occupied, , drop = FALSE] / exposure_total[occupied]
  segment_curves <- vector("list", length(samples))
  contracted_curves <- vector("list", length(samples))
  names(segment_curves) <- names(contracted_curves) <- names(samples)
  marginal_count_errors <- numeric(length(samples))
  names(marginal_count_errors) <- names(samples)
  for (partition in names(samples)) {
    segment_curve <- segment_weights %*% full$samples[[partition]]
    contracted_curve <- contracted_weights %*% samples[[partition]]
    marginal_count_errors[partition] <- max(abs(
      colSums(segment_curve * exposure_total) - colSums(counts[[partition]])))
    if (marginal_count_errors[partition] > tolerance) stop("Time bins did not conserve projected counts.")
    segment_curve[!occupied, ] <- NA_real_
    contracted_curve[!occupied, ] <- NA_real_
    segment_curves[[partition]] <- segment_curve
    contracted_curves[[partition]] <- contracted_curve
  }
  mapping <- data.frame(original_edge_id = seq_len(n_edges),
    pruned_edge_id = original_to_pruned,
    disposition = ifelse(is.na(original_to_pruned), "dropped_lineage", "retained_path"))
  mapping$disposition[excluded_stem_edges] <- "stem_above_retained_MRCA"
  checks <- list(
    maximum_path_length_error = max(abs(branches$older_age_Ma - branches$younger_age_Ma - branches$length_Myr)),
    total_retained_lineage_Myr = sum(tree$edge.length[retained_edges]),
    retained_minus_pruned_length = sum(tree$edge.length[retained_edges]) - sum(pruned_tree$edge.length),
    bin_exposure_minus_retained_length = sum(exposure_total) - sum(tree$edge.length[retained_edges]),
    maximum_binned_exposure_difference = exposure_difference,
    maximum_count_errors = count_errors,
    maximum_marginal_count_errors = marginal_count_errors
  )
  list(tree = pruned_tree, branches = branches, samples = samples, counts = counts,
    character_ids = full$character_ids, units = full$units, uncertainty = full$uncertainty,
    retained_labels = pruned_tree$tip.label, original_root_age_Ma = full_height,
    retained_root_age_Ma = retained_root_age, retained_root_original_node = mrca_node,
    excluded_stem_edges = excluded_stem_edges, original_to_pruned = mapping,
    edge_paths = edge_paths, retained_original_edges = retained_edges,
    marginal = list(breaks_Ma = breaks, exposure_lineage_Myr = exposure_total,
      original_segment_samples = segment_curves, contracted_branch_samples = contracted_curves,
      original_segment_exposure = segment_exposure, contracted_branch_exposure = contracted_exposure),
    checks = checks)
}
