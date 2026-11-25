# Run from comparisons/haika-wings after retrieving Roihu fit results.
# Rscript --vanilla scripts/05-compare-rates.R [--final]
# No fitting is performed here. Missing/failed fits are explicitly recorded.
library(ape)
source("scripts/project-ontophylo.R")

arguments <- commandArgs(trailingOnly = TRUE)
if (length(setdiff(arguments, "--final")) > 0L) stop("Only --final is supported.")
final_required <- "--final" %in% arguments
dir.create("results", showWarnings = FALSE)
datasets <- readRDS("data/prepared-datasets.rds")
manifest <- read.csv("data/fit-manifest.csv", stringsAsFactors = FALSE)
full_ontophylo <- readRDS("data/ontophylo-wing-branch-rate-samples.rds")
tip_metadata <- read.csv("data/tree-tip-metadata.csv", stringsAsFactors = FALSE)
partitions <- c("fore_wing", "fore_wing_sclerites", "all_wings", "hind_wing")
stopifnot(!anyDuplicated(manifest$task_id), all(partitions %in% names(full_ontophylo$samples)))
input_md5 <- tools::md5sum(c("data/prepared-datasets.rds", "data/fit-manifest.csv"))

# Small numerical helpers are shared by branch and temporal summaries.
sample_summary <- function(samples) {
  answer <- matrix(NA_real_, nrow(samples), 4L,
    dimnames = list(NULL, c("mean", "lower", "median", "upper")))
  for (row_index in seq_len(nrow(samples))) {
    values <- samples[row_index, is.finite(samples[row_index, ])]
    if (length(values) > 0L) {
      answer[row_index, ] <- c(mean(values), quantile(values, c(0.025, 0.5, 0.975)))
    }
  }
  as.data.frame(answer)
}

rank_draws <- function(samples) {
  ranks <- matrix(NA_real_, nrow(samples), ncol(samples))
  if (nrow(samples) < 3L) return(ranks)
  for (draw_index in seq_len(ncol(samples))) {
    values <- samples[, draw_index]
    if (all(is.finite(values))) {
      centered <- rank(values, ties.method = "average")
      centered <- centered - mean(centered)
      magnitude <- sqrt(sum(centered^2))
      if (magnitude > 0) ranks[, draw_index] <- centered / magnitude
    }
  }
  ranks
}

rank_association <- function(haika, ontophylo_mean, standardized_draw_ranks) {
  rho <- NA_real_
  draws <- rep(NA_real_, ncol(standardized_draw_ranks))
  valid <- is.finite(haika) & is.finite(ontophylo_mean)
  if (length(haika) >= 3L && all(valid) && length(unique(haika)) > 1L) {
    haika_rank <- rank(haika, ties.method = "average")
    centered <- haika_rank - mean(haika_rank)
    centered <- centered / sqrt(sum(centered^2))
    draws <- as.numeric(crossprod(centered, standardized_draw_ranks))
    if (length(unique(ontophylo_mean)) > 1L) {
      rho <- cor(haika, ontophylo_mean, method = "spearman")
    }
  }
  finite_draws <- draws[is.finite(draws)]
  interval <- rep(NA_real_, 3L)
  if (length(finite_draws) > 0L) interval <- quantile(finite_draws, c(0.025, 0.5, 0.975))
  list(rho = rho, draws = draws, lower = interval[1], median = interval[2],
    upper = interval[3], n_valid_draws = length(finite_draws))
}

# Ontophylo projections are cached by actual retained taxon set, so datasets
# that differ only in neural-network seed share exactly the same projection.
projections <- list()
projection_ids <- character(length(datasets))
names(projection_ids) <- names(datasets)
projection_tip_keys <- character()
for (dataset_id in names(datasets)) {
  dataset <- datasets[[dataset_id]]
  taxon_key <- paste(sort(dataset$tree$tip.label), collapse = "|")
  projection_index <- match(taxon_key, projection_tip_keys)
  if (is.na(projection_index)) {
    projection_index <- length(projections) + 1L
    projection_tip_keys[projection_index] <- taxon_key
    projection_id <- paste0("projection_", projection_index)
    projection <- project_ontophylo(full_ontophylo, dataset$tree$tip.label)
    projections[[projection_id]] <- projection
  }
  projection_id <- paste0("projection_", projection_index)
  projection_ids[dataset_id] <- projection_id
  projection <- projections[[projection_id]]
  if (!identical(dataset$tree$tip.label, projection$tree$tip.label) ||
      !identical(dataset$tree$edge, projection$tree$edge) ||
      max(abs(dataset$tree$edge.length - projection$tree$edge.length)) > 1e-10) {
    stop("Prepared dataset edge order differs from its checked projection: ", dataset_id)
  }
}

# Validate every expected fit before using its values.
fit_status <- manifest
fit_status$status <- "missing"
fit_status$success <- FALSE
fit_status$message <- "Expected result file is absent."
fit_status$file <- sprintf("results/fits/task-%04d.rds", manifest$task_id)
fit_status$convergence <- NA_integer_
fit_status$elapsed_seconds <- NA_real_
fit_status$dataset_lambda_complete <- FALSE
fits <- vector("list", nrow(manifest))
for (manifest_index in seq_len(nrow(manifest))) {
  file <- fit_status$file[manifest_index]
  if (!file.exists(file)) next
  fit <- tryCatch(readRDS(file), error = function(error) error)
  if (inherits(fit, "error")) {
    fit_status$status[manifest_index] <- "unreadable"
    fit_status$message[manifest_index] <- conditionMessage(fit)
    next
  }
  validation_message <- ""
  if (!isTRUE(fit$success)) {
    fit_status$status[manifest_index] <- "failed"
    fit_status$message[manifest_index] <- if (is.null(fit$error)) "Fit success flag is false." else fit$error
    next
  }
  expected <- manifest[manifest_index, , drop = FALSE]
  required_manifest <- c("task_id", "dataset", "axis", "lambda", "seed")
  if (is.null(fit$manifest) || !all(required_manifest %in% names(fit$manifest))) {
    validation_message <- "Stored manifest is missing required fields."
  } else {
    for (field in required_manifest) {
      if (!identical(as.character(fit$manifest[[field]]), as.character(expected[[field]]))) {
        validation_message <- paste("Stored manifest differs for", field)
      }
    }
  }
  dataset_tree <- datasets[[expected$dataset]]$tree
  if (is.null(fit$tree) || !identical(fit$tree$edge, dataset_tree$edge) ||
      !identical(fit$tree$tip.label, dataset_tree$tip.label) ||
      length(fit$tree$edge.length) != length(dataset_tree$edge.length) ||
      max(abs(fit$tree$edge.length - dataset_tree$edge.length)) > 1e-10) {
    validation_message <- "Fit tree does not match prepared tree in exact edge order."
  }
  if (length(fit$edge_sigma2) != nrow(dataset_tree$edge) ||
      any(!is.finite(fit$edge_sigma2)) || any(fit$edge_sigma2 <= 0)) {
    validation_message <- "Branch BM rates are missing, invalid or nonpositive."
  }
  if (!is.null(fit$input_md5)) {
    for (input_file in names(input_md5)) {
      if (!input_file %in% names(fit$input_md5) ||
          !identical(unname(fit$input_md5[input_file]), unname(input_md5[input_file]))) {
        validation_message <- paste("Input checksum mismatch:", input_file)
      }
    }
  }
  if (!identical(validation_message, "")) {
    fit_status$status[manifest_index] <- "invalid"
    fit_status$message[manifest_index] <- validation_message
    next
  }
  fit_status$status[manifest_index] <- "success"
  fit_status$success[manifest_index] <- TRUE
  fit_status$message[manifest_index] <- "Validated tree, rates, manifest and available checksums."
  if (!is.null(fit$convergence)) fit_status$convergence[manifest_index] <- fit$convergence
  if (!is.null(fit$task_summary$elapsed_seconds)) {
    fit_status$elapsed_seconds[manifest_index] <- fit$task_summary$elapsed_seconds
  }
  fits[[manifest_index]] <- fit
}
complete <- all(fit_status$success)
groups <- unique(manifest[, c("dataset", "lambda")])
for (group_index in seq_len(nrow(groups))) {
  group_rows <- which(manifest$dataset == groups$dataset[group_index] & manifest$lambda == groups$lambda[group_index])
  required_axes <- colnames(datasets[[groups$dataset[group_index]]]$traits)
  group_complete <- all(fit_status$success[group_rows]) &&
    !anyDuplicated(manifest$axis[group_rows]) && setequal(manifest$axis[group_rows], required_axes)
  fit_status$dataset_lambda_complete[group_rows] <- group_complete
}
fit_status$comparison_complete <- complete
write.csv(fit_status, "results/fit-status.csv", row.names = FALSE)

# Prepare normalized map draws, summaries and rank matrices once per tree.
onto_cache <- list()
for (projection_id in names(projections)) {
  projection <- projections[[projection_id]]
  branches <- projection$branches
  root <- setdiff(projection$tree$edge[, 1], projection$tree$edge[, 2])
  branches$root_edge <- branches$parent %in% root
  branches$n_descendant_families <- integer(nrow(branches))
  for (edge_index in seq_len(nrow(branches))) {
    labels <- strsplit(branches$descendant_tips[edge_index], "|", fixed = TRUE)[[1]]
    families <- tip_metadata$Family_Sharkey[match(labels, tip_metadata$analysis_label)]
    if (anyNA(families)) stop("Missing family metadata in a projected branch.")
    branches$n_descendant_families[edge_index] <- length(unique(families))
  }
  branches$interfamily_edge <- branches$n_descendant_families >= 2L
  parent_edge <- match(branches$parent, branches$child)
  branches$n_parent_descendant_families <- branches$n_descendant_families[parent_edge]
  all_families <- tip_metadata$Family_Sharkey[match(projection$tree$tip.label, tip_metadata$analysis_label)]
  branches$n_parent_descendant_families[is.na(parent_edge)] <- length(unique(all_families))
  branches$within_single_family <- branches$n_descendant_families == 1L &
    branches$n_parent_descendant_families == 1L
  branches$family_stem_edge <- branches$n_descendant_families == 1L &
    branches$n_parent_descendant_families >= 2L
  edge_subsets <- list(all_edges = seq_len(nrow(branches)),
    exclude_root_edges = which(!branches$root_edge),
    interfamily_edges = which(branches$interfamily_edge))
  bin_breaks <- projection$marginal$breaks_Ma
  bin_data <- data.frame(younger_age_Ma = head(bin_breaks, -1),
    older_age_Ma = tail(bin_breaks, -1),
    exposure_lineage_Myr = projection$marginal$exposure_lineage_Myr)
  bin_data$width_Myr <- bin_data$older_age_Ma - bin_data$younger_age_Ma
  occupied <- which(bin_data$exposure_lineage_Myr > 0)
  # Do not give the 0.0259-Myr oldest remainder the weight of a full bin.
  # All bins remain in the exported curves and exact exposure checks.
  full_width <- which(bin_data$exposure_lineage_Myr > 0 & bin_data$width_Myr >= 5)
  temporal_subsets <- list(occupied_bins = occupied, full_width_bins = full_width)
  partition_cache <- list()
  for (partition in partitions) {
    raw <- projection$samples[[partition]]
    normalization <- colSums(raw * branches$length_Myr) / sum(branches$length_Myr)
    valid_draw <- is.finite(normalization) & normalization > 0
    relative <- matrix(NA_real_, nrow(raw), ncol(raw), dimnames = dimnames(raw))
    relative[, valid_draw] <- sweep(raw[, valid_draw, drop = FALSE], 2, normalization[valid_draw], "/")
    contracted_raw <- projection$marginal$contracted_branch_samples[[partition]]
    segment_raw <- projection$marginal$original_segment_samples[[partition]]
    contracted_relative <- matrix(NA_real_, nrow(contracted_raw), ncol(contracted_raw))
    segment_relative <- contracted_relative
    contracted_relative[, valid_draw] <- sweep(contracted_raw[, valid_draw, drop = FALSE], 2, normalization[valid_draw], "/")
    segment_relative[, valid_draw] <- sweep(segment_raw[, valid_draw, drop = FALSE], 2, normalization[valid_draw], "/")
    branch_ranks <- list()
    for (scope in names(edge_subsets)) branch_ranks[[scope]] <- rank_draws(relative[edge_subsets[[scope]], , drop = FALSE])
    temporal_ranks <- list()
    for (scope in names(temporal_subsets)) {
      rows <- temporal_subsets[[scope]]
      temporal_ranks[[paste0("contracted_", scope)]] <- rank_draws(contracted_relative[rows, , drop = FALSE])
      temporal_ranks[[paste0("original_segments_", scope)]] <- rank_draws(segment_relative[rows, , drop = FALSE])
    }
    partition_cache[[partition]] <- list(raw_summary = sample_summary(raw),
      relative_summary = sample_summary(relative), relative_samples = relative,
      normalization = normalization, n_valid_maps = sum(valid_draw),
      contracted_raw_summary = sample_summary(contracted_raw),
      contracted_relative_summary = sample_summary(contracted_relative),
      segment_raw_summary = sample_summary(segment_raw),
      segment_relative_summary = sample_summary(segment_relative),
      branch_ranks = branch_ranks, temporal_ranks = temporal_ranks)
  }
  onto_cache[[projection_id]] <- list(branches = branches, edge_subsets = edge_subsets,
    bins = bin_data, temporal_subsets = temporal_subsets, partitions = partition_cache)
}

# Raw per-axis BM rates are summed before normalization for TOTAL_VARIANCE.
# This is the fitted trace of retained-PC diffusion, an exploratory summary;
# it is not an independent-axis uncertainty calculation.
haika_estimates <- list()
for (manifest_index in which(fit_status$success)) {
  row <- manifest[manifest_index, ]
  estimate_id <- paste(row$dataset, row$axis, row$lambda, sep = "::")
  haika_estimates[[estimate_id]] <- list(dataset = row$dataset, axis = row$axis,
    lambda = row$lambda, task_ids = row$task_id, n_axes_combined = 1L,
    group_complete = fit_status$dataset_lambda_complete[manifest_index],
    edge_sigma2 = fits[[manifest_index]]$edge_sigma2,
    node_sigma2 = fits[[manifest_index]]$node_sigma2)
}
for (group_index in seq_len(nrow(groups))) {
  dataset_id <- groups$dataset[group_index]
  penalty <- groups$lambda[group_index]
  rows <- which(manifest$dataset == dataset_id & manifest$lambda == penalty)
  if (!all(fit_status$dataset_lambda_complete[rows])) next
  total_variance <- rep(0, nrow(datasets[[dataset_id]]$tree$edge))
  for (manifest_index in rows) total_variance <- total_variance + fits[[manifest_index]]$edge_sigma2
  estimate_id <- paste(dataset_id, "TOTAL_VARIANCE", penalty, sep = "::")
  haika_estimates[[estimate_id]] <- list(dataset = dataset_id, axis = "TOTAL_VARIANCE",
    lambda = penalty, task_ids = manifest$task_id[rows], n_axes_combined = length(rows),
    group_complete = TRUE, edge_sigma2 = total_variance, node_sigma2 = NULL)
}

branch_tables <- list()
marginal_tables <- list()
association_tables <- list()
association_draws <- list()
table_index <- 0L
association_index <- 0L
for (estimate_id in names(haika_estimates)) {
  estimate <- haika_estimates[[estimate_id]]
  dataset <- datasets[[estimate$dataset]]
  projection_id <- projection_ids[estimate$dataset]
  projection <- projections[[projection_id]]
  cache <- onto_cache[[projection_id]]
  haika_raw <- as.numeric(estimate$edge_sigma2)
  haika_normalization <- sum(haika_raw * cache$branches$length_Myr) / sum(cache$branches$length_Myr)
  haika_relative <- haika_raw / haika_normalization
  exposure <- projection$marginal$contracted_branch_exposure
  haika_marginal_raw <- as.numeric(exposure %*% haika_raw) / cache$bins$exposure_lineage_Myr
  haika_marginal_raw[cache$bins$exposure_lineage_Myr == 0] <- NA_real_
  haika_marginal_relative <- haika_marginal_raw / haika_normalization
  haika_estimates[[estimate_id]]$normalization <- haika_normalization
  haika_estimates[[estimate_id]]$relative <- haika_relative
  haika_estimates[[estimate_id]]$marginal_raw <- haika_marginal_raw
  haika_estimates[[estimate_id]]$marginal_relative <- haika_marginal_relative
  for (partition in partitions) {
    onto <- cache$partitions[[partition]]
    common <- data.frame(dataset = estimate$dataset, source = dataset$metadata$source,
      subset = dataset$metadata$subset, matching_level = dataset$metadata$matching_level,
      axis = estimate$axis, lambda = estimate$lambda, n_axes_combined = estimate$n_axes_combined,
      task_ids = paste(estimate$task_ids, collapse = "|"),
      ontophylo_partition = partition, n_ontophylo_characters = length(projection$character_ids[[partition]]),
      n_ontophylo_maps = ncol(projection$samples[[partition]]),
      n_valid_normalized_maps = onto$n_valid_maps,
      haika_normalization = haika_normalization,
      dataset_lambda_complete = estimate$group_complete, comparison_complete = complete,
      haika_units = "PC score squared per Myr", ontophylo_units = "summed character transitions per Myr",
      relative_units = "rate / lineage-duration-weighted treewide mean",
      uncertainty = "Ontophylo map variation only; conditional on HAIKA point fit, fixed models, tree and proxy mapping",
      stringsAsFactors = FALSE)
    branch_table <- cbind(common[rep(1L, nrow(cache$branches)), , drop = FALSE], cache$branches,
      haika_raw = haika_raw, haika_relative = haika_relative)
    for (statistic in c("mean", "lower", "median", "upper")) {
      branch_table[[paste0("ontophylo_raw_", statistic)]] <- onto$raw_summary[[statistic]]
      branch_table[[paste0("ontophylo_relative_", statistic)]] <- onto$relative_summary[[statistic]]
    }
    marginal_table <- cbind(common[rep(1L, nrow(cache$bins)), , drop = FALSE], cache$bins,
      haika_raw = haika_marginal_raw, haika_relative = haika_marginal_relative)
    for (statistic in c("mean", "lower", "median", "upper")) {
      marginal_table[[paste0("ontophylo_contracted_raw_", statistic)]] <- onto$contracted_raw_summary[[statistic]]
      marginal_table[[paste0("ontophylo_contracted_relative_", statistic)]] <- onto$contracted_relative_summary[[statistic]]
      marginal_table[[paste0("ontophylo_original_segments_raw_", statistic)]] <- onto$segment_raw_summary[[statistic]]
      marginal_table[[paste0("ontophylo_original_segments_relative_", statistic)]] <- onto$segment_relative_summary[[statistic]]
    }
    table_index <- table_index + 1L
    branch_tables[[table_index]] <- branch_table
    marginal_tables[[table_index]] <- marginal_table

    # All correlations are descriptive. No independent-edge/time-bin P values.
    scopes <- c(names(cache$edge_subsets), paste0("contracted_", names(cache$temporal_subsets)),
      paste0("original_segments_", names(cache$temporal_subsets)))
    for (scope in scopes) {
      domain <- "branches"
      if (scope %in% names(cache$edge_subsets)) {
        rows <- cache$edge_subsets[[scope]]
        haika_values <- haika_relative[rows]
        onto_mean <- onto$relative_summary$mean[rows]
        ranks <- onto$branch_ranks[[scope]]
      } else {
        domain <- "time_bins"
        temporal_scope <- sub("^(contracted_|original_segments_)", "", scope)
        rows <- cache$temporal_subsets[[temporal_scope]]
        haika_values <- haika_marginal_relative[rows]
        if (startsWith(scope, "contracted_")) {
          onto_mean <- onto$contracted_relative_summary$mean[rows]
        } else {
          onto_mean <- onto$segment_relative_summary$mean[rows]
        }
        ranks <- onto$temporal_ranks[[scope]]
      }
      association <- rank_association(haika_values, onto_mean, ranks)
      output_scope <- scope
      if (scope == "contracted_full_width_bins") output_scope <- "time_full_width_bins"
      if (scope == "original_segments_full_width_bins") output_scope <- "original_segments_time_full_width_bins"
      scope_definition <- "All induced-tree edges."
      if (scope == "exclude_root_edges") scope_definition <- "All induced-tree edges excluding branches directly descending from the retained MRCA."
      if (scope == "interfamily_edges") scope_definition <- "Edges whose child clade contains at least two Family_Sharkey categories; excludes stems entering single families and terminal branches."
      if (domain == "time_bins") scope_definition <- "Occupied 10-Myr grid bins, including the short oldest remainder; equal weight of bins after lineage-exposure weighting within bins."
      if (grepl("full_width_bins$", scope)) scope_definition <- "Occupied time bins at least 5 Myr wide; excludes the short oldest remainder; equal weight of bins after lineage-exposure weighting within bins."
      association_index <- association_index + 1L
      association_tables[[association_index]] <- cbind(common, domain = domain, scope = output_scope,
        scope_definition = scope_definition,
        n_units = length(rows), spearman_rho = association$rho,
        map_rho_lower = association$lower, map_rho_median = association$median,
        map_rho_upper = association$upper, n_finite_map_correlations = association$n_valid_draws,
        interpretation = "Descriptive association; correlated branches/bins; no sampling-test P value")
      draw_id <- paste(estimate_id, partition, output_scope, sep = "::")
      association_draws[[draw_id]] <- association$draws
    }
  }
}

branch_output <- data.frame(dataset = character(), comparison_complete = logical())
marginal_output <- branch_output
association_output <- branch_output
if (length(branch_tables) > 0L) branch_output <- do.call(rbind, branch_tables)
if (length(marginal_tables) > 0L) marginal_output <- do.call(rbind, marginal_tables)
if (length(association_tables) > 0L) association_output <- do.call(rbind, association_tables)
write.csv(branch_output, "results/branch-comparison.csv", row.names = FALSE)
write.csv(marginal_output, "results/marginal-comparison.csv", row.names = FALSE)
write.csv(association_output, "results/association-summary.csv", row.names = FALSE)

# Keep reusable samples/projections, excluding cached rank matrices.
for (projection_id in names(onto_cache)) {
  for (partition in partitions) {
    onto_cache[[projection_id]]$partitions[[partition]]$branch_ranks <- NULL
    onto_cache[[projection_id]]$partitions[[partition]]$temporal_ranks <- NULL
  }
}
comparison_objects <- list(complete = complete, final_required = final_required,
  expected_fits = nrow(manifest), successful_fits = sum(fit_status$success),
  fit_status = fit_status, dataset_projection_ids = projection_ids,
  ontophylo_projections = projections, ontophylo_summaries = onto_cache,
  haika_estimates = haika_estimates, association_draws = association_draws,
  normalization = "Each HAIKA point fit and each ontophylo map divided by its own lineage-duration-weighted treewide mean",
  total_variance = "Sum raw sigma2 over all retained PCs before normalization; exploratory fitted trace, without independent-axis uncertainty",
  temporal_comparison = "Exact 10-Myr lineage-exposure bins of branch-constant rates on the same induced tree; original-segment ontophylo sensitivity kept separately",
  uncertainty = "Ontophylo reconstruction-map variation conditional on fixed HAIKA point fits, models, tree, matching and data selection; not joint rate uncertainty",
  input_md5 = input_md5, session_info = sessionInfo())
comparison_objects$datasets <- list()
for (dataset_id in names(datasets)) {
  dataset_estimates <- names(haika_estimates)[startsWith(names(haika_estimates), paste0(dataset_id, "::"))]
  comparison_objects$datasets[[dataset_id]] <- list(
    tree = datasets[[dataset_id]]$tree, metadata = datasets[[dataset_id]]$metadata,
    projection_id = unname(projection_ids[dataset_id]), estimate_ids = dataset_estimates)
}
saveRDS(comparison_objects, "results/comparison-objects.rds", compress = "xz")
cat("Successful fits:", sum(fit_status$success), "of", nrow(manifest),
  "; comparison complete:", complete, "\n")
cat("Branch rows:", nrow(branch_output), "; marginal rows:", nrow(marginal_output),
  "; association rows:", nrow(association_output), "\n")
if (final_required && !complete) {
  stop("Final comparison requires every manifest fit to be successful and valid. See results/fit-status.csv.")
}
