# Run on Roihu from comparisons/haika-wings (or its copied project root).
# Rscript --vanilla scripts/09-bound-sensitivity.R
# Eighteen sensitivity fits: two datasets, three PCs, offsets 5 and 15,
# plus a prespecified strong-smoothing case with lambda = original root age.
# The six original offset-10 fits are read without modification.
project_library <- Sys.getenv("HAIKA_RLIB", "/projappl/project_2004996/haika-wings/Rlib")
.libPaths(c(project_library, "/projappl/project_2004996/phylo-sampling/Rlib",
  "/projappl/project_2004996/hymenoptera-node-pattern-bias/Rlib", .libPaths()))
library(ape)
source("scripts/bm-model.R")
source("scripts/project-ontophylo.R")
if (!"lower_bound_offset" %in% names(formals(fit_bm_bounded))) {
  stop("The fitting helper must expose lower_bound_offset, defaulting to 10.")
}
dir.create("results/bounds", recursive = TRUE, showWarnings = FALSE)
datasets <- readRDS("data/prepared-datasets.rds")
manifest <- read.csv("data/fit-manifest.csv", stringsAsFactors = FALSE)
full_ontophylo <- readRDS("data/ontophylo-wing-branch-rate-samples.rds")
tips <- read.csv("data/tree-tip-metadata.csv", stringsAsFactors = FALSE)
dataset_ids <- c("seed3407_full77", "seed3407_family48_first")
original_root_age <- max(ape::node.depth.edgelength(full_ontophylo$tree))
scenarios <- data.frame(scenario = c("primary_reference", "lower_offset5",
  "lower_offset15", "strong_smoothing_root_age"),
  lower_bound_offset = c(10, 5, 15, 10), lambda = c(1, 1, 1, original_root_age),
  refit = c(FALSE, TRUE, TRUE, TRUE), stringsAsFactors = FALSE)
fit_results <- list()
axis_rows <- list()
start_rows <- list()
summary_rows <- list()
total_profiles <- list()
axis_index <- 0L
summary_index <- 0L
n_refits <- 0L
any_failed <- FALSE

spearman <- function(first, second) {
  keep <- is.finite(first) & is.finite(second)
  if (sum(keep) < 3L || length(unique(first[keep])) < 2L || length(unique(second[keep])) < 2L) {
    return(NA_real_)
  }
  cor(first[keep], second[keep], method = "spearman")
}

for (dataset_id in dataset_ids) {
  dataset <- datasets[[dataset_id]]
  axes <- colnames(dataset$traits)
  stopifnot(length(axes) == 3L)
  tree <- dataset$tree
  projection <- project_ontophylo(full_ontophylo, tree$tip.label)
  stopifnot(identical(tree$tip.label, projection$tree$tip.label),
    identical(tree$edge, projection$tree$edge),
    max(abs(tree$edge.length - projection$tree$edge.length)) < 1e-10)
  branches <- projection$branches
  root <- setdiff(tree$edge[, 1], tree$edge[, 2])
  n_families <- integer(nrow(branches))
  for (edge_index in seq_len(nrow(branches))) {
    labels <- strsplit(branches$descendant_tips[edge_index], "|", fixed = TRUE)[[1]]
    families <- tips$Family_Sharkey[match(labels, tips$analysis_label)]
    stopifnot(!anyNA(families))
    n_families[edge_index] <- length(unique(families))
  }
  edge_subsets <- list(all_edges = seq_len(nrow(branches)),
    exclude_root_edges = which(!branches$parent %in% root),
    interfamily_edges = which(n_families >= 2L))
  bin_breaks <- projection$marginal$breaks_Ma
  bin_widths <- diff(bin_breaks)
  bin_exposure <- projection$marginal$exposure_lineage_Myr
  full_width_bins <- which(bin_exposure > 0 & bin_widths >= 5)

  # Match the primary collector's point summaries exactly: normalize every
  # ontophylo map, then take its across-map mean on branches and in time bins.
  onto_raw <- projection$samples$fore_wing
  onto_normalization <- colSums(onto_raw * tree$edge.length) / sum(tree$edge.length)
  stopifnot(all(is.finite(onto_normalization)), all(onto_normalization > 0))
  onto_relative <- sweep(onto_raw, 2, onto_normalization, "/")
  onto_mean <- rowMeans(onto_relative)
  onto_time_raw <- projection$marginal$contracted_branch_samples$fore_wing
  onto_time_relative <- sweep(onto_time_raw, 2, onto_normalization, "/")
  onto_time_mean <- rowMeans(onto_time_relative)
  dataset_profiles <- list()

  for (scenario_index in seq_len(nrow(scenarios))) {
    scenario <- scenarios$scenario[scenario_index]
    offset <- scenarios$lower_bound_offset[scenario_index]
    penalty <- scenarios$lambda[scenario_index]
    group_fits <- list()
    group_diagnostics <- list()
    for (axis in axes) {
      manifest_row <- manifest[manifest$dataset == dataset_id &
        manifest$axis == axis & manifest$lambda == 1, , drop = FALSE]
      stopifnot(nrow(manifest_row) == 1L)
      trait <- dataset$traits[, axis]
      names(trait) <- rownames(dataset$traits)
      fit_id <- paste(dataset_id, scenario, axis, sep = "::")
      started <- Sys.time()
      elapsed_started <- proc.time()[["elapsed"]]
      error_message <- ""
      if (!scenarios$refit[scenario_index]) {
        fit_file <- sprintf("results/fits/task-%04d.rds", manifest_row$task_id)
        fit <- readRDS(fit_file)
        if (!isTRUE(fit$success) || !identical(fit$tree$edge, tree$edge) ||
            !identical(fit$tree$tip.label, tree$tip.label) ||
            max(abs(fit$tree$edge.length - tree$edge.length)) > 1e-10 ||
            as.character(fit$manifest$dataset) != dataset_id ||
            as.character(fit$manifest$axis) != axis || fit$manifest$lambda != 1) {
          stop("Primary comparison fit is absent, failed or mismatched: ", fit_file)
        }
        origin <- "Existing primary fit; not refitted"
      } else {
        n_refits <- n_refits + 1L
        cat("Refit", n_refits, "of 18:", dataset_id, axis,
          "lower offset", offset, "lambda", penalty, "seed", manifest_row$seed, "\n")
        flush.console()
        fit <- tryCatch(fit_bm_bounded(tree, trait, lambda = penalty, n_starts = 3L,
          seed = manifest_row$seed, maxit = 2000L, lower_bound_offset = offset),
          error = function(error) error)
        if (inherits(fit, "error")) {
          error_message <- conditionMessage(fit)
          fit <- list(success = FALSE, error = error_message)
        }
        origin <- "New sensitivity fit, same seed and starts as primary"
      }
      elapsed <- proc.time()[["elapsed"]] - elapsed_started
      success <- isTRUE(fit$success)
      bound_center <- log(var(trait) / max(ape::node.depth.edgelength(tree)))
      if (success) {
        if (abs((bound_center - fit$lower_log_bound) - offset) > 1e-8 ||
            abs((fit$upper_log_bound - bound_center) - 10) > 1e-8) {
          success <- FALSE
          error_message <- "Returned numerical bounds do not match requested lower offset and fixed upper +10."
        }
      }
      if (!success) {
        any_failed <- TRUE
        if (error_message == "") {
          error_message <- if (is.null(fit$error)) "Fit failed convergence or likelihood/objective validation." else fit$error
        }
      }
      gradient <- rep(NA_real_, length(tree$tip.label) + tree$Nnode)
      n_lower <- n_upper <- n_positive_lower <- n_negative_upper <- NA_integer_
      lower_grad_min <- lower_grad_median <- lower_grad_max <- interior_max_abs <- NA_real_
      best_start <- convergence <- NA_integer_
      n_converged_starts <- NA_integer_
      lower_hits_all_starts <- upper_hits_all_starts <- NA_integer_
      log_likelihood <- objective <- lower_bound <- upper_bound <- NA_real_
      likelihood_passed <- objective_passed <- FALSE
      if (success) {
        cached_objective <- make_bm_objective(tree, trait[tree$tip.label], lambda = penalty)
        log_node_rates <- log(fit$node_sigma2)
        gradient <- cached_objective$gr(log_node_rates)
        names(gradient) <- names(fit$node_sigma2)
        at_lower <- log_node_rates - fit$lower_log_bound < 1e-4
        at_upper <- fit$upper_log_bound - log_node_rates < 1e-4
        n_lower <- sum(at_lower)
        n_upper <- sum(at_upper)
        n_positive_lower <- sum(gradient[at_lower] > 1e-5)
        n_negative_upper <- sum(gradient[at_upper] < -1e-5)
        if (any(at_lower)) {
          lower_grad_min <- min(gradient[at_lower])
          lower_grad_median <- median(gradient[at_lower])
          lower_grad_max <- max(gradient[at_lower])
        }
        if (any(!at_lower & !at_upper)) interior_max_abs <- max(abs(gradient[!at_lower & !at_upper]))
        best_start <- fit$best_start
        convergence <- fit$convergence
        n_converged_starts <- sum(fit$diagnostics$convergence == 0, na.rm = TRUE)
        lower_hits_all_starts <- sum(fit$diagnostics$n_near_lower_bound, na.rm = TRUE)
        upper_hits_all_starts <- sum(fit$diagnostics$n_near_upper_bound, na.rm = TRUE)
        log_likelihood <- fit$logLik
        objective <- fit$penalized_objective
        lower_bound <- fit$lower_log_bound
        upper_bound <- fit$upper_log_bound
        likelihood_passed <- isTRUE(fit$likelihood_check$passed)
        objective_passed <- isTRUE(fit$objective_check$passed)
      }
      fit$bound_sensitivity <- list(dataset = dataset_id, axis = axis, scenario = scenario,
        lower_bound_offset = offset, upper_bound_offset = 10, lambda = penalty,
        seed = manifest_row$seed, original_task_id = manifest_row$task_id,
        origin = origin, success = success, error = error_message,
        objective_gradient_at_selected_node_rates = gradient,
        gradient_interpretation = "Objective minimized: positive derivative at lower bound favors smaller unconstrained node rates",
        started_utc = format(started, tz = "UTC", usetz = TRUE), elapsed_seconds = elapsed)
      fit_results[[fit_id]] <- fit
      group_fits[[axis]] <- fit
      axis_index <- axis_index + 1L
      axis_row <- data.frame(dataset = dataset_id, subset = dataset$metadata$subset,
        axis = axis, scenario = scenario, lambda = penalty, lower_bound_offset = offset, upper_bound_offset = 10,
        primary_task_id = manifest_row$task_id, seed = manifest_row$seed, origin = origin,
        success = success, error = error_message, selected_start = best_start,
        convergence = convergence, n_converged_starts = n_converged_starts,
        selected_lower_bound_hits = n_lower, selected_upper_bound_hits = n_upper,
        lower_bound_hits_all_starts = lower_hits_all_starts, upper_bound_hits_all_starts = upper_hits_all_starts,
        n_lower_gradients_above_1e_minus5 = n_positive_lower,
        n_upper_gradients_below_minus1e_minus5 = n_negative_upper,
        lower_gradient_min = lower_grad_min, lower_gradient_median = lower_grad_median,
        lower_gradient_max = lower_grad_max, interior_max_abs_gradient = interior_max_abs,
        lower_log_bound = lower_bound, upper_log_bound = upper_bound,
        logLik = log_likelihood, penalized_objective = objective,
        likelihood_check_passed = likelihood_passed, objective_check_passed = objective_passed,
        elapsed_seconds = elapsed, stringsAsFactors = FALSE)
      axis_rows[[axis_index]] <- axis_row
      group_diagnostics[[axis]] <- axis_row
      if (!is.null(fit$diagnostics)) {
        start_diagnostics <- fit$diagnostics
        start_diagnostics$dataset <- dataset_id
        start_diagnostics$axis <- axis
        start_diagnostics$lambda <- penalty
        start_diagnostics$scenario <- scenario
        start_diagnostics$lower_bound_offset <- offset
        start_diagnostics$origin <- origin
        # Select common start fields because primary and sensitivity wrappers
        # may attach different additional metadata to their diagnostics.
        fields <- c("dataset", "axis", "scenario", "lambda", "lower_bound_offset", "origin",
          "start", "objective", "convergence", "elapsed_seconds",
          "n_near_lower_bound", "n_near_upper_bound", "message")
        start_rows[[length(start_rows) + 1L]] <- start_diagnostics[, fields, drop = FALSE]
      }
    }

    diagnostics <- do.call(rbind, group_diagnostics)
    group_success <- all(diagnostics$success)
    summary <- data.frame(dataset = dataset_id, subset = dataset$metadata$subset,
      matching_level = dataset$metadata$matching_level, axis = "TOTAL_VARIANCE", scenario = scenario,
      lambda = penalty, lower_bound_offset = offset, upper_bound_offset = 10,
      n_axes = length(axes), n_successful_axes = sum(diagnostics$success),
      success = group_success, n_tips = length(tree$tip.label), n_edges = nrow(tree$edge),
      n_interfamily_edges = length(edge_subsets$interfamily_edges),
      n_temporal_full_width_bins = length(full_width_bins),
      selected_lower_bound_hits = sum(diagnostics$selected_lower_bound_hits),
      selected_upper_bound_hits = sum(diagnostics$selected_upper_bound_hits),
      lower_bound_hits_all_starts = sum(diagnostics$lower_bound_hits_all_starts),
      lower_gradients_positive = sum(diagnostics$n_lower_gradients_above_1e_minus5),
      all_likelihood_checks_passed = all(diagnostics$likelihood_check_passed),
      all_objective_checks_passed = all(diagnostics$objective_check_passed),
      total_variance_normalization = NA_real_, stringsAsFactors = FALSE)
    for (scope in names(edge_subsets)) {
      summary[[paste0("rate_agreement_default_", scope)]] <- NA_real_
      summary[[paste0("forewing_rho_", scope)]] <- NA_real_
      summary[[paste0("forewing_rho_default_", scope)]] <- NA_real_
      summary[[paste0("forewing_rho_change_", scope)]] <- NA_real_
    }
    summary$rate_agreement_default_time_full_width_bins <- NA_real_
    summary$forewing_rho_time_full_width_bins <- NA_real_
    summary$forewing_rho_default_time_full_width_bins <- NA_real_
    summary$forewing_rho_change_time_full_width_bins <- NA_real_
    if (group_success) {
      total_raw <- rep(0, nrow(tree$edge))
      for (axis in axes) total_raw <- total_raw + group_fits[[axis]]$edge_sigma2
      normalization <- sum(total_raw * tree$edge.length) / sum(tree$edge.length)
      total_relative <- total_raw / normalization
      time_raw <- as.numeric(projection$marginal$contracted_branch_exposure %*% total_raw) / bin_exposure
      time_raw[bin_exposure == 0] <- NA_real_
      time_relative <- time_raw / normalization
      profile <- list(tree = tree, raw_total_sigma2 = total_raw, relative_total_sigma2 = total_relative,
        normalization = normalization, marginal_raw = time_raw, marginal_relative = time_relative,
        breaks_Ma = bin_breaks, exposure_lineage_Myr = bin_exposure,
        original_root_age_Ma = projection$original_root_age_Ma,
        retained_root_age_Ma = projection$retained_root_age_Ma)
      dataset_profiles[[scenario]] <- profile
      baseline <- dataset_profiles[["primary_reference"]]
      if (is.null(baseline)) stop("The default-bound baseline must be available before comparisons.")
      summary$total_variance_normalization <- normalization
      for (scope in names(edge_subsets)) {
        selected <- edge_subsets[[scope]]
        summary[[paste0("rate_agreement_default_", scope)]] <- spearman(total_relative[selected], baseline$relative_total_sigma2[selected])
        current_rho <- spearman(total_relative[selected], onto_mean[selected])
        baseline_rho <- spearman(baseline$relative_total_sigma2[selected], onto_mean[selected])
        summary[[paste0("forewing_rho_", scope)]] <- current_rho
        summary[[paste0("forewing_rho_default_", scope)]] <- baseline_rho
        summary[[paste0("forewing_rho_change_", scope)]] <- current_rho - baseline_rho
      }
      summary$rate_agreement_default_time_full_width_bins <- spearman(time_relative[full_width_bins], baseline$marginal_relative[full_width_bins])
      summary$forewing_rho_time_full_width_bins <- spearman(time_relative[full_width_bins], onto_time_mean[full_width_bins])
      summary$forewing_rho_default_time_full_width_bins <- spearman(baseline$marginal_relative[full_width_bins], onto_time_mean[full_width_bins])
      summary$forewing_rho_change_time_full_width_bins <- summary$forewing_rho_time_full_width_bins - summary$forewing_rho_default_time_full_width_bins
    }
    summary_index <- summary_index + 1L
    summary_rows[[summary_index]] <- summary
  }
  total_profiles[[dataset_id]] <- list(profiles = dataset_profiles,
    ontophylo_forewing_relative_mean = onto_mean,
    ontophylo_forewing_marginal_relative_mean = onto_time_mean,
    branch_subsets = edge_subsets, time_full_width_bins = full_width_bins)
}

stopifnot(n_refits == 18L)
axis_diagnostics <- do.call(rbind, axis_rows)
start_diagnostics <- do.call(rbind, start_rows)
bound_summary <- do.call(rbind, summary_rows)
write.csv(axis_diagnostics, "results/bound-axis-diagnostics.csv", row.names = FALSE)
write.csv(start_diagnostics, "results/bound-start-diagnostics.csv", row.names = FALSE)
write.csv(bound_summary, "results/bound-sensitivity.csv", row.names = FALSE)
input_files <- c("data/prepared-datasets.rds", "data/fit-manifest.csv",
  "data/ontophylo-wing-branch-rate-samples.rds", "data/tree-tip-metadata.csv",
  "scripts/bm-model.R", "scripts/bm-cached-objective.R", "scripts/project-ontophylo.R",
  "scripts/09-bound-sensitivity.R")
input_md5 <- tools::md5sum(input_files)
saveRDS(list(complete = !any_failed, n_new_fits = n_refits, scenarios = scenarios, fits = fit_results,
  total_profiles = total_profiles, summary = bound_summary, axis_diagnostics = axis_diagnostics,
  start_diagnostics = start_diagnostics, input_md5 = input_md5,
  normalization = "Sum raw retained-PC branch diffusion rates, then divide by duration-weighted treewide mean",
  association = "Descriptive Spearman against mean individually normalized forewing ontophylo map draws",
  interfamily_scope = "Child clade contains at least two Family_Sharkey categories; family stems excluded",
  time_scope = "Branch-constant exact exposure bins; only occupied bins with width at least 5 Myr enter correlations",
  bound_gradient = "Gradient of minimized objective: positive at lower bound favors a smaller unconstrained rate; nodewise flags use1e-5 threshold",
  smoothing_sensitivity = "Prespecified lambda = original tree root age, offset10; log-rate diffusion SD approximately1 over original root-to-tip duration under penalty's time scale; not selected by correlation",
  session_info = sessionInfo()), "results/bounds/fits.rds", compress = "xz")
print(bound_summary)
if (any_failed) stop("One or more sensitivity fits failed. Diagnostics and partial results were preserved.")
cat("Completed 18 new fits; all eight bound/smoothing dataset summaries valid.\n")
