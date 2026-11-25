# Run from the Hymenoptera_paper repository root.
# This extracts existing reconstructions; it does not refit any model.
library(ape)

output_dir <- "comparisons/haika-wings"
dir.create(file.path(output_dir, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "provenance"), recursive = TRUE, showWarnings = FALSE)

# Read the two label sets without altering the original files.
tree_molecular <- readRDS("data/hym_tree.RDS")
analysis_data <- new.env()
load("data_out/paramo_stm_adult_final.RDA", envir = analysis_data)
partition_data <- new.env()
load("data_out/paramo_amalg_adult.RDA", envir = partition_data)
tree <- analysis_data$hym_tree
tip_data <- as.data.frame(analysis_data$hym_data)
stopifnot(identical(tree$edge, tree_molecular$edge))
stopifnot(identical(tree$edge.length, tree_molecular$edge.length))
stopifnot(identical(tree$tip.label, tip_data$Genus_Sharkey))
stopifnot(identical(tree_molecular$tip.label, tip_data$Blaimer_taxa))
stopifnot(is.ultrametric(tree), !anyDuplicated(tree$tip.label))
tip_data$tip_id <- seq_along(tree$tip.label)
tip_data$analysis_label <- tree$tip.label
tip_data$molecular_species <- sub("^([^_]+_[^_]+).*$", "\\1", tip_data$Blaimer_taxa)
tip_data$morphology_genus_differs_from_molecular <- tip_data$Genus_Sharkey != tip_data$Genus_Blaimer
write.csv(tip_data, file.path(output_dir, "data/tree-tip-metadata.csv"), row.names = FALSE)
saveRDS(tree, file.path(output_dir, "data/tree-77-morphology-labels.rds"))
saveRDS(tree_molecular, file.path(output_dir, "data/tree-77-molecular-labels.rds"))
write.tree(tree, file.path(output_dir, "data/tree-77-morphology-labels.nwk"))
write.tree(tree_molecular, file.path(output_dir, "data/tree-77-molecular-labels.nwk"))

# Branch identifiers preserve the original edge order; descendant sets give
# stable clade identities even if another program reorders the tree.
n_edges <- nrow(tree$edge)
n_tips <- length(tree$tip.label)
node_depth <- node.depth.edgelength(tree)
tree_height <- max(node_depth[seq_len(n_tips)])
branch_data <- data.frame(
  edge_id = seq_len(n_edges),
  branch_id = paste0("BR_", seq_len(n_edges)),
  parent = tree$edge[, 1], child = tree$edge[, 2],
  length_Myr = tree$edge.length,
  older_age_Ma = tree_height - node_depth[tree$edge[, 1]],
  younger_age_Ma = tree_height - node_depth[tree$edge[, 2]],
  terminal = tree$edge[, 2] <= n_tips,
  descendant_tips = character(n_edges),
  descendant_molecular_taxa = character(n_edges),
  n_descendant_tips = integer(n_edges)
)
for (edge_index in seq_len(n_edges)) {
  descendants <- tree$edge[edge_index, 2]
  frontier <- descendants
  while (length(frontier) > 0L) {
    child_nodes <- tree$edge[tree$edge[, 1] %in% frontier, 2]
    descendants <- c(descendants, child_nodes)
    frontier <- child_nodes
  }
  descendant_tips <- descendants[descendants <= n_tips]
  branch_data$descendant_tips[edge_index] <- paste(sort(tree$tip.label[descendant_tips]), collapse = "|")
  branch_data$descendant_molecular_taxa[edge_index] <- paste(sort(tree_molecular$tip.label[descendant_tips]), collapse = "|")
  branch_data$n_descendant_tips[edge_index] <- length(descendant_tips)
}
write.csv(branch_data, file.path(output_dir, "data/tree-branch-geometry.csv"), row.names = FALSE)

# The original wing partitions include one muscle character for each wing.
forewing_ids <- partition_data$anat_ent[["fore wing"]]
hindwing_ids <- partition_data$anat_ent[["hind wing"]]
wing_ids <- partition_data$body_reg$wings
stopifnot(length(intersect(forewing_ids, hindwing_ids)) == 0L)
stopifnot(setequal(wing_ids, c(forewing_ids, hindwing_ids)))
raw_annotations <- as.data.frame(readRDS("data/char_annot.RDS"))
wing_annotations <- raw_annotations[match(wing_ids, raw_annotations$char_id), ]
stopifnot(identical(wing_annotations$char_id, wing_ids))
write.csv(wing_annotations, file.path(output_dir, "data/ontophylo-wing-character-metadata.csv"), row.names = FALSE)
wing_tip_characters <- as.data.frame(analysis_data$hym_adult_final[, c("taxa", wing_ids)])
write.csv(wing_tip_characters, file.path(output_dir, "data/ontophylo-wing-tip-characters.csv"), row.names = FALSE)
saveRDS(list(tree = tree, tip_data = tip_data, characters = wing_tip_characters,
  fitted_models = analysis_data$fitted_models_final[wing_ids],
  root_prior_method = "maddfitz",
  fitting_source = "R/helpers.R: fit_models; corHMM ER/SYM/ARD compared by AICc for independent characters; custom Q for dependent characters",
  stochastic_mapping_source = "R/helpers.R: paramo_stm; phytools::make.simmap with fixed fitted Q, root_prior, tips_prior"),
  file.path(output_dir, "data/ontophylo-wing-fitted-models.rds"))

# Re-extract count samples to check the stored branch rates and enable
# character-level sensitivity analyses. Each map has fixed Q/root priors.
stored_rates <- readRDS("data_out/br_rates_all.RDS")
n_maps <- ncol(stored_rates[["fore wing"]])
count_samples <- array(0L, dim = c(n_edges, n_maps, length(wing_ids)),
  dimnames = list(branch_data$branch_id, paste0("S_", seq_len(n_maps)), wing_ids))
model_metadata <- data.frame(char_id = wing_ids, model = character(length(wing_ids)),
  n_states = integer(length(wing_ids)), Q_min_offdiagonal = numeric(length(wing_ids)),
  Q_max_offdiagonal = numeric(length(wing_ids)))
for (character_index in seq_along(wing_ids)) {
  character_id <- wing_ids[character_index]
  maps <- readRDS(paste0("stmaps/", character_id, ".RDS"))
  stopifnot(length(maps) == n_maps)
  for (map_index in seq_len(n_maps)) {
    one_map <- maps[[map_index]]
    stopifnot(identical(one_map$edge, tree$edge))
    stopifnot(identical(one_map$edge.length, tree$edge.length))
    stopifnot(identical(one_map$tip.label, tree$tip.label))
    count_samples[, map_index, character_index] <- lengths(one_map$maps) - 1L
  }
  fit <- analysis_data$fitted_models_final[[character_id]]
  model_metadata$model[character_index] <- as.character(fit$best_model$model)[1]
  model_metadata$n_states[character_index] <- nrow(fit$Q_matrix)
  offdiagonal <- fit$Q_matrix[row(fit$Q_matrix) != col(fit$Q_matrix)]
  model_metadata$Q_min_offdiagonal[character_index] <- min(offdiagonal)
  model_metadata$Q_max_offdiagonal[character_index] <- max(offdiagonal)
}
saveRDS(count_samples, file.path(output_dir, "data/ontophylo-wing-character-count-samples.rds"))
write.csv(model_metadata, file.path(output_dir, "data/ontophylo-wing-model-metadata.csv"), row.names = FALSE)

partition_ids <- list(
  fore_wing = forewing_ids,
  hind_wing = hindwing_ids,
  all_wings = wing_ids,
  fore_wing_sclerites = wing_annotations$char_id[wing_annotations$label == "fore wing" & wing_annotations$type == "sclerite"],
  hind_wing_sclerites = wing_annotations$char_id[wing_annotations$label == "hind wing" & wing_annotations$type == "sclerite"]
)
rate_samples <- vector("list", length(partition_ids))
names(rate_samples) <- names(partition_ids)
branch_summary_list <- vector("list", length(partition_ids))
for (partition_index in seq_along(partition_ids)) {
  ids <- partition_ids[[partition_index]]
  summed_counts <- matrix(0, n_edges, n_maps, dimnames = dimnames(count_samples)[1:2])
  for (character_id in ids) {
    summed_counts <- summed_counts + count_samples[, , character_id]
  }
  one_rate <- summed_counts / tree$edge.length
  rate_samples[[partition_index]] <- one_rate
  interval <- t(apply(one_rate, 1, quantile, probs = c(0.025, 0.5, 0.975)))
  one_summary <- branch_data
  one_summary$partition <- names(partition_ids)[partition_index]
  one_summary$n_characters <- length(ids)
  one_summary$n_maps <- n_maps
  one_summary$mean <- rowMeans(one_rate)
  one_summary$lower <- interval[, 1]
  one_summary$median <- interval[, 2]
  one_summary$upper <- interval[, 3]
  one_summary$mean_per_character <- one_summary$mean / length(ids)
  one_summary$zero_count_fraction <- rowMeans(summed_counts == 0)
  branch_summary_list[[partition_index]] <- one_summary
}
fore_max_error <- max(abs(rate_samples$fore_wing - stored_rates[["fore wing"]]))
hind_max_error <- max(abs(rate_samples$hind_wing - stored_rates[["hind wing"]]))
stopifnot(fore_max_error < 1e-12, hind_max_error < 1e-12)
branch_summary <- do.call(rbind, branch_summary_list)
write.csv(branch_summary, file.path(output_dir, "data/ontophylo-wing-branch-summary.csv"), row.names = FALSE)
saveRDS(list(tree = tree, branches = branch_data, character_ids = partition_ids,
  units = "Summed reconstructed character-state transitions per Myr",
  uncertainty = "Across stochastic maps conditional on fixed fitted Q matrices, root/tip priors and tree",
  samples = rate_samples), file.path(output_dir, "data/ontophylo-wing-branch-rate-samples.rds"))

# Pointwise marginal curves: equal weight for every lineage crossing that age.
# Node boundaries use younger <= t < older; root uses its outgoing branches.
ages <- seq(0, tree_height, length.out = 1000)
active_matrix <- matrix(0, length(ages), n_edges)
for (age_index in seq_along(ages)) {
  age <- ages[age_index]
  active <- branch_data$younger_age_Ma <= age + 1e-10 & branch_data$older_age_Ma > age + 1e-10
  # The stored tree differs among terminal depths by at most 1.1e-7 Myr.
  # All 77 extant terminal branches therefore represent age zero.
  if (age_index == 1L) active <- branch_data$terminal
  if (age_index == length(ages)) active <- tree$edge[, 1] == n_tips + 1L
  active_matrix[age_index, active] <- 1 / sum(active)
}
curve_samples <- vector("list", length(rate_samples))
names(curve_samples) <- names(rate_samples)
curve_summaries <- vector("list", length(rate_samples))
for (partition_index in seq_along(rate_samples)) {
  curves <- active_matrix %*% rate_samples[[partition_index]]
  curve_samples[[partition_index]] <- curves
  interval <- t(apply(curves, 1, quantile, probs = c(0.025, 0.5, 0.975)))
  curve_summaries[[partition_index]] <- data.frame(
    age_Ma = ages, partition = names(rate_samples)[partition_index],
    n_lineages = rowSums(active_matrix > 0), mean = rowMeans(curves),
    lower = interval[, 1], median = interval[, 2], upper = interval[, 3])
}
write.csv(do.call(rbind, curve_summaries), file.path(output_dir, "data/ontophylo-wing-marginal-pointwise.csv"), row.names = FALSE)
saveRDS(list(age_Ma = ages, samples = curve_samples, weights = active_matrix),
  file.path(output_dir, "data/ontophylo-wing-marginal-pointwise-samples.rds"))

# Exact 10-Myr bins weight each branch by its overlap with the time interval.
# This is total expected changes divided by total lineage-time exposure.
breaks <- c(seq(0, floor(tree_height / 10) * 10, by = 10), tree_height)
exposure <- matrix(0, length(breaks) - 1L, n_edges)
for (bin_index in seq_len(nrow(exposure))) {
  exposure[bin_index, ] <- pmax(0,
    pmin(branch_data$older_age_Ma, breaks[bin_index + 1L]) -
    pmax(branch_data$younger_age_Ma, breaks[bin_index]))
}
bin_weights <- exposure / rowSums(exposure)
bin_samples <- vector("list", length(rate_samples))
names(bin_samples) <- names(rate_samples)
bin_summaries <- vector("list", length(rate_samples))
for (partition_index in seq_along(rate_samples)) {
  bins <- bin_weights %*% rate_samples[[partition_index]]
  bin_samples[[partition_index]] <- bins
  interval <- t(apply(bins, 1, quantile, probs = c(0.025, 0.5, 0.975)))
  bin_summaries[[partition_index]] <- data.frame(
    younger_age_Ma = head(breaks, -1), older_age_Ma = tail(breaks, -1),
    partition = names(rate_samples)[partition_index], exposure_lineage_Myr = rowSums(exposure),
    mean = rowMeans(bins), lower = interval[, 1], median = interval[, 2], upper = interval[, 3])
}
write.csv(do.call(rbind, bin_summaries), file.path(output_dir, "data/ontophylo-wing-marginal-binned-10Myr.csv"), row.names = FALSE)
saveRDS(list(breaks_Ma = breaks, samples = bin_samples, exposure_lineage_Myr = exposure),
  file.path(output_dir, "data/ontophylo-wing-marginal-binned-10Myr-samples.rds"))

# Provenance and checks. Newick writers may reorder edges when read again;
# use the saved RDS and clade identifiers for rate matching.
source_files <- c("data/hym_tree.RDS", "data/taxon_data.RDS", "data/char_annot.RDS",
  "data_out/paramo_stm_adult_final.RDA", "data_out/paramo_amalg_adult.RDA",
  "data_out/br_rates_all.RDS", "make-branch-specific-rates.Rmd",
  "make-adult-marginal-rate.Rmd", "R/utils-branch-ama.R", "R/utils-ST.R", "R/helpers.R",
  paste0("stmaps/", wing_ids, ".RDS"))
hashes <- tools::md5sum(source_files)
write.csv(data.frame(source = names(hashes), md5 = unname(hashes)),
  file.path(output_dir, "provenance/ontophylo-input-checksums.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output_dir, "provenance/ontophylo-extraction-session.txt"))
git_head <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
git_remote <- system2("git", c("remote", "get-url", "origin"), stdout = TRUE)
git_status <- system2("git", c("status", "--short"), stdout = TRUE)
writeLines(c(paste("Repository:", git_remote), paste("HEAD:", git_head), "Working-tree status at extraction:", git_status),
  file.path(output_dir, "provenance/ontophylo-repository-state.txt"))
checks <- data.frame(check = c("tips", "edges", "maps", "height_Ma", "total_lineage_Myr",
  "morphology_genus_substitutions", "forewing_stored_rate_max_error", "hindwing_stored_rate_max_error",
  "minimum_active_lineages", "exposure_minus_tree_length"),
  value = c(n_tips, n_edges, n_maps, tree_height, sum(tree$edge.length),
  sum(tip_data$morphology_genus_differs_from_molecular), fore_max_error, hind_max_error,
  min(rowSums(active_matrix > 0)), sum(exposure) - sum(tree$edge.length)))
write.csv(checks, file.path(output_dir, "provenance/ontophylo-extraction-checks.csv"), row.names = FALSE)
print(checks)
