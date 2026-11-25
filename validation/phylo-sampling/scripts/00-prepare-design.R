# Run from validation/phylo-sampling. Tree files and outputs are relative to here.
library(ape)
n_trials <- 100L
proportions <- c(0.05, 0.10, 0.30, 0.50, 0.80)
epoch_lower <- c(0, 80, 160)
epoch_upper <- c(80, 160, 280)
simulation_seed <- 20261120L

full_tree <- reorder.phylo(read.tree("inputs/HYMtopC_bamm_hisse.phy"), "cladewise")
taxa <- readRDS("inputs/taxon_data.RDS")
study_tree <- reorder.phylo(keep.tip(full_tree, as.character(taxa$Blaimer_taxa)), "cladewise")
saved_study <- readRDS("inputs/hym_tree.RDS")
stopifnot(Ntip(full_tree) == 765L, Ntip(study_tree) == 77L,
          is.ultrametric(full_tree), is.ultrametric(study_tree),
          isTRUE(all.equal.phylo(study_tree, saved_study)))
full_tree$node.label <- NULL
study_tree$node.label <- NULL
apocrita_labels <- as.character(taxa$Blaimer_taxa[taxa$Tag != "NonApo"])
apocrita_node <- getMRCA(full_tree, apocrita_labels)
apocrita_tips <- extract.clade(full_tree, apocrita_node)$tip.label
stopifnot(length(apocrita_tips) == 736L, length(apocrita_labels) == 61L,
          !any(taxa$Blaimer_taxa[taxa$Tag == "NonApo"] %in% apocrita_tips))

# For cladewise-ordered trees, reversing the edges visits children before parents.
descendant_matrix <- function(tree) {
  n_tips <- Ntip(tree)
  answer <- matrix(FALSE, n_tips + tree$Nnode, n_tips)
  answer[cbind(seq_len(n_tips), seq_len(n_tips))] <- TRUE
  for (edge in rev(seq_len(nrow(tree$edge)))) {
    parent <- tree$edge[edge, 1]
    child <- tree$edge[edge, 2]
    answer[parent, ] <- answer[parent, ] | answer[child, ]
  }
  stopifnot(all(answer[n_tips + 1L, ]))
  answer
}
full_descendants <- descendant_matrix(full_tree)
study_descendants <- descendant_matrix(study_tree)
study_positions <- match(study_tree$tip.label, full_tree$tip.label)
root_age <- max(node.depth.edgelength(full_tree))
full_age <- root_age - node.depth.edgelength(full_tree)
study_age <- root_age - node.depth.edgelength(study_tree)
non_apocrita <- !full_tree$tip.label %in% apocrita_tips
non_apocrita_study <- !study_tree$tip.label %in% apocrita_labels

# Only branches wholly inside Apocrita can receive imposed changes.
# Exclude its crown node: the incoming stem is outside the crown clade.
inside_apocrita <- rowSums(full_descendants[, non_apocrita, drop = FALSE]) == 0L
internal_nodes <- seq.int(Ntip(full_tree) + 1L, Ntip(full_tree) + full_tree$Nnode)
eligible <- internal_nodes[inside_apocrita[internal_nodes] & internal_nodes != apocrita_node]
nodes <- data.frame(node = eligible, age_ma = full_age[eligible],
                    character = sprintf("N%04d", eligible),
                    edge = match(eligible, full_tree$edge[, 2]),
                    n_full_ones = rowSums(full_descendants[eligible, , drop = FALSE]),
                    n_study_ones = rowSums(full_descendants[eligible, study_positions, drop = FALSE]))
nodes$epoch <- findInterval(nodes$age_ma, c(0, 80, 160, 280))
nodes$variable_study <- nodes$n_study_ones > 0L & nodes$n_study_ones < Ntip(study_tree)
stopifnot(nrow(nodes) == 734L, all(nodes$epoch %in% 1:3),
          all(inside_apocrita[full_tree$edge[nodes$edge, 1]]))

# Exact correspondence between full branches and merged study-tree paths.
study_signatures <- character(nrow(study_tree$edge))
for (edge in seq_len(nrow(study_tree$edge))) {
  tips <- which(study_descendants[study_tree$edge[edge, 2], ])
  study_signatures[edge] <- paste(tips, collapse = ",")
}
full_to_study <- rep(NA_integer_, nrow(full_tree$edge))
for (edge in seq_len(nrow(full_tree$edge))) {
  tips <- which(full_descendants[full_tree$edge[edge, 2], study_positions])
  if (length(tips) > 0L) {
    full_to_study[edge] <- match(paste(tips, collapse = ","), study_signatures)
    stopifnot(!is.na(full_to_study[edge]))
  }
}
for (edge in seq_len(nrow(study_tree$edge))) {
  stopifnot(abs(sum(full_tree$edge.length[full_to_study %in% edge]) -
                study_tree$edge.length[edge]) < 1e-7)
}
nodes$study_edge <- full_to_study[nodes$edge]
nodes$node_preserved <- FALSE
retained <- which(nodes$variable_study)
study_child <- study_tree$edge[nodes$study_edge[retained], 2]
nodes$node_preserved[retained] <- abs(study_age[study_child] - nodes$age_ma[retained]) < 1e-7

full_characters <- matrix("0", Ntip(full_tree), nrow(nodes),
                          dimnames = list(full_tree$tip.label, nodes$character))
for (i in seq_len(nrow(nodes))) {
  full_characters[full_descendants[nodes$node[i], ], i] <- "1"
}
study_characters <- full_characters[study_tree$tip.label, , drop = FALSE]
stopifnot(all(full_characters[non_apocrita, ] == "0"),
          all(study_characters[non_apocrita_study, ] == "0"))

epochs <- data.frame(epoch = 1:3, lower_ma = epoch_lower, upper_ma = epoch_upper,
                     n_nodes = 0L, n_variable_study = 0L, n_nodes_preserved = 0L)
for (i in 1:3) {
  in_epoch <- nodes$epoch == i
  epochs$n_nodes[i] <- sum(in_epoch)
  epochs$n_variable_study[i] <- sum(nodes$variable_study[in_epoch])
  epochs$n_nodes_preserved[i] <- sum(nodes$node_preserved[in_epoch])
}

# Exact branch exposure inside each requested epoch.
branch_exposure <- function(tree, ages) {
  exposure <- matrix(0, nrow(tree$edge), 3L)
  for (i in 1:3) {
    exposure[, i] <- pmax(0, pmin(ages[tree$edge[, 1]], epoch_upper[i]) -
                            pmax(ages[tree$edge[, 2]], epoch_lower[i]))
  }
  exposure
}
full_exposure <- branch_exposure(full_tree, full_age)
study_exposure <- branch_exposure(study_tree, study_age)
full_apo_edges <- inside_apocrita[full_tree$edge[, 1]]
inside_study_apo <- rowSums(study_descendants[, non_apocrita_study, drop = FALSE]) == 0L
study_apo_edges <- inside_study_apo[study_tree$edge[, 1]]

dir.create("design", showWarnings = FALSE)
dir.create("simulated", showWarnings = FALSE)
conditions <- data.frame()
manifest <- data.frame()
condition <- 0L
task <- 0L
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
for (epoch in 1:3) {
  candidates <- which(nodes$epoch == epoch)
  for (proportion in proportions) {
    condition <- condition + 1L
    n_characters <- max(1L, as.integer(floor(length(candidates) * proportion + 0.5)))
    conditions <- rbind(conditions, data.frame(condition = condition, epoch = epoch,
      proportion = proportion, n_eligible = length(candidates), n_characters = n_characters,
      realized_proportion = n_characters / length(candidates)))
    for (trial in seq_len(n_trials)) {
      task <- task + 1L
      seed <- simulation_seed + condition * 10000L + trial
      set.seed(seed)
      chosen <- candidates[sample.int(length(candidates), n_characters, replace = FALSE)]
      selected <- nodes[chosen, , drop = FALSE]
      rownames(selected) <- NULL
      dataset <- list(task = task, condition = condition, epoch = epoch, trial = trial,
        seed = seed, proportion = proportion, selected_nodes = selected)
      filename <- sprintf("simulated/task-%04d.rds", task)
      if (file.exists(filename)) stopifnot(identical(readRDS(filename), dataset))
      if (!file.exists(filename)) saveRDS(dataset, filename)
      manifest <- rbind(manifest, data.frame(task = task, condition = condition,
        epoch = epoch, proportion = proportion, trial = trial, seed = seed,
        n_characters = n_characters, n_variable_study = sum(selected$variable_study),
        file = filename))
    }
  }
}
design <- list(full_tree = full_tree, study_tree = study_tree,
  full_age = full_age, study_age = study_age, nodes = nodes, epochs = epochs,
  conditions = conditions, n_trials = n_trials, root_age = root_age,
  apocrita_node = apocrita_node, apocrita_crown_age = full_age[apocrita_node],
  full_to_study = full_to_study, full_exposure = full_exposure, study_exposure = study_exposure,
  full_apo_edges = full_apo_edges, study_apo_edges = study_apo_edges,
  age_grid = seq(0, 280, by = 1), seed = simulation_seed)
saveRDS(design, "design/design.rds")
saveRDS(list(full = full_characters, study = study_characters), "design/node-characters.rds")
write.csv(manifest, "design/task-manifest.csv", row.names = FALSE)
write.csv(conditions, "design/conditions.csv", row.names = FALSE)
write.csv(epochs, "design/epoch-availability.csv", row.names = FALSE)
write.csv(nodes, "design/eligible-nodes.csv", row.names = FALSE)
write.csv(data.frame(taxon = full_tree$tip.label, apocrita = !non_apocrita,
                     study = full_tree$tip.label %in% study_tree$tip.label),
          "design/taxon-selection.csv", row.names = FALSE)
capture.output(sessionInfo(), file = "provenance/design-session.txt")
print(epochs)
print(conditions)
cat("Prepared", nrow(manifest), "paired trials;", sum(manifest$n_characters),
    "full-tree character instances.\n")

