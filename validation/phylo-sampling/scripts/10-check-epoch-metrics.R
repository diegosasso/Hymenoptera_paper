# Check temporal metrics independently from saved branch rates and tree ages.
# This supplements the raw-map, interval and marginal-curve audit in script 07.
design <- readRDS("design/design.rds")
manifest <- read.csv("design/task-manifest.csv")
maximum_error <- 0
checked_groups <- 0L
for (task in manifest$task) {
  result <- readRDS(sprintf("results/task-%04d/result.rds", task))
  lower_age <- design$epochs$lower_ma[result$epoch]
  upper_age <- design$epochs$upper_ma[result$epoch]
  for (group_name in names(result$groups)) {
    group <- result$groups[[group_name]]
    tree_name <- if (group_name == "study_all") "study" else "full"
    tree <- design[[paste0(tree_name, "_tree")]]
    ages <- design[[paste0(tree_name, "_age")]]
    parent_age <- ages[tree$edge[, 1]]
    child_age <- ages[tree$edge[, 2]]
    overlap <- pmax(0, pmin(parent_age, upper_age) - pmax(child_age, lower_age))
    stopifnot(all(overlap <= tree$edge.length + 1e-8), sum(overlap) > 0)
    counts <- group$rate_samples * tree$edge.length
    total_changes <- colSums(counts)
    epoch_changes <- colSums(group$rate_samples * overlap)
    epoch_rate <- epoch_changes / sum(overlap)
    epoch_fraction <- epoch_changes / total_changes
    mean_age <- colSums(counts * (parent_age + child_age) / 2) / total_changes
    origins <- group$branch_summary$origin_count
    reference_rate <- sum(origins * overlap / tree$edge.length) / sum(overlap)
    reference_fraction <- sum(origins * overlap / tree$edge.length) / sum(origins)
    reference_age <- sum(origins * (parent_age + child_age) / 2) / sum(origins)
    errors <- c(total_changes - group$metrics$changes,
      epoch_rate - group$metrics$target_epoch_rate,
      epoch_fraction - group$metrics$target_epoch_fraction,
      mean_age - group$metrics$mean_change_age,
      reference_rate - group$reference$target_epoch_rate,
      reference_fraction - group$reference$target_epoch_fraction,
      reference_age - group$reference$mean_change_age,
      colSums(abs(counts - origins)) / sum(origins) - group$metrics$branch_count_error)
    stopifnot(all(is.finite(errors)))
    maximum_error <- max(maximum_error, abs(errors))
    checked_groups <- checked_groups + 1L
  }
}
stopifnot(nrow(manifest) == 1500L, checked_groups == 4500L, maximum_error < 1e-9)
lines <- c("PASS: epoch metrics and single-origin references independently verified.",
  paste("Paired trials:", nrow(manifest), "; tree/character groups:", checked_groups),
  paste("Maximum numerical disagreement:", maximum_error),
  "Branch overlap with each target epoch was computed directly from node ages.",
  "Checked changes, epoch rates/fractions, mean ages, references and branch-count errors.")
writeLines(lines, "provenance/epoch-metric-verification.txt")
cat(lines, sep = "\n")
