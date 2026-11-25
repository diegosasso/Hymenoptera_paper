# Run from comparisons/haika-wings.
# Family proxies retain the original 77-tip topology; they are not species data.
library(ape)

tree <- readRDS("data/tree-77-morphology-labels.rds")
tips <- read.csv("data/tree-tip-metadata.csv", stringsAsFactors = FALSE)
crosswalk <- read.csv("data/haika-family-proxy-crosswalk.csv", stringsAsFactors = FALSE)
source_rows <- read.csv("data/haika-family-proxy-source-rows.csv", stringsAsFactors = FALSE)
stopifnot(identical(tree$tip.label, crosswalk$analysis_label))
stopifnot(identical(tree$tip.label, tips$analysis_label))

# A deterministic family-representative subset preserves original branch ages.
# Repeat the subset with the last representative to check exemplar dependence.
family_names <- sort(unique(crosswalk$selected_family))
family_first <- character(length(family_names))
family_last <- character(length(family_names))
for (family_index in seq_along(family_names)) {
  members <- sort(crosswalk$analysis_label[crosswalk$selected_family == family_names[family_index]])
  family_first[family_index] <- members[1]
  family_last[family_index] <- members[length(members)]
}
representatives <- data.frame(family = family_names, first = family_first,
                              last = family_last, stringsAsFactors = FALSE)
write.csv(representatives, "data/family-representatives.csv", row.names = FALSE)

# Each source is kept in its own feature space. PCA is fitted to family means
# with each family weighted once, not to 77 repeated family proxy observations.
source_files <- c(
  seed3407 = "data/haika-case2-seed3407-features64.csv",
  seed46567 = "data/haika-case2-seed46567-features64.csv",
  seed91507 = "data/haika-case2-seed91507-features64.csv",
  seed92779 = "data/haika-case2-seed92779-features64.csv",
  seed97081 = "data/haika-case2-seed97081-features64.csv",
  released28 = "data/haika-case2-features.csv"
)
datasets <- list()
pca_models <- list()
variance_rows <- list()
summary_rows <- list()
row_index <- 0L
for (source_index in seq_along(source_files)) {
  source_id <- names(source_files)[source_index]
  source_data <- read.csv(source_files[source_index], stringsAsFactors = FALSE)
  feature_columns <- startsWith(names(source_data), "feature_")
  image_features <- as.matrix(source_data[, feature_columns, drop = FALSE])
  rownames(image_features) <- source_data$taxon_label
  stopifnot(!anyDuplicated(rownames(image_features)), all(is.finite(image_features)))
  family_features <- matrix(NA_real_, nrow = length(family_names), ncol = ncol(image_features),
                            dimnames = list(family_names, colnames(image_features)))
  for (family_index in seq_along(family_names)) {
    selected <- source_rows$selected_family == family_names[family_index]
    source_labels <- unique(source_rows$source_S8_label[selected])
    if (source_id == "released28") {
      source_labels <- unique(source_rows$source_S3_label[selected])
    }
    stopifnot(length(source_labels) > 0, all(source_labels %in% rownames(image_features)))
    family_features[family_index, ] <- colMeans(image_features[source_labels, , drop = FALSE])
  }
  variable_features <- apply(family_features, 2, var) > 0
  family_features <- family_features[, variable_features, drop = FALSE]
  pca <- prcomp(family_features, center = TRUE, scale. = FALSE)
  variance_fraction <- pca$sdev^2 / sum(pca$sdev^2)
  retained_axes <- max(2L, which(cumsum(variance_fraction) >= 0.90)[1])
  stopifnot(retained_axes <= 10)
  pca_models[[source_id]] <- list(pca = pca, family_features = family_features,
                                 retained_axes = retained_axes,
                                 rule = "At least two axes and at least 90 percent variance")
  variance_rows[[source_id]] <- data.frame(source = source_id,
    axis = paste0("PC", seq_along(variance_fraction)), variance_fraction = variance_fraction,
    cumulative_fraction = cumsum(variance_fraction), retained = seq_along(variance_fraction) <= retained_axes)
  tip_traits <- pca$x[crosswalk$selected_family, seq_len(retained_axes), drop = FALSE]
  rownames(tip_traits) <- tree$tip.label
  subset_names <- c("full77", "family48_first")
  if (source_id == "seed3407") subset_names <- c(subset_names, "family48_last")
  for (subset_id in subset_names) {
    selected_tips <- tree$tip.label
    if (subset_id == "family48_first") selected_tips <- family_first
    if (subset_id == "family48_last") selected_tips <- family_last
    subset_tree <- keep.tip(tree, selected_tips)
    subset_traits <- tip_traits[subset_tree$tip.label, , drop = FALSE]
    dataset_id <- paste(source_id, subset_id, sep = "_")
    datasets[[dataset_id]] <- list(tree = subset_tree, traits = subset_traits,
      metadata = list(source = source_id, source_file = unname(source_files[source_index]),
        subset = subset_id, matching_level = "family_proxy",
        original_tip_count = 77L, observed_family_profiles = length(family_names),
        retained_axes = retained_axes, variance_retained = sum(variance_fraction[seq_len(retained_axes)]),
        biological_interpretation = "New PCA of image-derived family means; not published shape/venation axes"))
    row_index <- row_index + 1L
    summary_rows[[row_index]] <- data.frame(dataset = dataset_id, source = source_id,
      subset = subset_id, n_tips = length(subset_tree$tip.label),
      n_families = length(unique(crosswalk$selected_family[match(subset_tree$tip.label, crosswalk$analysis_label)])),
      n_distinct_profiles = nrow(unique(subset_traits)), n_axes = retained_axes,
      variance_retained = sum(variance_fraction[seq_len(retained_axes)]),
      stringsAsFactors = FALSE)
    trait_output <- data.frame(analysis_label = subset_tree$tip.label, subset_traits)
    write.csv(trait_output, paste0("data/traits-", dataset_id, ".csv"), row.names = FALSE)
  }
}
saveRDS(datasets, "data/prepared-datasets.rds", compress = "xz")
saveRDS(pca_models, "data/family-pca-models.rds", compress = "xz")
write.csv(do.call(rbind, variance_rows), "data/pca-variance.csv", row.names = FALSE)
write.csv(do.call(rbind, summary_rows), "data/dataset-summary.csv", row.names = FALSE)

# Primary seed 3407 gets the full prespecified smoothing grid. Other network
# initializations and the released 28-column matrix are sensitivity at lambda=1.
manifest_rows <- list()
task_id <- 0L
for (dataset_id in names(datasets)) {
  dataset <- datasets[[dataset_id]]
  penalties <- 1
  if (dataset$metadata$source == "seed3407") penalties <- c(0.1, 1, 10)
  for (axis in colnames(dataset$traits)) {
    for (penalty in penalties) {
      task_id <- task_id + 1L
      manifest_rows[[task_id]] <- data.frame(task_id = task_id, dataset = dataset_id,
        axis = axis, lambda = penalty, seed = 200000L + task_id, stringsAsFactors = FALSE)
    }
  }
}
manifest <- do.call(rbind, manifest_rows)
write.csv(manifest, "data/fit-manifest.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "provenance/preparation-session.txt")
print(do.call(rbind, summary_rows))
cat("Prepared", nrow(manifest), "bounded fits.\n")
