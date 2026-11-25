# Import the completed unrestricted shuffles; no simulation is repeated.
# This is a provenance/import script. Once figure-data.rds is present, the
# figure is independent of the older results folder. run-shuffling.R can
# regenerate the same data directly from the repository's empirical inputs.

analysis_dir <- "."
if (file.exists("validation/neff-shuffling/prepare-source-data.R")) {
  analysis_dir <- "validation/neff-shuffling"
}
analysis_dir <- normalizePath(analysis_dir)
repo_dir <- normalizePath(file.path(analysis_dir, "../.."))
source_dir <- file.path(analysis_dir, "../neff-recovery/empirical-branch-shuffle/results")
output_dir <- file.path(analysis_dir, "source-data")
dir.create(output_dir, showWarnings = FALSE)
dir.create(file.path(analysis_dir, "provenance"), showWarnings = FALSE)

original <- readRDS(file.path(source_dir, "analysis.rds"))
shuffles <- readRDS(file.path(source_dir, "unrestricted-shuffles.rds"))
manifest <- read.csv(file.path(source_dir, "input-checksums.csv"))
current_md5 <- tools::md5sum(file.path(repo_dir, manifest$path))
stopifnot(identical(unname(current_md5), manifest$md5),
          shuffles$null_type == "unrestricted", shuffles$n_shuffles == 9999,
          shuffles$seed == 20261123,
          identical(shuffles$group_neff, original$null_groups$unrestricted),
          identical(shuffles$binned_neff, original$null_bins$unrestricted))

figure_data <- list(slice_age = original$slice_age, slice_group = original$group,
                    observed_neff = original$observed_group[original$group],
                    shuffle_neff = shuffles$group_neff, observed_bins = original$bins,
                    shuffle_bins = shuffles$binned_neff, bin_weights = original$bin_weights,
                    units = original$units, n_maps = 1000L, n_shuffles = shuffles$n_shuffles,
                    seed = shuffles$seed, n_tips = length(original$tree$tip.label),
                    n_edges = nrow(original$tree$edge))
saveRDS(figure_data, file.path(output_dir, "figure-data.rds"))
files <- c("input-checksums.csv", "accelerated-branch-counts.csv", "validation.csv")
for (file in files) {
  copied <- file.copy(file.path(source_dir, file), file.path(output_dir, file), overwrite = TRUE)
  stopifnot(copied)
}
files <- c("analysis.rds", "unrestricted-shuffles.rds", "null-correlations.csv")
source_md5 <- tools::md5sum(file.path(source_dir, files))
write.csv(data.frame(source = paste0("validation/neff-recovery/empirical-branch-shuffle/results/", files),
                      md5 = unname(source_md5)),
          file.path(analysis_dir, "provenance/import-checksums.csv"), row.names = FALSE)
writeLines(c("Imported the complete, previously calculated unrestricted null on 2026-11-23.",
             "All 9,999 shuffles retained; no trials selected, replaced or regenerated.",
             "No length-constrained results are included in this publication package.",
             "Original empirical input checksums verified before importing."),
           file.path(analysis_dir, "provenance/import.txt"))
cat("Imported all", shuffles$n_shuffles, "unrestricted shuffles.\n")
