# Conditional homogeneous-BM null with the same family-proxy assignment.
# Run on Roihu: Rscript --vanilla scripts/06-null-comparison.R REPLICATE
# This tests association beyond shared tree geometry and proxy duplication.
arguments <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.integer(arguments[1])
stopifnot(length(replicate_id) == 1, is.finite(replicate_id), replicate_id >= 1)
.libPaths(c("/projappl/project_2004996/phylo-sampling/Rlib",
            "/projappl/project_2004996/hymenoptera-node-pattern-bias/Rlib", .libPaths()))
library(ape)
source("scripts/bm-model.R")
source("scripts/project-ontophylo.R")
datasets <- readRDS("data/prepared-datasets.rds")
crosswalk <- read.csv("data/haika-family-proxy-crosswalk.csv", stringsAsFactors = FALSE)
onto <- readRDS("data/ontophylo-wing-branch-rate-samples.rds")
reference <- datasets$seed3407_family48_first
reference_tree <- reference$tree
traits <- reference$traits
stopifnot(identical(rownames(traits), reference_tree$tip.label))

# Estimate the constant multivariate BM covariance with a fitted root mean.
# Cross-PC evolutionary covariance is retained in simulation, rather than
# treating ordinary PCA axes as phylogenetically independent replicates.
covariance <- vcv(reference_tree)
covariance <- covariance[reference_tree$tip.label, reference_tree$tip.label]
inverse <- solve(covariance)
root_mean <- colSums(inverse %*% traits) / sum(inverse)
centered <- sweep(traits, 2, root_mean, "-")
trait_covariance <- crossprod(centered, inverse %*% centered) / nrow(traits)
stopifnot(all(eigen(trait_covariance, symmetric = TRUE, only.values = TRUE)$values > 0))
set.seed(910000L + replicate_id)
noise <- matrix(rnorm(nrow(traits) * ncol(traits)), nrow(traits), ncol(traits))
simulated <- t(chol(covariance)) %*% noise %*% chol(trait_covariance)
simulated <- sweep(simulated, 2, root_mean, "+")
rownames(simulated) <- rownames(traits)
colnames(simulated) <- colnames(traits)
reference_families <- crosswalk$selected_family[match(rownames(simulated), crosswalk$analysis_label)]
stopifnot(!anyNA(reference_families), !anyDuplicated(reference_families))
rownames(simulated) <- reference_families

output <- list(replicate = replicate_id, generating_model = "Homogeneous multivariate BM on 48-family representative tree, then duplicate profiles to 77 tips",
               trait_covariance = trait_covariance, source = "seed3407", lambda = 1,
               results = list(), success = TRUE)
for (dataset_id in c("seed3407_full77", "seed3407_family48_first")) {
  dataset <- datasets[[dataset_id]]
  current_families <- crosswalk$selected_family[match(dataset$tree$tip.label, crosswalk$analysis_label)]
  stopifnot(!anyNA(current_families), all(current_families %in% rownames(simulated)))
  simulated_traits <- simulated[current_families, , drop = FALSE]
  rownames(simulated_traits) <- dataset$tree$tip.label
  fits <- list()
  axis_rates <- matrix(NA_real_, nrow(dataset$tree$edge), ncol(simulated_traits),
                       dimnames = list(NULL, colnames(simulated_traits)))
  for (axis_index in seq_len(ncol(simulated_traits))) {
    x <- simulated_traits[, axis_index]
    names(x) <- rownames(simulated_traits)
    fitted <- fit_bm_bounded(dataset$tree, x, lambda = 1, n_starts = 3L,
                             seed = 920000L + 100L * replicate_id + axis_index,
                             maxit = 2000L)
    fits[[axis_index]] <- list(success = fitted$success, diagnostics = fitted$diagnostics)
    if (isTRUE(fitted$success)) axis_rates[, axis_index] <- fitted$edge_sigma2
    if (!isTRUE(fitted$success)) output$success <- FALSE
  }
  output$results[[dataset_id]] <- list(tree = dataset$tree, axis_rates = axis_rates,
                                       diagnostics = fits)
}
output$session_info <- sessionInfo()
output$input_md5 <- tools::md5sum(c("data/prepared-datasets.rds", "scripts/bm-model.R",
                                   "scripts/bm-cached-objective.R",
                                   "data/haika-family-proxy-crosswalk.csv",
                                   "scripts/06-null-comparison.R"))
dir.create("results/null", recursive = TRUE, showWarnings = FALSE)
saveRDS(output, sprintf("results/null/rep-%04d.rds", replicate_id))
cat("Null replicate", replicate_id, "success", output$success, "\n")
if (!isTRUE(output$success)) quit(save = "no", status = 1L)
