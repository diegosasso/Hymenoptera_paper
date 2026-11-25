# Extended Data: Neff+ versus family-level PyRate and full-tree BAMM.
# Run with Source in RStudio. Edit the working directory if needed.
# Uses saved reconstructions; no character fitting or BAMM analysis is rerun.
setwd(path.expand("~/GitHub/Hymenoptera_paper"))
library(ape)
library(BAMMtools)

output_dir <- "validation/neff-diversification-uncertainty"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
rebuild_inputs <- FALSE  # Set TRUE after changing the metric or time bins below.
cache_file <- file.path(output_dir, "binned-trajectories.rds")

# 1. Input files and common 20-Myr bins ----
pyrate_files <- c(
  "Te_Ts-RJMCMC4_0__BDCS12_BDS_mcmc.log", "Te_Ts-RJMCMC4_1__BDCS12_BDS_mcmc.log",
  "Te_Ts-RJMCMC4_2__BDCS12_BDS_mcmc.log", "Te_Ts-RJMCMC4_3__BDCS12_BDS_mcmc.log",
  "Te_Ts-RJMCMC4_4__BDCS12_BDS_mcmc.log", "Te_Ts-RJMCMC4_5__BDCS12_BDS_mcmc.log",
  "Te_Ts-RJMCMC4_6__BDCS12_BDS_mcmc.log", "Te_Ts-RJMCMC4_7__BDCS12_BDS_mcmc.log",
  "Te_Ts-RJMCMC4_8__BDCS12_BDS_mcmc.log", "Te_Ts-RJMCMC4_9__BDCS12_BDS_mcmc.log")
pyrate_paths <- file.path("data_pyrate/BDCS/BDCS-fam-ssing-20", pyrate_files)
input_files <- c("data/hym_tree.RDS", "data_out/br_rates_all.RDS",
  "data_out/neff-per-replicate.RDS", "data_out/pyrate-dive.RDS",
  "data/bamm.rds", "make-BAMM/data/bamm_data.RDA", pyrate_paths)
stopifnot(all(file.exists(input_files)))
input_hashes <- tools::md5sum(input_files)
bins <- data.frame(bin = 1:13, older = seq(260, 20, -20), younger = seq(240, 0, -20))
bins$age <- (bins$older + bins$younger) / 2

# Reuse the expensive extraction only when all input files still match.
if (file.exists(cache_file) && !rebuild_inputs) {
  inputs <- readRDS(cache_file)
  stopifnot(identical(inputs$input_hashes, input_hashes), identical(inputs$bins, bins))
  cat("Using saved binned trajectories.\n")
} else {
  # 2. Neff+ for each of the 1000 existing branch-rate map replicates ----
  # This is the empirical >= 95th-percentile rule, across 15 adult units.
  tree <- readRDS("data/hym_tree.RDS")
  branch_rates <- readRDS("data_out/br_rates_all.RDS")
  regions <- names(branch_rates)[names(branch_rates) != "larval"]
  n_maps <- ncol(branch_rates[[1]])
  stopifnot(length(regions) == 15, n_maps == 1000, Ntip(tree) == 77)
  depth <- node.depth.edgelength(tree)
  edge_start <- depth[tree$edge[, 1]]
  edge_end <- depth[tree$edge[, 2]]
  times <- seq(0, max(depth), length.out = 1000)
  times[1000] <- times[1000] - 1e-4
  ages <- rev(times)  # Preserve the original manuscript's age convention.
  active_edges <- vector("list", 1000)
  slice_keys <- character(1000)
  for (i in 1:1000) {
    active_edges[[i]] <- which(times[i] >= edge_start & times[i] <= edge_end)
    slice_keys[i] <- paste(active_edges[[i]], collapse = ",")
  }
  # Adjacent slices with the same branches have identical Neff; evaluate once.
  slice_group <- match(slice_keys, unique(slice_keys))
  first_slice <- match(unique(slice_keys), slice_keys)
  neff_slices <- matrix(NA_real_, 1000, n_maps)
  for (j in 1:n_maps) {
    patterns <- numeric(nrow(tree$edge))
    for (k in seq_along(regions)) {
      rates <- branch_rates[[regions[k]]][, j]
      cutoff <- quantile(rates, 0.95, names = FALSE)
      patterns <- patterns + (rates >= cutoff) * 2^(k - 1)
    }
    neff <- numeric(length(first_slice))
    for (g in seq_along(first_slice)) {
      states <- patterns[active_edges[[first_slice[g]]]]
      accelerated <- states[states != 0]
      if (length(accelerated) > 0) {
        probabilities <- as.numeric(table(accelerated)) / length(accelerated)
        entropy <- -sum(probabilities * log(probabilities))
        frequency <- length(accelerated) / length(states)
        neff[g] <- frequency * exp(entropy)
      }
    }
    neff_slices[, j] <- neff[slice_group]
  }
  neff_draws <- matrix(NA_real_, 13, n_maps)
  for (i in 1:13) {
    keep <- ages >= bins$younger[i] & ages < bins$older[i]
    neff_draws[i, ] <- colMeans(neff_slices[keep, , drop = FALSE])
  }

  # 3. PyRate family-level net diversification, lambda - mu ----
  # Retain the last 200 rows per file, exactly as in make-pyrate-data.Rmd.
  pyrate_draws <- matrix(NA_real_, 13, 2000)
  pyrate_index <- data.frame()
  for (j in 1:10) {
    log_data <- read.delim(pyrate_paths[j], check.names = FALSE)
    stopifnot(nrow(log_data) >= 200)
    rows <- (nrow(log_data) - 199):nrow(log_data)
    retained <- log_data[rows, ]
    stopifnot(min(retained$root_age) > 260)
    lambda <- as.matrix(retained[, paste0("lambda_", 1:13)])
    mu <- as.matrix(retained[, paste0("mu_", 1:13)])
    columns <- ((j - 1) * 200 + 1):(j * 200)
    pyrate_draws[, columns] <- t(lambda - mu)
    pyrate_index <- rbind(pyrate_index, data.frame(draw = columns, file = pyrate_paths[j],
      row = rows, iteration = retained$it))
  }

  # 4. BAMM net diversification for every retained posterior configuration ----
  # Use the full 765-tip tree, as in the existing data/bamm.rds summary.
  # This extraction takes a few minutes; the binned result is cached afterwards.
  cat("Extracting all BAMM trajectories on the full tree...\n")
  bamm_environment <- new.env()
  load("make-BAMM/data/bamm_data.RDA", envir = bamm_environment)
  shifts <- bamm_environment$shifts
  stopifnot(length(shifts$tip.label) == 765, length(shifts$eventData) == 18001)
  rtt <- getRateThroughTimeMatrix(shifts, nslices = 280)
  bamm_net <- rtt$lambda - rtt$mu  # Keep lambda and mu from the same sample.
  bamm_ages <- rev(unname(rtt$times))  # Matrix columns run from root to present.
  bamm_draws <- matrix(NA_real_, 13, nrow(bamm_net))
  for (i in 1:13) {
    keep <- bamm_ages >= bins$younger[i] & bamm_ages < bins$older[i]
    bamm_draws[i, ] <- rowMeans(bamm_net[, keep, drop = FALSE])
  }
  bamm_saved <- readRDS("data/bamm.rds")
  bamm_error <- max(abs(rev(colMeans(bamm_net)) - bamm_saved$rate))
  stopifnot(bamm_error < 1e-10)
  colnames(neff_draws) <- paste0("map_", 1:ncol(neff_draws))
  colnames(pyrate_draws) <- paste0("pyrate_", 1:ncol(pyrate_draws))
  colnames(bamm_draws) <- paste0("bamm_", 1:ncol(bamm_draws))
  inputs <- list(bins = bins, neff = neff_draws, PyRate = pyrate_draws,
    BAMM = bamm_draws, pyrate_index = pyrate_index, bamm_mean_error = bamm_error,
    input_hashes = input_hashes)
  saveRDS(inputs, cache_file)
  rm(bamm_environment, shifts, rtt, bamm_net, branch_rates, neff_slices)
  gc()
}

# 5. Verify the extracted means against the original manuscript summaries ----
neff_draws <- inputs$neff
diversification <- list(PyRate = inputs$PyRate, BAMM = inputs$BAMM)
neff_saved <- readRDS("data_out/neff-per-replicate.RDS")
pyrate_saved <- readRDS("data_out/pyrate-dive.RDS")
neff_reference <- numeric(13)
for (i in 1:13) {
  keep <- neff_saved$time >= bins$younger[i] & neff_saved$time < bins$older[i]
  neff_reference[i] <- mean(neff_saved$mean_entropy[keep])
}
neff_error <- max(abs(rowMeans(neff_draws) - neff_reference))
pyrate_error <- max(abs(rowMeans(inputs$PyRate) - pyrate_saved$mean[match(bins$bin, pyrate_saved$bin)]))
validation <- data.frame(check = c("Neff mean maximum error", "PyRate mean maximum error",
  "BAMM mean maximum error"), value = c(neff_error, pyrate_error, inputs$bamm_mean_error))
stopifnot(all(validation$value < 1e-10))

# 6. Pearson correlation for ALL pairs of complete trajectories ----
# cor() compares columns: 2,000,000 PyRate and 18,001,000 BAMM pairs.
# The 13 bins stay together. These combinations are not independent experiments.
# Two-sided p-values use df = 13 - 2; they do not correct temporal autocorrelation.
results <- data.frame()
r_histograms <- list()
p_histograms <- list()
histogram_data <- data.frame()
for (method in names(diversification)) {
  rates <- diversification[[method]]
  stopifnot(all(is.finite(neff_draws)), all(is.finite(rates)))
  r <- cor(neff_draws, rates)
  p <- 2 * pt(-abs(r) * sqrt(11 / (1 - r^2)), df = 11)
  stopifnot(all(is.finite(r)), all(is.finite(p)), all(p > 0 & p <= 1))
  mean_test <- cor.test(rowMeans(neff_draws), rowMeans(rates))
  r_interval <- quantile(r, c(0.025, 0.5, 0.975), names = FALSE)
  p_interval <- quantile(p, c(0.025, 0.5, 0.975), names = FALSE)
  results <- rbind(results, data.frame(method = method, neff_maps = ncol(neff_draws),
    rate_draws = ncol(rates), pairs = length(r), mean_curve_r = unname(mean_test$estimate),
    mean_curve_p = mean_test$p.value, r_lower = r_interval[1], r_median = r_interval[2],
    r_upper = r_interval[3], p_lower = p_interval[1], p_median = p_interval[2],
    p_upper = p_interval[3], positive_percent = 100 * mean(r > 0),
    p_less_005_n = sum(p < 0.05), p_less_005_percent = 100 * mean(p < 0.05),
    positive_p_less_005_n = sum(r > 0 & p < 0.05),
    negative_p_less_005_n = sum(r < 0 & p < 0.05)))
  r_histograms[[method]] <- hist(as.vector(r), breaks = seq(-1, 1, 0.04), plot = FALSE)
  p_min <- floor(log10(min(p, mean_test$p.value)))
  p_histograms[[method]] <- hist(log10(as.vector(p)), breaks = seq(p_min, 0, 0.25), plot = FALSE)
  for (statistic in c("r", "log10_p")) {
    h <- r_histograms[[method]]
    if (statistic == "log10_p") h <- p_histograms[[method]]
    stopifnot(sum(h$counts) == length(r))
    histogram_data <- rbind(histogram_data, data.frame(method = method, statistic = statistic,
      lower = head(h$breaks, -1), upper = tail(h$breaks, -1), count = h$counts,
      percent = 100 * h$counts / length(r)))
  }
  # Check the matrix formula against direct cor.test calls, including extremes.
  set.seed(20261124)
  pairs <- unique(c(1, length(r), which.min(r), which.max(r), sample(length(r), 100)))
  index <- arrayInd(pairs, dim(r))
  r_error <- p_error <- 0
  for (j in seq_along(pairs)) {
    test <- cor.test(neff_draws[, index[j, 1]], rates[, index[j, 2]])
    r_error <- max(r_error, abs(unname(test$estimate) - r[pairs[j]]))
    p_error <- max(p_error, abs(test$p.value - p[pairs[j]]))
  }
  stopifnot(r_error < 1e-12, p_error < 1e-12)
  validation <- rbind(validation, data.frame(check = paste(method,
    c("direct cor.test r maximum error", "direct cor.test p maximum error")),
    value = c(r_error, p_error)))
  cat(method, ":", length(r), "complete trajectory pairs evaluated.\n")
}
write.csv(results, file.path(output_dir, "correlation-summary.csv"), row.names = FALSE)
write.csv(histogram_data, file.path(output_dir, "figure-source-data.csv"), row.names = FALSE)
write.csv(validation, file.path(output_dir, "validation.csv"), row.names = FALSE)
write.csv(inputs$pyrate_index, file.path(output_dir, "pyrate-draw-index.csv"), row.names = FALSE)

# 7. One four-panel figure, in the project's readable Arial style ----
# The PDF is editable; the TIFF is RGB at 300 dpi, 180 x 145 mm.
blue <- "#0072B2"
grey <- "#444444"
orange <- "#E69F00"
titles <- c("PyRate (families)", "BAMM")
p_axis_min <- min(histogram_data$lower[histogram_data$statistic == "log10_p"])
p_ticks <- seq(p_axis_min, 0, by = 2)
if (tail(p_ticks, 1) != 0) p_ticks <- c(p_ticks, 0)
for (format in c("pdf", "png", "tiff")) {
  filename <- file.path(output_dir, paste0("extended-data-neff-diversification.", format))
  if (format == "pdf") cairo_pdf(filename, width = 180 / 25.4, height = 145 / 25.4,
                                 family = "Arial", pointsize = 10)
  if (format == "png") png(filename, width = 180, height = 145, units = "mm",
                            res = 300, type = "cairo", family = "Arial", pointsize = 10)
  if (format == "tiff") tiff(filename, width = 180, height = 145, units = "mm",
    res = 300, type = "cairo", family = "Arial", pointsize = 10, compression = "lzw")
  par(mfrow = c(2, 2), mar = c(3.6, 4.1, 2.4, 0.8), oma = c(2, 0, 0, 0),
      cex = 1, cex.axis = 0.9, cex.lab = 1, cex.main = 1, font.main = 1,
      mgp = c(2.5, 0.7, 0), tcl = -0.2, las = 1, bty = "l", lwd = 2 / 3)
  for (i in 1:2) {
    method <- results$method[i]
    h <- r_histograms[[method]]
    h$counts <- 100 * h$counts / results$pairs[i]
    top <- max(h$counts) * 1.4
    plot(h, freq = TRUE, col = adjustcolor(blue, 0.75), border = NA,
      xlim = c(-1, 1), ylim = c(0, top), main = titles[i],
      xlab = expression("Pearson " * italic(r)), ylab = "Reconstruction pairs (%)")
    segments(results$mean_curve_r[i], 0, results$mean_curve_r[i], top * 0.74,
             col = grey, lty = 2, lwd = 1)
    text(-0.98, top * 0.94, sprintf("Median %.2f", results$r_median[i]), adj = 0, cex = 0.9)
    text(-0.98, top * 0.84, sprintf("95%% range %.2f to %.2f",
      results$r_lower[i], results$r_upper[i]), adj = 0, cex = 0.9)
    mtext(c("a", "c")[i], side = 3, line = 0.8, adj = -0.16, font = 2, cex = 1.1)
    h <- p_histograms[[method]]
    h$counts <- 100 * h$counts / results$pairs[i]
    top <- max(h$counts) * 1.25
    plot(h, freq = TRUE, col = adjustcolor(blue, 0.75), border = NA, xaxt = "n",
      xlim = c(p_axis_min, 0), ylim = c(0, top), main = titles[i],
      xlab = expression("Nominal Pearson " * italic(p)), ylab = "Reconstruction pairs (%)")
    axis(1, at = p_ticks, labels = parse(text = paste0("10^", p_ticks)))
    segments(log10(results$mean_curve_p[i]), 0, log10(results$mean_curve_p[i]), top * 0.82,
             col = grey, lty = 2, lwd = 1)
    segments(log10(0.05), 0, log10(0.05), top * 0.82, col = orange, lty = 3, lwd = 1.2)
    text(p_axis_min + 0.1, top * 0.94, sprintf("p < 0.05: %.2f%%",
      results$p_less_005_percent[i]), adj = 0, cex = 0.9)
    mtext(c("b", "d")[i], side = 3, line = 0.8, adj = -0.16, font = 2, cex = 1.1)
  }
  par(xpd = NA)
  legend(x = grconvertX(0.5, "ndc", "user"), y = grconvertY(0.02, "ndc", "user"),
    legend = c("Mean curves", "p = 0.05"), col = c(grey, orange), lty = c(2, 3),
    lwd = c(1, 1.2), horiz = TRUE, xjust = 0.5, yjust = 0, bty = "n", cex = 0.9)
  dev.off()
}

# 8. Provenance: verify that original inputs were not changed ----
stopifnot(identical(tools::md5sum(input_files), input_hashes))
write.csv(data.frame(file = input_files, md5 = unname(input_hashes)),
          file.path(output_dir, "input-checksums.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output_dir, "session-info.txt"))
print(results, row.names = FALSE)
cat("Figure, source data and cached trajectories:", output_dir, "\n")
