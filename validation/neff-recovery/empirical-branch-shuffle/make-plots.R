# Run after run-shuffles.R, from this directory or the repository root.
library(ggplot2)
library(gridExtra)
library(grid)

analysis_dir <- "."
if (file.exists("validation/neff-recovery/empirical-branch-shuffle/make-plots.R")) {
  analysis_dir <- "validation/neff-recovery/empirical-branch-shuffle"
}
analysis_dir <- normalizePath(analysis_dir)
output_dir <- Sys.getenv("NEFF_OUTPUT_DIR", file.path(analysis_dir, "results"))
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, showWarnings = FALSE)
source(file.path(analysis_dir, "../scripts/nature-style.R"))
analysis <- readRDS(file.path(output_dir, "analysis.rds"))
null <- read.csv(file.path(output_dir, "null-correlations.csv"))
observed <- analysis$bins
raw_summary <- analysis$summary[analysis$summary$statistic == "raw", ]

# Pointwise intervals across shuffled mean curves, not across ASR samples.
trajectory_data <- data.frame()
for (null_type in names(analysis$null_bins)) {
  values <- analysis$null_bins[[null_type]]
  curve <- observed
  curve$null_type <- null_type
  curve$null_median <- NA_real_
  curve$null_lower <- NA_real_
  curve$null_upper <- NA_real_
  for (bin in 1:nrow(curve)) {
    curve$null_median[bin] <- median(values[bin, ])
    curve$null_lower[bin] <- quantile(values[bin, ], 0.025)
    curve$null_upper[bin] <- quantile(values[bin, ], 0.975)
  }
  trajectory_data <- rbind(trajectory_data, curve)
}
write.csv(trajectory_data, file.path(output_dir, "figure-trajectories.csv"),
          row.names = FALSE)

# Common plotting functions change presentation only.
plot_trajectory <- function(null_type, title, tag) {
  curve <- trajectory_data[trajectory_data$null_type == null_type, ]
  ggplot(curve, aes(x = age)) +
    geom_ribbon(aes(ymin = null_lower, ymax = null_upper), fill = "grey85") +
    geom_line(aes(y = null_median, colour = "Shuffle median"),
              linewidth = 0.45, linetype = "dashed") +
    geom_line(aes(y = neff, colour = "Observed"), linewidth = 0.55) +
    geom_point(aes(y = neff, colour = "Observed"), size = 1.3) +
    scale_colour_manual(values = c("Observed" = unname(nature_colours["blue"]),
                                   "Shuffle median" = "grey35")) +
    scale_x_reverse(limits = c(260, 0), breaks = seq(250, 0, -50)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Age (Ma)", y = expression(N[eff]^"+"), colour = NULL,
         title = title, tag = tag) +
    theme_nature() +
    theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
          legend.key.width = unit(6, "mm"))
}
plot_correlation <- function(null_type, tag) {
  values <- null[null$null_type == null_type & null$statistic == "raw", ]
  result <- raw_summary[raw_summary$null_type == null_type, ]
  annotation <- paste0("Observed r = ", sprintf("%.3f", result$observed_r),
                       "\nPermutation P = ",
                       format(result$permutation_p_positive, scientific = FALSE))
  ggplot(values, aes(x = r)) +
    geom_histogram(breaks = seq(-1, 1, 0.04), fill = "grey60", colour = "white",
                   linewidth = 0.12) +
    geom_vline(xintercept = result$observed_r, colour = nature_colours["blue"],
                linewidth = 0.65) +
    annotate("text", x = -0.95, y = Inf, label = annotation,
              hjust = 0, vjust = 1.35, family = nature_font, size = nature_text_mm) +
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(x = expression("Pearson " * italic(r)), y = "Number of shuffles", tag = tag) +
    theme_nature()
}
p1 <- plot_trajectory("unrestricted", "Unrestricted shuffle", "a")
p2 <- plot_trajectory("length_quartile", "Within length quartiles", "b")
p3 <- plot_correlation("unrestricted", "c")
p4 <- plot_correlation("length_quartile", "d")

# Measure the grid on the same font-aware device used for export.
measurement <- tempfile(fileext = ".pdf")
cairo_pdf(measurement, family = nature_font, width = 183 / 25.4, height = 165 / 25.4)
print(p1)
combined <- arrangeGrob(p1, p2, p3, p4, ncol = 2, heights = c(1.15, 1))
dev.off()
unlink(measurement)
save_nature(combined, file.path(figure_dir, "branch-shuffle-test.pdf"),
             width_mm = 183, height_mm = 165)

# Empirical context: three separate axes avoid arbitrary double-axis scaling.
p_neff <- ggplot(observed, aes(age, neff)) +
  geom_line(colour = nature_colours["blue"], linewidth = 0.55) +
  geom_point(colour = nature_colours["blue"], size = 1.4) +
  scale_x_reverse(limits = c(260, 0), breaks = seq(250, 0, -50)) +
  labs(x = NULL, y = expression(N[eff]^"+"), tag = "a") + theme_nature()
p_div <- ggplot(observed, aes(age, pyrate)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_line(colour = nature_colours["orange"], linewidth = 0.55) +
  geom_point(colour = nature_colours["orange"], size = 1.4) +
  scale_x_reverse(limits = c(260, 0), breaks = seq(250, 0, -50)) +
  labs(x = NULL, y = "Net diversification\n(per Myr)", tag = "b") + theme_nature()
p_lineages <- ggplot(observed, aes(age, lineages)) +
  geom_line(colour = nature_colours["charcoal"], linewidth = 0.55) +
  geom_point(colour = nature_colours["charcoal"], size = 1.4) +
  scale_x_reverse(limits = c(260, 0), breaks = seq(250, 0, -50)) +
  labs(x = "Age (Ma)", y = "Branches per slice\n(bin mean)", tag = "c") + theme_nature()
cairo_pdf(measurement, family = nature_font, width = 110 / 25.4, height = 165 / 25.4)
print(p_div)
grobs <- list(ggplotGrob(p_neff), ggplotGrob(p_div), ggplotGrob(p_lineages))
widths <- unit.pmax(grobs[[1]]$widths, grobs[[2]]$widths, grobs[[3]]$widths)
for (i in 1:3) grobs[[i]]$widths <- widths
context <- arrangeGrob(grobs = grobs, ncol = 1)
dev.off()
unlink(measurement)
save_nature(context, file.path(figure_dir, "empirical-trajectories.pdf"),
             width_mm = 110, height_mm = 165)

captions <- c(
  "Branch-label randomization on the empirical tree.",
  paste0("The 77-tip tree, its branch lengths and 15 adult anatomical units were held fixed. ",
    "Accelerated-branch profiles were permuted independently within each unit, either ",
    "across all 152 branches (a, c) or within branch-length quartiles (b, d). Each profile ",
    "contains the same 1,000 reconstruction samples; their unit-specific counts are preserved. ",
    "There were ", format(analysis$n_shuffles, big.mark = ","), " permutations per null model. ",
    "Panels a and b show observed mean Neff+ (blue) and the median (dashed) and pointwise ",
    "central 95% interval (grey) of randomized mean trajectories in thirteen 20-Myr bins. ",
    "The vertical scales differ. Panels c and d show null Pearson correlations with the ",
    "fixed family-level PyRate mean net-diversification trajectory; blue lines mark the ",
    "observed correlation. Permutation P values are upper-tail probabilities for a positive ",
    "association, calculated as (1 + number of null correlations at least as large as ",
    "observed)/(1 + number of permutations)."),
  "", "Empirical temporal context.",
  paste0("Mean reconstructed Neff+ (a), family-level PyRate mean net diversification, ",
    "lambda minus mu (b), and the mean number of sampled-tree branches intersecting a ",
    "slice (c), in the same thirteen complete 20-Myr bins from 260 to 0 Ma. Neff+ is ",
    "averaged across 1,000 reconstruction samples. Branch counts describe the fixed ",
    "77-tip sampled phylogeny, not total historical lineage richness.")
)
writeLines(captions, file.path(output_dir, "figure-captions.txt"))

# A compact numerical report alongside the publication figures.
report <- c(
  "Empirical Neff+ branch-label randomization", "",
  paste0("77 tips; 152 branches; 15 adult units; 1,000 reconstruction samples; ",
         analysis$n_shuffles, " permutations per null model; seed ", analysis$seed, "."),
  "Thirteen complete 20-Myr bins, 260-0 Ma; fixed family PyRate mean net diversification.",
  sprintf("Observed Pearson r = %.6f; nominal two-sided P = %.6g.",
          raw_summary$observed_r[1], raw_summary$observed_nominal_p[1]),
  "")
for (i in 1:nrow(analysis$summary)) {
  result <- analysis$summary[i, ]
  report <- c(report,
    paste(result$null_type, "/", result$statistic),
    sprintf("  Observed r: %.6f; null median: %.6f; null central 95%%: [%.6f, %.6f].",
            result$observed_r, result$null_r_median, result$null_r_q025, result$null_r_q975),
    sprintf("  Null r >= observed: %d/%d; positive-tail permutation P: %.4f; absolute-r P: %.4f.",
            result$n_at_least_observed, result$n_shuffles,
            result$permutation_p_positive, result$permutation_p_absolute),
    sprintf("  Nominal Pearson P < 0.05: %.2f%% (either sign); %.2f%% (positive r only).",
            result$nominal_p_below_005_percent, result$positive_nominal_p_below_005_percent),
    "")
}
report <- c(report,
  sprintf("Sampled-tree branch count vs PyRate net diversification: r = %.6f.",
          cor(observed$lineages, observed$pyrate)),
  sprintf("Maximum difference from saved Neff curve: %.6g.", analysis$curve_error),
  sprintf("Maximum C++ vs plain-R difference: %.6g.", analysis$kernel_error),
  "",
  "Interpretation:",
  "Read the observed correlation against each empirical null distribution, not only against zero.",
  "Positive correlations can be common after preserving the association between branch length and acceleration.",
  "An unusually strong observed correlation supports temporal alignment beyond the specified null, not causality.",
  "These tests condition on the reconstructed branch-rate profiles and fixed PyRate mean curve.",
  "Branch-length quartiles approximate duration matching; they do not preserve every detection bias or clade structure.",
  "See README.md for the exact threshold convention and its existing zero-threshold edge case.")
writeLines(report, file.path(output_dir, "RESULTS.txt"))
