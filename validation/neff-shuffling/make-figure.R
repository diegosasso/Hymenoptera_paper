# Recreate the publication figure from the supplied unrestricted shuffles.
# In RStudio, set the working directory to this folder (or the repository root).
# No simulation or ancestral-state reconstruction is run by this script.

library(ggplot2)
library(grid)
library(gridExtra)

analysis_dir <- "."
if (file.exists("validation/neff-shuffling/make-figure.R")) {
  analysis_dir <- "validation/neff-shuffling"
}
analysis_dir <- normalizePath(analysis_dir)
source_dir <- file.path(analysis_dir, "source-data")
figure_dir <- file.path(analysis_dir, "figures")
dir.create(figure_dir, showWarnings = FALSE)
source(file.path(analysis_dir, "scripts/nature-style.R"))
data <- readRDS(file.path(source_dir, "figure-data.rds"))
n_shuffles <- data$n_shuffles
stopifnot(n_shuffles == 9999, data$n_tips == 77, data$n_edges == 152,
          length(data$units) == 15, data$n_maps == 1000,
          ncol(data$shuffle_neff) == n_shuffles,
          identical(dim(data$shuffle_bins), c(13L, n_shuffles)))

# Recalculate every correlation and nominal p value from the saved bin means.
observed <- cor.test(data$observed_bins$neff, data$observed_bins$pyrate,
                     method = "pearson")
observed_r <- unname(observed$estimate)
observed_p <- observed$p.value
statistics <- data.frame(shuffle = 1:n_shuffles, r = NA_real_,
                         pearson_p_two_sided = NA_real_)
for (shuffle in 1:n_shuffles) {
  test <- cor.test(data$shuffle_bins[, shuffle], data$observed_bins$pyrate,
                   method = "pearson")
  statistics$r[shuffle] <- unname(test$estimate)
  statistics$pearson_p_two_sided[shuffle] <- test$p.value
}
statistics$minus_log10_p <- -log10(statistics$pearson_p_two_sided)
n_exceed <- sum(statistics$r >= observed_r)
permutation_p <- (n_exceed + 1) / (n_shuffles + 1)
n_significant <- sum(statistics$pearson_p_two_sided < 0.05)
n_positive_significant <- sum(statistics$r > 0 & statistics$pearson_p_two_sided < 0.05)
n_negative_significant <- sum(statistics$r < 0 & statistics$pearson_p_two_sided < 0.05)
summary <- data.frame(n_shuffles = n_shuffles, observed_r = observed_r,
                      observed_pearson_p_two_sided = observed_p,
                      median_null_r = median(statistics$r),
                      lower_null_r = quantile(statistics$r, 0.025, names = FALSE),
                      upper_null_r = quantile(statistics$r, 0.975, names = FALSE),
                      maximum_null_r = max(statistics$r),
                      n_at_least_observed = n_exceed,
                      permutation_p_positive = permutation_p,
                      n_nominal_p_below_005 = n_significant,
                      percent_nominal_p_below_005 = 100 * n_significant / n_shuffles,
                      n_positive_and_nominal_p_below_005 = n_positive_significant,
                      n_negative_and_nominal_p_below_005 = n_negative_significant)
write.csv(statistics, file.path(source_dir, "correlations.csv"), row.names = FALSE)
write.csv(summary, file.path(source_dir, "summary.csv"), row.names = FALSE)
write.csv(data$observed_bins, file.path(source_dir, "observed-bins.csv"), row.names = FALSE)

# Expand the group-wise means onto the full original time grid.
shuffle_mean <- rowMeans(data$shuffle_neff)
trajectory_summary <- data.frame(age_Ma = data$slice_age, slice_group = data$slice_group,
                                 observed_neff = data$observed_neff,
                                 shuffle_mean = shuffle_mean[data$slice_group])
write.csv(trajectory_summary, file.path(source_dir, "trajectory-summary.csv"), row.names = FALSE)

# Every shuffle is plotted. Remove only redundant vertices on flat segments:
# these 144 vertices draw exactly the same lines as all 1,000 original slices.
changes <- which(diff(data$slice_group) != 0)
vertices <- sort(unique(c(1L, length(data$slice_group), changes, changes + 1L)))
curve_values <- data$shuffle_neff[data$slice_group[vertices], , drop = FALSE]
curves <- data.frame(shuffle = rep(1:n_shuffles, each = length(vertices)),
                     age_Ma = rep(data$slice_age[vertices], n_shuffles),
                     neff = as.vector(curve_values))
stopifnot(length(unique(curves$shuffle)) == n_shuffles)
vertex_error <- 0
for (shuffle in c(1, 2, 100, 1000, n_shuffles)) {
  expanded <- approx(data$slice_age[vertices], curve_values[, shuffle],
                      xout = data$slice_age)$y
  original <- data$shuffle_neff[data$slice_group, shuffle]
  vertex_error <- max(vertex_error, max(abs(expanded - original)))
}
stopifnot(vertex_error < 1e-12)

# Panel a: all randomized trajectories, their arithmetic mean and the observed curve.
trajectory_colours <- c("All shuffles" = "#777777", "Shuffle mean" = "#E69F00",
                         "Observed" = "#009E73")
p_time <- ggplot() +
  geom_line(data = curves, aes(age_Ma, neff, group = shuffle, colour = "All shuffles"),
              linewidth = 0.12, alpha = 0.015) +
  geom_line(data = trajectory_summary,
              aes(age_Ma, shuffle_mean, colour = "Shuffle mean"), linewidth = 0.75) +
  geom_line(data = trajectory_summary,
              aes(age_Ma, observed_neff, colour = "Observed"), linewidth = 0.65) +
  scale_colour_manual(values = trajectory_colours,
                       breaks = c("All shuffles", "Shuffle mean", "Observed")) +
  guides(colour = guide_legend(override.aes = list(alpha = c(0.6, 1, 1),
                                                   linewidth = c(0.35, 0.75, 0.65)))) +
  scale_x_reverse(limits = c(max(data$slice_age), 0), breaks = seq(250, 0, -50),
                    expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(0, NA), breaks = seq(0, 20, 5),
                      expand = expansion(mult = c(0, 0.04))) +
  labs(x = "Age (Ma)", y = expression(N[eff]^"+"), colour = NULL, tag = "a") +
  theme_nature() +
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
          legend.key.width = unit(8, "mm"))

# Panel b: correlations, with the positive-tail empirical permutation p value.
r_annotation <- sprintf("Observed r = %.3f\nPermutation p = %.4f\n%d / %s as large as observed",
                          observed_r, permutation_p, n_exceed,
                          format(n_shuffles, big.mark = ","))
p_r <- ggplot(statistics, aes(r)) +
  geom_histogram(aes(y = after_stat(count / sum(count) * 100)),
                   breaks = seq(-1, 1, 0.04), fill = "grey60", colour = "white",
                   linewidth = 0.1) +
  geom_vline(xintercept = observed_r, colour = nature_colours["green"], linewidth = 0.65) +
  annotate("text", x = -0.98, y = Inf, label = r_annotation, hjust = 0, vjust = 1.25,
              family = nature_font, size = nature_text_mm) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.4))) +
  labs(x = expression("Pearson " * italic(r)), y = "Shuffles (%)", tag = "b") +
  theme_nature()

# Panel c: nominal two-sided p values, distinguished from the permutation p.
p_annotation <- sprintf("p < 0.05: %.2f%% (all r < 0)\nPositive and p < 0.05: %.0f%%",
                          100 * n_significant / n_shuffles,
                          100 * n_positive_significant / n_shuffles)
stopifnot(n_significant == n_negative_significant)
p_p <- ggplot(statistics, aes(minus_log10_p)) +
  geom_histogram(aes(y = after_stat(count / sum(count) * 100)),
                   breaks = seq(0, 6.6, 0.1), fill = "grey60", colour = "white",
                   linewidth = 0.1) +
  geom_vline(xintercept = -log10(0.05), colour = "grey30", linetype = "dashed",
                linewidth = 0.45) +
  geom_vline(xintercept = -log10(observed_p), colour = nature_colours["green"],
                linewidth = 0.65) +
  annotate("label", x = 0.1, y = Inf, label = p_annotation, hjust = 0, vjust = 1.4,
              family = nature_font, size = nature_text_mm, fill = "white",
              linewidth = 0, label.padding = unit(0.15, "lines"),
              label.r = unit(0, "lines")) +
  scale_x_continuous(limits = c(0, 6.6), breaks = 0:6) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.4))) +
  labs(x = expression(-log[10]("Pearson " * italic(p))), y = "Shuffles (%)", tag = "c") +
  theme_nature()

# Assemble with embedded Arial text at the final physical size.
measurement <- tempfile(fileext = ".pdf")
cairo_pdf(measurement, family = nature_font, width = 183 / 25.4, height = 170 / 25.4)
# A small plot initializes font metrics without drawing the trajectory cloud twice.
print(p_r)
figure <- arrangeGrob(p_time, p_r, p_p,
                       layout_matrix = rbind(c(1, 1), c(2, 3)), heights = c(1.1, 1))
dev.off()
unlink(measurement)
save_nature(figure, file.path(figure_dir, "neff-shuffling.pdf"),
             width_mm = 183, height_mm = 170)

write.csv(data.frame(check = c("plotted_shuffles", "time_grid_points", "drawn_vertices_per_curve",
                               "vertex_compression_error", "mean_uses_all_shuffles"),
                      value = c(n_shuffles, length(data$slice_age), length(vertices),
                                 vertex_error, ncol(data$shuffle_neff))),
          file.path(analysis_dir, "provenance/figure-validation.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(analysis_dir, "provenance/figure-session-info.txt"))
print(summary, digits = 6, row.names = FALSE)
