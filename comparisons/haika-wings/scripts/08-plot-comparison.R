# Run from comparisons/haika-wings after scripts/05-compare-rates.R completes.
# Before the first PDF-authoring run, call the bundled PDF artifact marker for
# four PDF outputs. Plotting does not fit models or alter rate normalization.
library(ape)
library(ggplot2)
library(grid)
library(gridExtra)

# Readable project style, adapted from the nature-figures skill template.
nature_font <- "Arial"
nature_base_pt <- 10
nature_tick_pt <- 9
nature_tag_pt <- 11
nature_text_mm <- nature_tick_pt / ggplot2::.pt
nature_curve_mm <- 0.75 / ggplot2::.pt
font_match <- systemfonts::match_fonts(nature_font)
if (!grepl("Arial", font_match$path[1], ignore.case = TRUE)) {
  stop("Arial is unavailable. Choose and document a project-wide replacement before rendering.")
}
theme_nature <- function() {
  theme_classic(base_size = nature_base_pt, base_family = nature_font) +
    theme(text = element_text(colour = "#222222"),
      axis.title = element_text(size = nature_base_pt),
      axis.text = element_text(size = nature_tick_pt, colour = "#222222"),
      axis.line = element_line(linewidth = 0.5 / ggplot2::.pt),
      axis.ticks = element_line(linewidth = 0.5 / ggplot2::.pt),
      axis.ticks.length = unit(1, "mm"),
      legend.title = element_blank(), legend.text = element_text(size = nature_tick_pt),
      legend.key.size = unit(3.5, "mm"), legend.background = element_blank(),
      strip.background = element_blank(), strip.text = element_text(size = nature_base_pt),
      plot.title = element_text(size = nature_base_pt),
      plot.tag = element_text(size = nature_tag_pt, face = "bold"),
      plot.margin = margin(3, 3, 3, 3, unit = "mm"),
      plot.background = element_rect(fill = "white", colour = NA))
}
save_nature <- function(plot, filename, width_mm, height_mm) {
  ggsave(filename, plot, device = cairo_pdf, family = nature_font,
         width = width_mm, height = height_mm, units = "mm", bg = "white")
  ggsave(sub("[.]pdf$", ".png", filename), plot, device = ragg::agg_png,
         width = width_mm, height = height_mm, units = "mm", dpi = 450, bg = "white")
}
plot_legend <- function(plot) {
  # Called while the surrounding Cairo measurement device is active. cowplot's
  # default legend extraction opens a separate PostScript device and loses Arial.
  plot_grob <- ggplotGrob(plot)
  legend_index <- which(plot_grob$layout$name == "guide-box-bottom")
  stopifnot(length(legend_index) == 1)
  plot_grob$grobs[[legend_index]]
}

# Only completed fits are eligible for the comparison figures.
branches <- read.csv("results/branch-comparison.csv", stringsAsFactors = FALSE)
marginal <- read.csv("results/marginal-comparison.csv", stringsAsFactors = FALSE)
associations <- read.csv("results/association-summary.csv", stringsAsFactors = FALSE)
fit_status <- read.csv("results/fit-status.csv", stringsAsFactors = FALSE)
stopifnot(nrow(branches) > 0, nrow(marginal) > 0, nrow(associations) > 0)
stopifnot(all(fit_status$comparison_complete), all(branches$dataset_lambda_complete),
          all(marginal$dataset_lambda_complete), all(associations$dataset_lambda_complete))
stopifnot(all(is.finite(branches$haika_relative)), all(branches$haika_relative >= 0),
          all(is.finite(branches$ontophylo_relative_mean)),
          all(branches$ontophylo_relative_mean >= 0))
dir.create("figures", showWarnings = FALSE)

primary_branches <- branches[branches$source == "seed3407" & branches$lambda == 1 &
  branches$axis == "TOTAL_VARIANCE" & branches$ontophylo_partition == "fore_wing" &
  branches$subset %in% c("full77", "family48_first"), ]
primary_curves <- marginal[marginal$source == "seed3407" & marginal$lambda == 1 &
  marginal$axis == "TOTAL_VARIANCE" & marginal$ontophylo_partition == "fore_wing" &
  marginal$subset %in% c("full77", "family48_first") & marginal$width_Myr >= 5 &
  marginal$exposure_lineage_Myr > 0, ]
stopifnot(sum(primary_branches$subset == "full77") == 152,
          sum(primary_branches$subset == "family48_first") == 94)
primary_branches$branch_class <- "Multiple families"
primary_branches$branch_class[primary_branches$family_stem_edge] <- "Family stem / tip"
primary_branches$branch_class[primary_branches$within_single_family] <- "Within family"
primary_branches$branch_class <- factor(primary_branches$branch_class,
  levels = c("Multiple families", "Family stem / tip", "Within family"))
branch_colours <- c("Multiple families" = "#0072B2", "Family stem / tip" = "#9B6EA7",
                    "Within family" = "#888888")
branch_shapes <- c("Multiple families" = 16, "Family stem / tip" = 17, "Within family" = 1)
method_colours <- c("Ontophylo forewing" = "#0072B2", "HAIKA family proxies" = "#D55E00")
method_lines <- c("Ontophylo forewing" = "solid", "HAIKA family proxies" = "22")
rate_limit <- max(primary_branches$haika_relative, primary_branches$ontophylo_relative_mean)
rate_breaks <- c(0, 1, 3, 10, 30, 100, 300, 1000)
rate_breaks <- rate_breaks[rate_breaks <= rate_limit]

# One point per branch. log(1 + rate) retains genuine zero estimates.
scatter77 <- ggplot(primary_branches[primary_branches$subset == "full77", ],
    aes(x = ontophylo_relative_mean, y = haika_relative)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, colour = "#BBBBBB") +
  geom_point(aes(colour = branch_class, shape = branch_class), size = 1.7, alpha = 0.8) +
  scale_colour_manual(values = branch_colours, drop = FALSE) +
  scale_shape_manual(values = branch_shapes, drop = FALSE) +
  scale_x_continuous(trans = "log1p", limits = c(0, rate_limit), breaks = rate_breaks,
                     labels = as.character(rate_breaks)) +
  scale_y_continuous(trans = "log1p", limits = c(0, rate_limit), breaks = rate_breaks,
                     labels = as.character(rate_breaks)) +
  labs(x = "Ontophylo relative rate", y = "HAIKA relative rate",
       title = "77-tip family alignment", tag = "a") +
  theme_nature() + theme(legend.position = "bottom")
scatter48 <- scatter77 + primary_branches[primary_branches$subset == "family48_first", ]
scatter48 <- scatter48 + labs(title = "One tip per family (48 tips)", tag = "b")

# Expand each bin to its two boundaries: lines/ribbons are piecewise constant,
# without linear interpolation between bin midpoints.
primary_curves <- primary_curves[order(primary_curves$subset, -primary_curves$older_age_Ma), ]
curve_steps <- primary_curves[rep(seq_len(nrow(primary_curves)), each = 2), ]
curve_steps$age_Ma <- as.vector(t(as.matrix(primary_curves[, c("older_age_Ma", "younger_age_Ma")])))
curve_max <- max(curve_steps$haika_relative, curve_steps$ontophylo_contracted_relative_upper,
                 na.rm = TRUE)
curve77 <- ggplot(curve_steps[curve_steps$subset == "full77", ], aes(x = age_Ma)) +
  geom_hline(yintercept = 1, colour = "#DDDDDD", linewidth = 0.25) +
  geom_ribbon(aes(ymin = ontophylo_contracted_relative_lower,
                  ymax = ontophylo_contracted_relative_upper),
              fill = method_colours["Ontophylo forewing"], alpha = 0.18) +
  geom_path(aes(y = ontophylo_contracted_relative_mean,
                colour = "Ontophylo forewing", linetype = "Ontophylo forewing"),
            linewidth = nature_curve_mm) +
  geom_path(aes(y = haika_relative, colour = "HAIKA family proxies",
                linetype = "HAIKA family proxies"), linewidth = nature_curve_mm) +
  scale_colour_manual(name = "Method", values = method_colours, breaks = names(method_colours)) +
  scale_linetype_manual(name = "Method", values = method_lines, breaks = names(method_colours)) +
  scale_x_reverse(breaks = c(250, 200, 150, 100, 50, 0), limits = c(280, 0), expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0, curve_max * 1.03), expand = c(0, 0)) +
  labs(x = "Time before present (Ma)", y = "Marginal relative rate",
       title = "77-tip family alignment", tag = "c") +
  theme_nature() + theme(legend.position = "bottom")
curve48 <- curve77 + curve_steps[curve_steps$subset == "family48_first", ]
curve48 <- curve48 + labs(title = "One tip per family (48 tips)", tag = "d")

# Measure and assemble on a Cairo device so Arial metrics are respected.
measurement_pdf <- tempfile(fileext = ".pdf")
cairo_pdf(measurement_pdf, width = 183 / 25.4, height = 170 / 25.4, family = nature_font)
branch_legend <- plot_legend(scatter77 + guides(colour = guide_legend(nrow = 1)))
method_legend <- plot_legend(curve77)
main_figure <- arrangeGrob(
  scatter77 + theme(legend.position = "none"), scatter48 + theme(legend.position = "none"),
  branch_legend, curve77 + theme(legend.position = "none"), curve48 + theme(legend.position = "none"),
  method_legend, layout_matrix = rbind(c(1, 2), c(3, 3), c(4, 5), c(6, 6)),
  heights = c(1, 0.13, 1, 0.13))
dev.off()
unlink(measurement_pdf)
save_nature(main_figure, "figures/haika-ontophylo-comparison.pdf", 183, 170)

# Sensitivity correlations: intervals reflect only conditional ontophylo map
# variation. They are not confidence intervals for HAIKA or proxy uncertainty.
association_keep <- associations$axis == "TOTAL_VARIANCE" &
  associations$ontophylo_partition == "fore_wing" &
  ((associations$domain == "branches" & associations$scope == "exclude_root_edges") |
   (associations$domain == "time_bins" & associations$scope == "time_full_width_bins"))
correlations <- associations[association_keep, ]
correlations$domain_label <- ifelse(correlations$domain == "branches", "Branches", "Time bins")
correlations$subset_label <- ifelse(correlations$subset == "full77", "77 tips",
  ifelse(correlations$subset == "family48_first", "48 tips, first", "48 tips, last"))

# Additional smoothing diagnostic, separate from the original 51-fit grid.
# Reuse the completed total-PC profiles and the collector's normalized maps.
# The point estimates use the mean normalized ontophylo rate; whiskers use
# correlations against individual maps, exactly as in the original panels.
bounds <- readRDS("results/bounds/fits.rds")
comparison_objects <- readRDS("results/comparison-objects.rds")
strong_summary <- bounds$summary[bounds$summary$scenario == "strong_smoothing_root_age", ]
stopifnot(bounds$complete, comparison_objects$complete, nrow(strong_summary) == 2L,
          all(strong_summary$success), all(strong_summary$n_successful_axes == 3L))
strong_penalty <- unique(strong_summary$lambda)
stopifnot(length(strong_penalty) == 1L)
strong_correlations <- correlations[FALSE, ]
strong_curves <- marginal[FALSE, ]
for (dataset_id in c("seed3407_full77", "seed3407_family48_first")) {
  profiles <- bounds$total_profiles[[dataset_id]]
  diagnostic <- profiles$profiles$strong_smoothing_root_age
  baseline <- profiles$profiles$primary_reference
  projection_id <- comparison_objects$dataset_projection_ids[[dataset_id]]
  projection <- comparison_objects$ontophylo_projections[[projection_id]]
  onto <- comparison_objects$ontophylo_summaries[[projection_id]]$partitions$fore_wing
  stopifnot(identical(diagnostic$tree$edge, projection$tree$edge),
            identical(diagnostic$tree$tip.label, projection$tree$tip.label))
  onto_time <- sweep(projection$marginal$contracted_branch_samples$fore_wing,
                     2, onto$normalization, "/")
  for (domain in c("branches", "time_bins")) {
    row <- correlations[correlations$dataset == dataset_id &
      correlations$lambda == 1 & correlations$domain == domain, ]
    stopifnot(nrow(row) == 1L)
    if (domain == "branches") {
      selected <- profiles$branch_subsets$exclude_root_edges
      haika_values <- diagnostic$relative_total_sigma2[selected]
      baseline_values <- baseline$relative_total_sigma2[selected]
      onto_samples <- onto$relative_samples[selected, , drop = FALSE]
      expected_rho <- strong_summary$forewing_rho_exclude_root_edges[
        strong_summary$dataset == dataset_id]
    } else {
      selected <- profiles$time_full_width_bins
      stopifnot(length(selected) == 28L)
      haika_values <- diagnostic$marginal_relative[selected]
      baseline_values <- baseline$marginal_relative[selected]
      onto_samples <- onto_time[selected, , drop = FALSE]
      expected_rho <- strong_summary$forewing_rho_time_full_width_bins[
        strong_summary$dataset == dataset_id]
    }
    onto_mean <- rowMeans(onto_samples)
    baseline_point <- cor(baseline_values, onto_mean, method = "spearman")
    baseline_draws <- as.numeric(cor(baseline_values, onto_samples, method = "spearman"))
    baseline_interval <- quantile(baseline_draws, c(0.025, 0.5, 0.975))
    stopifnot(abs(baseline_point - row$spearman_rho) < 1e-12,
      max(abs(baseline_interval - unlist(row[, c("map_rho_lower", "map_rho_median",
                                                "map_rho_upper")]))) < 1e-12)
    draw_correlations <- as.numeric(cor(haika_values, onto_samples, method = "spearman"))
    map_interval <- quantile(draw_correlations, c(0.025, 0.5, 0.975))
    row$lambda <- strong_penalty
    row$spearman_rho <- cor(haika_values, onto_mean, method = "spearman")
    row$map_rho_lower <- map_interval[1]
    row$map_rho_median <- map_interval[2]
    row$map_rho_upper <- map_interval[3]
    row$n_finite_map_correlations <- sum(is.finite(draw_correlations))
    stopifnot(abs(row$spearman_rho - expected_rho) < 1e-12)
    strong_correlations <- rbind(strong_correlations, row)
  }
  rows <- marginal[marginal$dataset == dataset_id & marginal$lambda == 1 &
    marginal$axis == "TOTAL_VARIANCE" & marginal$ontophylo_partition == "fore_wing", ]
  rows <- rows[order(rows$younger_age_Ma), ]
  stopifnot(nrow(rows) == length(diagnostic$marginal_relative),
    max(abs(rows$younger_age_Ma - head(diagnostic$breaks_Ma, -1))) < 1e-10,
    max(abs(rows$haika_relative - baseline$marginal_relative)) < 1e-10)
  rows$lambda <- strong_penalty
  rows$haika_raw <- diagnostic$marginal_raw
  rows$haika_relative <- diagnostic$marginal_relative
  strong_curves <- rbind(strong_curves, rows)
}
subset_colours <- c("77 tips" = "#0072B2", "48 tips, first" = "#D55E00", "48 tips, last" = "#009E73")
subset_shapes <- c("77 tips" = 16, "48 tips, first" = 17, "48 tips, last" = 15)
source_levels <- c("seed3407", "seed46567", "seed91507", "seed92779", "seed97081", "released28")
source_labels <- c("3407", "46567", "91507", "92779", "97081", "Released 28D")
seed_correlations <- correlations[correlations$lambda == 1 &
  correlations$subset %in% c("full77", "family48_first"), ]
seed_correlations$row_position <- 7 - match(seed_correlations$source, source_levels)
seed_correlations$row_position <- seed_correlations$row_position +
  ifelse(seed_correlations$subset == "full77", 0.12, -0.12)
seed_plot <- ggplot(seed_correlations, aes(y = row_position, colour = subset_label, shape = subset_label)) +
  geom_vline(xintercept = 0, colour = "#DDDDDD", linewidth = 0.3) +
  geom_segment(aes(x = map_rho_lower, xend = map_rho_upper, yend = row_position), linewidth = 0.4) +
  geom_point(aes(x = spearman_rho), size = 2) +
  scale_colour_manual(values = subset_colours) + scale_shape_manual(values = subset_shapes) +
  scale_x_continuous(limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1)) +
  scale_y_continuous(breaks = 6:1, labels = source_labels, limits = c(0.5, 6.5)) +
  facet_grid(. ~ domain_label) +
  labs(x = "Spearman correlation", y = "ANN initialization", tag = "a") +
  theme_nature() + theme(legend.position = "bottom")
penalty_correlations <- correlations[correlations$source == "seed3407", ]
penalty_correlations <- rbind(penalty_correlations, strong_correlations)
penalty_correlations$penalty_position <- log10(penalty_correlations$lambda)
penalty_correlations$penalty_position <- penalty_correlations$penalty_position +
  ifelse(penalty_correlations$subset == "full77", -0.08,
    ifelse(penalty_correlations$subset == "family48_first", 0, 0.08))
penalty_plot <- ggplot(penalty_correlations,
    aes(x = penalty_position, colour = subset_label, shape = subset_label, group = subset_label)) +
  geom_hline(yintercept = 0, colour = "#DDDDDD", linewidth = 0.3) +
  geom_segment(aes(xend = penalty_position, y = map_rho_lower, yend = map_rho_upper), linewidth = 0.4) +
  geom_line(aes(y = spearman_rho), linewidth = 0.3) +
  geom_point(aes(y = spearman_rho), size = 2) +
  scale_colour_manual(values = subset_colours) + scale_shape_manual(values = subset_shapes) +
  scale_x_continuous(breaks = log10(c(0.1, 1, 10, strong_penalty)),
                     labels = c("0.1", "1", "10", "280.03")) +
  scale_y_continuous(limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1)) +
  facet_grid(. ~ domain_label) +
  labs(x = "Rate-smoothing penalty", y = "Spearman correlation", tag = "b") +
  theme_nature() + theme(legend.position = "bottom")
measurement_pdf <- tempfile(fileext = ".pdf")
cairo_pdf(measurement_pdf, width = 183 / 25.4, height = 160 / 25.4, family = nature_font)
correlation_legend <- plot_legend(penalty_plot)
correlation_figure <- arrangeGrob(seed_plot + theme(legend.position = "none"),
  penalty_plot + theme(legend.position = "none"), correlation_legend,
  ncol = 1, heights = c(1.15, 1, 0.13))
dev.off()
unlink(measurement_pdf)
save_nature(correlation_figure, "figures/haika-ontophylo-sensitivity-correlations.pdf", 183, 160)

# Separate sensitivity curves show individual fitted trajectories, never a
# shaded interval assembled from seeds or smoothing choices.
sensitivity_curves <- marginal[marginal$axis == "TOTAL_VARIANCE" &
  marginal$ontophylo_partition == "fore_wing" & marginal$width_Myr >= 5 &
  marginal$exposure_lineage_Myr > 0 & marginal$subset %in% c("full77", "family48_first"), ]
sensitivity_curves <- rbind(sensitivity_curves,
  strong_curves[strong_curves$width_Myr >= 5 & strong_curves$exposure_lineage_Myr > 0, ])
sensitivity_curves <- sensitivity_curves[order(sensitivity_curves$source, sensitivity_curves$subset,
  sensitivity_curves$lambda, -sensitivity_curves$older_age_Ma), ]
sensitivity_steps <- sensitivity_curves[rep(seq_len(nrow(sensitivity_curves)), each = 2), ]
sensitivity_steps$age_Ma <- as.vector(t(as.matrix(sensitivity_curves[, c("older_age_Ma", "younger_age_Ma")])))
sensitivity_steps$subset_label <- factor(
  ifelse(sensitivity_steps$subset == "full77", "77 tips", "48 tips"),
  levels = c("77 tips", "48 tips"))
seed_steps <- sensitivity_steps[sensitivity_steps$lambda == 1, ]
seed_steps$source_label <- factor(seed_steps$source, levels = source_levels, labels = source_labels)
seed_colours <- c("3407" = "#0072B2", "46567" = "#D55E00", "91507" = "#009E73",
                  "92779" = "#CC79A7", "97081" = "#E69F00", "Released 28D" = "#444444")
seed_lines <- c("3407" = "solid", "46567" = "22", "91507" = "42",
                "92779" = "13", "97081" = "44", "Released 28D" = "longdash")
seed_curve_plot <- ggplot(seed_steps, aes(x = age_Ma, y = haika_relative,
    colour = source_label, linetype = source_label, group = source_label)) +
  geom_hline(yintercept = 1, colour = "#DDDDDD", linewidth = 0.25) +
  geom_path(linewidth = nature_curve_mm) +
  scale_colour_manual(values = seed_colours) + scale_linetype_manual(values = seed_lines) +
  scale_x_reverse(breaks = c(250, 200, 150, 100, 50, 0), limits = c(280, 0), expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.04))) +
  facet_grid(. ~ subset_label) +
  labs(x = "Time before present (Ma)", y = "HAIKA marginal relative rate", tag = "a") +
  theme_nature() + theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))
penalty_steps <- sensitivity_steps[sensitivity_steps$source == "seed3407", ]
penalty_steps$penalty_label <- factor(penalty_steps$lambda,
  levels = c(0.1, 1, 10, strong_penalty), labels = c("0.1", "1", "10", "280.03"))
penalty_curve_plot <- ggplot(penalty_steps, aes(x = age_Ma, y = haika_relative,
    colour = penalty_label, linetype = penalty_label, group = penalty_label)) +
  geom_hline(yintercept = 1, colour = "#DDDDDD", linewidth = 0.25) +
  geom_path(linewidth = nature_curve_mm) +
  scale_colour_manual(values = c("0.1" = "#009E73", "1" = "#0072B2",
                                 "10" = "#D55E00", "280.03" = "#CC79A7")) +
  scale_linetype_manual(values = c("0.1" = "22", "1" = "solid",
                                   "10" = "42", "280.03" = "longdash")) +
  scale_x_reverse(breaks = c(250, 200, 150, 100, 50, 0), limits = c(280, 0), expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.04))) +
  facet_grid(. ~ subset_label) +
  labs(x = "Time before present (Ma)", y = "HAIKA marginal relative rate", tag = "b") +
  theme_nature() + theme(legend.position = "bottom", legend.title = element_text(size = nature_tick_pt)) +
  guides(colour = guide_legend(title = "Penalty"), linetype = guide_legend(title = "Penalty"))
measurement_pdf <- tempfile(fileext = ".pdf")
cairo_pdf(measurement_pdf, width = 183 / 25.4, height = 170 / 25.4, family = nature_font)
curve_sensitivity_figure <- arrangeGrob(seed_curve_plot, penalty_curve_plot, ncol = 1)
dev.off()
unlink(measurement_pdf)
save_nature(curve_sensitivity_figure, "figures/haika-sensitivity-trajectories.pdf", 183, 170)

# Large paired trees retain 9 pt genus labels. These are inspection/supplementary
# figures, deliberately taller than the main-figure artwork profile.
tree <- readRDS("data/tree-77-morphology-labels.rds")
tree_rates <- primary_branches[primary_branches$subset == "full77", ]
tree_rates <- tree_rates[match(seq_len(nrow(tree$edge)), tree_rates$edge_id), ]
stopifnot(identical(as.integer(tree_rates$parent), tree$edge[, 1]),
          identical(as.integer(tree_rates$child), tree$edge[, 2]))
n_tips <- length(tree$tip.label)
root_node <- setdiff(tree$edge[, 1], tree$edge[, 2])[1]
stack <- root_node
tip_order <- integer(0)
while (length(stack) > 0) {
  node <- stack[1]
  stack <- stack[-1]
  if (node <= n_tips) {
    tip_order <- c(tip_order, node)
  } else {
    stack <- c(tree$edge[tree$edge[, 1] == node, 2], stack)
  }
}
node_y <- rep(NA_real_, n_tips + tree$Nnode)
node_y[tip_order] <- rev(seq_len(n_tips))
pending <- unique(tree$edge[, 1])
while (length(pending) > 0) {
  resolved <- integer(0)
  for (node in pending) {
    children <- tree$edge[tree$edge[, 1] == node, 2]
    if (all(is.finite(node_y[children]))) {
      node_y[node] <- mean(node_y[children])
      resolved <- c(resolved, node)
    }
  }
  stopifnot(length(resolved) > 0)
  pending <- setdiff(pending, resolved)
}
depths <- node.depth.edgelength(tree)
root_age <- max(depths[seq_len(n_tips)])
node_age <- root_age - depths
tree_segments <- data.frame(parent_age = node_age[tree$edge[, 1]],
  child_age = node_age[tree$edge[, 2]], y = node_y[tree$edge[, 2]],
  parent_y = node_y[tree$edge[, 1]],
  ontophylo = tree_rates$ontophylo_relative_mean, haika = tree_rates$haika_relative)
tree_labels <- data.frame(age = -7, y = node_y[seq_len(n_tips)], label = tree$tip.label)
tree_rate_max <- max(tree_segments$ontophylo, tree_segments$haika)
tree_onto <- ggplot(tree_segments) +
  geom_segment(aes(x = parent_age, xend = parent_age, y = parent_y, yend = y),
               colour = "#AAAAAA", linewidth = 0.3) +
  geom_segment(aes(x = parent_age, xend = child_age, y = y, yend = y, colour = ontophylo),
               linewidth = 0.75) +
  geom_text(data = tree_labels, aes(x = age, y = y, label = label),
            inherit.aes = FALSE, hjust = 0, family = nature_font,
            fontface = "italic", size = nature_text_mm) +
  scale_colour_viridis_c(option = "C", trans = "log1p", limits = c(0, tree_rate_max),
                         name = "Relative rate", breaks = rate_breaks) +
  scale_x_reverse(limits = c(root_age + 2, -105), breaks = c(250, 200, 150, 100, 50, 0),
                  expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, n_tips + 1), expand = c(0, 0)) +
  labs(x = "Time before present (Ma)", y = NULL, title = "Ontophylo forewing", tag = "a") +
  theme_nature() +
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom", legend.title = element_text(size = nature_tick_pt)) +
  guides(colour = guide_colourbar(barwidth = unit(70, "mm"), barheight = unit(3, "mm"),
                                  title.position = "top"))
tree_haika <- tree_onto
tree_haika$layers[[2]] <- geom_segment(aes(x = parent_age, xend = child_age,
  y = y, yend = y, colour = haika), linewidth = 0.75)
tree_haika <- tree_haika + labs(title = "HAIKA family proxies", tag = "b")
measurement_pdf <- tempfile(fileext = ".pdf")
cairo_pdf(measurement_pdf, width = 330 / 25.4, height = 330 / 25.4, family = nature_font)
tree_legend <- plot_legend(tree_onto)
paired_trees <- arrangeGrob(tree_onto + theme(legend.position = "none"),
  tree_haika + theme(legend.position = "none"), tree_legend,
  layout_matrix = rbind(c(1, 2), c(3, 3)), heights = c(1, 0.045))
dev.off()
unlink(measurement_pdf)
save_nature(paired_trees, "figures/haika-ontophylo-paired-rate-trees.pdf", 330, 330)

writeLines(c(
  "Comparison: 77-tip family-proxy alignment and one representative per family (48 tips).",
  "HAIKA is the sum of reconstructed variance rates across retained PCs, seed 3407, penalty 1.",
  "Ontophylo is the forewing partition. Each rate is normalized by its own lineage-duration-weighted treewide mean.",
  "Scatter axes use log(1+relative rate) scaling to retain zero estimates. The diagonal denotes equal relative rates.",
  "Within-family branches have parent and child descendants in one family; multiple-family branches subtend at least two families; remaining branches are family stems or terminals.",
  "Time curves are lineage-exposure-weighted means in 10 Myr bins; the oldest incomplete 0.0259 Myr bin is omitted.",
  "Ontophylo rates are projected onto the compared tree and made constant along contracted paths for matching temporal resolution.",
  "Blue ribbons show central 95 percent conditional ontophylo map intervals. HAIKA lines are point estimates, without confidence bands.",
  "Sensitivity correlations exclude the two root-adjacent edges for branch summaries and use full-width occupied temporal bins.",
  "The penalty panels additionally show a strong-smoothing diagnostic at lambda = 280.02588667, the original root age in Myr. This corresponds to approximately unit log-rate diffusion SD over the root height. It is separate from the original 51-fit grid and was fitted for the 77-tip and first 48-tip alignments only.",
  "Correlation points use the mean ontophylo rate. Whiskers summarize correlations for individual map draws and are not confidence intervals around the plotted point; a point can therefore lie outside its whiskers.",
  "These whiskers condition on each HAIKA point fit and do not quantify HAIKA or proxy uncertainty.",
  "Seed and penalty trajectories are separate fitted scenarios, not confidence intervals.",
  "ANN coordinates and newly computed PCs are not equated across seeds or identified as the publication's shape/venation axes.",
  "Paired trees use a shared relative-rate colour scale and identical tip ordering; text remains editable in PDF.",
  "The 77 tips contain 48 family profiles, and within-family variation is not measured by the family-proxy analysis."
), "figures/figure-captions.txt")
writeLines(capture.output(sessionInfo()), "figures/plotting-session.txt")
