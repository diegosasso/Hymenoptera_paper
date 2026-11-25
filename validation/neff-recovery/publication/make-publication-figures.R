# Run from validation/neff-recovery, or open neff-recovery.Rproj in RStudio.
# source("publication/make-publication-figures.R")
# Only saved results are read. No simulation, fitting or mapping is repeated.
library(ape)
library(ggplot2)
library(gridExtra)
source("scripts/nature-style.R")

figure_dir <- "publication/figures"
data_dir <- "publication/source-data"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("publication/tmp", recursive = TRUE, showWarnings = FALSE)
blue <- "#0072B2"                 # Reconstructed histories, as in the existing figures.
orange <- "#D55E00"               # Elevated generating speciation.
grey <- "#666666"
scenario_labels <- c(Associated = "Correlated", Unassociated = "Independent")

# 1. Select all reconstructed trajectories under the primary strict cutoff ----
associated <- read.csv("summary/neff-curves.csv")
unassociated <- read.csv("null/summary/neff-curves.csv")
associated <- associated[associated$source == "reconstructed" & associated$cutoff == "strict", ]
unassociated <- unassociated[unassociated$source == "reconstructed" & unassociated$cutoff == "strict", ]
associated$scenario <- "Associated"
unassociated$scenario <- "Unassociated"
curves <- rbind(associated, unassociated)
curves$scenario <- factor(curves$scenario, levels = c("Associated", "Unassociated"))
stopifnot(nrow(associated) == 30000L, nrow(unassociated) == 30000L,
  identical(sort(unique(associated$trial)), 1:100),
  identical(sort(unique(unassociated$trial)), 1:100),
  all(is.finite(curves$neff)))
for (scenario in levels(curves$scenario)) {
  for (trial in 1:100) {
    d <- curves[curves$scenario == scenario & curves$trial == trial, ]
    stopifnot(nrow(d) == 300L, length(unique(d$age)) == 300L)
  }
}
means <- aggregate(neff ~ scenario + age, data = curves, FUN = mean)
write.csv(curves, file.path(data_dir, "reconstructed-neff-all-trials.csv"), row.names = FALSE)
write.csv(means, file.path(data_dir, "reconstructed-neff-means.csv"), row.names = FALSE)

# 2. One prespecified example tree, with exact epoch boundaries ----
simulation <- readRDS("simulated/trial-001/simulation.rds")
tree <- simulation$tree
branch <- simulation$branch_design
stopifnot(min(branch$age_young) >= -1e-8, max(branch$age_old) <= 300 + 1e-8)
node_y <- numeric(Ntip(tree) + tree$Nnode)
tip_order <- tree$edge[tree$edge[, 2] <= Ntip(tree), 2]
node_y[tip_order] <- seq_len(Ntip(tree))
depth <- node.depth.edgelength(tree)
connectors <- data.frame()
for (node in rev(unique(tree$edge[, 1]))) {
  children <- tree$edge[tree$edge[, 1] == node, 2]
  node_y[node] <- mean(node_y[children])
  connectors <- rbind(connectors, data.frame(age = 300 - depth[node],
    ymin = min(node_y[children]), ymax = max(node_y[children])))
}
segments <- data.frame()
for (edge in seq_len(nrow(tree$edge))) {
  breaks <- c(branch$age_old[edge], branch$age_young[edge])
  inside <- c(80, 180)[c(80, 180) > branch$age_young[edge] &
                         c(80, 180) < branch$age_old[edge]]
  breaks <- sort(unique(c(breaks, inside)))
  for (j in seq_len(length(breaks) - 1L)) {
    midpoint <- mean(breaks[c(j, j + 1L)])
    lambda <- if (midpoint >= 80 && midpoint < 180) .028 else .010
    segments <- rbind(segments, data.frame(edge, old = breaks[j + 1L], young = breaks[j],
      y = node_y[tree$edge[edge, 2]], lambda))
  }
}
write.csv(segments, file.path(data_dir, "example-tree-segments.csv"), row.names = FALSE)
lambda_steps <- data.frame(age = c(300, 180, 180, 80, 80, 0),
  lambda = c(.010, .010, .028, .028, .010, .010))
write.csv(lambda_steps, file.path(data_dir, "generating-lambda.csv"), row.names = FALSE)

# 3. Figure 1: common time axis and identical scales for the two Neff panels ----
compact_theme <- theme_nature() + theme(
  plot.margin = margin(1.5, 3, 1, 3, unit = "mm"),
  plot.subtitle = element_text(size = 9),
  axis.text.x = element_blank(), axis.ticks.x = element_blank(),
  axis.line.x = element_blank(), legend.position = "none")
age_scale <- scale_x_reverse(limits = c(300, 0), breaks = c(300, 200, 100, 0),
  expand = expansion(mult = .01), oob = scales::squish)
neff_scale <- scale_y_continuous(limits = c(0, 8.2), breaks = c(0, 4, 8),
  expand = expansion(mult = c(.015, .02)))

p_tree <- ggplot() +
  geom_segment(data = connectors, aes(x = age, xend = age, y = ymin, yend = ymax),
    colour = "grey70", linewidth = .15) +
  geom_segment(data = segments, aes(x = old, xend = young, y = y, yend = y,
    colour = factor(lambda)), linewidth = .20) +
  scale_colour_manual(values = c("0.01" = grey, "0.028" = orange)) + age_scale +
  labs(x = NULL, y = NULL, tag = "a") +
  compact_theme + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.line.y = element_blank())

p_lambda <- ggplot(lambda_steps, aes(age, lambda)) +
  annotate("rect", xmin = 80, xmax = 180, ymin = -Inf, ymax = Inf,
    fill = "#E69F00", alpha = .10) +
  geom_path(colour = grey, linewidth = .55) +
  annotate("segment", x = 180, xend = 80, y = .028, yend = .028,
    colour = orange, linewidth = .65) + age_scale +
  scale_y_continuous(limits = c(.008, .030), breaks = c(.010, .028),
    labels = c("0.010", "0.028"), expand = expansion(mult = 0)) +
  labs(x = NULL, y = expression(lambda~(Myr^{-1})), tag = "b") + compact_theme

neff_plots <- list()
for (scenario in c("Associated", "Unassociated")) {
  d <- curves[curves$scenario == scenario, ]
  average <- means[means$scenario == scenario, ]
  tag <- if (scenario == "Associated") "c" else "d"
  p <- ggplot(d, aes(age, neff)) +
    annotate("rect", xmin = 80, xmax = 180, ymin = -Inf, ymax = Inf,
      fill = "#E69F00", alpha = .10) +
    geom_line(aes(group = trial), colour = blue, alpha = .14, linewidth = .18) +
    geom_line(data = average, colour = blue, linewidth = .75) +
    age_scale + neff_scale +
    labs(x = NULL, y = expression(N[eff]^"+"), tag = tag,
      subtitle = unname(scenario_labels[scenario])) + compact_theme
  if (scenario == "Unassociated") {
    p <- p + labs(x = "Age (Ma)") + theme(axis.text.x = element_text(size = 9),
      axis.ticks.x = element_line(linewidth = nature_axis_mm),
      axis.line.x = element_line(linewidth = nature_axis_mm))
  }
  neff_plots[[scenario]] <- p
}

# Measure text on Cairo before aligning axes; the PDFs remain vector artwork.
cairo_pdf("publication/tmp/measure-figure-1.pdf", width = 183/25.4, height = 170/25.4,
  family = nature_font)
plots <- list(p_tree, p_lambda, neff_plots$Associated, neff_plots$Unassociated)
grobs <- vector("list", 4)
for (j in 1:4) {
  print(plots[[j]])
  grobs[[j]] <- ggplotGrob(plots[[j]])
}
common_widths <- grobs[[1]]$widths
for (j in 2:4) common_widths <- grid::unit.pmax(common_widths, grobs[[j]]$widths)
for (j in 1:4) grobs[[j]]$widths <- common_widths
figure_1 <- arrangeGrob(grobs = grobs, ncol = 1, heights = c(2.1, .9, 1.6, 1.9))
save_nature(figure_1, file.path(figure_dir, "figure-1-neff-and-diversification.pdf"),
  width_mm = 183, height_mm = 170, preview = FALSE)
dev.off()

# 4. Figure 2: every reconstructed correlation and its nominal p-value ----
correlations <- read.csv("null/summary/comparison-correlations.csv")
correlations <- correlations[correlations$source == "reconstructed" &
                               correlations$cutoff == "strict", ]
correlations$scenario[correlations$scenario == "Root ER control"] <- "Unassociated"
correlations$scenario <- factor(correlations$scenario, levels = c("Associated", "Unassociated"))
correlations <- correlations[order(correlations$scenario, correlations$trial), ]
stopifnot(nrow(correlations) == 200L, all(is.finite(correlations$Pearson_r)),
  all(is.finite(correlations$nominal_p)), all(correlations$n_bins == 15L))

# Horizontal offsets are for visibility only; r and p are never jittered.
set.seed(20261123)
offset <- runif(100, -.20, .20)
correlations$x <- as.integer(correlations$scenario) + offset[correlations$trial]
write.csv(correlations, file.path(data_dir, "reconstructed-correlations-and-p-values.csv"), row.names = FALSE)
counts <- data.frame()
for (scenario in levels(correlations$scenario)) {
  d <- correlations[correlations$scenario == scenario, ]
  stopifnot(nrow(d) == 100L, identical(sort(d$trial), 1:100))
  counts <- rbind(counts, data.frame(scenario, x = as.integer(d$scenario[1]),
    n = nrow(d), significant = sum(d$nominal_p < .05),
    median_r = median(d$Pearson_r), median_p = median(d$nominal_p)))
}
counts$label <- paste0(counts$significant, "/", counts$n, " below 0.05")
write.csv(counts, file.path(data_dir, "correlation-summary.csv"), row.names = FALSE)
scenario_axis <- scale_x_continuous(breaks = 1:2,
  labels = unname(scenario_labels[levels(correlations$scenario)]),
  limits = c(.5, 2.5), expand = expansion(mult = 0))
dot_theme <- theme_nature() + theme(plot.margin = margin(2, 4, 2, 3, unit = "mm"))

p_r <- ggplot(correlations, aes(x = x, y = Pearson_r)) +
  geom_hline(yintercept = 0, colour = "grey65", linetype = 3, linewidth = .35) +
  geom_boxplot(aes(x = as.integer(scenario), group = scenario), width = .48,
    outlier.shape = NA, colour = "#444444", fill = "#E1F1FC", linewidth = .35) +
  geom_point(colour = blue, size = 1.1, alpha = .55) +
  scenario_axis + scale_y_continuous(limits = c(-.10, .90), breaks = seq(0, .8, .2)) +
  labs(x = NULL, y = "Pearson correlation (r)", tag = "a") +
  dot_theme

p_p <- ggplot(correlations, aes(x = x, y = nominal_p)) +
  geom_hline(yintercept = .05, colour = orange, linetype = 2, linewidth = .5) +
  geom_boxplot(aes(x = as.integer(scenario), group = scenario), width = .48,
    outlier.shape = NA, colour = "#444444", fill = "#E1F1FC", linewidth = .35) +
  geom_point(colour = blue, size = 1.1, alpha = .55) +
  geom_text(data = counts, aes(x = x, y = .000055, label = label), inherit.aes = FALSE,
    family = nature_font, size = nature_text_mm, colour = "#222222") +
  scenario_axis +
  scale_y_log10(limits = c(.00004, 1.1), breaks = c(.0001, .001, .01, .05, .1, 1),
    labels = c("0.0001", "0.001", "0.01", "0.05", "0.1", "1")) +
  labs(x = NULL, y = "Nominal p-value (log scale)", tag = "b") + dot_theme

cairo_pdf("publication/tmp/measure-figure-2.pdf", width = 183/25.4, height = 112/25.4,
  family = nature_font)
print(p_r)
g_r <- ggplotGrob(p_r)
print(p_p)
g_p <- ggplotGrob(p_p)
common_heights <- grid::unit.pmax(g_r$heights, g_p$heights)
g_r$heights <- common_heights
g_p$heights <- common_heights
figure_2 <- arrangeGrob(g_r, g_p, ncol = 2)
save_nature(figure_2, file.path(figure_dir, "figure-2-correlations-and-p-values.pdf"),
  width_mm = 183, height_mm = 112, preview = FALSE)
dev.off()

# Render the authoritative PDFs so PNG previews retain exactly the same layout.
for (name in c("figure-1-neff-and-diversification", "figure-2-correlations-and-p-values")) {
  pdf <- file.path(figure_dir, paste0(name, ".pdf"))
  png_prefix <- file.path(figure_dir, name)
  status <- system2("pdftoppm", c("-f", "1", "-singlefile", "-r", "450", "-png", pdf, png_prefix))
  stopifnot(status == 0L)
}
capture.output(sessionInfo(), file = "publication/source-data/plotting-session.txt")
cat("Two publication figures saved; all 200 reconstructed trials retained.\n")
