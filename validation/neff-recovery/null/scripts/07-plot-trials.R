# Full diagnostics for all trials; detailed unit plots for prespecified trial 001.
library(ape)
library(ggplot2)
library(gridExtra)
source("scripts/nature-style.R")
arguments <- commandArgs(trailingOnly = TRUE)
trial_ids <- if (length(arguments)) as.integer(arguments) else 1:100
figure_dir <- Sys.getenv("NEFF_FIGURE_DIR", "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
source_labels <- c(truth = "Simulation truth", reconstructed = "ER reconstruction")
source_colours <- c("Simulation truth" = "#444444", "ER reconstruction" = "#0072B2")
# A measurement device is kept open for Arial text layout between PNG exports.
cairo_pdf("tmp/measure-trials.pdf", width = 183/25.4, height = 235/25.4, family = nature_font)
measurement_device <- dev.cur()
cairo_pdf(file.path(figure_dir, "all-trial-diagnostics.pdf"), width = 183/25.4, height = 235/25.4,
          family = nature_font, onefile = TRUE)
atlas_device <- dev.cur()
for (trial in trial_ids) {
  sim_dir <- sprintf("simulated/trial-%03d", trial)
  out_dir <- sprintf("results/trial-%03d", trial)
  s <- readRDS(file.path(sim_dir, "simulation.rds"))
  metric <- readRDS(file.path(out_dir, "metrics.rds"))
  curves <- read.csv(file.path(out_dir, "neff-curves.csv"))
  corr <- read.csv(file.path(out_dir, "correlations.csv"))
  tree <- s$tree
  branch <- s$branch_design
  # Display-only protection against terminal ages such as -6e-14 Ma.
  stopifnot(min(branch$age_young) >= -1e-8, max(branch$age_old) <= 300 + 1e-8)
  n_edges <- nrow(tree$edge)
  node_y <- numeric(Ntip(tree) + tree$Nnode)
  tips_in_order <- tree$edge[tree$edge[, 2] <= Ntip(tree), 2]
  node_y[tips_in_order] <- seq_len(Ntip(tree))
  connectors <- data.frame()
  depth <- node.depth.edgelength(tree)
  for (node in rev(unique(tree$edge[, 1]))) {
    children <- tree$edge[tree$edge[, 1] == node, 2]
    node_y[node] <- mean(node_y[children])
    connectors <- rbind(connectors, data.frame(age = 300 - depth[node],
      ymin = min(node_y[children]), ymax = max(node_y[children])))
  }
  segments <- data.frame()
  for (edge in seq_len(n_edges)) {
    boundaries <- sort(unique(c(branch$age_old[edge], branch$age_young[edge],
      c(80, 180)[c(80, 180) > branch$age_young[edge] & c(80, 180) < branch$age_old[edge]])))
    for (j in seq_len(length(boundaries) - 1L)) {
      midpoint <- mean(boundaries[c(j, j+1)])
      lambda <- if (midpoint >= 80 && midpoint < 180) .028 else .010
      segments <- rbind(segments, data.frame(edge, old = boundaries[j+1],
        young = boundaries[j], y = node_y[tree$edge[edge, 2]], lambda))
    }
  }
  p_tree <- ggplot() +
    geom_segment(data = connectors, aes(x = age, xend = age, y = ymin, yend = ymax),
      colour = "grey70", linewidth = .2) +
    geom_segment(data = segments, aes(x = old, xend = young, y = y, yend = y,
      colour = factor(lambda)), linewidth = .35) +
    scale_colour_manual(values = c("0.01" = "#444444", "0.028" = "#D55E00"),
      name = expression(lambda~(Myr^{-1}))) + scale_x_reverse(limits = c(300, 0), oob = scales::squish) +
    labs(x = NULL, y = NULL, tag = "a", subtitle = sprintf("Root ER trial %03d; %d tips; no focal branches", trial, Ntip(tree))) +
    theme_nature() + theme(legend.position = "top", axis.text.y = element_blank(),
      axis.ticks.y = element_blank(), axis.line.y = element_blank())
  p_div <- ggplot(s$diversification, aes(age, rate)) +
    geom_line(colour = "#D55E00", linewidth = .6) + scale_x_reverse(limits = c(300, 0), oob = scales::squish) +
    labs(x = NULL, y = expression(lambda-mu~(Myr^{-1})), tag = "b",
      subtitle = "Known generating marginal diversification") + theme_nature()
  marg <- data.frame()
  for (source_name in intersect(names(source_labels), names(metric$rate_results))) {
    rates <- metric$rate_results[[source_name]]
    total <- rates$groups[[1]]$marginal_samples * 0
    for (unit in names(rates$groups)) total <- total + rates$groups[[unit]]$marginal_samples
    d <- data.frame(age = s$diversification$age, source = unname(source_labels[source_name]),
      mean = rowMeans(total), lower = 0, upper = 0)
    for (i in seq_len(nrow(total))) {
      bounds <- quantile(total[i, ], c(.025, .975), names = FALSE)
      d$lower[i] <- bounds[1]; d$upper[i] <- bounds[2]
    }
    marg <- rbind(marg, d)
  }
  p_marg <- ggplot(marg, aes(age, mean, colour = source, fill = source)) +
    geom_hline(yintercept = 200 * .003, linetype = 3, colour = "#009E73", linewidth = .5) +
    annotate("rect", xmin = 80, xmax = 180, ymin = -Inf, ymax = Inf, fill = "#E69F00", alpha = .1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .18, colour = NA) +
    geom_line(linewidth = .55) + scale_x_reverse(limits = c(300, 0), oob = scales::squish) +
    scale_colour_manual(values = source_colours) + scale_fill_manual(values = source_colours) +
    labs(x = NULL, y = expression(Changes~Myr^{-1}), tag = "c",
      subtitle = "Marginal character rate; dotted: generating 200q = 0.6", colour = NULL, fill = NULL) +
    theme_nature() + theme(legend.position = "none")
  x <- curves[curves$cutoff == "strict" & curves$source %in% names(source_labels), ]
  x$source <- factor(source_labels[x$source], levels = source_labels)
  cc <- corr[corr$cutoff == "strict" & corr$source == "reconstructed", ]
  text <- "Reconstruction unavailable; see trial diagnostics"
  if (nrow(cc)) text <- sprintf("Neff vs known diversification: r = %.3f, nominal p = %.3g", cc$Pearson_r, cc$nominal_p)
  p_neff <- ggplot(x, aes(age, neff, colour = source, fill = source)) +
    annotate("rect", xmin = 80, xmax = 180, ymin = -Inf, ymax = Inf, fill = "#E69F00", alpha = .1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .18, colour = NA) + geom_line(linewidth = .55) +
    scale_x_reverse(limits = c(300, 0), oob = scales::squish) + scale_colour_manual(values = source_colours) +
    scale_fill_manual(values = source_colours) + labs(x = "Age (Ma)", y = expression(N[eff]^"+"),
      colour = NULL, fill = NULL, tag = "d", subtitle = text) +
    theme_nature() + theme(legend.position = "bottom")
  dev.set(measurement_device)
  plots <- list(p_tree, p_div, p_marg, p_neff)
  grobs <- vector("list", 4)
  for (i in 1:4) { print(plots[[i]]); grobs[[i]] <- ggplotGrob(plots[[i]]) }
  common_widths <- grobs[[1]]$widths
  for (i in 2:4) common_widths <- grid::unit.pmax(common_widths, grobs[[i]]$widths)
  for (i in 1:4) grobs[[i]]$widths <- common_widths
  page <- arrangeGrob(grobs = grobs, ncol = 1, heights = c(3.3, 1.4, 1.7, 2.2))
  dev.set(atlas_device)
  grid::grid.newpage(); grid::grid.draw(page)
  if (trial == 1L) {
    ggsave(file.path(figure_dir, "trial-001-overview.pdf"), page, device = cairo_pdf, family = nature_font,
      width = 183, height = 235, units = "mm", bg = "white")
    ggsave(file.path(figure_dir, "trial-001-overview.png"), page, device = grDevices::png, type = "cairo",
      width = 183, height = 235, units = "mm", dpi = 300, bg = "white")
    # Save the coordinates for detailed unit-by-unit reconstruction plots.
    saveRDS(list(simulation = s, metrics = metric, node_y = node_y, connectors = connectors),
            "tmp/trial-001-plot-data.rds")
  }
}
dev.off(atlas_device)
dev.off(measurement_device)
unlink("tmp/measure-trials.pdf")

# Detailed acceleration and Neff alignment for trial 001, chosen before results.
a <- readRDS("tmp/trial-001-plot-data.rds")
s <- a$simulation
metric <- a$metrics
node_y <- a$node_y
connectors <- a$connectors
ages <- s$diversification$age
n_edges <- nrow(s$tree$edge)
unit_fractions <- branch_acceleration <- neff <- data.frame()
for (source_name in intersect(names(source_labels), names(metric$strict))) {
  strict <- metric$strict[[source_name]]
  for (u in seq_along(s$groups)) {
    unit <- names(s$groups)[u]
    mask <- matrix(strict$regimes$accelerated[, u, ], nrow = n_edges)
    branch_acceleration <- rbind(branch_acceleration, data.frame(source = unname(source_labels[source_name]),
      unit, old = s$branch_design$age_old, young = s$branch_design$age_young,
      y = node_y[s$tree$edge[, 2]], frequency = rowMeans(mask)))
    for (i in seq_along(ages)) {
      active <- which(s$branch_design$age_old >= ages[i] & s$branch_design$age_young < ages[i])
      fraction <- colMeans(mask[active, , drop = FALSE])
      bounds <- quantile(fraction, c(.025, .975), names = FALSE)
      unit_fractions <- rbind(unit_fractions, data.frame(source = unname(source_labels[source_name]),
        unit, unit_number = u, age = ages[i], mean = mean(fraction), lower = bounds[1], upper = bounds[2]))
    }
  }
  d <- strict$summary[strict$summary$metric == "Neff_plus", ]
  neff <- rbind(neff, data.frame(source = unname(source_labels[source_name]), age = ages,
    mean = d$mean, lower = d$lower, upper = d$upper))
}
write.csv(unit_fractions, "summary/trial-001-unit-fractions.csv", row.names = FALSE)
write.csv(branch_acceleration, "summary/trial-001-branch-acceleration.csv", row.names = FALSE)
cairo_pdf("tmp/measure-units.pdf", width = 183/25.4, height = 235/25.4, family = nature_font)
measurement_device <- dev.cur()
cairo_pdf(file.path(figure_dir, "trial-001-unit-acceleration.pdf"), width = 183/25.4, height = 235/25.4,
          family = nature_font, onefile = TRUE)
atlas_device <- dev.cur()
for (u in seq_along(s$groups)) {
  unit <- names(s$groups)[u]
  b <- branch_acceleration[branch_acceleration$unit == unit, ]; b$panel <- "Branches"
  f <- unit_fractions[unit_fractions$unit == unit, ]; f$panel <- "Unit fraction+"
  n <- neff; n$panel <- "Global Neff+"
  curve <- rbind(f[,c("source","age","mean","lower","upper","panel")], n)
  ctree <- connectors; ctree$panel <- "Branches"
  levels_panel <- c("Branches", "Unit fraction+", "Global Neff+")
  b$panel <- factor(b$panel, levels = levels_panel)
  curve$panel <- factor(curve$panel, levels = levels_panel)
  ctree$panel <- factor(ctree$panel, levels = levels_panel)
  b$source <- factor(b$source, levels = source_labels)
  curve$source <- factor(curve$source, levels = source_labels)
  p <- ggplot() +
    geom_vline(xintercept = c(80, 180), colour = "grey70", linetype = 3, linewidth = .3) +
    geom_segment(data = ctree, aes(x = age, xend = age, y = ymin, yend = ymax),
      colour = "grey80", linewidth = .15) +
    geom_segment(data = b, aes(x = old, xend = young, y = y, yend = y, colour = frequency), linewidth = .4) +
    geom_ribbon(data = curve, aes(age, ymin = lower, ymax = upper), fill = "#b6d5eb", alpha = .5) +
    geom_line(data = curve, aes(age, mean), colour = "#444444", linewidth = .55) +
    geom_line(data = curve[curve$source == "ER reconstruction", ], aes(age, mean), colour = "#0072B2", linewidth = .55) +
    scale_colour_gradient(low = "#c8c8c8", high = "#0072B2", limits = c(0, 1),
      breaks = c(0, .5, 1), name = "True acceleration / fraction of ASR maps") +
    guides(colour = guide_colourbar(title.position = "top", title.hjust = .5,
      barwidth = grid::unit(45, "mm"))) +
    scale_x_reverse(limits = c(300, 0), oob = scales::squish) + facet_grid(panel ~ source, scales = "free_y") +
    labs(x = "Age (Ma)", y = NULL, subtitle = paste("Root ER trial 001 —", unit),
      caption = "Acceleration: rate > the unit’s 95th percentile in each map. No focal branches.\nBands: central 95% conditional mapping intervals.\nNeff uses joint patterns across all ten units.") +
    theme_nature() + theme(legend.position = "top", strip.text = element_text(size = 9))
  dev.set(measurement_device); print(p)
  g <- ggplotGrob(p)
  panel_rows <- sort(unique(g$layout$t[grepl("^panel", g$layout$name)]))
  g$heights[panel_rows] <- grid::unit(c(3.5, 1, 1.3), "null")
  dev.set(atlas_device); grid::grid.newpage(); grid::grid.draw(g)
  if (u == 1L) {
    ggsave(file.path(figure_dir, "trial-001-unit-01.png"), g, device = grDevices::png, type = "cairo",
      width = 183, height = 235, units = "mm", dpi = 300, bg = "white")
  }
}
dev.off(atlas_device); dev.off(measurement_device)
unlink("tmp/measure-units.pdf")
