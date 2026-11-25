# Run from validation/phylo-sampling. Uses saved results; no inference is rerun.
# Script 13 integrates marginal curves uniformly over time within each epoch.
library(ggplot2)
source("scripts/figure-style.R")

# Each row summarizes 100 paired trials for one epoch and node proportion.
rates <- read.csv("summary/epoch-average-audit-summary.csv")
stopifnot(nrow(rates) == 15L, all(rates$n_trials == 100L),
          all(is.finite(rates$time_ratio_mean)),
          all(rates$time_ratio_lower >= 0),
          all(rates$time_ratio_upper <= 1.6))

epoch_names <- c("0–80 Ma", "80–160 Ma", "160–280 Ma")
rates$epoch_label <- factor(rates$epoch, levels = 1:3, labels = epoch_names)

# Show the five tested proportions as equally spaced categories.
rates$condition_tick <- factor(rates$condition, levels = rates$condition)
# Tick rates are full-tree marginal means, averaged over the 100 trials.
rates$axis_label <- sprintf("%.4f (%d%%)", rates$full_time_mean,
                            round(100 * rates$proportion))
tick_labels <- setNames(rates$axis_label, rates$condition)
rates$ratio_label <- sprintf("%.2f", rates$time_ratio_mean)
epoch_colours <- c("0–80 Ma" = "#0072B2", "80–160 Ma" = "#D55E00",
                   "160–280 Ma" = "#009E73")

# Separate epochs so their empirical 95% trial ranges cannot overlap.
plot <- ggplot(rates, aes(condition_tick, time_ratio_mean, colour = epoch_label)) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "#444444",
             linewidth = 0.45) +
  geom_errorbar(aes(ymin = time_ratio_lower, ymax = time_ratio_upper),
                width = 0.16, linewidth = 0.55) +
  geom_line(aes(group = 1), linewidth = 0.55) +
  geom_point(size = 2.1) +
  geom_text(aes(label = ratio_label), hjust = -0.25, vjust = -0.5,
             colour = "#222222", family = nature_font, size = nature_text_mm) +
  facet_wrap(~epoch_label, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = epoch_colours, guide = "none") +
  scale_x_discrete(labels = tick_labels, expand = expansion(add = c(0.55, 0.95))) +
  scale_y_continuous(limits = c(0, 1.6), breaks = seq(0, 1.5, 0.5),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(x = "Mean full-tree marginal rate (changes per Myr)\nNodes selected (%)",
       y = "Study / full marginal-rate ratio") +
  theme_nature() +
  theme(strip.text = element_text(size = 10, margin = margin(b = 7)),
        panel.spacing.x = grid::unit(3, "mm"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_text(margin = margin(t = 7)),
        axis.title.y = element_text(margin = margin(r = 7)))

save_nature(plot, "figures/epoch-rate-ratios.pdf", width_mm = 183, height_mm = 90)
