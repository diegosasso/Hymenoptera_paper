# Compare the paired scenarios without changing either completed analysis.
library(ggplot2)
source("scripts/nature-style.R")
associated_correlations <- read.csv("../summary/correlations.csv")
null_correlations <- read.csv("summary/correlations.csv")
associated_correlations$scenario <- "Associated"
null_correlations$scenario <- "Root ER control"
correlations <- rbind(associated_correlations, null_correlations)
write.csv(correlations, "summary/comparison-correlations.csv", row.names = FALSE)

# Recompute every count directly from trial-level correlations.
percentages <- data.frame()
for (scenario in c("Associated", "Root ER control")) {
  reference <- if (scenario == "Associated") "../summary/p-value-percentages.csv" else "summary/p-value-percentages.csv"
  checked <- read.csv(reference)
  for (source in c("truth", "truth_retained", "reconstructed")) {
    for (cutoff in c("strict", "inclusive")) {
      d <- correlations[correlations$scenario == scenario & correlations$source == source &
                          correlations$cutoff == cutoff, ]
      finite <- is.finite(d$nominal_p) & is.finite(d$Pearson_r)
      for (alpha in c(.05, .01, .001)) {
        for (direction in c("two_sided_any_sign", "positive_r")) {
          detected <- finite & d$nominal_p < alpha
          if (direction == "positive_r") detected <- detected & d$Pearson_r > 0
          count <- sum(detected)
          interval <- binom.test(count, 100)$conf.int * 100
          original <- checked[checked$source == source & checked$cutoff == cutoff &
                                checked$alpha == alpha & checked$direction == direction, ]
          stopifnot(nrow(original) == 1L, original$n_detected == count,
                    original$n_estimable == sum(finite))
          percentages <- rbind(percentages, data.frame(scenario, source, cutoff, alpha,
            direction, n_planned = 100L, n_estimable = sum(finite),
            n_detected = count, percent = count, lower = interval[1], upper = interval[2]))
        }
      }
    }
  }
}
write.csv(percentages, "summary/comparison-percentages.csv", row.names = FALSE)

# A paired table records the scenario difference on each unchanged tree.
positive <- associated_correlations[associated_correlations$cutoff == "strict", ]
negative <- null_correlations[null_correlations$cutoff == "strict", ]
paired <- merge(positive, negative, by = c("trial", "source"),
  suffixes = c("_associated", "_null"))
paired$r_difference <- paired$Pearson_r_associated - paired$Pearson_r_null
write.csv(paired, "summary/paired-scenario-correlations.csv", row.names = FALSE)

source_labels <- c(truth = "Simulation truth", reconstructed = "ER reconstruction")
scenario_colours <- c("Associated" = "#0072B2", "Root ER control" = "#D55E00")
d <- percentages[percentages$cutoff == "strict" & percentages$source %in% names(source_labels), ]
d$source <- factor(source_labels[d$source], levels = source_labels)
d$scenario <- factor(d$scenario, levels = names(scenario_colours))
d$direction <- factor(d$direction, levels = c("two_sided_any_sign", "positive_r"),
  labels = c("Either sign", "Positive correlation"))
d$threshold <- factor(d$alpha, levels = c(.05, .01, .001), labels = c("0.05", "0.01", "0.001"))
d$label <- paste0(d$n_detected, "%")
p <- ggplot(d, aes(threshold, percent, colour = scenario, shape = scenario, group = scenario)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .12, linewidth = .4,
    position = position_dodge(width = .55)) +
  geom_point(size = 2, position = position_dodge(width = .55)) +
  geom_text(aes(y = upper + 4, label = label), family = nature_font, size = nature_text_mm,
    position = position_dodge(width = .55), vjust = 0, show.legend = FALSE) +
  facet_grid(direction ~ source) + scale_colour_manual(values = scenario_colours) +
  scale_y_continuous(limits = c(0, 112), breaks = seq(0, 100, 25)) +
  labs(x = "Nominal two-sided p-value threshold", y = "Trials below threshold (%)",
    colour = NULL, shape = NULL,
    caption = "100 trials per scenario, paired by tree. Bars: exact 95% binomial intervals.\nControl: root ER evolution independent of diversification; no temporal correction.") +
  theme_nature() + theme(legend.position = "top")
save_nature(p, "figures/scenario-p-value-percentages.pdf", 183, 170)

c <- correlations[correlations$cutoff == "strict" & correlations$source %in% names(source_labels), ]
c$source <- factor(source_labels[c$source], levels = source_labels)
c$scenario <- factor(c$scenario, levels = names(scenario_colours))
p <- ggplot(c, aes(nominal_p, fill = scenario)) +
  geom_histogram(breaks = seq(0, 1, .05), colour = "white", linewidth = .2) +
  geom_vline(xintercept = .05, linetype = 2, linewidth = .5, colour = "black") +
  facet_grid(scenario ~ source) + scale_fill_manual(values = scenario_colours) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
  labs(x = "Nominal two-sided Pearson p-value", y = "Trials",
    caption = "All trial-level correlations; strict acceleration cutoff. Dashed line: p = 0.05.") +
  theme_nature() + theme(legend.position = "none")
save_nature(p, "figures/scenario-p-value-distribution.pdf", 183, 170)

positive_trajectory <- read.csv("../summary/neff-trajectory-summary.csv")
negative_trajectory <- read.csv("summary/neff-trajectory-summary.csv")
positive_trajectory$scenario <- "Associated"
negative_trajectory$scenario <- "Root ER control"
trajectory <- rbind(positive_trajectory, negative_trajectory)
write.csv(trajectory, "summary/comparison-trajectories.csv", row.names = FALSE)
t <- trajectory[trajectory$cutoff == "strict" & trajectory$source %in% names(source_labels), ]
t$source <- factor(source_labels[t$source], levels = source_labels)
t$scenario <- factor(t$scenario, levels = names(scenario_colours))
source_colours <- c("Simulation truth" = "#444444", "ER reconstruction" = "#0072B2")
p <- ggplot(t, aes(age, mean, colour = source, fill = source)) +
  annotate("rect", xmin = 80, xmax = 180, ymin = -Inf, ymax = Inf, fill = "#E69F00", alpha = .1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), colour = NA, alpha = .13) +
  geom_line(linewidth = .6) + facet_wrap(~scenario, ncol = 1, scales = "free_y") +
  scale_x_reverse(limits = c(300, 0)) + scale_colour_manual(values = source_colours) +
  scale_fill_manual(values = source_colours) +
  labs(x = "Age (Ma)", y = expression(N[eff]^"+"), colour = NULL, fill = NULL,
    caption = "Lines: means across trials; bands: central 95% range of trial means.\nOrange interval: elevated generating diversification. Vertical scales differ.") +
  theme_nature() + theme(legend.position = "top")
save_nature(p, "figures/scenario-neff-through-time.pdf", 183, 170)

sources <- data.frame(quantity = c("Tree", "Diversification", "Terminal characters",
  "True histories and Neff", "Fitted Q", "ASR maps", "Reconstructed rates and Neff", "Pearson p-values"),
  source = c("Copied unchanged from the paired associated simulation",
    "Known generating lambda minus mu; not re-estimated",
    "New root-state-0 homogeneous ER simulation; q01=q10=0.003",
    "Exact new simulated histories, summarized by the same Neff formula",
    "Estimated from terminal states and tree, ER with maddfitz",
    "100 conditional histories per variable character; not simulation truth",
    "Computed from reconstructed maps; one metric per map before averaging",
    "Computed from 15 matched 20-Myr bins; nominal, two-sided"))
write.csv(sources, "summary/metric-sources.csv", row.names = FALSE)
print(d[, c("scenario", "source", "direction", "threshold", "n_detected", "n_estimable")])
