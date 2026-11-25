# Batch figures use only collected outputs. No inference occurs here.
library(ggplot2)
library(gridExtra)
source("scripts/nature-style.R")
curves <- read.csv("summary/neff-curves.csv")
trajectory <- read.csv("summary/neff-trajectory-summary.csv")
correlations <- read.csv("summary/correlations.csv")
percentages <- read.csv("summary/p-value-percentages.csv")
recovery <- read.csv("summary/recovery.csv")
trials <- read.csv("summary/trials.csv")
labels <- c(truth = "Simulation truth", reconstructed = "ER reconstruction")
colours <- c("Simulation truth" = "#444444", "ER reconstruction" = "#0072B2")
selected <- trajectory$cutoff == "strict" & trajectory$source %in% names(labels)
t <- trajectory[selected, ]
t$source <- factor(labels[t$source], levels = labels)
p_curves <- ggplot(t, aes(age, mean, colour = source, fill = source)) +
  annotate("rect", xmin = 80, xmax = 180, ymin = -Inf, ymax = Inf,
    fill = "#E69F00", alpha = .10) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .13, colour = NA) +
  geom_line(linewidth = .6) + scale_x_reverse(limits = c(300, 0)) +
  scale_colour_manual(values = colours) + scale_fill_manual(values = colours) +
  labs(x = "Age (Ma)", y = expression(N[eff]^"+"), colour = NULL, fill = NULL,
    caption = "Lines: means across trials. Bands: central 95% range of trial means.\nOrange interval: elevated generating diversification (180–80 Ma).") +
  theme_nature() + theme(legend.position = "top")
save_nature(p_curves, "figures/neff-through-time.pdf", 183, 125)

c <- correlations[correlations$cutoff == "strict" & correlations$source %in% names(labels), ]
c$source <- factor(labels[c$source], levels = labels)
p_limit <- max(.1, ceiling(max(c$nominal_p, na.rm = TRUE) * 10) / 10)
p_hist <- ggplot(c, aes(nominal_p, fill = source)) +
  geom_histogram(breaks = seq(0, p_limit, .01), colour = "white", linewidth = .2) +
  geom_vline(xintercept = .05, linetype = 2, colour = "#D55E00", linewidth = .5) +
  facet_wrap(~ source, ncol = 1) + scale_fill_manual(values = colours) +
  scale_x_continuous(limits = c(0, p_limit), breaks = seq(0, p_limit, .2)) +
  labs(x = "Nominal two-sided Pearson p-value", y = "Trials",
    caption = "One correlation per trial, using 15 matched 20-Myr bins. Dashed line: p = 0.05.\nThese tests do not adjust for temporal autocorrelation.") +
  theme_nature() + theme(legend.position = "none")
save_nature(p_hist, "figures/p-value-distribution.pdf", 183, 150)

d <- percentages[percentages$source %in% names(labels) & percentages$cutoff == "strict" &
                   percentages$direction %in% c("positive_r", "two_sided_any_sign"), ]
d$source <- factor(labels[d$source], levels = labels)
d$threshold <- factor(d$alpha, levels = c(.05, .01, .001),
                      labels = c("p < 0.05", "p < 0.01", "p < 0.001"))
d$direction <- factor(d$direction, levels = c("two_sided_any_sign", "positive_r"),
  labels = c("Either sign", "Positive correlation"))
d$label <- paste0(d$percent_all_trials, "%")
p_detect <- ggplot(d, aes(threshold, percent_all_trials, colour = source, group = source)) +
  geom_errorbar(aes(ymin = lower_all_trials, ymax = upper_all_trials),
    position = position_dodge(width = .65), width = .1, linewidth = .4) +
  geom_point(position = position_dodge(width = .65), size = 2) +
  geom_text(aes(y = upper_all_trials + 4, label = label), position = position_dodge(width = .65),
    vjust = 0, family = nature_font, size = nature_text_mm, show.legend = FALSE) +
  facet_wrap(~direction, ncol = 2) +
  scale_colour_manual(values = colours) + scale_y_continuous(limits = c(0, 110), breaks = seq(0,100,20)) +
  labs(x = NULL, y = "Trials below the p-value threshold (%)", colour = NULL,
    caption = "Denominator: all 100 prespecified trials. Bars: exact 95% binomial intervals.\nFalse association under homogeneous root ER evolution; nominal Pearson tests.") +
  theme_nature() + theme(legend.position = "top")
save_nature(p_detect, "figures/p-value-percentages.pdf", 183, 125)

true <- correlations[correlations$source == "truth" & correlations$cutoff == "strict", ]
est <- correlations[correlations$source == "reconstructed" & correlations$cutoff == "strict", ]
paired <- merge(true, est, by = "trial", suffixes = c("_truth", "_reconstructed"))
r <- recovery[recovery$reference == "truth", ]
p_correlation <- ggplot(paired, aes(Pearson_r_truth, Pearson_r_reconstructed)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey65", linetype = 2) +
  geom_point(size = 1.6, colour = "#0072B2", alpha = .7) + coord_equal() +
  labs(x = "True Neff–diversification r", y = "Reconstructed Neff–diversification r", tag = "a") + theme_nature()
p_agreement <- ggplot(r, aes(trajectory_r)) +
  geom_histogram(bins = 15, fill = "#0072B2", colour = "white", linewidth = .2) +
  labs(x = "Neff trajectory correlation", y = "Trials", tag = "b") + theme_nature()
p_rmse <- ggplot(r, aes(RMSE)) +
  geom_histogram(bins = 15, fill = "#0072B2", colour = "white", linewidth = .2) +
  labs(x = "Neff trajectory RMSE", y = "Trials", tag = "c") + theme_nature()
p_peaks <- ggplot(r, aes(true_peak_age, reconstructed_peak_age)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey65", linetype = 2) +
  geom_point(size = 1.6, colour = "#0072B2", alpha = .7) + coord_equal() +
  labs(x = "True peak age (Ma)", y = "Reconstructed peak age (Ma)", tag = "d") + theme_nature()
measure <- "tmp/measure-batch.pdf"
cairo_pdf(measure, width = 183/25.4, height = 170/25.4, family = nature_font)
plots <- list(p_correlation, p_agreement, p_rmse, p_peaks)
grobs <- vector("list", length(plots))
for (i in seq_along(plots)) { print(plots[[i]]); grobs[[i]] <- ggplotGrob(plots[[i]]) }
assembled <- arrangeGrob(grobs = grobs, ncol = 2)
save_nature(assembled, "figures/recovery-summary.pdf", 183, 170)
dev.off()
unlink(measure)

# Display every trial, ordered by its prespecified ID, not by its p-value.
all <- curves[curves$cutoff == "strict" & curves$source %in% names(labels), ]
all$source <- factor(labels[all$source], levels = labels)
cairo_pdf("figures/all-trial-neff-curves.pdf", width = 183/25.4, height = 235/25.4,
          family = nature_font, onefile = TRUE)
for (page in 0:9) {
  ids <- page*10 + 1:10
  x <- all[all$trial %in% ids, ]
  x$panel <- sprintf("Trial %03d", x$trial)
  p <- ggplot(x, aes(age, neff, colour = source)) +
    annotate("rect", xmin = 80, xmax = 180, ymin = -Inf, ymax = Inf,
      fill = "#E69F00", alpha = .10) + geom_line(linewidth = .45) +
    facet_wrap(~panel, ncol = 2, scales = "free_y") + scale_x_reverse(limits = c(300, 0)) +
    scale_colour_manual(values = colours) +
    labs(x = "Age (Ma)", y = expression(N[eff]^"+"), colour = NULL) +
    theme_nature() + theme(legend.position = "top")
  print(p)
}
dev.off()
