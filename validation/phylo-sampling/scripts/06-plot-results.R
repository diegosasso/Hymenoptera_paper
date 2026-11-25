# Publication plots use empirical ranges ACROSS trial means, not pooled map CIs.
library(ggplot2)
library(grid)
source("scripts/figure-style.R")
design <- readRDS("design/design.rds")
trials <- read.csv("summary/paired-trials.csv")
metrics <- read.csv("summary/metric-trials.csv")
curves <- read.csv("summary/marginal-summary.csv")
conditions <- design$conditions
figure_dir <- "figures"
is_pilot <- "--pilot" %in% commandArgs(trailingOnly = TRUE)
if (is_pilot) {
  figure_dir <- "figures/pilot"
} else {
  stopifnot(nrow(trials) == 1500L, length(unique(trials$task)) == 1500L,
            all(table(trials$condition) == 100L))
}
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
tree_colours <- c("Blaimer (765 tips)" = "#222222", "Study (77 tips)" = "#E69F00")
tree_lines <- c("Blaimer (765 tips)" = "solid", "Study (77 tips)" = "dashed")

# Repeat step boundaries so interval ribbons align with the stepped mean.
step_bands <- function(data) {
  answer <- data[FALSE, ]
  for (panel in unique(data$panel)) {
    for (tree in unique(data$tree)) {
      one <- data[data$panel == panel & data$tree == tree, , drop = FALSE]
      if (!nrow(one)) next
      one <- one[order(one$age, decreasing = TRUE), ]
      rows <- c(rep(seq_len(nrow(one)-1L), each=2), nrow(one))
      band <- one[rows, ]
      band$age <- c(one$age[1], rep(one$age[-1], each=2))
      answer <- rbind(answer, band)
    }
  }
  answer
}
for (comparison in c("all", "matched")) {
  groups <- c("full_all", "study_all")
  if (comparison == "matched") groups <- c("full_matched", "study_matched")
  filename <- file.path(figure_dir, paste0("marginal-rates-", comparison, ".pdf"))
  cairo_pdf(filename, width=183/25.4, height=230/25.4,
            family=nature_font, pointsize=nature_base_pt, onefile=TRUE)
  for (epoch in 1:3) {
    condition_ids <- conditions$condition[conditions$epoch == epoch]
    condition_ids <- intersect(condition_ids, unique(curves$condition))
    if (!length(condition_ids)) next
    data <- curves[curves$condition %in% condition_ids & curves$group %in% groups, ]
    data$tree <- ifelse(data$group == groups[1], "Blaimer (765 tips)", "Study (77 tips)")
    panel_labels <- character(length(condition_ids))
    for (i in seq_along(condition_ids)) {
      condition <- condition_ids[i]
      n <- unique(data$n_trials[data$condition == condition])
      stopifnot(length(n) == 1L)
      panel_labels[i] <- sprintf("%d%%: %d initial %s; %d %s",
        round(100*conditions$proportion[condition]),
        conditions$n_characters[condition],
        if (conditions$n_characters[condition] == 1L) "character" else "characters",
        n, if (n == 1L) "trial" else "trials")
    }
    data$panel <- factor(panel_labels[match(data$condition, condition_ids)],
                         levels=panel_labels)
    bands <- step_bands(data)
    # Publication explanations are in manuscript/figure-captions.tex.
    caption <- NULL
    if (is_pilot) caption <- "Screening only; across-trial ranges require the production run."
    plot <- ggplot(data, aes(age, mean, colour=tree, fill=tree)) +
      annotate("rect", xmin=design$epochs$lower_ma[epoch],
        xmax=design$epochs$upper_ma[epoch], ymin=-Inf, ymax=Inf,
        fill="#F2F2F2", colour=NA) +
      geom_ribbon(data=bands, aes(ymin=lower, ymax=upper), alpha=.18, colour=NA) +
      geom_step(aes(linetype=tree), linewidth=.55) +
      facet_wrap(~panel, ncol=1, scales="free_y") +
      scale_x_reverse(limits=c(280,0), breaks=c(280,240,160,80,0)) +
      scale_y_continuous(expand=expansion(mult=c(0,.08))) +
      scale_colour_manual(values=tree_colours) +
      scale_fill_manual(values=tree_colours) +
      scale_linetype_manual(values=tree_lines) +
      labs(x="Age (Ma)", y="Changes per Myr (combined character)",
        colour=NULL, fill=NULL, linetype=NULL,
        subtitle=sprintf("Origins at %d-%d Ma",
          design$epochs$lower_ma[epoch], design$epochs$upper_ma[epoch]),
        caption=caption) +
      theme_nature() + theme(legend.position="top",
        strip.text=element_text(hjust=0), panel.spacing=unit(2.5,"mm"))
    print(plot)
    ggsave(file.path(figure_dir, sprintf("marginal-rates-%s-epoch-%d.png",comparison,epoch)),
      plot, width=183, height=230, units="mm", dpi=300, bg="white")
  }
  dev.off()
}

# Use --marginal-only when revising just the marginal-rate figures.
if (!("--marginal-only" %in% commandArgs(trailingOnly = TRUE))) {
  # Four compact diagnostic panels, with identical epoch colours throughout.
  epoch_labels <- c("0-80 Ma","80-160 Ma","160-280 Ma")
  epoch_colours <- c("0-80 Ma"="#0072B2","80-160 Ma"="#D55E00","160-280 Ma"="#009E73")
  diagnostic <- data.frame()
  for (condition in conditions$condition) {
    one <- trials[trials$condition==condition,,drop=FALSE]
    if (!nrow(one)) next
    selected_metrics <- metrics[metrics$condition==condition & metrics$group=="study_all",]
    values <- list(
      "a  Characters remaining variable (%)"=100*one$retained_fraction,
      "b  Study / full epoch-rate ratio"=one$study_full_rate_ratio,
      "c  Study changes in target epoch (%)"=100*selected_metrics$mean[
        selected_metrics$metric=="target_epoch_fraction"],
      "d  Study changes on origin branch (%)"=100*selected_metrics$mean[
        selected_metrics$metric=="correct_origin_fraction"])
    for (label in names(values)) {
      x <- values[[label]]
      x <- x[is.finite(x)]
      if (!length(x)) next
      diagnostic <- rbind(diagnostic,data.frame(panel=label,condition=condition,
        epoch=epoch_labels[conditions$epoch[condition]],
        percent=100*conditions$proportion[condition],n_trials=length(x),
        mean=mean(x), lower=quantile(x,.025,names=FALSE),
        upper=quantile(x,.975,names=FALSE)))
    }
  }
  diagnostic$epoch <- factor(diagnostic$epoch,levels=epoch_labels)
  reference <- data.frame(panel=c("a  Characters remaining variable (%)",
    "b  Study / full epoch-rate ratio","c  Study changes in target epoch (%)",
    "d  Study changes on origin branch (%)"),value=c(100,1,100,100))
  recovery_caption <- paste0("Means and empirical central 95% ranges across trials.\n",
    "Temporal and branch localization use trials with surviving variable characters.")
  if (is_pilot) recovery_caption <- "Screening only; across-trial ranges require the production run."
  recovery <- ggplot(diagnostic,aes(percent,mean,colour=epoch,group=epoch)) +
    geom_hline(data=reference,aes(yintercept=value),linetype=3,colour="#999999",linewidth=.4) +
    geom_linerange(aes(ymin=lower,ymax=upper),linewidth=.45,alpha=.7) +
    geom_line(linewidth=.55) + geom_point(size=1.5) +
    facet_wrap(~panel,ncol=2,scales="free_y") +
    scale_colour_manual(values=epoch_colours) +
    scale_x_continuous(breaks=c(5,10,30,50,80)) +
    scale_y_continuous(expand=expansion(mult=c(.04,.08))) +
    labs(x="Nodes selected (%)",y=NULL,colour=NULL,
         caption=recovery_caption) +
    theme_nature() + theme(legend.position="top",
      strip.text=element_text(hjust=0),panel.spacing=unit(5,"mm"))
  ggsave(file.path(figure_dir,"recovery-summary.pdf"),recovery,device=cairo_pdf,
         width=183,height=155,units="mm",family=nature_font)
  ggsave(file.path(figure_dir,"recovery-summary.png"),recovery,
         width=183,height=155,units="mm",dpi=300,bg="white")
  write.csv(diagnostic,"summary/figure-diagnostics.csv",row.names=FALSE)

}
