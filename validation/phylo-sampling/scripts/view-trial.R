# Open phylo-sampling.Rproj in RStudio. Choose one completed trial below.
# This reads saved results only; it does not fit models or generate maps.
library(ggplot2)
source("scripts/figure-style.R")
task <- 1L
design <- readRDS("design/design.rds")
simulation <- readRDS(sprintf("simulated/task-%04d.rds", task))
result <- readRDS(sprintf("results/task-%04d/result.rds", task))
all_characters <- readRDS("design/node-characters.rds")
ids <- simulation$selected_nodes$character
characters_765 <- all_characters$full[, ids, drop=FALSE]
characters_77 <- all_characters$study[, ids, drop=FALSE]
print(simulation$selected_nodes)
print(c(simulated=result$n_characters, surviving=result$n_retained))
branches_765 <- result$groups$full_all$branch_summary
branches_77 <- result$groups$study_all$branch_summary

full <- result$groups$full_all$marginal_summary
study <- result$groups$study_all$marginal_summary
full$tree <- "Blaimer (765 tips)"
study$tree <- "Study (77 tips)"
all_label <- sprintf("All simulated: %d full, %d study",
                     result$n_characters, result$n_retained)
matched_label <- sprintf("Same surviving: %d characters on each tree", result$n_retained)
full$comparison <- study$comparison <- all_label
curves <- rbind(full,study)
if (!is.null(result$groups$full_matched)) {
  matched_full <- result$groups$full_matched$marginal_summary
  matched_full$tree <- "Blaimer (765 tips)"
  matched_full$comparison <- matched_label
  matched_study <- study
  matched_study$comparison <- matched_label
  curves <- rbind(curves,matched_full,matched_study)
}
curves$comparison <- factor(curves$comparison, levels=c(all_label, matched_label))
bands <- curves[FALSE,]
for(panel in unique(curves$comparison)) {
  for(tree in unique(curves$tree)) {
    one <- curves[curves$comparison==panel & curves$tree==tree,,drop=FALSE]
    one <- one[order(one$time,decreasing=TRUE),]
    rows <- c(rep(seq_len(nrow(one)-1L),each=2),nrow(one))
    band <- one[rows,]
    band$time <- c(one$time[1],rep(one$time[-1],each=2))
    bands <- rbind(bands,band)
  }
}
colours <- c("Blaimer (765 tips)"="#222222","Study (77 tips)"="#E69F00")
line_types <- c("Blaimer (765 tips)"="solid","Study (77 tips)"="dashed")
plot <- ggplot(curves,aes(time,mean,colour=tree,fill=tree)) +
  annotate("rect",xmin=design$epochs$lower_ma[result$epoch],
    xmax=design$epochs$upper_ma[result$epoch],ymin=-Inf,ymax=Inf,
    fill="#F2F2F2",colour=NA) +
  geom_ribbon(data=bands,aes(ymin=lower,ymax=upper),alpha=.18,colour=NA) +
  geom_step(aes(linetype=tree),linewidth=.65) +
  facet_wrap(~comparison,ncol=1) +
  scale_x_reverse(limits=c(280,0),breaks=c(280,240,160,80,0)) +
  scale_colour_manual(values=colours) + scale_fill_manual(values=colours) +
  scale_linetype_manual(values=line_types) +
  labs(x="Age (Ma)",y="Changes per Myr (combined character)",colour=NULL,fill=NULL,linetype=NULL,
    caption="Bands: pointwise 95% conditional mapping intervals for this trial.") +
  theme_nature() + theme(legend.position="top")
print(plot)

