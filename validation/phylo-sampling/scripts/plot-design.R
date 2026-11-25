# Publication design illustration: predetermined first trial at 50% in each epoch.
library(ape)
source("scripts/figure-style.R")
design <- readRDS("design/design.rds")
tree <- design$full_tree
tasks <- c(301L, 801L, 1301L)
dir.create("figures", showWarnings = FALSE)
draw_design <- function() {
  par(mfrow = c(3, 1), mar = c(2.5, 0.5, 2.4, 0.5),
      oma = c(0.7, 0, 0.6, 0), family = nature_font, cex = 0.9)
  for (panel in 1:3) {
    trial <- readRDS(sprintf("simulated/task-%04d.rds", tasks[panel]))
    selected <- trial$selected_nodes
    branch_colours <- rep("#E8E8E8", nrow(tree$edge))
    branch_colours[design$full_apo_edges] <- "#C7C7C7"
    plot(tree, show.tip.label = FALSE, edge.color = branch_colours, edge.width = 0.3,
         no.margin = FALSE, x.lim = c(0, design$root_age * 1.02))
    positions <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    node_colours <- rep("#666666", nrow(selected))
    node_colours[selected$variable_study] <- "#E69F00"
    points(positions$xx[selected$node], positions$yy[selected$node],
           pch = 16, col = node_colours, cex = 0.6)
    ages <- c(0, 80, 160, 280)
    axis(1, at = design$root_age - ages, labels = ages, cex.axis = 1)
    mtext("Age (Ma)", side = 1, line = 1.5, cex = 1)
    label <- sprintf("%s  %d-%d Ma; 50%% of nodes (%d characters)",
      letters[panel], design$epochs$lower_ma[panel], design$epochs$upper_ma[panel],
      nrow(selected))
    title(label, line = 0.5, cex.main = 1.05, font.main = 1, adj = 0)
  }
  par(mfrow = c(1, 1), fig = c(0, 1, 0, 1), new = TRUE,
      mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
  plot.new()
  legend("topright", inset = c(0.01, 0), bty = "n", horiz = TRUE, cex = 0.9,
         legend = c("Variable in study", "Lost from study"),
         pch = 16, col = c("#E69F00", "#666666"))
}
cairo_pdf("figures/simulation-design.pdf", width = 183/25.4, height = 170/25.4,
          family = nature_font, pointsize = nature_base_pt)
draw_design()
dev.off()
png("figures/simulation-design.png", width = 183, height = 170, units = "mm",
    res = 300, type = "cairo", pointsize = nature_base_pt)
draw_design()
dev.off()

