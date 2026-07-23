
# Load packages.
library(tidyverse)
library(phytools)
library(BAMMtools)

## ORGANIZE ORIGINAL BAMM DATA (Blaimer et al., 2023) ##

# Import original Blaimer's data.
hym_tree_full <- read.tree(file = "make-BAMM/data/HYMtopC_bamm_hisse.phy")

# Check tree.
plot.phylo(hym_tree_full, show.tip.label = F)
dev.off()

# Import event data.
e_data <- read.csv("make-BAMM/data/HYMtopC_bamm_event_data.txt")

# Get mismatching labels (left side).
wrong_lbs <- unique(e_data$leftchild[!(e_data$leftchild %in% hym_tree_full$tip.label)])

# Check correct labels.
hym_tree_full$tip.label[grep(hym_tree_full$tip.label, pattern = "Chalcimerus_borceai")]
hym_tree_full$tip.label[grep(hym_tree_full$tip.label, pattern = "Idarnes")]
hym_tree_full$tip.label[grep(hym_tree_full$tip.label, pattern = "Aperilampus_sp")]
hym_tree_full$tip.label[grep(hym_tree_full$tip.label, pattern = "Stilbula")]
hym_tree_full$tip.label[grep(hym_tree_full$tip.label, pattern = "Monacon_sp")]

# Replace wrong label.
e_data$leftchild[e_data$leftchild %in% wrong_lbs[[1]]] <- "Chalcimerus_borceai_JRAS01317_0301"
e_data$leftchild[e_data$leftchild %in% wrong_lbs[[2]]] <- "Idarnes_incertus_JRAS06107_0389"
e_data$leftchild[e_data$leftchild %in% wrong_lbs[[3]]] <- "Aperilampus_sp_JRAS07798_0101"
e_data$leftchild[e_data$leftchild %in% wrong_lbs[[4]]] <- "Stilbula_cyniformis_JRAS06279_0199"
e_data$leftchild[e_data$leftchild %in% wrong_lbs[[5]]] <- "Monacon_sp_JRAS07775_0189"

# Check again.
any(e_data$leftchild[!(e_data$leftchild %in% hym_tree_full$tip.label)])

# Get mismatching labels (right side).
wrong_lbs <- unique(e_data$rightchild[!(e_data$rightchild %in% hym_tree_full$tip.label)])
wrong_lbs <- wrong_lbs[-1]

# Check correct labels.
hym_tree_full$tip.label[grep(hym_tree_full$tip.label, pattern = "Aperilampus_sp")]
hym_tree_full$tip.label[grep(hym_tree_full$tip.label, pattern = "Idarnes")]
hym_tree_full$tip.label[grep(hym_tree_full$tip.label, pattern = "Chalcimerus_borceai")]
hym_tree_full$tip.label[grep(hym_tree_full$tip.label, pattern = "Stilbula")]
hym_tree_full$tip.label[grep(hym_tree_full$tip.label, pattern = "Monacon_sp")]

# Replace wrong label.
e_data$rightchild[e_data$rightchild %in% wrong_lbs[[1]]] <- "Aperilampus_sp_JRAS07798_0101"
e_data$rightchild[e_data$rightchild %in% wrong_lbs[[2]]] <- "Idarnes_incertus_JRAS06107_0389"
e_data$rightchild[e_data$rightchild %in% wrong_lbs[[3]]] <- "Chalcimerus_borceai_JRAS01317_0301"
e_data$rightchild[e_data$rightchild %in% wrong_lbs[[4]]] <- "Stilbula_cyniformis_JRAS06279_0199"
e_data$rightchild[e_data$rightchild %in% wrong_lbs[[5]]] <- "Monacon_sp_JRAS07775_0189"

# Check again.
any(!is.na(e_data$rightchild[!e_data$rightchild %in% hym_tree_full$tip.label]))

# Get number of shifts for the full tree.
shifts <- getEventData(hym_tree_full, eventdata = e_data, burnin = 0.1)

# Save data.
save(shifts, file = "make-BAMM/data/bamm_data.RDA")

#------------------------------#

## MAKE BAMM DATA FOR DOWNSTREAM ANALYSES ##

# Import raw BAMM data (again).
load("make-BAMM/data/bamm_data.RDA")

# Extract relevant data.
rtt <- getRateThroughTimeMatrix(shifts, nslices = 280)

dim(rtt$lambda)

lam <- apply(rtt$lambda, 2, mean)
mu <- apply(rtt$mu, 2, mean)
dd <- lam - mu
dd <- rev(dd)
ti <- rtt$times
names(ti) <- NULL

bamm <- cbind(time = ti, rate = dd)
bamm <- as_tibble(bamm)

# Save RDS object.
saveRDS(bamm, 'data/bamm.RDS')

#------------------------------#

## PLOT RATE-THROUGH-TIME CURVES ##

# Create folder to store images.
dir.create("figures/BAMM", recursive = TRUE)

# Load data.
load("data_out/paramo_stm_adult_final.RDA")
load("make-BAMM/data/bamm_data.RDA")

# Get data for the sub-tree containing the sampling from the present study.
div_shifts <- subtreeBAMM(shifts, tips = hym_data$Blaimer_taxa)

# Save PDF.
pdf("figures/BAMM/BAMM_rates.pdf", width = 7.2, height = 3.5)

par(mar = c(4.2, 4.6, 1.0, 0.8), las = 1, bty = "l", cex.lab = 1.1, cex.axis = 0.9)

plotRateThroughTime(div_shifts, intervals = NULL, smooth = FALSE, ylim = c(0, 0.5),
                    ratetype = "netdiv", avgCol = "green", axis.labels = FALSE)

plotRateThroughTime(div_shifts, intervals = seq(from = 0, 0.9, by = 0.01), smooth = FALSE, ylim = c(0, 0.5),
                    ratetype = "auto", intervalCol = "blue", avgCol = "blue", add = TRUE)

plotRateThroughTime(div_shifts, intervals = seq(from = 0, 0.9, by = 0.01), smooth = FALSE,
                    ratetype = "extinction", intervalCol = "red", avgCol = "red", add = TRUE)

title(xlab = "Time (Ma)", ylab = "Lineage rates (BAMM)")

dev.off()
