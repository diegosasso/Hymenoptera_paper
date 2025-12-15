
# Load packages.
library(tidyverse)
library(phytools)
library(BAMMtools)



#-------------------------------------------------------------------------#

###############################
### DIVERSIFICATION FIGURES ###
###############################

#-------------------------------------------------------------------------#

# Load data.
load("make-BAMM/data/paramo_stm_final.RDA")
# NOTE: The RDA file is a copy of the original file from the analysis conditioned on data (folder: "1_empirical_analyses").
load("make-BAMM/data/bamm_data.RDA")
# NOTE: in the original Blaimer's et al. study, the diversification rates were estimated using BAMM.

## GET DIVERSIFICATION DATA ##

# Get data from best shift configuration.
shifts <- getBestShiftConfiguration(shifts, expectedNumberOfShifts = 1, threshold = 5)

# Get data for the sub-tree containing the sampling from the present study.
div_shifts <- subtreeBAMM(shifts, tips = hym_data$Blaimer_taxa)

# Replace tip labels.
div_br$tip.label <- hym_data$Genus_Sharkey

# Get tip label colors matching ggplot2.
gg_cols <- scales::hue_pal()(4)

tip_lb <- rep(gg_cols[[3]], length(hym_tree$tip.label))
tip_lb[hym_data$Tag == "AcuEva"] <- gg_cols[[1]]
tip_lb[hym_data$Tag == "Ichneu"] <- gg_cols[[2]]
tip_lb[hym_data$Tag == "Procto"] <- gg_cols[[4]]

# Create folder to store figures.
suppressWarnings(dir.create("figures"))

#--------------------#

#-----------#
# PLOT TREE #
#-----------#

########
# FULL #
########

# NOTE: If desired, the code below can be used to plot the original Blaimer's et al. tree.

# Net diversification.
# Save a figure.
#png(paste0("figures/netdiv_baam_tree_full.png"), units = "in", width = 8, height = 8, bg = "white", res = 300)
#
#par(mfrow = c(1,1), mai = c(0,0,0,0))
#
#plot.bammdata(shifts, legend = TRUE, lwd = 2, method =  "phylogram", 
#              pal = "Spectral", breaksmethod = "jenks", spex = "netdiv")
#
#title(main = "BAMM - NetDiv (Blaimer et al. 2023)")
#
#dev.off()

#----------#

# Speciation.
# Save a figure.
#png(paste0("figures/speciation_baam_tree_full.png"), units = "in", width = 8, height = 8, bg = "white", res = 300)
#
#par(mfrow = c(1,1), mai = c(0,0,0,0))
#
#plot.bammdata(shifts, legend = TRUE, lwd = 2, method =  "phylogram", 
#              pal = "Spectral", breaksmethod = "jenks", spex = "s")
#
#title(main = "BAMM - Speciation (Blaimer et al. 2023)")
#
#dev.off()

#----------#

# Extinction.
# Save a figure.
#png(paste0("figures/extinction_baam_tree_full.png"), units = "in", width = 8, height = 8, bg = "white", res = 300)
#
#par(mfrow = c(1,1), mai = c(0,0,0,0))
#
#plot.bammdata(shifts, legend = TRUE, lwd = 2, method =  "phylogram", 
#              pal = "Spectral", breaksmethod = "jenks", spex = "e")
#
#title(main = "BAMM - Extinction (Blaimer et al. 2023)")
#
#dev.off()

#-----------#

#--------------------#

##########
# PRUNED #
##########

# Net diversification.
# Save a figure.
png(paste0("make-BAMM/figures/netdiv_baam_tree_pruned.png"), units = "in", width = 8, height = 8, bg = "white", res = 300)

par(mfrow = c(1,1), mai = c(0,0,0,0))

plot.bammdata(div_shifts, legend = TRUE, lwd = 5, method =  "phylogram", 
              pal = "Spectral", breaksmethod = "jenks", spex = "netdiv")

tiplabels(pch = 19, col = tip_lb, cex = 0.8)

title(main = "BAMM - NetDiv (Blaimer et al. 2023)")

dev.off()

#----------#

# Speciation.
# Save a figure.
png(paste0("make-BAMM/figures/speciation_baam_tree_pruned.png"), units = "in", width = 8, height = 8, bg = "white", res = 300)

par(mfrow = c(1,1), mai = c(0,0,0,0))

plot.bammdata(div_shifts, legend = TRUE, lwd = 5, method =  "phylogram", 
              pal = "Spectral", breaksmethod = "jenks", spex = "s")

tiplabels(pch = 19, col = tip_lb, cex = 0.8)

title(main = "BAMM - Speciation (Blaimer et al. 2023)")

dev.off()

#----------#

# Extinction.
# Save a figure.
png(paste0("make-BAMM/figures/extinction_baam_tree_pruned.png"), units = "in", width = 8, height = 8, bg = "white", res = 300)

par(mfrow = c(1,1), mai = c(0,0,0,0))

plot.bammdata(div_shifts, legend = TRUE, lwd = 5, method =  "phylogram", 
              pal = "Spectral", breaksmethod = "jenks", spex = "e")

tiplabels(pch = 19, col = tip_lb, cex = 0.8)

title(main = "BAMM - Extinction (Blaimer et al. 2023)")

dev.off()

#-----------#

#--------------------#

# NOTE: exported pngs used to make supplementary Figure S29 of the manuscript.

#----------------------------------------#

#-------------------------------#
# PLOT RATE-THROUGH-TIME CURVES #
#-------------------------------#

# Import raw BAMM data (again).
load("make-BAMM/data/bamm_data.RDA")

########
# FULL #
########

# NOTE: If desired, the code below can be used to plot the curves for the original Blaimer's et al. tree.

# PNG #
#png(paste0("curves_baam_tree_full.png"), units = "in", width = 16, height = 8, bg = "white", res = 300)
#
#plotRateThroughTime(shifts, intervals = NULL, smooth = FALSE, ylim = c(0, 0.5),
#                    ratetype = "netdiv", avgCol = "green")
#
#plotRateThroughTime(div_shifts, intervals = seq(from = 0, 0.9, by = 0.01), smooth = FALSE, ylim = c(0, 0.5),
#                    ratetype = "auto", intervalCol = "blue", avgCol = "blue", add = TRUE)
#
#plotRateThroughTime(div_shifts, intervals = seq(from = 0, 0.9, by = 0.01), smooth = FALSE, 
#                    ratetype = "extinction", intervalCol = "red", avgCol = "red", add = TRUE)
#
#abline(v = c(201, 237), col = "black", lty = 2, lwd = 2)
#
#dev.off()

#--------------------#

##########
# PRUNED #
##########

# Get data for the sub-tree containing the sampling from the present study.
div_shifts <- subtreeBAMM(shifts, tips = hym_data$Blaimer_taxa)



# PNG #
png(paste0("make-BAMM/curves_baam_tree_pruned.png"), units = "in", width = 16, height = 8, bg = "white", res = 300)

plotRateThroughTime(div_shifts, intervals = NULL, smooth = FALSE, ylim = c(0, 0.5),
                    ratetype = "netdiv", avgCol = "green")

plotRateThroughTime(div_shifts, intervals = seq(from = 0, 0.9, by = 0.01), smooth = FALSE, ylim = c(0, 0.5),
                    ratetype = "auto", intervalCol = "blue", avgCol = "blue", add = TRUE)

plotRateThroughTime(div_shifts, intervals = seq(from = 0, 0.9, by = 0.01), smooth = FALSE, 
                    ratetype = "extinction", intervalCol = "red", avgCol = "red", add = TRUE)

abline(v = c(201, 237), col = "black", lty = 2, lwd = 2)

dev.off()

#--------------------#

# NOTE: exported png used to make main Figure 1D of the manuscript.

#-------------------------------------------------------------------------#

###################################
### EXPLANATORY DIAGRAM FIGURES ###
###################################

#-------------------------------------------------------------------------#

# SET SIMULATION FUNCTIONS #

###--------------------------------------------------###

# BM function.
simu_bm <- function(n, sigma, x0) {
  
  # Simulate n values from a normal distribution.
  x <- rnorm(n = n, sd = sqrt(sig2) )
  
  # Get vector of values.
  x <- cumsum(c(x0, x))
  
  # Return results.
  return(x)
  
}

###--------------------------------------------------###

# OU function.
simu_ou <- function(n, alpha, theta, sig2, z0) {
  
  # Set initial vetors.
  dt <- 1
  x <- numeric(n)
  y <- numeric(n)
  x[1] <- z0
  y[1] <- alpha * (theta - z0) * dt
  
  for (i in 2:n) {
    
    # Simulate n - 1 values.
    x[i] <- x[i - 1] + alpha * (theta - x[i - 1]) * dt + sig2 * sqrt(dt) * rnorm(1)
    y[i] <- alpha * (theta - x[i - 1]) * dt
    
  }
  
  # Return results.
  return( list(x, y) )
  
}

###--------------------------------------------------###

# SIMULATIONS #

###--------------------------------------------------###

##################
### BM PROCESS ###
##################

# Set seed.
set.seed(42)

###----------###

# Set some parameters.
# Time bins.
n <- 100
# Diffusion rate.
sig2 <- 0.1
# Initial value.
x0 <- 0
# Number of simulations.
n_sim = 100

###----------###

# Simulate the BM process.
bm_sim <- numeric()
for (i in 1:n_sim) { bm_sim <- c(bm_sim, simu_bm(n, sig2, x0)) }

###--------------------------------------------------###

##################
### OU PROCESS ###
##################

# Set seed.
set.seed(42)

###----------###

# Set some parameters.
# Time bins.
n <- 100
# Diffusion rate.
sig2 <- 0.1
# Initial value.
z0 <- 0
# Number of simulations.
n_sim = 100
# Adaptation rate.
alpha <- 0.1
# Optimum value.
theta <- 5

###----------###

# Set initial vectors.
x <- numeric()
y <- numeric()

# Simulate the OU process.
for (i in 1:n_sim) { 
  
  ou_sim <- simu_ou(n = n, alpha = alpha, theta = theta, sig2 = sig2, z0 = z0)
  
  x <- c(x, ou_sim[[1]])
  y <- c(y, ou_sim[[2]])
  
}

###--------------------------------------------------###

# MAKE PLOTS #

###--------------------------------------------------###

##################
### BM PROCESS ###
##################

# Create a tibble for plotting.
plot_bm <- tibble(time = rep(0:n, n_sim), x = bm_sim, group = paste0("G",rep(1:n_sim, each = (n + 1))) )

# Plot.
ggplot(plot_bm, aes(x = time, y = x, group = group)) + geom_line(color = "grey", linewidth = 1.0) +
  
  labs(title = "Brownian Motion (BM)", x = "Time", y = "Trait value") + 
  
  ylim(c(-10, 10)) + 
  
  theme_minimal()

# Save figure.
ggsave(paste0("figures/simulation_bm.png"), units = "in", width = 8, height = 8, bg = "white")
dev.off()

# NOTE: exported pngs used to make main Figure 3E of the manuscript.

#-------------------------------------------------------------------------#

##################
### OU PROCESS ###
##################

# Create a tibble for plotting.
plot_ou <- tibble(time = rep(1:n, n_sim), x = x, y = y, group = paste0("G",rep(1:n_sim, each = (n))) )
plot_ou <- plot_ou %>% mutate(regime = rep(c(rep("dir", 20), rep("stab", 80)), n_sim))

# Plot.
ggplot(plot_ou, aes(x = time, y = x, group = group, color = regime)) + geom_line(show.legend = F, linewidth = 1.0) + 
  
  #scale_color_gradient(low = "purple", high = "orange") + 
  
  scale_color_manual(values = set_names(c("orange", "purple"), c("dir", "stab") ) ) + 
  
  labs(title = "Ornstein-Uhlenbeck (OU)", x = "Time", y = "Trait value") + 
  
  ylim(c(0, 8)) + 
  
  theme_minimal()

# Save figure.
ggsave(paste0("figures/simulation_ou.png"), units = "in", width = 8, height = 8, bg = "white")
dev.off()

# NOTE: exported pngs used to make main Figure 3E of the manuscript.

#-------------------------------------------------------------------------#
