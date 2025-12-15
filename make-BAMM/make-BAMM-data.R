# Load packages.
library(tidyverse)
library(phytools)
library(BAMMtools)

# Import raw BAMM data (again).
load("make-BAMM/data/bamm_data.RDA")

########
#
# Make bamm_fun to get diversification rate for downstream analyses
#
########


plotRateThroughTime(shifts, intervals = NULL, smooth = FALSE, ylim = c(0, 0.5),
                    ratetype = "netdiv", avgCol = "green")

rtt <- getRateThroughTimeMatrix(shifts, nslices = 280)

dim(rtt$lambda)

lam <- apply(rtt$lambda, 2, mean)
mu <- apply(rtt$mu, 2, mean)
dd <- lam-mu
dd <- rev(dd)
ti <- rtt$times
names(ti) <- NULL
plot(ti, dd, xlim=c(280, 0), type='l')

bamm <- cbind(ti, dd)

#----
make_bamm_fun <- function(time, divers) {
  if (length(time) != length(divers)) {
    stop("time and divers must have the same length")
  }
  # splinefun gives a smooth curve and linearly extrapolates
  f <- splinefun(time, divers, method = "natural")
  return(f)
}

# make function
bamm_fun <- make_bamm_fun(ti, dd)

# test
plot(ti, dd, xlim=c(280, 0), type='l')
lines(seq(0, 280, 0.1), bamm_fun(seq(0, 280, 0.1)), col='red')

saveRDS(bamm_fun, 'data/bamm_fun.rds')

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
