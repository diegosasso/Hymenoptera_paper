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

## Journal-style PDF figure: net diversification rates vs time



pdf("figures/BAMM/BAMM_rates.pdf", width = 7.2, height = 3.5)

par(
  mar = c(4.2, 4.6, 1.0, 0.8),
  las = 1,
  bty = "l",
  cex.lab = 1.1,
  cex.axis = 0.9
)

plot(
  ti, dd,
  xlim = c(280, 0),
  type = "l",
  lwd = 2,
  col='green',
  xlab = "Time (Ma)",
  ylab = "Net diversification rate (BAMM)"
)

dev.off()



bamm <- cbind(time=ti, rate=dd)
bamm <- as_tibble(bamm)

saveRDS(bamm, 'data/bamm.rds')
