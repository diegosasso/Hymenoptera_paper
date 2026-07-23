
library(ape)
library(phytools)
library(ontophylo)
library(dplyr)
library(parallel)
library(viridis)
library(ggplot2)
library(abind)
source("R/utils-ST.R")
source('R/region_class.R')



tree <-readRDS("data/hym_tree.RDS")
max(nodeHeights(tree))

rates <-readRDS("data_out/br_rates_all_ind.RDS")
rates <- do.call(cbind,rates)
rates <- rates[,-16]
rates
rates
dim(rates)

enr <- get_enrichment_mask(rates, lower_q = 0.05, upper_q = 0.95, direction = "over", return_type = "binary")
neff <- plot_enrichment_entropy(tree, enr, type = "Neff_scaled", n_points = 1000 )

neff$plot
neff$data

plot(neff$data$time, neff$data$entropy, type='l', xlim = c(280,0),lwd=3,)
lines(neff$data$time, neff$data$H_plus_exp*neff$data$prop_enriched, col='red', lty=2)

plot(neff$data$time, neff$data$entropy, type='l', xlim = c(280,0))
lines(neff$data$time, neff$data$prop_enriched*2.6, col='red')
lines(neff$data$time, neff$data$H_plus_exp*.5, col='blue')

plot(neff$data$time, neff$data$H_plus, col='red', type='l', lwd=3, xlim = c(280,0))
lines(neff$data$time, neff$data$H_plus_exp*.34, col='blue')

plot(neff$data$time, neff$data$N_enriched, type = 'l', xlim = c(280,0))
lines(neff$data$time, neff$data$N_unique_enriched, col='blue', lwd=3, lty=2)



plot(neff$data$time, neff$data$H_plus, col='black', type='l', lwd=3, xlim = c(280,0))
lines(neff$data$time, neff$data$H_plus_exp*.34, col='blue')
lines(neff$data$time, neff$data$N_enriched*.2, col='green', lwd=3, lty=3)
lines(neff$data$time, neff$data$N_unique_enriched*.2, col='blue', lwd=3, lty=2)

# Contrubution
plot(neff$data$time, log(neff$data$entropy), type='l', xlim = c(280,0))
lines(neff$data$time, log(neff$data$prop_enriched), col='red')
lines(neff$data$time, neff$data$H_plus, col='blue')

