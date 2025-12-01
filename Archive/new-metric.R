
library(ape)
library(phytools)
library(ontophylo)
library(dplyr)
library(parallel)
library(viridis)
library(ggplot2)
library(deeptime)
library(mgcv)
source("R/utils-ST.R")
source('R/region_class.R')



tree <-readRDS("data/hym_tree.RDS")

# read mean rates
rates <-readRDS("data_out/br_rates_all_ind.RDS")
rates <- do.call(cbind,rates)
# drop larva
rates <- rates[,-16]
class(rates)
dim(rates)

#-------------


enr <- get_enrichment_mask(rates, lower_q = 0.05, upper_q = 0.95, direction = "over", return_type = "binary")
enr.branch <- apply(enr, 1, sum)
enr.boolean <- (enr.branch > 0)*1
hist(enr.branch)

plot_branch_colored_tree(tree, enr.boolean, palette = viridis(15), 
                         title = "all", legend_title = 'Cols')


out <- plot_enrichment_entropy(tree, enr, type = "Neff_scaled",n_points = 1000)
out


branch_entropy <- function(x) {
  # x = numeric or logical vector of 0/1 values (enriched / not-enriched)
  x <- as.numeric(x)
  k <- length(x)
  p1 <- sum(x) / k
  p0 <- 1 - p1
  
  # entropy with 0 * log(0) = 0
  safe_log <- function(z) ifelse(z == 0, 0, log(z))
  H <- - (p1 * safe_log(p1) + p0 * safe_log(p0))
  
  return(H)
}

branch_entropy(enr[1,])
apply(enr, 1, function(x) branch_entropy(x) )

branch_entropy(c(1,1,0))
