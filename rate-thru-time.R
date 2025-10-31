library(ape)
library(phytools)
library(ontophylo)
library(dplyr)
library(parallel)
library(viridis)
library(ggplot2)
library(abind)
library(dplyr)
library(purrr)
source("R/utils-ST.R")
source('R/region_class.R')

# read replicates
tree <-readRDS("data/hym_tree.RDS")
rates.multi <- readRDS("data_out/br_rates_all.RDS")

class(rates.multi)
length(rates.multi)
dim(rates.multi$cranium)

# names of body regions
regions <- names(rates.multi)
n_reps  <- ncol(rates.multi[[1]])   # 1000
n_branches <- nrow(rates.multi[[1]]) # 152

# build new list: one element per replicate
replicate_list <- lapply(1:n_reps, function(j) {
  # extract column j from each body region and bind
  mat <- sapply(rates.multi, function(region_mat) region_mat[, j])
  # ensure it's a matrix with rows = branches, cols = body regions
  mat <- matrix(mat, nrow = n_branches, ncol = length(regions))
  colnames(mat) <- regions
  mat
})

# drop larva
replicate_list <- lapply(replicate_list, function(x) x[, -16, drop = FALSE])

# check dims
length(replicate_list)   # 1000
dim(replicate_list[[1]]) # 152 x 15


tb <- replicate_list[[1]]
branch_rate <- apply(tb, 1, sum)

one_map_rate <- branch_rate_thru_time(tree, branch_rate, n_points = 100, time_points=NULL)
plot(one_map_rate, type='l')

branch_rate_thru_time_Multi <- function(tree, replicate_list, n_points = 100, time_points=NULL) {
  tb <- replicate_list[[1]]
  branch_rate <- apply(tb, 1, sum)
  one_map_rate <- branch_rate_thru_time(tree, branch_rate, n_points = n_points, time_points=time_points)
}

rates_time <- branch_rate_thru_time_Multi(tree, replicate_list, n_points = 1000, time_points = NULL)
plot(rates_time$time, rates_time$mean, type = 'l', xlim=c(280,0))

branch_rate_thru_time_Multi <- function(tree, replicate_list, n_points = 100, time_points = NULL) {
  
  # Compute rate-through-time for each replicate
  rate_list <- lapply(replicate_list, function(mat) {
    branch_rate <- rowMeans(mat)  # average across body regions per branch
    one_map_rate <- branch_rate_thru_time(tree, branch_rate,
                                          n_points = n_points,
                                          time_points = time_points)
    one_map_rate
  })
  
  # Combine all replicates
  all_rates <- dplyr::bind_rows(rate_list, .id = "replicate")
  
  # Summarise mean, median, and 95% CI (2.5%–97.5%)
  summary_rates <- all_rates %>%
    dplyr::group_by(time) %>%
    dplyr::summarise(
      mean   = mean(rate, na.rm = TRUE),
      median = median(rate, na.rm = TRUE),
      lower  = quantile(rate, 0.025, na.rm = TRUE),
      upper  = quantile(rate, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(summary_rates)
}


branch_rate_thru_time <- function(tree, branch_rate, n_points = 100, time_points=NULL) {

  # branch times (start, end for each edge)
  H <- nodeHeights(tree)
  max_height <- max(H)
  if (is.null(time_points)){
    times <- seq(0, max_height, length.out = n_points)
    times[n_points] <- times[n_points] - 1e-4
  } else { # time points provided by user
    times <- time_points
  }
  
  
  rate_vals <- numeric(length(times))
  # i=10
  for (i in seq_along(times)) {
    t <- times[i]
    edge_ids <- which(t >= H[,1] & t <= H[,2])  # edges active at time t
    if (length(edge_ids) > 0) {
      rates <- branch_rate[edge_ids]
      rate_vals[i] <- mean(rates)
      
    } else {
      rate_vals[i] <- NA
    }
  }
  
  df <- data.frame(time = rev(times), rate = rate_vals)
  
  return(df)
}
