
# Load packages.
library(phytools)
library(ontophylo)
library(tidyverse)

# Load data.
load("data_out/paramo_stm_adult_final.RDA")

#------------------------------------------------------------#

### Define fome helper functions.

# Get time segments for each state.
get_time_seg <- function(stm) {
  
  # Get starting node heights.
  H <- nodeHeights(stm)[,1]
  
  # Get timings.
  br_sums <- lapply(stm$maps, function(x) cumsum(x) )
  times <- mapply(x = H, y = br_sums, function(x,y) unname(c(x, x+y)) )
  
  df <- vector(mode = "list", length = length(times))
  # Reorganize as data.frame.
  for (i in 1:length(df)) {
    
    df[[i]] <- data.frame("br_id" = i, "start" = head(times[[i]],-1), "end" = tail(times[[i]],-1), "state" = names(stm$maps[[i]]))
    
  }
  
  # Return results.
  return(do.call(rbind, df))
  
}

#----------#

frechet_hamming_mean_var <- function(strings, weights = NULL) {
  stopifnot(is.character(strings), length(strings) >= 1)
  
  # Check same length
  Ls <- nchar(strings)
  if (length(unique(Ls)) != 1)
    stop("All strings must have the same length.")
  L <- Ls[1]
  n <- length(strings)
  
  # Default weights
  if (is.null(weights)) weights <- rep(1, n)
  
  # Convert to character matrix
  mat <- do.call(rbind, strsplit(strings, ""))
  storage.mode(mat) <- "character"
  
  # Weighted consensus (majority rule)
  consensus <- character(L)
  for (j in seq_len(L)) {
    tab <- tapply(weights, mat[, j], sum)
    consensus[j] <- names(which.max(tab))
  }
  
  # Consensus string
  consensus_str <- paste(consensus, collapse = "")
  
  # Compute distances of each string from consensus
  dist_vec <- colSums(t(mat) != consensus)
  
  # Mean and variance
  #mean_dist <- mean(dist_vec)
  var_dist  <- var(dist_vec)
  
  var_dist
  
  # list(
  #   consensus     = consensus_str,
  #   mean_distance = mean_dist,
  #   var_distance  = var_dist
  # )
}


# Make disparity data.
# n_points = 1000
make_dispar_data <- function(char_list, n_points = 100, time_points = NULL) {

  # Get max time.
  tmax <- max(char_list[[1]]$end)
    
  # Make time points to sample from.
  if (is.null(time_points)){
    
    times <- seq(0, tmax, length.out = n_points)
    times[n_points] <- times[n_points] - 1e-4
    
  } else { # time points provided by user
    
    times <- time_points
    
  }
  
  # Initialize a vector to store results for all time points.
  dispar_ls <- vector(mode = "list", length = n_points)
  
  # Loop over all time points.
  # k=1
  for(k in 1:length(times)) {
    
    # Get time slice.
    t <- times[k]
    
    # Get all states/lineages represented in a time slice.
    states_time <- lapply(char_list, function(x) x[which(t >= x$start & t <= x$end),c(1,4)])
    
    # Amalgamate states.
    st <- apply(do.call(cbind,lapply(states_time, function(x) x$state )),1,paste0, collapse = "")
    
    # Calculate Hamming distances.
    # hm <- stringdist::stringdistmatrix(st, st, method = "hamming")
    # hm <- mean(hm[upper.tri(hm)])
    hm <- stringdist::stringdistmatrix(st, st, method = "hamming")
    hm <- var(hm[upper.tri(hm)])
    if (is.na(hm)){
      hm <- 0
    }
    
    # Frechet mean
    #hm <-frechet_hamming_mean(st)
    #hm <-frechet_hamming_mean_var(st)
    
    dispar_ls[[k]] <- tibble("Time" = t, "Disparity" = hm, "N_lineages" = length(st))
    
  }
  
  # Return results.
  return(do.call(bind_rows, dispar_ls) %>% mutate(Time = rev(Time)))
  
}



#------------------------------------------------------------#

library(parallel)

n_samples <- 100
CH_IDs    <- hym_annot_final$char_id

# Use all but one core
mc <- max(1, detectCores() - 1)


disp_data_all <- mclapply(1:n_samples, mc.cores = mc, FUN = function(i) {
  # cat(sprintf("\nWorking on map: %d %s\n", i, Sys.time()))
  
  char_list <- setNames(
    lapply(CH_IDs, function(ch) {
      stm <- readRDS(file.path("stmaps", paste0(ch, ".RDS")))
      get_time_seg(stm[[i]])
    }),
    CH_IDs
  )
  
  make_dispar_data(char_list, n_points = 1000)
})

# i=1
# char_list <- setNames(
#   lapply(CH_IDs[1:10], function(ch) {
#     stm <- readRDS(file.path("stmaps", paste0(ch, ".RDS")))
#     get_time_seg(stm[[i]])
#   }),
#   CH_IDs[1:10]
# )
# length(char_list)
#------------------------------------------------------------#
# # Set some parameters.
# n_samples = 1000
# 
# # Get character IDs.
# CH_IDs <- hym_annot_final$char_id
# 
# # Initialize list to store results from all samples.
# disp_data_all <- vector(mode = "list", length = n_samples)
# disp_data_all
# 
# i = 1
# # Loop over all stochastic maps. (NOT SETTING THE FULL LOOP YET)
# for (i in 1:100) {
#   
#   # Initialize list to store results.
#   char_list <- setNames(vector(mode = "list", length = length(CH_IDs)), CH_IDs)
#   
#   cat(paste0("\n", "Working on map: ", i, " ", Sys.time(), "\n"))
#   # Get timing of states for all characters.
#   for (j in 1:length(char_list)) { 
#     
#     # Import stochastic map.
#     stm <- readRDS(paste0("stmaps/", CH_IDs[[j]], ".RDS"))
#     
#     char_list[[j]] <- get_time_seg(stm[[i]]) 
#     
#   }
#   cat(paste0("\n", "Working on map: ", j, " ", Sys.time(), "\n"))
#   
#   # Make disparity data for one sample of stochastic maps.
#   disp_data <- make_dispar_data(char_list, n_points = 1000)
#   
#   # Store results.
#   disp_data_all[[i]] <- disp_data
#   
# }
# 
# # Calculate CI. Build final tibble.
# disp_data_all <- disp_data_all[1:50]

# Calculate CI. Build final tibble.
tb_dispar_summ <- tibble("Time" = apply(sapply(disp_data_all, function(x) x$Time ), 1, median),
                         "Median_disparity" = apply(sapply(disp_data_all, function(x) x$Disparity ), 1, median),
                         "Mean_disparity" = apply(sapply(disp_data_all, function(x) x$Disparity ), 1, mean),
                         "Low_CI_disparity" = apply(sapply(disp_data_all, function(x) x$Disparity ), 1, function(y) quantile(y, 0.025, na.rm = TRUE)),
                         "High_CI_disparity" = apply(sapply(disp_data_all, function(x) x$Disparity ), 1, function(y) quantile(y, 0.975, na.rm = TRUE)),
                         "Median_norm_disparity" = apply(sapply(disp_data_all, function(x) (x$Disparity/x$N_lineages) ), 1, median),
                         "Mean_norm_disparity" = apply(sapply(disp_data_all, function(x) (x$Disparity/x$N_lineages) ), 1, mean),
                         "Low_CI_norm_disparity" = apply(sapply(disp_data_all, function(x) (x$Disparity/x$N_lineages) ), 1, function(y) quantile(y, 0.025, na.rm = TRUE)),
                         "High_CI_norm_disparity" = apply(sapply(disp_data_all, function(x) (x$Disparity/x$N_lineages) ), 1, function(y) quantile(y, 0.975, na.rm = TRUE)),
                         "N_Lineages" = apply(sapply(disp_data_all, function(x) x$N_lineages ), 1, mean) )
tb_dispar_summ

# Save temporary data object.
#saveRDS(tb_dispar_summ, "data_out/data_disparity_frechet.RDS")
saveRDS(tb_dispar_summ, "data_out/data_disparity_var-hm.RDS")

#------------------------------------------------------------#

plot(tb_dispar_summ$Time, tb_dispar_summ$Mean_disparity, type='l', xlim=c(280, 0))
plot(tb_dispar_summ$Time, sqrt(tb_dispar_summ$Mean_disparity), type='l', xlim=c(280, 0))

plot(tb_dispar_summ$Time, tb_dispar_summ$Mean_norm_disparity, type='l', xlim=c(280, 0))

plot(tb_dispar_summ$Time, tb_dispar_summ$N_Lineages, type='l', xlim=c(280, 0))

tb_dispar_summ <- tb_dispar_summ %>%
  arrange(Time) %>%                # make sure time is increasing (or decreasing)
  mutate(
    dDispar = c(NA, diff(Mean_disparity)),      # change in disparity
    dTime   = c(NA, diff(Time)),                # change in time
    deriv   = dDispar / dTime                   # derivative: ΔDispar / ΔTime
  )

plot(tb_dispar_summ$Time, -tb_dispar_summ$deriv, type='l', xlim=c(280, 0))

### QUICK PLOT JUST TO CHECK.

# Make base plot.
base_plot <- tb_dispar_summ %>% ggplot(aes(x = Time, y =  Mean_norm_disparity)) + 

  geom_ribbon(data = tb_dispar_summ, aes(x = Time, ymin = Low_CI_norm_disparity, ymax = High_CI_norm_disparity), fill = "blue", alpha = 0.2) +  
    
  geom_line(data = tb_dispar_summ, aes(x = Time, y = Mean_norm_disparity), color = "black", linewidth = 0.8) +
  
  # Custom scale. Adjust if necessary.
  scale_x_reverse(breaks = seq(0, 280, by = 50)) +

  theme_minimal(base_size = 14) + 
  
  labs(title = "Mean Disparity Over Time", x = "Time", y = "Normalized Disparity") + 
  
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank())

# Save plot.
base_plot

# Save plot.
#ggsave("disparity_through_time_NEW.png", width = 12, height = 8, dpi = 300, bg = "white")
