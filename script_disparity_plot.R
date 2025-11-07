
# Load packages.
library(phytools)
library(ontophylo)
library(tidyverse)

# Load data.
# load("data_out/paramo_stm_adult_final.RDA")
# load("data_out/paramo_amalg_adult.RDA")

# Set some parameters.
nmaps = 1000
n_samples = 100
res = 1000
idx_s = seq(from = 1, to = nmaps, by = n_samples)
idx_e = tail(seq(from = 0, to = nmaps, by = n_samples), -1)

#--------------------------------------------------#

### INITIAL ORGANIZATION ###

# Get global timings.
hym_tree <-readRDS("data/hym_tree.RDS")
hym_tree_discr <- discr_Simmap(hym_tree, res)
node_heights <- nodeHeights(hym_tree_discr)[,1]
#nodeHeights(hym_tree)[,1]
br_sums <- lapply(hym_tree_discr$maps, function(x) cumsum(x) )
times <- mapply(x = node_heights, y = br_sums, function(x,y) round((x + y),5) )
times

times[1]

# Create a folder to store all temporary maps.
dir.create("stmaps_disp")

# Loop over all map samples.
for (i in 1:length(idx_s)) {
  
  # Import first stochastic map.
  stm <- readRDS(paste0("stmaps/CH1.RDS"))
  
  # Get a sample of stochastic maps.
  stm <- stm[idx_s[i]:idx_e[i]]
  
  # Discretize first stochastic map.
  stm_discr <- discr_Simmap_all(stm, res)
  
  # Extract all maps.
  br_maps <- lapply(stm_discr, function(x) lapply(x$maps, function(y) names(y) ) )
  
  # Set starting map.
  smaps <- br_maps
  
  # Get character IDs.
  CH_IDs <- hym_annot_final$char_id
  
  #--------------------------------------------------#
  
  # Loop over all remaining characters.
  for (j in 2:length(CH_IDs)) {
    
    # Import stochastic map.
    stm <- readRDS(paste0("stmaps/", CH_IDs[[j]], ".RDS"))
    
    # Get a sample of stochastic maps.
    stm <- stm[idx_s[i]:idx_e[i]]
    
    # Discretize stochastic map.
    stm_discr <- discr_Simmap_all(stm, res)
    
    # Extract all maps.
    br_maps <- lapply(stm_discr, function(x) lapply(x$maps, function(y) names(y) ) )
    
    cat(paste0("\n", "Working on: ", CH_IDs[[j]], " ", Sys.time(), "\n"))
    
    # Concatenate all branch maps for a new character across all stochastic maps. 
    for (k in 1:n_samples) {
      
      smaps[[k]] <- mapply(x = smaps[[k]], y = br_maps[[k]], function(x,y) rbind(x,y) )
      
    }
    
  }
  
  # Merge all individual states as single string.
  smaps_merg <- lapply(smaps, function(x) lapply(x, function(y) apply(y,2,paste0, collapse = "") ) )
  
  # Save temporary data object.
  saveRDS(smaps_merg, paste0("stmaps_disp/smaps_disp_", i, ".RDS"))
  
  # Clean workspace.
  rm(j, k, stm, stm_discr, smaps)
  
}

#------------------------------------------

# Re-Import all map samples.
k=1
for (k in 1:(nmaps/n_samples)) {
  
  cat(paste0("\n", "Working on sample: ", k, " ", Sys.time(), "\n"))
  
  # Import raw disparity information.
  smaps_merg <- readRDS(paste0("stmaps_disp/smaps_disp_", k, ".RDS"))
  length(smaps_merg)
  
  # Add timings to all amalgamated states.
  smaps_final <- vector(mode = "list", length = n_samples)
  i=1
  for (i in 1:n_samples) { 
    smaps_final[[i]] <- mapply(x = times, y = smaps_merg[[i]], function(x,y) setNames(x,y) ) 
    }
  
  # Pre-allocate list to store all Hamming distances.
  hm_dist_all <- vector(mode = "list", length = n_samples)
  
  # Loop over all stochastic maps.
  j=1
  for (j in 1:n_samples) {
    
    cat(paste0("\n", "Working on Map: ", j, " ", Sys.time(), "\n"))
    
    # Extract all bin groups for a given amalgamated stochastic map.
    t_bins <- split(sort(unlist(smaps_final[[j]])), f = as.factor(sort(unlist(smaps_final[[j]]))))
    
    # Remove singletons bins (due to approximations when discretizing maps).
    t_bins_clean <- t_bins[sapply(t_bins, function(x) length(x) > 1 )]
    
    # Calculate Hamming distances among all amalgamated states for each time bin.
    hm_dist_all[[j]] <- sapply(t_bins_clean, function(x) stringdist::stringdistmatrix(names(x), names(x), method = "hamming") %>% mean )
    
  }
  
  # Get number of lineages present at each time bin.
  ltt <- sapply(t_bins_clean, length)

  # Build tibble for plotting (already adding normalized disparity).
  tb_dispar <- lapply(hm_dist_all, function(x) tibble(Time = c(0, as.numeric(names(x))), 
                                                      Disparity = c(0, round(unname(x),2)),
                                                      Norm_Disparity = c(0, round(unname(x)/ltt,2)),
                                                      N_Lineages = c(1,ltt) ) )
  # Save temporary data object.
  saveRDS(tb_dispar, paste0("data_out/data_disp_", k, ".RDS"))
  
}

i = 1
# Re-import all tibbles.
tb_dispar <- list()
for (i in 1:(nmaps/n_samples)) { tb_dispar <- c(tb_dispar, readRDS(paste0("data_out/data_disp_", i, ".RDS"))) }

# Revert time axis. Calculate CI. Build final tibble.
tb_dispar_summ <- tibble("Time" = apply(sapply(tb_dispar, function(x) rev(x$Time) ), 1, median),
                         "Median_disparity" = apply(sapply(tb_dispar, function(x) x$Disparity ), 1, median),
                         "Mean_disparity" = apply(sapply(tb_dispar, function(x) x$Disparity ), 1, mean),
                         "Low_CI_disparity" = apply(sapply(tb_dispar, function(x) x$Disparity ), 1, function(y) quantile(y, 0.025, na.rm = TRUE)),
                         "High_CI_disparity" = apply(sapply(tb_dispar, function(x) x$Disparity ), 1, function(y) quantile(y, 0.975, na.rm = TRUE)),
                         "Median_norm_disparity" = apply(sapply(tb_dispar, function(x) x$Norm_Disparity ), 1, median),
                         "Mean_norm_disparity" = apply(sapply(tb_dispar, function(x) x$Norm_Disparity ), 1, mean),
                         "Low_CI_norm_disparity" = apply(sapply(tb_dispar, function(x) x$Norm_Disparity ), 1, function(y) quantile(y, 0.025, na.rm = TRUE)),
                         "High_CI_norm_disparity" = apply(sapply(tb_dispar, function(x) x$Norm_Disparity ), 1, function(y) quantile(y, 0.975, na.rm = TRUE)),
                         "N_Lineages" = apply(sapply(tb_dispar, function(x) rev(x$N_Lineages) ), 1, mean) )

# Save temporary data object.
saveRDS(tb_dispar_summ, "data_out/data_disparity_final.RDS")
tb_dispar_summ <- readRDS("data_out/data_disparity_final.RDS")
tb_dispar_summ

#------------------------------#
# NORMALIZED DISPARITY #
#------------------------------#

# Make base plot.
base_plot <- tb_dispar_summ %>% ggplot(aes(x = Time, y =  Mean_norm_disparity)) + 
  
  geom_ribbon(data = tb_dispar_summ, aes(x = Time, ymin = Low_CI_norm_disparity, ymax = High_CI_norm_disparity), fill = "blue", alpha = 0.2) +
  
  geom_line(data = tb_dispar_summ, aes(x = Time, y = Mean_norm_disparity), color = "black", linewidth = 0.8) +
  
  # Custom scale. Adjust if necessary.
  scale_x_reverse(breaks = seq(0, 280, by = 50)) +
  #scale_x_continuous(breaks = c(seq(0, max(tb_dispar_summ$Time), by = 50), max(tb_dispar_summ$Time)),
  #                   labels = c(round(seq(max(tb_dispar_summ$Time), 0, by = -50),0), 0)) +
  
  theme_minimal(base_size = 14) + 
  
  labs(title = "Mean Disparity Over Time (splines-smooth)", x = "Time", y = "Normalized Disparity") + 
  
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank())

# Save plot.
base_plot

ggsave("disparity_through_time_plot.png", width = 12, height = 8, dpi = 300, bg = "white")
