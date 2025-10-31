
# Load packages.
library(phytools)
library(ontophylo)
library(tidyverse)

# Load data.
load("data_out/paramo_stm_adult_final.RDA")
load("data_out/paramo_amalg_adult.RDA")

# Set some parameters.
res = 200
n_samples = 100

#--------------------------------------------------#

### INITIAL ORGANIZATION ###

# Get global timings.
hym_tree_discr <- discr_Simmap(hym_tree, res)
node_heights <- nodeHeights(hym_tree_discr)[,1]
br_sums <- lapply(hym_tree_discr$maps, function(x) cumsum(x) )
times <- mapply(x = node_heights, y = br_sums, function(x,y) round((x + y),4) )

# Import first stochastic map.
stm <- readRDS(paste0("stmaps/CH1.RDS"))

# Get a sample of stochastic maps.
stm <- stm[seq((length(stm)/n_samples), length(stm), (length(stm)/n_samples))]

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
for (i in 2:length(CH_IDs)) {
  
  # Import stochastic map.
  stm <- readRDS(paste0("stmaps/", CH_IDs[[i]], ".RDS"))
  
  # Get a sample of stochastic maps.
  stm <- stm[seq((length(stm)/n_samples), length(stm), (length(stm)/n_samples))]
  
  # Discretize stochastic map.
  stm_discr <- discr_Simmap_all(stm, res)
  
  # Extract all maps.
  br_maps <- lapply(stm_discr, function(x) lapply(x$maps, function(y) names(y) ) )
  
  cat(paste0("\n", "Working on: ", CH_IDs[[i]], " ", Sys.time(), "\n"))
  
  # Concatenate all branch maps for a new character across all stochastic maps. 
  for (j in 1:n_samples) {
    
    smaps[[j]] <- mapply(x = smaps[[j]], y = br_maps[[j]], function(x,y) rbind(x,y) )
    
  }
  
}

# Merge all individual states as single string.
smaps_merg <- lapply(smaps, function(x) lapply(x, function(y) apply(y,2,paste0, collapse = "") ) )

# Save temporary data object.
#saveRDS(smaps_merg, "data_out/smaps_disparity.RDS")
#smaps_merg <- readRDS("data_out/smaps_disparity.RDS")

# OBS.: alternatively one can simply import the amalgamated stochastic map and extract the mappings!

# Clean workspace.
rm(i, j, stm, stm_discr, smaps)

# Plot tree to check branch and node labels.
plot.phylo(hym_tree, show.tip.label = F)
edgelabels(frame = "none", col = "red", cex = 0.6)
nodelabels(frame = "none", col = "blue", cex = 0.6)

# Add timings to all amalgamated states.
smaps_final <- vector(mode = "list", length = n_samples)
for (i in 1:n_samples) { smaps_final[[i]] <- mapply(x = times, y = smaps_merg[[i]], function(x,y) setNames(x,y) ) }

# Pre-allocate list to store all Hamming distances.
hm_dist_all <- vector(mode = "list", length = n_samples)

# Loop over all stochastic maps.
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

# Quick checks.
apply(do.call(cbind, hm_dist_all),1,summary)

# Build tibble for plotting (already adding normalized disparity).
#tb_dispar <- lapply(hm_dist_all, function(x) tibble(Time = c(0, as.numeric(names(x))), 
#                                                    Disparity = c(0, round(unname(x),2)),
#                                                    Norm_Disparity = c(0, round(unname(x)/ltt,2))) )

# Build tibble for plotting.
tb_dispar <- lapply(hm_dist_all, function(x) tibble(Time = c(0, as.numeric(names(x))), 
                                                    Disparity = c(0, round(unname(x),2))) )
tb_dispar <- do.call(bind_rows, tb_dispar)
tb_dispar

# Save temporary data object.
#saveRDS(tb_dispar, "data_out/data_disparity.RDS")
#tb_dispar <- readRDS("data_out/data_disparity.RDS")

# Summarize data and add CI.
tb_dispar_summ <- tb_dispar %>% group_by(Time) %>% 
  
  summarise(Mean_Disparity = mean(Disparity), se = sd(Disparity)/sqrt(n()), .groups = "drop") %>% 
  
  mutate(CI_lower = Mean_Disparity - 1.96 * se, CI_upper = Mean_Disparity + 1.96 * se)

# Add column normalized by number of lineages.
tb_dispar_summ <- tb_dispar_summ %>% mutate(Norm_Disparity = Mean_Disparity/c(1,ltt))
tb_dispar_summ %>% mutate(N_lineages = c(1, ltt))



#------------------------------##------------------------------##------------------------------##------------------------------#
# Save temporary data object.
#saveRDS(tb_dispar_summ, "data_out/data_disparity_final.RDS")
tb_dispar_summ <- readRDS("data_out/data_disparity_final.RDS")
tb_dispar_summ

plot(tb_dispar_summ$Time, tb_dispar_summ$Norm_Disparity, type='l')

max(tb_dispar_summ$CI_upper)
plot(tb_dispar_summ$Time, tb_dispar_summ$Norm_Disparity, type='l', ylim=c(0, 108))
lines(tb_dispar_summ$Time, tb_dispar_summ$CI_upper, col='red')
lines(tb_dispar_summ$Time, tb_dispar_summ$CI_lower, col='blue')
#------------------------------#
# NORMALIZED DISPARITY #
#------------------------------#

# Make base plot.
base_plot <- tb_dispar_summ %>% ggplot(aes(x = Time, y =  Norm_Disparity)) + 
  
  #geom_smooth(method = "loess", linewidth = 0.1, color = "#B5EAD7", se = TRUE, fill = "#FFDAC1", alpha = 0.4) +
  #geom_smooth(method = "lm", formula = y ~ x + 0, linewidth = 0.1, color = "#B5EAD7", se = TRUE, fill = "#FFDAC1", alpha = 0.4) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), linewidth = 0.1, color = "#B5EAD7", se = TRUE, fill = "#FFDAC1", alpha = 0.4) +
  
  #geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), fill = "#FFDAC1", alpha = 0.8) + 
  #geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), fill = "black", alpha = 0.8) +
  
  #geom_line(linewidth = 1.2, color = "#B5EAD7") + 
  
  # Custom scale. Adjust if necessary.
  scale_x_continuous(breaks = c(seq(0, max(tb_dispar_summ$Time), by = 50), max(tb_dispar_summ$Time)),
                     labels = c(round(seq(max(tb_dispar_summ$Time), 0, by = -50),0), 0)) +
  
  theme_minimal(base_size = 14) + 
  
  labs(title = "Mean Disparity Over Time (splines-smooth)", x = "Time", y = "Normalized Disparity") + 
  
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank())

# Save plot.
base_plot
#ggsave("disparity_through_time_plot.png", width = 12, height = 8, dpi = 300, bg = "white")
