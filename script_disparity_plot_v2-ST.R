
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

frechet_hamming_mean <- function(strings, weights = NULL) {
  stopifnot(is.character(strings), length(strings) >= 1)
  
  # All strings must be same length
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
  
  # Weighted majority vote per position
  consensus <- character(L)
  for (j in seq_len(L)) {
    tab <- tapply(weights, mat[, j], sum)
    consensus[j] <- names(which.max(tab))
  }
  
  # Compute mean Hamming distance from consensus
  consensus_str <- paste(consensus, collapse = "")
  dist_vec <- colSums(t(mat) != consensus)
  mean_dist <- mean(dist_vec)
  
  mean_dist
  # list(
  #   consensus = consensus_str,
  #   mean_distance = mean_dist
  # )
}

#char_list <- lapply(maps1, function(x) get_time_seg(x))
#disp <- make_dispar_data(char_list, n_points = 100, time_points = NULL) 

# Make disparity data.
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
    
    # Frechet mean
    hm <-frechet_hamming_mean(st)
    
    dispar_ls[[k]] <- tibble("Time" = t, "Disparity" = hm, "N_lineages" = length(st))
    
  }
  
  # Return results.
  return(do.call(bind_rows, dispar_ls) %>% mutate(Time = rev(Time)))
  
}

#------------------------------------------------------------#

# Loop over all characters and get timings of states.
#for (j in 1:length(char_list)) {
#  
#  cat(paste0("\n", "Working on: ", CH_IDs[[j]], " ", Sys.time(), "\n"))
#  
#  # Import stochastic map.
#  stm <- readRDS(paste0("stmaps/", CH_IDs[[j]], ".RDS"))
#  
#  # Get timing of states.
#  char_list[[j]] <- lapply(stm, get_time_seg)
#  
#}

# Set some parameters.
n_samples = 1000

# Initialize list to store results from all samples.
disp_data_all <- vector(mode = "list", length = length(n_samples))

i = 1
# Loop over all stochastic maps. (NOT SETTING THE FULL LOOP YET)
#for (i in 1:n_samples) {}

# Get character IDs.
CH_IDs <- hym_annot_final$char_id

# Initialize list to store results.
char_list <- setNames(vector(mode = "list", length = length(CH_IDs)), CH_IDs)

cat(paste0("\n", "Working on map: ", j, " ", Sys.time(), "\n"))
# Get timing of states for all characters.
j=2
for (j in 1:length(char_list)) { 
#for (j in 1:10) { 
  
  # Import stochastic map.
  stm <- readRDS(paste0("stmaps/", CH_IDs[[j]], ".RDS"))
  
  char_list[[j]] <- get_time_seg(stm[[i]]) 
  
}
cat(paste0("\n", "Working on map: ", j, " ", Sys.time(), "\n"))

# Make disparity data for one sample of stochastic maps.
disp_data <- make_dispar_data(char_list, n_points = 100)


# Store results.
#disp_data)all[[i]] <- disp_data


plot(disp_data$Time, disp_data$Disparity, type='l', xlim=c(280,0))
plot(disp_data$Time, disp_data$Disparity/disp_data$N_lineages, type='l', xlim=c(280,0))
#lines(disp$Time, disp$N_lineages, col='red')
#------------------------------------------------------------#

### QUICK PLOT JUST TO CHECK.

# Make base plot.
base_plot <- disp_data %>% ggplot(aes(x = Time, y =  Disparity/N_lineages)) + 
  
  geom_line(data = disp_data, aes(x = Time, y = Disparity/N_lineages), color = "black", linewidth = 0.8) +
  
  # Custom scale. Adjust if necessary.
  scale_x_reverse(breaks = seq(0, 280, by = 50)) +

  theme_minimal(base_size = 14) + 
  
  labs(title = "Mean Disparity Over Time (splines-smooth)", x = "Time", y = "Normalized Disparity") + 
  
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank())

# Save plot.
base_plot

# Save plot.
ggsave("disparity_through_time_NEW.png", width = 12, height = 8, dpi = 300, bg = "white")
