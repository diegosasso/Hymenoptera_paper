# Load libraries
library(phytools)
source('R/utils-branch-ama.R')



# 1. Simulate a random tree with 50 tips
tree <- pbtree(n = 5, scale = 1)
plot(tree)

# 2. Define a transition rate matrix (Q) for a binary trait
Q <- matrix(c(-1, 1,
              1, -1), 2, 2)

# 3. Simulate three discrete characters
char1 <- sim.history(tree, Q)$states
char2 <- sim.history(tree, Q)$states
char3 <- sim.history(tree, Q)$states


# 4. Fit a model for each character and generate stochastic maps
maps1 <- make.simmap(tree, char1, model = "SYM", nsim = 100)
maps2 <- make.simmap(tree, char2, model = "SYM", nsim = 5)
maps3 <- make.simmap(tree, char3, model = "SYM", nsim = 5)

plot(maps1[[2]])
edgelabels()

#-----
stm <- maps1[[2]]
plot(stm)
edgelabels()
get_time_seg(stm)

char_list <- lapply(maps1, function(x) get_time_seg(x))

char_list <- lapply(maps1[2], function(x) get_time_seg(x))
disp <- make_dispar_data(char_list, n_points = 100, time_points = NULL) 

plot(disp$Time, disp$Disparity, type='l')
lines(disp$Time, disp$N_lineages, col='red')
