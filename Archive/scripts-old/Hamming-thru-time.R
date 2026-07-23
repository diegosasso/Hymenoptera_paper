library(rphenoscate)
library(expm)
library(stringdist)

Q1 <- initQ(c(0, 1), c(1,1))

Q <- amaSMM(Q1,Q1)

Q
start.state <- 1
expm(Q)


# Function to compute transition matrix for a time-varying rate multiplier
P_time_varying <- function(Q, rate_fun, t0, t1) {
  # Integrate rate over [t0, t1]
  R <- integrate(rate_fun, lower = t0, upper = t1)$value
  
  # Matrix exponential of Q scaled by integrated rate
  P <- expm(Q * R)
  return(P)
}


hamming_var_from_P_exact <- function(P, start_state = 1) {
  # Ensure matrix and names
  stopifnot(is.matrix(P))
  states <- colnames(P)
  probs  <- P[start_state, ]
  
  # Compute pairwise Hamming distance matrix
  D <- stringdistmatrix(states, states, method = "hamming")
  D <- as.matrix(D)
  
  # Expected Hamming distance
  E_D  <- sum(outer(probs, probs) * D)
  
  # Expected squared Hamming distance
  E_D2 <- sum(outer(probs, probs) * D^2)
  
  # Variance
  Var_D <- E_D2 - E_D^2
  
  Var_D
  # list(
  #   mean_hamming = E_D,
  #   var_hamming  = Var_D
  # )
}

# Example -----------------------------------------------------------


# Define a time-varying rate function, e.g., an exponential decay
# rate_fun <- function(t) {
#   0.55 + 0.45 * sin(2 * pi * 3 * t / 5)  # 3 full oscillations between 0.1 and 1
# }
rate_fun <- function(t) {
  0.1 + 0.9 * sin(pi * t / 5)
}
x <- seq(0, 30, .1)
y <- rate_fun(x)
plot(x,y, type='l')
 
# Compute transition probabilities from time 0 to 200
P <- P_time_varying(Q, rate_fun, 0, .1)
P
hamming_var_from_P_exact(P, start_state = 1)

library(purrr)

# --- Compute variance over time ---
time_points <- seq(0, 30, by = 0.01)

var_results <- map_dfr(time_points, function(t) {
  P <- P_time_varying(Q*.01, rate_fun, 0, t)
  res <- hamming_var_from_P_exact(P, start_state = 1)
  tibble(time = t, var=res)
})

# inspect
head(var_results)
plot(var_results$time, var_results$var, type='l', ylim=c(-0.12, 0.12))
lines(x,y*.1, col='red')
