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
rate_whole_branch <- apply(rates, 1, sum)
#colnames(rates)
edge_times <- edge_times_from_tips(tree)
head(edge_times)

# rates <- rates[,c(3:8)]
# rates
# ------ RND
# rates_shuffled <- apply(rates, 2, sample)
# dim(rates_shuffled) <- dim(rates)
# rownames(rates_shuffled) <- rownames(rates)
# colnames(rates_shuffled) <- colnames(rates)
# rates_shuffled
# rates <- rates_shuffled
# ----

n_perm <- 500
rates
edge_times

null_slopes <- numeric(n_perm)

for (i in seq_len(n_perm)) {
  rates_shuffled <- rates
  # shuffle each column independently
  for (j in seq_len(ncol(rates_shuffled))) {
    rates_shuffled[, j] <- sample(rates_shuffled[, j])
  }
  #rates_shuffled
  cor_null_df <- cor_vs_time_summary(rates_shuffled, edge_times,
                                     window_size = 100,
                                     step_size = 1,
                                     trim_quantiles = c(0.01, 0.99),
                                     log_transform = TRUE,
                                     abs_cor = TRUE)
  cor_null_df <- cor_null_df %>%
    filter(!is.na(mean_cor))
  
  #null_slopes[i] <- 
}

#------
source("R/utils-ST-test.R")

cor_time_df <- cor_vs_time_summary(rates, edge_times,
                                   window_size = 100,
                                   step_size = 1,
                                   trim_quantiles = c(0.01, 0.99),
                                   log_transform = TRUE,
                                   abs_cor = TRUE)
head(cor_time_df)
cor_time_df_clean <- cor_time_df %>%
  filter(!is.na(mean_cor))

m_lin <- lm(mean_cor ~ time_center, data = cor_time_df_clean,
            weights = n_edges)
summary(m_lin)

m_lin <- lm(log(mean_cor) ~ time_center, data = cor_time_df_clean, weights = n_edges)
summary(m_lin)

m_poly <- lm(mean_cor ~ poly(time_center, 2), data = cor_time_df_clean, weights = n_edges)
summary(m_poly)

plot(cor_time_df_clean$time_center, cor_time_df_clean$mean_cor, type='l', xlim=c(280,0), ylim=c(0, 0.7))
lines(cor_null_df$time_center, cor_null_df$mean_cor, col='red')
plot(cor_time_df_clean$mean_cor, cor_null_df$mean_cor)

plot(cor_time_df_clean$time_center, cor_time_df_clean$mean_cor/cor_null_df$time_center, xlim=c(280,0), ylim=c(0, 0.01))

model_corrected <- lm(cor_time_df_clean$mean_cor ~ cor_null_df$mean_cor)
summary(model_corrected)
residuals_corrected <- resid(model_corrected)

cor_time_df_clean$null_corrected <- resid(model_corrected)
mean(cor_time_df_clean$null_corrected)
var(cor_time_df_clean$null_corrected)

ggplot(cor_time_df_clean, aes(x = time_center, y = null_corrected)) +
  geom_line() +
  geom_smooth(span = 0.3, se = FALSE, color = "black") +
  scale_x_reverse() +
  theme_classic() +
  labs(y = "Null-corrected mean correlation", x = "Time (Ma)")

#-----
plot(cor_time_df_clean$time_center, log(cor_time_df_clean$mean_cor), type='l', xlim=c(280,0))
plot(log(cor_time_df_clean$time_center), (cor_time_df_clean$mean_cor), type='l')

ggplot(cor_time_df_clean,
       aes(x = time_center, y = mean_cor)) +
  geom_point(aes(size = n_edges), alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_x_reverse() +  # geological time, if older → larger
  theme_classic() +
  labs(x = "Time (Ma)",
       y = "Mean |correlation| between body regions",
       size = "Edges\nin window")


slope_obs <- get_slope(rates, edge_times, window_size = 30, step_size = 2)
slope_obs

#-----
# Ensure we only ever deal with a numeric matrix
rates_mat <- as.matrix(rates)

if (!is.numeric(rates_mat)) {
  stop("rates must be numeric (all columns).")
}


n_perm <- 500
null_slopes <- numeric(n_perm)

for (i in seq_len(n_perm)) {
  rates_shuffled <- rates_mat
  # shuffle each column independently
  for (j in seq_len(ncol(rates_shuffled))) {
    rates_shuffled[, j] <- sample(rates_shuffled[, j])
  }
  
  null_slopes[i] <- get_slope(
    rates_shuffled,
    edge_times,
    window_size = 30,
    step_size = 2
  )
}
#-----

p_emp <- mean(null_slopes <= slope_obs)  # for negative trend
p_emp

p_emp_two <- mean(abs(null_slopes) >= abs(slope_obs))
p_emp_two

ggplot(data.frame(slope = null_slopes), aes(x = slope)) +
  geom_histogram(bins = 40, fill = "grey80", color = "black") +
  geom_vline(xintercept = slope_obs, color = "red", linewidth = 1) +
  theme_classic() +
  labs(x = "Slope of mean |cor| ~ time (null)",
       y = "Count",
       title = paste("Observed slope =", round(slope_obs, 4)))

#------- 

cor_results <- rolling_rate_correlation(rates, edge_times,
                                        window_size = 40,
                                        step_size = 1,
                                        trim_quantiles = NULL,
                                        log_transform = F,
                                        use.abs.value=TRUE)

df_avg <- summarize_cor_bins(cor_results)
df_avg
lines(df_avg$time, df_avg$mean_cor, col='red')

ggplot(df_avg, aes(x = time, y = median_cor)) +
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = q05, ymax = q95),
              alpha = 0.2, fill = "grey") +
  scale_x_reverse(
    breaks = seq(0, 300, by = 25),
    expand = c(0, 0)
  ) +
  theme_classic() +
  labs(y = "Median pairwise correlation", x = "Time (Ma)")


#------------------------ OLD
#rates <- rates[,c(3:8)]

#------ RND
# rates_shuffled <- apply(rates, 2, sample)
# dim(rates_shuffled) <- dim(rates)
# rownames(rates_shuffled) <- rownames(rates)
# colnames(rates_shuffled) <- colnames(rates)
# rates_shuffled
# rates <- rates_shuffled
#----





#bins_20 <- seq(0, 280, by = 100)
bins_20 <- c(0, 100, 200, 250, 280)

cor_bins <- fixed_bin_correlation(
  rates       = rates,
  edge_times  = edge_times,
  breaks      = bins_20,
  trim_quantiles = NULL, #c(0.01, 0.99),
  log_transform  = TRUE,
  use.abs.value  = TRUE
)

df_avg <- summarize_cor_bins(cor_bins)
df_avg

ggplot(df_avg, aes(x = time, y = median_cor)) +
  geom_point(color = "black") +
  geom_ribbon(aes(ymin = q05, ymax = q95),
              alpha = 0.2, fill = "grey") +
  scale_x_reverse(
    breaks = seq(0, 300, by = 25),
    expand = c(0, 0)
  ) +
  theme_classic() +
  labs(y = "Median pairwise correlation", x = "Time (Ma)")



