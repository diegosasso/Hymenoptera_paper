
xx <- get_n_branches(tree)
plot(xx$time, xx$edge_len, type='l', xlim=c(280,0))
#----
adult <- apply(rates[,1:15], 1, sum)
rate_ad_lar <- cbind(adult, larva=rates[,16])
cor_pairwise(rate_ad_lar, 
             method = "pearson", 
             log_transform = T,
             log_const = 1e-6,
             trim_quantiles = c(0.01, 0.99))


corr <- cor_pairwise(rates, 
                       method = "pearson", 
                       log_transform = T,
                       log_const = 1e-6,
                       trim_quantiles = c(0.01, 0.99))

corr <- cor_pairwise(rates, 
                     method = "pearson", 
                     log_transform = T,
                     log_const = 1e-6,
                     trim_quantiles = NULL)

# cor_test <- cor_trimmed(rates, 
#                            method = "pearson", 
#                            log_transform = T, 
#                            log_const = 1e-6,
#                            trim_quantiles = c(0.01, 0.99),
#                            use = "pairwise.complete.obs")

corr$cor_sig
cor_test
#-----------------
library(ggplot2)
library(reshape2)

# Take absolute correlations and set diag to 0
mat <- abs(corr$cor_sig)
diag(mat) <- 0

# Preserve order
region_order <- colnames(mat)

# Melt for ggplot
df_heat <- melt(mat, varnames = c("Region1", "Region2"), value.name = "Correlation")

# Add a column indicating diagonal cells
df_heat$is_diag <- as.character(df_heat$Region1) == as.character(df_heat$Region2)

# Keep same order for both axes
df_heat$Region1 <- factor(df_heat$Region1, levels = region_order)
df_heat$Region2 <- factor(df_heat$Region2, levels = rev(region_order))  # flipped y-axis for symmetry

# Define maximum correlation
max_val <- max(df_heat$Correlation, na.rm = TRUE)

# Plot heatmap
ggplot(df_heat, aes(x = Region1, y = Region2, fill = Correlation)) +
  geom_tile(data = subset(df_heat, !is_diag),
            color = "grey85", linewidth = 0.3) +
  geom_tile(data = subset(df_heat, is_diag),
            fill = "grey", color = "grey85", linewidth = 0.3) +
  scale_fill_gradient(
    low = "white", high = "#b2182b",
    limits = c(0, max_val),
    name = "Abs. correlation"
  ) +
  scale_x_discrete(position = "top") +
  coord_fixed() +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 10)
  )

#-------- Cor vs Rate
ra <- abs(corr$cor_raw)
diag(ra) <- 0
apply(ra, 1, mean)

expand_pairwise_matrix <- function(mat, value_name = "value") {
  # Ensure it's a matrix
  mat <- as.matrix(mat)
  
  # Get the upper triangle indices (excluding diagonal)
  upper_idx <- which(upper.tri(mat), arr.ind = TRUE)
  
  # Build a data frame
  df <- data.frame(
    BR1 = rownames(mat)[upper_idx[, 1]],
    BR2 = colnames(mat)[upper_idx[, 2]],
    value = mat[upper_idx]
  )
  
  # Rename the value column
  names(df)[3] <- value_name
  
  return(df)
}
ra <- abs(corr$cor_raw)
tb <- expand_pairwise_matrix(ra)
(16^2-16)/2

load('data/char_num.RDA')
#anat_ent_state
# log of states per character per BR
# log_states <- lapply(anat_ent_state, function(x) sum(log(get_len(x))) ) %>% unlist
# log_states
tot_states <- lapply(anat_ent_state, function(x) sum(get_len(x)) ) %>% unlist
tot_states

rate_total <- apply(rates, 2, sum)
rate_total_norm <- rate_total/tot_states
# rate_total_norm <- exp(log(rate_total)-log_states)
rate_total_norm


#----

#cor_mat <- abs(corr$cor_sig)  # or raw correlations if signed values matter
cor_mat <- abs(corr$cor_raw)



# remove larva
cor_mat <-cor_mat[-16, -16]
diag(cor_mat) <- 0
mean(cor_mat)

group_traits <- c(
  "pronotum", "propectus", "mesonotum", 
  "mesopectus", "metanotum", "metapectal-propodeal complex"
)

# group_traits <- c(
#   "fore wing", "hind wing"
# )

# Extract all trait names
all_traits <- colnames(cor_mat)

# 1. Within-group correlations
within_vals <- cor_mat[group_traits, group_traits]
within_vals <- within_vals[upper.tri(within_vals, diag = FALSE)]

# 2. Between-group correlations
other_traits <- setdiff(all_traits, group_traits)
between_vals <- cor_mat[group_traits, other_traits]
between_vals <- as.vector(between_vals)

# 3. Compare with a statistical test
t_res <- t.test(within_vals, between_vals, alternative = "greater")
t_res
# 4. Summary
cat("Mean within-group correlation:", mean(within_vals, na.rm = TRUE), "\n")
cat("Mean between-group correlation:", mean(between_vals, na.rm = TRUE), "\n")
cat("T-test p-value (greater):", signif(t_res$p.value, 4), "\n")


#----------- Permutation test
permutation_test_modularity <- function(cor_mat, group_traits, n_perm = 10000, abs_vals = TRUE) {
  # cor_mat: correlation matrix (square, symmetric)
  # group_traits: vector of column names for focal module
  # n_perm: number of permutations
  # abs_vals: use absolute correlations (TRUE) or signed (FALSE)
  
  if (abs_vals) cor_mat <- abs(cor_mat)
  
  all_traits <- colnames(cor_mat)
  other_traits <- setdiff(all_traits, group_traits)
  
  # observed mean difference
  within_vals  <- cor_mat[group_traits, group_traits][upper.tri(cor_mat[group_traits, group_traits])]
  between_vals <- as.vector(cor_mat[group_traits, other_traits])
  delta_obs <- mean(within_vals, na.rm = TRUE) - mean(between_vals, na.rm = TRUE)
  
  # permutation loop
  deltas <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    # randomly assign same number of traits to "group"
    perm_group <- sample(all_traits, length(group_traits))
    perm_other <- setdiff(all_traits, perm_group)
    
    wv <- cor_mat[perm_group, perm_group][upper.tri(cor_mat[perm_group, perm_group])]
    bv <- as.vector(cor_mat[perm_group, perm_other])
    deltas[i] <- mean(wv, na.rm = TRUE) - mean(bv, na.rm = TRUE)
  }
  
  # empirical one-tailed p-value
  p_emp <- mean(deltas >= delta_obs)
  
  # return summary
  list(
    delta_obs = delta_obs,
    p_value = p_emp,
    mean_within = mean(within_vals, na.rm = TRUE),
    mean_between = mean(between_vals, na.rm = TRUE),
    deltas = deltas
  )
}

group_traits <- c(
  "pronotum", "propectus", "mesonotum",
  "mesopectus", "metanotum", "metapectal-propodeal complex"
)

res_perm <- permutation_test_modularity(corr$cor_sig, group_traits, n_perm = 10000)

cat("Observed Δ =", round(res_perm$delta_obs, 3), "\n")
cat("p-value =", signif(res_perm$p_value, 3), "\n")
cat("Within mean =", round(res_perm$mean_within, 3),
    "Between mean =", round(res_perm$mean_between, 3), "\n")
