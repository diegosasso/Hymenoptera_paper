
summary_df <- readRDS(file='data_out/neff-per-replicate.RDS')
summary_df_old <- readRDS(file='data_out/neff-per-replicate-old.RDS')

head(summary_df_old)
head(summary_df)

tail(summary_df_old)
tail(summary_df)

#---
mat <- replicate_list[[1]]
quantile(mat[,1], probs = c(0.05, 0.95), na.rm = TRUE)

q <- map(replicate_list, function(mat) {
  out <- mat[,1]
  out
})

q_un <- unlist(q)
quantile(q_un, probs = c(0.05, 0.95), na.rm = TRUE)
hist(q_un, breaks = 40)
lapply(q, function(x) x[1])


n_eff_list <- map(replicate_list, function(mat) {
  enr <- get_enrichment_mask(mat, lower_q = 0.05, upper_q = 0.95,
                             direction = "over", return_type = "binary")
  out <- plot_enrichment_entropy(tree, enr, type = "Neff_scaled",
                                 n_points = 1000)$data
  out
})

n_eff_df <- bind_rows(n_eff_list, .id = "replicate")
head(n_eff_df)