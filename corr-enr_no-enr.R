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

enr.boolean
rates_enr <- rates[enr.boolean==1,]
rates_no <- rates[enr.boolean==0,]

#----------------
# corr <- cor_pairwise(rates_enr,
#                      method = "pearson",
#                      log_transform = T,
#                      log_const = 1e-6,
#                      trim_quantiles = c(0.01, 0.99))

corr <- cor_pairwise(rates_no,
                     method = "pearson",
                     log_transform = T,
                     log_const = 1e-6,
                     trim_quantiles = c(0.01, 0.99))

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
p <- ggplot(df_heat, aes(x = Region1, y = Region2, fill = Correlation)) +
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

p
# ggsave("figures/cor/cor-pairwise-enriched.pdf", p, width = 8, height = 8, dpi = 300, units = "in")
# ggsave("figures/cor/cor-pairwise-no-enriched.pdf", p, width = 8, height = 8, dpi = 300, units = "in")

corr_enr <- cor_pairwise(rates_enr,
                     method = "pearson",
                     log_transform = T,
                     log_const = 1e-6,
                     trim_quantiles = c(0.01, 0.99))

corr_no <- cor_pairwise(rates_no, 
                     method = "pearson", 
                     log_transform = T,
                     log_const = 1e-6,
                     trim_quantiles = c(0.01, 0.99))


corr_enr$cor_sig - corr_no$cor_sig


