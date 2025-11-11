library(dplyr)
library(purrr)
library(tibble)

# rates: matrix/data.frame [branches x BRs]
# edge_times: data.frame with columns start, end (and maybe mid)
# window_size: width of time window (e.g. 30 My)
# step_size: step between window centers (e.g. 2 My)
# abs_cor: whether to use |correlation| or signed correlation

cor_vs_time_summary <- function(rates,
                                edge_times,
                                window_size = 30,
                                step_size = 2,
                                trim_quantiles = NULL,
                                log_transform = FALSE,
                                abs_cor = TRUE) {
  rates <- as.matrix(rates)
  
  # Add midpoints if not present
  if (!"mid" %in% names(edge_times)) {
    edge_times$mid <- (edge_times$start + edge_times$end) / 2
  }
  
  max_t <- max(edge_times$mid)
  min_t <- min(edge_times$mid)
  
  centers <- seq(min_t, max_t, by = step_size)
  
  out <- map_dfr(centers, function(c) {
    lower <- c - window_size / 2
    upper <- c + window_size / 2
    
    in_win <- edge_times$mid >= lower & edge_times$mid <= upper
    sub_rates <- rates[in_win, , drop = FALSE]
    n_edges <- sum(in_win)
    
    if (n_edges > 3) {
      # cor_mat <- cor(sub_rates, use = "pairwise.complete.obs", method = "pearson")
      cor_mat <- cor_trimmed(
        sub_rates,
        method = "pearson",
        log_transform = log_transform,
        log_const = 1e-6,
        trim_quantiles = trim_quantiles, #c(0.01, 0.99),
        use = "pairwise.complete.obs"
      )
      if (abs_cor) cor_mat <- abs(cor_mat)
      vals <- cor_mat[upper.tri(cor_mat)]
      
      tibble(
        time_center = c,
        n_edges = n_edges,
        mean_cor = mean(vals, na.rm = TRUE),
        median_cor = median(vals, na.rm = TRUE),
        q25 = quantile(vals, 0.25, na.rm = TRUE),
        q75 = quantile(vals, 0.75, na.rm = TRUE)
      )
    } else {
      tibble(
        time_center = c,
        n_edges = n_edges,
        mean_cor = NA_real_,
        median_cor = NA_real_,
        q25 = NA_real_,
        q75 = NA_real_
      )
    }
  })
  
  out <- out %>% arrange(time_center)
  out
}


get_slope <- function(rates, edge_times,
                      window_size = 30, step_size = 2,
                      abs_cor = TRUE) {
  df <- cor_vs_time_summary(rates, edge_times,
                            window_size = window_size,
                            step_size = step_size,
                            abs_cor = abs_cor)
  df <- df %>% filter(!is.na(mean_cor))
  if (nrow(df) < 3) return(NA_real_)
  
  m <- lm(mean_cor ~ time_center, data = df, weights = n_edges)
  coef(m)["time_center"]
}