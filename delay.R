cor_pairwise <- function(rates, 
                        method = "pearson", 
                        log_transform = FALSE, 
                        log_const = 1e-6,
                        trim_quantiles = NULL) {
  # rates: matrix (branches × body regions)
  # method: correlation method ("pearson", "spearman", etc.)
  # log_transform: apply log(x + log_const)
  # trim_quantiles: e.g. c(0.05, 0.95) removes outliers per column
  
  stopifnot(is.matrix(rates) || is.data.frame(rates))
  mat <- as.matrix(rates)
  p <- ncol(mat)
  
  # optional log transform
  if (log_transform) {
    mat <- log(mat + log_const)
  }
  
  # optional trimming per column
  if (!is.null(trim_quantiles)) {
    for (j in seq_len(p)) {
      col_vals <- mat[, j]
      qs <- quantile(col_vals, probs = trim_quantiles, na.rm = TRUE)
      keep <- col_vals >= qs[1] & col_vals <= qs[2]
      mat[!keep, j] <- NA
    }
  }
  
  # prepare result matrices
  cor_mat <- matrix(NA_real_, nrow = p, ncol = p,
                    dimnames = list(colnames(mat), colnames(mat)))
  p_mat <- cor_mat
  
  # compute pairwise correlations
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      x <- mat[, i]
      y <- mat[, j]
      keep <- is.finite(x) & is.finite(y)
      
      if (sum(keep) >= 3) {
        test <- suppressWarnings(cor.test(x[keep], y[keep], method = method, use = "pairwise.complete.obs"))
        cor_mat[i, j] <- cor_mat[j, i] <- unname(test$estimate)
        p_mat[i, j]   <- p_mat[j, i]   <- test$p.value
      }
    }
  }
  diag(cor_mat) <- 1
  diag(p_mat) <- 0
  
  # Bonferroni correction
  m_tests <- p * (p - 1) / 2
  p_bonf <- p.adjust(p_mat[upper.tri(p_mat)], method = "bonferroni")
  
  p_bonf_mat <- matrix(NA_real_, nrow = p, ncol = p,
                       dimnames = dimnames(p_mat))
  p_bonf_mat[upper.tri(p_bonf_mat)] <- p_bonf
  p_bonf_mat[lower.tri(p_bonf_mat)] <- t(p_bonf_mat)[lower.tri(p_bonf_mat)]
  diag(p_bonf_mat) <- 0
  
  list(
    cor   = cor_mat,
    p_raw = p_mat,
    p_bonf = p_bonf_mat
  )
}