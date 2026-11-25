# Check the four examples in Extended Data Figure 4 and one modified example.
# Binary columns follow the figure: wings, thorax, head.
# Run from validation/neff-shuffling or from the repository root.

analysis_dir <- "."
if (file.exists("validation/neff-shuffling/scripts/check-worked-example.R")) {
  analysis_dir <- "validation/neff-shuffling"
}
analysis_dir <- normalizePath(analysis_dir)
repo_dir <- normalizePath(file.path(analysis_dir, "../.."))

scenario_1 <- rbind(c(1, 1, 1), c(1, 1, 1), c(1, 1, 1))
scenario_2 <- rbind(c(1, 0, 0), c(1, 0, 0), c(1, 0, 0))
scenario_3 <- rbind(c(0, 0, 1), c(0, 1, 0), c(1, 0, 0))
scenario_4 <- rbind(c(0, 0, 0), c(0, 0, 0), c(0, 0, 0))
modified_3 <- rbind(c(0, 0, 1), c(0, 1, 0), c(0, 0, 0))
examples <- list(scenario_1, scenario_2, scenario_3, scenario_4, modified_3)
names(examples) <- c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4",
                     "Scenario 3 with third lineage unaccelerated")

# Load only the two existing metric definitions, avoiding unrelated code.
source_code <- parse(file.path(repo_dir, "R/utils-ST.R"))
for (expression in source_code) {
  if (!is.call(expression) || !identical(expression[[1]], as.name("<-"))) next
  if (identical(expression[[2]], as.name("diversity_metric")) ||
      identical(expression[[2]], as.name("relative_null_metric"))) {
    eval(expression)
  }
}

results <- data.frame(scenario = names(examples), patterns = NA_character_,
                      f_plus = NA_real_, H_plus = NA_real_, Neff_plus = NA_real_,
                      Pc = NA_real_, Neff_code_error = NA_real_, Pc_code_error = NA_real_)
for (example in 1:length(examples)) {
  matrix <- examples[[example]]
  N <- nrow(matrix)
  L <- ncol(matrix)
  a <- rowSums(matrix)
  patterns <- character(N)
  for (branch in 1:N) patterns[branch] <- paste(matrix[branch, ], collapse = "")

  positive <- a > 0
  f_plus <- sum(positive) / N
  H_plus <- NA_real_
  Neff_plus <- 0
  if (any(positive)) {
    p <- as.numeric(table(patterns[positive])) / sum(positive)
    H_plus <- -sum(p * log(p))
    Neff_plus <- f_plus * exp(H_plus)
  }
  Pc <- mean(a / L)

  results$patterns[example] <- paste(patterns, collapse = ", ")
  results$f_plus[example] <- f_plus
  results$H_plus[example] <- H_plus
  results$Neff_plus[example] <- Neff_plus
  results$Pc[example] <- Pc
  results$Neff_code_error[example] <- abs(Neff_plus - diversity_metric(patterns, "Neff_scaled"))
  results$Pc_code_error[example] <- abs(Pc - relative_null_metric(patterns, "scaled"))
}

stopifnot(max(abs(results$Neff_plus - c(1, 1, 3, 0, 4/3))) < 1e-12,
          max(abs(results$Pc - c(1, 1/3, 1/3, 0, 2/9))) < 1e-12,
          max(results$Neff_code_error) < 1e-12, max(results$Pc_code_error) < 1e-12)
write.csv(results, file.path(analysis_dir, "source-data/worked-example-values.csv"),
          row.names = FALSE, na = "not evaluated")
print(results, digits = 6, row.names = FALSE)
