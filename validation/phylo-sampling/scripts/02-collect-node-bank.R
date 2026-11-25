source("scripts/load-packages.R")
design <- readRDS("design/design.rds")
n_chunks <- ceiling(nrow(design$nodes) / 40)
for (tree_name in c("full", "study")) {
  combined <- readRDS(sprintf("fits/chunk-01/%s.rds", tree_name))
  combined$characters <- list()
  combined$status <- combined$status[FALSE, ]
  for (chunk in seq_len(n_chunks)) {
    one <- readRDS(sprintf("fits/chunk-%02d/%s.rds", chunk, tree_name))
    stopifnot(identical(one$tree, combined$tree),
              identical(one$settings, combined$settings),
              identical(one$versions, combined$versions))
    combined$characters <- c(combined$characters, one$characters)
    combined$status <- rbind(combined$status, one$status)
  }
  stopifnot(identical(names(combined$characters), design$nodes$character),
            !any(combined$status$status == "failed"),
            combined$settings$root_p == "maddfitz")
  combined$status$q <- NA_real_
  combined$status$boundary <- FALSE
  for (id in combined$status$character[combined$status$status == "ok"]) {
    row <- match(id, combined$status$character)
    fit <- combined$characters[[id]]
    stopifnot(fit$root_p == "maddfitz",
              ontophylo:::.asr_mapping_inputs(fit, combined$tree)$pi == "fitzjohn")
    combined$status$q[row] <- fit$Q_matrix["0", "1"]
    combined$status$boundary[row] <- fit$diagnostics$boundary
  }
  saveRDS(combined, sprintf("design/node-fits-%s.rds", tree_name))
  write.csv(combined$status, sprintf("design/node-fits-%s.csv", tree_name), row.names = FALSE)
  print(table(combined$status$status))
  cat("Boundary fits:",sum(combined$status$boundary),"\n")
}

