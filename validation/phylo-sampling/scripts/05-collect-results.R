# Collect complete tasks only; refuse to call an incomplete collection final.
arguments <- commandArgs(trailingOnly = TRUE)
require_all <- "--final" %in% arguments
design <- readRDS("design/design.rds")
manifest <- read.csv("design/task-manifest.csv")
dir.create("summary", showWarnings = FALSE)
trial_rows <- list()
curve_rows <- list()
metric_rows <- list()
completed <- integer()
missing <- integer()
row_curve <- 0L
row_metric <- 0L
for (task in manifest$task) {
  filename <- sprintf("results/task-%04d/result.rds", task)
  if (!file.exists(filename)) { missing <- c(missing, task); next }
  x <- readRDS(filename)
  completed <- c(completed, task)
  stopifnot(x$task == task, x$settings$root_p == "maddfitz",
            x$settings$normalize == "none", x$settings$n_maps == 100L,
            !any(x$diagnostics$full$mapping_status$status == "failed"),
            !any(x$diagnostics$study$mapping_status$status == "failed"))
  full <- x$groups$full_all
  study <- x$groups$study_all
  paired <- data.frame(task = task, condition = x$condition, epoch = x$epoch,
    proportion = x$proportion, trial = x$trial,
    n_characters = x$n_characters, n_retained = x$n_retained,
    retained_fraction = x$n_retained/x$n_characters,
    full_epoch_rate = mean(full$metrics$target_epoch_rate),
    study_epoch_rate = mean(study$metrics$target_epoch_rate),
    full_reference_epoch_rate = full$reference$target_epoch_rate,
    study_reference_epoch_rate = study$reference$target_epoch_rate,
    full_rate_vs_reference = mean(full$metrics$target_epoch_rate) /
                             full$reference$target_epoch_rate,
    study_rate_vs_reference = NA_real_,
    study_full_reference_rate_ratio = study$reference$target_epoch_rate /
                                      full$reference$target_epoch_rate,
    full_planted_mean_age = mean(x$selected_nodes$age_ma),
    study_planted_mean_age = NA_real_,
    study_full_rate_ratio = mean(study$metrics$target_epoch_rate) /
                            mean(full$metrics$target_epoch_rate),
    study_full_change_ratio = mean(study$metrics$changes) / mean(full$metrics$changes),
    study_full_age_difference = NA_real_,
    full_inference_age_error = mean(full$metrics$mean_change_age) -
                               full$reference$mean_change_age,
    study_inference_age_error = NA_real_, matched_rate_ratio = NA_real_,
    elapsed_seconds = x$elapsed_seconds)
  if (x$n_retained > 0L) {
    paired$study_rate_vs_reference <- mean(study$metrics$target_epoch_rate) /
                                      study$reference$target_epoch_rate
    paired$study_planted_mean_age <- mean(x$selected_nodes$age_ma[
                                          x$selected_nodes$variable_study])
    paired$study_full_age_difference <- mean(study$metrics$mean_change_age) -
                                        mean(full$metrics$mean_change_age)
    paired$study_inference_age_error <- mean(study$metrics$mean_change_age) -
                                        study$reference$mean_change_age
    paired$matched_rate_ratio <- mean(study$metrics$target_epoch_rate) /
                                 mean(x$groups$full_matched$metrics$target_epoch_rate)
  }
  trial_rows[[length(trial_rows)+1L]] <- paired
  for (group_name in names(x$groups)) {
    g <- x$groups[[group_name]]
    row_curve <- row_curve + 1L
    curve_rows[[row_curve]] <- data.frame(task = task, condition = x$condition,
      epoch = x$epoch, proportion = x$proportion, trial = x$trial, group = group_name,
      n_characters = g$n_characters, g$marginal_summary)
    if (group_name == "study_all" && x$n_retained > 0L) {
      matched_study <- curve_rows[[row_curve]]
      matched_study$group <- "study_matched"
      row_curve <- row_curve + 1L
      curve_rows[[row_curve]] <- matched_study
    }
    for (metric in c("changes","target_epoch_rate","target_epoch_fraction",
                    "mean_change_age","correct_origin_fraction",
                    "outside_apocrita_fraction","branch_count_error")) {
      values <- g$metrics[[metric]]
      values <- values[is.finite(values)]
      estimate <- lower <- upper <- NA_real_
      if (length(values)) {
        estimate <- mean(values)
        lower <- quantile(values,0.025,names=FALSE)
        upper <- quantile(values,0.975,names=FALSE)
      }
      row_metric <- row_metric+1L
      metric_rows[[row_metric]] <- data.frame(task = task, condition = x$condition,
        epoch = x$epoch, proportion = x$proportion, trial = x$trial,
        group = group_name, metric = metric, mean = estimate, lower = lower, upper = upper)
    }
  }
}
stopifnot(length(completed)>0L)
trials <- do.call(rbind,trial_rows)
curves <- do.call(rbind,curve_rows)
metrics <- do.call(rbind,metric_rows)
write.csv(trials,"summary/paired-trials.csv",row.names=FALSE)
saveRDS(curves,"summary/marginal-trials.rds")
write.csv(metrics,"summary/metric-trials.csv",row.names=FALSE)
write.csv(data.frame(task=missing),"summary/missing-tasks.csv",row.names=FALSE)

# Empirical ranges across trials describe replicate variation, not a CI for a mean.
curve_summary <- list()
counter <- 0L
for(condition in design$conditions$condition) {
  for(group in c("full_all","study_all","full_matched","study_matched")) {
    selected <- curves[curves$condition==condition & curves$group==group,,drop=FALSE]
    if(!nrow(selected)) next
    for(age in design$age_grid) {
      one <- selected[selected$time==age,,drop=FALSE]
      counter <- counter+1L
      curve_summary[[counter]] <- data.frame(condition=condition, group=group, age=age,
        n_trials=nrow(one), mean=mean(one$mean),
        lower=quantile(one$mean,.025,names=FALSE),
        upper=quantile(one$mean,.975,names=FALSE),
        origin_reference=mean(one$origin_reference))
    }
  }
}
write.csv(do.call(rbind,curve_summary),"summary/marginal-summary.csv",row.names=FALSE)
writeLines(c(paste("Completed:",length(completed),"of",nrow(manifest),"paired trials."),
             paste("Missing:",length(missing)),
             "Across-trial intervals are empirical central 95% ranges."),
           "summary/collection-status.txt")
cat(readLines("summary/collection-status.txt"),sep="\n")
if(require_all && length(missing)>0L) stop("Final collection requires all 1500 trials.")

