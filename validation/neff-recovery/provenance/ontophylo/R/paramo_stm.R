#' Sample stochastic maps from independently fitted character models
#'
#' Map each successful character with [phytools::make.simmap()] using its fitted
#' Q, observations and root setting. Save batches to disk or return them in memory.
#'
#' @param fitted_models An `ontophylo_fit` returned by [fit_models()].
#' @param n_stm Number of maps per character; the same replicate IDs are retained
#'   across characters for downstream amalgamation.
#' @param batch_size Maximum maps generated in one batch. Changing batch size
#'   does not change histories for a fixed seed and package versions.
#' @param workers Number of concurrent character jobs (PSOCK, also on Windows).
#' @param seed Nonnegative integer. Independent character-specific streams are
#'   reproducible across scheduling, batches and resumption; caller RNG is restored.
#' @param outdir Optional output directory. Disk runs return paths rather than
#'   collecting all histories in RAM. Use [read_stm()] to load selected characters.
#' @param resume Continue matching checkpoints, including partial characters.
#'   Increasing `n_stm` extends a run; decreasing it below saved counts is rejected.
#' @param store `"maps"`, `"branch_counts"` or `"both"`. Counts retain each
#'   branch and replicate, but cannot reproduce transition times or stacked maps.
#' @param timeout Seconds per mapping batch, or Inf. R.utils provides a
#'   cooperative limit, not a hard deadline for native code.
#' @param on_error Record failures and warn (`"continue"`) or stop after collecting
#'   jobs (`"stop"`). Successful checkpoints are preserved in either case.
#' @return An `ontophylo_stm` list containing `tree`, per-character results,
#'   `status`, settings and versions. In-memory character results contain `maps`
#'   and/or `branch_counts`; disk results contain `files`. Failed or skipped fits
#'   remain visible. [read_stm()] refuses incomplete selections.
#' @details These are conditional histories given the selected model, estimated
#'   Q and tree. They do not include parameter, model-selection or tree uncertainty.
#'   Combining replicate m across independently fitted characters samples their
#'   product distribution; dependent characters require a joint model first.
#'   Forbidden Q entries remain zero. Root settings are translated to phytools:
#'   equal, stationary (`estimated`), or FitzJohn (`fitzjohn`). Numeric probabilities
#'   retain their fitted state order. Existing checkpoints are never reused under
#'   different observations, Q, root settings, seeds, storage modes or versions.
#' @seealso [fit_models()], [read_stm()], [branch_rates()]
#' @examples
#' tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
#' x <- c(a = "0", b = "0", c = "1", d = "1")
#' fits <- fit_models(tree, x, models = "ER", root_p = "equal",
#'                    criterion = "AIC", nstarts = 0)
#' mapped <- paramo_stm(fits, n_stm = 5, seed = 9, store = "both")
#' rates <- branch_rates(read_stm(mapped))
#' @export
paramo_stm <- function(fitted_models, n_stm = 100L, batch_size = 25L,
                       workers = 1L, seed = 1L, outdir = NULL, resume = FALSE,
                       store = c("maps", "branch_counts", "both"), timeout = Inf,
                       on_error = c("continue", "stop")) {
  if (!inherits(fitted_models, "ontophylo_fit")) {
    stop("fitted_models must be the result of fit_models().", call. = FALSE)
  }
  n_stm <- .asr_integer(n_stm, "n_stm")
  batch_size <- .asr_integer(batch_size, "batch_size")
  workers <- .asr_integer(workers, "workers")
  seed <- .asr_integer(seed, "seed", 0L)
  resume <- .asr_flag(resume, "resume")
  store <- match.arg(store)
  on_error <- match.arg(on_error)
  if (!is.numeric(timeout) || length(timeout) != 1L || is.na(timeout) || timeout <= 0) {
    stop("timeout must be positive or Inf.", call. = FALSE)
  }
  outdir <- .asr_directory(outdir)
  if (resume && is.null(outdir)) stop("resume requires outdir.", call. = FALSE)
  tree <- .asr_tree(fitted_models$tree)
  versions <- .asr_versions()
  characters <- fitted_models$characters
  .asr_names(names(characters), "Character names")
  jobs <- vector("list", length(characters))
  for (i in seq_along(characters)) {
    id <- names(characters)[i]
    job <- list(id = id, tree = tree, fit = characters[[i]], n_stm = n_stm,
                  batch_size = batch_size, seed = seed, store = store,
                  timeout = timeout, resume = resume, directory = NULL)
    job$fingerprint <- .asr_hash(list(tree, id, characters[[i]]$Q_matrix,
      characters[[i]]$tips_prior, characters[[i]]$root_p,
      characters[[i]]$index_matrix, seed, store, versions))
    if (!is.null(outdir)) job$directory <- file.path(outdir, .asr_hash(id))
    jobs[[i]] <- job
  }
  results <- .asr_run_jobs(jobs, .asr_map_job, workers)
  names(results) <- names(characters)
  status <- data.frame(character = names(characters), status = "", n_maps = 0L,
                       reason = "", stringsAsFactors = FALSE)
  for (i in seq_along(results)) {
    status$status[i] <- results[[i]]$status
    status$n_maps[i] <- results[[i]]$n_maps
    status$reason[i] <- results[[i]]$reason
  }
  result <- list(tree = tree, characters = results, status = status,
                   settings = list(n_stm = n_stm, store = store, seed = seed),
                   versions = versions)
  class(result) <- "ontophylo_stm"
  .asr_check_failures(status, on_error)
  result
}

.asr_map_job <- function(job) {
  old_rng <- .asr_rng_save()
  on.exit(.asr_rng_restore(old_rng), add = TRUE)
  if (job$fit$status != "ok") {
    return(list(id = job$id, status = "skipped", n_maps = 0L,
                 reason = paste("Fit unavailable:", job$fit$reason)))
  }
  tryCatch(.asr_map_character(job), error = .asr_map_error)
}

.asr_map_error <- function(e) {
  list(status = "failed", n_maps = 0L, reason = conditionMessage(e))
}

.asr_mapping_inputs <- function(fit, tree) {
  Q <- fit$Q_matrix
  if (!is.matrix(Q) || !is.numeric(Q) || nrow(Q) < 2L || ncol(Q) != nrow(Q) ||
      any(!is.finite(Q)) || any(Q[row(Q) != col(Q)] < 0) ||
      any(abs(rowSums(Q)) > 1e-8 * max(1, max(abs(Q))))) stop("Invalid fitted Q.")
  .asr_names(rownames(Q), "Q state names")
  if (!identical(rownames(Q), colnames(Q))) stop("Q state order does not match.")
  tip <- fit$tips_prior
  if (!is.matrix(tip) || !is.numeric(tip) || any(!is.finite(tip)) ||
      any(tip < 0) || any(rowSums(tip) <= 0) ||
      !identical(rownames(tip), tree$tip.label) ||
      !identical(colnames(tip), colnames(Q))) stop("Tip likelihoods do not match Q/tree.")
  root <- .asr_root(fit$root_p, colnames(Q))
  if (identical(root, "yang")) root <- "estimated"
  if (identical(root, "maddfitz")) root <- "fitzjohn"
  list(Q = Q, x = tip, pi = root)
}

.asr_counts <- function(maps, tree) {
  counts <- matrix(0L, nrow(tree$edge), length(maps))
  rownames(counts) <- paste0("branch_", seq_len(nrow(tree$edge)))
  for (i in seq_along(maps)) {
    if (!identical(maps[[i]]$edge, tree$edge) ||
        !identical(maps[[i]]$tip.label, tree$tip.label) ||
        !isTRUE(all.equal(maps[[i]]$edge.length, tree$edge.length))) {
      stop("Mapping engine changed the reference tree or edge order.")
    }
    for (b in seq_len(nrow(tree$edge))) {
      labels <- names(maps[[i]]$maps[[b]])
      if (length(labels) > 1L) counts[b, i] <- sum(labels[-1L] != labels[-length(labels)])
    }
  }
  counts
}

.asr_map_character <- function(job) {
  inputs <- .asr_mapping_inputs(job$fit, job$tree)
  .asr_rng_start(job$seed, job$id, "mapping")
  completed <- 0L
  files <- character()
  warnings <- character()
  maps <- NULL
  counts <- NULL
  if (job$store %in% c("maps", "both")) maps <- list()
  if (job$store %in% c("branch_counts", "both")) {
    counts <- matrix(0L, nrow(job$tree$edge), 0L,
                       dimnames = list(paste0("branch_", seq_len(nrow(job$tree$edge))), NULL))
  }
  if (!is.null(job$directory)) {
    if (!dir.exists(job$directory)) dir.create(job$directory, recursive = TRUE)
    metadata_file <- file.path(job$directory, "metadata.rds")
    if (file.exists(metadata_file)) {
      if (!job$resume) stop("Mapping checkpoint exists; use resume=TRUE or a new outdir.")
      metadata <- readRDS(metadata_file)
      if (!identical(metadata$fingerprint, job$fingerprint)) stop("Mapping checkpoint mismatch.")
      files <- list.files(job$directory, "^chunk_[0-9]+_[0-9]+[.]rds$", full.names = TRUE)
      files <- sort(files)
      for (file in files) {
        chunk <- readRDS(file)
        if (!identical(chunk$fingerprint, job$fingerprint) ||
            chunk$first != completed + 1L || chunk$last < chunk$first) {
          stop("Noncontiguous or incompatible mapping checkpoint.")
        }
        completed <- chunk$last
        assign(".Random.seed", chunk$rng, envir = .GlobalEnv)
        warnings <- unique(c(warnings, chunk$warnings))
      }
      if (completed > job$n_stm) stop("Requested n_stm is smaller than saved map count.")
    } else {
      if (length(list.files(job$directory, all.files = FALSE))) {
        stop("Mapping directory is not empty and has no metadata.")
      }
      .asr_save(list(id = job$id, fingerprint = job$fingerprint), metadata_file)
    }
  }
  while (completed < job$n_stm) {
    first <- completed + 1L
    last <- min(job$n_stm, completed + job$batch_size)
    args <- list(tree = job$tree, x = inputs$x, Q = inputs$Q,
                  pi = inputs$pi, model = job$fit$index_matrix,
                  nsim = last - first + 1L, message = FALSE, tol = 0)
    attempt <- .asr_capture(phytools::make.simmap, args, job$timeout)
    warnings <- unique(c(warnings, attempt$warnings))
    if (!is.null(attempt$error)) {
      return(list(id = job$id, status = "failed", reason = attempt$error,
                    n_maps = completed, files = files, maps = maps,
                    branch_counts = counts, warnings = warnings))
    }
    batch <- attempt$value
    if (inherits(batch, "simmap")) batch <- list(batch)
    names(batch) <- paste0("map_", first:last)
    batch_counts <- .asr_counts(batch, job$tree)
    colnames(batch_counts) <- names(batch)
    chunk <- list(fingerprint = job$fingerprint, first = first, last = last,
                    rng = get(".Random.seed", envir = .GlobalEnv),
                    warnings = attempt$warnings)
    if (job$store %in% c("maps", "both")) chunk$maps <- batch
    if (job$store %in% c("branch_counts", "both")) chunk$branch_counts <- batch_counts
    if (is.null(job$directory)) {
      if (!is.null(chunk$maps)) maps <- c(maps, batch)
      if (!is.null(chunk$branch_counts)) counts <- cbind(counts, batch_counts)
    } else {
      file <- file.path(job$directory, sprintf("chunk_%010d_%010d.rds", first, last))
      .asr_save(chunk, file)
      files <- c(files, file)
    }
    completed <- last
  }
  if (!is.null(job$directory)) { maps <- NULL; counts <- NULL }
  if (!is.null(maps)) class(maps) <- c("multiSimmap", "multiPhylo")
  list(id = job$id, status = "ok", reason = "", n_maps = completed,
         files = files, maps = maps, branch_counts = counts, warnings = warnings)
}

#' Read selected stochastic-map outputs
#'
#' @param x An `ontophylo_stm` object from [paramo_stm()].
#' @param characters Character IDs to load; NULL selects all. Failed/skipped
#'   characters must be excluded explicitly before loading successful subsets.
#' @param type `"maps"` returns a named list of `multiSimmap` objects suitable for
#'   [branch_rates()] and PARAMO. `"branch_counts"` returns one branch-by-replicate
#'   matrix per character. Counts can also be extracted from stored maps.
#' @return A named list in the requested character order. This operation loads
#'   the selected outputs into memory, so select smaller groups for large runs.
#' @export
read_stm <- function(x, characters = NULL, type = c("maps", "branch_counts")) {
  if (!inherits(x, "ontophylo_stm")) stop("x must be a paramo_stm result.")
  type <- match.arg(type)
  if (is.null(characters)) characters <- names(x$characters)
  .asr_names(characters, "Selected character IDs")
  if (!all(characters %in% names(x$characters))) stop("Unknown character IDs.")
  result <- vector("list", length(characters))
  names(result) <- characters
  for (i in seq_along(characters)) {
    item <- x$characters[[characters[i]]]
    if (item$status != "ok" || item$n_maps != x$settings$n_stm) {
      stop("Incomplete character: ", characters[i], ". Select successful characters explicitly.")
    }
    value <- item[[type]]
    if (!is.null(item$files) && length(item$files)) {
      value <- NULL
      n <- 0L
      for (file in item$files) {
        chunk <- readRDS(file)
        if (chunk$first != n + 1L) stop("Noncontiguous mapping files.")
        part <- chunk[[type]]
        if (is.null(part) && type == "branch_counts" && !is.null(chunk$maps)) {
          part <- .asr_counts(chunk$maps, x$tree)
          colnames(part) <- names(chunk$maps)
        }
        if (is.null(part)) stop("Requested output was not stored: ", type)
        if (type == "maps") value <- c(value, part) else value <- cbind(value, part)
        n <- chunk$last
      }
      if (n != item$n_maps) stop("Stored map count does not match result.")
    }
    if (is.null(value) && type == "branch_counts" && !is.null(item$maps)) {
      value <- .asr_counts(item$maps, x$tree)
      colnames(value) <- names(item$maps)
    }
    if (is.null(value)) stop("Requested output was not stored: ", type)
    if (type == "maps") class(value) <- c("multiSimmap", "multiPhylo")
    result[[i]] <- value
  }
  result
}
