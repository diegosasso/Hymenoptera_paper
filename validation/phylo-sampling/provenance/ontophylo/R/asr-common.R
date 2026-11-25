# Shared validation, reproducibility and disk storage for character inference.

.asr_integer <- function(x, name, minimum = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != floor(x) || x < minimum || x > .Machine$integer.max) {
    stop(name, " must be an integer >= ", minimum, ".", call. = FALSE)
  }
  as.integer(x)
}

.asr_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  x
}

.asr_tree <- function(tree) {
  .check_branch_tree(tree)
  if (!ape::is.rooted(tree) || !ape::is.binary(tree)) {
    stop("Supply a rooted, bifurcating tree; resolve polytomies explicitly.", call. = FALSE)
  }
  if (anyNA(tree$tip.label) || anyDuplicated(tree$tip.label) ||
      any(tree$tip.label == "")) {
    stop("Tree tip labels must be nonempty and unique.", call. = FALSE)
  }
  # corHMM otherwise changes these lengths without changing the caller's tree.
  if (any(tree$edge.length <= 1e-5)) {
    stop("corHMM changes branches <= 1e-5. Rescale tree units explicitly first.",
         call. = FALSE)
  }
  tree <- ape::reorder.phylo(tree, "cladewise")
  tree$node.label <- NULL
  tree$root.edge <- NULL
  tree
}

.asr_names <- function(x, what) {
  if (is.null(x) || anyNA(x) || any(x == "") || anyDuplicated(x)) {
    stop(what, " must be nonempty and unique.", call. = FALSE)
  }
}

.asr_data <- function(data, tree) {
  if (is.atomic(data) && is.null(dim(data))) {
    data <- matrix(as.character(data), ncol = 1L,
                   dimnames = list(names(data), "character_1"))
  }
  if ((!is.matrix(data) && !is.data.frame(data)) || ncol(data) < 1L) {
    stop("data must be a named vector or a taxa-by-character matrix/data frame.",
         call. = FALSE)
  }
  .asr_names(rownames(data), "Taxon names")
  .asr_names(colnames(data), "Character names")
  if (!setequal(rownames(data), tree$tip.label)) {
    stop("Data taxon names must match tree tip labels exactly.", call. = FALSE)
  }
  data <- data[tree$tip.label, , drop = FALSE]
  result <- vector("list", ncol(data))
  names(result) <- colnames(data)
  for (i in seq_along(result)) result[[i]] <- as.character(data[, i])
  result
}

.asr_versions <- function() {
  c(ontophylo = as.character(utils::packageVersion("ontophylo")),
    corHMM = as.character(utils::packageVersion("corHMM")),
    phytools = as.character(utils::packageVersion("phytools")),
    R = as.character(getRversion()), schema = "1")
}

.asr_hash <- function(x) digest::digest(x, algo = "sha256")

.asr_rng_save <- function() {
  seed <- NULL
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    seed <- get(".Random.seed", envir = .GlobalEnv)
  }
  list(kind = RNGkind(), seed = seed)
}

.asr_rng_restore <- function(old) {
  do.call(RNGkind, as.list(old$kind))
  if (is.null(old$seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  } else assign(".Random.seed", old$seed, envir = .GlobalEnv)
}

.asr_rng_start <- function(seed, id, phase) {
  # Deterministic character-specific L'Ecuyer states, independent of job order.
  key <- .asr_hash(list(seed, id, phase))
  state <- integer(6L)
  for (i in seq_len(6L)) {
    state[i] <- strtoi(substr(key, 7L * i - 6L, 7L * i), base = 16L) + 1L
  }
  RNGkind("L'Ecuyer-CMRG")
  assign(".Random.seed", c(10407L, state), envir = .GlobalEnv)
}

.asr_save <- function(object, path) {
  tmp <- tempfile(".pending-", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(object, tmp)
  if (!file.rename(tmp, path)) stop("Could not save ", path, call. = FALSE)
}

.asr_directory <- function(outdir) {
  if (is.null(outdir)) return(NULL)
  if (!is.character(outdir) || length(outdir) != 1L || is.na(outdir) || outdir == "") {
    stop("outdir must be NULL or one directory name.", call. = FALSE)
  }
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  normalizePath(outdir, mustWork = TRUE)
}

.asr_capture <- function(fun, args, timeout = Inf) {
  warnings <- character()
  error <- NULL
  value <- NULL
  # Warnings are retained in the result instead of being discarded.
  record_warning <- function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
  record_error <- function(e) {
    error <<- conditionMessage(e)
    NULL
  }
  output <- utils::capture.output({
    value <- tryCatch(withCallingHandlers(
      R.utils::withTimeout(do.call(fun, args), timeout = timeout,
                           onTimeout = "error"), warning = record_warning),
      error = record_error)
  })
  list(value = value, warnings = unique(warnings), error = error)
}

.asr_init_worker <- function(functions, libpaths) {
  .libPaths(libpaths)
  # Explicit helpers also make PSOCK workers use the current development code.
  for (name in names(functions)) {
    fun <- functions[[name]]
    environment(fun) <- .GlobalEnv
    assign(name, fun, envir = .GlobalEnv)
  }
  NULL
}

.asr_run_jobs <- function(jobs, fun, workers) {
  old_rng <- .asr_rng_save()
  on.exit(.asr_rng_restore(old_rng), add = TRUE)
  workers <- min(workers, length(jobs))
  if (workers <= 1L) {
    result <- vector("list", length(jobs))
    for (i in seq_along(jobs)) result[[i]] <- fun(jobs[[i]])
    return(result)
  }
  helpers <- list()
  env <- environment(.asr_run_jobs)
  helper_names <- ls(env, all.names = TRUE, pattern = "^\\.asr_")
  for (name in helper_names) helpers[[name]] <- get(name, envir = env)
  cluster <- parallel::makePSOCKcluster(workers)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterCall(cluster, .asr_init_worker, helpers, .libPaths())
  # Jobs refer only to names explicitly installed in each worker above.
  parallel::parLapplyLB(cluster, jobs, fun = .asr_worker_call,
                       function_name = deparse(substitute(fun)))
}

.asr_worker_call <- function(job, function_name) {
  fun <- get(function_name, envir = .GlobalEnv)
  fun(job)
}

.asr_check_failures <- function(status, on_error) {
  bad <- status$status != "ok"
  if (any(bad)) {
    message <- paste(sum(bad), "character(s) failed or were skipped:",
                     paste(status$character[bad], collapse = ", "))
    if (on_error == "stop") stop(message, call. = FALSE)
    warning(message, "; inspect $status for details.", call. = FALSE)
  }
}
