#' Fit candidate Mk models independently to discrete characters
#'
#' Fit selected models with corHMM, retain their comparisons, and prepare a
#' consistent tree, state encoding and root setting for [paramo_stm()].
#'
#' @param tree A rooted, bifurcating `phylo` with branch lengths > 1e-5.
#' @param data Named character vector or taxa-by-character matrix/data frame.
#'   Row names must match tree tips exactly. `?` and `NA` mean unknown; `&`
#'   separates ambiguous states. Inapplicable `-` is rejected by default.
#' @param models Candidate names selected from `"ER"`, `"SYM"`, `"ARD"`, and
#'   names in `custom_models`. Binary ER and SYM are fitted only once.
#' @param custom_models Named list of parameter-index matrices. Rows are source
#'   states, columns are destination states. Positive integers label shared free
#'   rates; zero/NA off-diagonal entries forbid transitions. Dimnames must match
#'   character states. These are not numeric Q estimates.
#' @param model_overrides Named list of candidate-name vectors, by character ID.
#' @param root_p `"equal"` (`"flat"` is an alias), `"yang"` (stationary),
#'   `"maddfitz"`, or a named probability vector. A named list supplies one
#'   setting per character. The same setting is used in fitting and mapping.
#' @param criterion `"AICc"` or `"AIC"`. Undefined AICc values are ineligible;
#'   no automatic switch to AIC is made.
#' @param nstarts Number of corHMM random restarts per candidate. Each character
#'   runs serial restarts, avoiding nested parallelism.
#' @param ancestral `"none"` or `"marginal"`; reconstruct nodes only for the
#'   selected model, with its fitted rates held fixed.
#' @param workers Number of independent character jobs run concurrently.
#' @param seed Nonnegative integer seed. Character-specific streams are stable
#'   under reordering, parallel execution and resumption. Caller RNG is restored.
#' @param outdir Optional checkpoint directory, containing one RDS per character.
#' @param resume Reuse matching checkpoints. Mismatching inputs/settings fail.
#' @param on_error `"continue"` records character failures and warns at the end;
#'   `"stop"` raises an error after collecting jobs and saving completed work.
#' @param control Named list: `lower.bound`, `upper.bound`, `ip`, `opts`,
#'   `lewis.asc.bias`, `timeout`, `keep_fits`, `inapplicable`, `constant`.
#'   Defaults are 1e-9, 100, NULL, NULL, FALSE, Inf, FALSE, `"error"`, `"skip"`.
#'   `inapplicable = "missing"` explicitly treats `-` as unknown.
#'   `constant = "error"` makes constant observations a failed job rather than
#'   a skipped job. Neither option assumes a constant character has zero changes.
#'   Timeout uses R.utils and is cooperative, not a hard native-code deadline.
#' @return An `ontophylo_fit` list: `tree`, named `characters`, `status`,
#'   `settings`, and `versions`. Each successful character includes `model_table`,
#'   `selected_model`, `Q_matrix`, `index_matrix`, `tips_prior`, `root_p`,
#'   `state_labels`, optional `ancestral`, and diagnostics. Optimizer convergence
#'   is marked `"not_reported"` when corHMM does not expose it. Candidate failure
#'   makes a comparison incomplete even if another model can be selected.
#' @details Each input column is fitted independently, not as a joint corHMM
#'   model. Dependent characters must first be explicitly encoded as a composite
#'   character and supplied with an appropriate custom index matrix. This first
#'   implementation uses one observed-state rate category, not hidden rates.
#'   All states must occur in an observation or an ambiguity set. All-missing
#'   and constant characters are reported separately. Missing observations do
#'   not add states to the state space. AICc uses corHMM's tip-count convention.
#'   Large rates and boundary estimates are flagged, never filtered by refitting
#'   until a smaller estimate is obtained. Rates are per tree time unit.
#' @seealso [paramo_stm()], [read_stm()], [branch_rates()]
#' @examples
#' tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
#' x <- c(a = "0", b = "0", c = "1", d = "1")
#' fits <- fit_models(tree, x, models = "ER", root_p = "equal",
#'                    criterion = "AIC", nstarts = 0)
#' fits$characters[[1]]$Q_matrix
#' @export
fit_models <- function(tree, data, models = c("ER", "SYM", "ARD"),
                       custom_models = list(), model_overrides = list(),
                       root_p = "maddfitz", criterion = c("AICc", "AIC"),
                       nstarts = 5L, ancestral = c("none", "marginal"),
                       workers = 1L, seed = 1L, outdir = NULL, resume = FALSE,
                       on_error = c("continue", "stop"), control = list()) {
  tree <- .asr_tree(tree)
  data <- .asr_data(data, tree)
  criterion <- match.arg(criterion)
  ancestral <- match.arg(ancestral)
  on_error <- match.arg(on_error)
  nstarts <- .asr_integer(nstarts, "nstarts", 0L)
  workers <- .asr_integer(workers, "workers")
  seed <- .asr_integer(seed, "seed", 0L)
  resume <- .asr_flag(resume, "resume")
  control <- .asr_fit_control(control)
  if (!is.list(custom_models) || !is.list(model_overrides)) {
    stop("custom_models and model_overrides must be named lists.", call. = FALSE)
  }
  if (length(custom_models)) .asr_names(names(custom_models), "Custom model names")
  if (any(names(custom_models) %in% c("ER", "SYM", "ARD"))) {
    stop("Custom model names must not replace ER, SYM or ARD.", call. = FALSE)
  }
  if (length(model_overrides)) {
    .asr_names(names(model_overrides), "Override names")
    if (!all(names(model_overrides) %in% names(data))) {
      stop("Unknown character in model_overrides.", call. = FALSE)
    }
  }
  if (is.list(root_p)) {
    .asr_names(names(root_p), "Root-setting character names")
    if (!setequal(names(root_p), names(data))) {
      stop("A root_p list must name every character exactly once.", call. = FALSE)
    }
  }
  outdir <- .asr_directory(outdir)
  if (resume && is.null(outdir)) stop("resume requires outdir.", call. = FALSE)
  versions <- .asr_versions()
  jobs <- vector("list", length(data))
  for (i in seq_along(data)) {
    id <- names(data)[i]
    candidates <- models
    if (id %in% names(model_overrides)) candidates <- model_overrides[[id]]
    if (!is.character(candidates) || length(candidates) == 0L ||
        anyNA(candidates) || anyDuplicated(candidates) ||
        !all(candidates %in% c("ER", "SYM", "ARD", names(custom_models)))) {
      stop("Invalid or duplicated candidate models for ", id, ".", call. = FALSE)
    }
    root <- root_p
    if (is.list(root_p)) root <- root_p[[id]]
    job <- list(id = id, tree = tree, x = data[[i]], models = candidates,
                custom_models = custom_models, root_p = root, criterion = criterion,
                nstarts = nstarts, ancestral = ancestral, seed = seed,
                control = control, versions = versions)
    job$fingerprint <- .asr_hash(job)
    job$file <- NULL
    if (!is.null(outdir)) job$file <- file.path(outdir, paste0(.asr_hash(id), ".rds"))
    job$resume <- resume
    jobs[[i]] <- job
  }
  characters <- .asr_run_jobs(jobs, .asr_fit_job, workers)
  names(characters) <- names(data)
  status <- data.frame(character = names(data), status = "", reason = "",
                       stringsAsFactors = FALSE)
  for (i in seq_along(characters)) {
    status$status[i] <- characters[[i]]$status
    status$reason[i] <- characters[[i]]$reason
  }
  result <- list(tree = tree, characters = characters, status = status,
                 settings = list(models = models, criterion = criterion,
                   root_p = root_p, nstarts = nstarts, seed = seed), versions = versions)
  class(result) <- "ontophylo_fit"
  .asr_check_failures(status, on_error)
  incomplete <- status$status == "ok" & status$reason != ""
  if (any(incomplete)) {
    warning("Incomplete model comparisons: ", paste(status$character[incomplete], collapse = ", "),
            "; inspect model_table before interpretation.", call. = FALSE)
  }
  result
}

.asr_fit_control <- function(control) {
  defaults <- list(lower.bound = 1e-9, upper.bound = 100, ip = NULL, opts = NULL,
                    lewis.asc.bias = FALSE, timeout = Inf, keep_fits = FALSE,
                    inapplicable = "error", constant = "skip")
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || anyDuplicated(names(control)) ||
       !all(names(control) %in% names(defaults))))) {
    stop("Unknown or unnamed fit control option.", call. = FALSE)
  }
  for (name in names(control)) defaults[name] <- control[name]
  for (name in c("lower.bound", "upper.bound")) {
    x <- defaults[[name]]
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
      stop(name, " must be finite and positive.", call. = FALSE)
    }
  }
  if (defaults$lower.bound >= defaults$upper.bound) stop("Invalid rate bounds.")
  .asr_flag(defaults$lewis.asc.bias, "lewis.asc.bias")
  .asr_flag(defaults$keep_fits, "keep_fits")
  defaults$inapplicable <- match.arg(defaults$inapplicable, c("error", "missing"))
  defaults$constant <- match.arg(defaults$constant, c("skip", "error"))
  if (!is.numeric(defaults$timeout) || length(defaults$timeout) != 1L ||
      is.na(defaults$timeout) || defaults$timeout <= 0) stop("Invalid timeout.")
  defaults
}

.asr_encode <- function(x, tree, inapplicable) {
  x <- trimws(x)
  x[is.na(x)] <- "?"
  if (any(x == "-")) {
    if (inapplicable == "error") stop("Inapplicable '-' requires explicit recoding.")
    x[x == "-"] <- "?"
  }
  tokens <- strsplit(x, "&", fixed = TRUE)
  states <- character()
  for (i in seq_along(tokens)) {
    if (x[i] == "?") next
    if (x[i] == "" || startsWith(x[i], "&") || endsWith(x[i], "&") ||
        any(tokens[[i]] %in% c("", "?", "-"))) stop("Malformed ambiguity code.")
    states <- union(states, tokens[[i]])
  }
  states <- sort(states)
  if (length(states) < 2L) return(list(states = states))
  k <- length(states)
  # Zero-padded codes prevent corHMM's lexical sort from putting 10 before 2.
  codes <- sprintf(paste0("%0", nchar(k), "d"), seq_len(k))
  tip <- matrix(0, length(x), k, dimnames = list(tree$tip.label, states))
  encoded <- character(length(x))
  for (i in seq_along(x)) {
    if (x[i] == "?") {
      tip[i, ] <- 1
      encoded[i] <- "?"
    } else {
      allowed <- match(tokens[[i]], states)
      tip[i, allowed] <- 1
      encoded[i] <- paste(codes[sort(unique(allowed))], collapse = "&")
    }
  }
  list(states = states, tips = tip, codes = codes,
       data = data.frame(taxon = tree$tip.label, state = encoded))
}

.asr_root <- function(root, states) {
  if (is.character(root) && length(root) == 1L && !is.na(root)) {
    if (root == "flat") root <- "equal"
    if (!root %in% c("equal", "yang", "maddfitz")) stop("Unknown root_p setting.")
    return(root)
  }
  if (!is.numeric(root) || length(root) != length(states) ||
      any(!is.finite(root)) || any(root < 0) || sum(root) <= 0 ||
      !setequal(names(root), states) || anyDuplicated(names(root))) {
    stop("Numeric root_p must give named, nonnegative probabilities for every state.")
  }
  root[states] / sum(root)
}

.asr_index <- function(model, custom_models, states) {
  k <- length(states)
  index <- matrix(0L, k, k, dimnames = list(states, states))
  if (model == "ER") index[row(index) != col(index)] <- 1L
  if (model == "ARD") index[row(index) != col(index)] <- seq_len(k * (k - 1L))
  if (model == "SYM") {
    index[lower.tri(index)] <- seq_len(k * (k - 1L) / 2L)
    index <- index + t(index)
  }
  if (!model %in% c("ER", "SYM", "ARD")) {
    index <- custom_models[[model]]
    if (!is.matrix(index) || !is.numeric(index) ||
        !identical(dim(index), c(k, k)) ||
        !setequal(rownames(index), states) || !setequal(colnames(index), states) ||
        anyDuplicated(rownames(index)) || anyDuplicated(colnames(index))) {
      stop("Custom index matrix must have exactly the character state names.")
    }
    index <- index[states, states, drop = FALSE]
    index[is.na(index)] <- 0
    if (any(!is.finite(index)) || any(index < 0) || any(index != floor(index))) {
      stop("Custom index matrices require nonnegative integer parameter indices.")
    }
    if (any(diag(index) != 0)) stop("Index matrix diagonal must be zero or NA.")
  }
  parameters <- sort(unique(index[index > 0]))
  if (!length(parameters)) stop("Model must contain at least one free rate.")
  # Canonicalize labels, also allowing comparison of equivalent models.
  positive <- which(index > 0)
  index[positive] <- match(index[positive], unique(index[positive]))
  index
}

.asr_fit_job <- function(job) {
  old_rng <- .asr_rng_save()
  on.exit(.asr_rng_restore(old_rng), add = TRUE)
  if (!is.null(job$file) && file.exists(job$file)) {
    saved <- readRDS(job$file)
    if (!job$resume) stop("Checkpoint exists; use resume=TRUE or a new outdir.")
    if (!identical(saved$fingerprint, job$fingerprint)) stop("Fit checkpoint mismatch.")
    if (saved$status %in% c("ok", "skipped")) return(saved)
  }
  result <- tryCatch(.asr_fit_character(job), error = .asr_fit_error)
  result$id <- job$id
  result$fingerprint <- job$fingerprint
  if (!is.null(job$file)) .asr_save(result, job$file)
  result
}

.asr_fit_error <- function(e) list(status = "failed", reason = conditionMessage(e))

.asr_fit_character <- function(job) {
  encoded <- .asr_encode(job$x, job$tree, job$control$inapplicable)
  states <- encoded$states
  if (length(states) < 2L) {
    reason <- "Constant observations: insufficient observed state space for fitting."
    if (!length(states)) reason <- "All observations are missing."
    status <- "skipped"
    if (job$control$constant == "error") status <- "failed"
    return(list(status = status, reason = reason))
  }
  root <- .asr_root(job$root_p, states)
  engine_root <- root
  if (identical(root, "equal")) engine_root <- rep(1 / length(states), length(states))
  n <- length(job$models)
  tab <- data.frame(model = job$models, equivalent_to = NA_character_,
                    logLik = NA_real_, npar = NA_integer_, AIC = NA_real_,
                    AICc = NA_real_, delta = NA_real_, weight = NA_real_,
                    status = "failed", convergence = "not_reported",
                    boundary = FALSE, warnings = "", error = "",
                    stringsAsFactors = FALSE)
  indices <- fits <- vector("list", n)
  for (i in seq_len(n)) {
    model <- job$models[i]
    attempt <- tryCatch(.asr_index(model, job$custom_models, states),
                        error = .asr_fit_error)
    if (is.list(attempt)) {
      tab$error[i] <- attempt$reason
      next
    }
    indices[[i]] <- attempt
    same <- NA_integer_
    if (i > 1L) {
      for (j in seq_len(i - 1L)) {
        if (identical(indices[[i]], indices[[j]])) { same <- j; break }
      }
    }
    if (!is.na(same)) {
      tab$equivalent_to[i] <- job$models[same]
      tab$status[i] <- "equivalent"
      next
    }
    .asr_rng_start(job$seed, job$id, model)
    args <- list(phy = job$tree, data = encoded$data, rate.cat = 1,
                  rate.mat = indices[[i]], node.states = "none",
                  root.p = engine_root, nstarts = job$nstarts, n.cores = 1,
                  collapse = FALSE, lower.bound = job$control$lower.bound,
                  upper.bound = job$control$upper.bound, ip = job$control$ip,
                  opts = job$control$opts, lewis.asc.bias = job$control$lewis.asc.bias)
    attempt <- .asr_capture(corHMM::corHMM, args, job$control$timeout)
    tab$warnings[i] <- paste(attempt$warnings, collapse = "; ")
    if (!is.null(attempt$error)) { tab$error[i] <- attempt$error; next }
    fit <- attempt$value
    if (length(fit$loglik) != 1L || !is.finite(fit$loglik) ||
        !is.matrix(fit$solution) || !identical(dim(fit$solution), dim(indices[[i]])) ||
        any(!is.finite(fit$solution[indices[[i]] > 0]))) {
      tab$error[i] <- "Backend did not return finite likelihood and rates."
      next
    }
    fits[[i]] <- fit
    tab$status[i] <- "ok"
    tab$logLik[i] <- fit$loglik
    tab$npar[i] <- max(indices[[i]])
    tab$AIC[i] <- -2 * fit$loglik + 2 * tab$npar[i]
    denominator <- length(job$tree$tip.label) - tab$npar[i] - 1
    if (denominator > 0) {
      tab$AICc[i] <- tab$AIC[i] + 2 * tab$npar[i] * (tab$npar[i] + 1) / denominator
    }
    rates <- fit$solution[indices[[i]] > 0]
    tab$boundary[i] <- any(rates <= job$control$lower.bound * 1.01 |
                            rates >= job$control$upper.bound * 0.99)
  }
  eligible <- which(tab$status == "ok" & is.finite(tab[[job$criterion]]))
  if (!length(eligible)) {
    return(list(status = "failed", reason = paste("No valid model under", job$criterion),
                 model_table = tab))
  }
  values <- tab[[job$criterion]][eligible]
  tab$delta[eligible] <- values - min(values)
  weights <- exp(-0.5 * tab$delta[eligible])
  tab$weight[eligible] <- weights / sum(weights)
  best <- eligible[which.min(values)]
  fit <- fits[[best]]
  Q <- fit$solution
  Q[indices[[best]] == 0] <- 0
  diag(Q) <- -rowSums(Q)
  dimnames(Q) <- list(states, states)
  node_states <- NULL
  if (job$ancestral == "marginal") {
    rates <- numeric(tab$npar[best])
    for (i in seq_along(rates)) rates[i] <- fit$solution[indices[[best]] == i][1]
    args <- list(phy = job$tree, data = encoded$data, rate.cat = 1,
                  rate.mat = indices[[best]], p = rates, node.states = "marginal",
                  root.p = engine_root, collapse = FALSE,
                  lewis.asc.bias = job$control$lewis.asc.bias)
    attempt <- .asr_capture(corHMM::corHMM, args, job$control$timeout)
    if (!is.null(attempt$error)) stop("Selected-model ASR failed: ", attempt$error)
    node_states <- attempt$value$states
    dimnames(node_states) <- list(length(job$tree$tip.label) + seq_len(job$tree$Nnode), states)
  }
  complete <- !any(tab$status == "failed") &&
    all(is.finite(tab[[job$criterion]][tab$status == "ok"]))
  result <- list(status = "ok", reason = "", model_table = tab,
                   comparison_complete = complete, selected_model = job$models[best],
                   Q_matrix = Q, index_matrix = indices[[best]],
                   tips_prior = encoded$tips, root_p = root, state_labels = states,
                   ancestral = node_states,
                   diagnostics = list(boundary = tab$boundary[best],
                     convergence = "not_reported", lewis.asc.bias = job$control$lewis.asc.bias))
  if (!complete) result$reason <- "Incomplete model comparison; inspect model_table."
  if (job$control$keep_fits) { names(fits) <- job$models; result$fits <- fits }
  result
}
