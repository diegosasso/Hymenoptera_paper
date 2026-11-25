# Run from this project: Rscript --vanilla scripts/01-simulate.R TRIAL
source("scripts/load-packages.R")
source("scripts/tree-simulator.R")
source("scripts/birth-history.R")
arguments <- commandArgs(trailingOnly = TRUE)
trial <- as.integer(arguments[1])
stopifnot(trial >= 1L, trial <= 100L)
seed <- 202611210L + trial * 10L
tree_age <- 300                  # Arbitrary time units, displayed as Ma.
epoch_boundaries <- c(120, 220)   # Forward time: high lambda from 180 to 80 Ma (100 Myr).
lambda_low <- 0.010
lambda_high <- 0.028
lambda <- c(lambda_low, lambda_high, lambda_low)
mu <- 0                         # All epochs; no extinction or pruning.
tip_range <- c(120L, 250L)
n_units <- 10L
characters_per_unit <- 20L
min_descendant_tips <- 6L
min_outside_tips <- 5L          # Tree-design gate only; no character rejection.
overlap_width <- 6               # All positive focal edges share >=6 Ma.
min_high_fraction <- 0.5         # At least half each positive edge is high-lambda.
post_birth_q <- 0.0001         # Fixed q01 = q10 after the forced origin.
character_seed <- 420261120L + trial * 10L
n_maps <- 100L                    # Production mapping depth.
n_starts <- 3L
workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))                   # Serial is portable; use 2+ if local sockets work.
max_tree_attempts <- 1000L
run_dir <- sprintf("simulated/trial-%03d", trial)
if (file.exists(file.path(run_dir, "simulation.rds"))) {
  cat("Keeping simulated trial", trial, "\n")
  quit(save = "no", status = 0)
}
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(diff(epoch_boundaries) >= 100)
set.seed(seed)

# 2. Draw a tree satisfying ONLY prespecified size and branch-design criteria ----
# Expected tip count before conditioning is 2*exp(sum(lambda*epoch_duration)).
# No tree is accepted/rejected using Neff, a correlation, or an inference result.
expected_tips <- 2 * exp(sum(lambda * diff(c(0, epoch_boundaries, tree_age))))
cat("Unconditioned expected tip count:", expected_tips, "\n")
accepted <- FALSE
tree_attempts <- data.frame()
for (attempt in seq_len(max_tree_attempts)) {
  tree <- simulate_epoch_tree(tree_age, epoch_boundaries, lambda, tip_range[2])
  if (is.null(tree)) {
    tree_attempts <- rbind(tree_attempts, data.frame(attempt, tips = NA, reason = "oversized"))
    next
  }
  if (Ntip(tree) < tip_range[1] || min(tree$edge.length) <= 1e-5) {
    tree_attempts <- rbind(tree_attempts, data.frame(attempt, tips = Ntip(tree), reason = "size/short edge"))
    next
  }
  depth <- node.depth.edgelength(tree)
  edge_start <- depth[tree$edge[, 1]]
  edge_end <- depth[tree$edge[, 2]]
  n_edges <- nrow(tree$edge)
  descendants <- vector("list", n_edges)
  n_descendants <- integer(n_edges)
  for (edge in seq_len(n_edges)) {
    child <- tree$edge[edge, 2]
    tips <- child
    if (child > Ntip(tree)) {
      tips <- match(extract.clade(tree, child)$tip.label, tree$tip.label)
    }
    descendants[[edge]] <- tips
    n_descendants[edge] <- length(tips)
  }
  eligible <- which(n_descendants >= min_descendant_tips &
                     n_descendants <= Ntip(tree) - min_outside_tips)
  high_overlap <- pmax(0, pmin(edge_end, epoch_boundaries[2]) -
                           pmax(edge_start, epoch_boundaries[1]))
  high_fraction <- high_overlap / tree$edge.length
  centers <- seq(epoch_boundaries[1] + overlap_width / 2,
                  epoch_boundaries[2] - overlap_width / 2, by = 1)
  candidate_pools <- vector("list", length(centers))
  good_centers <- integer()
  for (i in seq_along(centers)) {
    active <- which(edge_start <= centers[i] & edge_end > centers[i])
    pool <- eligible[edge_start[eligible] <= centers[i] - overlap_width / 2 &
                       edge_end[eligible] >= centers[i] + overlap_width / 2 &
                       high_fraction[eligible] >= min_high_fraction]
    candidate_pools[[i]] <- pool
    if (length(active) >= 20L && length(pool) >= n_units) good_centers <- c(good_centers, i)
  }
  if (!length(good_centers)) {
    tree_attempts <- rbind(tree_attempts, data.frame(attempt, tips = Ntip(tree), reason = "insufficient overlapping edges"))
    next
  }
  chosen_center <- good_centers[sample.int(length(good_centers), 1)]
  focal_positive <- sample(candidate_pools[[chosen_center]], n_units, replace = FALSE)
  anchor_time <- centers[chosen_center]
  tree_attempts <- rbind(tree_attempts, data.frame(attempt, tips = Ntip(tree), reason = "accepted"))
  accepted <- TRUE
  break
}
write.csv(tree_attempts, file.path(run_dir, "tree-attempts.csv"), row.names = FALSE)
if (!accepted) stop("No eligible tree. Inspect tree-attempts.csv; do not condition on correlation.")
stopifnot(is.ultrametric(tree), is.binary(tree),
          sum(edge_start <= anchor_time & edge_end > anchor_time) >= 20L,
          all(edge_start[focal_positive] <= anchor_time - overlap_width / 2),
          all(edge_end[focal_positive] >= anchor_time + overlap_width / 2))
cat("Accepted proposal", attempt, ":", Ntip(tree), "tips; overlap centered at age",
    tree_age - anchor_time, "Ma.\n")

# 3. Known diversification histories and make.bisse.td() ----
# With the requested global middle-epoch pulse, lambda0=lambda1 within epochs.
# Thus this is time-dependent pure birth, a state-independent BiSSE special case.
# The 0/1/0 epoch marker is deterministic, not a separately evolving trait.
# make.bisse.td is a likelihood constructor, NOT a forward-tree simulator.
# All tips are in the final low regime. q01=q10=0 below supplies a constant dummy
# state to the likelihood; the epoch marker is not treated as a stochastic map.
tip_regime <- setNames(rep(0L, Ntip(tree)), tree$tip.label)
bisse_likelihood <- diversitree::make.bisse.td(tree, tip_regime, n.epoch = 3,
                                             strict = FALSE)
# Parameters ordered from present backward: boundaries, then 6 per epoch.
bisse_parameters <- c(rev(tree_age - epoch_boundaries),
                      lambda_low, lambda_low, 0, 0, 0, 0,
                      lambda_high, lambda_high, 0, 0, 0, 0,
                      lambda_low, lambda_low, 0, 0, 0, 0)
names(bisse_parameters) <- diversitree::argnames(bisse_likelihood)
bisse_loglik <- bisse_likelihood(bisse_parameters, root = diversitree::ROOT.GIVEN,
                                root.p = c(1, 0), condition.surv = FALSE)
stopifnot(is.finite(bisse_loglik))
write.csv(data.frame(parameter = names(bisse_parameters), value = bisse_parameters),
          file.path(run_dir, "bisse-parameters.csv"), row.names = FALSE)

# Equal 1-Ma slice midpoints, then 15 complete 20-Ma bins for Pearson tests.
# These are TRUE marginal lambda-mu rates, not node counts or fitted SSE rates.
ages <- seq(tree_age - 0.5, 0.5, by = -1)
forward_times <- tree_age - ages
slice_edges <- vector("list", length(ages))
diversification <- data.frame(age = ages, time_bin = floor(ages / 20),
                               n_lineages = NA_integer_, regime = NA_integer_, rate = NA_real_)
for (i in seq_along(ages)) {
  slice_edges[[i]] <- which(edge_start <= forward_times[i] & edge_end > forward_times[i])
  regime <- as.integer(forward_times[i] >= epoch_boundaries[1] &
                         forward_times[i] < epoch_boundaries[2])
  state_each_lineage <- rep(regime, length(slice_edges[[i]]))
  lambda_each_lineage <- ifelse(state_each_lineage == 1, lambda_high, lambda_low)
  diversification$n_lineages[i] <- length(slice_edges[[i]])
  diversification$regime[i] <- regime
  diversification$rate[i] <- mean(lambda_each_lineage - mu)
}
branch_design <- data.frame(edge = seq_len(n_edges), parent = tree$edge[, 1],
  child = tree$edge[, 2], age_old = tree_age - edge_start, age_young = tree_age - edge_end,
  duration = tree$edge.length, descendant_tips = n_descendants,
  high_epoch_fraction = high_fraction,
  mean_lambda = lambda_low + high_fraction * (lambda_high - lambda_low))
write.tree(tree, file.path(run_dir, "tree.phy"))
write.csv(branch_design, file.path(run_dir, "branch-design.csv"), row.names = FALSE)
write.csv(diversification, file.path(run_dir, "true-diversification.csv"), row.names = FALSE)

# 4. Force one origin per character, then let its descendants evolve slowly ----
# Origins occur uniformly within focal-edge portions inside the high-lambda epoch.
# A forced origin changes 0 -> 1. Subsequent flips have q01 = q10 = post_birth_q.
# Evolution continues after the high-lambda epoch ends. All outside tips stay 0.
# Characters are drawn ONCE, including rare constant outcomes; no ascertainment.
set.seed(character_seed)
unit_names <- sprintf("unit_%02d", seq_len(n_units))
groups <- setNames(vector("list", n_units), unit_names)
for (unit in seq_len(n_units)) groups[[unit]] <- sprintf("U%02d_C%02d", unit, seq_len(characters_per_unit))
character_names <- unlist(groups, use.names = FALSE)
observations <- matrix("0", Ntip(tree), length(character_names),
  dimnames = list(tree$tip.label, character_names))
true_maps <- setNames(vector("list", length(character_names)), character_names)
true_counts <- matrix(0L, n_edges, length(character_names), dimnames = list(NULL, character_names))
character_checks <- birth_events <- unit_design <- data.frame()
for (unit in seq_len(n_units)) {
  focal <- focal_positive[unit]
  focal_tips <- descendants[[focal]]
  in_subtree <- logical(n_edges)
  for (edge in seq_len(n_edges)) in_subtree[edge] <- all(descendants[[edge]] %in% focal_tips)
  birth_lower <- max(0, epoch_boundaries[1] - edge_start[focal])
  birth_upper <- min(tree$edge.length[focal], epoch_boundaries[2] - edge_start[focal])
  stopifnot(birth_upper > birth_lower)
  unit_design <- rbind(unit_design, data.frame(unit = unit_names[unit], focal_edge = focal,
    descendants = length(focal_tips), focal_duration = tree$edge.length[focal],
    origin_age_old = tree_age - edge_start[focal] - birth_lower,
    origin_age_young = tree_age - edge_start[focal] - birth_upper,
    expected_extra_changes = post_birth_q * (sum(tree$edge.length[in_subtree]) - mean(c(birth_lower, birth_upper)))))
  for (character in groups[[unit]]) {
    offset <- runif(1, birth_lower + 1e-10, birth_upper - 1e-10)
    draw <- simulate_birth_history(tree, focal, offset, post_birth_q, in_subtree)
    stopifnot(all(draw$tip_states[-focal_tips] == "0"), sum(draw$events$type == "birth") == 1L)
    observations[, character] <- draw$tip_states
    true_maps[[character]] <- list(draw$history)
    true_counts[, character] <- draw$counts
    events <- draw$events
    events$age <- tree_age - edge_start[events$edge] - events$offset
    birth_events <- rbind(birth_events, data.frame(unit = unit_names[unit], character, events))
    exposure <- sum(tree$edge.length[in_subtree]) - offset
    character_checks <- rbind(character_checks, data.frame(scenario = "associated",
      unit = unit_names[unit], character, focal_edge = focal, proposals = 1L,
      birth_age = tree_age - edge_start[focal] - offset,
      exposure_after_birth = exposure, expected_extra_changes = post_birth_q * exposure,
      tips_0 = sum(draw$tip_states == "0"), tips_1 = sum(draw$tip_states == "1"),
      focal_changes = draw$counts[focal], extra_changes = sum(draw$counts) - 1L,
      total_changes = sum(draw$counts)))
  }
}
settings <- list(trial = trial, seed = seed, character_seed = character_seed, lambda = lambda, mu = mu,
  tree_age = tree_age, epoch_boundaries = epoch_boundaries, tip_range = tip_range,
  n_units = n_units, characters_per_unit = characters_per_unit,
  overlap_width = overlap_width, min_high_fraction = min_high_fraction,
  character_model = "forced_birth_slow_er", post_birth_q = post_birth_q,
  n_maps = n_maps, n_starts = n_starts, workers = workers)
simulation <- list(tree = tree, branch_design = branch_design, diversification = diversification,
  anchor_age = tree_age - anchor_time, focal_edges = list(associated = focal_positive),
  groups = groups, bisse_parameters = bisse_parameters, bisse_loglik = bisse_loglik,
  scenarios = list(associated = list(observations = observations, true_maps = true_maps,
    true_counts = true_counts, checks = character_checks, events = birth_events,
    unit_design = unit_design)))
saveRDS(settings, file.path(run_dir, "settings.rds"))
saveRDS(simulation, file.path(run_dir, "simulation.rds"))
write.csv(observations, file.path(run_dir, "associated-tip-states.csv"))
write.csv(character_checks, file.path(run_dir, "associated-character-checks.csv"), row.names = FALSE)
write.csv(birth_events, file.path(run_dir, "character-events.csv"), row.names = FALSE)
write.csv(unit_design, file.path(run_dir, "unit-design.csv"), row.names = FALSE)
cat("Origins:", sum(birth_events$type == "birth"), "; additional changes:", sum(character_checks$extra_changes), "\n")
cat("Characters retaining both tip states:", sum(character_checks$tips_0 > 0 & character_checks$tips_1 > 0), "of", length(character_names), "\n")

capture.output(sessionInfo(), file = file.path(run_dir, "session-info.txt"))
writeLines("complete", file.path(run_dir, "simulation-complete.txt"))
