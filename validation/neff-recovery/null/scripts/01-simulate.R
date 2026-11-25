# Homogeneous root ER control on the exact corresponding associated tree.
# Run from null/: Rscript --vanilla scripts/01-simulate.R TRIAL
source("scripts/load-packages.R")
source("scripts/tree-simulator.R")
arguments <- commandArgs(trailingOnly = TRUE)
trial <- as.integer(arguments[1])
stopifnot(trial >= 1L, trial <= 100L)
manifest <- read.csv("design/task-manifest.csv")
seed_row <- manifest[manifest$trial == trial, ]
source_dir <- sprintf("../simulated/trial-%03d", trial)
run_dir <- sprintf("simulated/trial-%03d", trial)
if (file.exists(file.path(run_dir, "simulation.rds"))) {
  cat("Keeping simulated root ER trial", trial, "\n")
  quit(save = "no", status = 0)
}
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
source_file <- file.path(source_dir, "simulation.rds")
stopifnot(unname(tools::md5sum(source_file)) == seed_row$source_simulation_md5)
associated <- readRDS(source_file)
tree <- associated$tree
groups <- associated$groups
characters <- unlist(groups, use.names = FALSE)
n_edges <- nrow(tree$edge)
homogeneous_q <- 0.003             # Preserved from the original root ER pilot.
root_state <- 0L
character_seed <- seed_row$character_seed
stopifnot(homogeneous_q == seed_row$homogeneous_q)
settings <- list(trial = trial, scenario = "root_er", tree_age = 300,
  lambda = c(.010, .028, .010), mu = 0, homogeneous_q = homogeneous_q,
  root_state = root_state, n_units = 10L, characters_per_unit = 20L,
  character_seed = character_seed, n_maps = 100L, n_starts = 3L,
  source_simulation_md5 = seed_row$source_simulation_md5,
  source_tree_sha256 = seed_row$source_tree_sha256)

# Every character is drawn once. No focal branches, clade restriction, forced
# origins, tip-balance filtering, or rejection of constant terminal states.
set.seed(character_seed)
q01 <- q10 <- rep(homogeneous_q, n_edges)
observations <- matrix("0", Ntip(tree), length(characters),
  dimnames = list(tree$tip.label, characters))
true_maps <- setNames(vector("list", length(characters)), characters)
true_counts <- matrix(0L, n_edges, length(characters),
  dimnames = list(NULL, characters))
checks <- events <- data.frame()
for (unit in names(groups)) {
  for (character in groups[[unit]]) {
    draw <- simulate_binary_history(tree, q01, q10)
    observations[, character] <- draw$tip_states
    true_maps[[character]] <- list(draw$history)
    true_counts[, character] <- draw$counts
    checks <- rbind(checks, data.frame(scenario = "root_er", unit, character,
      proposals = 1L, tips_0 = sum(draw$tip_states == "0"),
      tips_1 = sum(draw$tip_states == "1"), total_changes = sum(draw$counts)))
    for (edge in which(draw$counts > 0L)) {
      dwell <- draw$history$maps[[edge]]
      event_index <- seq_len(length(dwell) - 1L)
      offset <- cumsum(dwell)[event_index]
      events <- rbind(events, data.frame(unit, character, edge, offset,
        from = names(dwell)[event_index], to = names(dwell)[event_index + 1L],
        type = "ER_transition", age = associated$branch_design$age_old[edge] - offset))
    }
  }
}
Q <- matrix(c(-homogeneous_q, homogeneous_q, homogeneous_q, -homogeneous_q),
  2, 2, dimnames = list(c("0", "1"), c("0", "1")))
simulation <- list(tree = tree, groups = groups,
  branch_design = associated$branch_design,
  diversification = associated$diversification,
  focal_edges = list(root_er = integer()),
  scenarios = list(root_er = list(observations = observations,
    true_maps = true_maps, true_counts = true_counts, checks = checks,
    events = events, Q = Q, root_state = root_state)))
saveRDS(settings, file.path(run_dir, "settings.rds"))
write.csv(checks, file.path(run_dir, "root-er-character-checks.csv"), row.names = FALSE)
write.csv(observations, file.path(run_dir, "root-er-tip-states.csv"))
write.csv(events, file.path(run_dir, "character-events.csv"), row.names = FALSE)
write.csv(simulation$branch_design, file.path(run_dir, "branch-design.csv"), row.names = FALSE)
write.csv(simulation$diversification, file.path(run_dir, "true-diversification.csv"), row.names = FALSE)
file.copy(file.path(source_dir, "tree-attempts.csv"), file.path(run_dir, "tree-attempts.csv"))
file.copy(file.path(source_dir, "tree.phy"), file.path(run_dir, "tree.phy"))
capture.output(sessionInfo(), file = file.path(run_dir, "session-info.txt"))
saveRDS(simulation, file.path(run_dir, "simulation.rds.tmp"))
stopifnot(file.rename(file.path(run_dir, "simulation.rds.tmp"), file.path(run_dir, "simulation.rds")))
writeLines("complete", file.path(run_dir, "simulation-complete.txt"))
cat("Root ER trial", trial, ":", nrow(checks), "characters;", sum(true_counts), "changes.\n")
