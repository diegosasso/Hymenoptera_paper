# Two simulators used by neff-recovery-one-trial.R.
# No inference or file writes happen in these functions.

# Exact pure-birth simulation with piecewise-constant, global birth rates.
# Start with two lineages at a crown split; stop at a fixed age, not tip count.
# max_tips only aborts oversized proposals, which are discarded in full.
simulate_epoch_tree <- function(duration, boundaries, lambda, max_tips) {
  parents <- c(1L, 1L)
  children <- c(2L, 3L)
  starts <- c(0, 0)
  ends <- c(NA_real_, NA_real_)
  alive <- c(1L, 2L)
  next_node <- 4L
  time <- 0
  epoch_ends <- c(boundaries, duration)
  for (epoch in seq_along(lambda)) {
    while (time < epoch_ends[epoch]) {
      waiting <- rexp(1, rate = length(alive) * lambda[epoch])
      if (time + waiting >= epoch_ends[epoch]) {
        time <- epoch_ends[epoch]
        break
      }
      time <- time + waiting
      split_edge <- alive[sample.int(length(alive), 1)]
      ends[split_edge] <- time
      new_edges <- length(parents) + 1:2
      parents <- c(parents, rep(children[split_edge], 2))
      children <- c(children, next_node, next_node + 1L)
      starts <- c(starts, time, time)
      ends <- c(ends, NA_real_, NA_real_)
      next_node <- next_node + 2L
      alive <- c(alive[alive != split_edge], new_edges)
      if (length(alive) > max_tips) return(NULL)
    }
  }
  ends[alive] <- duration
  tip_nodes <- children[alive]
  internal_nodes <- c(1L, children[!seq_along(children) %in% alive])
  node_order <- c(tip_nodes, internal_nodes)
  tree <- list(edge = cbind(match(parents, node_order), match(children, node_order)),
               edge.length = ends - starts, Nnode = length(internal_nodes),
               tip.label = sprintf("tip_%03d", seq_along(tip_nodes)))
  class(tree) <- "phylo"
  storage.mode(tree$edge) <- "integer"
  ape::reorder.phylo(tree, "cladewise")
}

# Exact binary CTMC on every branch, inheriting the parent's end state.
# Q = [-q01, q01; q10, -q10]. Zero rates freeze the character.
# The root is 0. With q=0 outside a focal subtree, those tips remain 0.
simulate_binary_history <- function(tree, q01, q10) {
  n_tip <- length(tree$tip.label)
  states <- integer(n_tip + tree$Nnode)
  branch_maps <- vector("list", nrow(tree$edge))
  mapped_edge <- matrix(0, nrow(tree$edge), 2,
                         dimnames = list(NULL, c("0", "1")))
  counts <- integer(nrow(tree$edge))
  for (edge in seq_len(nrow(tree$edge))) {
    state <- states[tree$edge[edge, 1]]
    time <- 0
    dwell <- numeric()
    labels <- character()
    while (time < tree$edge.length[edge]) {
      rate <- if (state == 0) q01[edge] else q10[edge]
      wait <- Inf
      if (rate > 0) wait <- rexp(1, rate)
      remaining <- tree$edge.length[edge] - time
      dwell <- c(dwell, min(wait, remaining))
      labels <- c(labels, as.character(state))
      if (wait >= remaining) break
      time <- time + wait
      state <- 1L - state
      counts[edge] <- counts[edge] + 1L
    }
    states[tree$edge[edge, 2]] <- state
    branch_maps[[edge]] <- setNames(dwell, labels)
    mapped_edge[edge, 1] <- sum(dwell[labels == "0"])
    mapped_edge[edge, 2] <- sum(dwell[labels == "1"])
  }
  history <- tree
  history$maps <- branch_maps
  history$mapped.edge <- mapped_edge
  class(history) <- c("simmap", "phylo")
  tip_states <- setNames(as.character(states[seq_len(n_tip)]), tree$tip.label)
  list(tip_states = tip_states, history = history, counts = counts)
}
