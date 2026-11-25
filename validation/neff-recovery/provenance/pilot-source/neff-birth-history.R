# One forced origin, followed by slow ER evolution in its descendant lineage.
# All ancestors and branches outside the focal clade remain 0.
simulate_birth_history <- function(tree, focal, birth_offset, q, in_subtree) {
  n_tip <- length(tree$tip.label)
  states <- integer(n_tip + tree$Nnode)
  maps <- vector("list", nrow(tree$edge))
  mapped_edge <- matrix(0, nrow(tree$edge), 2, dimnames = list(NULL, c("0", "1")))
  counts <- integer(nrow(tree$edge))
  events <- data.frame(edge = integer(), offset = numeric(), from = integer(),
                       to = integer(), type = character())
  for (edge in seq_len(nrow(tree$edge))) {
    state <- states[tree$edge[edge, 1]]
    elapsed <- 0
    dwell <- numeric()
    labels <- character()
    if (edge == focal) {
      stopifnot(state == 0L, birth_offset > 0, birth_offset < tree$edge.length[edge])
      dwell <- birth_offset
      labels <- "0"
      elapsed <- birth_offset
      state <- 1L
      counts[edge] <- 1L
      events <- rbind(events, data.frame(edge, offset = elapsed, from = 0L, to = 1L, type = "birth"))
    }
    while (elapsed < tree$edge.length[edge]) {
      wait <- Inf
      if (in_subtree[edge] && q > 0) wait <- rexp(1, rate = q)
      remaining <- tree$edge.length[edge] - elapsed
      dwell <- c(dwell, min(wait, remaining))
      labels <- c(labels, as.character(state))
      if (wait >= remaining) break
      elapsed <- elapsed + wait
      previous <- state
      state <- 1L - state
      counts[edge] <- counts[edge] + 1L
      events <- rbind(events, data.frame(edge, offset = elapsed, from = previous,
        to = state, type = "post_birth"))
    }
    states[tree$edge[edge, 2]] <- state
    maps[[edge]] <- setNames(dwell, labels)
    mapped_edge[edge, 1] <- sum(dwell[labels == "0"])
    mapped_edge[edge, 2] <- sum(dwell[labels == "1"])
  }
  history <- tree
  history$maps <- maps
  history$mapped.edge <- mapped_edge
  class(history) <- c("simmap", "phylo")
  list(history = history, counts = counts, events = events,
       tip_states = setNames(as.character(states[seq_len(n_tip)]), tree$tip.label))
}
