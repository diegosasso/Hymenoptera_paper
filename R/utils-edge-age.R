# Returns a table with:
# - edge_id: numbering used by edgelabels() on a plotted tree
# - avg_edge_age: average of parent and child node ages (present = 0)
edge_average_age_table <- function(tree, use_plot_order = TRUE) {
  if (!inherits(tree, "phylo")) stop("'tree' must be an object of class 'phylo'.")
  if (is.null(tree$edge.length)) stop("'tree' must have edge lengths (a time tree).")
  
  n_tip <- length(tree$tip.label)
  bt <- branching.times(tree)
  
  node_age <- function(node_id) {
    if (node_id <= n_tip) return(0)
    age <- unname(bt[as.character(node_id)])
    if (is.na(age)) stop("Could not find age for internal node: ", node_id)
    age
  }
  
  if (use_plot_order) {
    # Build the same edge order used by edgelabels() by plotting once.
    plot(tree, show.tip.label = FALSE)
    lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    edge_mat <- lp$edge
    edge_id <- seq_len(nrow(edge_mat))
  } else {
    edge_mat <- tree$edge
    edge_id <- seq_len(nrow(edge_mat))
  }
  
  parent_age <- vapply(edge_mat[, 1], node_age, numeric(1))
  child_age  <- vapply(edge_mat[, 2], node_age, numeric(1))
  
  tibble(
    edge_id = edge_id,
    avg_edge_age = (parent_age + child_age) / 2
  )
}

# # Example: table for your main time tree (edge IDs match edgelabels numbering)
# edge_age_tbl <- edge_average_age_table(tree)
# edge_age_tbl