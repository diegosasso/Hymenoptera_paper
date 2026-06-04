# install.packages(c("igraph", "ggraph", "ggplot2"))
library(igraph)
library(ggraph)
library(ggplot2)



# ── palette (matches alluvial_base.R) ─────────────────────────────────────────
group_pal <- c(
  head     = "#4e79a7",
  mesosoma = "#59a14f",
  metasoma = "#e15759",
  leg      = "#f28e2b",
  wing     = "#b07aa1",
  body     = "#2d2d2d"
)

"#76b7b2"   # teal
"#edc948"  # mustard yellow
"#9c755f"  # brown
"#ff9da7"  # pink

# ── edges ──────────────────────────────────────────────────────────────────────
# "metasoma_g" is the internal Level-2 group node; avoids name collision with
# the Level-1 leaf also called "metasoma" (cleaned from "metasoma gen.").
# Edge order controls left-to-right leaf order in the dendrogram.
edges <- rbind(
  data.frame(from = "body",       to = c("head", "mesosoma", "metasoma_g", "leg", "wing")),
  data.frame(from = "head",       to = c("cranium", "mouthparts")),
  data.frame(from = "mesosoma",   to = c("pronotum", "propectus", "mesonotum",
                                          "mesopectus", "metanotum", "mp. complex")),
  data.frame(from = "metasoma_g", to = c("metasoma", "genitalia")),
  data.frame(from = "leg",        to = c("fore leg", "mid leg", "hind leg")),
  data.frame(from = "wing",       to = c("fore wing", "hind wing")),
  # Level-1 → Level-0: each region connects to its own bottom copy
  data.frame(
    from = c("cranium", "mouthparts",
             "pronotum", "propectus", "mesonotum", "mesopectus", "metanotum", "mp. complex",
             "metasoma", "genitalia",
             "fore leg", "mid leg", "hind leg",
             "fore wing", "hind wing"),
    to   = paste0("lf_", c("cranium", "mouthparts",
                            "pronotum", "propectus", "mesonotum", "mesopectus", "metanotum", "mp. complex",
                            "metasoma", "genitalia",
                            "fore leg", "mid leg", "hind leg",
                            "fore wing", "hind wing"))
  )
)

# ── nodes ──────────────────────────────────────────────────────────────────────
nodes <- data.frame(
  name = c(
    "body",
    "head",      "mesosoma",    "metasoma_g",       "leg",            "wing",
    "cranium",   "mouthparts",
    "pronotum",  "propectus",   "mesonotum",        "mesopectus",
    "metanotum", "mp. complex",
    "metasoma",  "genitalia",
    "fore leg",  "mid leg",     "hind leg",
    "fore wing", "hind wing",
    # Level-0 leaves (lf_ prefix for unique igraph names, same display labels)
    paste0("lf_", c("cranium", "mouthparts",
                    "pronotum", "propectus", "mesonotum", "mesopectus", "metanotum", "mp. complex",
                    "metasoma", "genitalia",
                    "fore leg", "mid leg", "hind leg",
                    "fore wing", "hind wing"))
  ),
  label = c(
    "body",
    "head",      "mesosoma",    "metasoma",         "leg",            "wing",
    "cranium",   "mouthparts",
    "pronotum",  "propectus",   "mesonotum",        "mesopectus",
    "metanotum", "mp. complex",
    "metasoma",  "genitalia",
    "fore leg",  "mid leg",     "hind leg",
    "fore wing", "hind wing",
    # same display labels as Level-1
    c("cranium", "mouthparts",
      "pronotum", "propectus", "mesonotum", "mesopectus", "metanotum", "mp. complex",
      "metasoma", "genitalia",
      "fore leg", "mid leg", "hind leg",
      "fore wing", "hind wing")
  ),
  group = c(
    "body",
    "head",      "mesosoma",    "metasoma",         "leg",            "wing",
    "head",      "head",
    "mesosoma",  "mesosoma",    "mesosoma",         "mesosoma",
    "mesosoma",  "mesosoma",
    "metasoma",  "metasoma",
    "leg",       "leg",         "leg",
    "wing",      "wing",
    # Level-0 inherits same groups as Level-1
    c("head",     "head",
      "mesosoma", "mesosoma", "mesosoma", "mesosoma", "mesosoma", "mesosoma",
      "metasoma", "metasoma",
      "leg",      "leg",      "leg",
      "wing",     "wing")
  ),
  stringsAsFactors = FALSE
)

# ── graph ──────────────────────────────────────────────────────────────────────
g <- graph_from_data_frame(edges, vertices = nodes, directed = TRUE)

# ══ OPTIONS ═══════════════════════════════════════════════════════════════════
# Switch between node styles:  "circle"  or  "rect"
node_style <- "rect"

# Circle options
circle_size_internal <- 2   # dot diameter for root + Level-2 nodes
circle_size_leaf     <- 2    # dot diameter for leaf nodes
circle_label_size    <- 1    # text size inside circles

# Rect options
rect_corner_radius   <- 0.3  # unit(r, "lines"); 0 = sharp rectangle
rect_label_size      <- 1  # text size inside rectangles
# ══════════════════════════════════════════════════════════════════════════════

# ── node layers ───────────────────────────────────────────────────────────────
# Labels are always OUTSIDE the nodes in black.
# Internal-node labels appear above; leaf labels appear below at 45°.
if (node_style == "circle") {
  node_geoms <- list(
    geom_node_point(
      data = function(x) x[!x$leaf, ],
      aes(fill = group), shape = 21, size = circle_size_internal,
      colour = "white", stroke = .5
    ),
    geom_node_point(
      data = function(x) x[x$leaf, ],
      aes(fill = group), shape = 21, size = circle_size_leaf,
      colour = "white", stroke = .5
    ),
    geom_node_text(
      data = function(x) x[!x$leaf, ],
      aes(label = label), colour = "black",
      size = circle_label_size, fontface = "bold",
      nudge_y = 0.12, vjust = 0, hjust = 0.5
    ),
    geom_node_text(
      data = function(x) x[x$leaf, ],
      aes(label = label), colour = "black",
      size = circle_label_size, angle = 45, hjust = 1,
      nudge_y = -0.08
    )
  )
} else {
  # "rect" uses shape 22 (filled square); increase circle_size_internal/leaf
  # to make them bigger. rect_corner_radius is unused in this mode.
  node_geoms <- list(
    geom_node_point(
      data = function(x) x[!x$leaf, ],
      aes(fill = group), shape = 22, size = circle_size_internal,
      colour = "white", stroke = .5
    ),
    geom_node_point(
      data = function(x) x[x$leaf, ],
      aes(fill = group), shape = 22, size = circle_size_leaf,
      colour = "white", stroke = .5
    ),
    geom_node_text(
      data = function(x) x[!x$leaf, ],
      aes(label = label), colour = "black",
      size = rect_label_size, fontface = "bold",
      nudge_y = 0.12, vjust = 0, hjust = 0.5
    ),
    geom_node_text(
      data = function(x) x[x$leaf, ],
      aes(label = label), colour = "black",
      size = rect_label_size, angle = 45, hjust = 1,
      nudge_y = -0.08
    )
  )
}

# ── plot ───────────────────────────────────────────────────────────────────────
# Vertical layout: root at top, leaves at bottom.
# Edges body→Level2 share one neutral colour; Level2→Level1 edges are coloured
# by the child's anatomical group.
body_edge_col <- "#aaaaaa"   # colour for the 5 root→group branches

p <- ggraph(g, layout = "dendrogram") +
  # root → Level-2 edges: one neutral colour
  geom_edge_diagonal(
    aes(filter = node1.group == "body"),
    colour = body_edge_col, width = .5, alpha = 1
  ) +
  # Level-2 → Level-1 edges: coloured by child's group
  geom_edge_diagonal(
    aes(colour = node2.group, filter = node1.group != "body"),
    width = .5, alpha = 1
  ) +
  node_geoms +
  scale_fill_manual(values = group_pal,  guide = "none") +
  scale_colour_manual(values = group_pal, guide = "none") +
  scale_edge_colour_manual(values = group_pal, guide = "none") +
  coord_cartesian(clip = "off") +   # allow leaf labels to extend below panel
  theme_void() +
  theme(plot.margin = margin(40, 20, 120, 20))  # large bottom margin for labels

print(p)

# ── save ──────────────────────────────────────────────────────────────────────
ggsave("alluvial-plot/anatomy_tree.pdf", p, width = 1.5, height = 5,
        device = cairo_pdf)
# ggsave("../figures/anatomy_tree.svg", p, width = 14, height = 8)
