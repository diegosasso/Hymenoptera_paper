library(dplyr)
library(stringr)
library(ape)

#----- Matrix
mt <- readRDS("data/adult_matrix.RDS")
tree <-readRDS("data/hym_tree.RDS")
mt
mt$taxa

#---------

# Copy the matrix
df <- mt

# Convert polymorphisms "1&2" into PHYLIP-friendly "?"
df <- df %>%
  mutate(across(-taxa, ~ str_replace_all(.x, "([0-9])&([0-9])", "?")))

# Replace "-" with "?"
df <- df %>%
  mutate(across(-taxa, ~ str_replace_all(.x, "-", "?")))

# Build PHYLIP lines
ntax  <- nrow(df)
nchar <- ncol(df) - 1

# PHYLIP requires fixed-width taxon names (10 chars is standard)
df$taxa_fmt <- str_pad(df$taxa, width = 10, side = "right")

# Paste characters for each taxon
df$chars <- apply(df[ , -c(1, ncol(df))], 1, paste0, collapse = "")

# Build final text
phylip_lines <- c(
  paste(ntax, nchar),
  # paste(df$taxa_fmt, df$chars)
  paste(tree$tip.label, df$chars)
)

# Write to file
writeLines(phylip_lines, "ml-morpho/iqtree/morphology.phy")

# --- SAVE tree to phylip
# 1. Remove branch lengths
tree_noBL <- tree
tree_noBL$edge.length <- NULL
# 2. Write to PHYLIP-compatible file
write.tree(tree_noBL, file = "ml-morpho/iqtree/phylo.tre")


#----- RUN IQTREE first to remove invariable sites and then final run

# iqtree -s morphology.phy -m MK+G+ASC -g phylo.tre --prefix ml-morpho
# iqtree -s ml-morpho.varsites.phy -m MK+G+ASC -g phylo.tre --prefix ml-morpho


