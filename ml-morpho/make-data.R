library(dplyr)
library(stringr)
library(ape)

#----- Matrix
# Original data
mt <- readRDS("data/adult_matrix.RDS")
tree <-readRDS("data/hym_tree.RDS")
mt
mt$taxa


# Our data
load("data/paramo_stm_adult_final.RDA")
#hym_adult_final$taxa
mt1 <- hym_adult_final
mt1

# get all individual chars from our data
ch <- names(mt1)[-1]
out <- unlist(strsplit(ch, "_"))
out
length(out)


#---------

# Copy the matrix
df <- mt

# select from Original data, only those individual chars used by us
df <- df %>%
  select(taxa, any_of(out))

# # count state freq
# uniq_counts <- df %>%
#   summarise(across(-taxa, ~ {
#     v <- suppressWarnings(as.numeric(.x))
#     length(unique(v[!is.na(v)]))
#   })) %>% unlist
# table(uniq_counts)


# Convert polymorphisms "1&2" into PHYLIP-friendly "?"
# df <- df %>%
#   mutate(across(-taxa, ~ str_replace_all(.x, "([0-9])&([0-9])", "?")))
df <- df %>%
  mutate(across(-taxa, ~ str_replace_all(.x, "\\b[0-9]+(&[0-9]+)+\\b", "?")))


# Replace "-" with "?"
df <- df %>%
  mutate(across(-taxa, ~ str_replace_all(.x, "-", "?")))

# filter out constant chars
df <- df %>%
  select(
    taxa,
    where(function(x) {
      vals <- suppressWarnings(as.numeric(x))
      vals <- vals[!is.na(vals)]
      length(unique(vals)) > 1
    })
  )


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



