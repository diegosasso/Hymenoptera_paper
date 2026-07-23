library(alluvial)

# ── data ──────────────────────────────────────────────────────────────────────
al <- read.csv("alluvial_plot_data.csv",
               header = TRUE, stringsAsFactors = FALSE,
               na.strings = c("", NA))

# clean display labels
al$Level_1[al$Level_1 == "met.-prop. complex"] <- "mp complex"
al$Level_1[al$Level_1 == "metasoma gen."]      <- "metasoma"
al$Level_3 <- "body"

# ── ordering ──────────────────────────────────────────────────────────────────
# alluvial stacks bottom-to-top, so list groups bottom-first → head ends at top
level2_order <- c("wing", "leg", "metasoma", "mesosoma", "head")

level1_order <- c(
  "fore wing", "hind wing",                                                       # wing  (bottom)
  "fore leg",  "mid leg",  "hind leg",                                             # leg
  "metasoma",  "female genitalia",                                                 # metasoma
  "pronotum",  "propectus", "mesonotum", "mesopectus", "metanotum", "mp complex",  # mesosoma
  "cranium",   "mouthparts"                                                        # head  (top)
)

al$Level_1 <- factor(al$Level_1, levels = level1_order)
al$Level_2 <- factor(al$Level_2, levels = level2_order)
al <- al[order(al$Level_1), ]

# ── inter-group gap rows ───────────────────────────────────────────────────────
# Hidden rows inserted after each group (except head).
# They reserve blank space at Level_1 without drawing a visible flow.
gap_freq <- 0.6   # height of gap relative to 1 region; tweak to taste

gap_rows <- data.frame(
  Level_1 = c("__g1__", "__g2__", "__g3__",   "__g4__",   "__g5__"),
  Level_2 = c("wing",   "leg",    "metasoma", "mesosoma", "head"),
  Level_3 = "body",
  stringsAsFactors = FALSE
)

# narrow al to just the 3 axis columns before rbind so gap_rows matches
al3 <- al[, c("Level_1", "Level_2", "Level_3")]

al_plot <- rbind(
  al3[al3$Level_2 == "wing",     , drop = FALSE],  gap_rows[1, ],
  al3[al3$Level_2 == "leg",      , drop = FALSE],  gap_rows[2, ],
  al3[al3$Level_2 == "metasoma", , drop = FALSE],  gap_rows[3, ],
  al3[al3$Level_2 == "mesosoma", , drop = FALSE],  gap_rows[4, ],
  al3[al3$Level_2 == "head",     , drop = FALSE],  gap_rows[5, ]
)

level1_full <- c(
  "fore wing",  "hind wing",  "__g1__",
  "fore leg",   "mid leg",    "hind leg",   "__g2__",
  "metasoma",   "female genitalia",          "__g3__",
  "pronotum",   "propectus",  "mesonotum",  "mesopectus", "metanotum", "mp complex", "__g4__",
  "cranium",    "mouthparts", "__g5__"
)
al_plot$Level_1 <- factor(al_plot$Level_1, levels = level1_full)
al_plot$Level_2 <- factor(al_plot$Level_2, levels = level2_order)

hide_vec <- grepl("^__g[0-9]+__$", as.character(al_plot$Level_1))

# Equal Level_2 blocks: assign freq = target / group_size so every body-part
# block has the same height. Each group also gets one gap row (gap_freq), so
# Level_2 totals are equal across all five groups.
# Trade-off: Level_1 blocks are no longer equal (mesosoma regions are smaller).
group_sizes    <- c(wing = 2, leg = 3, metasoma = 2, mesosoma = 6, head = 2)
target_per_group <- 3   # total real-row freq per group; tweak to taste
freq_real <- target_per_group / group_sizes[as.character(al_plot$Level_2)]
freq_vec  <- ifelse(hide_vec, gap_freq, freq_real)

# ── colours ───────────────────────────────────────────────────────────────────
group_pal <- c(
  "head"     = "#4e79a7",   # slate blue
  "mesosoma" = "#59a14f",   # forest green
  "metasoma" = "#e15759",   # coral red
  "leg"      = "#f28e2b",   # amber
  "wing"     = "#b07aa1"    # muted purple
)
flow_col <- group_pal[as.character(al_plot$Level_2)]
flow_col[is.na(flow_col)] <- "#ffffff"   # gap rows: irrelevant (hidden)

# ── plot ──────────────────────────────────────────────────────────────────────
do_plot <- function() {
  par(mar = c(1, 1, 1, 1))
  alluvial(
    al_plot[, c("Level_1", "Level_2", "Level_3")],
    freq        = freq_vec,
    col         = flow_col,
    border      = flow_col,
    alpha       = 0.75,
    gap.width   = 0.05,
    xw          = 0.25,
    blocks      = TRUE,
    cex         = 0.85,
    hide        = hide_vec,
    axis_labels = c("", "", "")
  )
}

do_plot()

# ── save ──────────────────────────────────────────────────────────────────────
# svg("alluvial_anatomy_base.svg", width = 14, height = 6)
# do_plot()
# dev.off()
