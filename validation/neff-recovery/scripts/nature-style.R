# Shared publication style. Source this file before building a figure.
# Functions do not change data, scales, fitted models or global ggplot defaults.

## Typography and physical sizes ----
nature_font <- "Arial"
# Readable is the project default. Change to "submission" only for compact artwork.
nature_profile <- "readable"
nature_base_pt <- 10
nature_tick_pt <- 9
nature_label_pt <- 9
nature_table_pt <- 9
nature_tag_pt <- 11
if (nature_profile == "submission") {
  nature_base_pt <- 7
  nature_tick_pt <- 6.5
  nature_label_pt <- 6
  nature_table_pt <- 7
  nature_tag_pt <- 8
}
nature_text_mm <- nature_label_pt / ggplot2::.pt
nature_axis_mm <- 0.5 / ggplot2::.pt
nature_curve_mm <- 0.75 / ggplot2::.pt
nature_point_mm <- 1.5

## Palette ----
nature_colours <- c(
  blue = "#0072B2", pale_blue = "#E1F1FC", charcoal = "#444444",
  text = "#222222", grey = "#888888", orange = "#E69F00",
  green = "#009E73", purple = "#CC79A7", vermillion = "#D55E00",
  sky_blue = "#56B4E9"
)

## Shared ggplot theme ----
theme_nature <- function() {
  ggplot2::theme_classic(base_size = nature_base_pt, base_family = nature_font) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = nature_colours["text"]),
      axis.title = ggplot2::element_text(size = nature_base_pt),
      axis.text = ggplot2::element_text(size = nature_tick_pt, colour = nature_colours["text"]),
      axis.line = ggplot2::element_line(linewidth = nature_axis_mm, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = nature_axis_mm, colour = "black"),
      axis.ticks.length = grid::unit(1, "mm"),
      legend.title = ggplot2::element_text(size = nature_tick_pt),
      legend.text = ggplot2::element_text(size = nature_tick_pt),
      legend.key.size = grid::unit(3, "mm"),
      legend.background = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = nature_base_pt, face = "plain"),
      plot.title = ggplot2::element_text(size = nature_base_pt, face = "plain"),
      plot.subtitle = ggplot2::element_text(size = nature_tick_pt),
      plot.caption = ggplot2::element_text(size = nature_label_pt),
      plot.tag = ggplot2::element_text(size = nature_tag_pt, face = "bold"),
      plot.margin = ggplot2::margin(4, 4, 4, 4, unit = "mm"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

## Export a main-figure PDF and preview at final size ----
save_nature <- function(plot, filename, width_mm = 183, height_mm = 120,
                        preview = TRUE) {
  if (!grepl("[.]pdf$", filename)) stop("Use a .pdf filename.")
  if (height_mm > 170) warning("Main-figure artwork normally needs height <= 170 mm.")
  ggplot2::ggsave(filename, plot, device = grDevices::cairo_pdf,
                  family = nature_font, width = width_mm, height = height_mm,
                  units = "mm", bg = "white")
  if (preview) {
    preview_file <- sub("[.]pdf$", ".png", filename)
    ggplot2::ggsave(preview_file, plot, device = grDevices::png, type = "cairo",
                    width = width_mm, height = height_mm, units = "mm",
                    dpi = 450, bg = "white")
  }
}
