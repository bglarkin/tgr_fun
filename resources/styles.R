# Graphics styles for TGR repository

# Color palettes
# hcl_wizard()
# hcl_color_picker()
# Field types: corn, restored, remnant:
ft_pal <- c("#EABB38", "#96DC51", "#078F78")

# Shapes and sizes
lw <- 0.4
sm_size <- 3.1
lg_size <- 4.6
yrtx_size <- 2.1

# Ordination style
theme_ord <-
  theme_bw() +
  theme(
    plot.margin = margin(t = 0, r = 0, b = 1, l = 2, unit = "mm"),
    axis.title.x = element_text(
      size = 10, face = 1,
      margin = margin(t = 2, r = 0, b = 0, l = 0, unit = "mm")
    ),
    axis.title.y = element_text(
      size = 10, face = 1,
      margin = margin(t = 0, r = 1, b = 0, l = 0, unit = "mm")
    ),
    axis.text.x = element_text(
      size = 8, face = 1,
      margin = margin(t = 1, r = 0, b = 0, l = 0, unit = "mm")
    ),
    axis.text.y = element_text(
      size = 8, face = 1,
      margin = margin(t = 0, r = 1, b = 0, l = 0, unit = "mm")
    ),
    axis.ticks.length = unit(-1.4, "mm"),
    legend.text = element_text(size = 8, face = 1),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )

# Correlation plots
theme_cor <-
  theme_classic() +
  theme(
    plot.margin = margin(t = 0, r = 0, b = 1, l = 1, unit = "mm"),
    axis.title.x = element_text(
      size = 10, face = 1,
      margin = margin(t = 1, r = 0, b = 0, l = 0, unit = "mm")
    ),
    axis.title.y = element_text(
      size = 10, face = 1,
      margin = margin(t = 0, r = 1, b = 0, l = 0, unit = "mm")
    ),
    axis.text.x = element_text(
      size = 8, face = 1,
      margin = margin(t = 2, r = 0, b = 0, l = 0, unit = "mm")
    ),
    axis.text.y = element_text(
      size = 8, face = 1,
      margin = margin(t = 0, r = 1, b = 0, l = 0, unit = "mm")
    ),
    axis.ticks.length.y = unit(-2, "mm"),
    axis.ticks.length.x = unit(0, "mm"),
    legend.text = element_text(size = 10, face = 1),
    panel.grid = element_blank()
  )

# Correlation plots in facets
theme_corf <-
  theme_bw() +
  theme(
    plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "mm"),
    axis.title.x = element_text(
      size = 12, face = 1,
      margin = margin(t = 2, r = 0, b = 0, l = 0, unit = "mm")
    ),
    axis.title.y = element_text(
      size = 12, face = 1,
      margin = margin(t = 0, r = 2, b = 0, l = 0, unit = "mm")
    ),
    axis.text.x = element_text(
      size = 10, face = 1,
      margin = margin(t = 2, r = 0, b = 0, l = 0, unit = "mm")
    ),
    axis.text.y = element_text(
      size = 10, face = 1,
      margin = margin(t = 0, r = 2, b = 0, l = 0, unit = "mm")
    ),
    axis.ticks.length = unit(-2, "mm"),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = 1)
  )