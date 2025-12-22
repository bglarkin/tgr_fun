#' ---
#' title: "Site locations and pairwise distances"
#' author: "Beau Larkin\n"
#' date: "Last updated: `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   github_document:
#'     toc: true
#'     toc_depth: 2
#' ---
#'
#' # Description
#' 
#' - Calculate inter-site pairwise distances and summarize
#' - Produce continental and area views in a two-panel map
#' 
#' # Package and library installation
packages_needed = c(
    "tidyverse", "colorspace", "knitr", "conflicted", "grid", "gridExtra",
    "geosphere", "ggrepel", "sf", "rnaturalearth","rnaturalearthdata", 
    "rnaturalearthhires", "ggspatial", "maps", "ggpubr", "ggpmisc",
    "osmdata", "cowplot"
)

to_install <- setdiff(packages_needed, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)
invisible(lapply(packages_needed, library, character.only = TRUE))

#' ## Root path function
root_path <- function(...) rprojroot::find_rstudio_root_file(...)

#+ conflicts,message=FALSE
conflicts_prefer(
  dplyr::filter,
  dplyr::select,
  purrr::map,
  ggpp::annotate
)
#' 
#+ graphics_styles
source(root_path("resources", "styles.R"))
#' 

#' 
#' # Functions
#' Executed from a separate script to save lines here; to view the function navigate to 
#' `functions.R` in the code folder, accessible from the root dir of the repo.
# Functions ———————— ####
source(root_path("code", "functions.R"))

#' # Sites
sites <-
  read_csv(root_path("clean_data/sites.csv"), show_col_types = FALSE) %>% 
  mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant")))
#' Calculate region locations
region_locs <- sites %>%
  group_by(region) %>%
  summarize(long_cen = mean(long),
              lat_cen = mean(lat),
              .groups = "drop")
#' ## Site types
#' How many sites are in each field type?
field_types <- sites %>% 
  group_by(region, field_type) %>% 
  summarize(n = n(), yr_min = min(yr_since), yr_max = max(yr_since), .groups = "drop") %>% 
  pivot_wider(names_from = field_type, values_from = n:yr_max) %>% 
  select(region, starts_with("n"), yr_min_restored, yr_max_restored)
kable(field_types,
  format = "pandoc",
  caption = "Count of sites by type in each area & age of restored fields:\nBM = Blue Mounds, FG = Faville Grove,\nFL = Fermilab, LP = Lake Petite")
#'
#' ## Site pairwise distances
#' Convert site coordinates to distances using the haversine method (package geosphere). 
field_dist <- as.dist(distm(sites[, c("long", "lat")], fun = distHaversine))
#' Stats, pairwise distances in regions
rbind(
    split(sites, sites$region) %>% 
        map(\(df) as.dist(distm(df[, c("long", "lat")], fun = distHaversine))) %>% 
        map(
            \(df) data.frame(
                minimum = min(df),
                median = median(df),
                maximum = max(df)
            ) %>% map( ~ round(. / 1000, 1))
        ) %>% 
        bind_rows(.id = "region"),
    data.frame(
        minimum = min(field_dist),
        median = median(field_dist),
        maximum = max(field_dist)
    ) %>% map( ~ round(. / 1000, 1)) %>% 
        bind_rows() %>% 
        mutate(region = "All") %>% 
        select(region, everything())
) %>% 
    kable(format = "pandoc", caption = "Summary of intra-region and overall pairwise distances (km)")
#' Use function `reg_dist_stats()` to produce pairwise summaries in regions and with field types
reg_ft_stats <- 
  list(
    BM = reg_dist_stats(field_dist, sites, "BM"),
    FG = reg_dist_stats(field_dist, sites, "FG"),
    FL = reg_dist_stats(field_dist, sites, "FL"),
    LP = reg_dist_stats(field_dist, sites, "LP")
  ) %>% 
  bind_rows(.id = "region") %>% 
  arrange(median_dist)
kable(reg_ft_stats, format = "pandoc", caption = "Summary of intra-region and field type pairwise distances (km)")

ft_stats <- 
  reg_ft_stats %>% 
  group_by(group_pair) %>% 
  summarize(median_dist = round(median(median_dist), 1), .groups = "drop") %>% 
  arrange(median_dist)
kable(ft_stats, format = "pandoc", caption = "Summary of field type pairwise distances")

#' 
#' # Maps
#' ## Map data
#' Create two panel map with regional and area views.
#' ### Continental map data
#' Retrieve layer for US states and Canadian provinces.
cont <- ne_states(country = c("United States of America", "Canada"),
                  returnclass = "sf")
#' Retrieve metadata for populated places (cities)
cities <- ne_download( 
    scale = 50,
    type = "populated_places",
    category = "cultural",
    returnclass = "sf"
)
#' Continental mapping objects
cont_box <- st_bbox(c(
    xmin = -95,
    ymin = 39.25,
    xmax = -81,
    ymax = 49.5
), crs = 4326)
sf_use_s2(FALSE)
#+ cont_crop,message=FALSE,warning=FALSE
cont_crop <- st_crop(cont, cont_box)
#+ cities_crop,message=FALSE,warning=FALSE
cities_crop <- st_crop(cities, cont_box)
#' 
#' ### Area map data
area_cities <- ne_download(
    scale = 10,
    type = "populated_places",
    category = "cultural",
    returnclass = "sf"
)
counties <- st_as_sf(maps::map("county", plot = FALSE, fill = TRUE))
#' Area mapping objects
area_box <- st_bbox(c(
    xmin = -90.3,
    ymin = 41.4,
    xmax = -87.4,
    ymax = 43.6
), crs = 4326)
sf_use_s2(FALSE)
#+ area_crop,message=FALSE,warning=FALSE
area_crop <- st_crop(cont, area_box)
#+ area_cities_crop,message=FALSE,warning=FALSE
area_cities_crop <- st_crop(area_cities, area_box)
#+ counties_crop,message=FALSE,warning=FALSE
counties_crop <- st_crop(st_transform(counties, 4326), area_box)
#+ area_roads,message=FALSE,warning=FALSE
# area_roads <- get_osm_roads(area_box, density = 2) # leave commented unless missing from env
#' 
#' ### Site map data
sites_sf <- st_as_sf(sites, coords = c("long", "lat"), crs = 4326, remove = FALSE)
#' Per-region bboxes with user-defined buffer distance from map border. Uses 
#' function `bbox_buffer_km()`.
bb_BM <- bbox_buffer_km(sites_sf %>% filter(region == "BM"), buffer_km = 5)
bb_FG <- bbox_buffer_km(sites_sf %>% filter(region == "FG"), buffer_km = 0.3)
bb_FL <- bbox_buffer_km(sites_sf %>% filter(region == "FL"), buffer_km = 1)
bb_LP <- bbox_buffer_km(sites_sf %>% filter(region == "LP"), buffer_km = 0.2)
#' Retrieve roads data
# Don't execute if roads data are in the local env to save time
# rd_BM = get_osm_roads(bb_BM, density=4)
# rd_FG = get_osm_roads(bb_FG, density=8)
# rd_FL = get_osm_roads(bb_FL, density=8)
# rd_LP = get_osm_roads(bb_LP, density=8)
# site_roads <- list(rd_BM = rd_BM, rd_FG = rd_FG, rd_FL= rd_FL, rd_LP = rd_LP)
#' 
#' 
#' ### Map styles
state_lab_size <- 2.4
state_lab_col <- "darkslateblue"
city_lab_size <- 1.8
city_pt_size <- 1.2
city_col <- "grey35"
panel_lab_x <- 0.02
panel_lab_y <- 0.98

#'
#' ## Continental map
cont_map <-
  ggplot() +
  geom_sf(
    data = cont_crop,
    fill = "ivory",
    color = "black",
    size = 0.3
  ) +
  geom_rect(
    aes(
      xmin = -90.3,
      ymin = 41.5,
      xmax = -87.4,
      ymax = 43.45
    ),
    fill = NA,
    color = "indianred",
    linewidth = 0.6
  ) +
  geom_text_repel(
    data = cont_crop %>% filter(name %in% c("Minnesota", "Iowa", "Ontario", "Michigan",
                                            "Ohio", "Illinois", "Wisconsin","Indiana")),
    aes(x = longitude, y = latitude, label = name),
    size = state_lab_size,
    color = state_lab_col,
    segment.color = NA
  ) +
  geom_point(
    data = cities_crop,
    aes(x = LONGITUDE, y = LATITUDE),
    color = city_col,
    size = city_pt_size
  ) +
  geom_text_repel(
    data = cities_crop,
    aes(x = LONGITUDE, y = LATITUDE, label = NAME),
    size = city_lab_size,
    color = city_col,
    nudge_y = 0
  ) +
  annotation_scale(
    location = "bl",
    width_hint = 0.4,
    height = unit(0.15, "cm")
  ) +
  geom_label_npc(
    aes(npcx = panel_lab_x, npcy = panel_lab_y, label = "A"),
    hjust = "left",
    vjust = "top",
    size = 3,
    fontface = "bold",
    label.r = unit(0.3, "mm"),
    label.size = 0.4,
    label.padding = unit(c(0.4, 0.3, 0.15, 0.3), "lines")
  ) +
  coord_sf(
    xlim = c(cont_box$xmin, cont_box$xmax),
    ylim = c(cont_box$ymin, cont_box$ymax),
    expand = FALSE
  ) +
  theme_void() +
  theme(
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = "aliceblue", color = "black", linewidth = 0.5)
  )

#'
#' ## Area map
area_map <-
  ggplot() +
  geom_sf(
    data = counties,
    fill = "ivory",
    color = "gray80") +
  geom_sf(
    data = area_crop,
    fill = NA,
    color = "black",
    size = 0.3
  ) +
  geom_sf(
    data = area_roads,
    color = "grey70",
    linewidth = 0.3
  ) +
  geom_point(
    data = region_locs,
    aes(x = long_cen, y = lat_cen),
    color = "indianred",
    size = 3
  ) +
  geom_label_repel(
    data = region_locs,
    aes(x = long_cen, y = lat_cen, label = region),
    color = "indianred",
    fill = "snow",
    size = 3,
    fontface = "bold",
    label.r = unit(0.3, "mm"),
    label.size = 0.4,
    label.padding = unit(1.3, "mm"),
    nudge_x = c(-0.25, 0.25, -0.25, 0.10),
    nudge_y = c(0.06, -0.05, 0.02, 0.15)
  ) +
  geom_text_repel(
    data = area_crop %>% filter(name %in% c("Wisconsin", "Illinois")),
    aes(x = longitude, y = latitude, label = name),
    size = state_lab_size,
    color = state_lab_col,
    nudge_x = c(0.8, 0.1),
    nudge_y = c(-1, 1.8),
    segment.color = NA
  ) +
  geom_point(
    data = area_cities_crop,
    aes(x = LONGITUDE, y = LATITUDE),
    color = city_col,
    size = city_pt_size
  ) +
  geom_text_repel(
    data = area_cities_crop,
    aes(x = LONGITUDE, y = LATITUDE, label = NAME),
    size = city_lab_size,
    color = city_col,
    nudge_y = 0
  ) +
  annotation_scale(
    location = "bl",
    width_hint = 0.4,
    height = unit(0.15, "cm")
  ) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    height = unit(0.75, "cm"),
    width = unit(0.75, "cm"),
    pad_x = unit(0.2, "cm"),
    pad_y = unit(0.25, "cm"),
    style = north_arrow_fancy_orienteering()
  ) +
  geom_label_npc(
    aes(npcx = panel_lab_x, npcy = panel_lab_y, label = "B"),
    hjust = "left",
    vjust = "top",
    size = 3,
    fontface = "bold",
    label.r = unit(0.3, "mm"),
    label.size = 0.4,
    label.padding = unit(c(0.4, 0.3, 0.15, 0.3), "lines")
  ) +
  coord_sf(
    xlim = c(area_box$xmin, area_box$xmax),
    ylim = c(area_box$ymin, area_box$ymax),
    expand = FALSE
  ) +
  theme_void() +
  theme(
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = "aliceblue", color = "black", linewidth = 0.5)
  )

#'
#' ## Site maps
#' 
#' Use `nudge_coords` to prevent overlap of points where sites lie close together. 
#' Nudge values are in meters. Positive values nudge in the indicated directions of east
#' or north (e.g., nudge_e_m is east in meters); negative values nudge to the west 
#' or south.
sites_plot <-
  nudge_coords(sites %>%
                 mutate(nudge_e_m = case_when(
                     field_name %in% c("FLRP4") ~ 220,
                     field_name %in% c("FLRSP3") ~ -140,
                     field_name %in% c("FLRP5") ~ -180,
                     field_name %in% c("FLREM1") ~ -60,
                     TRUE ~ 0),
                   nudge_n_m = case_when(
                     field_name %in% c("MBRP1", "PHRP1", "MHRP2") ~ -1900,
                     field_name %in% c("FLRSP1") ~ -140,
                     field_name %in% c("FLREM1") ~ 60,
                     TRUE ~ 0)
                 ))
#' 
#' Produce individual region panels using `make_zoom_map()`. 
#+ map_BM,message=FALSE,warning=FALSE
map_BM <- make_zoom_map(bb_BM, panel_tag = "BM", road_data = site_roads$rd_BM)
#+ map_FG,message=FALSE,warning=FALSE
map_FG <- make_zoom_map(bb_FG, panel_tag = "FG", road_data = site_roads$rd_FG)
#+ map_FL,message=FALSE,warning=FALSE
map_FL <- make_zoom_map(bb_FL, panel_tag = "FL", road_data = site_roads$rd_FL)
#+ map_LP,message=FALSE,warning=FALSE
map_LP <- make_zoom_map(bb_LP, panel_tag = "LP", road_data = site_roads$rd_LP)
#' 
#+ sites_grid
region_zoom_grid <- ggarrange(
  map_BM, NULL, map_FG, NULL, map_FL, NULL, map_LP,
  nrow = 1, align = "h", widths = c(rep(c(1, 0.02), 3), 1)
)
#' 
#' ## Legend and citation grobs
legend_plot <-
  ggplot(sites_plot, aes(x = long_plot, y = lat_plot, fill = field_type)) +
  geom_point(
    shape = 21,
    size = sm_size,
    stroke = lw,
    color = "black"
  ) +
  scale_fill_manual(
    values = ft_pal,
    name = "Field type",
    breaks = levels(sites_plot$field_type)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_void() +
  theme(
    legend.position      = "bottom",
    legend.box           = "horizontal",
    legend.key.height    = unit(4, "mm"),
    legend.key.width     = unit(8, "mm"),
    legend.text          = element_text(size = 8),
    legend.title        = element_text(size = 8),
    plot.margin          = margin(0, 0, 0, 0)
  )
legend_grob <- cowplot::get_legend(legend_plot)

credits_grob <- textGrob(
  "Roads: © OpenStreetMap contributors via {osmextract}\nBoundaries & Names: Natural Earth via {rnaturalearth}",
  x = unit(1, "npc") - unit(3, "mm"),
  y = 0.5,
  hjust = 1,
  vjust = 0.5,
  gp = gpar(cex = 0.5, col = "grey20")
)

footer_row <- arrangeGrob(
  grobs   = list(legend_grob, credits_grob),
  ncol    = 2,
  widths  = c(1, 1),
  heights = unit(1, "null")
)
footer_panel <- as_ggplot(footer_row)
#' 
#' ## Final map 
#' Balance row heights to fill white space
fhs <- c(0.56, 0.01, 0.44)
#' Arrange
maps_fig <- ggarrange(
  ggarrange(
    ggarrange(
      cont_map, NULL, area_map, 
      nrow = 1, ncol = 3, widths = c(1, 0.02, 1)
      ),
    NULL,
    region_zoom_grid,
    nrow = 3, heights = fhs
  ),
  footer_panel,
  ncol = 1, 
  heights = c(1, 0.06),
  widths = c(1, 0.8)
)
#+ tgr_map,message=FALSE,fig.height=6.8,fig.width=7
maps_fig
#' 
ggsave(root_path("figs/fig1.svg"), plot = maps_fig, device = "svg",
    width = 6.5, height = (6.5 * fhs[3] / fhs[1]) + 1.2, units = "in",
    dpi = 600)


