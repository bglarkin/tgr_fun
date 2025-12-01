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
    "tidyverse", "colorspace", "knitr", "conflicted", "gridExtra",
    "geosphere", "ggrepel", "sf", "rnaturalearth","rnaturalearthdata", 
    "rnaturalearthhires", "ggspatial", "maps", "ggpubr", "grid", "ggpmisc",
    "osmdata", "digest"
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
cont_crop <- st_crop(cont, cont_box)
cities_crop <- st_crop(cities, cont_box)
#' 
#' ### Area map data
#' Retrieve populated places (cities)
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
area_crop <- st_crop(cont, area_box)
area_cities_crop <- st_crop(area_cities, area_box)
counties_crop <- st_crop(st_transform(counties, 4326), area_box)






#' ### Map styles
state_lab_size <- 2.4
city_lab_size <- 1.8
city_pt_size <- 1.2
city_col <- "grey40"
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
        size = 0.5
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
        data = cont_crop %>% filter(
            name %in% c(
                "Minnesota",
                "Iowa",
                "Ontario",
                "Michigan",
                "Ohio",
                "Illinois",
                "Wisconsin",
                "Indiana"
            )
        ),
        aes(x = longitude, y = latitude, label = name),
        size = state_lab_size,
        alpha = 0.5,
        color = "darkblue",
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
    annotation_scale(location = "bl", 
                     width_hint = 0.4,
                     height = unit(0.15, "cm")) +
    annotation_north_arrow(
        location = "bl",
        which_north = "true",
        height = unit(0.75, "cm"),
        width = unit(0.75, "cm"),
        pad_x = unit(0.25, "cm"),
        pad_y = unit(0.55, "cm"),
        style = north_arrow_fancy_orienteering()
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
    theme(#aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "aliceblue"))

#'
#' ## Area map
area_map <-
    ggplot() +
    geom_sf(
        data = area_crop,
        fill = "ivory",
        color = "black",
        size = 0.5
    ) +
    geom_sf(data = counties,
            fill = NA,
            color = "gray80") +
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
        label.size = 0.6,
        label.padding = unit(1.3, "mm"),
        nudge_x = c(-0.25, 0.25, -0.25, 0.10),
        nudge_y = c(0.06, -0.05, 0.02, 0.15)
    ) +
    geom_text_repel(
        data = area_crop %>%
            filter(name %in% c("Wisconsin", "Illinois")),
        aes(x = longitude, y = latitude, label = name),
        size = state_lab_size,
        alpha = 0.5,
        color = "darkblue",
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
    annotation_scale(location = "bl", 
                     width_hint = 0.4,
                     height = unit(0.15, "cm")) +
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
    theme(#aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "aliceblue"))





# --- Per-region zoom maps (points only; no labels) ---------------------------



# Keep sites as sf (retain long/lat columns for easy aes())
sites_sf <- st_as_sf(sites, coords = c("long", "lat"), crs = 4326, remove = FALSE)

# Buffer a bbox by meters (robust for labeling/jitter padding)
bbox_buffer_km <- function(pts_sf, buffer_km = 20) {
  albers <- 5070 # NAD83 / Conus Albers
  pts_sf %>%
    st_transform(albers) %>%
    st_bbox() %>%
    st_as_sfc(crs = albers) %>%
    st_buffer(dist = buffer_km * 1000) %>%
    st_transform(4326) %>%
    st_bbox()
}




# Convert per-row meter nudges to degrees and add *_plot columns
# Expect these columns in `sites` (you add them manually where needed):
#   nudge_e_m = meters to move EAST  (positive = east, negative = west)
#   nudge_n_m = meters to move NORTH (positive = north, negative = south)
nudge_coords <- function(sites, lon_col = "long", lat_col = "lat",
                         east_col = "nudge_e_m", north_col = "nudge_n_m") {
  stopifnot(all(c(lon_col, lat_col) %in% names(sites)))
  # If nudge columns are missing, treat as zeros
  if (!(east_col  %in% names(sites))) sites[[east_col]]  <- 0
  if (!(north_col %in% names(sites))) sites[[north_col]] <- 0
  
  # Vectorized meters → degrees (WGS84 approximation)
  m_per_deg_lat <- 111320                      # ~ meters per degree latitude
  m_per_deg_lon <- m_per_deg_lat * cospi(sites[[lat_col]] / 180)
  
  dx_deg <- sites[[east_col]]  / m_per_deg_lon
  dy_deg <- sites[[north_col]] / m_per_deg_lat
  
  sites %>%
    mutate(
      long_plot = .data[[lon_col]] + ifelse(is.finite(dx_deg), dx_deg, 0),
      lat_plot  = .data[[lat_col]] + ifelse(is.finite(dy_deg), dy_deg, 0)
    )
}


sites_plot <-
  nudge_coords(
    sites %>%
      mutate(
        nudge_e_m = case_when(
          field_name %in% c("FLRP4") ~  100,   # move east
          field_name %in% c("FLRSP3") ~ -30,   # move west
          field_name %in% c("FLRP5") ~ -80,   # move west
          TRUE ~ 0
        ),
        nudge_n_m = case_when(
          field_name %in% c("MBRP1", "PHRP1", "MHRP2") ~ -1300,   # move south
          field_name %in% c("FLRSP1") ~ -80,   # move south
          TRUE ~ 0
        )
      )
  )



get_osm_roads <- function(bb, density = 8) {
  
  # validate density
  types <- c("motorway","trunk","primary","secondary",
             "tertiary","unclassified","residential","service")
  if (!is.numeric(density) || length(density) != 1L || !is.finite(density)) {
    stop("`density` must be a single finite number in 1:8")
  }
  density <- as.integer(max(1L, min(length(types), density)))
  vals <- types[seq_len(density)]
  
  # build layer
  res <- 
    opq(bbox = bb) %>% 
    add_osm_feature(key = "highway", value = vals) %>% 
    osmdata_sf()
  
  return(st_crop(res$osm_lines, st_as_sfc(bb)))
  
}
  

# Build a zoom map over a bbox; points only, fill by field_type
make_zoom_map <- function(bb,
                          panel_tag = NULL,
                          legend = FALSE,
                          show_counties = TRUE,
                          road_density = 8) {
  
  crop_states   <- st_crop(cont, bb)
  crop_counties <- if (show_counties) st_crop(st_transform(counties, 4326), bb) else NULL
  roads         <- get_osm_roads(bb, density=road_density)
  
  pts <- sites_sf %>%
    filter(long >= bb["xmin"], long <= bb["xmax"],
           lat  >= bb["ymin"], lat  <= bb["ymax"]) %>%
    st_drop_geometry()
  
  # Labels only where yr_since is available (restored sites)
  pts_lab <- sites_plot %>%
    filter(!is.na(yr_since)) %>%
    mutate(lbl = as.character(round(yr_since, 0)))
  
  g <- ggplot() +
    geom_sf(data = crop_states, fill = "ivory", color = "black", linewidth = 0.5) +
    { if (!is.null(crop_counties)) geom_sf(data = crop_counties, fill = NA, color = "gray85", linewidth = 0.3) } +
    { if (!is.null(roads)) geom_sf(data = roads, color = "grey60", linewidth = 0.35, alpha = 0.7) } +
    geom_point(
      data = sites_plot,
      aes(x = long_plot, y = lat_plot, fill = field_type),
      shape = 21, size = sm_size, stroke = lw, color = "black"
    ) +
    geom_text(
      data = pts_lab,
      aes(x = long_plot, y = lat_plot, label = lbl),
      size = yrtx_size, family = "sans", fontface = 2, color = "black"
    ) +
    scale_fill_manual(values = ft_pal) +
    annotation_scale(location = "bl", width_hint = 0.35, height = grid::unit(0.15, "cm")) +
    geom_label_npc(
      aes(npcx = 0.02, npcy = 0.98, label = panel_tag %||% ""),
      hjust = "left", vjust = "top",
      size = 3, fontface = "bold",
      label.r = unit(0.3, "mm"),
      label.size = 0.4,
      label.padding = unit(c(0.4, 0.3, 0.15, 0.3), "lines")
    ) +
    coord_sf(
      xlim = c(bb["xmin"], bb["xmax"]),
      ylim = c(bb["ymin"], bb["ymax"]),
      expand = FALSE
    ) +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "aliceblue"),
      legend.position  = if (legend) "right" else "none"
    )
  
  g
}

# Per-region bboxes (tweak buffer_km to taste)
bb_BM <- bbox_buffer_km(sites_sf %>% filter(region == "BM"), buffer_km = 5)
bb_FG <- bbox_buffer_km(sites_sf %>% filter(region == "FG"), buffer_km = 0.3)
bb_FL <- bbox_buffer_km(sites_sf %>% filter(region == "FL"), buffer_km = 1)
bb_LP <- bbox_buffer_km(sites_sf %>% filter(region == "LP"), buffer_km = 0.2)

# Build four zoom panels (letters continue A/B from your top row)
map_BM <- make_zoom_map(bb_BM, panel_tag = "BM", show_counties = FALSE, road_density = 4)
map_FG <- make_zoom_map(bb_FG, panel_tag = "FG", show_counties = FALSE)
map_FL <- make_zoom_map(bb_FL, panel_tag = "FL", show_counties = FALSE, road_density = 7)
map_LP <- make_zoom_map(bb_LP, panel_tag = "LP", show_counties = FALSE)

# Assemble the 2x2 zoom grid
region_zoom_grid <- ggarrange(
  map_BM, NULL, map_FG, NULL, map_FL, NULL, map_LP,
  nrow = 1, align = "h", widths = c(rep(c(1, 0.02), 3), 1)
)

# Combine with your existing two-panel (A/B) figure
maps_fig_all <- ggarrange(
  ggarrange(cont_map, NULL, area_map, nrow = 1, ncol = 3, widths = c(1, 0.02, 1)),
  region_zoom_grid,
  nrow = 2, heights = c(0.57,0.43)
)

maps_fig_all

# Save if desired
ggsave(
  root_path("figs/fig1_plus_regions.png"),
  maps_fig_all, width = 6.5, height = 6.6, units = "in", dpi = 600
)











#' 
#' ## Final map 
maps_fig <-
    ggarrange(
        cont_map,
        NULL,
        area_map,
        nrow = 1,
        ncol = 3,
        align = "h",
        widths = c(1, 0.02, 1)
    )
#+ tgr_map,message=FALSE,fig.height=3.25,fig.width=6.5,fig.align='center'
maps_fig
#' 
ggsave(
    root_path("figs/fig1.png"),
    width = 6.5,
    height = 3.25,
    units = "in",
    dpi = 600
)









# 1) Confirm your bbox is sane (numbers in degrees, xmin < xmax, etc.)
bb_FL

# 2) Try a single pull verbosely
rtest <- get_roads_osmextract(bb_FL, density = 8,
                              extracts = c("north-america/us/illinois"),
                              cache_dir = root_path("cache","osm"),
                              quiet = FALSE)
rtest; if (!is.null(rtest)) plot(sf::st_geometry(rtest))
