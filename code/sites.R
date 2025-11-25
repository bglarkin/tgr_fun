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
    "rnaturalearthhires", "ggspatial", "maps", "ggpubr", "grid", "ggpmisc"
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
    glimpse() # Any reason to not recode field_type to factor here?
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
#' ### Regional map data
#' Retrieve layer for US states and Canadian provinces.
regions <- ne_states(country = c("United States of America", "Canada"),
                     returnclass = "sf")
#' Retrieve metadata for populated places (cities)
cities <- ne_download(
    scale = 50,
    type = "populated_places",
    category = "cultural",
    returnclass = "sf"
)
#' Regional mapping objects
regional_box <- st_bbox(c(
    xmin = -95,
    ymin = 39.25,
    xmax = -81,
    ymax = 49.5
), crs = 4326)
sf_use_s2(FALSE)
regions_crop <- st_crop(regions, regional_box)
cities_crop <- st_crop(cities, regional_box)
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
area_crop <- st_crop(regions, area_box)
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
#' ## Regional map
regional_map <-
    ggplot() +
    geom_sf(
        data = regions_crop,
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
        data = regions_crop %>% filter(
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
        xlim = c(regional_box$xmin, regional_box$xmax),
        ylim = c(regional_box$ymin, regional_box$ymax),
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
#' 
#' ## Final map 
maps_fig <-
    ggarrange(
        regional_map,
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