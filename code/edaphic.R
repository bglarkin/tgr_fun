#' ---
#' title: "Soil properties"
#' author: "Beau Larkin\n"
#' date: "Last updated: `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   github_document:
#'     toc: true
#'     toc_depth: 2
#' ---
#'
#' # Description
#' Soil nutrients were analyzed by [Ward Laboratories, Inc.](https://www.wardlab.com/services/soil-health-analysis/), 
#' analysis methods available in local files or at the link included here.
#' Soil organic matter is in percent determined by the loss-on-ignition method.
#' Soil pH is in a log scale as is typical, and all the other minerals are in parts per million. 
#' This may need to be converted to $mg*kg^{-1}$ or other unit. 
#' 
#' This script provides a quick overview of the soil abiotic property data and tests differences among field types
#' based on soil properties.
#' 
#' # Packages and libraries
packages_needed <- c("tidyverse", "knitr", "vegan", "patchwork", "conflicted", 
                     "permute", "geosphere", "ape")

to_install <- setdiff(packages_needed, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)
invisible(lapply(packages_needed, library, character.only = TRUE))

#' ## Root path function
root_path <- function(...) rprojroot::find_rstudio_root_file(...)

#+ conflicts,message=FALSE
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("diversity", "vegan")
#' 
#+ graphics_styles
source(root_path("resources", "styles.R"))

#+ functions
source(root_path("code", "functions.R"))

#' 
#' #' # Data
#' ## Site metadata and design
sites <- read_csv(root_path("clean_data/sites.csv"), show_col_types = FALSE) %>% 
  mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant")))
#' 
#' ## Soil properties
soil <- read_csv(root_path("clean_data/soil.csv"), show_col_types = FALSE)[-c(26:27), ]
soil_units <- read_csv(root_path("clean_data/soil_units.csv"), show_col_types = FALSE)
#' 
#' # Results
#' ## Quantities in field types
soil_ft_avg <- 
    soil %>% 
    left_join(sites %>% select(field_name, field_type), by = join_by(field_name)) %>% 
    select(-field_key) %>% 
    pivot_longer(pH:Na, names_to = "soil_property", values_to = "qty") %>% 
    group_by(field_type, soil_property) %>% 
    summarize(avg_qty = mean(qty), .groups = "drop") %>% 
    pivot_wider(names_from = "field_type", values_from = "avg_qty") %>% 
    left_join(soil_units, by = join_by(soil_property)) %>% 
    select(soil_property, units, everything()) %>% 
  rowwise() %>% mutate(cv = sd(c_across(corn:remnant)) / mean(c_across(corn:remnant)),
                       across(where(is.numeric), ~ round(.x, 2))) %>% 
  arrange(-cv)

#' ## PCA ordination, variable correlations, and PERMANOVA
soil_z <- decostand(data.frame(soil[, -1], row.names = 1), "standardize")
soil_pca <- rda(soil_z)
summary(soil_pca)
#' Axes 1 and 2 explain 52% of the variation in sites. Axes 1 through 6 account for 91%. 
#' 
#' ## Soil variable loadings and correlations
#' Which soil properties explain the most variation among sites?
site_sco <- scores(soil_pca, display = "sites", choices = c(1,2))
soil_cor <- 
    data.frame(cor(soil_z, site_sco)) %>% 
    mutate(PCA_correlation = sqrt(PC1^2 + PC2^2)) %>% 
    arrange(-PCA_correlation) %>% 
    rownames_to_column(var = "soil_property") %>% 
  mutate(across(where(is.numeric), ~ round(.x, 2)))

#' Use the variable correlations to sort the soil property averages in a table 
#' highlighting field types:
soil_ft_avg %>% 
    left_join(soil_cor %>% select(soil_property, PCA_cor = PCA_correlation), by = join_by(soil_property)) %>% 
  rowwise() %>% 
    mutate(cv = sd(c_across(corn:remnant)) / mean(c_across(corn:remnant)), 
           across(where(is.numeric), ~ round(.x, 2))) %>% 
  arrange(-cv) %>% 
    kable(format = "pandoc")

#' Axis 1 & 2 eigenvalue proportions
eig_prop <- round(summary(soil_pca)$cont$importance[2, 1:2] * 100, 1)
soil_ord_scores <-
    data.frame(site_sco) %>%
    rownames_to_column(var = "field_name") %>%
    left_join(sites, by = join_by(field_name))

#' 
#' ### PERMANOVA on field type
d_soil = dist(soil_z, method = "euclidean")
mva_soil <- soilperm(d = d_soil, env = sites)
#+ soil_ord_results
mva_soil$dispersion_test
mva_soil$permanova
mva_soil$pairwise_contrasts[c(1,3,2), c(1,2,4,3,7,8)] %>% 
  arrange(group1, desc(group2)) %>% 
  kable(format = "pandoc", caption = "Pairwise permanova contrasts")
#' ### Plotting and Fig S2
soil_ord_ft_centers <- soil_ord_scores %>%
  group_by(field_type) %>%
  summarize(across(starts_with("PC"), list(mean = mean, ci_l = ci_l, ci_u = ci_u), .names = "{.fn}_{.col}"), .groups = "drop") %>%
  mutate(across(c(ci_l_PC1, ci_u_PC1), ~ mean_PC1 + .x),
         across(c(ci_l_PC2, ci_u_PC2), ~ mean_PC2 + .x))
soil_ord_ftypes <-
  ggplot(soil_ord_scores, aes(x = PC1, y = PC2)) +
  geom_linerange(data = soil_ord_ft_centers, aes(x = mean_PC1, y = mean_PC2, xmin = ci_l_PC1, xmax = ci_u_PC1), linewidth = lw) +
  geom_linerange(data = soil_ord_ft_centers, aes(x = mean_PC1, y = mean_PC2, ymin = ci_l_PC2, ymax = ci_u_PC2), linewidth = lw) +
  geom_point(data = soil_ord_ft_centers, 
             aes(x = mean_PC1, y = mean_PC2, fill = field_type), 
             size = lg_size, stroke = lw, shape = 21, show.legend = c(fill = FALSE)) +
  geom_point(aes(fill = field_type), size = sm_size, shape = 21, stroke = lw, show.legend = c(fill = TRUE)) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "serif", fontface = 2, color = "black") +
  scale_fill_manual(name = "Field type", values = ft_pal) +
  xlab(paste0("PCA 1 (", eig_prop[1], "%)")) +
  ylab(paste0("PCA 2 (", eig_prop[2], "%)")) +
  theme_ord +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme(legend.title = element_text(size = 9, face = 1),
        legend.text = element_text(size = 8, face = 1))
#+ figS4,warning=FALSE,fig.height=3.5,fig.width=6.5
ggsave(root_path("figs", "figS4.svg"), plot = soil_ord_ftypes, device = svglite::svglite,
       width = 5.25, height = 4.25, units = "in")

