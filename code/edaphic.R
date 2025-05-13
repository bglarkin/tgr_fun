#' ---
#' title: "Soil properties"
#' author: "Beau Larkin\n"
#' date: "Last updated: `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   github_document:
#'     toc: true
#'     toc_depth: 3
#'     df_print: paged
#'     fig_width: 8
#'     fig_height: 7
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
packages_needed = c("tidyverse", "vegan", "knitr")
packages_installed = packages_needed %in% rownames(installed.packages())
#+ packages,message=FALSE
if (any(!packages_installed)) {
    install.packages(packages_needed[!packages_installed])
}
#+ libraries,message=FALSE
for (i in 1:length(packages_needed)) {
    library(packages_needed[i], character.only = T)
}

#' ## Root path function
root_path <- function(...) rprojroot::find_rstudio_root_file(...)

#+ conflicts,message=FALSE
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("diversity", "vegan")
#' 
#+ graphics_styles
source(root_path("resources", "styles.R"))

#' 
#' # Functions
#' ## Confidence intervals
#' Calculate upper and lower confidence intervals with alpha=0.05
#+ ci_function
ci_u <- function(x) {(sd(x) / sqrt(length(x))) * qnorm(0.975)}
ci_l <- function(x) {(sd(x) / sqrt(length(x))) * qnorm(0.025)}
#' 

#' 
#' #' # Data
#' ## Site metadata and design
sites <- read_csv(root_path("clean_data/sites.csv"), show_col_types = FALSE) %>% 
    mutate(field_type = factor(field_type, ordered = TRUE, levels = c("corn", "restored", "remnant")))
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
    select(soil_property, units, everything())
    



#' ## PCA ordination, variable correlations, and PERMANOVA

soil_z <- decostand(data.frame(soil[, -1], row.names = 1), "standardize")
soil_pca <- rda(soil_z)
summary(soil_pca)
plot(soil_pca)
#' Axes 1 and 2 explain 52% of the variation in sites. Axes 1 through 6 account for 91%. 
#' 
#' ## Soil variable loadings and correlations
#' Which soil properties explain the most variation among sites?


site_sco <- scores(soil_pca, display = "sites", choices = c(1,2))

soil_cor <- 
    data.frame(cor(soil_z, site_sco)) %>% 
    mutate(PCA_correlation = sqrt(PC1^2 + PC2^2)) %>% 
    arrange(-PCA_correlation) %>% 
    rownames_to_column(var = "soil_property")

#' Use the variable correlations to sort the soil property averages in a table 
#' highlighting field types:
soil_ft_avg %>% 
    left_join(soil_cor %>% select(soil_property, PCA_correlation), by = join_by(soil_property)) %>% 
    arrange(-PCA_correlation) %>% 
    mutate(across(where(is.numeric), ~ round(.x, 2))) %>% 
    select(PCA_correlation, everything()) %>% 
    kable(format = "pandoc")


#' Axis 1 & 2 eigenvalue proportions
eig_prop <- round(summary(soil_pca)$cont$importance[2, 1:2] * 100, 1)


soil_ord_scores <-
    site_sco %>%
    data.frame() %>% 
    rownames_to_column(var = "field_name") %>%
    left_join(sites, by = join_by(field_name))


soilperm <- function(clust_vec, clust_var) {
    #' ### Permanova
    soil_d <- dist(site_sco)
    
    # Test multivariate dispersions
    soil_disper <- betadisper(soil_d, clust_vec, bias.adjust = TRUE)
    soil_mvdisper <- permutest(soil_disper, pairwise = TRUE, permutations = 1999)
    # Global PERMANOVA
    soil_gl_permtest <- adonis2(
        as.formula(paste("soil_d ~", clust_var)),
        data = soil_ord_scores,
        permutations = 1999,
        by = "terms")
    # Pairwise PERMANOVA
    group_var <- as.character(clust_vec)
    groups <- as.data.frame(t(combn(unique(group_var), m = 2)))
    soil_contrasts <- data.frame(
        group1 = groups$V1,
        group2 = groups$V2,
        R2 = NA,
        F_value = NA,
        df1 = NA,
        df2 = NA,
        p_value = NA
    )
    for (i in seq(nrow(soil_contrasts))) {
        group_subset <-
            group_var == soil_contrasts$group1[i] |
            group_var == soil_contrasts$group2[i]
        contrast_d <- as.matrix(soil_d)[group_subset, group_subset]
        fit <- adonis2(
            contrast_d ~ group_var[group_subset],
            permutations = 1999,
            by = "terms")
        # Prepare contrasts table
        soil_contrasts$R2[i] <- round(fit[grep("group_var", rownames(fit)), "R2"], digits = 3)
        soil_contrasts$F_value[i] <- round(fit[grep("group_var", rownames(fit)), "F"], digits = 3)
        soil_contrasts$df1[i] <- fit[grep("group_var", rownames(fit)), "Df"]
        soil_contrasts$df2[i] <- fit[grep("Residual", rownames(fit)), "Df"]
        soil_contrasts$p_value[i] <- fit[grep("group_var", rownames(fit)), 5]
    }
    soil_contrasts$p_value_adj <- p.adjust(soil_contrasts$p_value, method = "fdr") %>% round(., 4)
    
    out = list(
        mvdisper = soil_mvdisper,
        gl_permtest = soil_gl_permtest,
        contrasts = soil_contrasts
    )
    
    return(out)
}


soilperm_field_type <- soilperm(soil_ord_scores$field_type, "field_type")
soilperm_field_type$gl_permtest
ggplot(soil_ord_scores, aes(x = PC1, y = PC2)) +
    geom_polygon(data = soil_ord_scores %>% group_by(field_type) %>% slice(chull(PC1, PC2)),
                 aes(group = field_type), fill = "transparent", color = "black") +
    geom_point(aes(fill = field_type, shape = region), size = sm_size, stroke = lw) +
    scale_fill_manual(values = c("gray", "black", "white")) +
    scale_shape_manual(values = c(22:25)) +
    labs(
        x = paste0("Axis 1 (", eig_prop[1], "%)"),
        y = paste0("Axis 2 (", eig_prop[2], "%)")) +
    theme_ord +
    theme(legend.position = "none")


soilperm_region <- soilperm(soil_ord_scores$region, "region")
soilperm_region$gl_permtest
ggplot(soil_ord_scores, aes(x = PC1, y = PC2)) +
    geom_polygon(data = soil_ord_scores %>% group_by(region) %>% slice(chull(PC1, PC2)),
                 aes(group = region), fill = "transparent", color = "black") +
    geom_point(aes(fill = field_type, shape = region), size = sm_size, stroke = lw) +
    scale_fill_manual(values = c("gray", "black", "white")) +
    scale_shape_manual(values = c(22:25)) +
    labs(
        x = paste0("Axis 1 (", eig_prop[1], "%)"),
        y = paste0("Axis 2 (", eig_prop[2], "%)")) +
    theme_ord +
    theme(legend.position = "none")






#' Plotting results: 
soil_ord_centers <- soil_ord_scores %>% 
    group_by(field_type) %>% 
    summarize(across(starts_with("PC"), list(mean = mean, ci_l = ci_l, ci_u = ci_u), .names = "{.fn}_{.col}"), .groups = "drop") %>% 
    mutate(across(c(ci_l_PC1, ci_u_PC1), ~ mean_PC1 + .x),
           across(c(ci_l_PC2, ci_u_PC2), ~ mean_PC2 + .x))
soil_ord <- 
    ggplot(soil_ord_scores, aes(x = PC1, y = PC2)) +
    geom_point(aes(fill = field_type, shape = region), size = sm_size, stroke = lw) +
    geom_text(aes(label = yr_since), size = yrtx_size, family = "serif", fontface = 2, color = "white") +
    geom_linerange(data = soil_ord_centers, aes(x = mean_PC1, y = mean_PC2, xmin = ci_l_PC1, xmax = ci_u_PC1), linewidth = lw) +
    geom_linerange(data = soil_ord_centers, aes(x = mean_PC1, y = mean_PC2, ymin = ci_l_PC2, ymax = ci_u_PC2), linewidth = lw) +
    geom_point(data = soil_ord_centers, aes(x = mean_PC1, y = mean_PC2, fill = field_type), size = lg_size, stroke = lw, shape = 21) +
    scale_fill_manual(values = c("gray", "black", "white")) +
    scale_shape_manual(values = c(22:25)) +
    labs(
        x = paste0("Axis 1 (", eig_prop[1], "%)"),
        y = paste0("Axis 2 (", eig_prop[2], "%)")) +
    theme_ord +
    guides(fill = guide_legend(position = "inside"),
           shape = guide_legend(position = "inside")) +
    theme(legend.justification = c(0.03, 0.98),
          legend.box = "horizontal")
#+ soil_ord_plot,warning=FALSE,fig.height=5,fig.width=5,fig.align='center',echo=FALSE
soil_ord
