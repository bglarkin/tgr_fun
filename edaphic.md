Soil properties
================
Beau Larkin

Last updated: 13 February, 2026

- [Description](#description)
- [Packages and libraries](#packages-and-libraries)
  - [Root path function](#root-path-function)
  - [Soil properties](#soil-properties)
  - [Distance-based MEM](#distance-based-mem)
- [Results](#results)
  - [Averages in field types](#averages-in-field-types)
  - [Boxplot displays](#boxplot-displays)
  - [PCA ordination, variable correlations, and
    PERMANOVA](#pca-ordination-variable-correlations-and-permanova)
  - [Test spatial structure on soil
    data](#test-spatial-structure-on-soil-data)
  - [Soil variable loadings and
    correlations](#soil-variable-loadings-and-correlations)

# Description

Soil nutrients were analyzed by [Ward Laboratories,
Inc.](https://www.wardlab.com/services/soil-health-analysis/), analysis
methods available in local files or at the link included here. Soil
organic matter is in percent determined by the loss-on-ignition method.
Soil pH is in a log scale as is typical, and all the other minerals are
in parts per million. This may need to be converted to $mg*kg^{-1}$ or
other unit.

This script provides a quick overview of the soil abiotic property data
and tests differences among field types based on soil properties.

# Packages and libraries

``` r
packages_needed <- c("tidyverse", "knitr", "vegan", "patchwork", "conflicted", 
                     "permute", "geosphere", "ape", "adespatial")

to_install <- setdiff(packages_needed, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)
invisible(lapply(packages_needed, library, character.only = TRUE))
```

## Root path function

``` r
root_path <- function(...) rprojroot::find_rstudio_root_file(...)
```

``` r
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("diversity", "vegan")
```

``` r
source(root_path("resources", "styles.R"))
```

``` r
source(root_path("code", "functions.R"))
```

\#’ \# Data \## Site metadata and design

``` r
sites <- read_csv(root_path("clean_data/sites.csv"), show_col_types = FALSE) %>% 
  mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant")))
```

## Soil properties

``` r
soil <- read_csv(root_path("clean_data/soil.csv"), show_col_types = FALSE)[-c(26:27), ]
soil_units <- read_csv(root_path("clean_data/soil_units.csv"), show_col_types = FALSE)
```

## Distance-based MEM

``` r
coord_tbl <- sites %>% select(long, lat) %>% as.matrix()
rownames(coord_tbl) <- sites$field_name
mem <- dbmem(coord_tbl) %>% as.data.frame()
setequal(sites$field_name, rownames(mem))
```

    ## [1] TRUE

# Results

## Averages in field types

``` r
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
  rowwise() %>% 
  mutate(
    cv = sd(c_across(corn:remnant)) / mean(c_across(corn:remnant)), across(where(is.numeric), ~ round(.x, 2))
    ) %>% 
  arrange(-cv)
```

## Boxplot displays

``` r
soil_p_main <- 
  soil %>% 
  pivot_longer(pH:Na, names_to = "soil_property", values_to = "value") %>% 
  left_join(sites %>% select(field_name, field_type), by = join_by(field_name)) %>% 
  left_join(soil_ft_avg %>% select(soil_property, cv, units), by = join_by(soil_property)) %>% 
  mutate(facet_labs = paste0(soil_property, " (", units, ")"),
         facet_labs = fct_reorder(as.factor(facet_labs), -cv)) %>% 
  ggplot(aes(x = field_type, y = value)) +
  facet_wrap(vars(facet_labs), ncol = 4, scales = "free_y") +
  labs(x = NULL, y = NULL) +
  geom_boxplot(aes(fill = field_type)) +
  scale_fill_manual(name = "Field type", values = ft_pal) +
  theme_corf +
  theme(legend.position = "none")
soil_p_legend <- soil_p_main + theme(legend.position = "right")
```

``` r
ggsave(root_path("figs", "figS4.svg"), plot = soil_p_main, device = svglite::svglite,
       width = 19, height = 15, units = "cm")
ggsave(root_path("figs", "figS4_legend.svg"), plot = soil_p_legend, device = svglite::svglite,
       width = 19, height = 15, units = "cm")
```

## PCA ordination, variable correlations, and PERMANOVA

``` r
soil_z <- decostand(data.frame(soil[, -1], row.names = 1), "standardize")
soil_pca <- rda(soil_z)
summary(soil_pca)
```

    ## 
    ## Call:
    ## rda(X = soil_z) 
    ## 
    ## Partitioning of variance:
    ##               Inertia Proportion
    ## Total              13          1
    ## Unconstrained      13          1
    ## 
    ## Eigenvalues, and their contribution to the variance 
    ## 
    ## Importance of components:
    ##                          PC1    PC2    PC3    PC4     PC5     PC6     PC7     PC8     PC9     PC10     PC11    PC12     PC13
    ## Eigenvalue            4.4790 2.3301 1.7896 1.4132 1.16188 0.70120 0.38467 0.27747 0.20308 0.129286 0.069183 0.04186 0.019432
    ## Proportion Explained  0.3445 0.1792 0.1377 0.1087 0.08938 0.05394 0.02959 0.02134 0.01562 0.009945 0.005322 0.00322 0.001495
    ## Cumulative Proportion 0.3445 0.5238 0.6614 0.7701 0.85952 0.91346 0.94305 0.96440 0.98002 0.989963 0.995285 0.99851 1.000000

Axes 1 and 2 explain 52% of the variation in sites. Axes 1 through 6
account for 91%.

## Test spatial structure on soil data

Using db-MEM

``` r
setequal(rownames(soil_z), rownames(mem))
```

    ## [1] TRUE

``` r
forward.sel(soil_z, mem, alpha = 0.05, nperm = 1999)
```

    ## Testing variable 1
    ## Testing variable 2
    ## Testing variable 3
    ## Procedure stopped (alpha criteria): pvalue for variable 3 is 0.263000 (> 0.050000)

    ##   variables order        R2     R2Cum  AdjR2Cum        F pvalue
    ## 1      MEM3     3 0.1621929 0.1621929 0.1257666 4.452622  5e-04
    ## 2      MEM1     1 0.1373308 0.2995238 0.2358441 4.313178  5e-04

``` r
soil_mem_rda <- rda(soil_z, mem[, c(1,3)])
round(RsquareAdj(soil_mem_rda)$adj.r.squared, 3)
```

    ## [1] 0.236

``` r
anova(soil_mem_rda, permutations = 1999) %>% 
  as.data.frame() %>% 
  mutate(p.adj = p.adjust(`Pr(>F)`, "fdr")) %>% 
  kable(, format = "pandoc")
```

|          |  Df | Variance |        F | Pr(\>F) | p.adj |
|----------|----:|---------:|---------:|--------:|------:|
| Model    |   2 | 3.893809 | 4.703603 |   5e-04 | 5e-04 |
| Residual |  22 | 9.106191 |       NA |      NA |    NA |

MEM3 and MEM1 explain 23.6%

## Soil variable loadings and correlations

Which soil properties explain the most variation among sites?

``` r
site_sco <- scores(soil_pca, display = "sites", choices = c(1,2))
soil_cor <- 
    data.frame(cor(soil_z, site_sco)) %>% 
    mutate(PCA_correlation = sqrt(PC1^2 + PC2^2)) %>% 
    arrange(-PCA_correlation) %>% 
    rownames_to_column(var = "soil_property") %>% 
  mutate(across(where(is.numeric), ~ round(.x, 2)))
```

Use the variable correlations to sort the soil property averages in a
table highlighting field types:

``` r
soil_ft_avg %>% 
    left_join(soil_cor %>% select(soil_property, PCA_cor = PCA_correlation), by = join_by(soil_property)) %>% 
  rowwise() %>% 
    mutate(cv = sd(c_across(corn:remnant)) / mean(c_across(corn:remnant)), 
           across(where(is.numeric), ~ round(.x, 2))) %>% 
  arrange(-cv) %>% 
    kable(format = "pandoc")
```

| soil_property | units                |    corn | restored | remnant |   cv | PCA_cor |
|:--------------|:---------------------|--------:|---------:|--------:|-----:|--------:|
| P             | mg/L (Mehlich P-III) |   64.40 |     8.06 |    5.50 | 1.28 |    0.95 |
| NO3           | mg/L                 |   21.54 |     5.07 |    4.38 | 0.94 |    0.85 |
| K             | mg/L                 |  214.40 |   111.62 |   96.00 | 0.46 |    0.57 |
| OM            | % LOI                |    4.68 |     5.26 |    7.28 | 0.24 |    0.95 |
| Ca            | mg/L                 | 2803.20 |  1992.75 | 2856.50 | 0.19 |    0.91 |
| Zn            | mg/L                 |    2.72 |     3.60 |    2.61 | 0.18 |    0.40 |
| SO4           | mg/L                 |   21.20 |    16.69 |   16.00 | 0.16 |    0.79 |
| Cu            | mg/L                 |    2.90 |     2.73 |    2.15 | 0.15 |    0.33 |
| Mn            | mg/L                 |   15.42 |    20.36 |   16.70 | 0.15 |    0.79 |
| Fe            | mg/L                 |   47.34 |    50.16 |   55.92 | 0.09 |    0.27 |
| Na            | mg/L                 |   15.00 |    13.31 |   13.75 | 0.06 |    0.52 |
| Mg            | mg/L                 |  562.40 |   556.56 |  512.75 | 0.05 |    0.87 |
| pH            | NULL                 |    6.88 |     6.42 |    6.68 | 0.03 |    0.72 |

Axis 1 & 2 eigenvalue proportions

``` r
eig_prop <- round(summary(soil_pca)$cont$importance[2, 1:2] * 100, 1)
soil_ord_scores <-
    data.frame(site_sco) %>%
    rownames_to_column(var = "field_name") %>%
    left_join(sites, by = join_by(field_name))
```

### PERMANOVA on field type

``` r
d_soil = dist(soil_z, method = "euclidean")
mva_soil <- soilperm(d = d_soil, env = cbind(sites, mem), covar = c("MEM3", "MEM1"))
```

``` r
mva_soil$dispersion_test
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## Response: Distances
    ##           Df Sum Sq Mean Sq      F N.Perm Pr(>F)
    ## Groups     2  0.066 0.03278 0.0157   1999 0.9805
    ## Residuals 22 45.980 2.08998                     
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##             corn remnant restored
    ## corn             0.99150   0.8810
    ## remnant  0.99277           0.8885
    ## restored 0.87268 0.89386

``` r
mva_soil$permanova
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## adonis2(formula = perm_form, data = env, permutations = nperm, by = "terms")
    ##            Df SumOfSqs      R2      F Pr(>F)    
    ## MEM3        1   50.604 0.16219 6.1548  5e-04 ***
    ## MEM1        1   42.847 0.13733 5.2113  5e-04 ***
    ## field_type  2   54.110 0.17343 3.2906  5e-04 ***
    ## Residual   20  164.439 0.52705                  
    ## Total      24  312.000 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
mva_soil$pairwise_contrasts[c(1,3,2), c(1,2,4,3,7,8)] %>% 
  arrange(group1, desc(group2)) %>% 
  kable(format = "pandoc", caption = "Pairwise permanova contrasts")
```

| group1  | group2   | F_value |    R2 | p_value | p_value_adj |
|:--------|:---------|--------:|------:|--------:|------------:|
| corn    | restored |   5.249 | 0.165 |  0.0005 |      0.0015 |
| corn    | remnant  |   4.627 | 0.277 |  0.0010 |      0.0015 |
| remnant | restored |   0.577 | 0.023 |  0.7470 |      0.7470 |

Pairwise permanova contrasts

### Plotting and Fig S2

``` r
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
```

``` r
ggsave(root_path("figs", "figS5.svg"), plot = soil_ord_ftypes, device = svglite::svglite,
       width = 5.25, height = 4.25, units = "in")
```
