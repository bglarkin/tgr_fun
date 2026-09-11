Soil properties
================
Beau Larkin

Last updated: 11 September, 2026

- [Description](#description)
- [Packages and libraries](#packages-and-libraries)
  - [Root path function](#root-path-function)
  - [Soil properties](#soil-properties)
  - [Distance-based MEM](#distance-based-mem)
- [Results](#results)
  - [Averages in field types](#averages-in-field-types)
  - [Boxplot displays](#boxplot-displays)
  - [Test differences among field
    types](#test-differences-among-field-types)
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
                     "permute", "geosphere", "ape", "adespatial", "broom")

to_install <- setdiff(packages_needed, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)
invisible(lapply(packages_needed, library, character.only = TRUE))
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.2.1     ✔ readr     2.2.0
    ## ✔ forcats   1.0.1     ✔ stringr   1.6.0
    ## ✔ ggplot2   4.0.3     ✔ tibble    3.3.1
    ## ✔ lubridate 1.9.5     ✔ tidyr     1.3.2
    ## ✔ purrr     1.2.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
    ## Loading required package: permute
    ## 
    ## 
    ## Attaching package: 'ape'
    ## 
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     where
    ## 
    ## 
    ## Registered S3 methods overwritten by 'adegraphics':
    ##   method         from
    ##   biplot.dudi    ade4
    ##   kplot.foucart  ade4
    ##   kplot.mcoa     ade4
    ##   kplot.mfa      ade4
    ##   kplot.pta      ade4
    ##   kplot.sepan    ade4
    ##   kplot.statis   ade4
    ##   scatter.coa    ade4
    ##   scatter.dudi   ade4
    ##   scatter.nipals ade4
    ##   scatter.pco    ade4
    ##   score.acm      ade4
    ##   score.mix      ade4
    ##   score.pca      ade4
    ##   screeplot.dudi ade4
    ## 
    ## Registered S3 method overwritten by 'spdep':
    ##   method   from
    ##   plot.mst ape 
    ## 
    ## Registered S3 method overwritten by 'adespatial':
    ##   method          from       
    ##   plot.multispati adegraphics

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

\#’ \# Data \## Site metadata and design Identify plots that need to be
collapsed into single replicate for the biofuel plots. Average location
data from collapsed plots.

``` r
biofuel_plots <- c("FLRSP1", "FLRSP2", "FLRSP3")
sites <- read_csv(root_path("clean_data/sites.csv"), show_col_types = FALSE) %>% 
  mutate(
    biofuel = field_name %in% biofuel_plots,
    field_key = if_else(biofuel, 12, field_key),
    field_name = if_else(biofuel, "FLRSP1", field_name),
    field_code = if_else(biofuel, "FL-6", field_code),
    field_type = factor(field_type, levels = c("corn", "restored", "remnant"))
  ) %>% 
  group_by(field_key, field_name, field_code, field_type, region, yr_restore, yr_since) %>% 
  summarize(across(where(is.numeric), mean), .groups = "drop")
```

## Soil properties

Remove rows 26-27 which were old fields and not applicable here. FLRSP1,
2, and 3 are replicate control plots within a single biofuel experiment,
not independent restored fields. Collapse to FLRSP1.

``` r
soil <- read_csv(root_path("clean_data/soil.csv"), show_col_types = FALSE)[-c(26:27), ] %>% 
  select(-field_key) %>% 
  mutate(field_name = if_else(field_name %in% biofuel_plots, "FLRSP1", field_name)) %>% 
  group_by(field_name) %>% 
  summarize(across(where(is.numeric), mean), .groups = "drop")
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
  geom_boxplot(aes(fill = field_type), shape = 21) +
  scale_fill_manual(name = "Field type", values = ft_pal) +
  theme_corf +
  theme(legend.position = "none")
```

``` r
ggsave(root_path("figs", "figS6.svg"), plot = soil_p_main, 
       device = svglite::svglite, fix_text_size = FALSE, 
       width = 19, height = 20, units = "cm")
```

## Test differences among field types

Use Kruskal-Wallis tests with FDR corrected p values

``` r
soil_kw_data <- 
  soil %>% 
  pivot_longer(pH:Na, names_to = "soil_property", values_to = "value") %>% 
  left_join(sites %>% select(field_name, field_type), by = join_by(field_name))
split(soil_kw_data, soil_kw_data$soil_property) %>% 
  map(\(df) kruskal.test(df$value, df$field_type) %>% 
        tidy()) %>% 
  bind_rows(.id = "property") %>% 
  mutate(p.adj = p.adjust(p.value, "fdr"), 
         across(where(is.numeric), ~ round(.x, 4))) %>% 
  select(property, kw_stat = statistic, p.val = p.value, p.adj) %>% 
  arrange(p.val) %>% 
  kable(format = "pandoc", caption = "Kruskal-Wallis rank sum test results on soil properties across field types.\nDf=2, FDR correction used.")
```

| property | kw_stat |  p.val |  p.adj |
|:---------|--------:|-------:|-------:|
| P        | 12.0440 | 0.0024 | 0.0210 |
| NO3      | 11.4670 | 0.0032 | 0.0210 |
| K        |  7.4989 | 0.0235 | 0.1020 |
| Na       |  2.9518 | 0.2286 | 0.5842 |
| Fe       |  2.6944 | 0.2600 | 0.5842 |
| Zn       |  2.2565 | 0.3236 | 0.5842 |
| SOM      |  2.1520 | 0.3409 | 0.5842 |
| Ca       |  1.8388 | 0.3988 | 0.5842 |
| Cu       |  1.8104 | 0.4045 | 0.5842 |
| pH       |  1.3833 | 0.5007 | 0.6265 |
| SO4      |  1.2694 | 0.5301 | 0.6265 |
| Mn       |  0.3990 | 0.8192 | 0.8874 |
| Mg       |  0.0875 | 0.9572 | 0.9572 |

Kruskal-Wallis rank sum test results on soil properties across field
types. Df=2, FDR correction used.

## PCA ordination, variable correlations, and PERMANOVA

``` r
soil_z <- decostand(data.frame(soil, row.names = 1), "standardize")
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
    ##                          PC1    PC2    PC3    PC4     PC5     PC6     PC7
    ## Eigenvalue            4.5726 2.3551 1.8160 1.4454 1.05915 0.71450 0.37842
    ## Proportion Explained  0.3517 0.1812 0.1397 0.1112 0.08147 0.05496 0.02911
    ## Cumulative Proportion 0.3517 0.5329 0.6726 0.7838 0.86525 0.92021 0.94932
    ##                           PC8     PC9     PC10     PC11     PC12     PC13
    ## Eigenvalue            0.25333 0.16652 0.116394 0.065487 0.038287 0.018816
    ## Proportion Explained  0.01949 0.01281 0.008953 0.005037 0.002945 0.001447
    ## Cumulative Proportion 0.96881 0.98162 0.990570 0.995607 0.998553 1.000000

Axes 1 and 2 explain 53% of the variation in sites. Axes 1 through 6
account for 92%.

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
    ## Procedure stopped (alpha criteria): pvalue for variable 3 is 0.119500 (> 0.050000)

    ##   variables order        R2     R2Cum  AdjR2Cum        F pvalue
    ## 1      MEM3     3 0.1758088 0.1758088 0.1365616 4.479525  0.001
    ## 2      MEM1     1 0.1379889 0.3137977 0.2451775 4.021814  0.001

``` r
soil_mem_rda <- rda(soil_z, mem[, c(1,3)])
round(RsquareAdj(soil_mem_rda)$adj.r.squared, 3)
```

    ## [1] 0.245

``` r
anova(soil_mem_rda, permutations = 1999) %>% 
  as.data.frame() %>% 
  mutate(p.adj = p.adjust(`Pr(>F)`, "fdr")) %>% 
  kable(, format = "pandoc")
```

|          |  Df | Variance |        F | Pr(\>F) | p.adj |
|----------|----:|---------:|---------:|--------:|------:|
| Model    |   2 |  4.07937 | 4.572962 |   5e-04 | 5e-04 |
| Residual |  20 |  8.92063 |       NA |      NA |    NA |

MEM3 and MEM1 explain 24.5%

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
| P             | mg/L (Mehlich P-III) |   64.40 |     8.98 |    5.50 | 1.26 |    0.95 |
| NO3           | mg/L                 |   21.54 |     4.01 |    4.38 | 1.00 |    0.93 |
| K             | mg/L                 |  214.40 |   116.62 |   96.00 | 0.44 |    0.56 |
| SOM           | % LOI                |    4.68 |     5.19 |    7.28 | 0.24 |    0.95 |
| Zn            | mg/L                 |    2.72 |     3.66 |    2.61 | 0.19 |    0.40 |
| Ca            | mg/L                 | 2803.20 |  2029.86 | 2856.50 | 0.18 |    0.91 |
| Cu            | mg/L                 |    2.90 |     2.82 |    2.15 | 0.16 |    0.33 |
| Mn            | mg/L                 |   15.42 |    20.88 |   16.70 | 0.16 |    0.80 |
| SO4           | mg/L                 |   21.20 |    16.88 |   16.00 | 0.15 |    0.78 |
| Fe            | mg/L                 |   47.34 |    43.23 |   55.92 | 0.13 |    0.13 |
| Na            | mg/L                 |   15.00 |    13.12 |   13.75 | 0.07 |    0.56 |
| Mg            | mg/L                 |  562.40 |   567.31 |  512.75 | 0.06 |    0.88 |
| pH            | NULL                 |    6.88 |     6.49 |    6.68 | 0.03 |    0.71 |

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
    ##           Df Sum Sq Mean Sq     F N.Perm Pr(>F)
    ## Groups     2  0.091 0.04553 0.021   1999 0.9755
    ## Residuals 20 43.265 2.16323                    
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##             corn remnant restored
    ## corn             0.91350   0.8245
    ## remnant  0.91894           0.9760
    ## restored 0.82943 0.97691

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
    ## MEM3        1   50.281 0.17581 6.1188 0.0005 ***
    ## MEM1        1   39.465 0.13799 4.8026 0.0015 ** 
    ## field_type  2   48.339 0.16902 2.9413 0.0020 ** 
    ## Residual   18  147.914 0.51718                  
    ## Total      22  286.000 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
mva_soil$pairwise_contrasts[c(1,3,2), c(1,2,4,3,7,8)] %>% 
  arrange(group1, desc(group2)) %>% 
  kable(format = "pandoc", caption = "Pairwise permanova contrasts")
```

| group1  | group2   | F_value |    R2 | p_value | p_value_adj |
|:--------|:---------|--------:|------:|--------:|------------:|
| corn    | restored |   4.350 | 0.153 |  0.0010 |      0.0023 |
| corn    | remnant  |   4.195 | 0.272 |  0.0015 |      0.0023 |
| remnant | restored |   0.569 | 0.025 |  0.7550 |      0.7550 |

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
ggsave(root_path("figs", "figS7.svg"), plot = soil_ord_ftypes, 
       device = svglite::svglite, fix_text_size = FALSE,
       width = 5.25, height = 4.25, units = "in")
```
