Correlations: guilds and environment
================
Beau Larkin

Last updated: 16 May, 2025

- [Description](#description)
- [Packages and libraries](#packages-and-libraries)
  - [Root path function](#root-path-function)
- [Functions](#functions)
  - [Model distribution
    probabilities](#model-distribution-probabilities)
- [Data](#data)
  - [Site metadata and design](#site-metadata-and-design)
  - [Sites-species tables](#sites-species-tables)
  - [Microbial species metadata](#microbial-species-metadata)
  - [Plant data](#plant-data)
- [Data wrangling](#data-wrangling)
  - [Grass-forb index](#grass-forb-index)
  - [Whole soil fungi](#whole-soil-fungi)
  - [AMF](#amf)
- [AMF abundance in families](#amf-abundance-in-families)

# Description

Site‑average guild abundances were correlated with plant functional
groups. Sequence data are compositional, so abundances were centered via
the additive‑log‑ratio (ALR) method [Gloor
2017](https://www.frontiersin.org/articles/10.3389/fmicb.2017.02224) and
[Greenacre 2021](https://doi.org/10.3389/fmicb.2021.727398) before
differential tests. Plant cover (10×1 m quadrats, WI only; M. Healy
2016) was analysed without log‑transformation because values are not
instrument‑bounded and may exceed 100% [Gloor
2017](https://www.frontiersin.org/articles/10.3389/fmicb.2017.02224).
Trait data was obtained from the [TRY database](https://www.try-db.org)
([Kattge 2010](https://doi.org/10.1111/j.2041-210X.2010.00067.x),
[Kattge 2011](https://doi.org/10.1111/j.1365-2486.2011.02451.x)).

# Packages and libraries

``` r
packages_needed <- c(
  "tidyverse", "colorspace", "vegan", "knitr", "conflicted",
  "geosphere", "ape", "performance", "patchwork"
)
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
```

``` r
source(root_path("resources", "styles.R"))
```

# Functions

``` r
# Functions ———————— ####
```

## Model distribution probabilities

Probable distributions of response and residuals. Package performance
prints javascript which doesn’t render on github documents.

``` r
distribution_prob <- function(df) {
  print(
    performance::check_distribution(df) %>% 
      as.data.frame() %>% 
      select(Distribution, p_Residuals) %>% 
      arrange(-p_Residuals) %>% 
      slice_head(n = 3) %>% 
      kable(format = "pandoc"))
  print(
    performance::check_distribution(df) %>% 
      as.data.frame() %>% 
      select(Distribution, p_Response) %>% 
      arrange(-p_Response) %>% 
      slice_head(n = 3) %>% 
      kable(format = "pandoc"))
}
```

# Data

``` r
# Data ———————— ####
```

## Site metadata and design

``` r
sites <- read_csv(paste0(getwd(), "/clean_data/sites.csv"), show_col_types = FALSE) %>% 
  mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant")))
```

## Sites-species tables

List *spe* holds average sequence abundances from repeated samples in
each site. CSV files were produced in `sequence_data.R`.

``` r
spe <- list(
  its_avg     = read_csv(root_path("clean_data/spe_ITS_avg.csv"), show_col_types = FALSE),
  amf_avg     = read_csv(root_path("clean_data/spe_18S_avg.csv"), show_col_types = FALSE),
  amf_avg_uni = read_delim(root_path("otu_tables/18S/18S_avg_4unifrac.tsv"), show_col_types = FALSE)
)
```

## Microbial species metadata

``` r
spe_meta <- list(
  its = read_csv(root_path("clean_data/spe_ITS_metadata.csv"), show_col_types = FALSE) %>% 
    mutate(primary_lifestyle = case_when(str_detect(primary_lifestyle, "_saprotroph$") ~ "saprotroph", 
                                         TRUE ~ primary_lifestyle)),
  amf = read_csv(root_path("clean_data/spe_18S_metadata.csv"), show_col_types = FALSE)
) %>% 
  map(. %>% mutate(across(everything(), ~ replace_na(., "unidentified"))))
```

## Plant data

Abundance in functional groups and by species are only available from
Wisconsin sites. Only C4_grass and forbs are used. Others: C3_grass,
legume, and shrubTree were found previously to have high VIF in models
or were not chosen in forward selection.

``` r
pfg <- read_csv(root_path("clean_data", "plant_traits.csv"), show_col_types = FALSE) %>% 
  select(field_name, C4_grass, forb)
```

# Data wrangling

``` r
# Data wrangling ———————— ####
```

- C4 grass and forb cover are transformed into a single index using PCA
  in restored sites only.
- The OTU abundance tables must be wrangled to perform a log-ratio
  transformation, which reduces data skewness and compositionality bias
  (Aitchison 1986, Gloor et al. 2017). The transformation will be
  applied across guilds for whole soil fungi and families for AMF. Raw
  abundances are kept for plotting. Abundance data are also joined with
  site and env paramaters to facilitate downstream analyses.

## Grass-forb index

C4 grass and forb cover are highly correlated (*r* = -0.91) in restored
prairies. In models or constrained ordinations, they are collinear and
cannot be used simultaneously. An index of grass-forb cover is created
to solve this problem.

``` r
pfg_pca <- 
  pfg %>% 
  left_join(sites %>% select(field_name, field_type), by = join_by(field_name)) %>% 
  filter(field_type == "restored") %>% 
  select(-field_type) %>% 
  column_to_rownames(var = "field_name") %>% 
  decostand(method = "standardize") %>% 
  rda()
pfg_pca %>% summary() # 92% variation on first axis
```

    ## 
    ## Call:
    ## rda(X = .) 
    ## 
    ## Partitioning of variance:
    ##               Inertia Proportion
    ## Total               2          1
    ## Unconstrained       2          1
    ## 
    ## Eigenvalues, and their contribution to the variance 
    ## 
    ## Importance of components:
    ##                          PC1     PC2
    ## Eigenvalue            1.8482 0.15184
    ## Proportion Explained  0.9241 0.07592
    ## Cumulative Proportion 0.9241 1.00000

``` r
gf_index = scores(pfg_pca, choices = 1, display = "sites") %>% 
  data.frame() %>% 
  rename(gf_index = PC1) %>% 
  rownames_to_column(var = "field_name")
```

## Whole soil fungi

### Raw abundances

``` r
its_guab <- 
  spe$its_avg %>% 
  pivot_longer(starts_with("otu"), names_to = "otu_num", values_to = "abund") %>% 
  left_join(spe_meta$its %>% select(otu_num, primary_lifestyle), by = join_by(otu_num)) %>% 
  group_by(field_name, primary_lifestyle) %>% summarize(abund = sum(abund), .groups = "drop") %>% 
  arrange(field_name, -abund) %>% 
  pivot_wider(names_from = "primary_lifestyle", values_from = "abund") %>% 
  select(field_name, unidentified, saprotroph, plant_pathogen, everything())
its_guab_pfg <- 
  its_guab %>% 
  left_join(pfg, by = join_by(field_name)) %>% 
  left_join(gf_index, by = join_by(field_name)) %>% 
  left_join(sites %>% select(field_name, field_type, region, yr_since), by = join_by(field_name)) %>% 
  select(field_name, field_type, yr_since, region, everything())
```

### Log ratio transformed abundances

``` r
its_gulr <- 
  its_guab %>% 
  column_to_rownames(var = "field_name") %>% 
  decostand("rclr", MARGIN = 2) %>%
  # decostand("clr", MARGIN = 2, pseudocount = 0.2) %>%
  # decostand("alr", MARGIN = 2, reference = 1, pseudocount = 0.2) %>%
  rownames_to_column(var = "field_name") %>% 
  as_tibble()
its_gulr_pfg <- 
  its_gulr %>% 
  left_join(pfg, by = join_by(field_name)) %>% 
  left_join(gf_index, by = join_by(field_name)) %>% 
  left_join(sites %>% select(field_name, field_type, region, yr_since), by = join_by(field_name)) %>% 
  select(field_name, field_type, yr_since, region, everything())
```

## AMF

### Raw abundances

``` r
amf_fmab <- 
  spe$amf_avg %>% 
  pivot_longer(starts_with("otu"), names_to = "otu_num", values_to = "abund") %>% 
  left_join(spe_meta$amf %>% select(otu_num, family), by = join_by(otu_num)) %>% 
  group_by(field_name, family) %>% summarize(abund = sum(abund), .groups = "drop") %>% 
  arrange(field_name, -abund) %>% 
  pivot_wider(names_from = "family", values_from = "abund")
amf_fmab_pfg <- 
  amf_fmab %>% 
  left_join(pfg, by = join_by(field_name)) %>% 
  left_join(gf_index, by = join_by(field_name)) %>% 
  left_join(sites %>% select(field_name, field_type, region, yr_since), by = join_by(field_name)) %>% 
  select(field_name, field_type, yr_since, region, everything())
```

### Log ratio transformed abundances

``` r
amf_fmlr <- 
  amf_fmab %>% 
  column_to_rownames(var = "field_name") %>% 
  decostand("rclr", MARGIN = 2) %>% 
  # decostand("alr", MARGIN = 2, reference = 1, pseudocount = 0.2) %>% 
  rownames_to_column(var = "field_name") %>% 
  as_tibble()
amf_fmlr_pfg <- 
  amf_fmlr %>% 
  left_join(pfg, by = join_by(field_name)) %>% 
  left_join(gf_index, by = join_by(field_name)) %>% 
  left_join(sites %>% select(field_name, field_type, region, yr_since), by = join_by(field_name)) %>% 
  select(field_name, field_type, yr_since, region, everything())
```

# AMF abundance in families

``` r
# AMF abundance in families ———————— ####
```

Display raw abundances in a table but separate means with log ratio
transformed data

``` r
amf_fmab_ft <- 
  amf_fmab %>% 
  left_join(sites %>% select(field_name, field_type, yr_since), by = join_by(field_name)) %>% 
  pivot_longer(Glomeraceae:Gigasporaceae, names_to = "family", values_to = "abund") %>% 
  group_by(field_type, family) %>% 
  summarize(abund = mean(abund), .groups = "drop") %>% 
  pivot_wider(names_from = field_type, values_from = abund) %>% 
  mutate(total = rowSums(across(where(is.numeric))), across(where(is.numeric), ~ round(.x, 1))) %>% 
  arrange(-total)
kable(amf_fmab_ft, format = "pandoc", caption = "AMF abundance in families and field types")
```

| family               |   corn | restored | remnant |  total |
|:---------------------|-------:|---------:|--------:|-------:|
| Glomeraceae          | 3001.3 |   2924.2 |  2952.1 | 8877.5 |
| Claroideoglomeraceae |  151.2 |    455.5 |   367.4 |  974.1 |
| Paraglomeraceae      |  330.9 |    228.6 |    89.8 |  649.4 |
| Diversisporaceae     |  110.1 |     81.7 |    61.1 |  252.9 |
| Gigasporaceae        |    3.5 |     17.5 |    10.9 |   31.9 |

AMF abundance in families and field types

Test RCLR transformed abundances across field types for each family

``` r
glom_lm <- lm(Glomeraceae ~ field_type, data = amf_fmlr_pfg)
summary(glom_lm)
```

    ## 
    ## Call:
    ## lm(formula = Glomeraceae ~ field_type, data = amf_fmlr_pfg)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.46762 -0.12054  0.03023  0.18259  0.30843 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)         0.03598    0.10567   0.340    0.737
    ## field_typerestored -0.04618    0.12106  -0.381    0.707
    ## field_typeremnant  -0.04014    0.15850  -0.253    0.802
    ## 
    ## Residual standard error: 0.2363 on 22 degrees of freedom
    ## Multiple R-squared:  0.006636,   Adjusted R-squared:  -0.08367 
    ## F-statistic: 0.07349 on 2 and 22 DF,  p-value: 0.9294

NS

``` r
clar_lm <- lm(Claroideoglomeraceae ~ field_type, data = amf_fmlr_pfg)
summary(clar_lm)
```

    ## 
    ## Call:
    ## lm(formula = Claroideoglomeraceae ~ field_type, data = amf_fmlr_pfg)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.2116 -0.3921  0.0490  0.2913  1.3191 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)         -0.7689     0.2870  -2.679  0.01372 * 
    ## field_typerestored   1.0087     0.3288   3.067  0.00564 **
    ## field_typeremnant    0.7708     0.4305   1.790  0.08717 . 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.6418 on 22 degrees of freedom
    ## Multiple R-squared:  0.2996, Adjusted R-squared:  0.2359 
    ## F-statistic: 4.704 on 2 and 22 DF,  p-value: 0.01991

``` r
distribution_prob(clar_lm)
```

    ## 
    ## 
    ## Distribution    p_Residuals
    ## -------------  ------------
    ## normal              0.71875
    ## cauchy              0.18750
    ## gamma               0.09375
    ## 
    ## 
    ## Distribution    p_Response
    ## -------------  -----------
    ## normal             0.71875
    ## cauchy             0.18750
    ## gamma              0.06250

``` r
leveneTest(Claroideoglomeraceae ~ field_type, data = amf_fmlr_pfg) %>% as.data.frame() %>% kable(format = "pandoc")
```

|       |  Df |   F value |   Pr(\>F) |
|-------|----:|----------:|----------:|
| group |   2 | 0.2616302 | 0.7721571 |
|       |  22 |        NA |        NA |

``` r
TukeyHSD(aov(Claroideoglomeraceae ~ field_type, data = amf_fmlr_pfg))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = Claroideoglomeraceae ~ field_type, data = amf_fmlr_pfg)
    ## 
    ## $field_type
    ##                        diff        lwr      upr     p adj
    ## restored-corn     1.0086676  0.1826065 1.834729 0.0149598
    ## remnant-corn      0.7708261 -0.3107418 1.852394 0.1960301
    ## remnant-restored -0.2378415 -1.1391481 0.663465 0.7870933

Model R2_adj 0.24, p\<0.02

``` r
ggplot(amf_fmlr_pfg, aes(x = field_type, y = Paraglomeraceae)) + geom_boxplot()
```

![](resources/guild_envir_correlations_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
para_lm <- lm(Paraglomeraceae ~ field_type, data = amf_fmlr_pfg)
summary(para_lm)
```

    ## 
    ## Call:
    ## lm(formula = Paraglomeraceae ~ field_type, data = amf_fmlr_pfg)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -2.1910 -1.0760  0.2353  0.9508  2.2856 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)          0.3096     0.5922   0.523    0.606
    ## field_typerestored  -0.2269     0.6785  -0.334    0.741
    ## field_typeremnant   -1.0270     0.8884  -1.156    0.260
    ## 
    ## Residual standard error: 1.324 on 22 degrees of freedom
    ## Multiple R-squared:  0.06421,    Adjusted R-squared:  -0.02087 
    ## F-statistic: 0.7547 on 2 and 22 DF,  p-value: 0.4819

NS

``` r
dive_lm <- lm(Diversisporaceae ~ field_type, data = amf_fmlr_pfg)
summary(dive_lm)
```

    ## 
    ## Call:
    ## lm(formula = Diversisporaceae ~ field_type, data = amf_fmlr_pfg)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -3.0799 -0.6190  0.1989  0.8105  2.0707 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)          0.3472     0.5047   0.688    0.499
    ## field_typerestored  -0.4499     0.5782  -0.778    0.445
    ## field_typeremnant   -0.3706     0.7570  -0.490    0.629
    ## 
    ## Residual standard error: 1.129 on 22 degrees of freedom
    ## Multiple R-squared:  0.02687,    Adjusted R-squared:  -0.0616 
    ## F-statistic: 0.3037 on 2 and 22 DF,  p-value: 0.7411

NS

``` r
giga_lm <- lm(Gigasporaceae ~ field_type, data = amf_fmlr_pfg)
summary(giga_lm)
```

    ## 
    ## Call:
    ## lm(formula = Gigasporaceae ~ field_type, data = amf_fmlr_pfg)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -2.7186 -0.9071  0.3179  0.9320  2.2616 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)         -0.4378     0.5768  -0.759    0.456
    ## field_typerestored   0.5587     0.6608   0.846    0.407
    ## field_typeremnant    0.5012     0.8652   0.579    0.568
    ## 
    ## Residual standard error: 1.29 on 22 degrees of freedom
    ## Multiple R-squared:  0.03196,    Adjusted R-squared:  -0.05604 
    ## F-statistic: 0.3632 on 2 and 22 DF,  p-value: 0.6995

NS

``` r
# Pathogens and plants
# Formally incorporate into this script later...
ggplot(its_guab_pfg %>% filter(field_type == "restored", region != "FL"), aes(x = gf_index, y = plant_pathogen)) +
  geom_smooth(method = "lm") +
  geom_point()
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](resources/guild_envir_correlations_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
gf_patho_lm <- lm(plant_pathogen ~ gf_index, data = its_gulr_pfg %>% filter(field_type == "restored", region != "FL"))
```

Diagnostics

``` r
par(mfrow = c(2,2))
plot(gf_patho_lm) 
```

![](resources/guild_envir_correlations_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

Some residual structure, but a smooth qq fit. Minor leverage with point
4 pulls the slope to more level, risking a type II error rather than
type I. Model is still significant with point 4 removed.

``` r
distribution_prob(gf_patho_lm)
```

    ## 
    ## 
    ## Distribution    p_Residuals
    ## -------------  ------------
    ## normal              0.50000
    ## cauchy              0.18750
    ## beta                0.15625
    ## 
    ## 
    ## Distribution    p_Response
    ## -------------  -----------
    ## normal             0.56250
    ## cauchy             0.15625
    ## beta               0.09375

Response and residuals normal, no transformations warranted and linear
model appropriate.

``` r
summary(gf_patho_lm)
```

    ## 
    ## Call:
    ## lm(formula = plant_pathogen ~ gf_index, data = its_gulr_pfg %>% 
    ##     filter(field_type == "restored", region != "FL"))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.33142 -0.15670  0.03228  0.16276  0.33424 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.25081    0.07504   3.342 0.010196 *  
    ## gf_index     0.62189    0.11521   5.398 0.000647 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2373 on 8 degrees of freedom
    ## Multiple R-squared:  0.7846, Adjusted R-squared:  0.7577 
    ## F-statistic: 29.14 on 1 and 8 DF,  p-value: 0.0006475
