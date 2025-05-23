Correlations: guilds and environment
================
Beau Larkin

Last updated: 23 May, 2025

- [Description](#description)
- [Packages and libraries](#packages-and-libraries)
  - [Root path function](#root-path-function)
- [Functions](#functions)
  - [Model distribution
    probabilities](#model-distribution-probabilities)
  - [Perform Indicator Species
    Analysis](#perform-indicator-species-analysis)
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

## Perform Indicator Species Analysis

Function `inspan()` takes a combined species and sites data frame and
filters OTUs for indicators of field types.

``` r
inspan <- function(spe_sites, spe_meta, nperm=1999) {
  spe <- data.frame(
    spe_sites %>% select(field_key, starts_with("otu")),
    row.names = 1
  )
  grp = spe_sites$field_type
  mp <- multipatt(
    spe, grp, max.order = 1, 
    control = how(nperm = np))
  si <- mp$sign %>% 
    select(index, stat, p.value) %>% 
    mutate(field_type = case_when(index == 1 ~ "corn", 
                                  index == 2 ~ "restored", 
                                  index == 3 ~ "remnant")) %>% 
    filter(p.value < 0.05) %>% 
    rownames_to_column(var = "otu_num") %>%
    select(-index) %>% 
    as_tibble()
  A  <- mp$A %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "otu_num") %>% 
    pivot_longer(cols = corn:remnant, 
                 names_to = "field_type", 
                 values_to = "A")
  B <- mp$B %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "otu_num") %>% 
    pivot_longer(cols = corn:remnant, 
                 names_to = "field_type", 
                 values_to = "B")
  out <- 
    si %>% 
    left_join(A, by = join_by(otu_num, field_type)) %>% 
    left_join(B, by = join_by(otu_num, field_type)) %>% 
    left_join(meta %>% select(-otu_ID), by = join_by(otu_num)) %>% 
    select(otu_num, A, B, stat, p.value, 
           field_type, primary_lifestyle, everything()) %>% 
    arrange(field_type, -stat)
  
  return(out)
  
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
  amf_fmab_pfg %>% 
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
# Pathogen indicators

patho
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":["field_name"],"name":[1],"type":["chr"],"align":["left"]},{"label":["otu_1"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["otu_3"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["otu_7"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["otu_13"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["otu_16"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["otu_21"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["otu_23"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["otu_28"],"name":[9],"type":["dbl"],"align":["right"]},{"label":["otu_33"],"name":[10],"type":["dbl"],"align":["right"]},{"label":["otu_43"],"name":[11],"type":["dbl"],"align":["right"]},{"label":["otu_53"],"name":[12],"type":["dbl"],"align":["right"]},{"label":["otu_58"],"name":[13],"type":["dbl"],"align":["right"]},{"label":["otu_65"],"name":[14],"type":["dbl"],"align":["right"]},{"label":["otu_68"],"name":[15],"type":["dbl"],"align":["right"]},{"label":["otu_87"],"name":[16],"type":["dbl"],"align":["right"]},{"label":["otu_99"],"name":[17],"type":["dbl"],"align":["right"]},{"label":["otu_135"],"name":[18],"type":["dbl"],"align":["right"]},{"label":["otu_137"],"name":[19],"type":["dbl"],"align":["right"]},{"label":["otu_153"],"name":[20],"type":["dbl"],"align":["right"]},{"label":["otu_172"],"name":[21],"type":["dbl"],"align":["right"]},{"label":["otu_179"],"name":[22],"type":["dbl"],"align":["right"]},{"label":["otu_200"],"name":[23],"type":["dbl"],"align":["right"]},{"label":["otu_212"],"name":[24],"type":["dbl"],"align":["right"]},{"label":["otu_279"],"name":[25],"type":["dbl"],"align":["right"]},{"label":["otu_285"],"name":[26],"type":["dbl"],"align":["right"]},{"label":["otu_289"],"name":[27],"type":["dbl"],"align":["right"]},{"label":["otu_294"],"name":[28],"type":["dbl"],"align":["right"]},{"label":["otu_304"],"name":[29],"type":["dbl"],"align":["right"]},{"label":["otu_315"],"name":[30],"type":["dbl"],"align":["right"]},{"label":["otu_317"],"name":[31],"type":["dbl"],"align":["right"]},{"label":["otu_319"],"name":[32],"type":["dbl"],"align":["right"]},{"label":["otu_325"],"name":[33],"type":["dbl"],"align":["right"]},{"label":["otu_332"],"name":[34],"type":["dbl"],"align":["right"]},{"label":["otu_383"],"name":[35],"type":["dbl"],"align":["right"]},{"label":["otu_391"],"name":[36],"type":["dbl"],"align":["right"]},{"label":["otu_408"],"name":[37],"type":["dbl"],"align":["right"]},{"label":["otu_416"],"name":[38],"type":["dbl"],"align":["right"]},{"label":["otu_425"],"name":[39],"type":["dbl"],"align":["right"]},{"label":["otu_432"],"name":[40],"type":["dbl"],"align":["right"]},{"label":["otu_504"],"name":[41],"type":["dbl"],"align":["right"]},{"label":["otu_511"],"name":[42],"type":["dbl"],"align":["right"]},{"label":["otu_521"],"name":[43],"type":["dbl"],"align":["right"]},{"label":["otu_528"],"name":[44],"type":["dbl"],"align":["right"]},{"label":["otu_543"],"name":[45],"type":["dbl"],"align":["right"]},{"label":["otu_552"],"name":[46],"type":["dbl"],"align":["right"]},{"label":["otu_553"],"name":[47],"type":["dbl"],"align":["right"]},{"label":["otu_559"],"name":[48],"type":["dbl"],"align":["right"]},{"label":["otu_564"],"name":[49],"type":["dbl"],"align":["right"]},{"label":["otu_573"],"name":[50],"type":["dbl"],"align":["right"]},{"label":["otu_607"],"name":[51],"type":["dbl"],"align":["right"]},{"label":["otu_649"],"name":[52],"type":["dbl"],"align":["right"]},{"label":["otu_672"],"name":[53],"type":["dbl"],"align":["right"]},{"label":["otu_685"],"name":[54],"type":["dbl"],"align":["right"]},{"label":["otu_692"],"name":[55],"type":["dbl"],"align":["right"]},{"label":["otu_704"],"name":[56],"type":["dbl"],"align":["right"]},{"label":["otu_740"],"name":[57],"type":["dbl"],"align":["right"]},{"label":["otu_749"],"name":[58],"type":["dbl"],"align":["right"]},{"label":["otu_758"],"name":[59],"type":["dbl"],"align":["right"]},{"label":["otu_762"],"name":[60],"type":["dbl"],"align":["right"]},{"label":["otu_777"],"name":[61],"type":["dbl"],"align":["right"]},{"label":["otu_796"],"name":[62],"type":["dbl"],"align":["right"]},{"label":["otu_797"],"name":[63],"type":["dbl"],"align":["right"]},{"label":["otu_802"],"name":[64],"type":["dbl"],"align":["right"]},{"label":["otu_822"],"name":[65],"type":["dbl"],"align":["right"]},{"label":["otu_830"],"name":[66],"type":["dbl"],"align":["right"]},{"label":["otu_833"],"name":[67],"type":["dbl"],"align":["right"]},{"label":["otu_845"],"name":[68],"type":["dbl"],"align":["right"]},{"label":["otu_848"],"name":[69],"type":["dbl"],"align":["right"]},{"label":["otu_861"],"name":[70],"type":["dbl"],"align":["right"]},{"label":["otu_870"],"name":[71],"type":["dbl"],"align":["right"]},{"label":["otu_875"],"name":[72],"type":["dbl"],"align":["right"]},{"label":["otu_921"],"name":[73],"type":["dbl"],"align":["right"]},{"label":["otu_942"],"name":[74],"type":["dbl"],"align":["right"]},{"label":["otu_968"],"name":[75],"type":["dbl"],"align":["right"]},{"label":["otu_969"],"name":[76],"type":["dbl"],"align":["right"]},{"label":["otu_972"],"name":[77],"type":["dbl"],"align":["right"]},{"label":["otu_1000"],"name":[78],"type":["dbl"],"align":["right"]},{"label":["otu_1013"],"name":[79],"type":["dbl"],"align":["right"]},{"label":["otu_1038"],"name":[80],"type":["dbl"],"align":["right"]},{"label":["otu_1079"],"name":[81],"type":["dbl"],"align":["right"]},{"label":["otu_1099"],"name":[82],"type":["dbl"],"align":["right"]},{"label":["otu_1114"],"name":[83],"type":["dbl"],"align":["right"]},{"label":["otu_1151"],"name":[84],"type":["dbl"],"align":["right"]},{"label":["otu_1159"],"name":[85],"type":["dbl"],"align":["right"]},{"label":["otu_1172"],"name":[86],"type":["dbl"],"align":["right"]},{"label":["otu_1189"],"name":[87],"type":["dbl"],"align":["right"]},{"label":["otu_1193"],"name":[88],"type":["dbl"],"align":["right"]},{"label":["otu_1198"],"name":[89],"type":["dbl"],"align":["right"]},{"label":["otu_1211"],"name":[90],"type":["dbl"],"align":["right"]},{"label":["otu_1241"],"name":[91],"type":["dbl"],"align":["right"]},{"label":["otu_1265"],"name":[92],"type":["dbl"],"align":["right"]},{"label":["otu_1292"],"name":[93],"type":["dbl"],"align":["right"]},{"label":["otu_1311"],"name":[94],"type":["dbl"],"align":["right"]},{"label":["otu_1334"],"name":[95],"type":["dbl"],"align":["right"]},{"label":["otu_1351"],"name":[96],"type":["dbl"],"align":["right"]},{"label":["otu_1357"],"name":[97],"type":["dbl"],"align":["right"]},{"label":["otu_1388"],"name":[98],"type":["dbl"],"align":["right"]},{"label":["otu_1412"],"name":[99],"type":["dbl"],"align":["right"]},{"label":["otu_1414"],"name":[100],"type":["dbl"],"align":["right"]},{"label":["otu_1477"],"name":[101],"type":["dbl"],"align":["right"]},{"label":["otu_1496"],"name":[102],"type":["dbl"],"align":["right"]},{"label":["otu_1497"],"name":[103],"type":["dbl"],"align":["right"]},{"label":["otu_1508"],"name":[104],"type":["dbl"],"align":["right"]},{"label":["otu_1535"],"name":[105],"type":["dbl"],"align":["right"]},{"label":["otu_1541"],"name":[106],"type":["dbl"],"align":["right"]},{"label":["otu_1556"],"name":[107],"type":["dbl"],"align":["right"]},{"label":["otu_1582"],"name":[108],"type":["dbl"],"align":["right"]},{"label":["otu_1662"],"name":[109],"type":["dbl"],"align":["right"]},{"label":["otu_1669"],"name":[110],"type":["dbl"],"align":["right"]},{"label":["otu_1674"],"name":[111],"type":["dbl"],"align":["right"]},{"label":["otu_1716"],"name":[112],"type":["dbl"],"align":["right"]},{"label":["otu_1720"],"name":[113],"type":["dbl"],"align":["right"]},{"label":["otu_1746"],"name":[114],"type":["dbl"],"align":["right"]},{"label":["otu_1773"],"name":[115],"type":["dbl"],"align":["right"]},{"label":["otu_1789"],"name":[116],"type":["dbl"],"align":["right"]},{"label":["otu_1794"],"name":[117],"type":["dbl"],"align":["right"]},{"label":["otu_1834"],"name":[118],"type":["dbl"],"align":["right"]},{"label":["otu_1841"],"name":[119],"type":["dbl"],"align":["right"]},{"label":["otu_1910"],"name":[120],"type":["dbl"],"align":["right"]},{"label":["otu_1987"],"name":[121],"type":["dbl"],"align":["right"]},{"label":["otu_2004"],"name":[122],"type":["dbl"],"align":["right"]},{"label":["otu_2011"],"name":[123],"type":["dbl"],"align":["right"]},{"label":["otu_2034"],"name":[124],"type":["dbl"],"align":["right"]},{"label":["otu_2044"],"name":[125],"type":["dbl"],"align":["right"]},{"label":["otu_2048"],"name":[126],"type":["dbl"],"align":["right"]},{"label":["otu_2055"],"name":[127],"type":["dbl"],"align":["right"]},{"label":["otu_2060"],"name":[128],"type":["dbl"],"align":["right"]},{"label":["otu_2064"],"name":[129],"type":["dbl"],"align":["right"]},{"label":["otu_2101"],"name":[130],"type":["dbl"],"align":["right"]},{"label":["otu_2118"],"name":[131],"type":["dbl"],"align":["right"]},{"label":["otu_2149"],"name":[132],"type":["dbl"],"align":["right"]},{"label":["otu_2157"],"name":[133],"type":["dbl"],"align":["right"]},{"label":["otu_2166"],"name":[134],"type":["dbl"],"align":["right"]},{"label":["otu_2203"],"name":[135],"type":["dbl"],"align":["right"]},{"label":["otu_2238"],"name":[136],"type":["dbl"],"align":["right"]},{"label":["otu_2242"],"name":[137],"type":["dbl"],"align":["right"]},{"label":["otu_2264"],"name":[138],"type":["dbl"],"align":["right"]},{"label":["otu_2281"],"name":[139],"type":["dbl"],"align":["right"]},{"label":["otu_2302"],"name":[140],"type":["dbl"],"align":["right"]},{"label":["otu_2311"],"name":[141],"type":["dbl"],"align":["right"]},{"label":["otu_2354"],"name":[142],"type":["dbl"],"align":["right"]},{"label":["otu_2365"],"name":[143],"type":["dbl"],"align":["right"]},{"label":["otu_2385"],"name":[144],"type":["dbl"],"align":["right"]},{"label":["otu_2401"],"name":[145],"type":["dbl"],"align":["right"]},{"label":["otu_2403"],"name":[146],"type":["dbl"],"align":["right"]},{"label":["otu_2410"],"name":[147],"type":["dbl"],"align":["right"]},{"label":["otu_2414"],"name":[148],"type":["dbl"],"align":["right"]},{"label":["otu_2447"],"name":[149],"type":["dbl"],"align":["right"]},{"label":["otu_2460"],"name":[150],"type":["dbl"],"align":["right"]},{"label":["otu_2498"],"name":[151],"type":["dbl"],"align":["right"]},{"label":["otu_2547"],"name":[152],"type":["dbl"],"align":["right"]},{"label":["otu_2607"],"name":[153],"type":["dbl"],"align":["right"]},{"label":["otu_2613"],"name":[154],"type":["dbl"],"align":["right"]},{"label":["otu_2632"],"name":[155],"type":["dbl"],"align":["right"]},{"label":["otu_2644"],"name":[156],"type":["dbl"],"align":["right"]},{"label":["otu_2698"],"name":[157],"type":["dbl"],"align":["right"]},{"label":["otu_2718"],"name":[158],"type":["dbl"],"align":["right"]},{"label":["otu_2724"],"name":[159],"type":["dbl"],"align":["right"]},{"label":["otu_2725"],"name":[160],"type":["dbl"],"align":["right"]},{"label":["otu_2729"],"name":[161],"type":["dbl"],"align":["right"]},{"label":["otu_2732"],"name":[162],"type":["dbl"],"align":["right"]},{"label":["otu_2744"],"name":[163],"type":["dbl"],"align":["right"]},{"label":["otu_2757"],"name":[164],"type":["dbl"],"align":["right"]},{"label":["otu_2785"],"name":[165],"type":["dbl"],"align":["right"]},{"label":["otu_2804"],"name":[166],"type":["dbl"],"align":["right"]},{"label":["otu_2806"],"name":[167],"type":["dbl"],"align":["right"]},{"label":["otu_2815"],"name":[168],"type":["dbl"],"align":["right"]},{"label":["otu_2885"],"name":[169],"type":["dbl"],"align":["right"]},{"label":["otu_2914"],"name":[170],"type":["dbl"],"align":["right"]},{"label":["otu_2927"],"name":[171],"type":["dbl"],"align":["right"]},{"label":["otu_2948"],"name":[172],"type":["dbl"],"align":["right"]},{"label":["otu_2977"],"name":[173],"type":["dbl"],"align":["right"]},{"label":["otu_2981"],"name":[174],"type":["dbl"],"align":["right"]},{"label":["otu_3001"],"name":[175],"type":["dbl"],"align":["right"]},{"label":["otu_3032"],"name":[176],"type":["dbl"],"align":["right"]},{"label":["otu_3035"],"name":[177],"type":["dbl"],"align":["right"]},{"label":["otu_3065"],"name":[178],"type":["dbl"],"align":["right"]},{"label":["otu_3095"],"name":[179],"type":["dbl"],"align":["right"]},{"label":["otu_3110"],"name":[180],"type":["dbl"],"align":["right"]},{"label":["otu_3137"],"name":[181],"type":["dbl"],"align":["right"]},{"label":["otu_3139"],"name":[182],"type":["dbl"],"align":["right"]},{"label":["otu_3147"],"name":[183],"type":["dbl"],"align":["right"]},{"label":["otu_3153"],"name":[184],"type":["dbl"],"align":["right"]}],"data":[{"1":"BBRP1","2":"146.1000","3":"11.30000","4":"473.80000","5":"2.900000","6":"21.600000","7":"0.0000000","8":"49.10000","9":"0.0","10":"35.60000","11":"4.400000","12":"0.0000000","13":"7.5","14":"3.000000","15":"31.200000","16":"4.200000","17":"57.800000","18":"0.000000","19":"0.000000","20":"6.80000","21":"0.0","22":"10.700000","23":"0.0000000","24":"0.00000","25":"0.0","26":"0.0","27":"0.0","28":"0.0","29":"0.7","30":"0.0","31":"0.0","32":"0.0","33":"0.0","34":"0.800000","35":"0.00000","36":"0.0","37":"0.0000000","38":"0.000000","39":"0.000000","40":"0.0","41":"0.600000","42":"0.0","43":"0.0","44":"0","45":"0.0","46":"0","47":"0.000000","48":"0.0","49":"0.7000000","50":"0.0","51":"1.400000","52":"1.0","53":"0.000000","54":"0.0","55":"0.0","56":"0","57":"0.0","58":"0.000000","59":"0.0","60":"0.0","61":"0.0","62":"0.0000000","63":"0.000000","64":"0.0","65":"0.000000","66":"0.0","67":"0.3000000","68":"0.0","69":"0.000000","70":"0.000000","71":"0.0000000","72":"0.0","73":"0.0","74":"0.0","75":"0.0","76":"0.0","77":"0.0","78":"0.0000000","79":"1.0000000","80":"0.0","81":"0.0","82":"0.0","83":"0.0","84":"0.000000","85":"0.0","86":"0.0","87":"0.0","88":"0.6","89":"0.0","90":"0.0","91":"0.0","92":"0.0","93":"0.0","94":"0.0","95":"0.000000","96":"0.0","97":"1.4000000","98":"0.0","99":"0.0","100":"0.0","101":"0.0","102":"0","103":"0.0","104":"0.400000","105":"0.0","106":"0.0","107":"0.0","108":"0.0","109":"0.0","110":"0.0","111":"0.0","112":"0.2","113":"0.000000","114":"0.000000","115":"0.0","116":"0.000000","117":"0.0","118":"0.0","119":"0.0","120":"0.0","121":"0.0","122":"0.0","123":"0.0","124":"0.0","125":"0.0","126":"0.0","127":"0.0","128":"0.0","129":"0.000000","130":"0.0","131":"0.0","132":"0.0","133":"0.0","134":"0.0","135":"0.0","136":"0.0000000","137":"0.0","138":"0.0","139":"0.0","140":"0.0","141":"0","142":"0.0","143":"0.0","144":"0.0000000","145":"0.0","146":"0.0","147":"0.0","148":"0.0","149":"0.0","150":"0.0","151":"0.0","152":"0.0","153":"0.0","154":"0.0","155":"0.0","156":"0.0","157":"0.5","158":"0.0","159":"0.0","160":"0.0","161":"0.0","162":"0.0","163":"0.0","164":"0.0","165":"0.0","166":"0.0","167":"0.0","168":"0.0","169":"0.0","170":"0.0","171":"0.0","172":"0.3","173":"0.0","174":"0.0","175":"0.0","176":"0.0","177":"0.0","178":"0.0000000","179":"0.0","180":"0.0","181":"0.0","182":"0.0","183":"0.0","184":"0.0"},{"1":"ERRP1","2":"508.0000","3":"128.20000","4":"80.20000","5":"108.400000","6":"33.600000","7":"3.6000000","8":"65.90000","9":"2.4","10":"129.50000","11":"19.800000","12":"0.0000000","13":"0.0","14":"10.300000","15":"41.700000","16":"26.100000","17":"67.500000","18":"3.000000","19":"2.900000","20":"6.50000","21":"0.0","22":"147.900000","23":"0.0000000","24":"1.30000","25":"0.0","26":"0.0","27":"42.1","28":"1.6","29":"12.2","30":"13.0","31":"0.0","32":"0.0","33":"0.0","34":"106.100000","35":"0.00000","36":"0.4","37":"0.9000000","38":"1.900000","39":"2.100000","40":"1.1","41":"0.000000","42":"1.0","43":"0.0","44":"0","45":"34.8","46":"0","47":"0.000000","48":"4.9","49":"0.0000000","50":"3.5","51":"1.100000","52":"0.8","53":"0.000000","54":"0.0","55":"0.0","56":"0","57":"6.7","58":"0.000000","59":"0.0","60":"0.0","61":"0.0","62":"5.4000000","63":"0.400000","64":"0.0","65":"0.000000","66":"0.0","67":"0.0000000","68":"0.0","69":"0.000000","70":"0.000000","71":"0.0000000","72":"4.7","73":"0.0","74":"0.0","75":"0.0","76":"0.0","77":"0.0","78":"1.6000000","79":"0.0000000","80":"0.0","81":"0.0","82":"0.0","83":"0.0","84":"0.000000","85":"0.0","86":"0.0","87":"0.0","88":"0.3","89":"0.0","90":"0.0","91":"0.0","92":"3.4","93":"0.0","94":"0.0","95":"0.000000","96":"0.0","97":"2.8000000","98":"0.0","99":"0.0","100":"0.0","101":"0.0","102":"0","103":"1.3","104":"0.000000","105":"0.0","106":"0.0","107":"0.0","108":"0.0","109":"0.0","110":"0.7","111":"0.0","112":"0.0","113":"0.000000","114":"0.000000","115":"0.0","116":"0.000000","117":"0.0","118":"0.0","119":"0.0","120":"0.0","121":"0.0","122":"0.0","123":"0.0","124":"0.0","125":"0.0","126":"0.0","127":"0.7","128":"0.0","129":"0.000000","130":"0.0","131":"0.0","132":"0.0","133":"0.0","134":"0.0","135":"0.0","136":"0.0000000","137":"0.0","138":"0.0","139":"0.0","140":"0.0","141":"0","142":"0.0","143":"0.0","144":"0.0000000","145":"0.0","146":"0.0","147":"0.0","148":"0.0","149":"0.0","150":"0.0","151":"0.0","152":"0.0","153":"0.0","154":"0.0","155":"0.0","156":"0.6","157":"0.0","158":"0.0","159":"0.0","160":"0.0","161":"0.0","162":"0.0","163":"0.5","164":"0.0","165":"0.0","166":"0.0","167":"0.0","168":"0.0","169":"0.4","170":"0.0","171":"0.0","172":"0.0","173":"0.0","174":"0.0","175":"0.0","176":"0.0","177":"0.0","178":"0.0000000","179":"0.0","180":"0.0","181":"0.0","182":"0.0","183":"0.0","184":"0.0"},{"1":"FGC1","2":"63.3000","3":"61.70000","4":"48.70000","5":"311.100000","6":"115.300000","7":"300.5000000","8":"33.10000","9":"26.5","10":"23.40000","11":"35.200000","12":"0.0000000","13":"0.0","14":"11.000000","15":"13.000000","16":"4.200000","17":"0.300000","18":"0.000000","19":"10.200000","20":"0.00000","21":"60.7","22":"0.000000","23":"15.3000000","24":"12.40000","25":"0.0","26":"0.0","27":"0.0","28":"0.0","29":"0.0","30":"0.0","31":"0.0","32":"0.0","33":"36.0","34":"0.000000","35":"0.00000","36":"2.9","37":"0.3000000","38":"0.000000","39":"0.400000","40":"49.4","41":"0.300000","42":"15.3","43":"0.7","44":"0","45":"0.0","46":"0","47":"11.100000","48":"0.0","49":"0.0000000","50":"4.4","51":"0.000000","52":"0.4","53":"0.000000","54":"0.0","55":"0.0","56":"0","57":"2.4","58":"0.000000","59":"0.0","60":"0.0","61":"0.0","62":"18.6000000","63":"0.700000","64":"0.0","65":"0.000000","66":"0.0","67":"0.8000000","68":"0.0","69":"0.000000","70":"0.000000","71":"1.4000000","72":"0.0","73":"0.0","74":"0.0","75":"0.0","76":"0.0","77":"0.0","78":"1.4000000","79":"0.7000000","80":"0.0","81":"0.0","82":"0.0","83":"0.0","84":"0.000000","85":"4.2","86":"0.0","87":"0.0","88":"0.0","89":"0.0","90":"0.0","91":"0.0","92":"1.0","93":"0.0","94":"0.0","95":"0.000000","96":"0.0","97":"0.4000000","98":"0.0","99":"0.0","100":"0.0","101":"0.0","102":"0","103":"0.0","104":"0.000000","105":"0.0","106":"0.0","107":"0.0","108":"3.0","109":"0.0","110":"0.0","111":"0.0","112":"0.0","113":"0.000000","114":"0.000000","115":"0.0","116":"0.000000","117":"0.0","118":"0.0","119":"1.0","120":"0.0","121":"0.0","122":"0.0","123":"0.0","124":"0.0","125":"0.0","126":"0.0","127":"0.0","128":"0.0","129":"0.000000","130":"0.0","131":"0.0","132":"0.0","133":"0.0","134":"0.0","135":"0.0","136":"0.0000000","137":"0.0","138":"0.0","139":"0.0","140":"0.0","141":"0","142":"0.0","143":"0.0","144":"0.0000000","145":"0.0","146":"0.0","147":"0.2","148":"0.9","149":"0.0","150":"0.0","151":"0.0","152":"0.0","153":"0.0","154":"0.0","155":"0.0","156":"0.0","157":"0.0","158":"0.0","159":"0.0","160":"0.0","161":"0.0","162":"0.0","163":"0.0","164":"0.0","165":"0.0","166":"0.0","167":"0.0","168":"0.0","169":"0.0","170":"0.0","171":"0.0","172":"0.0","173":"0.0","174":"0.0","175":"0.0","176":"0.0","177":"0.0","178":"0.0000000","179":"0.0","180":"0.0","181":"0.0","182":"0.0","183":"0.0","184":"0.0"},{"1":"FGREM1","2":"259.2000","3":"0.50000","4":"228.40000","5":"0.400000","6":"0.000000","7":"0.0000000","8":"20.60000","9":"0.0","10":"11.20000","11":"0.000000","12":"8.9000000","13":"135.6","14":"5.600000","15":"13.000000","16":"2.600000","17":"20.500000","18":"61.100000","19":"43.400000","20":"91.40000","21":"0.0","22":"0.000000","23":"0.0000000","24":"0.70000","25":"0.0","26":"0.0","27":"7.9","28":"0.0","29":"0.0","30":"1.7","31":"0.0","32":"0.0","33":"0.0","34":"0.700000","35":"0.00000","36":"0.8","37":"0.0000000","38":"0.000000","39":"2.000000","40":"0.0","41":"55.300000","42":"3.4","43":"0.0","44":"0","45":"0.0","46":"0","47":"0.000000","48":"0.0","49":"0.0000000","50":"0.0","51":"0.000000","52":"2.1","53":"0.000000","54":"0.0","55":"0.0","56":"0","57":"0.0","58":"0.000000","59":"0.0","60":"0.0","61":"0.0","62":"0.0000000","63":"0.000000","64":"0.0","65":"0.000000","66":"0.0","67":"0.0000000","68":"15.2","69":"0.000000","70":"0.000000","71":"0.0000000","72":"0.0","73":"0.0","74":"11.7","75":"0.0","76":"12.0","77":"0.0","78":"0.0000000","79":"0.0000000","80":"0.0","81":"0.0","82":"0.0","83":"0.0","84":"0.000000","85":"0.0","86":"0.0","87":"0.0","88":"0.0","89":"0.0","90":"0.0","91":"0.0","92":"0.0","93":"0.0","94":"0.4","95":"0.000000","96":"0.0","97":"0.0000000","98":"0.0","99":"0.0","100":"0.0","101":"0.0","102":"0","103":"0.2","104":"0.000000","105":"0.0","106":"0.0","107":"0.0","108":"0.0","109":"0.0","110":"0.0","111":"0.0","112":"2.6","113":"0.000000","114":"0.000000","115":"0.0","116":"0.000000","117":"0.0","118":"0.0","119":"0.0","120":"0.0","121":"0.0","122":"0.0","123":"0.0","124":"0.0","125":"0.0","126":"0.0","127":"0.0","128":"0.0","129":"0.000000","130":"1.4","131":"0.0","132":"0.0","133":"1.3","134":"0.0","135":"0.0","136":"0.0000000","137":"0.0","138":"0.0","139":"0.0","140":"0.0","141":"0","142":"0.4","143":"0.0","144":"0.0000000","145":"0.0","146":"0.0","147":"0.0","148":"0.0","149":"0.0","150":"0.0","151":"0.5","152":"0.0","153":"0.4","154":"0.0","155":"0.0","156":"0.0","157":"0.0","158":"0.0","159":"0.0","160":"0.0","161":"0.0","162":"0.0","163":"0.0","164":"0.0","165":"0.0","166":"0.0","167":"0.0","168":"0.0","169":"0.0","170":"0.0","171":"0.4","172":"0.0","173":"0.0","174":"0.0","175":"0.0","176":"0.0","177":"0.0","178":"0.0000000","179":"0.0","180":"0.0","181":"0.0","182":"0.0","183":"0.0","184":"0.0"},{"1":"FGRP1","2":"78.8000","3":"67.30000","4":"649.70000","5":"28.100000","6":"98.300000","7":"18.2000000","8":"29.30000","9":"13.8","10":"19.60000","11":"43.700000","12":"237.8000000","13":"229.0","14":"9.900000","15":"17.000000","16":"3.100000","17":"0.600000","18":"0.000000","19":"23.000000","20":"1.10000","21":"3.7","22":"1.300000","23":"0.0000000","24":"1.70000","25":"0.0","26":"0.0","27":"0.0","28":"0.0","29":"0.0","30":"2.1","31":"78.9","32":"0.0","33":"0.0","34":"4.900000","35":"0.00000","36":"1.0","37":"0.0000000","38":"0.000000","39":"0.000000","40":"0.0","41":"0.500000","42":"2.3","43":"0.0","44":"0","45":"0.0","46":"0","47":"0.000000","48":"1.7","49":"0.0000000","50":"0.0","51":"0.000000","52":"0.0","53":"0.900000","54":"0.0","55":"0.0","56":"0","57":"0.0","58":"0.000000","59":"0.0","60":"0.0","61":"0.0","62":"0.0000000","63":"0.000000","64":"0.0","65":"0.000000","66":"0.5","67":"0.0000000","68":"0.0","69":"0.000000","70":"0.000000","71":"3.7000000","72":"0.0","73":"0.0","74":"0.0","75":"0.0","76":"0.8","77":"0.0","78":"0.0000000","79":"0.6000000","80":"0.0","81":"0.0","82":"0.0","83":"0.0","84":"0.000000","85":"0.0","86":"0.0","87":"0.0","88":"0.4","89":"0.6","90":"0.0","91":"0.0","92":"0.0","93":"0.0","94":"0.0","95":"0.000000","96":"0.0","97":"0.0000000","98":"0.0","99":"0.0","100":"0.0","101":"0.0","102":"0","103":"0.0","104":"0.000000","105":"0.0","106":"0.0","107":"0.0","108":"0.0","109":"1.2","110":"1.1","111":"0.0","112":"0.0","113":"0.000000","114":"0.000000","115":"0.0","116":"0.000000","117":"2.3","118":"0.0","119":"0.0","120":"0.0","121":"0.0","122":"0.0","123":"1.6","124":"0.0","125":"0.0","126":"0.0","127":"1.5","128":"0.0","129":"0.000000","130":"0.0","131":"0.0","132":"0.0","133":"0.0","134":"0.0","135":"0.0","136":"0.0000000","137":"0.0","138":"0.0","139":"0.0","140":"0.0","141":"0","142":"0.6","143":"0.0","144":"0.0000000","145":"0.0","146":"0.0","147":"0.0","148":"0.0","149":"0.0","150":"0.0","151":"0.0","152":"0.0","153":"0.0","154":"0.0","155":"0.0","156":"0.0","157":"0.0","158":"0.0","159":"0.0","160":"0.0","161":"0.0","162":"0.0","163":"0.0","164":"0.0","165":"0.0","166":"0.5","167":"0.5","168":"0.0","169":"0.0","170":"0.0","171":"0.0","172":"0.0","173":"0.0","174":"0.0","175":"0.0","176":"0.0","177":"0.0","178":"0.0000000","179":"0.0","180":"0.0","181":"0.0","182":"0.0","183":"0.0","184":"0.2"},{"1":"FLC1","2":"180.6000","3":"29.50000","4":"162.00000","5":"2.600000","6":"89.500000","7":"209.2000000","8":"72.50000","9":"23.7","10":"4.00000","11":"49.200000","12":"1.1000000","13":"0.3","14":"91.500000","15":"0.300000","16":"20.500000","17":"0.000000","18":"0.000000","19":"9.800000","20":"4.40000","21":"0.7","22":"0.000000","23":"17.4000000","24":"0.00000","25":"0.0","26":"0.0","27":"0.0","28":"0.0","29":"0.0","30":"0.0","31":"0.0","32":"0.0","33":"0.0","34":"0.000000","35":"0.00000","36":"1.1","37":"0.0000000","38":"0.000000","39":"0.000000","40":"0.0","41":"0.000000","42":"0.0","43":"0.0","44":"0","45":"0.0","46":"0","47":"34.900000","48":"0.0","49":"0.0000000","50":"0.6","51":"1.600000","52":"0.0","53":"0.000000","54":"0.0","55":"0.0","56":"0","57":"0.0","58":"0.000000","59":"0.5","60":"0.0","61":"0.7","62":"0.0000000","63":"2.100000","64":"0.0","65":"3.100000","66":"0.0","67":"0.0000000","68":"0.0","69":"0.000000","70":"0.000000","71":"0.5000000","72":"0.0","73":"0.0","74":"0.0","75":"0.0","76":"0.0","77":"0.0","78":"0.0000000","79":"1.6000000","80":"0.0","81":"0.0","82":"0.0","83":"0.0","84":"0.000000","85":"0.0","86":"0.0","87":"0.0","88":"0.0","89":"0.0","90":"0.0","91":"0.0","92":"0.0","93":"0.0","94":"0.0","95":"0.000000","96":"0.0","97":"0.0000000","98":"0.0","99":"0.0","100":"0.0","101":"0.0","102":"0","103":"0.0","104":"0.000000","105":"0.0","106":"0.0","107":"0.0","108":"0.3","109":"0.0","110":"0.0","111":"0.0","112":"0.0","113":"0.000000","114":"0.000000","115":"0.0","116":"0.000000","117":"0.0","118":"0.0","119":"0.4","120":"0.0","121":"0.0","122":"1.0","123":"0.0","124":"0.0","125":"0.0","126":"0.0","127":"0.0","128":"0.0","129":"0.000000","130":"0.0","131":"0.0","132":"0.0","133":"0.0","134":"0.0","135":"0.0","136":"0.0000000","137":"0.0","138":"0.0","139":"0.0","140":"0.0","141":"0","142":"0.0","143":"0.0","144":"0.0000000","145":"0.0","146":"0.0","147":"0.0","148":"0.0","149":"0.0","150":"0.0","151":"0.0","152":"0.0","153":"0.0","154":"0.0","155":"0.0","156":"0.0","157":"0.0","158":"0.0","159":"0.0","160":"0.0","161":"0.0","162":"0.0","163":"0.0","164":"0.0","165":"0.0","166":"0.0","167":"0.0","168":"0.0","169":"0.0","170":"0.0","171":"0.0","172":"0.0","173":"0.0","174":"0.0","175":"0.0","176":"0.0","177":"0.0","178":"0.0000000","179":"0.0","180":"0.0","181":"0.0","182":"0.0","183":"0.0","184":"0.0"},{"1":"FLC2","2":"206.5000","3":"182.10000","4":"15.90000","5":"111.200000","6":"117.500000","7":"85.4000000","8":"32.70000","9":"74.2","10":"38.00000","11":"30.900000","12":"0.0000000","13":"0.0","14":"20.100000","15":"2.700000","16":"9.000000","17":"0.000000","18":"0.700000","19":"9.000000","20":"0.00000","21":"61.2","22":"0.000000","23":"55.4000000","24":"20.60000","25":"0.0","26":"0.0","27":"2.9","28":"0.0","29":"0.0","30":"0.0","31":"0.0","32":"0.0","33":"9.5","34":"0.000000","35":"0.00000","36":"13.0","37":"0.5000000","38":"0.000000","39":"0.600000","40":"5.7","41":"0.000000","42":"9.5","43":"32.5","44":"0","45":"0.0","46":"0","47":"1.100000","48":"0.0","49":"0.0000000","50":"0.0","51":"0.000000","52":"0.5","53":"0.000000","54":"0.0","55":"0.0","56":"0","57":"0.8","58":"0.000000","59":"29.2","60":"0.0","61":"0.0","62":"3.9000000","63":"0.600000","64":"0.0","65":"0.000000","66":"0.0","67":"0.0000000","68":"0.0","69":"0.000000","70":"0.000000","71":"0.0000000","72":"0.0","73":"0.0","74":"0.0","75":"0.0","76":"0.0","77":"1.3","78":"3.0000000","79":"6.0000000","80":"0.0","81":"0.0","82":"0.0","83":"8.2","84":"0.000000","85":"1.1","86":"0.5","87":"0.0","88":"0.0","89":"7.1","90":"0.0","91":"0.0","92":"0.0","93":"0.0","94":"0.0","95":"0.000000","96":"4.3","97":"0.0000000","98":"0.0","99":"0.0","100":"0.0","101":"0.0","102":"0","103":"0.0","104":"0.000000","105":"0.0","106":"9.6","107":"0.0","108":"0.0","109":"0.0","110":"0.0","111":"0.0","112":"0.0","113":"0.000000","114":"0.000000","115":"0.0","116":"0.000000","117":"0.0","118":"0.0","119":"0.0","120":"1.9","121":"0.0","122":"0.2","123":"0.0","124":"0.0","125":"0.0","126":"0.0","127":"0.0","128":"0.0","129":"0.000000","130":"0.0","131":"0.0","132":"0.0","133":"0.0","134":"0.0","135":"0.0","136":"0.0000000","137":"0.4","138":"0.0","139":"0.0","140":"0.0","141":"0","142":"0.0","143":"0.0","144":"0.0000000","145":"0.0","146":"0.0","147":"0.0","148":"0.0","149":"0.0","150":"0.0","151":"0.0","152":"0.0","153":"0.0","154":"0.0","155":"0.3","156":"0.0","157":"0.0","158":"0.0","159":"0.0","160":"0.0","161":"0.0","162":"0.0","163":"0.0","164":"0.0","165":"0.0","166":"0.0","167":"0.0","168":"0.5","169":"0.0","170":"0.0","171":"0.0","172":"0.0","173":"0.0","174":"0.0","175":"0.0","176":"0.0","177":"0.0","178":"0.0000000","179":"0.0","180":"0.0","181":"0.0","182":"0.0","183":"0.0","184":"0.0"},{"1":"FLREM1","2":"224.3000","3":"106.70000","4":"230.20000","5":"20.500000","6":"19.000000","7":"12.1000000","8":"8.90000","9":"12.8","10":"74.60000","11":"22.400000","12":"0.0000000","13":"25.4","14":"7.500000","15":"44.000000","16":"3.000000","17":"22.300000","18":"13.200000","19":"17.400000","20":"7.10000","21":"0.0","22":"35.600000","23":"0.0000000","24":"4.20000","25":"0.0","26":"0.0","27":"7.4","28":"0.0","29":"17.8","30":"1.2","31":"0.0","32":"7.3","33":"0.0","34":"6.000000","35":"0.00000","36":"1.5","37":"2.7000000","38":"0.000000","39":"0.000000","40":"0.0","41":"0.300000","42":"13.1","43":"0.0","44":"0","45":"0.0","46":"0","47":"0.500000","48":"0.0","49":"0.9000000","50":"12.9","51":"0.000000","52":"0.0","53":"5.400000","54":"0.0","55":"8.0","56":"0","57":"1.3","58":"0.000000","59":"0.0","60":"0.0","61":"0.0","62":"0.0000000","63":"0.000000","64":"0.0","65":"6.200000","66":"0.0","67":"0.0000000","68":"0.0","69":"0.000000","70":"0.000000","71":"0.0000000","72":"0.0","73":"0.0","74":"0.0","75":"0.0","76":"0.0","77":"0.0","78":"0.0000000","79":"0.0000000","80":"0.0","81":"0.0","82":"0.0","83":"0.0","84":"0.000000","85":"0.0","86":"0.0","87":"0.0","88":"0.0","89":"0.0","90":"0.0","91":"0.0","92":"0.0","93":"0.0","94":"0.0","95":"0.000000","96":"0.0","97":"0.0000000","98":"0.0","99":"0.0","100":"0.0","101":"0.0","102":"0","103":"0.0","104":"0.000000","105":"0.0","106":"0.0","107":"0.0","108":"0.0","109":"0.0","110":"0.0","111":"0.0","112":"1.1","113":"0.000000","114":"0.000000","115":"0.0","116":"0.400000","117":"0.0","118":"0.0","119":"0.0","120":"0.0","121":"0.0","122":"0.0","123":"0.0","124":"0.0","125":"0.0","126":"0.0","127":"0.0","128":"0.6","129":"0.000000","130":"0.0","131":"0.0","132":"0.0","133":"0.0","134":"0.0","135":"0.0","136":"0.0000000","137":"0.0","138":"0.0","139":"0.0","140":"0.0","141":"0","142":"0.0","143":"0.0","144":"0.0000000","145":"0.0","146":"0.0","147":"0.5","148":"0.0","149":"0.0","150":"0.0","151":"0.0","152":"0.0","153":"0.0","154":"0.0","155":"0.0","156":"0.0","157":"0.0","158":"0.0","159":"0.0","160":"0.0","161":"0.0","162":"0.0","163":"0.0","164":"0.0","165":"0.0","166":"0.0","167":"0.0","168":"0.0","169":"0.0","170":"0.0","171":"0.0","172":"0.0","173":"0.0","174":"0.0","175":"0.0","176":"0.0","177":"0.0","178":"0.0000000","179":"0.0","180":"0.0","181":"0.0","182":"0.0","183":"0.0","184":"0.0"},{"1":"FLRP1","2":"328.2222","3":"160.66667","4":"155.22222","5":"41.444444","6":"7.666667","7":"0.2222222","8":"97.33333","9":"0.0","10":"41.11111","11":"1.444444","12":"1.4444444","13":"57.0","14":"5.666667","15":"11.000000","16":"6.000000","17":"2.888889","18":"10.888889","19":"5.555556","20":"35.22222","21":"0.0","22":"2.888889","23":"0.3333333","24":"12.55556","25":"0.0","26":"0.0","27":"0.0","28":"0.0","29":"0.0","30":"0.0","31":"0.0","32":"0.0","33":"0.0","34":"3.666667","35":"41.77778","36":"0.0","37":"0.0000000","38":"0.000000","39":"0.000000","40":"0.0","41":"24.000000","42":"0.0","43":"0.0","44":"0","45":"0.0","46":"0","47":"1.333333","48":"0.0","49":"0.0000000","50":"0.0","51":"4.777778","52":"0.0","53":"6.444444","54":"0.0","55":"0.0","56":"0","57":"0.0","58":"0.000000","59":"0.0","60":"0.0","61":"0.0","62":"0.4444444","63":"0.000000","64":"0.0","65":"1.111111","66":"0.0","67":"0.0000000","68":"0.0","69":"5.888889","70":"4.222222","71":"0.6666667","72":"0.0","73":"0.0","74":"0.0","75":"0.0","76":"0.0","77":"0.0","78":"0.0000000","79":"0.0000000","80":"0.0","81":"0.0","82":"0.0","83":"0.0","84":"0.000000","85":"0.0","86":"0.0","87":"0.0","88":"0.0","89":"0.0","90":"0.0","91":"0.0","92":"0.0","93":"0.0","94":"0.0","95":"0.000000","96":"0.0","97":"0.8888889","98":"0.0","99":"0.0","100":"0.0","101":"0.0","102":"0","103":"0.0","104":"0.000000","105":"0.0","106":"0.0","107":"0.0","108":"0.0","109":"0.0","110":"0.0","111":"0.0","112":"0.0","113":"2.888889","114":"2.555556","115":"0.0","116":"1.555556","117":"0.0","118":"0.0","119":"0.0","120":"0.0","121":"0.0","122":"0.0","123":"0.0","124":"0.0","125":"0.0","126":"0.0","127":"0.0","128":"0.0","129":"0.000000","130":"0.0","131":"0.0","132":"0.0","133":"0.0","134":"0.0","135":"0.0","136":"0.0000000","137":"0.0","138":"0.0","139":"0.0","140":"0.0","141":"0","142":"0.0","143":"0.0","144":"0.2222222","145":"0.0","146":"0.0","147":"0.0","148":"0.0","149":"0.0","150":"0.0","151":"0.0","152":"0.0","153":"0.0","154":"0.0","155":"0.0","156":"0.0","157":"0.0","158":"0.0","159":"0.0","160":"0.0","161":"0.0","162":"0.0","163":"0.0","164":"0.0","165":"0.0","166":"0.0","167":"0.0","168":"0.0","169":"0.0","170":"0.0","171":"0.0","172":"0.0","173":"0.0","174":"0.0","175":"0.0","176":"0.0","177":"0.0","178":"0.3333333","179":"0.0","180":"0.0","181":"0.0","182":"0.0","183":"0.0","184":"0.0"},{"1":"FLRP4","2":"134.9000","3":"116.80000","4":"212.70000","5":"13.800000","6":"1.500000","7":"0.2000000","8":"8.50000","9":"0.0","10":"50.20000","11":"0.000000","12":"12.9000000","13":"33.3","14":"9.500000","15":"22.900000","16":"0.700000","17":"10.600000","18":"20.100000","19":"3.000000","20":"15.70000","21":"0.0","22":"15.700000","23":"0.0000000","24":"0.00000","25":"0.0","26":"0.0","27":"1.2","28":"2.3","29":"0.5","30":"0.0","31":"0.0","32":"0.0","33":"0.0","34":"6.300000","35":"8.30000","36":"0.5","37":"0.0000000","38":"0.000000","39":"0.000000","40":"0.0","41":"0.300000","42":"0.0","43":"0.0","44":"0","45":"0.0","46":"0","47":"0.000000","48":"0.0","49":"0.0000000","50":"0.0","51":"1.400000","52":"0.0","53":"0.000000","54":"0.0","55":"0.9","56":"0","57":"0.0","58":"0.000000","59":"0.0","60":"0.0","61":"0.4","62":"0.0000000","63":"0.000000","64":"0.0","65":"0.000000","66":"0.0","67":"0.0000000","68":"0.0","69":"0.000000","70":"10.500000","71":"0.7000000","72":"0.0","73":"0.0","74":"0.0","75":"0.0","76":"0.0","77":"0.0","78":"0.0000000","79":"0.0000000","80":"0.0","81":"0.0","82":"0.0","83":"0.0","84":"0.000000","85":"0.0","86":"0.0","87":"0.0","88":"0.0","89":"0.0","90":"0.0","91":"0.8","92":"0.0","93":"0.0","94":"0.0","95":"0.000000","96":"0.0","97":"0.0000000","98":"0.0","99":"0.0","100":"0.0","101":"0.0","102":"0","103":"1.9","104":"0.400000","105":"0.0","106":"0.0","107":"0.0","108":"0.0","109":"0.0","110":"0.0","111":"0.0","112":"0.0","113":"0.000000","114":"2.400000","115":"0.0","116":"0.000000","117":"0.0","118":"0.0","119":"0.0","120":"0.0","121":"0.0","122":"0.0","123":"0.0","124":"0.0","125":"0.0","126":"0.0","127":"0.0","128":"0.8","129":"0.000000","130":"0.0","131":"0.0","132":"0.0","133":"0.0","134":"0.0","135":"0.0","136":"0.0000000","137":"0.0","138":"0.0","139":"0.0","140":"0.4","141":"0","142":"0.0","143":"0.9","144":"0.0000000","145":"0.0","146":"0.0","147":"0.0","148":"0.0","149":"0.0","150":"0.0","151":"0.0","152":"0.0","153":"0.6","154":"0.0","155":"0.0","156":"0.0","157":"0.0","158":"0.0","159":"0.0","160":"0.0","161":"0.0","162":"0.0","163":"0.0","164":"0.0","165":"0.0","166":"0.0","167":"0.0","168":"0.0","169":"0.0","170":"0.0","171":"0.0","172":"0.0","173":"0.0","174":"0.0","175":"0.0","176":"0.0","177":"0.0","178":"0.0000000","179":"0.0","180":"0.0","181":"0.0","182":"0.0","183":"0.0","184":"0.0"},{"1":"FLRP5","2":"175.1000","3":"205.50000","4":"78.40000","5":"62.000000","6":"11.400000","7":"0.0000000","8":"19.90000","9":"1.4","10":"44.00000","11":"0.000000","12":"9.3000000","13":"68.0","14":"5.100000","15":"11.000000","16":"5.000000","17":"12.000000","18":"35.600000","19":"4.100000","20":"17.10000","21":"0.0","22":"11.100000","23":"0.0000000","24":"14.70000","25":"0.0","26":"0.0","27":"0.0","28":"0.0","29":"0.0","30":"0.4","31":"0.0","32":"0.0","33":"0.0","34":"12.700000","35":"12.20000","36":"0.0","37":"0.5000000","38":"8.000000","39":"0.000000","40":"0.4","41":"6.900000","42":"0.7","43":"0.0","44":"0","45":"0.0","46":"0","47":"0.300000","48":"1.9","49":"0.0000000","50":"0.4","51":"3.800000","52":"0.0","53":"11.300000","54":"0.0","55":"14.0","56":"0","57":"0.0","58":"0.000000","59":"0.0","60":"0.0","61":"0.0","62":"0.0000000","63":"0.000000","64":"0.0","65":"1.400000","66":"0.0","67":"0.3000000","68":"0.0","69":"0.700000","70":"0.000000","71":"4.8000000","72":"0.0","73":"0.0","74":"0.0","75":"0.0","76":"0.0","77":"0.0","78":"2.8000000","79":"0.0000000","80":"0.0","81":"0.0","82":"0.0","83":"0.0","84":"0.000000","85":"0.0","86":"0.0","87":"0.0","88":"0.0","89":"0.0","90":"0.0","91":"0.0","92":"0.0","93":"0.0","94":"0.0","95":"0.000000","96":"0.0","97":"0.0000000","98":"0.0","99":"0.0","100":"0.0","101":"0.0","102":"0","103":"0.9","104":"0.000000","105":"0.0","106":"0.0","107":"0.0","108":"0.0","109":"1.0","110":"0.0","111":"0.0","112":"0.0","113":"0.000000","114":"0.000000","115":"2.3","116":"0.000000","117":"0.0","118":"0.0","119":"0.0","120":"0.4","121":"0.0","122":"0.0","123":"0.0","124":"0.0","125":"0.0","126":"0.0","127":"0.0","128":"0.0","129":"0.000000","130":"0.0","131":"0.0","132":"0.0","133":"0.0","134":"1.2","135":"0.0","136":"0.0000000","137":"0.0","138":"0.0","139":"0.0","140":"0.0","141":"0","142":"0.0","143":"0.0","144":"0.0000000","145":"0.0","146":"0.0","147":"0.0","148":"0.0","149":"0.0","150":"0.0","151":"0.0","152":"0.0","153":"0.0","154":"0.0","155":"0.0","156":"0.0","157":"0.0","158":"0.0","159":"0.0","160":"0.0","161":"0.0","162":"0.0","163":"0.0","164":"0.0","165":"0.0","166":"0.0","167":"0.0","168":"0.0","169":"0.0","170":"0.0","171":"0.0","172":"0.0","173":"0.0","174":"0.0","175":"0.0","176":"0.0","177":"0.0","178":"0.0000000","179":"0.0","180":"0.0","181":"0.0","182":"0.0","183":"0.0","184":"0.0"},{"1":"FLRSP1","2":"181.6000","3":"14.40000","4":"74.80000","5":"4.900000","6":"178.400000","7":"0.0000000","8":"14.00000","9":"0.0","10":"22.40000","11":"144.100000","12":"0.0000000","13":"0.0","14":"1.700000","15":"9.900000","16":"9.200000","17":"2.500000","18":"15.500000","19":"3.000000","20":"0.00000","21":"0.0","22":"7.500000","23":"0.0000000","24":"0.00000","25":"0.0","26":"0.0","27":"6.4","28":"0.0","29":"0.0","30":"0.0","31":"0.0","32":"0.0","33":"0.0","34":"14.700000","35":"0.40000","36":"0.0","37":"0.0000000","38":"1.800000","39":"22.700000","40":"0.0","41":"2.600000","42":"0.0","43":"0.0","44":"0","45":"0.0","46":"0","47":"0.200000","48":"2.0","49":"3.1000000","50":"0.0","51":"0.700000","52":"0.0","53":"0.000000","54":"0.0","55":"0.0","56":"0","57":"0.0","58":"0.000000","59":"0.0","60":"0.0","61":"0.0","62":"0.0000000","63":"0.000000","64":"0.0","65":"1.400000","66":"0.0","67":"0.6000000","68":"0.0","69":"0.000000","70":"0.000000","71":"0.0000000","72":"0.0","73":"0.0","74":"0.0","75":"0.0","76":"0.0","77":"1.2","78":"1.7000000","79":"0.0000000","80":"0.0","81":"0.0","82":"0.0","83":"0.0","84":"1.600000","85":"0.0","86":"0.0","87":"0.0","88":"0.0","89":"0.0","90":"0.0","91":"0.0","92":"0.0","93":"0.0","94":"0.0","95":"0.000000","96":"0.0","97":"0.0000000","98":"0.0","99":"0.0","100":"0.0","101":"0.0","102":"0","103":"0.0","104":"0.700000","105":"0.0","106":"0.0","107":"0.0","108":"0.0","109":"0.0","110":"0.0","111":"0.0","112":"0.0","113":"0.000000","114":"0.000000","115":"0.0","116":"0.000000","117":"0.0","118":"0.0","119":"0.0","120":"0.0","121":"0.0","122":"0.0","123":"0.0","124":"0.7","125":"0.0","126":"0.0","127":"0.0","128":"0.0","129":"0.400000","130":"0.0","131":"1.3","132":"0.0","133":"0.0","134":"0.0","135":"0.0","136":"0.0000000","137":"0.0","138":"0.0","139":"0.0","140":"0.0","141":"0","142":"0.0","143":"0.0","144":"0.0000000","145":"0.0","146":"0.0","147":"0.0","148":"0.0","149":"0.8","150":"0.0","151":"0.0","152":"0.0","153":"0.0","154":"0.0","155":"0.0","156":"0.0","157":"0.0","158":"0.0","159":"0.0","160":"0.0","161":"0.0","162":"0.0","163":"0.0","164":"0.0","165":"0.0","166":"0.0","167":"0.0","168":"0.0","169":"0.0","170":"0.0","171":"0.0","172":"0.0","173":"0.0","174":"0.0","175":"0.0","176":"0.0","177":"0.0","178":"0.0000000","179":"0.0","180":"0.0","181":"0.0","182":"0.0","183":"0.0","184":"0.0"},{"1":"FLRSP2","2":"288.6000","3":"63.00000","4":"174.00000","5":"4.900000","6":"47.900000","7":"0.2000000","8":"88.90000","9":"1.4","10":"145.80000","11":"29.200000","12":"0.0000000","13":"0.0","14":"1.700000","15":"1.300000","16":"16.000000","17":"1.400000","18":"6.200000","19":"5.900000","20":"0.00000","21":"0.3","22":"2.600000","23":"0.0000000","24":"0.00000","25":"0.0","26":"0.0","27":"12.1","28":"0.0","29":"5.2","30":"0.0","31":"0.0","32":"0.0","33":"0.0","34":"0.000000","35":"0.00000","36":"0.3","37":"0.0000000","38":"0.000000","39":"1.500000","40":"0.0","41":"1.300000","42":"0.0","43":"0.0","44":"0","45":"0.0","46":"0","47":"0.600000","48":"0.0","49":"6.0000000","50":"0.0","51":"2.100000","52":"0.0","53":"0.000000","54":"0.0","55":"0.0","56":"0","57":"0.0","58":"17.700000","59":"0.0","60":"0.0","61":"13.9","62":"0.0000000","63":"0.000000","64":"0.0","65":"5.900000","66":"0.0","67":"4.6000000","68":"0.0","69":"0.000000","70":"0.000000","71":"0.0000000","72":"0.0","73":"0.0","74":"0.0","75":"0.0","76":"0.0","77":"1.8","78":"0.0000000","79":"0.0000000","80":"0.0","81":"0.0","82":"0.0","83":"0.0","84":"1.700000","85":"0.0","86":"0.0","87":"0.0","88":"0.0","89":"0.0","90":"0.0","91":"0.0","92":"0.0","93":"0.0","94":"0.0","95":"4.100000","96":"0.0","97":"0.0000000","98":"0.0","99":"0.0","100":"0.0","101":"0.0","102":"0","103":"0.0","104":"0.000000","105":"0.0","106":"0.0","107":"0.0","108":"0.0","109":"0.0","110":"0.0","111":"0.0","112":"0.7","113":"0.000000","114":"0.000000","115":"0.0","116":"0.000000","117":"0.0","118":"0.0","119":"0.0","120":"0.0","121":"0.0","122":"0.0","123":"0.0","124":"0.0","125":"0.0","126":"0.0","127":"0.0","128":"0.0","129":"0.000000","130":"0.0","131":"0.0","132":"0.0","133":"0.0","134":"0.0","135":"0.0","136":"0.0000000","137":"0.0","138":"0.0","139":"0.0","140":"0.6","141":"0","142":"0.0","143":"0.0","144":"0.0000000","145":"0.0","146":"0.0","147":"0.0","148":"0.0","149":"0.0","150":"0.0","151":"0.0","152":"0.0","153":"0.0","154":"0.2","155":"0.0","156":"0.0","157":"0.0","158":"0.0","159":"0.0","160":"0.0","161":"0.0","162":"0.0","163":"0.0","164":"0.0","165":"0.0","166":"0.0","167":"0.0","168":"0.0","169":"0.0","170":"0.0","171":"0.0","172":"0.0","173":"0.0","174":"0.0","175":"0.0","176":"0.0","177":"0.0","178":"0.0000000","179":"0.0","180":"0.0","181":"0.0","182":"0.0","183":"0.0","184":"0.0"},{"1":"FLRSP3","2":"309.5556","3":"42.55556","4":"50.77778","5":"6.333333","6":"45.555556","7":"0.0000000","8":"21.00000","9":"0.0","10":"14.66667","11":"65.555556","12":"0.7777778","13":"0.0","14":"2.444444","15":"4.777778","16":"9.666667","17":"4.666667","18":"4.444444","19":"1.666667","20":"0.00000","21":"0.0","22":"13.444444","23":"0.0000000","24":"0.00000","25":"0.0","26":"0.0","27":"0.0","28":"0.0","29":"0.0","30":"0.0","31":"0.0","32":"0.0","33":"0.0","34":"13.666667","35":"0.00000","36":"0.0","37":"0.4444444","38":"2.666667","39":"8.444444","40":"0.0","41":"2.666667","42":"0.0","43":"0.0","44":"0","45":"0.0","46":"0","47":"2.000000","48":"0.0","49":"0.7777778","50":"0.0","51":"0.000000","52":"0.0","53":"0.000000","54":"0.0","55":"0.0","56":"0","57":"0.0","58":"2.111111","59":"0.0","60":"0.0","61":"0.0","62":"0.0000000","63":"3.777778","64":"0.0","65":"0.000000","66":"0.0","67":"0.7777778","68":"0.0","69":"0.000000","70":"0.000000","71":"0.0000000","72":"0.0","73":"0.0","74":"0.0","75":"0.0","76":"0.0","77":"8.0","78":"0.7777778","79":"0.5555556","80":"0.0","81":"0.0","82":"0.0","83":"0.0","84":"4.222222","85":"0.0","86":"0.0","87":"0.0","88":"0.0","89":"0.0","90":"0.0","91":"0.0","92":"0.0","93":"0.0","94":"0.0","95":"1.555556","96":"0.0","97":"0.0000000","98":"0.0","99":"0.0","100":"0.0","101":"0.0","102":"0","103":"0.0","104":"2.111111","105":"0.0","106":"0.0","107":"0.0","108":"0.0","109":"0.0","110":"0.0","111":"0.0","112":"0.0","113":"0.000000","114":"0.000000","115":"0.0","116":"0.000000","117":"0.0","118":"0.0","119":"0.0","120":"0.0","121":"0.0","122":"0.0","123":"0.0","124":"0.0","125":"0.0","126":"0.0","127":"0.0","128":"0.0","129":"1.111111","130":"0.0","131":"0.0","132":"0.0","133":"0.0","134":"0.0","135":"0.0","136":"0.5555556","137":"0.0","138":"0.0","139":"0.0","140":"0.0","141":"0","142":"0.0","143":"0.0","144":"0.0000000","145":"0.0","146":"0.0","147":"0.0","148":"0.0","149":"0.0","150":"0.0","151":"0.0","152":"0.0","153":"0.0","154":"0.0","155":"0.0","156":"0.0","157":"0.0","158":"0.0","159":"0.0","160":"0.0","161":"0.0","162":"0.0","163":"0.0","164":"0.0","165":"0.0","166":"0.0","167":"0.0","168":"0.0","169":"0.0","170":"0.0","171":"0.0","172":"0.0","173":"0.0","174":"0.0","175":"0.0","176":"0.0","177":"0.0","178":"0.0000000","179":"0.0","180":"0.0","181":"0.0","182":"0.0","183":"0.0","184":"0.0"},{"1":"KORP1","2":"120.0000","3":"6.10000","4":"215.40000","5":"1.800000","6":"18.400000","7":"0.6000000","8":"25.30000","9":"0.0","10":"105.50000","11":"9.500000","12":"0.0000000","13":"0.0","14":"0.900000","15":"1.000000","16":"29.800000","17":"10.500000","18":"41.500000","19":"10.000000","20":"2.90000","21":"0.0","22":"0.000000","23":"6.2000000","24":"0.00000","25":"1.4","26":"0.0","27":"0.0","28":"1.7","29":"0.0","30":"0.0","31":"0.0","32":"0.0","33":"0.0","34":"0.000000","35":"0.00000","36":"0.0","37":"0.0000000","38":"2.000000","39":"1.200000","40":"0.3","41":"12.400000","42":"0.0","43":"0.5","44":"0","45":"0.0","46":"0","47":"0.000000","48":"0.0","49":"0.5000000","50":"0.0","51":"3.300000","52":"0.0","53":"0.000000","54":"1.8","55":"0.0","56":"0","57":"0.0","58":"0.000000","59":"0.0","60":"0.0","61":"0.0","62":"0.0000000","63":"0.000000","64":"0.0","65":"0.000000","66":"0.0","67":"0.5000000","68":"0.0","69":"0.000000","70":"0.000000","71":"0.0000000","72":"0.0","73":"0.0","74":"0.0","75":"0.0","76":"0.0","77":"0.0","78":"0.0000000","79":"0.0000000","80":"15.4","81":"0.0","82":"0.0","83":"0.0","84":"1.500000","85":"0.0","86":"0.4","87":"0.0","88":"0.0","89":"0.0","90":"0.0","91":"3.4","92":"0.2","93":"0.0","94":"0.0","95":"0.000000","96":"0.0","97":"0.0000000","98":"0.0","99":"0.0","100":"0.0","101":"0.0","102":"0","103":"0.0","104":"0.000000","105":"0.0","106":"0.0","107":"0.0","108":"0.0","109":"0.0","110":"0.0","111":"0.0","112":"0.0","113":"0.000000","114":"0.000000","115":"0.0","116":"0.000000","117":"0.0","118":"0.4","119":"0.0","120":"0.0","121":"0.5","122":"0.0","123":"0.0","124":"0.0","125":"0.0","126":"0.0","127":"0.0","128":"0.0","129":"0.000000","130":"0.0","131":"0.0","132":"0.0","133":"0.0","134":"0.0","135":"0.0","136":"0.0000000","137":"0.0","138":"0.0","139":"0.0","140":"0.0","141":"1","142":"0.0","143":"0.0","144":"0.0000000","145":"0.0","146":"0.0","147":"0.0","148":"0.0","149":"0.0","150":"0.0","151":"0.2","152":"0.0","153":"0.0","154":"0.0","155":"0.0","156":"0.0","157":"0.0","158":"0.5","159":"0.0","160":"0.0","161":"0.0","162":"0.0","163":"0.0","164":"0.0","165":"0.0","166":"0.0","167":"0.0","168":"0.0","169":"0.0","170":"0.0","171":"0.0","172":"0.0","173":"0.3","174":"0.0","175":"0.0","176":"0.0","177":"0.0","178":"0.0000000","179":"0.0","180":"0.0","181":"0.0","182":"0.0","183":"0.0","184":"0.0"},{"1":"LPC1","2":"480.2000","3":"280.40000","4":"21.80000","5":"236.700000","6":"222.100000","7":"357.8000000","8":"12.30000","9":"175.8","10":"29.00000","11":"17.500000","12":"0.0000000","13":"0.0","14":"22.300000","15":"1.200000","16":"16.600000","17":"0.000000","18":"0.400000","19":"0.000000","20":"1.80000","21":"16.9","22":"0.000000","23":"12.7000000","24":"62.40000","25":"0.0","26":"0.0","27":"0.0","28":"0.0","29":"0.0","30":"0.0","31":"0.0","32":"0.3","33":"30.1","34":"0.000000","35":"0.00000","36":"11.5","37":"0.0000000","38":"0.000000","39":"0.000000","40":"11.7","41":"0.000000","42":"0.0","43":"1.3","44":"0","45":"0.0","46":"0","47":"2.800000","48":"8.0","49":"0.0000000","50":"1.5","51":"0.000000","52":"7.4","53":"0.000000","54":"0.0","55":"0.0","56":"0","57":"2.3","58":"0.000000","59":"0.0","60":"0.0","61":"0.0","62":"1.6000000","63":"4.400000","64":"0.0","65":"0.000000","66":"0.0","67":"0.0000000","68":"0.0","69":"0.900000","70":"0.000000","71":"0.0000000","72":"0.0","73":"0.0","74":"0.6","75":"0.0","76":"0.0","77":"0.0","78":"0.8000000","79":"0.0000000","80":"0.0","81":"0.0","82":"0.0","83":"0.0","84":"0.000000","85":"2.1","86":"0.0","87":"0.0","88":"0.0","89":"0.0","90":"0.0","91":"0.0","92":"0.0","93":"0.0","94":"0.0","95":"0.000000","96":"0.0","97":"0.0000000","98":"0.0","99":"0.0","100":"0.0","101":"0.0","102":"0","103":"0.9","104":"0.000000","105":"0.0","106":"0.0","107":"0.0","108":"0.0","109":"0.0","110":"0.3","111":"0.0","112":"0.0","113":"0.000000","114":"0.000000","115":"0.0","116":"0.000000","117":"0.0","118":"0.0","119":"0.2","120":"0.0","121":"0.0","122":"0.0","123":"0.0","124":"0.0","125":"0.0","126":"0.0","127":"0.0","128":"0.0","129":"0.000000","130":"0.0","131":"0.0","132":"0.0","133":"0.0","134":"0.0","135":"0.0","136":"0.6000000","137":"0.7","138":"0.0","139":"0.0","140":"0.0","141":"0","142":"0.0","143":"0.0","144":"0.0000000","145":"0.0","146":"0.0","147":"0.0","148":"0.0","149":"0.0","150":"0.0","151":"0.0","152":"0.0","153":"0.0","154":"0.0","155":"0.0","156":"0.0","157":"0.0","158":"0.0","159":"0.5","160":"0.5","161":"0.0","162":"0.0","163":"0.0","164":"0.0","165":"0.0","166":"0.0","167":"0.0","168":"0.0","169":"0.0","170":"0.0","171":"0.0","172":"0.0","173":"0.0","174":"0.3","175":"0.0","176":"0.0","177":"0.0","178":"0.0000000","179":"0.0","180":"0.0","181":"0.0","182":"0.0","183":"0.0","184":"0.0"},{"1":"LPREM1","2":"486.8000","3":"223.70000","4":"106.00000","5":"55.300000","6":"154.400000","7":"2.7000000","8":"16.50000","9":"76.7","10":"19.00000","11":"100.600000","12":"0.0000000","13":"0.0","14":"39.100000","15":"2.300000","16":"27.800000","17":"1.500000","18":"22.000000","19":"0.000000","20":"1.90000","21":"0.0","22":"0.000000","23":"0.0000000","24":"0.30000","25":"0.0","26":"0.0","27":"18.7","28":"0.0","29":"15.5","30":"0.0","31":"0.0","32":"12.4","33":"0.0","34":"8.800000","35":"0.00000","36":"0.0","37":"0.0000000","38":"33.200000","39":"9.700000","40":"0.0","41":"0.700000","42":"0.0","43":"0.0","44":"0","45":"0.0","46":"0","47":"1.900000","48":"0.0","49":"3.3000000","50":"0.0","51":"1.100000","52":"0.0","53":"0.000000","54":"0.0","55":"0.0","56":"0","57":"0.0","58":"0.000000","59":"0.0","60":"19.1","61":"0.0","62":"0.0000000","63":"0.000000","64":"0.0","65":"4.100000","66":"0.0","67":"2.3000000","68":"0.0","69":"8.000000","70":"0.000000","71":"0.0000000","72":"0.0","73":"0.0","74":"0.7","75":"0.0","76":"0.0","77":"0.0","78":"0.0000000","79":"0.0000000","80":"0.0","81":"0.0","82":"0.0","83":"0.0","84":"0.000000","85":"0.0","86":"0.0","87":"0.0","88":"0.0","89":"0.0","90":"0.0","91":"0.0","92":"0.0","93":"0.0","94":"0.0","95":"0.000000","96":"0.0","97":"0.0000000","98":"0.0","99":"0.0","100":"0.0","101":"0.0","102":"0","103":"0.0","104":"0.000000","105":"0.0","106":"0.0","107":"0.0","108":"0.0","109":"0.0","110":"0.0","111":"0.0","112":"0.0","113":"0.000000","114":"0.000000","115":"0.0","116":"0.000000","117":"0.0","118":"0.3","119":"0.0","120":"0.0","121":"0.0","122":"0.0","123":"0.0","124":"0.0","125":"0.0","126":"0.0","127":"0.0","128":"0.0","129":"0.000000","130":"0.0","131":"0.0","132":"0.0","133":"0.0","134":"0.0","135":"0.0","136":"0.0000000","137":"0.0","138":"0.0","139":"0.0","140":"0.0","141":"0","142":"0.9","143":"0.0","144":"0.0000000","145":"0.0","146":"0.0","147":"0.0","148":"0.0","149":"0.0","150":"0.8","151":"0.0","152":"0.0","153":"0.0","154":"0.0","155":"0.3","156":"0.0","157":"0.0","158":"0.0","159":"0.0","160":"0.0","161":"0.5","162":"0.5","163":"0.0","164":"0.0","165":"0.0","166":"0.0","167":"0.0","168":"0.0","169":"0.0","170":"0.0","171":"0.0","172":"0.0","173":"0.0","174":"0.0","175":"0.0","176":"0.0","177":"0.0","178":"0.0000000","179":"0.2","180":"0.0","181":"0.0","182":"0.0","183":"0.0","184":"0.0"},{"1":"LPRP1","2":"725.5000","3":"345.40000","4":"79.80000","5":"99.800000","6":"99.700000","7":"123.6000000","8":"59.40000","9":"108.0","10":"123.10000","11":"44.800000","12":"0.0000000","13":"0.0","14":"40.900000","15":"18.700000","16":"50.400000","17":"1.600000","18":"1.100000","19":"0.000000","20":"0.00000","21":"0.0","22":"12.200000","23":"0.0000000","24":"25.60000","25":"0.0","26":"91.9","27":"28.1","28":"0.0","29":"0.0","30":"2.6","31":"0.0","32":"36.2","33":"0.0","34":"5.400000","35":"0.00000","36":"13.0","37":"10.5000000","38":"0.000000","39":"0.200000","40":"0.0","41":"8.700000","42":"10.7","43":"0.0","44":"0","45":"0.0","46":"0","47":"2.100000","48":"0.6","49":"6.4000000","50":"6.7","51":"0.800000","52":"0.0","53":"0.000000","54":"0.0","55":"0.0","56":"0","57":"0.0","58":"0.000000","59":"0.0","60":"0.0","61":"0.0","62":"1.7000000","63":"0.000000","64":"0.0","65":"0.000000","66":"15.2","67":"0.0000000","68":"0.0","69":"0.000000","70":"0.000000","71":"1.4000000","72":"0.4","73":"0.0","74":"0.4","75":"0.0","76":"0.0","77":"0.0","78":"0.0000000","79":"0.0000000","80":"0.0","81":"4.2","82":"0.0","83":"0.0","84":"0.000000","85":"0.0","86":"0.0","87":"6.2","88":"0.0","89":"0.0","90":"0.0","91":"2.2","92":"0.0","93":"0.0","94":"0.0","95":"0.000000","96":"0.0","97":"0.0000000","98":"0.0","99":"0.0","100":"0.0","101":"0.0","102":"0","103":"0.0","104":"0.000000","105":"0.0","106":"0.0","107":"0.0","108":"0.0","109":"0.0","110":"0.0","111":"0.0","112":"0.0","113":"0.000000","114":"0.000000","115":"0.0","116":"0.000000","117":"0.0","118":"0.0","119":"0.0","120":"0.0","121":"0.0","122":"0.0","123":"0.0","124":"0.0","125":"0.0","126":"0.0","127":"0.7","128":"0.0","129":"0.000000","130":"0.0","131":"0.0","132":"0.0","133":"0.0","134":"0.0","135":"0.0","136":"0.0000000","137":"0.0","138":"0.0","139":"0.0","140":"0.0","141":"0","142":"0.0","143":"0.0","144":"0.0000000","145":"0.0","146":"0.0","147":"0.0","148":"0.0","149":"0.0","150":"0.0","151":"0.0","152":"0.2","153":"0.0","154":"0.0","155":"0.0","156":"0.0","157":"0.0","158":"0.0","159":"0.0","160":"0.0","161":"0.0","162":"0.0","163":"0.0","164":"0.0","165":"0.0","166":"0.0","167":"0.0","168":"0.0","169":"0.0","170":"0.0","171":"0.0","172":"0.0","173":"0.0","174":"0.0","175":"0.0","176":"0.0","177":"0.0","178":"0.0000000","179":"0.0","180":"0.0","181":"0.0","182":"0.0","183":"0.0","184":"0.0"},{"1":"LPRP2","2":"364.4000","3":"138.30000","4":"81.00000","5":"32.800000","6":"60.800000","7":"63.4000000","8":"34.30000","9":"23.0","10":"20.00000","11":"89.000000","12":"0.0000000","13":"0.0","14":"28.900000","15":"9.300000","16":"10.000000","17":"7.100000","18":"0.800000","19":"2.700000","20":"1.00000","21":"1.0","22":"8.800000","23":"1.3000000","24":"0.00000","25":"0.0","26":"1.6","27":"10.7","28":"98.5","29":"0.0","30":"1.0","31":"0.0","32":"22.5","33":"0.0","34":"7.200000","35":"0.00000","36":"2.1","37":"0.0000000","38":"0.000000","39":"0.300000","40":"0.0","41":"6.300000","42":"0.0","43":"0.0","44":"0","45":"0.0","46":"0","47":"1.000000","48":"4.1","49":"0.0000000","50":"3.9","51":"0.000000","52":"0.0","53":"0.000000","54":"0.0","55":"0.0","56":"0","57":"3.1","58":"0.000000","59":"0.0","60":"0.0","61":"0.0","62":"0.0000000","63":"6.200000","64":"0.0","65":"0.000000","66":"0.0","67":"3.3000000","68":"0.0","69":"0.000000","70":"0.000000","71":"0.7000000","72":"0.0","73":"0.0","74":"0.0","75":"0.0","76":"0.0","77":"0.0","78":"0.0000000","79":"0.9000000","80":"0.0","81":"0.0","82":"8.4","83":"0.0","84":"0.300000","85":"0.0","86":"0.0","87":"1.0","88":"0.0","89":"0.0","90":"0.0","91":"0.0","92":"0.0","93":"0.0","94":"0.0","95":"0.000000","96":"0.8","97":"0.0000000","98":"0.0","99":"3.1","100":"0.0","101":"0.0","102":"4","103":"0.0","104":"0.000000","105":"3.5","106":"0.0","107":"0.0","108":"0.0","109":"0.0","110":"0.0","111":"0.0","112":"0.0","113":"0.000000","114":"0.700000","115":"0.0","116":"0.000000","117":"0.0","118":"0.0","119":"0.3","120":"0.0","121":"0.0","122":"0.0","123":"0.0","124":"0.3","125":"0.0","126":"0.0","127":"0.0","128":"0.0","129":"0.000000","130":"0.0","131":"0.0","132":"0.9","133":"0.0","134":"0.0","135":"0.0","136":"0.0000000","137":"0.0","138":"1.1","139":"0.0","140":"0.0","141":"0","142":"0.0","143":"0.0","144":"0.0000000","145":"0.0","146":"0.0","147":"0.0","148":"0.0","149":"0.0","150":"0.0","151":"0.0","152":"0.0","153":"0.0","154":"0.0","155":"0.0","156":"0.0","157":"0.0","158":"0.0","159":"0.0","160":"0.0","161":"0.0","162":"0.0","163":"0.0","164":"0.0","165":"0.0","166":"0.0","167":"0.0","168":"0.0","169":"0.0","170":"0.0","171":"0.0","172":"0.0","173":"0.0","174":"0.0","175":"0.3","176":"0.0","177":"0.0","178":"0.0000000","179":"0.0","180":"0.2","181":"0.0","182":"0.0","183":"0.0","184":"0.0"},{"1":"MBREM1","2":"209.5000","3":"3.70000","4":"69.70000","5":"1.900000","6":"3.000000","7":"0.0000000","8":"8.00000","9":"0.7","10":"14.70000","11":"7.900000","12":"0.0000000","13":"4.5","14":"0.200000","15":"0.400000","16":"9.500000","17":"10.700000","18":"22.400000","19":"0.000000","20":"0.00000","21":"0.0","22":"0.000000","23":"1.5000000","24":"0.00000","25":"0.0","26":"2.9","27":"0.0","28":"0.4","29":"0.0","30":"0.0","31":"0.0","32":"0.2","33":"0.0","34":"0.000000","35":"0.00000","36":"0.0","37":"0.0000000","38":"0.800000","39":"0.400000","40":"0.0","41":"12.200000","42":"0.0","43":"0.0","44":"0","45":"0.0","46":"0","47":"0.000000","48":"0.0","49":"0.8000000","50":"0.0","51":"0.000000","52":"0.0","53":"1.600000","54":"5.9","55":"0.0","56":"0","57":"0.0","58":"0.000000","59":"0.0","60":"0.0","61":"0.0","62":"0.0000000","63":"0.000000","64":"0.0","65":"2.000000","66":"0.0","67":"0.0000000","68":"0.0","69":"0.000000","70":"0.000000","71":"0.0000000","72":"0.0","73":"0.2","74":"0.0","75":"0.0","76":"0.0","77":"0.0","78":"0.0000000","79":"0.0000000","80":"0.0","81":"0.0","82":"0.0","83":"0.0","84":"0.000000","85":"0.0","86":"5.0","87":"0.0","88":"0.0","89":"0.0","90":"0.0","91":"0.0","92":"0.0","93":"0.0","94":"0.0","95":"0.000000","96":"0.0","97":"0.0000000","98":"0.0","99":"0.0","100":"0.0","101":"0.0","102":"0","103":"1.9","104":"0.200000","105":"0.0","106":"0.0","107":"0.0","108":"0.0","109":"0.0","110":"0.0","111":"0.0","112":"0.0","113":"0.000000","114":"0.000000","115":"0.0","116":"0.000000","117":"0.0","118":"0.0","119":"0.0","120":"0.0","121":"0.0","122":"0.0","123":"0.0","124":"0.0","125":"1.5","126":"0.0","127":"0.0","128":"0.0","129":"0.000000","130":"0.0","131":"0.0","132":"0.0","133":"0.0","134":"0.0","135":"0.0","136":"0.0000000","137":"0.0","138":"0.0","139":"0.0","140":"0.0","141":"0","142":"0.0","143":"0.0","144":"0.5000000","145":"0.0","146":"0.0","147":"0.0","148":"0.0","149":"0.0","150":"0.0","151":"0.0","152":"0.0","153":"0.0","154":"0.4","155":"0.0","156":"0.0","157":"0.0","158":"0.0","159":"0.0","160":"0.0","161":"0.0","162":"0.0","163":"0.0","164":"0.8","165":"0.0","166":"0.0","167":"0.0","168":"0.0","169":"0.0","170":"0.0","171":"0.0","172":"0.0","173":"0.0","174":"0.0","175":"0.0","176":"0.0","177":"0.0","178":"0.0000000","179":"0.0","180":"0.0","181":"0.0","182":"0.0","183":"0.0","184":"0.0"},{"1":"MBRP1","2":"440.3000","3":"2.50000","4":"103.90000","5":"1.400000","6":"0.000000","7":"1.3000000","8":"17.50000","9":"0.4","10":"20.10000","11":"25.700000","12":"0.0000000","13":"0.3","14":"1.000000","15":"1.100000","16":"26.000000","17":"42.700000","18":"11.300000","19":"0.000000","20":"7.90000","21":"0.0","22":"0.000000","23":"6.7000000","24":"12.70000","25":"0.0","26":"2.6","27":"0.0","28":"6.4","29":"48.2","30":"0.0","31":"0.0","32":"0.0","33":"0.0","34":"0.000000","35":"0.00000","36":"0.0","37":"0.2000000","38":"0.900000","39":"10.800000","40":"0.0","41":"9.800000","42":"0.0","43":"0.0","44":"0","45":"0.0","46":"0","47":"0.000000","48":"0.0","49":"2.8000000","50":"0.0","51":"4.400000","52":"0.0","53":"11.400000","54":"15.4","55":"0.0","56":"0","57":"0.0","58":"0.000000","59":"0.0","60":"1.1","61":"2.3","62":"0.0000000","63":"0.000000","64":"0.0","65":"2.600000","66":"0.0","67":"2.4000000","68":"0.0","69":"0.000000","70":"0.000000","71":"0.0000000","72":"0.0","73":"12.2","74":"0.0","75":"0.0","76":"0.0","77":"0.0","78":"0.0000000","79":"0.8000000","80":"0.0","81":"1.1","82":"0.0","83":"0.0","84":"0.000000","85":"0.0","86":"0.7","87":"0.0","88":"0.9","89":"0.0","90":"0.0","91":"0.0","92":"0.0","93":"0.0","94":"0.5","95":"0.000000","96":"0.0","97":"0.0000000","98":"0.0","99":"0.0","100":"0.0","101":"3.9","102":"0","103":"2.5","104":"0.000000","105":"0.0","106":"0.0","107":"0.0","108":"0.0","109":"0.0","110":"0.0","111":"2.8","112":"0.7","113":"0.000000","114":"0.000000","115":"0.0","116":"0.500000","117":"0.0","118":"0.0","119":"0.0","120":"0.0","121":"1.1","122":"0.0","123":"0.0","124":"0.0","125":"0.0","126":"0.7","127":"0.0","128":"0.0","129":"0.000000","130":"0.0","131":"0.0","132":"0.0","133":"0.0","134":"0.0","135":"0.0","136":"0.0000000","137":"0.0","138":"0.0","139":"1.1","140":"0.0","141":"0","142":"0.0","143":"0.0","144":"0.0000000","145":"0.0","146":"0.0","147":"0.0","148":"0.0","149":"0.0","150":"0.0","151":"0.0","152":"0.0","153":"0.0","154":"0.0","155":"0.0","156":"0.0","157":"0.0","158":"0.0","159":"0.0","160":"0.0","161":"0.0","162":"0.0","163":"0.0","164":"0.0","165":"0.0","166":"0.0","167":"0.0","168":"0.0","169":"0.0","170":"0.0","171":"0.0","172":"0.0","173":"0.0","174":"0.0","175":"0.0","176":"0.0","177":"0.0","178":"0.0000000","179":"0.0","180":"0.0","181":"0.0","182":"0.0","183":"0.0","184":"0.0"},{"1":"MHRP1","2":"579.9000","3":"53.30000","4":"730.90000","5":"25.800000","6":"100.200000","7":"26.8000000","8":"83.00000","9":"2.0","10":"186.70000","11":"0.000000","12":"252.7000000","13":"0.0","14":"58.400000","15":"109.400000","16":"21.100000","17":"18.500000","18":"3.600000","19":"13.000000","20":"17.20000","21":"0.0","22":"40.100000","23":"0.0000000","24":"0.00000","25":"0.0","26":"0.0","27":"0.0","28":"1.7","29":"0.0","30":"0.0","31":"0.0","32":"0.0","33":"0.0","34":"3.500000","35":"0.00000","36":"0.0","37":"0.2000000","38":"1.700000","39":"2.500000","40":"0.0","41":"1.700000","42":"0.0","43":"0.0","44":"37","45":"0.0","46":"0","47":"0.000000","48":"2.1","49":"3.5000000","50":"0.9","51":"6.600000","52":"7.1","53":"0.600000","54":"0.0","55":"0.0","56":"22","57":"0.0","58":"0.000000","59":"0.0","60":"0.0","61":"0.7","62":"0.8000000","63":"1.300000","64":"0.0","65":"0.000000","66":"0.0","67":"1.2000000","68":"0.0","69":"0.000000","70":"0.000000","71":"6.7000000","72":"0.0","73":"0.0","74":"0.0","75":"8.3","76":"0.0","77":"0.0","78":"9.7000000","79":"0.0000000","80":"0.0","81":"0.7","82":"0.0","83":"0.0","84":"0.000000","85":"0.0","86":"0.7","87":"0.0","88":"3.9","89":"0.0","90":"0.0","91":"0.0","92":"8.6","93":"0.0","94":"5.0","95":"0.000000","96":"0.0","97":"0.0000000","98":"0.0","99":"0.0","100":"2.4","101":"0.0","102":"0","103":"0.0","104":"0.000000","105":"0.0","106":"0.0","107":"0.0","108":"0.0","109":"0.0","110":"0.0","111":"0.0","112":"1.4","113":"0.000000","114":"0.000000","115":"0.0","116":"0.000000","117":"0.0","118":"1.4","119":"0.0","120":"0.0","121":"0.0","122":"0.4","123":"0.0","124":"0.0","125":"0.0","126":"0.8","127":"0.0","128":"0.0","129":"0.000000","130":"0.0","131":"0.0","132":"0.0","133":"0.0","134":"0.8","135":"0.0","136":"0.0000000","137":"0.0","138":"0.0","139":"0.0","140":"0.0","141":"0","142":"0.0","143":"0.0","144":"0.2000000","145":"0.9","146":"0.9","147":"0.0","148":"0.0","149":"0.0","150":"0.0","151":"0.0","152":"0.0","153":"0.0","154":"0.0","155":"0.0","156":"0.0","157":"0.0","158":"0.0","159":"0.0","160":"0.0","161":"0.0","162":"0.0","163":"0.0","164":"0.0","165":"0.0","166":"0.0","167":"0.0","168":"0.0","169":"0.0","170":"0.0","171":"0.0","172":"0.0","173":"0.0","174":"0.0","175":"0.0","176":"0.0","177":"0.0","178":"0.0000000","179":"0.0","180":"0.0","181":"0.0","182":"0.0","183":"0.0","184":"0.0"},{"1":"MHRP2","2":"721.3000","3":"133.60000","4":"124.00000","5":"143.100000","6":"99.900000","7":"49.9000000","8":"90.70000","9":"365.0","10":"160.50000","11":"2.900000","12":"0.0000000","13":"0.0","14":"54.100000","15":"94.700000","16":"14.700000","17":"2.600000","18":"2.300000","19":"37.500000","20":"8.40000","21":"0.0","22":"0.000000","23":"0.0000000","24":"0.00000","25":"95.5","26":"0.0","27":"12.3","28":"0.0","29":"0.0","30":"40.9","31":"0.0","32":"0.0","33":"0.0","34":"23.700000","35":"0.00000","36":"4.6","37":"35.0000000","38":"0.000000","39":"0.800000","40":"0.0","41":"0.000000","42":"5.5","43":"3.4","44":"0","45":"0.0","46":"0","47":"4.200000","48":"1.2","49":"3.3000000","50":"0.0","51":"0.000000","52":"4.3","53":"0.000000","54":"0.0","55":"0.0","56":"0","57":"0.0","58":"0.000000","59":"0.0","60":"0.0","61":"0.6","62":"4.2000000","63":"1.300000","64":"16.8","65":"0.000000","66":"0.0","67":"0.0000000","68":"0.0","69":"0.000000","70":"0.000000","71":"0.3000000","72":"0.0","73":"0.0","74":"0.0","75":"2.6","76":"0.0","77":"0.0","78":"0.0000000","79":"0.0000000","80":"0.0","81":"2.7","82":"0.0","83":"0.0","84":"0.000000","85":"0.0","86":"0.6","87":"0.0","88":"0.6","89":"0.0","90":"6.9","91":"0.0","92":"0.0","93":"0.0","94":"3.9","95":"0.000000","96":"0.0","97":"0.0000000","98":"6.9","99":"0.0","100":"0.0","101":"0.0","102":"0","103":"0.5","104":"0.000000","105":"0.0","106":"0.0","107":"3.4","108":"0.0","109":"0.6","110":"0.0","111":"0.0","112":"0.0","113":"0.000000","114":"0.000000","115":"0.0","116":"0.000000","117":"0.0","118":"0.0","119":"0.0","120":"0.0","121":"0.0","122":"0.0","123":"0.0","124":"0.0","125":"0.0","126":"0.0","127":"0.0","128":"0.0","129":"0.000000","130":"0.0","131":"0.0","132":"1.3","133":"0.0","134":"0.0","135":"1.2","136":"0.0000000","137":"0.0","138":"0.0","139":"0.0","140":"0.0","141":"0","142":"0.0","143":"0.0","144":"0.0000000","145":"0.0","146":"0.0","147":"0.4","148":"0.0","149":"0.0","150":"0.0","151":"0.3","152":"0.5","153":"0.0","154":"0.0","155":"0.0","156":"0.0","157":"0.0","158":"0.0","159":"0.0","160":"0.0","161":"0.0","162":"0.0","163":"0.0","164":"0.0","165":"0.0","166":"0.0","167":"0.0","168":"0.0","169":"0.0","170":"0.4","171":"0.0","172":"0.0","173":"0.0","174":"0.0","175":"0.0","176":"0.3","177":"0.3","178":"0.0000000","179":"0.0","180":"0.0","181":"0.2","182":"0.2","183":"0.0","184":"0.0"},{"1":"PHC1","2":"139.0000","3":"77.00000","4":"29.90000","5":"142.100000","6":"54.600000","7":"121.8000000","8":"102.10000","9":"20.5","10":"51.60000","11":"0.000000","12":"0.0000000","13":"0.0","14":"6.300000","15":"6.100000","16":"7.900000","17":"0.000000","18":"0.000000","19":"2.300000","20":"0.00000","21":"12.5","22":"0.000000","23":"36.4000000","24":"0.00000","25":"0.0","26":"0.0","27":"0.0","28":"0.0","29":"0.0","30":"0.0","31":"0.0","32":"0.0","33":"0.5","34":"0.000000","35":"0.00000","36":"2.3","37":"1.8000000","38":"0.000000","39":"0.400000","40":"2.0","41":"0.000000","42":"0.4","43":"0.0","44":"0","45":"0.0","46":"34","47":"2.000000","48":"0.0","49":"0.0000000","50":"0.0","51":"0.000000","52":"0.3","53":"0.000000","54":"0.0","55":"0.0","56":"0","57":"2.7","58":"0.000000","59":"2.3","60":"0.0","61":"0.0","62":"0.6000000","63":"0.000000","64":"0.0","65":"0.000000","66":"0.0","67":"0.0000000","68":"0.0","69":"0.000000","70":"0.000000","71":"0.0000000","72":"0.0","73":"0.0","74":"0.0","75":"0.0","76":"0.0","77":"0.0","78":"0.0000000","79":"0.0000000","80":"0.0","81":"0.0","82":"0.0","83":"0.0","84":"0.000000","85":"1.2","86":"0.0","87":"0.0","88":"0.0","89":"0.0","90":"0.0","91":"0.0","92":"0.0","93":"5.9","94":"0.0","95":"0.000000","96":"0.0","97":"0.0000000","98":"0.0","99":"0.0","100":"0.0","101":"0.0","102":"0","103":"0.0","104":"0.000000","105":"0.0","106":"0.0","107":"0.0","108":"0.0","109":"0.0","110":"0.0","111":"0.0","112":"0.0","113":"0.000000","114":"0.000000","115":"0.0","116":"0.000000","117":"0.0","118":"0.0","119":"0.4","120":"0.0","121":"0.0","122":"0.0","123":"0.0","124":"0.0","125":"0.0","126":"0.0","127":"0.0","128":"0.0","129":"0.000000","130":"0.0","131":"0.0","132":"0.0","133":"0.0","134":"0.0","135":"0.0","136":"0.0000000","137":"0.0","138":"0.0","139":"0.0","140":"0.0","141":"0","142":"0.0","143":"0.0","144":"0.0000000","145":"0.0","146":"0.0","147":"0.0","148":"0.0","149":"0.0","150":"0.0","151":"0.0","152":"0.0","153":"0.0","154":"0.0","155":"0.0","156":"0.0","157":"0.0","158":"0.0","159":"0.0","160":"0.0","161":"0.0","162":"0.0","163":"0.0","164":"0.0","165":"0.0","166":"0.0","167":"0.0","168":"0.0","169":"0.0","170":"0.0","171":"0.0","172":"0.0","173":"0.0","174":"0.0","175":"0.0","176":"0.0","177":"0.0","178":"0.0000000","179":"0.0","180":"0.0","181":"0.0","182":"0.0","183":"0.0","184":"0.0"},{"1":"PHRP1","2":"339.2000","3":"237.20000","4":"520.70000","5":"27.400000","6":"67.000000","7":"24.8000000","8":"79.50000","9":"1.6","10":"196.20000","11":"1.600000","12":"22.6000000","13":"36.8","14":"42.300000","15":"166.400000","16":"4.700000","17":"92.200000","18":"3.100000","19":"8.500000","20":"0.00000","21":"4.3","22":"78.000000","23":"0.0000000","24":"1.20000","25":"0.0","26":"0.0","27":"4.3","28":"0.0","29":"29.0","30":"97.1","31":"0.0","32":"0.0","33":"0.0","34":"11.300000","35":"2.10000","36":"2.5","37":"0.0000000","38":"2.800000","39":"0.000000","40":"0.0","41":"0.000000","42":"0.0","43":"0.0","44":"0","45":"0.0","46":"0","47":"0.000000","48":"8.6","49":"0.0000000","50":"3.2","51":"0.800000","52":"2.2","53":"0.000000","54":"0.0","55":"0.0","56":"0","57":"0.6","58":"0.000000","59":"0.0","60":"0.0","61":"1.3","62":"0.0000000","63":"0.000000","64":"0.0","65":"0.000000","66":"0.0","67":"0.0000000","68":"0.0","69":"0.000000","70":"0.000000","71":"0.0000000","72":"8.9","73":"0.0","74":"0.0","75":"0.0","76":"0.0","77":"0.0","78":"0.0000000","79":"0.4000000","80":"0.0","81":"0.0","82":"0.0","83":"0.0","84":"0.000000","85":"0.0","86":"0.6","87":"0.0","88":"0.4","89":"0.0","90":"0.0","91":"0.0","92":"0.0","93":"0.0","94":"0.0","95":"0.000000","96":"0.0","97":"0.0000000","98":"0.0","99":"1.4","100":"2.1","101":"0.0","102":"0","103":"0.0","104":"1.100000","105":"0.0","106":"0.0","107":"0.0","108":"0.0","109":"0.0","110":"1.0","111":"0.0","112":"0.7","113":"0.000000","114":"0.000000","115":"0.0","116":"0.000000","117":"0.0","118":"0.0","119":"0.0","120":"0.0","121":"0.0","122":"0.0","123":"0.0","124":"1.2","125":"0.0","126":"0.0","127":"0.0","128":"0.0","129":"0.000000","130":"0.0","131":"0.0","132":"0.0","133":"0.0","134":"0.0","135":"0.0","136":"0.0000000","137":"0.0","138":"0.0","139":"0.0","140":"0.0","141":"0","142":"0.0","143":"0.0","144":"0.0000000","145":"0.0","146":"0.0","147":"0.0","148":"0.0","149":"0.0","150":"0.0","151":"0.0","152":"0.0","153":"0.0","154":"0.0","155":"0.0","156":"0.0","157":"0.0","158":"0.0","159":"0.0","160":"0.0","161":"0.0","162":"0.0","163":"0.0","164":"0.0","165":"0.5","166":"0.0","167":"0.0","168":"0.0","169":"0.0","170":"0.0","171":"0.0","172":"0.0","173":"0.0","174":"0.0","175":"0.0","176":"0.0","177":"0.0","178":"0.0000000","179":"0.0","180":"0.0","181":"0.0","182":"0.0","183":"0.2","184":"0.0"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

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
