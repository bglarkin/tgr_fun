Correlations: guilds and environment
================
Beau Larkin

Last updated: 15 May, 2025

- [Description](#description)
  - [Root path function](#root-path-function)
- [Data](#data)
  - [Site metadata and design](#site-metadata-and-design)
  - [Sites-species tables](#sites-species-tables)
  - [Microbial species metadata](#microbial-species-metadata)
  - [Plant data](#plant-data)
- [Data wrangling](#data-wrangling)
  - [Raw abundances](#raw-abundances)
  - [ALR-transformed abundances](#alr-transformed-abundances)
- [AMF abundance in families](#amf-abundance-in-families)

# Description

Site‑average guild abundances were correlated with plant functional
groups. Sequence data are compositional, so abundances were centered via
the additive‑log‑ratio (ALR) method [Gloor
2017](https://www.frontiersin.org/articles/10.3389/fmicb.2017.02224) and
[Greenacre 2021](https://doi.org/10.3389/fmicb.2021.727398) before
differential tests. Plant cover (10×1 m quadrats, WI only; M. Healy
2016) was analysed without log‑transformation because values are not
instrument‑bounded and may exceed 100% \[Gloor 2017\]. Trait data was
obtained from the [TRY database](https://www.try-db.org)  
([Kattge 2010](https://doi.org/10.1111/j.2041-210X.2010.00067.x),
[Kattge 2011](https://doi.org/10.1111/j.1365-2486.2011.02451.x)). \#
Packages and libraries

``` r
packages_needed <- c(
  "tidyverse", "colorspace", "vegan", "knitr", "conflicted",
  "grid", "gridExtra", "gridtext", "GGally", "geosphere", "ape"
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
  its = read_csv(root_path("clean_data/spe_ITS_metadata.csv"), show_col_types = FALSE),
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

The OTU tables must be wrangled to perform the ALR transformation in
guilds for whole soild fungi and families for AMF. Raw abundances are
kept for plotting. Abundance data are also joined with site and env
paramaters to facilitate downstream analyses.

## Raw abundances

### Whole soil fungi

``` r
its_guab <- 
  spe$its_avg %>% 
  pivot_longer(starts_with("otu"), names_to = "otu_num", values_to = "abund") %>% 
  left_join(spe_meta$its %>% select(otu_num, primary_lifestyle), by = join_by(otu_num)) %>% 
  group_by(field_name, primary_lifestyle) %>% summarize(abund = sum(abund), .groups = "drop") %>% 
  arrange(field_name, -abund) %>% 
  pivot_wider(names_from = "primary_lifestyle", values_from = "abund")
its_guab_pfg <- 
  its_guab %>% 
  left_join(pfg, by = join_by(field_name))
```

### AMF

``` r
amf_fmab <- 
  spe$amf_avg %>% 
  pivot_longer(starts_with("otu"), names_to = "otu_num", values_to = "abund") %>% 
  left_join(spe_meta$amf %>% select(otu_num, family), by = join_by(otu_num)) %>% 
  group_by(field_name, family) %>% summarize(abund = sum(abund), .groups = "drop") %>% 
  arrange(field_name, -abund) %>% 
  pivot_wider(names_from = "family", values_from = "abund") %>% 
  left_join(its_guab %>% select(field_name, unidentified), by = join_by(field_name)) %>% 
  select(field_name, unidentified, everything())
amf_fmab_pfg <- 
  amf_fmab %>% 
  left_join(pfg, by = join_by(field_name))
```

## ALR-transformed abundances

### Whole soil fungi

``` r
its_gulr <- 
  its_guab %>% 
  column_to_rownames(var = "field_name") %>% 
  decostand("alr", MARGIN = 2, reference = 1, pseudocount = 0.2) %>% 
  rownames_to_column(var = "field_name") %>% 
  left_join(sites %>% select(field_name, field_type, yr_since), by = join_by(field_name)) %>% 
  select(field_name, field_type, yr_since, everything()) %>% 
  as_tibble()
its_gulr_pfg <- 
  its_gulr %>% 
  left_join(pfg, by = join_by(field_name))
```

### AMF

``` r
amf_fmlr <- 
  amf_fmab %>% 
  column_to_rownames(var = "field_name") %>% 
  decostand("alr", MARGIN = 2, reference = 1, pseudocount = 0.2) %>% 
  rownames_to_column(var = "field_name") %>% 
  left_join(sites %>% select(field_name, field_type, yr_since), by = join_by(field_name)) %>% 
  select(field_name, field_type, yr_since, everything()) %>% 
  as_tibble()
amf_fmlr_pfg <- 
  amf_fmlr %>% 
  left_join(pfg, by = join_by(field_name))
```

# AMF abundance in families

Display raw abundances in field types and perform tests on
ALR-transformed data

``` r
amf_fmab_ft <- 
  amf_fmab %>% 
  left_join(sites %>% select(field_name, field_type, yr_since), by = join_by(field_name)) %>% 
  pivot_longer(Glomeraceae:Gigasporaceae, names_to = "family", values_to = "abund") %>% 
  group_by(field_type, family) %>% 
  summarize(abund = sum(abund), .groups = "drop") %>% 
  pivot_wider(names_from = field_type, values_from = abund) %>% 
  mutate(total = corn+restored+remnant, across(where(is.numeric), ~ round(.x, 1))) %>% 
  arrange(-total)
kable(amf_fmab_ft, format = "pandoc", caption = "AMF abundance in families and field types")
```

| family               |    corn | restored | remnant |   total |
|:---------------------|--------:|---------:|--------:|--------:|
| Glomeraceae          | 15006.3 |  46787.1 | 11808.2 | 73601.6 |
| Claroideoglomeraceae |   755.8 |   7287.6 |  1469.7 |  9513.1 |
| Paraglomeraceae      |  1654.7 |   3657.9 |   359.3 |  5671.9 |
| Diversisporaceae     |   550.3 |   1307.6 |   244.3 |  2102.2 |
| Gigasporaceae        |    17.5 |    279.3 |    43.6 |   340.4 |

AMF abundance in families and field types

``` r
# Correlate litter and saprotrophs
read_csv(root_path("clean_data", "plant_avg.csv"), show_col_types = FALSE) %>% 
  rename(field_name = SITE) %>% select(LITTER)
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":["LITTER"],"name":[1],"type":["dbl"],"align":["right"]}],"data":[{"1":"8.3"},{"1":"3.2"},{"1":"0.4"},{"1":"2.1"},{"1":"2.5"},{"1":"29.7"},{"1":"0.2"},{"1":"1.7"},{"1":"8.8"},{"1":"4.1"},{"1":"4.3"},{"1":"4.7"},{"1":"3.4"},{"1":"2.3"},{"1":"0.0"},{"1":"0.0"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>
