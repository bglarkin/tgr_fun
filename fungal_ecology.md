Results: Soil Fungal Communities
================
Beau Larkin

Last updated: 14 January, 2026

- [Description](#description)
- [Packages and libraries](#packages-and-libraries)
  - [Root path function](#root-path-function)
- [Data](#data)
  - [Site metadata and design](#site-metadata-and-design)
  - [Fatty Acids: Biomass](#fatty-acids-biomass)
  - [Sites-species tables](#sites-species-tables)
  - [Microbial species metadata](#microbial-species-metadata)
  - [Plant data](#plant-data)
  - [Soil properties](#soil-properties)
- [Species, environment, and
  metadata](#species-environment-and-metadata)
  - [Plant communities](#plant-communities)
  - [Plant functional groups](#plant-functional-groups)
  - [Whole soil fungi](#whole-soil-fungi)
- [Functions](#functions)
- [Whole Soil Fungi](#whole-soil-fungi-1)
  - [Diversity Indices](#diversity-indices)
  - [Abundance (PLFA biomass)](#abundance-plfa-biomass)
  - [Beta Diversity](#beta-diversity)
  - [Constrained Analysis](#constrained-analysis)
  - [Soil fungi and plant community
    correlations](#soil-fungi-and-plant-community-correlations)
- [Arbuscular mycorrhizal fungi](#arbuscular-mycorrhizal-fungi)
  - [Diversity Indices](#diversity-indices-1)
  - [Abundance (NLFA biomass)](#abundance-nlfa-biomass)
  - [Beta Diversity](#beta-diversity-1)
  - [AMF abundance in families](#amf-abundance-in-families)
  - [Indicator species analysis](#indicator-species-analysis)
  - [Constrained Analysis](#constrained-analysis-1)
  - [AM fungi and plant community
    correlations](#am-fungi-and-plant-community-correlations)
- [Putative plant pathogens](#putative-plant-pathogens)
  - [Diversity Indices](#diversity-indices-2)
  - [Abundance](#abundance)
  - [Beta Diversity](#beta-diversity-2)
  - [Pathogen Indicator Species](#pathogen-indicator-species)
  - [Constrained Analysis](#constrained-analysis-2)
  - [Pathogen and plant community
    correlations](#pathogen-and-plant-community-correlations)
- [Putative saprotrophs](#putative-saprotrophs)
  - [Diversity Indices](#diversity-indices-3)
  - [Abundance](#abundance-1)
  - [Beta Diversity](#beta-diversity-3)
  - [Saprotroph Indicator Species](#saprotroph-indicator-species)
  - [Constrained Analysis](#constrained-analysis-3)
  - [Saprotroph and plant community
    correlations](#saprotroph-and-plant-community-correlations)
- [Summary results](#summary-results)
  - [Richness summary](#richness-summary)
  - [Shannon diversity summary](#shannon-diversity-summary)
  - [Biomass (abundance) summary](#biomass-abundance-summary)
  - [Beta diversity summary](#beta-diversity-summary)
  - [Constrained analysis summary](#constrained-analysis-summary)

# Description

**Scope** – Biomass (PLFA/NLFA), OTU richness and diversity, and
β‑diversity of soil fungi across corn, restored, and remnant prairie
fields.

**Alpha diversity** – 97 %-OTUs (ITS & 18S); site means are replicates;
means separation model selection based on response and residuals
distributions; sequencing depth used as covariate per [Bálint
2015](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.13018) when
warranted; pairwise LSMs via *emmeans*.

**Beta diversity** – Workflow after [Song
2015](https://doi.org/10.1371/journal.pone.0127234):  
1. PCoA of Bray (ITS) or UNIFRAC (18S) distances 1. homogeneity test
diagnostics 1. PERMANOVA (+ pairwise)

Cartesian inter‑site distance enters models as a covariate per [Redondo
2020](https://doi.org/10.1093/femsec/fiaa082).

# Packages and libraries

``` r
# Libraries ———————— ####
```

``` r
packages_needed <- c(
  # Analysis
  "emmeans", "vegan", "phyloseq", "ape", "phangorn", "geosphere", 
  "car", "rlang", "rsq", "sandwich", "lmtest", "performance", "boot", 
  "indicspecies", "MASS", "DHARMa", "broom",
  # Scripting
  "rprojroot", "conflicted", "purrr", "knitr", "tidyverse", 
  # Graphics
  "colorspace", "grid", "gridExtra", "ggpubr", "patchwork"
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
conflicts_prefer(
  dplyr::filter(),
  dplyr::select(),
  dplyr::where(),
  vegan::diversity(),
  purrr::map()
)
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
sites <- read_csv(root_path("clean_data/sites.csv"), show_col_types = FALSE) %>% 
  mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant")))
```

### Wrangle site metadata

Intersite geographic distance will be used as a covariate in clustering.
Raw coordinates in data file aren’t distances; convert to distance
matrix and summarize with PCoA

``` r
field_dist <- as.dist(distm(sites[, c("long", "lat")], fun = distHaversine))
field_dist_pcoa <- pcoa(field_dist)
field_dist_pcoa$values[c(1,2), c(1,2)] %>% 
  kable(format = "pandoc")
```

|  Eigenvalues | Relative_eig |
|-------------:|-------------:|
| 146898426293 |    0.9053961 |
|  15349390146 |    0.0946047 |

First axis of geographic distance PCoA explains 91% of the variation
among sites.

``` r
sites$dist_axis_1 <- field_dist_pcoa$vectors[, 1]
```

## Fatty Acids: Biomass

Use only 18.2 for soil fungi

``` r
fa <- read_csv(root_path("clean_data/plfa.csv"), show_col_types = FALSE) %>% 
  rename(fungi_18.2 = fa_18.2) %>% 
  select(field_name, fungi_18.2, amf) %>%
  left_join(
    sites %>% select(field_name, field_type),
    by = join_by(field_name)
  )
```

## Sites-species tables

CSV files were produced in `sequence_data.R`. Average sequence abundance
at sites included here. Amf_avg_uni table is in species-samples format
to enable use of `Unifrac()` later.

``` r
its_avg     = read_csv(root_path("clean_data/spe_ITS_avg.csv"), show_col_types = FALSE)
amf_avg     = read_csv(root_path("clean_data/spe_18S_avg.csv"), show_col_types = FALSE)
```

## Microbial species metadata

``` r
its_meta = read_csv(root_path("clean_data/spe_ITS_metadata.csv"), show_col_types = FALSE) %>% 
  mutate(primary_lifestyle = case_when(str_detect(primary_lifestyle, "_saprotroph$") ~ "saprotroph",
                                       TRUE ~ primary_lifestyle),
         across(everything(), ~ replace_na(., "unidentified")))
amf_meta = read_csv(root_path("clean_data/spe_18S_metadata.csv"), show_col_types = FALSE) %>% 
  mutate(across(everything(), ~ replace_na(., "unidentified")))
```

### Wrangle additional species and metadata objects

``` r
# Species data wrangling ———————— ####
```

1.  Proportional species abundance corrected for site biomass
2.  Unifrac products for AMF
3.  Phyloseq products to process the Unifrac distance cleanly

``` r
its_avg_ma <- its_avg %>% # ma = sequence proportion of biomass
  rowwise() %>%
  mutate(total = sum(c_across(where(is.numeric))),
         across(starts_with("otu"), ~ if_else(total > 0, .x / total, 0))) %>%
  left_join(fa %>% select(-amf, -field_type), by = join_by(field_name)) %>%
  mutate(across(starts_with("otu"), ~ .x * fungi_18.2)) %>%
  select(field_name, starts_with("otu")) %>% 
  ungroup()
amf_avg_ma <- amf_avg %>% 
  rowwise() %>%
  mutate(total = sum(c_across(where(is.numeric))),
         across(starts_with("otu"), ~ if_else(total > 0, .x / total, 0))) %>%
  left_join(fa %>% select(-fungi_18.2, -field_type), by = join_by(field_name)) %>%
  mutate(across(starts_with("otu"), ~ .x * amf)) %>%
  select(field_name, starts_with("otu")) %>% 
  ungroup()
amf_avg_uni <- amf_avg %>%
  column_to_rownames("field_name") %>%
  t() %>% as.data.frame() %>% rownames_to_column("otu_num") %>%
  left_join(amf_meta %>% select(otu_num, otu_ID), by = "otu_num") %>%
  select(otu_ID, everything(), -otu_num) %>% 
  as_tibble()
amf_ps <- phyloseq(
  otu_table(data.frame(amf_avg_uni, row.names = 1), taxa_are_rows = TRUE),
  tax_table(as.matrix(data.frame(amf_meta, row.names = 2))),
  read.dna(root_path("otu_tables/18S/18S_sequences.fasta"), format = "fasta") %>%
    phyDat(type = "DNA") %>% dist.hamming() %>% NJ(),
  sample_data(sites %>% column_to_rownames(var = "field_name"))
)
```

## Plant data

Abundance in functional groups and by species are only available from
Wisconsin sites. Only C4_grass and forbs are used. Others: C3_grass,
legume, and shrubTree were found previously to have high VIF in models
or were not chosen in forward selection.

``` r
pfg <- read_csv(root_path("clean_data", "plant_traits.csv"), show_col_types = FALSE) 
```

## Soil properties

``` r
soil <- read_csv(root_path("clean_data/soil.csv"), show_col_types = FALSE)[-c(26:27), ]
```

# Species, environment, and metadata

``` r
# Metadata wrangling ———————— ####
```

## Plant communities

### Plant species richness

``` r
plant <- read_csv(root_path("clean_data/plant_avg.csv"), show_col_types = FALSE)
prich <- plant %>% 
  select(-BARESOIL, -LITTER, -ROSA, -SALIX) %>% # remove non-species entries
  rowwise() %>% 
  mutate(pl_rich = sum(c_across(where(is.numeric)) > 0),
         pl_shan = exp(diversity(c_across(where(is.numeric))))
         ) %>% 
  select(field_name = SITE, pl_rich, pl_shan) %>% 
  left_join(sites, by = join_by(field_name)) %>% 
  filter(field_type != "corn") %>% 
  ungroup()
with(prich[-c(3,6,9), ], cor.test(yr_since, pl_rich)) # Remnant fields don't have an age
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  yr_since and pl_rich
    ## t = -1.2324, df = 8, p-value = 0.2528
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.8222704  0.3075217
    ## sample estimates:
    ##       cor 
    ## -0.399447

Years since restoration isn’t obviously related to plant species
richness.

## Plant functional groups

C4 grass and forb cover are transformed into a single index using PCA in
restored sites only. Abundance data are also joined with site and env
paramaters to facilitate downstream analyses.

### Grass-forb index

C4 grass and forb cover are highly correlated (*r* = -0.91) in restored
prairies. In models or constrained ordinations, they are collinear and
cannot be used simultaneously. An index of grass-forb cover is created
to solve this problem.

``` r
pfg_pca <- 
  pfg %>%
  select(field_name, C3_grass:shrubTree) %>%
  rowwise() %>%
  mutate(total = sum(c_across(where(is.numeric))),
         across(C3_grass:shrubTree, ~ if_else(total > 0, .x / total, 0))) %>%
  ungroup() %>%
  select(field_name, C4_grass, forb) %>% 
  left_join(sites %>% select(field_name, field_type), by = join_by(field_name)) %>% 
  filter(field_type != "corn") %>% 
  select(-field_type) %>% 
  column_to_rownames(var = "field_name") %>% 
  rda()
pfg_pca %>% summary() # 92% variation on first axis
```

    ## 
    ## Call:
    ## rda(X = .) 
    ## 
    ## Partitioning of variance:
    ##               Inertia Proportion
    ## Total         0.07131          1
    ## Unconstrained 0.07131          1
    ## 
    ## Eigenvalues, and their contribution to the variance 
    ## 
    ## Importance of components:
    ##                           PC1     PC2
    ## Eigenvalue            0.06655 0.00476
    ## Proportion Explained  0.93324 0.06676
    ## Cumulative Proportion 0.93324 1.00000

Define the grass_forb index

``` r
gf_index = scores(pfg_pca, choices = 1, display = "sites") %>% 
  data.frame() %>% 
  rename(gf_index = PC1) %>% 
  rownames_to_column(var = "field_name")
```

Are field age and gf_index correlated?

``` r
gfi_yrs <- gf_index %>% 
  left_join(sites %>% select(field_name, yr_since), by = join_by(field_name)) %>% 
  arrange(-gf_index)
with(gfi_yrs, cor.test(yr_since, gf_index))
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  yr_since and gf_index
    ## t = -7.1535, df = 8, p-value = 9.676e-05
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.9836360 -0.7245743
    ## sample estimates:
    ##       cor 
    ## -0.929948

The relatively strong correlation suggests that different restoration
methods over time are still reflected in plant composition. Years since
restoration is highly related to plant community change.

Visualize grass forb gradient compared with plant composition, grass and
forb cover, and years since restoration.

``` r
fs8_xlab <- pfg %>% 
  left_join(sites, by = join_by(field_name)) %>% 
  filter(field_type != "corn") %>% 
  group_by(field_type) %>% 
  mutate(xlab = case_when(field_type == "restored" ~ paste("RS", 1:n(), sep = "-"),
                          field_type == "remnant"  ~ paste("RM", 1:n(), sep = "-"))
         ) %>% 
  ungroup() %>% 
  select(field_name, xlab)
plt_div <- 
  prich %>% 
  left_join(gf_index, by = join_by(field_name)) %>% 
  left_join(fs8_xlab, by = join_by(field_name)) %>% 
  select(field_name, xlab, gf_index, pl_rich, pl_shan) %>%
  pivot_longer(pl_rich:pl_shan, names_to = "var", values_to = "value") %>% 
  ggplot(aes(x = fct_reorder(xlab, gf_index), y = value, group = var)) +
  geom_col(aes(fill = var), position = position_dodge()) +
  labs(x = NULL, y = expression(atop("Alpha diversity", paste("(", italic(n), " species)")))) +
  scale_fill_discrete_qualitative(name = "Diversity index", palette = "Dynamic", 
                                  labels = c(expression("richness"), expression(paste("Shannon (", italic(e)^{italic(H)*"\u2032"}, ")")))) +
  theme_cor +
  theme(plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
pfg_comp <- 
  pfg %>% 
  select(field_name, C3_grass:shrubTree) %>% 
  rowwise() %>% 
  mutate((across(where(is.numeric), ~ .x / sum(c_across(where(is.numeric))))) * 100) %>% 
  ungroup() %>% 
  pivot_longer(C3_grass:shrubTree, names_to = "pfg", values_to = "pct_comp") %>% 
  left_join(sites, by = join_by(field_name)) %>% 
  left_join(gf_index, by = join_by(field_name)) %>% 
  filter(field_type != "corn") %>% 
  select(field_name, yr_since, gf_index, pfg, pct_comp) %>% 
  mutate(pfg = factor(pfg, levels = c("shrubTree", "legume", "C3_grass", "C4_grass", "forb"),
                      labels = c("shrub, tree", "legume", "grass (C3)", "grass (C4)", "forb"))) %>% 
  left_join(fs8_xlab, by = join_by(field_name))
pfg_comp_fig <- 
  ggplot(pfg_comp, aes(x = fct_reorder(xlab, gf_index), y = pct_comp, group = pfg)) +
  geom_col(aes(fill = pfg)) +
  labs(x = NULL, y = "Composition (%)") +
  scale_fill_manual(name = "Functional group", values = pfg_col,
                    labels = c(expression("shrub, tree"), expression("legume"), 
                               expression("grass ("*C[3]*")"), expression("grass ("*C[4]*")"), expression("forb"))) +
  theme_cor +
  theme(plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
pfg_pct <- 
  pfg %>% 
  select(field_name, C4_grass, forb) %>%
  pivot_longer(C4_grass:forb, names_to = "pfg", values_to = "pct_cvr") %>% 
  left_join(sites, by = join_by(field_name)) %>% 
  left_join(gf_index, by = join_by(field_name)) %>% 
  filter(field_type != "corn") %>% 
  select(field_name, yr_since, gf_index, pfg, pct_cvr)  %>% 
  mutate(pfg = factor(pfg, levels = c("C4_grass", "forb"),
                      labels = c("grass (C4)", "forb"))) %>% 
  left_join(fs8_xlab, by = join_by(field_name))
gf_pct_fig <- 
  ggplot(pfg_pct, aes(x = fct_reorder(xlab, gf_index), y = pct_cvr, group = pfg)) +
  geom_step(aes(color = pfg), linejoin = "round", lineend = "round") +
  geom_point(aes(color = pfg), shape = 21, size = 1.8, fill = "white", stroke = 0.9) +
  scale_color_manual(name = "Functional group", values = pfg_col[4:5], 
                     labels = c(expression("grass ("*C[4]*")"), expression("forb"))) +
  labs(x = NULL, y = "Cover (%)") +
  theme_cor +
  theme(plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
gfi_loc_fig <- 
  gfi_yrs %>% 
  left_join(sites, by = join_by(field_name)) %>% 
  ggplot(aes(x = gf_index, y = rep("PCA 1", nrow(gfi_yrs)))) + 
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.3, color = "gray20") +
  geom_point(aes(color = field_type), shape = 21, size = 1.8, fill = "white", stroke = 0.9) +
  labs(y = NULL, x = "Grass-forb index") +
  scale_color_manual(name = "Field type", values = ft_pal[2:3]) +
  theme_cor +
  theme(plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1.1),
        axis.text.y = element_text(hjust = 0))
```

#### Unified figure

``` r
pfg_pct_fig <- (plt_div / plot_spacer() / pfg_comp_fig / plot_spacer() / gf_pct_fig / plot_spacer() / gfi_loc_fig) +
  plot_layout(heights = c(1,0.01,1,0.01,0.7,0.01,0.2)) +
  plot_annotation(tag_levels = 'A') 
```

``` r
pfg_pct_fig
```

![](resources/fungal_ecology_files/figure-gfm/pfg_fig-1.png)<!-- -->

## Whole soil fungi

Wrangle data to produce proportional biomass in guilds for its and
families for amf

``` r
its_guild_ma <- # guild biomass (proportion of total biomass)
  its_avg_ma %>% 
  pivot_longer(starts_with("otu"), names_to = "otu_num", values_to = "abund") %>% 
  left_join(its_meta %>% select(otu_num, primary_lifestyle), by = join_by(otu_num)) %>% 
  group_by(field_name, primary_lifestyle) %>% summarize(abund = sum(abund), .groups = "drop") %>% 
  arrange(field_name, -abund) %>% 
  pivot_wider(names_from = "primary_lifestyle", values_from = "abund") %>% 
  select(field_name, patho_mass = plant_pathogen, sapro_mass = saprotroph) %>% 
  left_join(pfg, by = join_by(field_name)) %>% 
  left_join(gf_index, by = join_by(field_name)) %>% 
  left_join(sites %>% select(field_name, field_type, region, yr_since), by = join_by(field_name)) %>% 
  select(field_name, field_type, yr_since, region, everything())
```

Wrangle a second set to compare raw sequence abundances and proportion
of biomass values together

``` r
its_guild <- 
  its_avg %>% 
  pivot_longer(starts_with("otu"), names_to = "otu_num", values_to = "abund") %>% 
  left_join(its_meta %>% select(otu_num, primary_lifestyle), by = join_by(otu_num)) %>% 
  group_by(field_name, primary_lifestyle) %>% summarize(abund = sum(abund), .groups = "drop") %>% 
  arrange(field_name, -abund) %>% 
  pivot_wider(names_from = "primary_lifestyle", values_from = "abund") %>% 
  rowwise() %>% 
  mutate(fungi_abund = sum(c_across(where(is.numeric)))) %>% 
  select(field_name, patho_abund = plant_pathogen, sapro_abund = saprotroph, fungi_abund) %>% 
  left_join(fa %>% select(field_name, fungi_mass = fungi_18.2), by = join_by(field_name)) %>% 
  left_join(pfg, by = join_by(field_name)) %>% 
  left_join(gf_index, by = join_by(field_name)) %>% 
  left_join(sites %>% select(field_name, field_type, region, yr_since), by = join_by(field_name)) %>% 
  select(field_name, field_type, yr_since, region, everything()) %>% 
  ungroup()
```

#### AMF

``` r
amf_fam <- # sequence abundance in families
  amf_avg %>% 
  pivot_longer(starts_with("otu"), names_to = "otu_num", values_to = "abund") %>% 
  left_join(amf_meta %>% select(otu_num, family), by = join_by(otu_num)) %>% 
  group_by(field_name, family) %>% summarize(abund = sum(abund), .groups = "drop") %>% 
  arrange(field_name, -abund) %>% 
  pivot_wider(names_from = "family", values_from = "abund") %>% 
  rename_with(~ paste0(abbreviate(.x, minlength = 5, strict = TRUE), "_ab"),
              Glomeraceae:Ambisporaceae) %>% 
  left_join(pfg %>% select(field_name, C3_grass:shrubTree), by = join_by(field_name)) %>% 
  left_join(gf_index, by = join_by(field_name)) %>% 
  left_join(sites %>% select(field_name, field_type, region, yr_since), by = join_by(field_name)) %>% 
  select(field_name, field_type, yr_since, region, everything())
amf_fam_ma <- # family biomass (proportion of total biomass)
  amf_avg_ma %>% 
  pivot_longer(starts_with("otu"), names_to = "otu_num", values_to = "abund") %>% 
  left_join(amf_meta %>% select(otu_num, family), by = join_by(otu_num)) %>% 
  group_by(field_name, family) %>% summarize(abund = sum(abund), .groups = "drop") %>% 
  arrange(field_name, -abund) %>% 
  pivot_wider(names_from = "family", values_from = "abund") %>% 
  rename_with(~ paste0(abbreviate(.x, minlength = 5, strict = TRUE), "_mass"),
              Glomeraceae:Ambisporaceae) %>% 
  left_join(pfg %>% select(field_name, C3_grass:shrubTree), by = join_by(field_name)) %>% 
  left_join(gf_index, by = join_by(field_name)) %>% 
  left_join(sites %>% select(field_name, field_type, region, yr_since), by = join_by(field_name)) %>% 
  select(field_name, field_type, yr_since, region, everything())
```

# Functions

Executed from a separate script to save lines here; to view the function
navigate to `functions.R` in the code folder, accessible from the root
dir of the repo.

``` r
# Functions ———————— ####
source(root_path("code", "functions.R"))
```

# Whole Soil Fungi

``` r
# Whole soil fungi ———————— ####
```

## Diversity Indices

``` r
its_div <- calc_div(its_avg, sites) %>% 
  mutate(depth_csq = sqrt(depth) - mean(sqrt(depth)))
```

### Richness

Sequence depth square root transformed and centered. Negative binomial
model used to handle count data. Poisson model was overdispersed (not
shown).

Test interaction

``` r
its_rich_glm_i <- glm.nb(richness ~ depth_csq * field_type, data = its_div)
Anova(its_rich_glm_i, type = 3, test.statistic = "LR") # no interaction detected
```

    ## Analysis of Deviance Table (Type III tests)
    ## 
    ## Response: richness
    ##                      LR Chisq Df Pr(>Chisq)    
    ## depth_csq               0.383  1     0.5359    
    ## field_type             41.979  2  7.661e-10 ***
    ## depth_csq:field_type    2.890  2     0.2358    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Fit additive model

``` r
its_rich_glm <- glm.nb(richness ~ depth_csq + field_type, data = its_div)
```

Diagnostics

``` r
check_model(its_rich_glm)
```

![](resources/fungal_ecology_files/figure-gfm/its_rich_covar_diagnostics-1.png)<!-- -->

``` r
check_overdispersion(its_rich_glm)
```

    ## # Overdispersion test
    ## 
    ##  dispersion ratio = 1.047
    ##           p-value = 0.792

    ## No overdispersion detected.

``` r
check_collinearity(its_rich_glm)
```

    ## # Check for Multicollinearity
    ## 
    ## Low Correlation
    ## 
    ##        Term  VIF    VIF 95% CI adj. VIF Tolerance Tolerance 95% CI
    ##   depth_csq 1.06 [1.00, 15.14]     1.03      0.94     [0.07, 1.00]
    ##  field_type 1.06 [1.00, 15.14]     1.01      0.94     [0.07, 1.00]

Long tails, some midrange structure, no leverage points

``` r
distribution_prob(its_rich_glm)
```

    ## 
    ## 
    ## Distribution    p_Residuals
    ## -------------  ------------
    ## normal               0.5625
    ## cauchy               0.1875
    ## gamma                0.1250
    ## 
    ## 
    ## Distribution                  p_Response
    ## ---------------------------  -----------
    ## lognormal                         0.3750
    ## neg. binomial (zero-infl.)        0.3125
    ## beta-binomial                     0.1250

residuals distribution normal or long-tailed, response log

``` r
leveneTest(richness ~ field_type, data = its_div) %>% as.data.frame() %>% kable(format = "pandoc")
```

|       |  Df |   F value |   Pr(\>F) |
|-------|----:|----------:|----------:|
| group |   2 | 0.4909782 | 0.6185753 |
|       |  22 |        NA |        NA |

``` r
leveneTest(residuals(its_rich_glm) ~ its_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc")
```

|       |  Df |   F value |   Pr(\>F) |
|-------|----:|----------:|----------:|
| group |   2 | 0.1626475 | 0.8509033 |
|       |  22 |        NA |        NA |

Residuals/response distributions do not suggest the need for
transformation. Levene’s p \> 0.05 → fail to reject = variances can be
considered equal across groups.

Model results, group means, and post-hoc. Use Type II LR test of
variables due to unbalanced design.

``` r
Anova(its_rich_glm, type = 2, test.statistic = "LR")
```

    ## Analysis of Deviance Table (Type II tests)
    ## 
    ## Response: richness
    ##            LR Chisq Df Pr(>Chisq)    
    ## depth_csq     8.128  1   0.004359 ** 
    ## field_type   40.506  2  1.601e-09 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Sequence depth is significant, less so than field type. Proceed with
means separation by obtaining estimated marginal means for field type.

``` r
its_rich_em <- emmeans(its_rich_glm, ~ field_type, type = "response")
```

Results tables below show the emmeans summary of group means and
confidence intervals, with sequencing depth as a covariate, and the post
hoc contrast of richness among field types.

| field_type | response |       SE |  df | asymp.LCL | asymp.UCL |
|:-----------|---------:|---------:|----:|----------:|----------:|
| corn       | 391.0204 | 15.42835 | Inf |  361.9210 |  422.4594 |
| restored   | 499.1467 | 10.64842 | Inf |  478.7065 |  520.4597 |
| remnant    | 553.8016 | 23.85811 | Inf |  508.9603 |  602.5936 |

Confidence level used: 0.95

| contrast           |     ratio |        SE |  df | null |   z.ratio |   p.value |
|:-------------------|----------:|----------:|----:|-----:|----------:|----------:|
| corn / restored    | 0.7833777 | 0.0350650 | Inf |    1 | -5.454267 | 0.0000001 |
| corn / remnant     | 0.7060658 | 0.0414771 | Inf |    1 | -5.924805 | 0.0000000 |
| restored / remnant | 0.9013096 | 0.0436858 | Inf |    1 | -2.143762 | 0.0812207 |

P value adjustment: tukey method for comparing a family of 3 estimates

OTU richness in cornfields is significantly less than in restored or
remnant fields (p\<0.001), which don’t differ.

``` r
its_rich_fig <- 
  ggplot(summary(its_rich_em), aes(x = field_type, y = response)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = response, ymax = asymp.UCL), width = 0, linewidth = lw) +
  geom_text(aes(y = asymp.UCL, label = c("a", "b", "b")),  vjust = -1, family = "sans", size = 3.5) +
  labs(x = NULL, y = expression(atop("Richness", paste("(", italic(n), " OTUs)")))) +
  lims(y = c(0, 760)) +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1))
```

### Shannon’s diversity

Sequence depth square root transformed and centered

``` r
its_shan_lm <- lm(shannon ~ depth_csq + field_type, data = its_div)
```

Diagnostics

``` r
check_model(its_shan_lm)
```

![](resources/fungal_ecology_files/figure-gfm/its_shan_covar_diagnostics-1.png)<!-- -->

Some residual structure, no leverage points, no evidence for increasing
mean/var relationship.

``` r
distribution_prob(its_shan_lm)
```

    ## 
    ## 
    ## Distribution    p_Residuals
    ## -------------  ------------
    ## cauchy               0.5625
    ## normal               0.4375
    ## bernoulli            0.0000
    ## 
    ## 
    ## Distribution    p_Response
    ## -------------  -----------
    ## gamma              0.31250
    ## lognormal          0.31250
    ## chi                0.15625

residuals distribution most likely cauchy/normal; symmetric but long
tails, response gamma

``` r
leveneTest(shannon ~ field_type, data = its_div) %>% as.data.frame() %>% kable(format = "pandoc")
```

|       |  Df |  F value |   Pr(\>F) |
|-------|----:|---------:|----------:|
| group |   2 | 2.742412 | 0.0864225 |
|       |  22 |       NA |        NA |

``` r
leveneTest(residuals(its_shan_lm) ~ its_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc")
```

|       |  Df |  F value |   Pr(\>F) |
|-------|----:|---------:|----------:|
| group |   2 | 2.265973 | 0.1274059 |
|       |  22 |       NA |        NA |

Residuals distribution does not suggest the need for transformation.
Levene’s p \> 0.05 → fail to reject = variances can be considered equal.
Response more suspicious. Examine CV in groups to assess changes in
variance.

``` r
augment(its_shan_lm) %>%
  mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant"))) %>% 
  group_by(field_type) %>%
  summarise(
    mean_fitted = mean(.fitted),
    sd_resid    = sd(.resid),
    cv_resid    = sd_resid / mean_fitted
  ) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>% 
  kable(format = "pandoc", caption = "CV of residuals and fitted means in groups")
```

| field_type | mean_fitted | sd_resid | cv_resid |
|:-----------|------------:|---------:|---------:|
| corn       |       79.92 |    15.20 |     0.19 |
| restored   |      111.36 |    21.39 |     0.19 |
| remnant    |      120.89 |     7.76 |     0.06 |

CV of residuals and fitted means in groups

Residuals’ CV constant to declining. Relatively low Levene’s p value
likely due to unequal variance in restored and remnant despite similar
means. Unbalanced data and possible biological reality likely causing
this. No need for further transformation.

Model results, group means, and post-hoc. Type II SS used due to
unbalanced design.

``` r
Anova(its_shan_lm, type = 2)
```

    ## Anova Table (Type II tests)
    ## 
    ## Response: shannon
    ##            Sum Sq Df F value   Pr(>F)   
    ## depth_csq   328.6  1  0.8662 0.362577   
    ## field_type 4968.9  2  6.5489 0.006162 **
    ## Residuals  7966.9 21                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Sequence depth is not a significant predictor of Shannon diversity.
Proceed with means separation by obtaining estimated marginal means for
field type. Arithmetic means calculated in this case.

``` r
its_shan_em <- emmeans(its_shan_lm, ~ field_type, type = "response")
```

Results tables below show the emmeans summary of group means and
confidence intervals, with sequencing depth as a covariate, and the post
hoc contrast of richness among field types.

| field_type |    emmean |       SE |  df |  lower.CL |  upper.CL |
|:-----------|----------:|---------:|----:|----------:|----------:|
| corn       |  79.49764 | 8.722310 |  21 |  61.35861 |  97.63668 |
| restored   | 110.98223 | 4.886025 |  21 | 100.82118 | 121.14327 |
| remnant    | 122.91831 | 9.978907 |  21 | 102.16604 | 143.67059 |

Confidence level used: 0.95

| contrast           |  estimate |        SE |  df |   t.ratio |   p.value |
|:-------------------|----------:|----------:|----:|----------:|----------:|
| corn - restored    | -31.48458 |  9.979385 |  21 | -3.154962 | 0.0127116 |
| corn - remnant     | -43.42067 | 13.327491 |  21 | -3.257978 | 0.0100752 |
| restored - remnant | -11.93609 | 11.189520 |  21 | -1.066720 | 0.5445052 |

P value adjustment: tukey method for comparing a family of 3 estimates

Shannon’s diversity in cornfields is significantly less than in restored
or remnant fields, which don’t differ.

``` r
its_shan_fig <- 
  ggplot(summary(its_shan_em), aes(x = field_type, y = emmean)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = emmean, ymax = upper.CL), width = 0, linewidth = lw) +
  geom_text(aes(y = upper.CL, label = c("a", "b", "b")), vjust = -1, family = "sans", size = 3.5) +
  labs(x = "Field type", y = expression(Shannon~diversity~paste("(", italic(e)^italic(H), ")"))) +
  lims(y = c(0, 160)) +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
```

## Abundance (PLFA biomass)

``` r
plfa_lm <- lm(fungi_18.2 ~ field_type, data = fa)
par(mfrow = c(2,2))
plot(plfa_lm) 
```

![](resources/fungal_ecology_files/figure-gfm/unnamed-chunk-37-1.png)<!-- -->

variance differs slightly in groups. Tails on qq plot diverge, lots of
groups structure visible.

``` r
distribution_prob(plfa_lm)
```

    ## 
    ## 
    ## Distribution    p_Residuals
    ## -------------  ------------
    ## normal              0.62500
    ## cauchy              0.12500
    ## gamma               0.09375
    ## 
    ## 
    ## Distribution     p_Response
    ## --------------  -----------
    ## weibull             0.21875
    ## normal              0.18750
    ## beta-binomial       0.12500

Residuals distribution fits normal, response normal-ish

``` r
leveneTest(residuals(plfa_lm) ~ fa$field_type) %>% as.data.frame() %>% kable(format = "pandoc") # No covariate, response and residuals tests equivalent
```

|       |  Df |   F value |   Pr(\>F) |
|-------|----:|----------:|----------:|
| group |   2 | 0.9749963 | 0.3929075 |
|       |  22 |        NA |        NA |

Residuals distribution does not suggest the need for transformation.
Levene’s p \> 0.05 → fail to reject = variances can be considered equal.

Model results, group means, and post-hoc, with arithmetic means from
emmeans

``` r
anova(plfa_lm)
```

    ## Analysis of Variance Table
    ## 
    ## Response: fungi_18.2
    ##            Df Sum Sq Mean Sq F value Pr(>F)
    ## field_type  2 16.229  8.1146  2.3426 0.1196
    ## Residuals  22 76.207  3.4639

``` r
plfa_em <- emmeans(plfa_lm, ~ field_type, type = "response")
```

| field_type |   emmean |        SE |  df | lower.CL | upper.CL |
|:-----------|---------:|----------:|----:|---------:|---------:|
| corn       | 3.094661 | 0.8323383 |  22 | 1.368497 | 4.820825 |
| restored   | 5.129779 | 0.4652913 |  22 | 4.164824 | 6.094734 |
| remnant    | 5.011704 | 0.9305825 |  22 | 3.081794 | 6.941614 |

Confidence level used: 0.95

| contrast           |   estimate |        SE |  df |    t.ratio |   p.value |
|:-------------------|-----------:|----------:|----:|-----------:|----------:|
| corn - restored    | -2.0351183 | 0.9535633 |  22 | -2.1342246 | 0.1058898 |
| corn - remnant     | -1.9170436 | 1.2485075 |  22 | -1.5354682 | 0.2943078 |
| restored - remnant |  0.1180747 | 1.0404229 |  22 |  0.1134873 | 0.9929269 |

P value adjustment: tukey method for comparing a family of 3 estimates

``` r
plfa_fig <- 
  ggplot(summary(plfa_em), aes(x = field_type, y = emmean)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = emmean, ymax = upper.CL), width = 0, linewidth = lw) +
  labs(x = "Field type", y = expression(atop("Biomass", paste("(", nmol[PLFA], " × ", g[soil]^{-1}, ")")))) +
  # labs(x = "Field Type", y = "Biomass") +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1.1))
```

## Beta Diversity

1.  Using biomass-weighted relative abundance [Waller et
    al. 2023](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/1365-2745.14392)
2.  Using sequence-based relative abundance [McKnight et
    al. 2019](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13115)
3.  Contrast the two with procrustes

### Biomass-weighted relative abundance

``` r
d_its_ma <- its_avg_ma %>% 
  data.frame(row.names = 1) %>% 
  vegdist("bray")
mva_its_ma <- mva(d = d_its_ma, env = sites)
```

``` r
mva_its_ma$dispersion_test
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
    ## Groups     2 0.014652 0.0073262 2.3478   1999  0.116
    ## Residuals 22 0.068649 0.0031204                     
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##             corn remnant restored
    ## corn             0.18950   0.1705
    ## remnant  0.18914           0.1220
    ## restored 0.16490 0.11917

``` r
mva_its_ma$permanova
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## adonis2(formula = d ~ dist_axis_1 + field_type, data = env, permutations = nperm, by = "terms")
    ##             Df SumOfSqs      R2      F Pr(>F)    
    ## dist_axis_1  1   0.4368 0.06225 1.7169 0.0270 *  
    ## field_type   2   1.2378 0.17639 2.4327 0.0005 ***
    ## Residual    21   5.3426 0.76136                  
    ## Total       24   7.0171 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
mva_its_ma$pairwise_contrasts[c(1,3,2), c(1,2,4,3,8)] %>% 
  arrange(group1, desc(group2)) %>% 
  kable(format = "pandoc", caption = "Pairwise permanova contrasts")
```

| group1  | group2   | F_value |    R2 | p_value_adj |
|:--------|:---------|--------:|------:|------------:|
| corn    | restored |   3.809 | 0.160 |      0.0015 |
| corn    | remnant  |   2.662 | 0.268 |      0.0030 |
| remnant | restored |   0.980 | 0.050 |      0.4530 |

Pairwise permanova contrasts

No eignevalue correction was needed. Two relative eigenvalues exceeded
broken stick model. Based on the homogeneity of variance test, the null
hypothesis of equal variance among groups is accepted across all
clusters and in pairwise comparison of clusters (both p\>0.05),
supporting the application of a PERMANOVA test.

Clustering revealed that community variation was related to geographic
distance, the covariate in the model. With geographic distance accounted
for, the test variable ‘field type’ significantly explained variation in
fungal communities, with a post-hoc test revealing that communities in
corn fields differed from communities in restored and remnant fields.

Plotting results:

``` r
its_ma_ord_data <- mva_its_ma$ordination_scores %>% mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant")))
p_its_ma_centers <- its_ma_ord_data %>% 
  group_by(field_type) %>% 
  summarize(across(starts_with("Axis"), list(mean = mean, ci_l = ci_l, ci_u = ci_u), .names = "{.fn}_{.col}"), .groups = "drop") %>% 
  mutate(across(c(ci_l_Axis.1, ci_u_Axis.1), ~ mean_Axis.1 + .x),
         across(c(ci_l_Axis.2, ci_u_Axis.2), ~ mean_Axis.2 + .x))
its_ma_ord <- 
  ggplot(its_ma_ord_data, aes(x = Axis.1, y = Axis.2)) +
  geom_linerange(data = p_its_ma_centers, aes(x = mean_Axis.1, y = mean_Axis.2, xmin = ci_l_Axis.1, xmax = ci_u_Axis.1), linewidth = lw) +
  geom_linerange(data = p_its_ma_centers, aes(x = mean_Axis.1, y = mean_Axis.2, ymin = ci_l_Axis.2, ymax = ci_u_Axis.2), linewidth = lw) +
  geom_point(data = p_its_ma_centers, aes(x = mean_Axis.1, y = mean_Axis.2, fill = field_type), size = lg_size, stroke = lw, shape = 21) +
  geom_point(aes(fill = field_type), size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  scale_fill_manual(values = ft_pal) +
  labs(
    x = paste0("PCoA 1 (", mva_its_ma$axis_pct[1], "%)"),
    y = paste0("PCoA 2 (", mva_its_ma$axis_pct[2], "%)")) +
  theme_ord +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
```

#### Unified figure

``` r
fig2_ls <- (its_rich_fig / plot_spacer() / plfa_fig) +
  plot_layout(heights = c(1,0.01,1)) 
fig2 <- (fig2_ls | plot_spacer() | its_ma_ord) +
  plot_layout(widths = c(0.35, 0.01, 0.64)) +
  plot_annotation(tag_levels = 'A') 
```

``` r
fig2
```

![](resources/fungal_ecology_files/figure-gfm/fig2-1.png)<!-- -->

**Fig 2.** Whole-soil fungal communities in **corn**, **restored**, and
**remnant** prairie fields. **a** OTU richness and **b** fungal biomass
(nmol PLFA g soil^-1) are shown as columns with 95 % CIs; lowercase
letters mark significant pairwise differences. **c**
Principal-coordinate (PCoA) ordination of ITS-based (97 % OTU) community
distances: small points = sites, large circles = field-type centroids
(error bars = 95 % CI). Cornfields cluster apart from restored or
remnant prairies. Numbers in circles give years since restoration. Axis
labels show the percent variation explained.

### Sequence-based relative abundance

Comparison figure and stats for supplemental

``` r
d_its <- its_avg %>% 
  data.frame(row.names = 1) %>%
  decostand("total") %>%
  vegdist("bray")
mva_its <- mva(d = d_its, env = sites)
```

``` r
mva_its$dispersion_test
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)  
    ## Groups     2 0.018698 0.0093489 3.2104   1999 0.0595 .
    ## Residuals 22 0.064065 0.0029121                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##              corn  remnant restored
    ## corn              0.123000   0.0710
    ## remnant  0.126039            0.1385
    ## restored 0.068726 0.135570

``` r
mva_its$permanova
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## adonis2(formula = d ~ dist_axis_1 + field_type, data = env, permutations = nperm, by = "terms")
    ##             Df SumOfSqs      R2      F Pr(>F)    
    ## dist_axis_1  1   0.4225 0.06253 1.7391 0.0270 *  
    ## field_type   2   1.2321 0.18236 2.5358 0.0005 ***
    ## Residual    21   5.1017 0.75510                  
    ## Total       24   6.7563 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
mva_its$pairwise_contrasts[c(1,3,2), c(1,2,4,3,8)] %>% 
  arrange(group1, desc(group2)) %>% 
  kable(format = "pandoc", caption = "Pairwise permanova contrasts")
```

| group1  | group2   | F_value |    R2 | p_value_adj |
|:--------|:---------|--------:|------:|------------:|
| corn    | restored |   3.913 | 0.164 |      0.0015 |
| corn    | remnant  |   2.858 | 0.281 |      0.0030 |
| remnant | restored |   1.062 | 0.054 |      0.3285 |

Pairwise permanova contrasts

No eignevalue correction was needed. Two relative eigenvalues exceeded
broken stick model. Based on the homogeneity of variance test, the null
hypothesis of equal variance among groups is accepted across all
clusters and in pairwise comparison of clusters (both p\>0.05),
supporting the application of a PERMANOVA test.

Clustering revealed that community variation was related to geographic
distance, the covariate in the model. With geographic distance accounted
for, the test variable ‘field type’ significantly explained variation in
fungal communities, with a post-hoc test revealing that communities in
corn fields differed from communities in restored and remnant fields.

Plotting results:

``` r
its_ord_data <- mva_its$ordination_scores %>% mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant")))
p_its_centers <- its_ord_data %>% 
  group_by(field_type) %>% 
  summarize(across(starts_with("Axis"), list(mean = mean, ci_l = ci_l, ci_u = ci_u), .names = "{.fn}_{.col}"), .groups = "drop") %>% 
  mutate(across(c(ci_l_Axis.1, ci_u_Axis.1), ~ mean_Axis.1 + .x),
         across(c(ci_l_Axis.2, ci_u_Axis.2), ~ mean_Axis.2 + .x))
its_ord <- 
  ggplot(its_ord_data, aes(x = Axis.1, y = Axis.2)) +
  geom_linerange(data = p_its_centers, aes(x = mean_Axis.1, y = mean_Axis.2, xmin = ci_l_Axis.1, xmax = ci_u_Axis.1), linewidth = lw) +
  geom_linerange(data = p_its_centers, aes(x = mean_Axis.1, y = mean_Axis.2, ymin = ci_l_Axis.2, ymax = ci_u_Axis.2), linewidth = lw) +
  geom_point(data = p_its_centers, aes(x = mean_Axis.1, y = mean_Axis.2, fill = field_type), size = lg_size, stroke = lw, shape = 21) +
  geom_point(aes(fill = field_type), size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  labs(
    x = paste0("PCoA 1 (", mva_its$axis_pct[1], "%)"),
    y = paste0("PCoA 2 (", mva_its$axis_pct[2], "%)")) +
  scale_fill_manual(values = ft_pal) +
  theme_ord +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
```

#### Supplemental figure

``` r
its_shan_ord_sup <- (its_shan_fig | plot_spacer() | its_ord) +
  plot_layout(widths = c(0.45, 0.01, 0.55)) +
  plot_annotation(tag_levels = 'A') 
```

``` r
its_shan_ord_sup
```

![](resources/fungal_ecology_files/figure-gfm/its_shan_ord_sup-1.png)<!-- -->

#### Contrast community metrics

Procrustes test on PCoA values using axes with eigenvalues exceeding a
broken stick model

``` r
set.seed(20251111)
its_protest <- protest(
  pcoa(d_its)$vectors[, 1:2],
  pcoa(d_its_ma)$vectors[, 1:2],
  permutations = 1999
)
its_protest
```

    ## 
    ## Call:
    ## protest(X = pcoa(d_its)$vectors[, 1:2], Y = pcoa(d_its_ma)$vectors[,      1:2], permutations = 1999) 
    ## 
    ## Procrustes Sum of Squares (m12 squared):        0.02873 
    ## Correlation in a symmetric Procrustes rotation: 0.9855 
    ## Significance:  5e-04 
    ## 
    ## Permutation: free
    ## Number of permutations: 1999

Including biomass changes little. The spatial configurations of both
ordinations are highly correlated. $R^{2}=$ 0.97, p\<0.001.

## Constrained Analysis

Test explanatory variables for correlation with site ordination. Using
plant data, so the analysis is restricted to Wisconsin sites. Edaphic
variables are too numerous to include individually, so transform micro
nutrients using PCA. Forb and grass cover is highly collinear; use the
grass-forb index produced previously with PCA.

``` r
soil_micro_pca <- 
  soil %>% 
  left_join(sites %>% select(field_name, field_type, region), by = join_by(field_name)) %>% 
  filter(field_type != "corn", region != "FL") %>% 
  select(field_name, SO4, Zn, Fe, Mn, Cu, Ca, Mg, Na, -field_key, -field_type, -region) %>% 
  column_to_rownames(var = "field_name") %>% 
  decostand(method = "standardize") %>% 
  rda()
summary(soil_micro_pca) # 70% on first two axes
```

    ## 
    ## Call:
    ## rda(X = .) 
    ## 
    ## Partitioning of variance:
    ##               Inertia Proportion
    ## Total               8          1
    ## Unconstrained       8          1
    ## 
    ## Eigenvalues, and their contribution to the variance 
    ## 
    ## Importance of components:
    ##                          PC1    PC2    PC3     PC4    PC5     PC6     PC7      PC8
    ## Eigenvalue            3.1372 1.8654 1.3683 0.65097 0.5312 0.32857 0.09349 0.024912
    ## Proportion Explained  0.3921 0.2332 0.1710 0.08137 0.0664 0.04107 0.01169 0.003114
    ## Cumulative Proportion 0.3921 0.6253 0.7964 0.87773 0.9441 0.98520 0.99689 1.000000

``` r
soil_micro_index <- scores(soil_micro_pca, choices = c(1, 2), display = "sites") %>% 
  data.frame() %>% 
  rename(soil_micro_1 = PC1, soil_micro_2 = PC2) %>% 
  rownames_to_column(var = "field_name")
soil_macro <- 
  soil %>% 
  left_join(sites %>% select(field_name, field_type, region), by = join_by(field_name)) %>% 
  filter(field_type != "corn", region != "FL") %>% 
  select(-c(field_key, field_type, region, SO4, Zn, Fe, Mn, Cu, Ca, Mg, Na))
```

Assemble explanatory variables and begin iterative selection process.
Plant functional groups and traits not included here were eliminated in
previous forward selection procedures (not shown). Check the VIF for
each explanatory variable to test for collinearity if model overfitting
is detected. Then run forward selection in `dbrda()`.

``` r
env_vars <- sites %>% 
  filter(field_type != "corn", region != "FL") %>% 
  select(field_name, dist_axis_1) %>% # 90% on axis 1
  left_join(soil_micro_index, by = join_by(field_name)) %>% # 70% on first two axes
  left_join(soil_macro, by = join_by(field_name)) %>% 
  left_join(gf_index, by = join_by(field_name)) %>% # 92% on axis 1
  left_join(prich %>% select(field_name, pl_rich), by = join_by(field_name)) %>% # plant richness
  select(-starts_with("field_key"), -soil_micro_1, -K) %>% # soil_micro_1, K removed based on initial VIF check
  column_to_rownames(var = "field_name") %>% 
  as.data.frame()
env_cov <- env_vars[,"dist_axis_1", drop = TRUE]
env_expl <- env_vars[, setdiff(colnames(env_vars), "dist_axis_1"), drop = FALSE] %>% 
  decostand("standardize")
```

Check VIF

``` r
env_expl %>% 
  cor() %>% 
  solve() %>% 
  diag() %>% 
  sort()
```

    ##          NO3      pl_rich soil_micro_2     gf_index            P           pH           OM 
    ##     2.250651     2.586388     2.645543     2.898170     2.903741     5.045472     5.463773

OM, K, and soil_micro_1 with high VIF in initial VIF check. Removed
soil_micro_1 and K to maintain OM in the model. No overfitting detected
in full model; proceed with forward selection.

``` r
spe_its_wi_resto <- its_avg_ma %>% 
  filter(field_name %in% rownames(env_expl)) %>% 
  column_to_rownames(var = "field_name")
mod_null <- dbrda(spe_its_wi_resto ~ 1 + Condition(env_cov), data = env_expl, distance = "bray")
mod_full <- dbrda(spe_its_wi_resto ~ . + Condition(env_cov), data = env_expl, distance = "bray")
mod_step <- ordistep(mod_null, 
                     scope = formula(mod_full), 
                     direction = "forward", 
                     permutations = 1999, 
                     trace = FALSE)
```

### Constrained Analysis Results

``` r
mod_step
```

    ## 
    ## Call: dbrda(formula = spe_its_wi_resto ~ Condition(env_cov) + gf_index + pH + pl_rich, data = env_expl, distance = "bray")
    ## 
    ##               Inertia Proportion Rank
    ## Total          3.3581     1.0000     
    ## Conditional    0.3207     0.0955    1
    ## Constrained    1.3561     0.4038    3
    ## Unconstrained  1.6813     0.5007    8
    ## 
    ## Inertia is squared Bray distance
    ## 
    ## Eigenvalues for constrained axes:
    ## dbRDA1 dbRDA2 dbRDA3 
    ## 0.6286 0.3979 0.3296 
    ## 
    ## Eigenvalues for unconstrained axes:
    ##   MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8 
    ## 0.3433 0.3313 0.2409 0.1994 0.1838 0.1568 0.1267 0.0990

``` r
(mod_r2   <- RsquareAdj(mod_step, permutations = 1999))
```

    ## $r.squared
    ## [1] 0.403825
    ## 
    ## $adj.r.squared
    ## [1] 0.2357145

``` r
(mod_glax <- anova(mod_step, permutations = 1999))
```

    ## Permutation test for dbrda under reduced model
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## Model: dbrda(formula = spe_its_wi_resto ~ Condition(env_cov) + gf_index + pH + pl_rich, data = env_expl, distance = "bray")
    ##          Df SumOfSqs      F Pr(>F)    
    ## Model     3   1.3561 2.1508  5e-04 ***
    ## Residual  8   1.6813                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(mod_inax <- anova(mod_step, by = "axis", permutations = 1999))
```

    ## Permutation test for dbrda under reduced model
    ## Forward tests for axes
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## Model: dbrda(formula = spe_its_wi_resto ~ Condition(env_cov) + gf_index + pH + pl_rich, data = env_expl, distance = "bray")
    ##          Df SumOfSqs      F Pr(>F)    
    ## dbRDA1    1  0.62863 2.9911 0.0010 ***
    ## dbRDA2    1  0.39788 2.1299 0.0035 ** 
    ## dbRDA3    1  0.32956 1.9602 0.0035 ** 
    ## Residual  8  1.68130                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(mod_axpct <- round(100 * mod_step$CCA$eig / sum(mod_step$CCA$eig), 1))
```

    ## dbRDA1 dbRDA2 dbRDA3 
    ##   46.4   29.3   24.3

``` r
anova(mod_step, by = "margin", permutations = 1999) %>% 
  as.data.frame() %>% 
  mutate(p.adj = p.adjust(`Pr(>F)`, "fdr")) %>% 
  kable(, format = "pandoc")
```

|          |  Df |  SumOfSqs |        F | Pr(\>F) |   p.adj |
|----------|----:|----------:|---------:|--------:|--------:|
| gf_index |   1 | 0.4958334 | 2.359281 |  0.0025 | 0.00750 |
| pH       |   1 | 0.4435865 | 2.110679 |  0.0055 | 0.00825 |
| pl_rich  |   1 | 0.3673145 | 1.747760 |  0.0215 | 0.02150 |
| Residual |   8 | 1.6813037 |       NA |      NA |      NA |

Based on permutation tests with n=1999 permutations, the model shows a
significant correlation between the site ordination on fungal
communities and the selected explanatory variables (p\<0.001). The first
two constrained axes are also significant (p\<0.001, P\<0.01). The
selected variables explain $R^{2}_{\text{Adj}}$=21.3% of the community
variation. Selected explanatory variables are pH and the grass-forb
index; see table for individual p values and statistics.

Create the figure objects. Figure 6 will be produced with panels from
other groups, see code in the saprotrophs section.

``` r
mod_step_eig <- round(mod_step$CCA$eig * 100, 1)
mod_scor <- scores(
  mod_step,
  choices = c(1, 2),
  display = c("bp", "sites"),
  tidy = FALSE
)
mod_scor_site <- mod_scor$sites %>% 
  data.frame() %>%
  rownames_to_column(var = "field_name") %>% 
  left_join(sites, by = join_by(field_name))
mod_scor_bp <- bind_rows(
  mod_scor$biplot %>% 
    data.frame() %>% 
    rownames_to_column(var = "envvar") %>% 
    mutate(envlabs = c(">forb", "pH", "plant spp.")),
  data.frame(
    envvar = "gf_index",
    dbRDA1 = -mod_scor$biplot["gf_index", 1],
    dbRDA2 = -mod_scor$biplot["gf_index", 2],
    envlabs = ">grass")
) %>% 
  arrange(envvar, envlabs) %>% 
  mutate(
    origin = 0,
    m = dbRDA2 / dbRDA1, 
    d = sqrt(dbRDA1^2 + dbRDA2^2), 
    dadd = sqrt((max(dbRDA1)-min(dbRDA2))^2 + (max(dbRDA2)-min(dbRDA2))^2)*dadd_adj,
    labx = ((d+dadd)*cos(atan(m)))*(dbRDA1/abs(dbRDA1)), 
    laby = ((d+dadd)*sin(atan(m)))*(dbRDA1/abs(dbRDA1)))
```

## Soil fungi and plant community correlations

Data for these tests

``` r
fungi_resto <- its_div %>% 
  left_join(fa %>% select(field_name, fungi_mass = fungi_18.2), by = join_by(field_name)) %>% 
  left_join(sites, by = join_by(field_name, field_type)) %>% 
  left_join(gf_index, by = join_by(field_name)) %>% 
  left_join(prich %>% select(field_name, pl_rich, pl_shan), by = join_by(field_name)) %>% 
  filter(field_type != "corn", region != "FL") %>% 
  select(field_name, fungi_ab = depth, fungi_mass, gf_index, pl_rich, pl_shan)
```

### Plant alpha diversity and fungal biomass

Is plant richness related to pathogen mass?

``` r
fa_prich_lm <- lm(fungi_mass ~ pl_rich, data = fungi_resto)
summary(fa_prich_lm)
```

    ## 
    ## Call:
    ## lm(formula = fungi_mass ~ pl_rich, data = fungi_resto)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -2.4809 -0.6583  0.2743  0.7575  1.9430 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  7.84735    1.38756   5.656 0.000148 ***
    ## pl_rich     -0.07289    0.03368  -2.164 0.053335 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.436 on 11 degrees of freedom
    ## Multiple R-squared:  0.2986, Adjusted R-squared:  0.2348 
    ## F-statistic: 4.683 on 1 and 11 DF,  p-value: 0.05334

Fungal mass and plant richness are weakly correlated but driven by a
high-leverage point (not shown). When seq abund is a response and
log(mass) included as a covariate, no relationship is detected (not
shown).

Is plant diversity related to fungal mass?

``` r
fa_pshan_lm <- lm(fungi_mass ~ pl_shan, data = fungi_resto)
summary(fa_pshan_lm)
```

    ## 
    ## Call:
    ## lm(formula = fungi_mass ~ pl_shan, data = fungi_resto)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -2.0899 -1.4918  0.4519  0.7842  2.6614 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)   7.9113     1.8730   4.224  0.00143 **
    ## pl_shan      -0.2222     0.1378  -1.612  0.13517   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.542 on 11 degrees of freedom
    ## Multiple R-squared:  0.1912, Adjusted R-squared:  0.1176 
    ## F-statistic:   2.6 on 1 and 11 DF,  p-value: 0.1352

Fungal biomass and plant diversity are negatively related but the
correlation is not significant. It’s driven almost entirely by KORP (not
shown) and wouldn’t be close to significant otherwise, no further
testing warranted.

### Fungal biomass and grass/forb composition

Inspect simple linear relationship. Naïve model.

``` r
fuma_rest_m <- lm(fungi_mass ~ gf_index, data = fungi_resto)
```

``` r
check_model(fuma_rest_m)
```

![](resources/fungal_ecology_files/figure-gfm/cm1-1.png)<!-- -->

``` r
summary(fuma_rest_m)
```

    ## 
    ## Call:
    ## lm(formula = fungi_mass ~ gf_index, data = fungi_resto)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -2.7014 -1.6462  0.6455  1.0298  2.2321 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    4.971      0.442  11.246 2.26e-07 ***
    ## gf_index      -2.180      1.657  -1.315    0.215    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.594 on 11 degrees of freedom
    ## Multiple R-squared:  0.1359, Adjusted R-squared:  0.05736 
    ## F-statistic:  1.73 on 1 and 11 DF,  p-value: 0.2151

The relationship is poor and needs no further analysis

# Arbuscular mycorrhizal fungi

``` r
# AMF ———————— ####
```

## Diversity Indices

``` r
amf_div <- calc_div(amf_avg, sites) %>% 
  mutate(depth_csq = sqrt(depth) - mean(sqrt(depth)))
```

### Richness

Sequence depth square root transformed and centered. Negative binomial
model was underdispersed and failed to converge at default iterations;
use poisson glm instead.

Test interaction

``` r
amf_rich_glm_i <- glm(richness ~ depth_csq * field_type, data = amf_div, family = poisson(link = "log")) 
Anova(amf_rich_glm_i, type = 3, test.statistic = "LR") # interaction detected
```

    ## Analysis of Deviance Table (Type III tests)
    ## 
    ## Response: richness
    ##                      LR Chisq Df Pr(>Chisq)   
    ## depth_csq              4.4899  1   0.034096 * 
    ## field_type            10.4165  2   0.005471 **
    ## depth_csq:field_type   6.1037  2   0.047271 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
check_overdispersion(amf_rich_glm_i) # not overdispersed
```

    ## # Overdispersion test
    ## 
    ##        dispersion ratio =  0.659
    ##   Pearson's Chi-Squared = 12.519
    ##                 p-value =  0.862

    ## No overdispersion detected.

``` r
augment(amf_rich_glm_i) # corn site has cooks >0.9
```

    ## # A tibble: 25 × 9
    ##    richness depth_csq field_type .fitted   .resid   .hat .sigma    .cooksd .std.resid
    ##       <int>     <dbl> <fct>        <dbl>    <dbl>  <dbl>  <dbl>      <dbl>      <dbl>
    ##  1       59   -10.3   restored      4.00  0.573   0.232   0.824 0.0221        0.654  
    ##  2       47     9.30  restored      3.94 -0.624   0.171   0.823 0.0157       -0.685  
    ##  3       38    -0.862 corn          3.69 -0.357   0.214   0.833 0.00722      -0.403  
    ##  4       52     3.80  remnant       4.03 -0.577   0.366   0.821 0.0491       -0.724  
    ##  5       53     8.24  restored      3.94  0.197   0.147   0.837 0.00131       0.213  
    ##  6       60     3.69  corn          3.91  1.39    0.547   0.683 0.914         2.06   
    ##  7       33    -4.72  corn          3.51 -0.0975  0.560   0.838 0.00457      -0.147  
    ##  8       62     6.11  remnant       4.06  0.560   0.500   0.817 0.107         0.792  
    ##  9       53     1.52  restored      3.96  0.0425  0.0642  0.838 0.0000221     0.0439 
    ## 10       54    -5.89  restored      3.99  0.00740 0.120   0.838 0.00000142    0.00789
    ## # ℹ 15 more rows

``` r
check_collinearity(amf_rich_glm_i) # depth and field_type VIF > 26
```

    ## Model has interaction terms. VIFs might be inflated.
    ##   Try to center the variables used for the interaction, or check multicollinearity among predictors of a model without interaction terms.

    ## # Check for Multicollinearity
    ## 
    ## Low Correlation
    ## 
    ##        Term  VIF        VIF 95% CI adj. VIF Tolerance Tolerance 95% CI
    ##  field_type 1.02 [ 1.00, 2.55e+05]     1.00      0.98     [0.00, 1.00]
    ## 
    ## High Correlation
    ## 
    ##                  Term   VIF        VIF 95% CI adj. VIF Tolerance Tolerance 95% CI
    ##             depth_csq 26.03 [15.52,    44.15]     5.10      0.04     [0.02, 0.06]
    ##  depth_csq:field_type 26.25 [15.65,    44.51]     5.12      0.04     [0.02, 0.06]

An interaction was detected, but including it in the model leads to very
poor diagnostics. It’s driven by one site in corn with high leverage,
and it introduces high multicollinearity. Further, the outlier point
would tend to lead to a Type II error of inference, making it a
conservative choice to stick with the additive model.

Fit additive model

``` r
amf_rich_glm <- glm(richness ~ depth_csq + field_type, data = amf_div, family = poisson(link = "log")) 
```

Diagnostics

``` r
check_model(amf_rich_glm)
```

![](resources/fungal_ecology_files/figure-gfm/amf_rich_covar_diagnostics-1.png)<!-- -->

``` r
check_overdispersion(amf_rich_glm)
```

    ## # Overdispersion test
    ## 
    ##        dispersion ratio =  0.927
    ##   Pearson's Chi-Squared = 19.467
    ##                 p-value =  0.555

    ## No overdispersion detected.

``` r
check_collinearity(amf_rich_glm)
```

    ## # Check for Multicollinearity
    ## 
    ## Low Correlation
    ## 
    ##        Term  VIF       VIF 95% CI adj. VIF Tolerance Tolerance 95% CI
    ##   depth_csq 1.01 [1.00, 1.26e+07]     1.01      0.99     [0.00, 1.00]
    ##  field_type 1.01 [1.00, 1.26e+07]     1.00      0.99     [0.00, 1.00]

Long tails, some midrange structure, no leverage points, overdispersion,
or multicollinearity

``` r
distribution_prob(amf_rich_glm)
```

    ## 
    ## 
    ## Distribution    p_Residuals
    ## -------------  ------------
    ## normal              0.31250
    ## cauchy              0.28125
    ## exponential         0.15625
    ## 
    ## 
    ## Distribution                  p_Response
    ## ---------------------------  -----------
    ## beta-binomial                    0.46875
    ## neg. binomial (zero-infl.)       0.25000
    ## normal                           0.09375

residuals distribution normal or long-tailed, response count-distributed

``` r
leveneTest(richness ~ field_type, data = amf_div) %>% as.data.frame() %>% kable(format = "pandoc")
```

|       |  Df |   F value |   Pr(\>F) |
|-------|----:|----------:|----------:|
| group |   2 | 0.7049808 | 0.5049423 |
|       |  22 |        NA |        NA |

``` r
leveneTest(residuals(amf_rich_glm) ~ amf_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc")
```

|       |  Df |  F value |   Pr(\>F) |
|-------|----:|---------:|----------:|
| group |   2 | 1.001082 | 0.3836147 |
|       |  22 |       NA |        NA |

Residuals/response distributions do not suggest the need for
transformation. Levene’s p \> 0.05 → fail to reject = variances can be
considered equal across groups.

Model results, group means, and post-hoc. Use Type II LR test of
variables due to unbalanced design.

``` r
Anova(amf_rich_glm, type = 2, test.statistic = "LR")
```

    ## Analysis of Deviance Table (Type II tests)
    ## 
    ## Response: richness
    ##            LR Chisq Df Pr(>Chisq)   
    ## depth_csq    0.2733  1   0.601100   
    ## field_type  10.2189  2   0.006039 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Sequencing depth not a significant predictor of amf richness

``` r
amf_rich_em <- emmeans(amf_rich_glm, ~ field_type, type = "response")
```

Results tables below show the emmeans summary of estimated marginal
means and confidence intervals, and the post hoc contrast of richness
among field types. Main effect in model significant; pairwise contrast
warranted.

| field_type |     rate |       SE |  df | asymp.LCL | asymp.UCL |
|:-----------|---------:|---------:|----:|----------:|----------:|
| corn       | 41.82633 | 2.893593 | Inf |  36.52268 |  47.90015 |
| restored   | 52.80884 | 1.820856 | Inf |  49.35794 |  56.50100 |
| remnant    | 53.44030 | 3.678202 | Inf |  46.69627 |  61.15832 |

Confidence level used: 0.95

| contrast           |     ratio |        SE |  df | null |   z.ratio |   p.value |
|:-------------------|----------:|----------:|----:|-----:|----------:|----------:|
| corn / restored    | 0.7920328 | 0.0612532 | Inf |    1 | -3.014772 | 0.0072555 |
| corn / remnant     | 0.7826740 | 0.0763188 | Inf |    1 | -2.512954 | 0.0320995 |
| restored / remnant | 0.9881838 | 0.0762908 | Inf |    1 | -0.153965 | 0.9870163 |

P value adjustment: tukey method for comparing a family of 3 estimates

OTU richness in cornfields is significantly less than in restored or
remnant fields, which don’t differ. Plot the results

``` r
amf_rich_fig <- 
  ggplot(summary(amf_rich_em), aes(x = field_type, y = rate)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = rate, ymax = asymp.UCL), width = 0, linewidth = lw) +
  geom_text(aes(y = asymp.UCL, label = c("a", "b", "b")),  vjust = -1, family = "sans", size = 3.5) +
  labs(x = NULL, y = expression(atop("Richness", paste("(", italic(n), " OTUs)")))) +
  lims(y = c(0, 75)) +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1))
```

### Shannon diversity

Sequence depth square root transformed and centered

``` r
amf_shan_lm <- lm(shannon ~ depth_csq + field_type, data = amf_div)
```

Diagnostics

``` r
check_model(amf_shan_lm)
```

![](resources/fungal_ecology_files/figure-gfm/amf_shan_covar_diagnostics-1.png)<!-- -->

Variance appears somewhat non-constant in groups, qqplot fit is poor,
one leverage point (Cook’s \> 0.5), a cornfield with high richness. Mean
richness in corn fields is lowest; this outlier would make the pairwise
contrast less significant, possible Type II error which is more
acceptable.

``` r
distribution_prob(amf_shan_lm)
```

    ## 
    ## 
    ## Distribution    p_Residuals
    ## -------------  ------------
    ## normal              0.68750
    ## cauchy              0.18750
    ## chi                 0.03125
    ## 
    ## 
    ## Distribution    p_Response
    ## -------------  -----------
    ## normal             0.25000
    ## uniform            0.15625
    ## pareto             0.12500

Residuals/response distributions most likely normal.

``` r
leveneTest(shannon ~ field_type, data = amf_div) %>% as.data.frame() %>% kable(format = "pandoc")
```

|       |  Df |   F value |   Pr(\>F) |
|-------|----:|----------:|----------:|
| group |   2 | 0.2896534 | 0.7513336 |
|       |  22 |        NA |        NA |

``` r
leveneTest(residuals(amf_shan_lm) ~ amf_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc")
```

|       |  Df |   F value |   Pr(\>F) |
|-------|----:|----------:|----------:|
| group |   2 | 0.3022751 | 0.7421555 |
|       |  22 |        NA |        NA |

Residuals/response distributions do not suggest the need for
transformation. Levene’s p \> 0.05 → fail to reject = variances can be
considered equal.

Model results, group means, and post-hoc

``` r
Anova(amf_shan_lm, type = 2)
```

    ## Anova Table (Type II tests)
    ## 
    ## Response: shannon
    ##             Sum Sq Df F value    Pr(>F)    
    ## depth_csq    0.123  1  0.0127 0.9114841    
    ## field_type 258.303  2 13.2691 0.0001881 ***
    ## Residuals  204.398 21                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Sequencing depth not a significant predictor of Shannon diversity.
Produce arithmetic means in groups and post hoc contrasts

``` r
amf_shan_em <- emmeans(amf_shan_lm, ~ field_type, type = "response")
```

Results tables below show the emmeans summary of group means and
confidence intervals, with sequencing depth as a covariate, and the post
hoc contrast of richness among field types.

| field_type |   emmean |        SE |  df | lower.CL | upper.CL |
|:-----------|---------:|----------:|----:|---------:|---------:|
| corn       | 14.72048 | 1.3955076 |  21 | 11.81837 | 17.62260 |
| restored   | 21.52479 | 0.7814914 |  21 | 19.89959 | 23.14999 |
| remnant    | 24.84687 | 1.5681686 |  21 | 21.58568 | 28.10805 |

Confidence level used: 0.95

| contrast           |   estimate |       SE |  df |   t.ratio |   p.value |
|:-------------------|-----------:|---------:|----:|----------:|----------:|
| corn - restored    |  -6.804309 | 1.600293 |  21 | -4.251915 | 0.0009952 |
| corn - remnant     | -10.126384 | 2.097026 |  21 | -4.828927 | 0.0002549 |
| restored - remnant |  -3.322075 | 1.756599 |  21 | -1.891197 | 0.1659595 |

P value adjustment: tukey method for comparing a family of 3 estimates

Shannon’s diversity in cornfields is significantly less than in restored
or remnant fields, which don’t differ.

``` r
amf_shan_fig <- 
  ggplot(summary(amf_shan_em), aes(x = field_type, y = emmean)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = emmean, ymax = upper.CL), width = 0, linewidth = lw) +
  geom_text(aes(y = upper.CL, label = c("a", "b", "b")),  vjust = -1, family = "sans", size = 3.5) +
  labs(x = "Field type", y = expression(Shannon~diversity~paste("(", italic(e)^italic(H), ")"))) +
  lims(y = c(0, 32)) +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
```

## Abundance (NLFA biomass)

``` r
nlfa_lm <- lm(amf ~ field_type, data = fa)
```

Diagnostics

``` r
par(mfrow = c(2,2))
plot(nlfa_lm) # variance obviously not constant in groups
```

![](resources/fungal_ecology_files/figure-gfm/unnamed-chunk-65-1.png)<!-- -->

``` r
distribution_prob(nlfa_lm)
```

    ## 
    ## 
    ## Distribution    p_Residuals
    ## -------------  ------------
    ## cauchy              0.62500
    ## normal              0.15625
    ## chi                 0.06250
    ## 
    ## 
    ## Distribution    p_Response
    ## -------------  -----------
    ## gamma              0.40625
    ## chi                0.12500
    ## half-cauchy        0.12500

``` r
# response distribution gamma; resids likely normal
leveneTest(residuals(nlfa_lm) ~ fa$field_type) # No covariate, response and residuals tests equivalent
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  2  3.3372 0.05423 .
    ##       22                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Residuals distribution variance may not be equal in groups. Levene’s p =
0.054, close to rejecting the null of equal variance. Check CV in
groups.

``` r
fa %>%
  mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant"))) %>%
  group_by(field_type) %>%
  summarize(mean = mean(amf),
            cv = sd(amf) / mean) %>%
  mutate(across(mean:cv, ~ round(.x, 2))) %>%
  kable(format = "pandoc", caption = "Mean and CV relationship in groups")
```

| field_type |  mean |   cv |
|:-----------|------:|-----:|
| corn       |  3.79 | 0.26 |
| restored   | 32.27 | 0.54 |
| remnant    | 34.82 | 0.59 |

Mean and CV relationship in groups

CV increases with mean, suggesting \> proportional mean/variance
relationship. Determine best model choice of log-transformed response or
gamma glm. Log:

``` r
nlfa_lm_log <- lm(log(amf) ~ field_type, data = fa)
par(mfrow = c(2,2))
plot(nlfa_lm_log) # qqplot ok, one high leverage point in remnants
```

![](resources/fungal_ecology_files/figure-gfm/unnamed-chunk-67-1.png)<!-- -->

``` r
ncvTest(nlfa_lm_log) # p=0.16, null of constant variance not rejected
```

    ## Non-constant Variance Score Test 
    ## Variance formula: ~ fitted.values 
    ## Chisquare = 1.989836, Df = 1, p = 0.15836

Gamma glm:

``` r
nlfa_glm  <- glm(amf ~ field_type, family = Gamma(link = "log"), data = fa)
nlfa_glm_diag <- glm.diag(nlfa_glm)
glm.diag.plots(nlfa_glm, nlfa_glm_diag) # qqplot shows strong fit; no leverage >0.5
```

![](resources/fungal_ecology_files/figure-gfm/unnamed-chunk-68-1.png)<!-- -->

``` r
performance::check_overdispersion(nlfa_glm) # not detected
```

    ## # Overdispersion test
    ## 
    ##  dispersion ratio = 1.259
    ##           p-value =  0.48

    ## No overdispersion detected.

Gamma glm is the best choice; no high-leverage point

Model results, group means, and post-hoc

``` r
Anova(nlfa_glm, test.statistic = "LR") 
```

    ## Analysis of Deviance Table (Type II tests)
    ## 
    ## Response: amf
    ##            LR Chisq Df Pr(>Chisq)    
    ## field_type   46.369  2  8.533e-11 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
nlfa_em <- emmeans(nlfa_glm, ~ field_type, type = "response")
```

| field_type |  response |        SE |  df |  lower.CL | upper.CL |
|:-----------|----------:|----------:|----:|----------:|---------:|
| corn       |  3.789798 | 0.8571948 |  22 |  2.370816 |  6.05807 |
| restored   | 32.269772 | 4.0802271 |  22 | 24.826411 | 41.94477 |
| remnant    | 34.817071 | 8.8046212 |  22 | 20.607643 | 58.82422 |

Confidence level used: 0.95. Intervals are back-transformed from the log
scale

| contrast           |     ratio |        SE |  df | null |   t.ratio |   p.value |
|:-------------------|----------:|----------:|----:|-----:|----------:|----------:|
| corn / restored    | 0.1174411 | 0.0304322 |  22 |    1 | -8.265507 | 0.0000001 |
| corn / remnant     | 0.1088489 | 0.0369299 |  22 |    1 | -6.536823 | 0.0000041 |
| restored / remnant | 0.9268376 | 0.2620457 |  22 |    1 | -0.268725 | 0.9610443 |

P value adjustment: tukey method for comparing a family of 3 estimates.
Tests are performed on the log scale

``` r
nlfa_fig <-
  ggplot(summary(nlfa_em), aes(x = field_type, y = response)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = response, ymax = upper.CL), width = 0, linewidth = lw) +
  geom_text(aes(y = upper.CL, label = c("a", "b", "b")),  vjust = -1, family = "sans", size = 3.5) +
  labs(x = "Field type", y = expression(atop("Biomass", paste("(", nmol[NLFA], " × ", g[soil]^{-1}, ")")))) +
  # labs(x = "Field Type", y = "Biomass") +
  scale_fill_manual(values = ft_pal) +
  lims(y = c(0, 75)) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1.1))
```

## Beta Diversity

1.  Using biomass-weighted relative abundance and b-c distance (unifrac
    is scale invariant; it’s based on the proportion of reads on each
    descending branch, multiplying rows by any constant doesn’t change
    this).
2.  Using sequence-based relative abundance and unifrac distance to
    display phylogenetically aware information.
3.  Contrast the two with procrustes.

### Biomass-weighted relative abundance

``` r
d_amf_ma <- amf_avg_ma %>% 
  data.frame(row.names = 1) %>% 
  vegdist("bray")
mva_amf_ma <- mva(d = d_amf_ma, env = sites, corr = "lingoes")
```

``` r
mva_amf_ma$dispersion_test
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)
    ## Groups     2 0.020317 0.010158 0.9907   1999 0.3895
    ## Residuals 22 0.225579 0.010254                     
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##             corn remnant restored
    ## corn             0.41250   0.1685
    ## remnant  0.38741           0.6570
    ## restored 0.16747 0.65170

``` r
mva_amf_ma$permanova
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## adonis2(formula = d ~ dist_axis_1 + field_type, data = env, permutations = nperm, by = "terms")
    ##             Df SumOfSqs      R2      F Pr(>F)    
    ## dist_axis_1  1   0.2327 0.03986 1.3035 0.2215    
    ## field_type   2   1.8558 0.31794 5.1982 0.0005 ***
    ## Residual    21   3.7485 0.64220                  
    ## Total       24   5.8369 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
mva_amf_ma$pairwise_contrasts[c(1,3,2), c(1,2,4,3,8)] %>% 
  arrange(group1, desc(group2)) %>% 
  kable(format = "pandoc", caption = "Pairwise permanova contrasts")
```

| group1  | group2   | F_value |    R2 | p_value_adj |
|:--------|:---------|--------:|------:|------------:|
| corn    | restored |   9.902 | 0.334 |      0.0015 |
| corn    | remnant  |   6.251 | 0.464 |      0.0015 |
| remnant | restored |   0.430 | 0.023 |      0.9405 |

Pairwise permanova contrasts

Lingoes correction was applied to negative eignevalues. Three relative
eigenvalues exceeded broken stick model. Based on the homogeneity of
variance test, the null hypothesis of equal variance among groups is
accepted across all clusters and in pairwise comparison of clusters
(both p\>0.05), supporting the application of a PERMANOVA test.

Clustering revealed that community variation was not related to
geographic distance, the covariate in the model. With geographic
distance accounted for, the test variable ‘field type’ significantly
explained variation in fungal communities, with a post-hoc test
revealing that communities in corn fields differed from communities in
restored and remnant fields.

Plotting results:

``` r
amf_ma_ord_data <- mva_amf_ma$ordination_scores %>% mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant")))
p_amf_ma_centers <- amf_ma_ord_data %>% 
  group_by(field_type) %>% 
  summarize(across(starts_with("Axis"), list(mean = mean, ci_l = ci_l, ci_u = ci_u), .names = "{.fn}_{.col}"), .groups = "drop") %>% 
  mutate(across(c(ci_l_Axis.1, ci_u_Axis.1), ~ mean_Axis.1 + .x),
         across(c(ci_l_Axis.2, ci_u_Axis.2), ~ mean_Axis.2 + .x),
         across(ends_with("Axis.1"), ~ .x))
amf_ma_ord <- 
  ggplot(amf_ma_ord_data, aes(x = Axis.1, y = Axis.2)) + 
  geom_linerange(data = p_amf_ma_centers, aes(x = mean_Axis.1, y = mean_Axis.2, xmin = ci_l_Axis.1, xmax = ci_u_Axis.1), linewidth = lw) +
  geom_linerange(data = p_amf_ma_centers, aes(x = mean_Axis.1, y = mean_Axis.2, ymin = ci_l_Axis.2, ymax = ci_u_Axis.2), linewidth = lw) +
  geom_point(data = p_amf_ma_centers, aes(x = mean_Axis.1, y = mean_Axis.2, fill = field_type), size = lg_size, stroke = lw, shape = 21) +
  geom_point(aes(fill = field_type), size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  scale_fill_manual(values = ft_pal) +
  labs(
    x = paste0("PCoA 1 (", mva_amf_ma$axis_pct[1], "%)"),
    y = paste0("PCoA 2 (", mva_amf_ma$axis_pct[2], "%)")) +
  theme_ord +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
```

#### Unified figure

``` r
fig3_ls <- (amf_rich_fig / plot_spacer() / nlfa_fig) +
  plot_layout(heights = c(1,0.01,1)) 
fig3 <- (fig3_ls | plot_spacer() | amf_ma_ord) +
  plot_layout(widths = c(0.35, 0.01, 0.64)) +
  plot_annotation(tag_levels = 'A') 
```

``` r
fig3
```

![](resources/fungal_ecology_files/figure-gfm/fig3-1.png)<!-- -->

**Fig 3.** AMF communities in corn, restored, and remnant prairie
fields. OTU richness **a**; biomass **b** (nmol NLFA g soil^-1), 95 %
CI, letters = Tukey groups). PCoA of BC distances on proportion of
biomass abundance data (18S, 97 % OTUs) **c**: small points = sites,
large rings = field‑type centroids ±95 % CI. Numbers in points give
years since restoration. Axes show % variance. Corn clusters apart from
both prairie types.

### Sequence-based relative abundance, unifrac distance

``` r
d_amf <- UniFrac(amf_ps, weighted = TRUE, normalized = TRUE)
mva_amf <- mva(d = d_amf, env = sites, corr = "lingoes")
```

``` r
mva_amf$dispersion_test
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
    ## Groups     2 0.000418 0.0002089 0.0647   1999 0.9285
    ## Residuals 22 0.071014 0.0032279                     
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##             corn remnant restored
    ## corn             0.89350    0.853
    ## remnant  0.89942            0.702
    ## restored 0.85873 0.71820

``` r
mva_amf$permanova
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## adonis2(formula = d ~ dist_axis_1 + field_type, data = env, permutations = nperm, by = "terms")
    ##             Df SumOfSqs      R2      F Pr(>F)    
    ## dist_axis_1  1  0.04776 0.05566 1.6777 0.1135    
    ## field_type   2  0.21243 0.24757 3.7307 0.0005 ***
    ## Residual    21  0.59788 0.69677                  
    ## Total       24  0.85808 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
mva_amf$pairwise_contrasts[c(1,3,2), c(1,2,4,3,8)] %>% 
  arrange(group1, desc(group2)) %>% 
  kable(format = "pandoc", caption = "Pairwise permanova contrasts")
```

| group1  | group2   | F_value |    R2 | p_value_adj |
|:--------|:---------|--------:|------:|------------:|
| corn    | restored |   6.478 | 0.250 |      0.0008 |
| corn    | remnant  |   4.655 | 0.355 |      0.0008 |
| remnant | restored |   0.442 | 0.023 |      0.8705 |

Pairwise permanova contrasts

Lingoes eigenvalue correction was used. The first three relative
eigenvalues exceeded broken stick model. Based on the homogeneity of
variance test, the null hypothesis of equal variance among groups is
accepted across all clusters and in pairwise comparison of clusters
(both p\>0.05), supporting the application of a PERMANOVA test.

Clustering revealed that geographic distance among sites did not
significantly explain AMF community variation. With geographic distance
accounted for, the test variable field type significantly explained
variation in AMF communities, with a post-hoc test revealing that
communities in corn fields differed from communities in restored and
remnant fields.

Plotting the result:

``` r
amf_ord_data <- mva_amf$ordination_scores %>% mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant")))
p_amf_centers <- amf_ord_data %>% 
  group_by(field_type) %>% 
  summarize(across(starts_with("Axis"), list(mean = mean, ci_l = ci_l, ci_u = ci_u), .names = "{.fn}_{.col}"), .groups = "drop") %>% 
  mutate(across(c(ci_l_Axis.1, ci_u_Axis.1), ~ mean_Axis.1 + .x),
         across(c(ci_l_Axis.2, ci_u_Axis.2), ~ mean_Axis.2 + .x),
         across(ends_with("Axis.1"), ~ .x * -1)) # reversed for consistency
amf_ord <- 
  ggplot(amf_ord_data, aes(x = Axis.1 * -1, y = Axis.2)) + # reversed for consistency
  geom_linerange(data = p_amf_centers, aes(x = mean_Axis.1, y = mean_Axis.2, xmin = ci_l_Axis.1, xmax = ci_u_Axis.1), linewidth = lw) +
  geom_linerange(data = p_amf_centers, aes(x = mean_Axis.1, y = mean_Axis.2, ymin = ci_l_Axis.2, ymax = ci_u_Axis.2), linewidth = lw) +
  geom_point(data = p_amf_centers, aes(x = mean_Axis.1, y = mean_Axis.2, fill = field_type), size = lg_size, stroke = lw, shape = 21) +
  geom_point(aes(fill = field_type), size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  scale_fill_manual(values = ft_pal) +
  labs(
    x = paste0("PCoA 1 (", mva_amf$axis_pct[1], "%)"),
    y = paste0("PCoA 2 (", mva_amf$axis_pct[2], "%)")) +
  theme_ord +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
```

#### Supplemental figure

``` r
amf_shan_ord_sup <- (amf_shan_fig | plot_spacer() | amf_ord) +
  plot_layout(widths = c(0.45, 0.01, 0.55)) +
  plot_annotation(tag_levels = 'A') 
```

``` r
amf_shan_ord_sup
```

![](resources/fungal_ecology_files/figure-gfm/amf_shan_ord_sup-1.png)<!-- -->

#### Contrast community metrics

Procrustes test on PCoA values using axes with eigenvalues exceeding a
broken stick model

``` r
set.seed(20251111)
amf_protest <- protest(
  pcoa(d_amf, correction = "lingoes")$vectors[, 1:3],
  pcoa(d_amf_ma, correction = "lingoes")$vectors[, 1:3],
  permutations = 1999
)
amf_protest
```

    ## 
    ## Call:
    ## protest(X = pcoa(d_amf, correction = "lingoes")$vectors[, 1:3],      Y = pcoa(d_amf_ma, correction = "lingoes")$vectors[, 1:3],      permutations = 1999) 
    ## 
    ## Procrustes Sum of Squares (m12 squared):        0.3712 
    ## Correlation in a symmetric Procrustes rotation: 0.793 
    ## Significance:  5e-04 
    ## 
    ## Permutation: free
    ## Number of permutations: 1999

The null that these solutions are unrelated is rejected at p\<0.001.
However, the alignment isn’t perfect. Clearly, the low biomass in
cornfields is a driving difference in the biomass-aware ordination,
which, as a result, should possibly be preferred in this case.

## AMF abundance in families

Test proportion of biomass across field types for each family. A similar
contrast with sequence abundance only returns largely NS results,
highlighting the incorrect inference when using compositional data to
conduct site-differential analysis.

For each AMF family, fit a Gamma GLM with log link predicting
biomass-weighted relative abundance by field type. We controlled the
false discovery rate across the four family-level tests using a
Benjamini–Hochberg (BH) false discovery rate correction. For pairwise
contrasts, within models showing a significant field-type effect,
compare estimated marginal means with Tukey adjustment for all pairwise
contrasts (emmeans).

### Glomeraceae

``` r
glom_lm <- lm(Glmrc_mass ~ field_type, data = amf_fam_ma)
par(mfrow = c(2,2))
plot(glom_lm) # variance non-constant among groups
```

![](resources/fungal_ecology_files/figure-gfm/unnamed-chunk-72-1.png)<!-- -->

``` r
distribution_prob(glom_lm)
```

    ## 
    ## 
    ## Distribution    p_Residuals
    ## -------------  ------------
    ## normal              0.50000
    ## cauchy              0.46875
    ## gamma               0.03125
    ## 
    ## 
    ## Distribution    p_Response
    ## -------------  -----------
    ## chi                0.28125
    ## gamma              0.21875
    ## uniform            0.12500

residuals normal, response chi/gamma…mean variance relationship
suggested

``` r
leveneTest(residuals(glom_lm) ~ amf_fam_ma$field_type)
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  2  4.3054 0.02643 *
    ##       22                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Levene’s p\<0.05, null of equal variance rejected

``` r
amf_fam_ma %>% 
  mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant"))) %>%
  group_by(field_type) %>%
  summarize(mean = mean(Glmrc_mass),
            cv = sd(Glmrc_mass) / mean) %>%
  mutate(across(mean:cv, ~ round(.x, 2))) %>%
  kable(format = "pandoc", caption = "Mean and CV relationship in groups")
```

| field_type |  mean |   cv |
|:-----------|------:|-----:|
| corn       |  3.25 | 0.39 |
| restored   | 25.21 | 0.49 |
| remnant    | 29.20 | 0.57 |

Mean and CV relationship in groups

Corroborates the mean/variance relationship. Proceed with log
transformation or gamma glm

``` r
glom_lm_log   <- lm(log(Glmrc_mass) ~ field_type, data = amf_fam_ma)
par(mfrow = c(2,2))
plot(glom_lm_log) # some improvement, leverage point still exists
```

![](resources/fungal_ecology_files/figure-gfm/unnamed-chunk-75-1.png)<!-- -->

``` r
ncvTest(glom_lm_log) # p=0.29, null of constant variance not rejected, model fit is improved
```

    ## Non-constant Variance Score Test 
    ## Variance formula: ~ fitted.values 
    ## Chisquare = 1.129022, Df = 1, p = 0.28798

Gamma glm to reduce leverage point…

``` r
glom_glm <- glm(Glmrc_mass ~ field_type, family = Gamma(link = "log"), data = amf_fam_ma)
```

Diagnostics

``` r
glom_glm_diag <- glm.diag(glom_glm)
glm.diag.plots(glom_glm, glom_glm_diag) # qqplot shows strong fit; no leverage >0.5
```

![](resources/fungal_ecology_files/figure-gfm/unnamed-chunk-77-1.png)<!-- -->

``` r
check_model(glom_glm) # corroborates
```

![](resources/fungal_ecology_files/figure-gfm/cm3-1.png)<!-- -->

``` r
performance::check_overdispersion(glom_glm) # not detected
```

    ## # Overdispersion test
    ## 
    ##  dispersion ratio = 1.096
    ##           p-value = 0.696

    ## No overdispersion detected.

Gamma glm is the best choice; no high-leverage point

Model results, group means, and post-hoc

``` r
anova(glom_glm) # Decline in residual deviance worth the cost in df
```

    ## Analysis of Deviance Table
    ## 
    ## Model: Gamma, link: log
    ## 
    ## Response: Glmrc_mass
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ## 
    ##            Df Deviance Resid. Df Resid. Dev      F    Pr(>F)    
    ## NULL                          24    17.1887                     
    ## field_type  2   11.237        22     5.9519 23.499 3.463e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
glom_em <- emmeans(glom_glm, ~ field_type, type = "response")
```

### Claroideoglomeraceae

``` r
clar_lm <- lm(Clrdg_mass ~ field_type, data = amf_fam_ma)
par(mfrow = c(2,2))
plot(clar_lm) # variance non-constant among groups
```

![](resources/fungal_ecology_files/figure-gfm/unnamed-chunk-79-1.png)<!-- -->

``` r
distribution_prob(clar_lm)
```

    ## 
    ## 
    ## Distribution    p_Residuals
    ## -------------  ------------
    ## cauchy              0.68750
    ## normal              0.15625
    ## tweedie             0.09375
    ## 
    ## 
    ## Distribution    p_Response
    ## -------------  -----------
    ## half-cauchy        0.34375
    ## chi                0.21875
    ## weibull            0.15625

residuals and response very long tailed

``` r
leveneTest(residuals(clar_lm) ~ amf_fam_ma$field_type)
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  0.9056 0.4189
    ##       22

Levene’s p\>0.05, null of equal variance not rejected Look to gamma to
reduce outlier and remain consistent with this model set

``` r
clar_glm <- glm(Clrdg_mass ~ field_type, family = Gamma(link = "log"), data = amf_fam_ma)
```

Diagnostics

``` r
clar_glm_diag <- glm.diag(clar_glm)
glm.diag.plots(clar_glm, clar_glm_diag) # qqplot shows strong fit; no outlier point
```

![](resources/fungal_ecology_files/figure-gfm/unnamed-chunk-82-1.png)<!-- -->

``` r
check_model(clar_glm) # corroborates
```

![](resources/fungal_ecology_files/figure-gfm/cm4-1.png)<!-- -->

``` r
performance::check_overdispersion(clar_glm) # not detected
```

    ## # Overdispersion test
    ## 
    ##  dispersion ratio = 2.205
    ##           p-value =  0.08

    ## No overdispersion detected.

Gamma glm is the best choice; no high-leverage point

Model results, group means, and post-hoc

``` r
anova(clar_glm) # Decline in residual deviance worth the cost in df
```

    ## Analysis of Deviance Table
    ## 
    ## Model: Gamma, link: log
    ## 
    ## Response: Clrdg_mass
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ## 
    ##            Df Deviance Resid. Df Resid. Dev     F    Pr(>F)    
    ## NULL                          24     41.280                    
    ## field_type  2   23.413        22     17.867 10.01 0.0008104 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
clar_em <- emmeans(clar_glm, ~ field_type, type = "response")
```

### Paraglomeraceae

``` r
para_lm <- lm(Clrdg_mass ~ field_type, data = amf_fam_ma)
par(mfrow = c(2,2))
plot(para_lm) # variance non-constant among groups
```

![](resources/fungal_ecology_files/figure-gfm/unnamed-chunk-84-1.png)<!-- -->

``` r
distribution_prob(para_lm)
```

    ## 
    ## 
    ## Distribution    p_Residuals
    ## -------------  ------------
    ## cauchy              0.68750
    ## normal              0.15625
    ## tweedie             0.09375
    ## 
    ## 
    ## Distribution    p_Response
    ## -------------  -----------
    ## half-cauchy        0.34375
    ## chi                0.21875
    ## weibull            0.15625

residuals and response very long tailed

``` r
leveneTest(residuals(para_lm) ~ amf_fam_ma$field_type)
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  0.9056 0.4189
    ##       22

Levene’s p\>0.05, null of equal variance not rejected Look to gamma to
reduce outlier and remain consistent with this model set

``` r
para_glm <- glm(Prglm_mass ~ field_type, family = Gamma(link = "log"), data = amf_fam_ma)
```

Diagnostics

``` r
para_glm_diag <- glm.diag(para_glm)
glm.diag.plots(para_glm, para_glm_diag) # qqplot shows strong fit; no outlier point
```

![](resources/fungal_ecology_files/figure-gfm/unnamed-chunk-87-1.png)<!-- -->

``` r
check_model(para_glm) # corroborates
```

![](resources/fungal_ecology_files/figure-gfm/cm5-1.png)<!-- -->

``` r
performance::check_overdispersion(para_glm) # not detected
```

    ## # Overdispersion test
    ## 
    ##  dispersion ratio = 1.929
    ##           p-value = 0.216

    ## No overdispersion detected.

Gamma glm is the best choice; no high-leverage point

Model results, group means, and post-hoc

``` r
anova(para_glm) # field_type not significant
```

    ## Analysis of Deviance Table
    ## 
    ## Model: Gamma, link: log
    ## 
    ## Response: Prglm_mass
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ## 
    ##            Df Deviance Resid. Df Resid. Dev     F Pr(>F)
    ## NULL                          24     39.286             
    ## field_type  2   9.2546        22     30.031 2.534 0.1022

``` r
para_em <- emmeans(para_glm, ~ field_type, type = "response")
```

Pairwise significance not appropriate due to overall variable NS

### Diversisporaceae

``` r
diver_glm <- glm(Dvrss_mass ~ field_type, family = Gamma(link = "log"), data = amf_fam_ma)
```

Diagnostics

``` r
diver_glm_diag <- glm.diag(diver_glm)
glm.diag.plots(diver_glm, diver_glm_diag) # qqplot shows strong fit; no outlier point
```

![](resources/fungal_ecology_files/figure-gfm/unnamed-chunk-90-1.png)<!-- -->

``` r
check_model(diver_glm) # corroborates
```

![](resources/fungal_ecology_files/figure-gfm/cm6-1.png)<!-- -->

``` r
performance::check_overdispersion(diver_glm) # not detected
```

    ## # Overdispersion test
    ## 
    ##  dispersion ratio = 1.239
    ##           p-value =  0.44

    ## No overdispersion detected.

Gamma glm is the best choice; no high-leverage point Model results,
group means, and post-hoc

``` r
anova(diver_glm) # Decline in residual deviance worth the cost in df
```

    ## Analysis of Deviance Table
    ## 
    ## Model: Gamma, link: log
    ## 
    ## Response: Dvrss_mass
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ## 
    ##            Df Deviance Resid. Df Resid. Dev      F   Pr(>F)   
    ## NULL                          24     26.307                   
    ## field_type  2   8.6471        22     17.660 5.8721 0.009046 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
diver_em <- emmeans(diver_glm, ~ field_type, type = "response")
```

### Gigasporaceae and others

Zeroes in data; comparison not warranted

### Results table

Combine later with confidence intervals and corrected p values

``` r
amf_fam_diff <- 
  amf_fam_ma %>% 
  pivot_longer(Glmrc_mass:Dvrss_mass, names_to = "family", values_to = "mass") %>% 
  group_by(field_type, family) %>% 
  summarize(mass = mean(mass), .groups = "drop") %>% 
  pivot_wider(names_from = field_type, values_from = mass) %>% 
  mutate(total = rowSums(across(where(is.numeric))), across(where(is.numeric), ~ round(.x, 2))) %>% 
  arrange(-total) %>% select(-total) %>% 
  left_join(
    list(
      Glmrc_mass = glom_em,
      Clrdg_mass = clar_em,
      Prglm_mass = para_em,
      Dvrss_mass = diver_em
    ) %>% map(\(df) df %>% as.data.frame() %>% select(field_type, lower.CL, upper.CL) %>% 
                pivot_wider(names_from = field_type, values_from = c(lower.CL, upper.CL), names_glue = "{field_type}_{.value}")) %>% 
      bind_rows(.id = "family") %>% 
      mutate(across(where(is.numeric), ~ round(.x, 2))),
    by = join_by(family)
  ) %>% 
  mutate(corn_ci = paste(corn_lower.CL, corn_upper.CL, sep = "–"),
         restored_ci = paste(restored_lower.CL, restored_upper.CL, sep = "–"),
         remnant_ci = paste(remnant_lower.CL, remnant_upper.CL, sep = "–")) %>% 
  select(family, corn, corn_ci, restored, restored_ci, remnant, remnant_ci)
kable(amf_fam_diff, format = "pandoc", caption = "Table: AMF mass in families across field types")
```

| family     | corn | corn_ci   | restored | restored_ci | remnant | remnant_ci  |
|:-----------|-----:|:----------|---------:|:------------|--------:|:------------|
| Glmrc_mass | 3.25 | 2.07–5.12 |    25.21 | 19.56–32.48 |   29.20 | 17.58–48.47 |
| Clrdg_mass | 0.15 | 0.05–0.4  |     4.56 | 2.6–7.98    |    3.94 | 1.28–12.1   |
| Prglm_mass | 0.29 | 0.08–1.01 |     1.70 | 0.84–3.42   |    0.82 | 0.2–3.34    |
| Dvrss_mass | 0.10 | 0.04–0.22 |     0.55 | 0.35–0.86   |    0.65 | 0.26–1.57   |

Table: AMF mass in families across field types

Pairwise contrasts Paraglomeraceae model was NS; no pairwise indicated

``` r
list(glom = glom_em, clar = clar_em, diver = diver_em) %>% 
  map(\(df) pairs(df, adjust = "tukey"))
```

    ## $glom
    ##  contrast           ratio     SE df null t.ratio p.value
    ##  corn / restored    0.129 0.0323 22    1  -8.172  <.0001
    ##  corn / remnant     0.111 0.0366 22    1  -6.689  <.0001
    ##  restored / remnant 0.863 0.2360 22    1  -0.537  0.8540
    ## 
    ## P value adjustment: tukey method for comparing a family of 3 estimates 
    ## Tests are performed on the log scale 
    ## 
    ## $clar
    ##  contrast            ratio     SE df null t.ratio p.value
    ##  corn / restored    0.0322 0.0178 22    1  -6.203  <.0001
    ##  corn / remnant     0.0372 0.0270 22    1  -4.538  0.0005
    ##  restored / remnant 1.1558 0.6990 22    1   0.240  0.9689
    ## 
    ## P value adjustment: tukey method for comparing a family of 3 estimates 
    ## Tests are performed on the log scale 
    ## 
    ## $diver
    ##  contrast           ratio     SE df null t.ratio p.value
    ##  corn / restored    0.178 0.0781 22    1  -3.930  0.0020
    ##  corn / remnant     0.151 0.0872 22    1  -3.279  0.0092
    ##  restored / remnant 0.853 0.4090 22    1  -0.332  0.9412
    ## 
    ## P value adjustment: tukey method for comparing a family of 3 estimates 
    ## Tests are performed on the log scale

Model results

``` r
list(glom = glom_glm, clar = clar_glm, para = para_glm, diver = diver_glm) %>% 
  map(\(df) df %>% anova() %>% as.data.frame()) %>% 
  bind_rows(.id = "family") %>% 
  mutate(p.adj = p.adjust(`Pr(>F)`, "fdr") %>% round(., 5))
```

    ##                family Df  Deviance Resid. Df Resid. Dev         F       Pr(>F)   p.adj
    ## NULL...1         glom NA        NA        24   17.18866        NA           NA      NA
    ## field_type...2   glom  2 11.236779        22    5.95188 23.499356 3.462571e-06 0.00001
    ## NULL...3         clar NA        NA        24   41.27991        NA           NA      NA
    ## field_type...4   clar  2 23.412599        22   17.86731 10.009596 8.104473e-04 0.00162
    ## NULL...5         para NA        NA        24   39.28586        NA           NA      NA
    ## field_type...6   para  2  9.254643        22   30.03122  2.534027 1.022403e-01 0.10224
    ## NULL...7        diver NA        NA        24   26.30731        NA           NA      NA
    ## field_type...8  diver  2  8.647118        22   17.66019  5.872076 9.046193e-03 0.01206

## Indicator species analysis

``` r
amf_ind <- inspan(spe=amf_avg_ma, meta=amf_meta, guild=NULL, site_dat=sites)
amf_ind %>% 
  select(field_type, family, taxon, starts_with("corn"), starts_with("restor"), starts_with("rem")) %>% 
  filter(taxon != "unidentified") %>% 
  arrange(field_type) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% 
  kable(format = "pandoc", caption = "Indicator species analysis results with avg biomass")
```

| field_type | family | taxon | corn_avg | corn_ci | restored_avg | restored_ci | remnant_avg | remnant_ci |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|
| corn | Glomeraceae | Glomus sp. | 0.085 | 0.164 | 0.000 | 0.000 | 0.000 | 0.000 |
| corn | Glomeraceae | Glomus Glo2 | 0.012 | 0.020 | 0.000 | 0.000 | 0.000 | 0.000 |
| corn | Glomeraceae | Glomus ORVIN GLO3D | 0.003 | 0.004 | 0.000 | 0.000 | 0.000 | 0.000 |
| remnant | Glomeraceae | Glomus perpusillum | 0.000 | 0.000 | 0.060 | 0.058 | 0.677 | 1.301 |
| remnant | Glomeraceae | Glomus MO-G15 | 0.016 | 0.020 | 1.785 | 0.680 | 3.058 | 2.363 |
| remnant | Glomeraceae | Glomus Glo3 | 0.000 | 0.000 | 0.001 | 0.001 | 0.113 | 0.209 |
| remnant | Acaulosporaceae | Acaulospora Acau2 | 0.000 | 0.000 | 0.002 | 0.004 | 0.048 | 0.070 |
| restored | Glomeraceae | Glomus MO-G7 | 0.004 | 0.006 | 0.548 | 0.209 | 0.328 | 0.246 |

Indicator species analysis results with avg biomass

None of the identified species seem relevant for further discussion…

## Constrained Analysis

Env covars processed in the ITS section (see above)

``` r
spe_amf_wi_resto <- amf_avg_ma %>%
  filter(field_name %in% rownames(env_expl)) %>%
  column_to_rownames(var = "field_name")

amf_mod_null <- dbrda(spe_amf_wi_resto ~ 1 + Condition(env_cov), data = env_expl, distance = "bray")
amf_mod_full <- dbrda(spe_amf_wi_resto ~ . + Condition(env_cov), data = env_expl, distance = "bray")
amf_mod_step <- ordistep(amf_mod_null,
                         scope = formula(amf_mod_full),
                         direction = "forward",
                         permutations = 1999,
                         trace = FALSE)
```

### Constrained Analysis Results

``` r
amf_mod_step
```

    ## 
    ## Call: dbrda(formula = spe_amf_wi_resto ~ Condition(env_cov) + gf_index + pH, data = env_expl, distance = "bray")
    ## 
    ##               Inertia Proportion Rank
    ## Total         2.37652    1.00000     
    ## Conditional   0.21725    0.09141    1
    ## Constrained   1.01747    0.42814    2
    ## Unconstrained 1.14180    0.48045    9
    ## 
    ## Inertia is squared Bray distance
    ## 
    ## Eigenvalues for constrained axes:
    ## dbRDA1 dbRDA2 
    ## 0.7036 0.3139 
    ## 
    ## Eigenvalues for unconstrained axes:
    ##    MDS1    MDS2    MDS3    MDS4    MDS5    MDS6    MDS7    MDS8    MDS9 
    ## 0.30831 0.23208 0.18162 0.12011 0.10056 0.08747 0.05770 0.03574 0.01822

``` r
(amf_mod_r2   <- RsquareAdj(amf_mod_step, permutations = 1999))
```

    ## $r.squared
    ## [1] 0.4281361
    ## 
    ## $adj.r.squared
    ## [1] 0.3505848

``` r
(amf_mod_glax <- anova(amf_mod_step, permutations = 1999))
```

    ## Permutation test for dbrda under reduced model
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## Model: dbrda(formula = spe_amf_wi_resto ~ Condition(env_cov) + gf_index + pH, data = env_expl, distance = "bray")
    ##          Df SumOfSqs    F Pr(>F)    
    ## Model     2   1.0175 4.01  5e-04 ***
    ## Residual  9   1.1418                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(amf_mod_inax <- anova(amf_mod_step, by = "axis", permutations = 1999))
```

    ## Permutation test for dbrda under reduced model
    ## Forward tests for axes
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## Model: dbrda(formula = spe_amf_wi_resto ~ Condition(env_cov) + gf_index + pH, data = env_expl, distance = "bray")
    ##          Df SumOfSqs      F Pr(>F)    
    ## dbRDA1    1  0.70357 5.5457 0.0005 ***
    ## dbRDA2    1  0.31390 2.7492 0.0015 ** 
    ## Residual  9  1.14180                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(amf_mod_axpct <- round(100 * amf_mod_step$CCA$eig / sum(amf_mod_step$CCA$eig), 1))
```

    ## dbRDA1 dbRDA2 
    ##   69.1   30.9

``` r
amf_mod_step$anova %>% 
  as.data.frame() %>% 
  mutate(p.adj = p.adjust(`Pr(>F)`, "fdr")) %>% 
  kable(, format = "pandoc")
```

|             |  Df |       AIC |        F | Pr(\>F) |  p.adj |
|-------------|----:|----------:|---------:|--------:|-------:|
| \+ gf_index |   1 | 10.566795 | 4.027540 |  0.0020 | 0.0040 |
| \+ pH       |   1 |  8.683317 | 3.133299 |  0.0085 | 0.0085 |

Based on permutation tests with n=1999 permutations, after accounting
for inter-site pairwise distance as a covariate, the model shows a
significant correlation between the site ordination on fungal
communities and the selected explanatory variables (p\<0.001). The first
two constrained axes are also significant (p\<0.001, p\<0.02). The
selected variables explain $R^{2}_{\text{Adj}}$=35.1 of the community
variation. Selected explanatory variables are pH and the grass-forb
index; see table for individual p values and statistics.

#### AMF constrained figure

Produce figure objects. Code for multipanel fig 6 is shown in the
saprotroph section.

``` r
amf_mod_step_eig <- round(amf_mod_step$CCA$eig * 100, 1)
amf_mod_scor <- scores(
  amf_mod_step,
  choices = c(1, 2),
  display = c("bp", "sites"),
  tidy = FALSE
)
amf_mod_scor_site <- amf_mod_scor$sites %>%
  data.frame() %>%
  rownames_to_column(var = "field_name") %>%
  left_join(sites, by = join_by(field_name))
amf_mod_scor_bp <- bind_rows(
  amf_mod_scor$biplot %>%
    data.frame() %>%
    rownames_to_column(var = "envvar") %>%
    mutate(envlabs = c(">forb", "pH")),
  data.frame(
    envvar = "gf_index",
    dbRDA1 = -amf_mod_scor$biplot["gf_index", 1],
    dbRDA2 = -amf_mod_scor$biplot["gf_index", 2],
    envlabs = ">grass")
) %>% 
  arrange(envvar, envlabs) %>% 
  mutate(
    origin = 0,
    m = dbRDA2 / dbRDA1,
    d = sqrt(dbRDA1^2 + dbRDA2^2),
    dadd = sqrt((max(dbRDA1)-min(dbRDA2))^2 + (max(dbRDA2)-min(dbRDA2))^2)*dadd_adj,
    labx = ((d+dadd)*cos(atan(m)))*(dbRDA1/abs(dbRDA1)),
    laby = ((d+dadd)*sin(atan(m)))*(dbRDA1/abs(dbRDA1)))
```

## AM fungi and plant community correlations

Data for these tests

``` r
amf_resto <- amf_div %>% 
  left_join(fa %>% select(field_name, amf_mass = amf), by = join_by(field_name)) %>% 
  left_join(sites, by = join_by(field_name, field_type)) %>% 
  left_join(gf_index, by = join_by(field_name)) %>% 
  left_join(prich %>% select(field_name, pl_rich, pl_shan), by = join_by(field_name)) %>% 
  filter(field_type != "corn", region != "FL") %>% 
  select(field_name, amf_ab = depth, amf_mass, gf_index, pl_rich, pl_shan) 
```

### Plant richness and fungal biomass

Is plant richness related to am fungal mass?

``` r
amfa_prich_lm <- lm(amf_mass ~ pl_rich, data = amf_resto)
summary(amfa_prich_lm)
```

    ## 
    ## Call:
    ## lm(formula = amf_mass ~ pl_rich, data = amf_resto)
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -25.68 -10.68  -1.42  10.06  45.05 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept) 35.651114  18.693669   1.907   0.0829 .
    ## pl_rich     -0.003654   0.453789  -0.008   0.9937  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 19.34 on 11 degrees of freedom
    ## Multiple R-squared:  5.893e-06,  Adjusted R-squared:  -0.0909 
    ## F-statistic: 6.483e-05 on 1 and 11 DF,  p-value: 0.9937

AM fungal mass and plant richness are nearly perfectly unrelated.

Is plant diversity related to am fungal mass?

``` r
amfa_pshan_lm <- lm(amf_mass ~ pl_shan, data = amf_resto)
summary(amfa_pshan_lm)
```

    ## 
    ## Call:
    ## lm(formula = amf_mass ~ pl_shan, data = amf_resto)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -25.120 -10.432  -2.183  10.726  43.314 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)  23.5363    23.2042   1.014    0.332
    ## pl_shan       0.9046     1.7072   0.530    0.607
    ## 
    ## Residual standard error: 19.1 on 11 degrees of freedom
    ## Multiple R-squared:  0.02489,    Adjusted R-squared:  -0.06376 
    ## F-statistic: 0.2808 on 1 and 11 DF,  p-value: 0.6067

AM fungal biomass and plant diversity are positively related but only
weakly so, no further testing warranted.

### AM fungal biomass and grass/forb composition

Inspect simple linear relationship. Naïve model.

``` r
amma_rest_m <- lm(amf_mass ~ gf_index, data = amf_resto)
```

``` r
check_model(amma_rest_m)
```

![](resources/fungal_ecology_files/figure-gfm/cm7-1.png)<!-- -->

``` r
summary(amma_rest_m)
```

    ## 
    ## Call:
    ## lm(formula = amf_mass ~ gf_index, data = amf_resto)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -15.847  -9.842  -2.154   5.266  45.036 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   35.507      4.551   7.801 8.29e-06 ***
    ## gf_index      35.320     17.063   2.070   0.0628 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 16.41 on 11 degrees of freedom
    ## Multiple R-squared:  0.2803, Adjusted R-squared:  0.2149 
    ## F-statistic: 4.285 on 1 and 11 DF,  p-value: 0.06277

AM fungal mass increases with grass-forb index slightly with p value
just above 0.95 alpha cutoff. Account for mass and gf_index as
explanatory factors for sequence abundance.

``` r
summary(lm(amf_ab ~ amf_mass + gf_index, data = amf_resto))
```

    ## 
    ## Call:
    ## lm(formula = amf_ab ~ amf_mass + gf_index, data = amf_resto)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -775.66 -580.92  -24.34  458.86 1072.20 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  2466.09     477.64   5.163 0.000423 ***
    ## amf_mass       36.43      12.38   2.942 0.014723 *  
    ## gf_index      422.36     825.83   0.511 0.620149    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 673.8 on 10 degrees of freedom
    ## Multiple R-squared:  0.5936, Adjusted R-squared:  0.5124 
    ## F-statistic: 7.304 on 2 and 10 DF,  p-value: 0.01108

AM fungal mass predicts increase in amf sequences, but after accounting
for that, gf_index doesn’t explain any new variation. Logistic
regression unavailable as AM fungal sequence abundance aren’t
proportional of a total.

# Putative plant pathogens

``` r
# Putative plant pathogens ———————— ####
```

Retrieve pathogen sequence abundance

``` r
patho <- guildseq(its_avg, its_meta, "plant_pathogen")
```

## Diversity Indices

``` r
patho_div <- calc_div(patho, sites) %>% 
  mutate(depth_csq = sqrt(depth) - mean(sqrt(depth)))
```

### Richness

Sequence depth square root transformed and centered. Negative binomial
model was underdispersed and failed to converge at default iterations;
use poisson glm instead.

Test interaction

``` r
patho_rich_glm_i <- glm(richness ~ depth_csq * field_type, data = patho_div, family = poisson(link = "log")) 
Anova(patho_rich_glm_i, type = 3, test.statistic = "LR") # no interaction detected
```

    ## Analysis of Deviance Table (Type III tests)
    ## 
    ## Response: richness
    ##                      LR Chisq Df Pr(>Chisq)
    ## depth_csq             2.01799  1     0.1554
    ## field_type            1.81313  2     0.4039
    ## depth_csq:field_type  0.35638  2     0.8368

Fit additive model

``` r
patho_rich_glm <- glm(richness ~ depth_csq + field_type, data = patho_div, family = poisson(link = "log")) 
```

Diagnostics

``` r
check_model(patho_rich_glm)
```

![](resources/fungal_ecology_files/figure-gfm/patho_rich_covar_diagnostics-1.png)<!-- -->

``` r
check_overdispersion(patho_rich_glm)
```

    ## # Overdispersion test
    ## 
    ##        dispersion ratio =  0.665
    ##   Pearson's Chi-Squared = 13.957
    ##                 p-value =  0.871

    ## No overdispersion detected.

``` r
check_collinearity(patho_rich_glm)
```

    ## # Check for Multicollinearity
    ## 
    ## Low Correlation
    ## 
    ##        Term  VIF    VIF 95% CI adj. VIF Tolerance Tolerance 95% CI
    ##   depth_csq 1.05 [1.00, 26.49]     1.03      0.95     [0.04, 1.00]
    ##  field_type 1.05 [1.00, 26.49]     1.01      0.95     [0.04, 1.00]

Some midrange structure, no leverage points, overdispersion, or
multicollinearity

``` r
distribution_prob(patho_rich_glm)
```

    ## 
    ## 
    ## Distribution    p_Residuals
    ## -------------  ------------
    ## normal              0.71875
    ## cauchy              0.18750
    ## gamma               0.06250
    ## 
    ## 
    ## Distribution                  p_Response
    ## ---------------------------  -----------
    ## beta-binomial                    0.50000
    ## neg. binomial (zero-infl.)       0.25000
    ## chi                              0.09375

residuals distribution normal or long-tailed, response count-distributed

``` r
leveneTest(richness ~ field_type, data = patho_div) %>% as.data.frame() %>% kable(format = "pandoc")
```

|       |  Df |  F value |   Pr(\>F) |
|-------|----:|---------:|----------:|
| group |   2 | 1.455231 | 0.2549471 |
|       |  22 |       NA |        NA |

``` r
leveneTest(residuals(patho_rich_glm) ~ patho_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc")
```

|       |  Df |   F value |   Pr(\>F) |
|-------|----:|----------:|----------:|
| group |   2 | 0.9728451 | 0.3936847 |
|       |  22 |        NA |        NA |

Residuals/response distributions do not suggest the need for
transformation. Levene’s p \> 0.05 → fail to reject = variances can be
considered equal across groups.

Model results, group means, and post-hoc. Use Type II LR test of
variables due to unbalanced design.

``` r
Anova(patho_rich_glm, type = 2, test.statistic = "LR")
```

    ## Analysis of Deviance Table (Type II tests)
    ## 
    ## Response: richness
    ##            LR Chisq Df Pr(>Chisq)    
    ## depth_csq   11.5413  1  0.0006807 ***
    ## field_type   1.5352  2  0.4641180    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Sequence depth is highly significant; richness doesn’t vary in groups.

``` r
patho_div %>% 
  group_by(field_type) %>% 
  summarize(across(c(depth, richness), ~ round(mean(.x), 0))) %>% 
  kable(format = "pandoc", caption = "Average sequence depth and pathogen richness in field types")
```

| field_type | depth | richness |
|:-----------|------:|---------:|
| corn       |  1295 |       39 |
| restored   |  1324 |       42 |
| remnant    |   979 |       37 |

Average sequence depth and pathogen richness in field types

Depth is correlated with richness in field types. Differences in
richness are small and with depth variance removed first, this explains
why richness isn’t significantly different.

Calculate confidence intervals for figure. Arithmetic means calculated
in this case.

``` r
patho_rich_em <- emmeans(patho_rich_glm, ~ field_type, type = "response")
```

| field_type |     rate |       SE |  df | asymp.LCL | asymp.UCL |
|:-----------|---------:|---------:|----:|----------:|----------:|
| corn       | 38.01625 | 2.743835 | Inf |  33.00149 |  43.79303 |
| restored   | 41.73822 | 1.627226 | Inf |  38.66772 |  45.05253 |
| remnant    | 39.20495 | 3.257339 | Inf |  33.31340 |  46.13844 |

Confidence level used: 0.95

| contrast           |     ratio |        SE |  df | null |    z.ratio |   p.value |
|:-------------------|----------:|----------:|----:|-----:|-----------:|----------:|
| corn / restored    | 0.9108259 | 0.0743125 | Inf |    1 | -1.1448178 | 0.4864567 |
| corn / remnant     | 0.9696798 | 0.1073621 | Inf |    1 | -0.2780858 | 0.9582676 |
| restored / remnant | 1.0646159 | 0.0988045 | Inf |    1 |  0.6746647 | 0.7782483 |

P value adjustment: tukey method for comparing a family of 3 estimates

``` r
patho_rich_fig <- 
  ggplot(summary(patho_rich_em), aes(x = field_type, y = rate)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = rate, ymax = asymp.UCL), width = 0, linewidth = lw) +
  labs(x = NULL, y = expression(atop("Richness", paste("(", italic(n), " OTUs)")))) +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1))
```

### Shannon’s diversity

Sequence depth square root transformed and centered

``` r
patho_shan_lm <- lm(shannon ~ depth_csq + field_type, data = patho_div)
```

Diagnostics

``` r
check_model(patho_shan_lm)
```

![](resources/fungal_ecology_files/figure-gfm/patho_shan_covar_diagnostics-1.png)<!-- -->

``` r
distribution_prob(patho_shan_lm)
```

    ## 
    ## 
    ## Distribution    p_Residuals
    ## -------------  ------------
    ## normal              0.84375
    ## cauchy              0.12500
    ## pareto              0.03125
    ## 
    ## 
    ## Distribution    p_Response
    ## -------------  -----------
    ## normal             0.40625
    ## weibull            0.15625
    ## pareto             0.12500

residuals distribution most likely cauchy/normal; symmetric but long
tails response normal

``` r
leveneTest(shannon ~ field_type, data = patho_div) %>% as.data.frame() %>% kable(format = "pandoc")
```

|       |  Df |   F value |   Pr(\>F) |
|-------|----:|----------:|----------:|
| group |   2 | 0.0061893 | 0.9938316 |
|       |  22 |        NA |        NA |

``` r
leveneTest(residuals(patho_shan_lm) ~ patho_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc")
```

|       |  Df |   F value |   Pr(\>F) |
|-------|----:|----------:|----------:|
| group |   2 | 0.1904989 | 0.8278957 |
|       |  22 |        NA |        NA |

Residuals distribution does not suggest the need for further model
selection. Levene’s p \> 0.05 → fail to reject = variances can be
considered equal.

Model results, group means, and post-hoc

``` r
Anova(patho_shan_lm, type = 2)
```

    ## Anova Table (Type II tests)
    ## 
    ## Response: shannon
    ##             Sum Sq Df F value  Pr(>F)  
    ## depth_csq   21.506  1  3.7988 0.06477 .
    ## field_type  15.342  2  1.3551 0.27958  
    ## Residuals  118.884 21                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Neither predictor is significant

``` r
patho_shan_em <- emmeans(patho_shan_lm, ~ field_type, type = "response")
```

Results tables below show the emmeans summary of group means and
confidence intervals, with sequencing depth as a covariate, and the post
hoc contrast of richness among field types.

| field_type |   emmean |        SE |  df |  lower.CL | upper.CL |
|:-----------|---------:|----------:|----:|----------:|---------:|
| corn       | 12.37989 | 1.0655872 |  21 | 10.163883 | 14.59590 |
| restored   | 10.37509 | 0.5965357 |  21 |  9.134528 | 11.61566 |
| remnant    | 10.74690 | 1.2159866 |  21 |  8.218121 | 13.27569 |

Confidence level used: 0.95

| contrast           |   estimate |       SE |  df |    t.ratio |   p.value |
|:-------------------|-----------:|---------:|----:|-----------:|----------:|
| corn - restored    |  2.0048015 | 1.219093 |  21 |  1.6445023 | 0.2497492 |
| corn - remnant     |  1.6329897 | 1.625668 |  21 |  1.0045038 | 0.5822155 |
| restored - remnant | -0.3718118 | 1.362784 |  21 | -0.2728324 | 0.9598747 |

P value adjustment: tukey method for comparing a family of 3 estimates

``` r
patho_shan_fig <- 
  ggplot(summary(patho_shan_em), aes(x = field_type, y = emmean)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = emmean, ymax = upper.CL), width = 0, linewidth = lw) +
  labs(x = "Field type", y = expression(Shannon~diversity~paste("(", italic(e)^italic(H), ")"))) +
  # lims(y = c(0, 160)) +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
```

## Abundance

``` r
patho_ma_lm <- lm(patho_mass ~ field_type, data = its_guild_ma)
par(mfrow = c(2,2))
plot(patho_ma_lm) 
```

![](resources/fungal_ecology_files/figure-gfm/unnamed-chunk-115-1.png)<!-- -->

no serious violations observed

``` r
distribution_prob(patho_ma_lm)
```

    ## 
    ## 
    ## Distribution    p_Residuals
    ## -------------  ------------
    ## normal              0.53125
    ## cauchy              0.15625
    ## gamma               0.12500
    ## 
    ## 
    ## Distribution    p_Response
    ## -------------  -----------
    ## gamma              0.43750
    ## exponential        0.09375
    ## F                  0.09375

Residuals distribution fits normal, response gamma?

``` r
leveneTest(residuals(patho_ma_lm) ~ fa$field_type) %>% as.data.frame() %>% kable(format = "pandoc") 
```

|       |  Df |  F value |   Pr(\>F) |
|-------|----:|---------:|----------:|
| group |   2 | 2.415849 | 0.1125951 |
|       |  22 |       NA |        NA |

No covariate, response and residuals tests equivalent. Residuals
distribution does not suggest the need for transformation. Levene’s p \>
0.05 → fail to reject = variances can be considered equal.

Model results, group means, and post-hoc, with arithmetic means from
emmeans

``` r
anova(patho_ma_lm)
```

    ## Analysis of Variance Table
    ## 
    ## Response: patho_mass
    ##            Df Sum Sq Mean Sq F value Pr(>F)
    ## field_type  2 0.4790 0.23950  1.7121 0.2037
    ## Residuals  22 3.0774 0.13988

``` r
patho_ma_em <- emmeans(patho_ma_lm, ~ field_type, type = "response")
```

| field_type |    emmean |        SE |  df |  lower.CL |  upper.CL |
|:-----------|----------:|----------:|----:|----------:|----------:|
| corn       | 0.4420484 | 0.1672623 |  22 | 0.0951678 | 0.7889291 |
| restored   | 0.7950284 | 0.0935024 |  22 | 0.6011162 | 0.9889406 |
| remnant    | 0.6749910 | 0.1870049 |  22 | 0.2871666 | 1.0628154 |

Confidence level used: 0.95

| contrast           |   estimate |        SE |  df |    t.ratio |   p.value |
|:-------------------|-----------:|----------:|----:|-----------:|----------:|
| corn - restored    | -0.3529799 | 0.1916230 |  22 | -1.8420542 | 0.1795254 |
| corn - remnant     | -0.2329426 | 0.2508934 |  22 | -0.9284524 | 0.6284990 |
| restored - remnant |  0.1200374 | 0.2090778 |  22 |  0.5741277 | 0.8351761 |

P value adjustment: tukey method for comparing a family of 3 estimates

``` r
patho_ma_fig <- 
  ggplot(summary(patho_ma_em), aes(x = field_type, y = emmean)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = emmean, ymax = upper.CL), width = 0, linewidth = lw) +
  labs(x = "Field type", y = expression(atop("Biomass (scaled)", paste(bold(`(`), "(", nmol[PLFA], " × ", g[soil]^{-1}, ")", " × ", paste("(rel. abund)", bold(`)`)))))) +
  # labs(x = "Field Type", y = "Biomass (scaled)") +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1.1))
```

## Beta Diversity

Community distances handled similarly to previous

### Biomass-weighted relative abundance

``` r
patho_ma <- guildseq(its_avg_ma, its_meta, "plant_pathogen")
```

``` r
d_patho_ma <- patho_ma %>% 
  data.frame(row.names = 1) %>% 
  vegdist("bray")
mva_patho_ma <- mva(d = d_patho_ma, env = sites, corr = "lingoes")
```

``` r
mva_patho_ma$dispersion_test
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)
    ## Groups     2 0.034719 0.017359 1.6491   1999  0.226
    ## Residuals 22 0.231584 0.010526                     
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##              corn  remnant restored
    ## corn              0.310000    0.086
    ## remnant  0.310250             0.622
    ## restored 0.080028 0.602352

``` r
mva_patho_ma$permanova
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## adonis2(formula = d ~ dist_axis_1 + field_type, data = env, permutations = nperm, by = "terms")
    ##             Df SumOfSqs      R2      F Pr(>F)   
    ## dist_axis_1  1   0.2740 0.06235 1.7945 0.0875 . 
    ## field_type   2   0.9141 0.20802 2.9935 0.0015 **
    ## Residual    21   3.2061 0.72964                 
    ## Total       24   4.3942 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
mva_patho_ma$pairwise_contrasts[c(1,3,2), c(1,2,4,3,8)] %>% 
  arrange(group1, desc(group2)) %>% 
  kable(format = "pandoc", caption = "Pairwise permanova contrasts")
```

| group1  | group2   | F_value |    R2 | p_value_adj |
|:--------|:---------|--------:|------:|------------:|
| corn    | restored |   5.335 | 0.205 |      0.0015 |
| corn    | remnant  |   3.772 | 0.360 |      0.0060 |
| remnant | restored |   0.581 | 0.030 |      0.8100 |

Pairwise permanova contrasts

Lingoes correction was needed. Three axes were significant based on
broken stick test. Based on the homogeneity of variance test, the null
hypothesis of equal variance among groups is accepted across all
clusters and in pairwise comparison of clusters (both p\>0.05),
supporting the application of a PERMANOVA test.

An effect of geographic distance (covariate) on pathogen communities was
not supported. With geographic distance accounted for, the test variable
‘field type’ significantly explained variation in fungal communities,
with a post-hoc test revealing that communities in corn fields differed
from communities in restored and remnant fields.

Plotting results:

``` r
patho_ma_ord_data <- mva_patho_ma$ordination_scores %>% mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant")))
p_patho_ma_centers <- patho_ma_ord_data %>% 
  group_by(field_type) %>% 
  summarize(across(starts_with("Axis"), list(mean = mean, ci_l = ci_l, ci_u = ci_u), .names = "{.fn}_{.col}"), .groups = "drop") %>% 
  mutate(across(c(ci_l_Axis.1, ci_u_Axis.1), ~ mean_Axis.1 + .x),
         across(c(ci_l_Axis.2, ci_u_Axis.2), ~ mean_Axis.2 + .x),
         across(ends_with("Axis.1"), ~ .x * -1)) # reverse axis values to be consistent with other plots
patho_ma_ord <- 
  ggplot(patho_ma_ord_data, aes(x = Axis.1 * -1, y = Axis.2)) + # reverse axis
  geom_linerange(data = p_patho_ma_centers, aes(x = mean_Axis.1, y = mean_Axis.2, xmin = ci_l_Axis.1, xmax = ci_u_Axis.1), linewidth = lw) +
  geom_linerange(data = p_patho_ma_centers, aes(x = mean_Axis.1, y = mean_Axis.2, ymin = ci_l_Axis.2, ymax = ci_u_Axis.2), linewidth = lw) +
  geom_point(data = p_patho_ma_centers, aes(x = mean_Axis.1, y = mean_Axis.2, fill = field_type), size = lg_size, stroke = lw, shape = 21) +
  geom_point(aes(fill = field_type), size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  scale_fill_manual(values = ft_pal) +
  labs(
    x = paste0("PCoA 1 (", mva_patho_ma$axis_pct[1], "%)"),
    y = paste0("PCoA 2 (", mva_patho_ma$axis_pct[2], "%)")) +
  theme_ord +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1))
```

#### Unified figure

``` r
fig4_ls <- (patho_rich_fig / plot_spacer() / patho_ma_fig) +
  plot_layout(heights = c(1,0.01,1)) 
fig4 <- (fig4_ls | plot_spacer() | patho_ma_ord) +
  plot_layout(widths = c(0.35, 0.01, 0.64)) +
  plot_annotation(tag_levels = 'A') 
```

``` r
fig4
```

![](resources/fungal_ecology_files/figure-gfm/fig4-1.png)<!-- -->

**Fig 4.** Putative plant pathogen communities in **corn**,
**restored**, and **remnant** prairie fields. **a** OTU richness and
**b** biomass (nmol PLFA g soil^-1 \* (proportional sequence abundance))
are shown as columns with 95 % CIs. **c** Principal-coordinate (PCoA)
ordination of ITS-based (97 % OTU) community distances: small points =
sites, large circles = field-type centroids (error bars = 95 % CI).
Cornfields cluster apart from restored or remnant prairies (P \< 0.01).
Numbers in black circles give years since restoration. Axis labels show
the percent variation explained.

### Sequence-based relative abundance

``` r
d_patho <- patho %>% 
  data.frame(row.names = 1) %>% 
  decostand("total") %>%
  vegdist("bray")
mva_patho <- mva(d = d_patho, env = sites, corr = "lingoes")
```

Diagnostics/results

``` r
mva_patho$dispersion_test
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
    ## Groups     2 0.015563 0.0077814 1.4865   1999  0.251
    ## Residuals 22 0.115164 0.0052347                     
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##             corn remnant restored
    ## corn             0.31450   0.1035
    ## remnant  0.30990           0.7800
    ## restored 0.10559 0.77814

``` r
mva_patho$permanova
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## adonis2(formula = d ~ dist_axis_1 + field_type, data = env, permutations = nperm, by = "terms")
    ##             Df SumOfSqs      R2      F Pr(>F)    
    ## dist_axis_1  1   0.1805 0.05290 1.6072 0.1345    
    ## field_type   2   0.8732 0.25589 3.8872 0.0005 ***
    ## Residual    21   2.3587 0.69121                  
    ## Total       24   3.4124 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
mva_patho$pairwise_contrasts[c(1,3,2), c(1,2,4,3,8)] %>% 
  arrange(group1, desc(group2)) %>% 
  kable(format = "pandoc", caption = "Pairwise permanova contrasts")
```

| group1  | group2   | F_value |    R2 | p_value_adj |
|:--------|:---------|--------:|------:|------------:|
| corn    | restored |   6.418 | 0.246 |      0.0015 |
| corn    | remnant  |   5.690 | 0.453 |      0.0098 |
| remnant | restored |   0.768 | 0.040 |      0.6310 |

Pairwise permanova contrasts

Lingoes correction was needed. Three axes were significant based on a
broken stick test. Based on the homogeneity of variance test, the null
hypothesis of equal variance among groups is accepted across all
clusters and in pairwise comparison of clusters (both p\>0.05),
supporting the application of a PERMANOVA test. An effect of geographic
distance (covariate) on pathogen communities was not supported. With
geographic distance accounted for, the test variable ‘field type’
significantly explained variation in fungal communities, with a post-hoc
test revealing that communities in corn fields differed from communities
in restored and remnant fields.

Plot results

``` r
patho_ord_data <- mva_patho$ordination_scores %>% mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant")))
p_patho_centers <- patho_ord_data %>% 
  group_by(field_type) %>% 
  summarize(across(starts_with("Axis"), list(mean = mean, ci_l = ci_l, ci_u = ci_u), .names = "{.fn}_{.col}"), .groups = "drop") %>% 
  mutate(across(c(ci_l_Axis.1, ci_u_Axis.1), ~ mean_Axis.1 + .x),
         across(c(ci_l_Axis.2, ci_u_Axis.2), ~ mean_Axis.2 + .x))
patho_ord <- 
  ggplot(patho_ord_data, aes(x = Axis.1, y = Axis.2)) +
  geom_linerange(data = p_patho_centers, aes(x = mean_Axis.1, y = mean_Axis.2, xmin = ci_l_Axis.1, xmax = ci_u_Axis.1), linewidth = lw) +
  geom_linerange(data = p_patho_centers, aes(x = mean_Axis.1, y = mean_Axis.2, ymin = ci_l_Axis.2, ymax = ci_u_Axis.2), linewidth = lw) +
  geom_point(data = p_patho_centers, aes(x = mean_Axis.1, y = mean_Axis.2, fill = field_type), size = lg_size, stroke = lw, shape = 21) +
  geom_point(aes(fill = field_type), size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "serif", fontface = 2, color = "black") +
  scale_fill_manual(values = ft_pal) +
  labs(
    x = paste0("PCoA 1 (", mva_patho$axis_pct[1], "%)"),
    y = paste0("PCoA 2 (", mva_patho$axis_pct[2], "%)")) +
  theme_ord +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1))
```

#### Supplemental figure

``` r
patho_shan_ord_sup <- (patho_shan_fig | plot_spacer() | patho_ord) +
  plot_layout(widths = c(0.45, 0.01, 0.55)) +
  plot_annotation(tag_levels = 'A') 
```

``` r
patho_shan_ord_sup
```

![](resources/fungal_ecology_files/figure-gfm/patho_shan_ord_sup-1.png)<!-- -->

#### Contrast community metrics

Procrustes test on PCoA values using axes with eigenvalues exceeding a
broken stick model

``` r
set.seed(20251112)
patho_protest <- protest(
  pcoa(d_patho)$vectors[, 1:3],
  pcoa(d_patho_ma)$vectors[, 1:3],
  permutations = 1999
)
patho_protest
```

    ## 
    ## Call:
    ## protest(X = pcoa(d_patho)$vectors[, 1:3], Y = pcoa(d_patho_ma)$vectors[,      1:3], permutations = 1999) 
    ## 
    ## Procrustes Sum of Squares (m12 squared):        0.2189 
    ## Correlation in a symmetric Procrustes rotation: 0.8838 
    ## Significance:  5e-04 
    ## 
    ## Permutation: free
    ## Number of permutations: 1999

Including biomass changes little. The spatial configurations of both
ordinations are highly correlated. $R^{2}=$ 0.78, p\<0.001.

## Pathogen Indicator Species

Use as a tool to find species for discussion. Unbalanced design and bias
to agricultural soil research may make the indicator stats less
appropriate for other use. Using biomass-weighted relative abundance for
this differential analysis.

``` r
patho_ind <- inspan(its_avg_ma, its_meta, "plant_pathogen", sites)
patho_ind %>% 
  select(A, B, stat, p_val_adj, field_type, species, starts_with("corn"), starts_with("restor"), starts_with("rem")) %>% 
  filter(species != "unidentified") %>% 
  arrange(field_type, p_val_adj) %>% 
  mutate(across(A:p_val_adj, ~ round(.x, 3)),
         across(corn_avg:remnant_ci, ~ num(.x, notation = "sci"))) %>% 
  kable(format = "pandoc", caption = "Indicator species analysis results with biomass-aware relative abundances in field types")
```

| A | B | stat | p_val_adj | field_type | species | corn_avg | corn_ci | restored_avg | restored_ci | remnant_avg | remnant_ci |
|---:|---:|---:|---:|:---|:---|---:|---:|---:|---:|---:|---:|
| 0.975 | 1.00 | 0.987 | 0.037 | corn | Corynespora_cassiicola | 1.28e-2 | 1.13e-2 | 3.28e-4 | 4.00e-4 | 0 | 0 |
| 0.866 | 1.00 | 0.930 | 0.037 | corn | Setophoma_terrestris | 7.14e-2 | 3.77e-2 | 8.47e-3 | 6.09e-3 | 2.62e-3 | 3.95e-3 |
| 1.000 | 0.60 | 0.775 | 0.103 | corn | Pseudocoleophoma_polygonicola | 2.61e-3 | 4.41e-3 | 0 | 0 | 0 | 0 |
| 0.635 | 1.00 | 0.797 | 0.410 | corn | Plectosphaerella_cucumerina | 6.45e-2 | 4.86e-2 | 2.31e-2 | 1.32e-2 | 1.40e-2 | 1.86e-2 |
| 0.749 | 0.75 | 0.749 | 0.318 | remnant | Ustilago_nunavutica | 1.11e-5 | 2.17e-5 | 1.19e-3 | 1.65e-3 | 3.59e-3 | 4.32e-3 |
| 0.912 | 0.50 | 0.675 | 0.318 | remnant | Monosporascus_eutypoides | 0 | 0 | 2.30e-5 | 4.51e-5 | 2.38e-4 | 3.08e-4 |

Indicator species analysis results with biomass-aware relative
abundances in field types

## Constrained Analysis

Env covars processed in the ITS section (see above)

``` r
spe_patho_wi_resto <- patho_ma %>%
  filter(field_name %in% rownames(env_expl)) %>%
  column_to_rownames(var = "field_name")

patho_mod_null <- dbrda(spe_patho_wi_resto ~ 1 + Condition(env_cov), data = env_expl, distance = "bray")
patho_mod_full <- dbrda(spe_patho_wi_resto ~ . + Condition(env_cov), data = env_expl, distance = "bray")
patho_mod_step <- ordistep(patho_mod_null,
                           scope = formula(patho_mod_full),
                           direction = "forward",
                           permutations = 1999,
                           trace = FALSE)
```

### Constrained Analysis Results

``` r
patho_mod_step
```

    ## 
    ## Call: dbrda(formula = spe_patho_wi_resto ~ 1 + Condition(env_cov), data = env_expl, distance = "bray")
    ## 
    ##               Inertia Proportion Rank RealDims
    ## Total          2.1779     1.0000              
    ## Conditional    0.2519     0.1157    1         
    ## Unconstrained  1.9260     0.8843   11       10
    ## 
    ## Inertia is squared Bray distance
    ## 
    ## Eigenvalues for unconstrained axes:
    ##    MDS1    MDS2    MDS3    MDS4    MDS5    MDS6    MDS7    MDS8    MDS9   MDS10   iMDS1 
    ##  0.7992  0.5233  0.2020  0.1376  0.0802  0.0734  0.0462  0.0435  0.0221  0.0099 -0.0114

``` r
(patho_mod_r2   <- RsquareAdj(patho_mod_step, permutations = 1999))
```

    ## $r.squared
    ## numeric(0)
    ## 
    ## $adj.r.squared
    ## numeric(0)

``` r
(patho_mod_glax <- anova(patho_mod_step, permutations = 1999))
```

    ## No constrained component
    ## 
    ## Model: dbrda(formula = spe_patho_wi_resto ~ 1 + Condition(env_cov), data = env_expl, distance = "bray")
    ##          Df SumOfSqs  F Pr(>F)
    ## Model     0    0.000  0       
    ## Residual 12    1.926

``` r
(patho_mod_inax <- anova(patho_mod_step, by = "axis", permutations = 1999))
```

    ## No constrained component
    ## 
    ## Model: dbrda(formula = spe_patho_wi_resto ~ 1 + Condition(env_cov), data = env_expl, distance = "bray")
    ##          Df SumOfSqs  F Pr(>F)
    ## Model     0    0.000  0       
    ## Residual 12    1.926

``` r
(patho_mod_axpct <- round(100 * patho_mod_step$CCA$eig / sum(patho_mod_step$CCA$eig), 1))
```

    ## numeric(0)

Based on permutation tests with n=1999 permutations, after accounting
for inter-site pairwise distance as a covariate, the model shows no
significant correlation between pathogen community turnover and
explanatory variables.

#### Pathogen (un)constrained figure

Produce figure objects. In this case the panel will show a PCoA
ordination of sites based on pathogen communities. Code for multipanel
fig 6 is shown in the saprotroph section.

``` r
patho_mod_eig <- round(patho_mod_step$CA$eig / sum(patho_mod_step$CA$eig) * 100, 1)[1:2]
patho_mod_scor <- scores(
  patho_mod_step,
  choices = c(1, 2),
  display = c("sites"),
  tidy = FALSE
)
patho_mod_scor_site <- patho_mod_scor %>%
  data.frame() %>%
  rownames_to_column(var = "field_name") %>%
  left_join(sites, by = join_by(field_name))
```

## Pathogen and plant community correlations

Data for these tests

``` r
patho_resto <- its_guild %>% 
  filter(field_type != "corn", region != "FL") %>% 
  left_join(its_guild_ma %>% select(field_name, patho_mass), by = join_by(field_name)) %>% 
  left_join(prich %>% select(field_name, pl_rich, pl_shan), by = join_by(field_name)) %>% 
  mutate(
    patho_prop = patho_abund / fungi_abund, # no zeroes present...
    notpatho_abund = fungi_abund - patho_abund,
    fungi_mass_lc = as.numeric(scale(log(fungi_mass), center = TRUE, scale = FALSE))
  ) %>% 
  select(-sapro_abund, -c(annual:shrubTree)) 
```

### Plant richness and pathogen biomass

Is plant richness related to pathogen mass?

``` r
pathofa_prich_lm = lm(patho_mass ~ pl_rich, data = patho_resto)
summary(pathofa_prich_lm)
```

    ## 
    ## Call:
    ## lm(formula = patho_mass ~ pl_rich, data = patho_resto)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.59887 -0.36293  0.02093  0.28457  0.79546 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)  1.126277   0.426049   2.644   0.0228 *
    ## pl_rich     -0.006329   0.010342  -0.612   0.5530  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4409 on 11 degrees of freedom
    ## Multiple R-squared:  0.03292,    Adjusted R-squared:  -0.055 
    ## F-statistic: 0.3744 on 1 and 11 DF,  p-value: 0.553

Pathogen mass and plant richness aren’t correlated, though the direction
is negative. Relationship is weak enough that no further tests are
warranted.

Is plant diversity related to pathogen mass?

``` r
pathofa_pshan_lm = lm(patho_mass ~ pl_shan, data = patho_resto)
summary(pathofa_pshan_lm)
```

    ## 
    ## Call:
    ## lm(formula = patho_mass ~ pl_shan, data = patho_resto)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.72208 -0.33569  0.08258  0.30109  0.72374 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)  1.11070    0.53974   2.058   0.0641 .
    ## pl_shan     -0.01770    0.03971  -0.446   0.6645  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4443 on 11 degrees of freedom
    ## Multiple R-squared:  0.01773,    Adjusted R-squared:  -0.07157 
    ## F-statistic: 0.1986 on 1 and 11 DF,  p-value: 0.6645

NS

### Pathogen biomass and grass/forb composition

#### Inspect simple linear relationship. Naïve model.

``` r
patho_gf_lm <- lm(patho_mass ~ gf_index, data = patho_resto)
```

``` r
check_model(patho_gf_lm)
```

![](resources/fungal_ecology_files/figure-gfm/cm9-1.png)<!-- -->

``` r
summary(patho_gf_lm)
```

    ## 
    ## Call:
    ## lm(formula = patho_mass ~ gf_index, data = patho_resto)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.5889 -0.3432  0.1578  0.2245  0.5870 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   0.8765     0.1096   7.997 6.55e-06 ***
    ## gf_index      0.7301     0.4109   1.777    0.103    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3952 on 11 degrees of freedom
    ## Multiple R-squared:  0.223,  Adjusted R-squared:  0.1524 
    ## F-statistic: 3.157 on 1 and 11 DF,  p-value: 0.1032

The naïve model shows a positive relationship that isn’t significant.
Diagnostic reveals noisy fit and lots of structure. How much does fungal
mass vary across sites?

``` r
(fma_cv <- sd(patho_resto$fungi_mass) / mean(patho_resto$fungi_mass) * 100)
```

    ## [1] 33.02084

Decompose biomass and sequence abundance in a model to test changes in
each given the other.

#### Multiple, weighted logistic glm

Note on interpretation: exponentiated coefficients are interpreted as
odds ratios for pathogen dominance within the fungal community. Sequence
abundances are used as analytic weights so that sites with higher
sequencing depth contributed proportionally more information to the
likelihood.

``` r
patho_gf_glm <- glm(patho_prop ~ fungi_mass_lc + gf_index,
                    data = patho_resto, family = quasibinomial(link = "logit"),
                    weights = fungi_abund)
summary(patho_gf_glm) # dispersion parameter >117 justifies quasibinomial
```

    ## 
    ## Call:
    ## glm(formula = patho_prop ~ fungi_mass_lc + gf_index, family = quasibinomial(link = "logit"), 
    ##     data = patho_resto, weights = fungi_abund)
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   -1.60740    0.09286 -17.310 8.77e-09 ***
    ## fungi_mass_lc  0.57428    0.29444   1.950 0.079698 .  
    ## gf_index       1.90791    0.38447   4.962 0.000568 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for quasibinomial family taken to be 117.7646)
    ## 
    ##     Null deviance: 4421.9  on 12  degrees of freedom
    ## Residual deviance: 1201.0  on 10  degrees of freedom
    ## AIC: NA
    ## 
    ## Number of Fisher Scoring iterations: 4

Diagnostics

``` r
check_model(patho_gf_glm)
```

![](resources/fungal_ecology_files/figure-gfm/unnamed-chunk-133-1.png)<!-- -->

``` r
check_collinearity(patho_gf_glm)
```

    ## # Check for Multicollinearity
    ## 
    ## Low Correlation
    ## 
    ##           Term  VIF    VIF 95% CI adj. VIF Tolerance Tolerance 95% CI
    ##  fungi_mass_lc 1.08 [1.00, 12.78]     1.04      0.93     [0.08, 1.00]
    ##       gf_index 1.08 [1.00, 12.78]     1.04      0.93     [0.08, 1.00]

``` r
augment(patho_gf_glm)
```

    ## # A tibble: 13 × 10
    ##    patho_prop fungi_mass_lc  gf_index `(weights)` .fitted  .resid   .hat .sigma  .cooksd .std.resid
    ##         <dbl>         <dbl>     <dbl>       <dbl>   <dbl>   <dbl>  <dbl>  <dbl>    <dbl>      <dbl>
    ##  1     0.126        0.470   -0.150          6936.  -1.62   -8.91  0.229   11.0  0.0807      -0.935 
    ##  2     0.217        0.142    0.114          7578.  -1.31    0.886 0.109   11.5  0.000306     0.0865
    ##  3     0.140        0.110    0.167          7300.  -1.23  -18.7   0.116    9.46 0.132       -1.83  
    ##  4     0.207        0.0610  -0.000380       8145.  -1.57    8.15  0.0856  11.2  0.0202       0.785 
    ##  5     0.0906       0.468   -0.575          7226.  -2.44    3.12  0.323   11.5  0.0201       0.349 
    ##  6     0.168        0.320   -0.222          8810.  -1.85    8.44  0.213   11.1  0.0736       0.877 
    ##  7     0.266       -0.540    0.267          8343.  -1.41   15.3   0.356    9.65 0.612        1.75  
    ##  8     0.133       -0.511    0.120          8853.  -1.67   -6.69  0.290   11.2  0.0696      -0.731 
    ##  9     0.0549      -0.439   -0.351          7392.  -2.53   -6.50  0.269   11.3  0.0552      -0.700 
    ## 10     0.0837      -0.339   -0.137         10437.  -2.06   -9.74  0.255   10.9  0.113       -1.04  
    ## 11     0.273        0.212    0.168          9062.  -1.16    7.66  0.198   11.2  0.0528       0.788 
    ## 12     0.281        0.00271  0.450          8537.  -0.747  -8.14  0.450   11.0  0.272       -1.01  
    ## 13     0.257        0.0438   0.150          8324.  -1.30    9.18  0.107   11.1  0.0332       0.895

Long tails and low n showing structure. Moderate leverage at LPRP1: high
pathogens, high gf_index but very low biomass…this is evidence of the
noise that caused the naïve model to fail.

``` r
distribution_prob(patho_gf_glm)
```

    ## 
    ## 
    ## Distribution    p_Residuals
    ## -------------  ------------
    ## normal              0.68750
    ## cauchy              0.12500
    ## gamma               0.09375
    ## 
    ## 
    ## Distribution    p_Response
    ## -------------  -----------
    ## weibull            0.31250
    ## uniform            0.18750
    ## beta               0.15625

Residuals distribution normal

``` r
loocv_paglm_gfi <- map_dbl(seq_len(nrow(patho_resto)), function(i){
  exp(coef(glm(patho_prop ~ fungi_mass_lc + gf_index, 
           data = patho_resto[-i, ], 
           family = quasibinomial(link = "logit"),
           weights = fungi_abund))["gf_index"])
})
summary(loocv_paglm_gfi)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   5.860   6.301   6.584   6.827   7.272   9.095

``` r
(cv_paglm <- (sd(loocv_paglm_gfi) / mean(loocv_paglm_gfi) * 100) %>% round(., 1))
```

    ## [1] 12.6

Grass-forb index LOOCV variation of 12.6% on the back-transformed scale
shows that the influential points (LPRP1) and three other potential
outliers from the qq plot do not significantly affect fit. Sign and
magnitude of LOO slopes wouldn’t change inference.

``` r
loocv_paglm_fma <- map_dbl(seq_len(nrow(patho_resto)), function(i){
  exp(coef(glm(patho_prop ~ fungi_mass_lc + gf_index, 
               data = patho_resto[-i, ], 
               family = quasibinomial(link = "logit"),
               weights = fungi_abund))["fungi_mass_lc"])
})
summary(loocv_paglm_fma)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   1.577   1.652   1.742   1.797   1.894   2.457

``` r
(cv_paglm <- (sd(loocv_paglm_fma) / mean(loocv_paglm_fma) * 100) %>% round(., 1))
```

    ## [1] 12.9

Similarly, fungal mass LOOCV variation of 12.9% on the back-transformed
scale shows that the influential points (LPRP1) and three other
potential outliers from the qq plot do not significantly affect fit.
Sign and magnitude of LOO slopes wouldn’t change inference. The higher
variability here shows that the noise of PLFA variation is substantial.
View partial regression plots for consistency and save data for
plotting.

``` r
paglm_crpldata <- as.data.frame(crPlots(patho_gf_glm))
```

![](resources/fungal_ecology_files/figure-gfm/unnamed-chunk-137-1.png)<!-- -->

Noise in fungal mass data is obvious here. Fit of partial gf_index is
clean.

Partial R2 values

``` r
data.frame(
  term = c("fungal_mass", "gf_index"),
  partial_R2 = rsq.partial(patho_gf_glm, adj = TRUE)$partial.rsq
) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% 
  kable(format = "pandoc", caption = "Partial R2 from weighted logistic regression")
```

| term        | partial_R2 |
|:------------|-----------:|
| fungal_mass |      0.220 |
| gf_index    |      0.713 |

Partial R2 from weighted logistic regression

Model summary

``` r
patho_null_glm <- glm(patho_prop ~ 1,
                  data = patho_resto, family = quasibinomial(link = "logit"),
                  weights = fungi_abund)
anova(patho_null_glm, patho_gf_glm, test = "F")
```

    ## Analysis of Deviance Table
    ## 
    ## Model 1: patho_prop ~ 1
    ## Model 2: patho_prop ~ fungi_mass_lc + gf_index
    ##   Resid. Df Resid. Dev Df Deviance      F   Pr(>F)   
    ## 1        12     4421.9                               
    ## 2        10     1201.0  2   3220.9 13.675 0.001376 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Deviance explained

``` r
(patho_gf_glm_pr2 <- round(1-(summary(patho_gf_glm)$deviance / summary(patho_gf_glm)$null.deviance), 3))
```

    ## [1] 0.728

Summary of terms

``` r
tidy(patho_gf_glm) %>% 
  mutate(odds_ratio = exp(estimate), exp_std.error = exp(std.error),
         across(where(is.numeric), ~ round(.x, 3))) %>% 
  select(term, estimate, odds_ratio, std.error, exp_std.error, statistic, p.value) %>% 
  kable(format = "pandoc", caption = "Summary of terms from weighted logistic regression")
```

| term          | estimate | odds_ratio | std.error | exp_std.error | statistic | p.value |
|:--------------|---------:|-----------:|----------:|--------------:|----------:|--------:|
| (Intercept)   |   -1.607 |      0.200 |     0.093 |         1.097 |   -17.310 |   0.000 |
| fungi_mass_lc |    0.574 |      1.776 |     0.294 |         1.342 |     1.950 |   0.080 |
| gf_index      |    1.908 |      6.739 |     0.384 |         1.469 |     4.962 |   0.001 |

Summary of terms from weighted logistic regression

``` r
(patho_or_pct <- round((exp(coef(patho_gf_glm)[3])^0.1)-1, 3))
```

    ## gf_index 
    ##     0.21

Confidence intervals on the prediction scale, results on the increment
of 0.1 increase in grass-forb index desired due to scale of that
variable. Note: in the following output, percent predicted changes are
calculated *in excess* of 100% (e.g., 1.124 = 12.4%).

``` r
((exp(confint(patho_gf_glm))^0.1)-1) %>%
  as.data.frame() %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% 
  kable(format = "pandoc", caption = "95% confidence intervals, back transformed from the log-scale")
```

|               |  2.5 % | 97.5 % |
|---------------|-------:|-------:|
| (Intercept)   | -0.164 | -0.133 |
| fungi_mass_lc |  0.001 |  0.123 |
| gf_index      |  0.124 |  0.307 |

95% confidence intervals, back transformed from the log-scale

Create objects for plotting

``` r
paglm_med_fungi <- median(patho_resto$fungi_mass_lc, na.rm = TRUE)
paglm_med_abund <- median(patho_resto$fungi_abund, na.rm = TRUE) # Needed for weight context
paglm_newdat <- tibble(
  gf_index = seq(min(patho_resto$gf_index, na.rm = TRUE),
                 max(patho_resto$gf_index, na.rm = TRUE),
                 length.out = 200),
  fungi_mass_lc = paglm_med_fungi,
  fungi_abund = paglm_med_abund 
)
```

Predict on link scale, back-transform with plogis

``` r
paglm_pred <- predict(patho_gf_glm, newdata = paglm_newdat, type = "link", se.fit = TRUE) %>%
  as_tibble() %>%
  bind_cols(paglm_newdat) %>%
  mutate(
    fit_prob = plogis(fit),
    lwr_prob = plogis(fit - 1.96 * se.fit),
    upr_prob = plogis(fit + 1.96 * se.fit)
  )
```

``` r
fig7a <- 
  ggplot(paglm_pred, aes(x = gf_index, y = fit_prob)) +
  # geom_ribbon(aes(ymin = lwr_med, ymax = upr_med), fill = "gray90") +
  geom_line(color = "black", linewidth = lw) +
  geom_point(data = patho_resto, aes(x = gf_index, y = patho_prop, fill = field_type),
             size = sm_size, stroke = lw, shape = 21) +
  geom_text(data = patho_resto, aes(x = gf_index, y = patho_prop, label = yr_since), 
            size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  labs(
    x = "Grass–forb index",
    y = "Pathogen proportion (share of total sequences)",
    tag = "A"
  ) +
  scale_fill_manual(name = "Field type", values = ft_pal[2:3]) +
  theme_cor +
  theme(legend.position = c(0.03, 1),
        legend.justification = c(0, 1),
        legend.title = element_text(size = 9, face = 1),
        legend.text = element_text(size = 8, face = 1),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2),
        legend.key = element_rect(fill = "white"),
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1))
```

``` r
fig7b <- 
  cbind(paglm_crpldata, patho_resto %>% select(field_type, yr_since)) %>% 
  ggplot(aes(x = fungi_mass_lc.fungi_mass_lc, y = fungi_mass_lc.patho_prop)) +
  geom_smooth(method = "lm", color = "black", linewidth = lw, se = FALSE) + 
  geom_point(aes(fill = field_type), size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  labs(x = expression("Fungal biomass"~paste("(", nmol[PLFA], " × ", g[soil]^{-1}, ")")), y = NULL, tag = "B") +
  scale_fill_manual(values = ft_pal[2:3]) +
  scale_y_continuous(breaks = c(-0.5, 0, 0.5)) +
  scale_x_continuous(breaks = log(c(2.7, 3.8, 5.3, 7.5)) - mean(log(patho_resto$fungi_mass)), 
                     labels = c(2.7, 3.8, 5.3, 7.5)) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0.175, 1))
```

``` r
fig7c <- 
  cbind(paglm_crpldata, patho_resto %>% select(field_type, yr_since)) %>% 
  ggplot(aes(x = gf_index.gf_index, y = gf_index.patho_prop)) +
  geom_smooth(method = "lm", color = "black", linewidth = lw, se = FALSE) + 
  geom_point(aes(fill = field_type), size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  labs(x = "Grass–forb index", y = NULL, tag = "C") +
  scale_fill_manual(values = ft_pal[2:3]) +
  scale_y_continuous(breaks = c(-1, 0, 1)) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0.175, 1))
fig7yax_grob <- textGrob(
  "Pathogen proportion (partial residuals, logit scale)",
  x = 0.5,
  y = 0.5,
  hjust = 0.5,
  vjust = 0.5,
  rot = 90,
  gp = gpar(cex = 9/12)
)
```

``` r
fig7rh <- (fig7b / plot_spacer() / fig7c) +
  plot_layout(heights = c(0.5, 0.01,0.5))
fig7 <- (fig7a | plot_spacer() | (wrap_elements(full = fig7yax_grob) & theme(plot.tag = element_blank())) | fig7rh) +
  plot_layout(widths = c(0.62, 0.004, 0.02, 0.38))
```

``` r
fig7
```

![](resources/fungal_ecology_files/figure-gfm/fig7-1.png)<!-- -->

# Putative saprotrophs

``` r
# Putative saprotrophs ———————— ####
```

Retrieve pathogen sequence abundance

``` r
sapro <- guildseq(its_avg, its_meta, "saprotroph")
```

## Diversity Indices

``` r
sapro_div <- calc_div(sapro, sites) %>% 
  mutate(depth_csq = sqrt(depth) - mean(sqrt(depth)))
```

### Richness

Sequence depth square root transformed and centered. Poisson model was
overdispersed (not shown), use negative binomial instead.

Test interaction

``` r
sapro_rich_glm_i <- glm.nb(richness ~ depth_csq * field_type, data = sapro_div) 
Anova(sapro_rich_glm_i, type = 3, test.statistic = "LR") # interaction detected
```

    ## Analysis of Deviance Table (Type III tests)
    ## 
    ## Response: richness
    ##                      LR Chisq Df Pr(>Chisq)   
    ## depth_csq              2.4867  1   0.114816   
    ## field_type             3.3102  2   0.191071   
    ## depth_csq:field_type  12.2617  2   0.002175 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
check_model(sapro_rich_glm_i)
```

![](resources/fungal_ecology_files/figure-gfm/unnamed-chunk-143-1.png)<!-- -->

``` r
check_overdispersion(sapro_rich_glm_i) # not overdispersed
```

    ## # Overdispersion test
    ## 
    ##  dispersion ratio = 1.056
    ##           p-value = 0.792

    ## No overdispersion detected.

``` r
augment(sapro_rich_glm_i) # corn site has cooks >0.9
```

    ## # A tibble: 25 × 9
    ##    richness depth_csq field_type .fitted  .resid   .hat .sigma  .cooksd .std.resid
    ##       <int>     <dbl> <fct>        <dbl>   <dbl>  <dbl>  <dbl>    <dbl>      <dbl>
    ##  1      100    -7.23  restored      4.64 -0.306  0.339    1.16 0.0119      -0.376 
    ##  2      145     0.850 restored      4.82  1.76   0.0821   1.09 0.0533       1.83  
    ##  3      118    -0.233 corn          4.83 -0.622  0.776    1.13 0.972       -1.31  
    ##  4       96   -15.2   remnant       4.56  0.0285 0.615    1.17 0.000564     0.0460
    ##  5      106    -4.15  restored      4.71 -0.415  0.144    1.16 0.00553     -0.448 
    ##  6       93    12.5   corn          4.64 -1.05   0.392    1.12 0.186       -1.34  
    ##  7      106     8.88  corn          4.70 -0.333  0.214    1.16 0.00633     -0.376 
    ##  8      120    -7.19  remnant       4.75  0.382  0.260    1.16 0.0118       0.445 
    ##  9      111    -3.14  restored      4.73 -0.187  0.103    1.17 0.000747    -0.198 
    ## 10      116    -6.67  restored      4.65  1.04   0.296    1.13 0.113        1.24  
    ## # ℹ 15 more rows

``` r
check_collinearity(sapro_rich_glm_i) # depth and interaction VIF > 6
```

    ## # Check for Multicollinearity
    ## 
    ## Low Correlation
    ## 
    ##        Term  VIF    VIF 95% CI adj. VIF Tolerance Tolerance 95% CI
    ##  field_type 4.32 [2.80,  7.14]     1.44      0.23     [0.14, 0.36]
    ## 
    ## Moderate Correlation
    ## 
    ##                  Term  VIF    VIF 95% CI adj. VIF Tolerance Tolerance 95% CI
    ##             depth_csq 7.41 [4.61, 12.40]     2.72      0.13     [0.08, 0.22]
    ##  depth_csq:field_type 6.01 [3.78, 10.01]     2.45      0.17     [0.10, 0.26]

An interaction was detected, but including it in the model leads to very
poor diagnostics. It’s driven by one site in corn with high leverage,
and it introduces high multicollinearity.

Fit additive model

``` r
sapro_rich_glm <- glm.nb(richness ~ depth_csq + field_type, data = sapro_div) 
```

Diagnostics

``` r
check_model(sapro_rich_glm)
```

![](resources/fungal_ecology_files/figure-gfm/sapro_rich_covar_diagnostics-1.png)<!-- -->

``` r
check_overdispersion(sapro_rich_glm)
```

    ## # Overdispersion test
    ## 
    ##  dispersion ratio = 1.041
    ##           p-value = 0.824

    ## No overdispersion detected.

``` r
check_collinearity(sapro_rich_glm)
```

    ## # Check for Multicollinearity
    ## 
    ## Low Correlation
    ## 
    ##        Term  VIF   VIF 95% CI adj. VIF Tolerance Tolerance 95% CI
    ##   depth_csq 2.04 [1.43, 3.50]     1.43      0.49     [0.29, 0.70]
    ##  field_type 2.04 [1.43, 3.50]     1.19      0.49     [0.29, 0.70]

Long tails, some structure throughout, no leverage points,
overdispersion, or multicollinearity

``` r
distribution_prob(sapro_rich_glm)
```

    ## 
    ## 
    ## Distribution    p_Residuals
    ## -------------  ------------
    ## normal              0.56250
    ## cauchy              0.18750
    ## gamma               0.09375
    ## 
    ## 
    ## Distribution                  p_Response
    ## ---------------------------  -----------
    ## neg. binomial (zero-infl.)       0.37500
    ## beta-binomial                    0.34375
    ## weibull                          0.09375

residuals distribution normal or long-tailed, response count-distributed

``` r
leveneTest(richness ~ field_type, data = sapro_div) %>% as.data.frame() %>% kable(format = "pandoc")
```

|       |  Df |   F value |   Pr(\>F) |
|-------|----:|----------:|----------:|
| group |   2 | 0.3571473 | 0.7036516 |
|       |  22 |        NA |        NA |

``` r
leveneTest(residuals(sapro_rich_glm) ~ sapro_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc")
```

|       |  Df |   F value |   Pr(\>F) |
|-------|----:|----------:|----------:|
| group |   2 | 0.7215314 | 0.4971548 |
|       |  22 |        NA |        NA |

Residuals/response distributions do not suggest the need for
transformation. Levene’s p \> 0.05 → fail to reject = variances can be
considered equal across groups.

Model results, group means, and post-hoc. Use Type II LR test of
variables due to unbalanced design.

``` r
Anova(sapro_rich_glm, type = 2, test.statistic = "LR")
```

    ## Analysis of Deviance Table (Type II tests)
    ## 
    ## Response: richness
    ##            LR Chisq Df Pr(>Chisq)  
    ## depth_csq    5.8012  1    0.01602 *
    ## field_type   5.9477  2    0.05111 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Differences in richness are very close to significance. Calculate
confidence intervals for figure. Estimated marginal means calculated in
this case

``` r
sapro_rich_em <- emmeans(sapro_rich_glm, ~ field_type, type = "response")
```

| field_type | response |       SE |  df | asymp.LCL | asymp.UCL |
|:-----------|---------:|---------:|----:|----------:|----------:|
| corn       | 100.7910 | 7.063812 | Inf |  87.85495 |  115.6317 |
| restored   | 120.5198 | 3.598141 | Inf | 113.66998 |  127.7825 |
| remnant    | 129.3314 | 8.806868 | Inf | 113.17258 |  147.7974 |

Confidence level used: 0.95

Model NS; no post hoc comparison…

``` r
sapro_rich_fig <- 
  ggplot(summary(sapro_rich_em), aes(x = field_type, y = response)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = response, ymax = asymp.UCL), width = 0, linewidth = lw) +
  labs(x = NULL, y = expression(atop("Richness", paste("(", italic(n), " OTUs)")))) +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1))
```

### Shannon’s diversity

Sequence depth square root transformed and centered

``` r
sapro_shan_lm <- lm(shannon ~ depth_csq + field_type, data = sapro_div)
```

Diagnostics

``` r
check_model(sapro_shan_lm)
```

![](resources/fungal_ecology_files/figure-gfm/sapro_shan_covar_diagnostics-1.png)<!-- -->

``` r
distribution_prob(sapro_shan_lm)
```

    ## 
    ## 
    ## Distribution    p_Residuals
    ## -------------  ------------
    ## normal              0.84375
    ## cauchy              0.12500
    ## gamma               0.03125
    ## 
    ## 
    ## Distribution    p_Response
    ## -------------  -----------
    ## chi                 0.2500
    ## gamma               0.2500
    ## weibull             0.1875

residuals distribution most likely normal, qq fit good, no evidence of
mean/variance increase response non-normal, check variance in groups
though

``` r
leveneTest(shannon ~ field_type, data = sapro_div) %>% as.data.frame() %>% kable(format = "pandoc")
```

|       |  Df |   F value |   Pr(\>F) |
|-------|----:|----------:|----------:|
| group |   2 | 0.3384785 | 0.7165013 |
|       |  22 |        NA |        NA |

``` r
leveneTest(residuals(sapro_shan_lm) ~ sapro_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc")
```

|       |  Df |   F value |   Pr(\>F) |
|-------|----:|----------:|----------:|
| group |   2 | 0.0967895 | 0.9081315 |
|       |  22 |        NA |        NA |

Residuals distribution does not suggest the need for transformation.
Levene’s p \> 0.05 → fail to reject = variances can be considered equal.

Model results, group means, and post-hoc

``` r
Anova(sapro_shan_lm, type = 2)
```

    ## Anova Table (Type II tests)
    ## 
    ## Response: shannon
    ##             Sum Sq Df F value Pr(>F)
    ## depth_csq    10.49  1  0.1744 0.6804
    ## field_type  114.23  2  0.9497 0.4029
    ## Residuals  1263.01 21

Sequence depth is not a significant predictor of Shannon diversity, nor
field type

``` r
sapro_shan_em <- emmeans(sapro_shan_lm, ~ field_type, type = "response")
```

Results tables below show the emmeans summary of group means and
confidence intervals, with sequencing depth as a covariate, and the post
hoc contrast of richness among field types.

``` r
sapro_shan_fig <- 
  ggplot(summary(sapro_shan_em), aes(x = field_type, y = emmean)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = emmean, ymax = upper.CL), width = 0, linewidth = lw) +
  labs(x = "Field type", y = expression(Shannon~diversity~paste("(", italic(e)^italic(H), ")"))) +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
```

## Abundance

Proportional biomass

``` r
sapro_ma_lm <- lm(sapro_mass ~ field_type, data = its_guild_ma)
par(mfrow = c(2,2))
plot(sapro_ma_lm) 
```

![](resources/fungal_ecology_files/figure-gfm/unnamed-chunk-153-1.png)<!-- -->

Variance looks consistent, no leverage points, poor qq fit

``` r
distribution_prob(sapro_ma_lm)
```

    ## 
    ## 
    ## Distribution    p_Residuals
    ## -------------  ------------
    ## normal              0.40625
    ## cauchy              0.15625
    ## gamma               0.15625
    ## 
    ## 
    ## Distribution    p_Response
    ## -------------  -----------
    ## gamma               0.6875
    ## exponential         0.1250
    ## binomial            0.0625

Residuals distribution fits normal, so do residuals

``` r
leveneTest(residuals(sapro_ma_lm) ~ sapro_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc") 
```

|       |  Df |   F value |   Pr(\>F) |
|-------|----:|----------:|----------:|
| group |   2 | 0.1060768 | 0.8998127 |
|       |  22 |        NA |        NA |

No covariate; response and residuals tests equivalent Residuals
distribution does not suggest the need for transformation. Levene’s p \>
0.05 → fail to reject = variances can be considered equal (aka
homoscedastic by group).

Produce model results, group means, and post-hoc, with arithmetic means
from emmeans

``` r
anova(sapro_ma_lm)
```

    ## Analysis of Variance Table
    ## 
    ## Response: sapro_mass
    ##            Df Sum Sq Mean Sq F value Pr(>F)
    ## field_type  2 0.1960 0.09799  0.3756 0.6912
    ## Residuals  22 5.7401 0.26091

``` r
sapro_ma_em <- emmeans(sapro_ma_lm, ~ field_type, type = "response")
```

| field_type |   emmean |        SE |  df |  lower.CL | upper.CL |
|:-----------|---------:|----------:|----:|----------:|---------:|
| corn       | 1.084440 | 0.2284349 |  22 | 0.6106951 | 1.558185 |
| restored   | 1.255165 | 0.1276990 |  22 | 0.9903330 | 1.519996 |
| remnant    | 1.055516 | 0.2553980 |  22 | 0.5258529 | 1.585179 |

Confidence level used: 0.95

``` r
sapro_ma_fig <- 
  ggplot(summary(sapro_ma_em), aes(x = field_type, y = emmean)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = emmean, ymax = upper.CL), width = 0, linewidth = lw) +
  labs(x = "Field type", y = expression(atop("Biomass (scaled)", paste(bold(`(`), "(", nmol[PLFA], " × ", g[soil]^{-1}, ")", " × ", paste("(rel. abund)", bold(`)`)))))) +
  # labs(x = "Field Type", y = "Biomass (scaled)") +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1.1))
```

## Beta Diversity

Community distance handled similarly to previous

### Biomass-weighted relative abundance

``` r
sapro_ma <- guildseq(its_avg_ma, its_meta, "saprotroph")
```

``` r
d_sapro_ma <- sapro_ma %>% 
  data.frame(row.names = 1) %>% 
  vegdist("bray")
mva_sapro_ma <- mva(d = d_sapro_ma, env = sites)
```

``` r
mva_sapro_ma$dispersion_test
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
    ## Groups     2 0.016606 0.0083028 1.4423   1999  0.265
    ## Residuals 22 0.126644 0.0057566                     
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##              corn  remnant restored
    ## corn              0.396500   0.6995
    ## remnant  0.392824            0.0880
    ## restored 0.702162 0.094644

``` r
mva_sapro_ma$permanova
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## adonis2(formula = d ~ dist_axis_1 + field_type, data = env, permutations = nperm, by = "terms")
    ##             Df SumOfSqs      R2      F Pr(>F)    
    ## dist_axis_1  1   0.5385 0.07454 2.0226 0.0045 ** 
    ## field_type   2   1.0951 0.15157 2.0564 0.0005 ***
    ## Residual    21   5.5916 0.77389                  
    ## Total       24   7.2253 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
mva_sapro_ma$pairwise_contrasts[c(1,3,2), c(1,2,4,3,8)] %>% 
  arrange(group1, desc(group2)) %>% 
  kable(format = "pandoc", caption = "Pairwise permanova contrasts")
```

| group1  | group2   | F_value |    R2 | p_value_adj |
|:--------|:---------|--------:|------:|------------:|
| corn    | restored |   3.042 | 0.131 |      0.0015 |
| corn    | remnant  |   1.766 | 0.195 |      0.0030 |
| remnant | restored |   1.178 | 0.059 |      0.2135 |

Pairwise permanova contrasts

Lingoes correction was not necessary. Based on the homogeneity of
variance test, the null hypothesis of equal variance among groups is
accepted across all clusters and in pairwise comparison of clusters
(both p\>0.05), supporting the application of a PERMANOVA test.

An effect of geographic distance (covariate) on pathogen communities was
detected. With geographic distance accounted for, the test variable
‘field type’ significantly explained variation in fungal communities,
with a post-hoc test revealing that communities in corn fields differed
from communities in restored and remnant fields.

Plotting results:

``` r
sapro_ma_ord_data <- mva_sapro_ma$ordination_scores %>% mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant")))
p_sapro_ma_centers <- sapro_ma_ord_data %>% 
  group_by(field_type) %>% 
  summarize(across(starts_with("Axis"), list(mean = mean, ci_l = ci_l, ci_u = ci_u), .names = "{.fn}_{.col}"), .groups = "drop") %>% 
  mutate(across(c(ci_l_Axis.1, ci_u_Axis.1), ~ mean_Axis.1 + .x),
         across(c(ci_l_Axis.2, ci_u_Axis.2), ~ mean_Axis.2 + .x),
         across(ends_with("Axis.1"), ~ .x * -1)) # reverse axis values to be consistent with other plots
sapro_ma_ord <- 
  ggplot(sapro_ma_ord_data, aes(x = Axis.1 * -1, y = Axis.2)) + # reverse axis
  geom_linerange(data = p_sapro_ma_centers, aes(x = mean_Axis.1, y = mean_Axis.2, xmin = ci_l_Axis.1, xmax = ci_u_Axis.1), linewidth = lw) +
  geom_linerange(data = p_sapro_ma_centers, aes(x = mean_Axis.1, y = mean_Axis.2, ymin = ci_l_Axis.2, ymax = ci_u_Axis.2), linewidth = lw) +
  geom_point(data = p_sapro_ma_centers, aes(x = mean_Axis.1, y = mean_Axis.2, fill = field_type), size = lg_size, stroke = lw, shape = 21) +
  geom_point(aes(fill = field_type), size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  scale_fill_manual(values = ft_pal) +
  labs(
    x = paste0("PCoA 1 (", mva_sapro_ma$axis_pct[1], "%)"),
    y = paste0("PCoA 2 (", mva_sapro_ma$axis_pct[2], "%)")) +
  theme_ord +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1))
```

#### Unified figure

``` r
fig5_ls <- (sapro_rich_fig / plot_spacer() / sapro_ma_fig) +
  plot_layout(heights = c(1,0.01,1)) 
fig5 <- (fig5_ls | plot_spacer() | sapro_ma_ord) +
  plot_layout(widths = c(0.35, 0.01, 0.64)) +
  plot_annotation(tag_levels = 'A') 
```

``` r
fig5
```

![](resources/fungal_ecology_files/figure-gfm/fig5-1.png)<!-- -->

**Fig 5.** Putative plant pathogen communities in **corn**,
**restored**, and **remnant** prairie fields. **a** OTU richness and
**b** biomass (nmol PLFA g soil^-1 \* (proportional sequence abundance))
are shown as columns with 95 % CIs. **c** Principal-coordinate (PCoA)
ordination of ITS-based (97 % OTU) community distances: small points =
sites, large circles = field-type centroids (error bars = 95 % CI).
Cornfields cluster apart from restored or remnant prairies (P \< 0.01).
Numbers in black circles give years since restoration. Axis labels show
the percent variation explained. Colours/shading: corn = grey, restored
= black, remnant = white.

### Sequence-based relative abundance

``` r
d_sapro <- sapro %>%
  data.frame(row.names = 1) %>%
  decostand("total") %>%
  vegdist("bray")
mva_sapro <- mva(d = d_sapro, env = sites)
```

``` r
mva_sapro$dispersion_test
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## Response: Distances
    ##           Df  Sum Sq   Mean Sq     F N.Perm Pr(>F)
    ## Groups     2 0.01522 0.0076101 1.229   1999 0.3265
    ## Residuals 22 0.13623 0.0061922                    
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##             corn remnant restored
    ## corn             0.35300   0.9465
    ## remnant  0.34863           0.1100
    ## restored 0.94409 0.10783

``` r
mva_sapro$permanova
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## adonis2(formula = d ~ dist_axis_1 + field_type, data = env, permutations = nperm, by = "terms")
    ##             Df SumOfSqs      R2      F Pr(>F)    
    ## dist_axis_1  1   0.5499 0.07873 2.1652  3e-03 ** 
    ## field_type   2   1.1017 0.15773 2.1691  5e-04 ***
    ## Residual    21   5.3330 0.76354                  
    ## Total       24   6.9845 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
mva_sapro$pairwise_contrasts[c(1,3,2), c(1,2,4,3,8)] %>%
  arrange(group1, desc(group2)) %>% 
  kable(format = "pandoc", caption = "Pairwise permanova contrasts")
```

| group1  | group2   | F_value |    R2 | p_value_adj |
|:--------|:---------|--------:|------:|------------:|
| corn    | restored |   3.196 | 0.136 |      0.0015 |
| corn    | remnant  |   1.947 | 0.208 |      0.0037 |
| remnant | restored |   1.231 | 0.061 |      0.1845 |

Pairwise permanova contrasts

Lingoes correction was not necessary. Based on the homogeneity of
variance test, the null hypothesis of equal variance among groups is
accepted across all clusters and in pairwise comparison of clusters
(both p\>0.05), supporting the application of a PERMANOVA test.

An effect of geographic distance (covariate) on pathogen communities was
detected With geographic distance accounted for, the test variable
‘field type’ significantly explained variation in fungal communities,
with a post-hoc test revealing that communities in corn fields differed
from communities in restored and remnant fields.

Plotting results:

``` r
sapro_ord_data <- mva_sapro$ordination_scores %>% mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant")))
p_sapro_centers <- sapro_ord_data %>%
  group_by(field_type) %>%
  summarize(across(starts_with("Axis"), list(mean = mean, ci_l = ci_l, ci_u = ci_u), .names = "{.fn}_{.col}"), .groups = "drop") %>%
  mutate(across(c(ci_l_Axis.1, ci_u_Axis.1), ~ mean_Axis.1 + .x),
         across(c(ci_l_Axis.2, ci_u_Axis.2), ~ mean_Axis.2 + .x))
sapro_ord <-
  ggplot(sapro_ord_data, aes(x = Axis.1, y = Axis.2)) +
  geom_linerange(data = p_sapro_centers, aes(x = mean_Axis.1, y = mean_Axis.2, xmin = ci_l_Axis.1, xmax = ci_u_Axis.1), linewidth = lw) +
  geom_linerange(data = p_sapro_centers, aes(x = mean_Axis.1, y = mean_Axis.2, ymin = ci_l_Axis.2, ymax = ci_u_Axis.2), linewidth = lw) +
  geom_point(data = p_sapro_centers, aes(x = mean_Axis.1, y = mean_Axis.2, fill = field_type), size = lg_size, stroke = lw, shape = 21) +
  geom_point(aes(fill = field_type), size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  scale_fill_manual(values = ft_pal) +
  labs(
    x = paste0("PCoA 1 (", mva_sapro$axis_pct[1], "%)"),
    y = paste0("PCoA 2 (", mva_sapro$axis_pct[2], "%)")) +
  theme_ord +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1))
```

#### Supplemental figure

``` r
sapro_shan_ord_sup <- (sapro_shan_fig | plot_spacer() | sapro_ord) +
  plot_layout(widths = c(0.45, 0.01, 0.55)) +
  plot_annotation(tag_levels = 'A') 
```

``` r
sapro_shan_ord_sup
```

![](resources/fungal_ecology_files/figure-gfm/sapro_shan_ord_sup-1.png)<!-- -->

#### Contrast community metrics

Procrustes test on PCoA values using axes with eigenvalues exceeding a
broken stick model

``` r
set.seed(20251119)
sapro_protest <- protest(
  pcoa(d_sapro)$vectors[, 2],
  pcoa(d_sapro_ma)$vectors[, 2],
  permutations = 1999
)
sapro_protest
```

    ## 
    ## Call:
    ## protest(X = pcoa(d_sapro)$vectors[, 2], Y = pcoa(d_sapro_ma)$vectors[,      2], permutations = 1999) 
    ## 
    ## Procrustes Sum of Squares (m12 squared):        0.1409 
    ## Correlation in a symmetric Procrustes rotation: 0.9269 
    ## Significance:  5e-04 
    ## 
    ## Permutation: free
    ## Number of permutations: 1999

Including biomass changes little. The spatial configurations of both
ordinations are highly correlated. $R^{2}=$ 0.86, p\<0.001.

## Saprotroph Indicator Species

``` r
sapro_ind <- inspan(its_avg_ma, its_meta, "saprotroph", sites)
sapro_ind %>% 
  select(A, B, stat, p_val_adj, field_type, species, starts_with("corn"), starts_with("restor"), starts_with("rem")) %>% 
  filter(species != "unidentified") %>% 
  arrange(field_type, p_val_adj) %>% 
  mutate(across(A:p_val_adj, ~ round(.x, 3)),
         across(corn_avg:remnant_ci, ~ num(.x, notation = "sci"))) %>% 
  kable(format = "pandoc", caption = "Indicator species analysis results with biomass-aware relative abundances in field types")
```

| A | B | stat | p_val_adj | field_type | species | corn_avg | corn_ci | restored_avg | restored_ci | remnant_avg | remnant_ci |
|---:|---:|---:|---:|:---|:---|---:|---:|---:|---:|---:|---:|
| 1.000 | 1.00 | 1.000 | 0.180 | corn | Conocybe_apala | 1.22e-2 | 2.24e-2 | 0 | 0 | 0 | 0 |
| 0.995 | 1.00 | 0.998 | 0.240 | corn | Bolbitius_titubans | 2.91e-2 | 3.36e-2 | 1.36e-4 | 1.66e-4 | 0 | 0 |
| 0.930 | 0.80 | 0.863 | 0.450 | corn | Stachybotrys_limonispora | 2.31e-3 | 1.68e-3 | 1.73e-4 | 1.82e-4 | 0 | 0 |
| 0.707 | 1.00 | 0.841 | 0.450 | corn | Humicola_grisea | 6.37e-2 | 3.46e-2 | 2.04e-2 | 1.25e-2 | 5.96e-3 | 4.55e-3 |
| 1.000 | 0.60 | 0.775 | 0.450 | corn | Murispora_galii | 5.38e-3 | 7.81e-3 | 0 | 0 | 0 | 0 |
| 0.957 | 1.00 | 0.978 | 0.498 | corn | Tausonia_pullulans | 1.30e-1 | 1.99e-1 | 2.43e-3 | 3.18e-3 | 3.36e-3 | 2.59e-3 |
| 0.983 | 0.60 | 0.768 | 0.650 | corn | Cyphellophora_suttonii | 7.68e-3 | 1.18e-2 | 1.34e-4 | 2.03e-4 | 0 | 0 |
| 0.932 | 0.60 | 0.748 | 0.650 | corn | Cheilymenia_stercorea | 1.21e-3 | 1.52e-3 | 8.80e-5 | 1.42e-4 | 0 | 0 |
| 0.801 | 0.60 | 0.693 | 0.650 | corn | Preussia_terricola | 3.85e-3 | 5.24e-3 | 9.56e-4 | 1.87e-3 | 0 | 0 |
| 0.770 | 0.60 | 0.680 | 0.650 | corn | Schizothecium_inaequale | 4.28e-4 | 4.23e-4 | 1.28e-4 | 2.50e-4 | 0 | 0 |
| 0.903 | 0.75 | 0.823 | 0.650 | remnant | Coniochaeta_decumbens | 2.46e-5 | 4.82e-5 | 5.01e-4 | 3.97e-4 | 4.88e-3 | 5.44e-3 |
| 0.628 | 1.00 | 0.793 | 0.650 | remnant | Gliomastix_roseogrisea | 0 | 0 | 2.11e-3 | 1.54e-3 | 3.56e-3 | 3.10e-3 |
| 1.000 | 0.50 | 0.707 | 0.650 | remnant | Coprinellus_xanthothrix | 0 | 0 | 0 | 0 | 9.52e-2 | 1.86e-1 |
| 1.000 | 0.50 | 0.707 | 0.650 | remnant | Savoryella_paucispora | 0 | 0 | 0 | 0 | 3.00e-4 | 4.98e-4 |
| 0.851 | 0.50 | 0.652 | 0.650 | remnant | Minimelanolocus_asiaticus | 0 | 0 | 2.71e-4 | 4.08e-4 | 1.55e-3 | 2.35e-3 |

Indicator species analysis results with biomass-aware relative
abundances in field types

## Constrained Analysis

Env covars processed in the ITS section (see above)

``` r
spe_sapro_wi_resto <- sapro_ma %>%
  filter(field_name %in% rownames(env_expl)) %>%
  column_to_rownames(var = "field_name")

sapro_mod_null <- dbrda(spe_sapro_wi_resto ~ 1 + Condition(env_cov), data = env_expl, distance = "bray")
sapro_mod_full <- dbrda(spe_sapro_wi_resto ~ . + Condition(env_cov), data = env_expl, distance = "bray")
sapro_mod_step <- ordistep(sapro_mod_null,
                           scope = formula(sapro_mod_full),
                           direction = "forward",
                           permutations = 1999,
                           trace = FALSE)
```

### Constrained Analysis Results

``` r
sapro_mod_step
```

    ## 
    ## Call: dbrda(formula = spe_sapro_wi_resto ~ Condition(env_cov) + OM + gf_index + pl_rich + NO3, data = env_expl, distance = "bray")
    ## 
    ##               Inertia Proportion Rank
    ## Total          3.6087     1.0000     
    ## Conditional    0.3655     0.1013    1
    ## Constrained    1.7265     0.4784    4
    ## Unconstrained  1.5166     0.4203    7
    ## 
    ## Inertia is squared Bray distance
    ## 
    ## Eigenvalues for constrained axes:
    ## dbRDA1 dbRDA2 dbRDA3 dbRDA4 
    ## 0.6161 0.4542 0.3871 0.2690 
    ## 
    ## Eigenvalues for unconstrained axes:
    ##   MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7 
    ## 0.3541 0.3072 0.2463 0.2111 0.1730 0.1461 0.0789

``` r
(sapro_mod_r2   <- RsquareAdj(sapro_mod_step, permutations = 1999))
```

    ## $r.squared
    ## [1] 0.4784291
    ## 
    ## $adj.r.squared
    ## [1] 0.2599313

``` r
(sapro_mod_glax <- anova(sapro_mod_step, permutations = 1999))
```

    ## Permutation test for dbrda under reduced model
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## Model: dbrda(formula = spe_sapro_wi_resto ~ Condition(env_cov) + OM + gf_index + pl_rich + NO3, data = env_expl, distance = "bray")
    ##          Df SumOfSqs      F Pr(>F)    
    ## Model     4   1.7265 1.9921  5e-04 ***
    ## Residual  7   1.5166                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(sapro_mod_inax <- anova(sapro_mod_step, by = "axis", permutations = 1999))
```

    ## Permutation test for dbrda under reduced model
    ## Forward tests for axes
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## Model: dbrda(formula = spe_sapro_wi_resto ~ Condition(env_cov) + OM + gf_index + pl_rich + NO3, data = env_expl, distance = "bray")
    ##          Df SumOfSqs      F Pr(>F)    
    ## dbRDA1    1  0.61611 2.8436 0.0005 ***
    ## dbRDA2    1  0.45423 2.3960 0.0015 ** 
    ## dbRDA3    1  0.38713 2.2973 0.0015 ** 
    ## dbRDA4    1  0.26902 1.7738 0.0175 *  
    ## Residual  7  1.51664                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(sapro_mod_axpct <- round(100 * sapro_mod_step$CCA$eig / sum(sapro_mod_step$CCA$eig), 1))
```

    ## dbRDA1 dbRDA2 dbRDA3 dbRDA4 
    ##   35.7   26.3   22.4   15.6

``` r
sapro_mod_step$anova %>% 
  as.data.frame() %>% 
  mutate(p.adj = p.adjust(`Pr(>F)`, "fdr")) %>% 
  kable(, format = "pandoc")
```

|             |  Df |      AIC |        F | Pr(\>F) |     p.adj |
|-------------|----:|---------:|---------:|--------:|----------:|
| \+ OM       |   1 | 18.11061 | 1.792868 |   0.011 | 0.0146667 |
| \+ gf_index |   1 | 17.66275 | 1.864729 |   0.008 | 0.0146667 |
| \+ pl_rich  |   1 | 16.98962 | 1.826341 |   0.009 | 0.0146667 |
| \+ NO3      |   1 | 16.37391 | 1.560162 |   0.041 | 0.0410000 |

Based on permutation tests with n=1999 permutations, after accounting
for inter-site pairwise distance as a covariate, the model shows
correlations between the site ordination on saprotroph communities and
the selected explanatory variables (p\<0.001). The first four
constrained axes are also significant (p\<0.05). The selected variables
explain $R^{2}_{\text{Adj}}$= 26% of the community variation. Selected
explanatory variables are SOM, grass-forb index, plant richness, and
nitrate; see table for individual p values and statistics.

#### Saprotroph constrained figure

``` r
sapro_mod_step_eig <- round(sapro_mod_step$CCA$eig * 100, 1)
sapro_mod_scor <- scores(
  sapro_mod_step,
  choices = c(1, 2),
  display = c("bp", "sites"),
  tidy = FALSE
)
sapro_mod_scor_site <- sapro_mod_scor$sites %>%
  data.frame() %>%
  rownames_to_column(var = "field_name") %>%
  left_join(sites, by = join_by(field_name))
sapro_mod_scor_bp <- bind_rows(
  sapro_mod_scor$biplot %>%
    data.frame() %>%
    rownames_to_column(var = "envvar") %>%
    mutate(envlabs = c("'OM'", "'>forb'", "plant~spp.", 'NO[3]^"–"')),
  data.frame(
    envvar = "gf_index",
    dbRDA1 = -sapro_mod_scor$biplot["gf_index", 1],
    dbRDA2 = -sapro_mod_scor$biplot["gf_index", 2],
    envlabs = "'>grass'")
) %>% 
  arrange(envvar, envlabs) %>% 
  mutate(
    origin = 0,
    m = dbRDA2 / dbRDA1,
    d = sqrt(dbRDA1^2 + dbRDA2^2),
    dadd = sqrt((max(dbRDA1)-min(dbRDA2))^2 + (max(dbRDA2)-min(dbRDA2))^2)*dadd_adj,
    labx = ((d+dadd)*cos(atan(m)))*(dbRDA1/abs(dbRDA1)),
    laby = ((d+dadd)*sin(atan(m)))*(dbRDA1/abs(dbRDA1)))
```

## Saprotroph and plant community correlations

Data for these tests

``` r
sapro_resto <- its_guild %>% 
  filter(field_type != "corn", region != "FL") %>% 
  left_join(its_guild_ma %>% select(field_name, sapro_mass), by = join_by(field_name)) %>% 
  left_join(prich %>% select(field_name, pl_rich, pl_shan), by = join_by(field_name)) %>% 
  mutate(
    sapro_prop = sapro_abund / fungi_abund, # no zeroes present...
    notsapro_abund = fungi_abund - sapro_abund,
    fungi_mass_lc = as.numeric(scale(log(fungi_mass), center = TRUE, scale = FALSE))
  ) %>% 
  select(-patho_abund, -c(annual:shrubTree))
```

### Plant richness and saprotroph biomass

Is plant richness related to saprotroph mass?

``` r
saprofa_prich_lm <- lm(sapro_mass ~ pl_rich, data = sapro_resto)
summary(saprofa_prich_lm)
```

    ## 
    ## Call:
    ## lm(formula = sapro_mass ~ pl_rich, data = sapro_resto)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.56697 -0.32365  0.06416  0.24916  0.53019 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  2.422247   0.349522   6.930 2.49e-05 ***
    ## pl_rich     -0.030600   0.008485  -3.607  0.00412 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3617 on 11 degrees of freedom
    ## Multiple R-squared:  0.5418, Adjusted R-squared:  0.5001 
    ## F-statistic: 13.01 on 1 and 11 DF,  p-value: 0.004123

Strong negative relationship worthy of closer examination. It isn’t
driven by grass-forb index (see below), so this isn’t a confounding with
C4 grass abundance… KORP site appears to have some leverage (not shown),
conduct logistic model and diagnostics as before.

Multiple, weighted logistic glm

``` r
sapro_prich_glm <- glm(sapro_prop ~ fungi_mass_lc + pl_rich,
                       data = sapro_resto, family = quasibinomial(link = "logit"),
                       weights = fungi_abund)
summary(sapro_prich_glm) # dispersion parameter >62 justifies quasibinomial
```

    ## 
    ## Call:
    ## glm(formula = sapro_prop ~ fungi_mass_lc + pl_rich, family = quasibinomial(link = "logit"), 
    ##     data = sapro_resto, weights = fungi_abund)
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)   -0.517543   0.218103  -2.373   0.0391 *
    ## fungi_mass_lc -0.216797   0.186892  -1.160   0.2730  
    ## pl_rich       -0.015656   0.005417  -2.890   0.0161 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for quasibinomial family taken to be 62.46492)
    ## 
    ##     Null deviance: 1187.92  on 12  degrees of freedom
    ## Residual deviance:  663.24  on 10  degrees of freedom
    ## AIC: NA
    ## 
    ## Number of Fisher Scoring iterations: 4

Diagnostics

``` r
check_model(sapro_prich_glm)
```

![](resources/fungal_ecology_files/figure-gfm/unnamed-chunk-167-1.png)<!-- -->

``` r
check_collinearity(sapro_prich_glm)
```

    ## # Check for Multicollinearity
    ## 
    ## Low Correlation
    ## 
    ##           Term  VIF   VIF 95% CI adj. VIF Tolerance Tolerance 95% CI
    ##  fungi_mass_lc 1.28 [1.04, 2.87]     1.13      0.78     [0.35, 0.96]
    ##        pl_rich 1.28 [1.04, 2.87]     1.13      0.78     [0.35, 0.96]

``` r
augment(sapro_prich_glm)
```

    ## # A tibble: 13 × 10
    ##    sapro_prop fungi_mass_lc pl_rich `(weights)` .fitted   .resid   .hat .sigma     .cooksd .std.resid
    ##         <dbl>         <dbl>   <int>       <dbl>   <dbl>    <dbl>  <dbl>  <dbl>       <dbl>      <dbl>
    ##  1      0.217       0.470        31       6936.  -1.10   -6.15   0.199    8.27 0.0608        -0.870  
    ##  2      0.290       0.142        37       7578.  -1.13    9.08   0.0864   7.98 0.0473         1.20   
    ##  3      0.130       0.110        46       7300.  -1.26  -19.8    0.119    4.93 0.284         -2.67   
    ##  4      0.216       0.0610       46       8145.  -1.25   -1.51   0.119    8.57 0.00184       -0.203  
    ##  5      0.303       0.468        13       7226.  -0.823  -0.504  0.476    8.58 0.00234       -0.0880 
    ##  6      0.285       0.320        33       8810.  -1.10    7.78   0.165    8.10 0.0787         1.08   
    ##  7      0.275      -0.540        36       8343.  -0.964  -0.257  0.392    8.58 0.000374      -0.0417 
    ##  8      0.280      -0.511        37       8853.  -0.986   1.84   0.363    8.55 0.0163         0.292  
    ##  9      0.206      -0.439        57       7392.  -1.31   -1.27   0.231    8.57 0.00332       -0.183  
    ## 10      0.213      -0.339        55      10437.  -1.31    0.0161 0.264    8.58 0.000000675    0.00238
    ## 11      0.229       0.212        53       9062.  -1.39    6.98   0.314    8.11 0.180          1.07   
    ## 12      0.255       0.00271      42       8537.  -1.18    4.12   0.0859   8.46 0.00948        0.545  
    ## 13      0.259       0.0438       27       8324.  -0.950  -4.14   0.186    8.45 0.0253        -0.581

Long tails and low n showing moderate structure.

``` r
distribution_prob(sapro_prich_glm)
```

    ## 
    ## 
    ## Distribution    p_Residuals
    ## -------------  ------------
    ## cauchy              0.65625
    ## normal              0.31250
    ## pareto              0.03125
    ## 
    ## 
    ## Distribution    p_Response
    ## -------------  -----------
    ## weibull            0.31250
    ## beta               0.28125
    ## uniform            0.18750

Residuals distribution long-tailed / normal

``` r
loocv_saglm_gfi <- map_dbl(seq_len(nrow(sapro_resto)), function(i){
  exp(coef(glm(sapro_prop ~ fungi_mass_lc + pl_rich,
               data = sapro_resto[-i, ],
               family = quasibinomial(link = "logit"),
               weights = fungi_abund))["pl_rich"])
})
summary(loocv_saglm_gfi)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.9813  0.9842  0.9844  0.9844  0.9847  0.9875

``` r
(cv_saglm <- (sd(loocv_saglm_gfi) / mean(loocv_saglm_gfi) * 100) %>% round(., 1))
```

    ## [1] 0.1

Grass-forb index LOOCV variation of 0.1% on the back-transformed scale
shows remarkably consistent prediction.

``` r
loocv_saglm_fma <- map_dbl(seq_len(nrow(sapro_resto)), function(i){
  exp(coef(glm(sapro_prop ~ fungi_mass_lc + pl_rich,
               data = sapro_resto[-i, ],
               family = quasibinomial(link = "logit"),
               weights = fungi_abund))["fungi_mass_lc"])
})
summary(loocv_saglm_fma)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.7349  0.7953  0.8007  0.8051  0.8099  0.8803

``` r
(cv_saglm <- (sd(loocv_saglm_fma) / mean(loocv_saglm_fma) * 100) %>% round(., 1))
```

    ## [1] 4.6

Slope CV on fungal mass of 4.6% is low but shows greater noise with the
covariate than the test variable.

``` r
saglm_crpldata <- as.data.frame(crPlots(sapro_prich_glm))
```

![](resources/fungal_ecology_files/figure-gfm/unnamed-chunk-171-1.png)<!-- -->

Noise in fungal mass data is obvious here. Fit of partial gf_index is
clean.

Partial R2 values

``` r
data.frame(
  term = c("fungal_mass", "pl_rich"),
  partial_R2 = rsq.partial(sapro_prich_glm, adj = TRUE)$partial.rsq
) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% 
  kable(format = "pandoc", caption = "Partial R2 from weighted logistic regression")
```

| term        | partial_R2 |
|:------------|-----------:|
| fungal_mass |      0.123 |
| pl_rich     |      0.455 |

Partial R2 from weighted logistic regression

Model summary

``` r
sapro_null_glm <- glm(sapro_prop ~ 1,
                      data = sapro_resto, family = quasibinomial(link = "logit"),
                      weights = fungi_abund)
anova(sapro_null_glm, sapro_prich_glm, test = "F")
```

    ## Analysis of Deviance Table
    ## 
    ## Model 1: sapro_prop ~ 1
    ## Model 2: sapro_prop ~ fungi_mass_lc + pl_rich
    ##   Resid. Df Resid. Dev Df Deviance      F  Pr(>F)  
    ## 1        12    1187.92                             
    ## 2        10     663.24  2   524.68 4.1998 0.04742 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Deviance explained

``` r
(sapro_prich_glm_pr2 <- round(1-(summary(sapro_prich_glm)$deviance / summary(sapro_prich_glm)$null.deviance), 3))
```

    ## [1] 0.442

Summary of terms Odds ratio prediction and confidence intervals on the
prediction scale, results on the increment of an increase of 10 plant
species desired due to scale of that variable. Note: in the following
output, percent predicted changes are calculated *in excess* of 100%
(e.g., 0.085 = -15%).

``` r
tidy(sapro_prich_glm) %>% 
  mutate(odds_ratio = exp(estimate), exp_std.error = exp(std.error),
         across(where(is.numeric), ~ round(.x, 3))) %>% 
  select(term, estimate, odds_ratio, std.error, exp_std.error, statistic, p.value) %>% 
  kable(format = "pandoc", caption = "Summary of terms from weighted logistic regression")
```

| term          | estimate | odds_ratio | std.error | exp_std.error | statistic | p.value |
|:--------------|---------:|-----------:|----------:|--------------:|----------:|--------:|
| (Intercept)   |   -0.518 |      0.596 |     0.218 |         1.244 |    -2.373 |   0.039 |
| fungi_mass_lc |   -0.217 |      0.805 |     0.187 |         1.205 |    -1.160 |   0.273 |
| pl_rich       |   -0.016 |      0.984 |     0.005 |         1.005 |    -2.890 |   0.016 |

Summary of terms from weighted logistic regression

``` r
(sapro_or_pct <- round(1-(exp(coef(sapro_prich_glm)[3])^10), 3))
```

    ## pl_rich 
    ##   0.145

``` r
(1-(exp(confint(sapro_prich_glm))^10)) %>%
  as.data.frame() %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% 
  kable(format = "pandoc", caption = "95% confidence intervals, back transformed from the log-scale")
```

|               | 2.5 % | 97.5 % |
|---------------|------:|-------:|
| (Intercept)   | 1.000 |  0.603 |
| fungi_mass_lc | 0.997 | -3.547 |
| pl_rich       | 0.231 |  0.049 |

95% confidence intervals, back transformed from the log-scale

Create objects for plotting

``` r
saglm_med_fungi <- median(sapro_resto$fungi_mass_lc, na.rm = TRUE)
saglm_med_abund <- median(sapro_resto$fungi_abund, na.rm = TRUE) # Needed for weight context
saglm_newdat <- tibble(
  pl_rich = seq(min(sapro_resto$pl_rich, na.rm = TRUE),
                max(sapro_resto$pl_rich, na.rm = TRUE),
                length.out = 200),
  fungi_mass_lc = saglm_med_fungi,
  fungi_abund = saglm_med_abund
)
```

Predict on link scale, back-transform with plogis

``` r
saglm_pred <- predict(sapro_prich_glm, newdata = saglm_newdat, type = "link", se.fit = TRUE) %>%
  as_tibble() %>%
  bind_cols(saglm_newdat) %>%
  mutate(
    fit_prob = plogis(fit),
    lwr_prob = plogis(fit - 1.96 * se.fit),
    upr_prob = plogis(fit + 1.96 * se.fit)
  )
```

``` r
fig8a <-
  ggplot(saglm_pred, aes(x = pl_rich, y = fit_prob)) +
  # geom_ribbon(aes(ymin = lwr_med, ymax = upr_med), fill = "gray90") +
  geom_line(color = "black", linewidth = lw) +
  geom_point(data = sapro_resto, aes(x = pl_rich, y = sapro_prop, fill = field_type),
             size = sm_size, stroke = lw, shape = 21) +
  geom_text(data = sapro_resto, aes(x = pl_rich, y = sapro_prop, label = yr_since),
            size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  labs(
    x = expression(paste("Plant richness (", italic(n), " species)")),
    y = "Saprotroph proportion (share of total sequences)",
    tag = "A"
  ) +
  scale_fill_manual(name = "Field type", values = ft_pal[2:3]) +
  theme_cor +
  theme(legend.position = c(0.03, 0.25),
        legend.justification = c(0, 1),
        legend.title = element_text(size = 9, face = 1),
        legend.text = element_text(size = 8, face = 1),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2),
        legend.key = element_rect(fill = "white"),
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1))
```

``` r
fig8b <-
  cbind(saglm_crpldata, sapro_resto %>% select(field_type, yr_since)) %>%
  ggplot(aes(x = fungi_mass_lc.fungi_mass_lc, y = fungi_mass_lc.sapro_prop)) +
  geom_smooth(method = "lm", color = "black", linewidth = lw, se = FALSE) +
  geom_point(aes(fill = field_type), size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  labs(x = expression("Fungal biomass"~paste("(", nmol[PLFA], " × ", g[soil]^{-1}, ")")), 
       y = NULL, 
       tag = "B") +
  scale_fill_manual(values = ft_pal[2:3]) +
  scale_y_continuous(breaks = c(-0.5, 0, 0.5)) +
  scale_x_continuous(breaks = log(c(2.7, 3.8, 5.3, 7.5)) - mean(log(sapro_resto$fungi_mass)),
  labels = c(2.7, 3.8, 5.3, 7.5)) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0.025, 1))
```

``` r
fig8c <-
  cbind(saglm_crpldata, sapro_resto %>% select(field_type, yr_since)) %>%
  ggplot(aes(x = pl_rich.pl_rich, y = pl_rich.sapro_prop)) +
  geom_smooth(method = "lm", color = "black", linewidth = lw, se = FALSE) +
  geom_point(aes(fill = field_type), size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  labs(x = expression(paste("Plant richness (", italic(n), " species)")), 
       y = NULL, tag = "C") +
  scale_fill_manual(values = ft_pal[2:3]) +
  scale_y_continuous(breaks = c(-1, 0, 1)) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0.025, 1))
fig8yax_grob <- textGrob(
  "Saprotroph proportion (partial residuals, logit scale)",
  x = 0.5,
  y = 0.5,
  hjust = 0.5,
  vjust = 0.5,
  rot = 90,
  gp = gpar(cex = 9/12)
)
```

``` r
fig8rh <- (fig8b / plot_spacer() / fig8c) +
  plot_layout(heights = c(0.5, 0.01,0.5))
fig8 <- (fig8a | plot_spacer() | (wrap_elements(full = fig8yax_grob) & theme(plot.tag = element_blank())) | fig8rh) +
  plot_layout(widths = c(0.62, 0.004, 0.02, 0.38))
```

``` r
fig8
```

![](resources/fungal_ecology_files/figure-gfm/fig8-1.png)<!-- -->

### Plant diversity and saprotroph biomass

Is plant diversity related to saprotroph mass?

``` r
saprofa_pshan_lm <- lm(sapro_mass ~ pl_shan, data = sapro_resto)
summary(saprofa_pshan_lm)
```

    ## 
    ## Call:
    ## lm(formula = sapro_mass ~ pl_shan, data = sapro_resto)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.65416 -0.32568  0.03834  0.32797  0.71819 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  2.43532    0.52764   4.616 0.000746 ***
    ## pl_shan     -0.09224    0.03882  -2.376 0.036753 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4344 on 11 degrees of freedom
    ## Multiple R-squared:  0.3392, Adjusted R-squared:  0.2791 
    ## F-statistic: 5.646 on 1 and 11 DF,  p-value: 0.03675

Saprotroph mass and plant diversity may also deserve further
examination.

``` r
sapro_pshan_glm <- glm(sapro_prop ~ fungi_mass_lc + pl_shan,
                       data = sapro_resto, family = quasibinomial(link = "logit"),
                       weights = fungi_abund)
summary(sapro_pshan_glm)
```

    ## 
    ## Call:
    ## glm(formula = sapro_prop ~ fungi_mass_lc + pl_shan, family = quasibinomial(link = "logit"), 
    ##     data = sapro_resto, weights = fungi_abund)
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)   -0.75763    0.34121  -2.220   0.0507 .
    ## fungi_mass_lc -0.07564    0.23605  -0.320   0.7552  
    ## pl_shan       -0.02839    0.02533  -1.121   0.2886  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for quasibinomial family taken to be 101.3357)
    ## 
    ##     Null deviance: 1187.9  on 12  degrees of freedom
    ## Residual deviance: 1059.9  on 10  degrees of freedom
    ## AIC: NA
    ## 
    ## Number of Fisher Scoring iterations: 4

NS

### Saprotroph biomass and grass/forb composition

``` r
sama_rest_m <- lm(sapro_mass ~ gf_index, data = sapro_resto)
summary(sama_rest_m)
```

    ## 
    ## Call:
    ## lm(formula = sapro_mass ~ gf_index, data = sapro_resto)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.8700 -0.3270  0.1771  0.3471  0.6007 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   1.2147     0.1336    9.09  1.9e-06 ***
    ## gf_index     -0.7966     0.5010   -1.59     0.14    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4818 on 11 degrees of freedom
    ## Multiple R-squared:  0.1869, Adjusted R-squared:  0.113 
    ## F-statistic: 2.528 on 1 and 11 DF,  p-value: 0.1401

Likely not enough of a relationship to warrant further attention Examine
mass+richness models, similar to pathogen models

``` r
sapro_gf_glm <- glm(sapro_prop ~ fungi_mass_lc + gf_index,
                    data = sapro_resto, family = quasibinomial(link = "logit"),
                    weights = fungi_abund)
summary(sapro_gf_glm)
```

    ## 
    ## Call:
    ## glm(formula = sapro_prop ~ fungi_mass_lc + gf_index, family = quasibinomial(link = "logit"), 
    ##     data = sapro_resto, weights = fungi_abund)
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   -1.13286    0.07532 -15.041 3.41e-08 ***
    ## fungi_mass_lc  0.01934    0.23267   0.083    0.935    
    ## gf_index      -0.07304    0.29889  -0.244    0.812    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for quasibinomial family taken to be 111.5614)
    ## 
    ##     Null deviance: 1187.9  on 12  degrees of freedom
    ## Residual deviance: 1178.6  on 10  degrees of freedom
    ## AIC: NA
    ## 
    ## Number of Fisher Scoring iterations: 4

NS

### Plant species and saprotroph abundance

Field age and restored vs. remnant status don’t appear to relate with
saprotroph abundance. Do any particular plant species associate with
high or low saprotroph sites? Try an indicator species analysis to find
out.

``` r
sapro_groups <- setNames(ifelse(sapro_resto$sapro_prop > median(sapro_resto$sapro_prop), "hsap", "lsap"), 
                         sapro_resto$field_name)
plantsap_ind <- multipatt(
  x = plant %>% filter(SITE %in% names(sapro_groups)) %>% 
    column_to_rownames("SITE") %>% select(-BARESOIL, -LITTER, -ROSA, -SALIX),
  cluster = sapro_groups, 
  func = "IndVal.g", max.order = 1, control = how(nperm = 1999)
)
sapro_sig_codes <- plantsap_ind$sign %>%
  filter(p.value < 0.1) %>%
  rownames_to_column("CODE") %>%
  mutate(group = c("hsap", "lsap")[index]) %>%
  select(CODE, group, stat, p.value)
```

``` r
sapro_plant_summary <- plant %>%
  pivot_longer(-SITE, names_to = "CODE", values_to = "pct") %>%
  filter(SITE %in% names(sapro_groups), CODE %in% sapro_sig_codes$CODE) %>%
  mutate(grp = sapro_groups[SITE]) %>%
  group_by(grp, CODE) %>%
  summarize(mean_pct = mean(pct), .groups = "drop") %>%
  pivot_wider(names_from = grp, values_from = mean_pct, values_fill = 0) %>%
  left_join(sapro_sig_codes, by = "CODE") %>%
  left_join(read_csv(root_path("clean_data/plant_meta.csv"), show_col_types = FALSE), by = "CODE") %>%
  select(species = NAME, family = FAMILY, life_hist = LIFEHIST, plcvr_hsap = hsap, plcvr_lsap = lsap, stat, p.value) %>%
  mutate(p.adj = p.adjust(p.value, "fdr"), across(where(is.numeric), ~ round(.x, 3))) %>% 
  arrange(desc(plcvr_hsap)) # Sort by high-sapro indicators first
kable(sapro_plant_summary, format = "pandoc", caption = "Plant indicators in low or high saprotroph sites")
```

| species | family | life_hist | plcvr_hsap | plcvr_lsap | stat | p.value | p.adj |
|:---|:---|:---|---:|---:|---:|---:|---:|
| Solidago canadensis | Asteraceae | perennial | 6.483 | 2.657 | 0.842 | 0.028 | 0.078 |
| Artemisia ludoviciana | Asteraceae | perennial | 6.117 | 0.000 | 0.707 | 0.076 | 0.078 |
| Echinacea purpurea | Asteraceae | perennial | 3.750 | 0.000 | 0.707 | 0.064 | 0.078 |
| Symphyotrichum novae-angliae | Asteraceae | perennial | 3.083 | 0.729 | 0.821 | 0.061 | 0.078 |
| Heliopsis helianthoides | Asteraceae | perennial | 2.867 | 0.171 | 0.971 | 0.003 | 0.033 |
| Symphyotrichum laevis | Asteraceae | perennial | 2.100 | 0.200 | 0.872 | 0.060 | 0.078 |
| Cirsium arvense | Asteraceae | perennial | 2.067 | 0.000 | 0.707 | 0.068 | 0.078 |
| Elymus canadensis | Poaceae | perennial | 1.617 | 0.000 | 0.816 | 0.025 | 0.078 |
| Ambrosia artemisiifolia | Asteraceae | annual | 0.000 | 8.457 | 0.756 | 0.072 | 0.078 |
| Medicago lupulina | Fabaceae | biennial | 0.000 | 0.186 | 0.756 | 0.074 | 0.078 |
| Schizachyrium scoparium | Poaceae | perennial | 0.000 | 1.043 | 0.756 | 0.078 | 0.078 |

Plant indicators in low or high saprotroph sites

# Summary results

The following outputs display unified summary reports where needed for
tables, etc.

## Richness summary

``` r
list(
  its_rich_nb     = Anova(its_rich_glm, type = 2, test.statistic = "LR"),
  amf_rich_pois   = Anova(amf_rich_glm, type = 2, test.statistic = "LR"),
  patho_rich_pois = Anova(patho_rich_glm, type = 2, test.statistic = "LR"),
  sapro_rich_nb   = Anova(sapro_rich_glm, type = 2, test.statistic = "LR")
) %>% map(\(df) tidy(df)) %>% 
  bind_rows(.id = "guild_test") %>% 
  mutate(p.adj = if_else(term == "field_type", p.adjust(p.value, "fdr"), NA_real_),
         across(where(is.numeric), ~ round(.x, 3)),
         LRchisq_df = paste0(statistic, " (", df, ", 21)")) %>% 
  select(guild_test, term, LRchisq_df, p.value, p.adj) %>% 
  kable(format = "pandoc", caption = "Fungal OTU richness differences across field types accounting for sequencing depth.\nField type effects were evaluated using Type II Analysis of Deviance.\nP-values for field type were adjusted for multiple comparisons\nacross fungal groups using the Benjamini-Hochberg procedure.")
```

| guild_test      | term       | LRchisq_df     | p.value | p.adj |
|:----------------|:-----------|:---------------|--------:|------:|
| its_rich_nb     | depth_csq  | 8.128 (1, 21)  |   0.004 |    NA |
| its_rich_nb     | field_type | 40.506 (2, 21) |   0.000 | 0.000 |
| amf_rich_pois   | depth_csq  | 0.273 (1, 21)  |   0.601 |    NA |
| amf_rich_pois   | field_type | 10.219 (2, 21) |   0.006 | 0.012 |
| patho_rich_pois | depth_csq  | 11.541 (1, 21) |   0.001 |    NA |
| patho_rich_pois | field_type | 1.535 (2, 21)  |   0.464 | 0.530 |
| sapro_rich_nb   | depth_csq  | 5.801 (1, 21)  |   0.016 |    NA |
| sapro_rich_nb   | field_type | 5.948 (2, 21)  |   0.051 | 0.068 |

Fungal OTU richness differences across field types accounting for
sequencing depth. Field type effects were evaluated using Type II
Analysis of Deviance. P-values for field type were adjusted for multiple
comparisons across fungal groups using the Benjamini-Hochberg procedure.

## Shannon diversity summary

Fungal OTU Shannon Diversity differences across field types accounting
for sequencing depth. Field type effects were evaluated using Type II
Analysis of Deviance. P-values for field type were adjusted for multiple
comparisons across fungal groups using the Benjamini-Hochberg procedure.

``` r
list(
  its_shan_lm   = Anova(its_shan_lm, type = 2),
  amf_shan_lm   = Anova(amf_shan_lm, type = 2),
  patho_shan_lm = Anova(patho_shan_lm, type = 2),
  sapro_shan_lm = Anova(sapro_shan_lm, type = 2)
) %>% map(\(df) tidy(df)) %>% 
  bind_rows(.id = "guild_test") %>% 
  mutate(p.adj = if_else(term == "field_type", p.adjust(p.value, "fdr"), NA_real_),
         across(where(is.numeric), ~ round(.x, 3)),
         `F` = paste0(statistic, " (", df, ", 21)")) %>% 
  select(guild_test, term, `F`, p.value, p.adj) %>% 
  kable(format = "pandoc")
```

| guild_test    | term       | F              | p.value | p.adj |
|:--------------|:-----------|:---------------|--------:|------:|
| its_shan_lm   | depth_csq  | 0.866 (1, 21)  |   0.363 |    NA |
| its_shan_lm   | field_type | 6.549 (2, 21)  |   0.006 | 0.025 |
| its_shan_lm   | Residuals  | NA (21, 21)    |      NA |    NA |
| amf_shan_lm   | depth_csq  | 0.013 (1, 21)  |   0.911 |    NA |
| amf_shan_lm   | field_type | 13.269 (2, 21) |   0.000 | 0.002 |
| amf_shan_lm   | Residuals  | NA (21, 21)    |      NA |    NA |
| patho_shan_lm | depth_csq  | 3.799 (1, 21)  |   0.065 |    NA |
| patho_shan_lm | field_type | 1.355 (2, 21)  |   0.280 | 0.537 |
| patho_shan_lm | Residuals  | NA (21, 21)    |      NA |    NA |
| sapro_shan_lm | depth_csq  | 0.174 (1, 21)  |   0.680 |    NA |
| sapro_shan_lm | field_type | 0.95 (2, 21)   |   0.403 | 0.537 |
| sapro_shan_lm | Residuals  | NA (21, 21)    |      NA |    NA |

## Biomass (abundance) summary

Fungal biomass differences differences across field types. Field type
effects were evaluated using ANOVA or Analysis of Deviance. P-values for
field type were adjusted for multiple comparisons across fungal groups
using the Benjamini-Hochberg procedure.”

``` r
list(
  its_ma_lm   = anova(plfa_lm), 
  amf_ma_glm  = Anova(nlfa_glm, test.statistic = "LR"), 
  patho_ma_lm = anova(patho_ma_lm), 
  sapro_ma_lm = anova(sapro_ma_lm)
) %>% map(\(df) tidy(df) %>% select(term, statistic, df, p.value)) %>% 
  bind_rows(.id = "guild_test") %>% 
  mutate(p.adj = if_else(term == "field_type", p.adjust(p.value, "fdr"), NA_real_),
         across(where(is.numeric), ~ round(.x, 3)),
         `F` = paste0(statistic, " (", df, ", 21)")) %>% 
  select(guild_test, term, `F`, p.value, p.adj) %>% 
  kable(format = "pandoc")
```

| guild_test  | term       | F              | p.value | p.adj |
|:------------|:-----------|:---------------|--------:|------:|
| its_ma_lm   | field_type | 2.343 (2, 21)  |   0.120 | 0.239 |
| its_ma_lm   | Residuals  | NA (22, 21)    |      NA |    NA |
| amf_ma_glm  | field_type | 46.369 (2, 21) |   0.000 | 0.000 |
| patho_ma_lm | field_type | 1.712 (2, 21)  |   0.204 | 0.272 |
| patho_ma_lm | Residuals  | NA (22, 21)    |      NA |    NA |
| sapro_ma_lm | field_type | 0.376 (2, 21)  |   0.691 | 0.691 |
| sapro_ma_lm | Residuals  | NA (22, 21)    |      NA |    NA |

## Beta diversity summary

### Using abundance-scaled biomass

Fungal community differences differences among field types. Field type
effects were evaluated using Permanova. P-values for field type were
adjusted for multiple comparisons across fungal groups using the
Benjamini-Hochberg procedure.”

``` r
list(
  its_ma   = mva_its_ma$permanova,
  amf_ma   = mva_amf_ma$permanova,
  patho_ma = mva_patho_ma$permanova,
  sapro_ma = mva_sapro_ma$permanova
) %>% map(\(df) tidy(df) %>% select(term, pseudo_F = statistic, df, R2, p.value)) %>% 
  bind_rows(.id = "guild") %>% 
  mutate(p.adj = if_else(term == "field_type", p.adjust(p.value, "fdr"), NA_real_),
         across(where(is.numeric), ~ round(.x, 3)),
         `Pseudo_F_(df)` = paste0(pseudo_F, " (", df, ", 21)")) %>% 
  filter(term %in% c("dist_axis_1", "field_type")) %>% 
  select(guild, term, `Pseudo_F_(df)`, R2, p.value, p.adj) %>% 
  kable(format = "pandoc")
```

| guild    | term        | Pseudo_F\_(df) |    R2 | p.value | p.adj |
|:---------|:------------|:---------------|------:|--------:|------:|
| its_ma   | dist_axis_1 | 1.717 (1, 21)  | 0.062 |   0.027 |    NA |
| its_ma   | field_type  | 2.433 (2, 21)  | 0.176 |   0.000 | 0.001 |
| amf_ma   | dist_axis_1 | 1.303 (1, 21)  | 0.040 |   0.222 |    NA |
| amf_ma   | field_type  | 5.198 (2, 21)  | 0.318 |   0.000 | 0.001 |
| patho_ma | dist_axis_1 | 1.794 (1, 21)  | 0.062 |   0.088 |    NA |
| patho_ma | field_type  | 2.994 (2, 21)  | 0.208 |   0.002 | 0.003 |
| sapro_ma | dist_axis_1 | 2.023 (1, 21)  | 0.075 |   0.004 |    NA |
| sapro_ma | field_type  | 2.056 (2, 21)  | 0.152 |   0.000 | 0.001 |

### Using relative sequence abundance

Fungal community differences differences among field types. Field type
effects were evaluated using Permanova. P-values for field type were
adjusted for multiple comparisons across fungal groups using the
Benjamini-Hochberg procedure.”

``` r
list(
  its_ma   = mva_its$permanova,
  amf_ma   = mva_amf$permanova,
  patho_ma = mva_patho$permanova,
  sapro_ma = mva_sapro$permanova
) %>% map(\(df) tidy(df) %>% select(term, pseudo_F = statistic, df, R2, p.value)) %>% 
  bind_rows(.id = "guild") %>% 
  mutate(p.adj = if_else(term == "field_type", p.adjust(p.value, "fdr"), NA_real_),
         across(where(is.numeric), ~ round(.x, 3)),
         `Pseudo_F_(df)` = paste0(pseudo_F, " (", df, ", 21)")) %>% 
  filter(term %in% c("dist_axis_1", "field_type")) %>% 
  select(guild, term, `Pseudo_F_(df)`, R2, p.value, p.adj) %>% 
  kable(format = "pandoc")
```

| guild    | term        | Pseudo_F\_(df) |    R2 | p.value | p.adj |
|:---------|:------------|:---------------|------:|--------:|------:|
| its_ma   | dist_axis_1 | 1.739 (1, 21)  | 0.063 |   0.027 |    NA |
| its_ma   | field_type  | 2.536 (2, 21)  | 0.182 |   0.000 | 0.001 |
| amf_ma   | dist_axis_1 | 1.678 (1, 21)  | 0.056 |   0.114 |    NA |
| amf_ma   | field_type  | 3.731 (2, 21)  | 0.248 |   0.000 | 0.001 |
| patho_ma | dist_axis_1 | 1.607 (1, 21)  | 0.053 |   0.134 |    NA |
| patho_ma | field_type  | 3.887 (2, 21)  | 0.256 |   0.000 | 0.001 |
| sapro_ma | dist_axis_1 | 2.165 (1, 21)  | 0.079 |   0.003 |    NA |
| sapro_ma | field_type  | 2.169 (2, 21)  | 0.158 |   0.000 | 0.001 |

## Constrained analysis summary

Environmental drivers were identified via partial distance-based
Redundancy Analysis (db-RDA) using forward selection. Geographic
distance (PCoA Axis 1) was included as a conditional term to partial out
spatial effects. Radj2 represents the cumulative variance explained by
the final selected model. P-values are based on 1,999 permutations; Padj
reflects FDR correction within the guild.

### Adjusted R2

``` r
#= dbrda_r2
data.frame(
  guild = c("all_fungi", "amf", "saprotrophs"),
  r2    = round(c(mod_r2$adj.r.squared, amf_mod_r2$adj.r.squared, sapro_mod_r2$adj.r.squared), 3)
) %>% 
  kable(format = "pandoc")
```

| guild       |    r2 |
|:------------|------:|
| all_fungi   | 0.236 |
| amf         | 0.351 |
| saprotrophs | 0.260 |

``` r
rdf <- data.frame(
  guild = c("all_fungi", "amf", "saprotrophs"),
  rdf   = c(mod_inax["Residual", "Df"], amf_mod_inax["Residual", "Df"], sapro_mod_inax["Residual", "Df"])
)
```

### Selected constraining variables

``` r
list(
  all_fungi   = mod_step$anova,
  amf         = amf_mod_step$anova,
  saprotrophs = sapro_mod_step$anova
) %>% 
  map(\(df) df %>% 
        tidy() %>% 
        mutate(p.adj = p.adjust(p.value, "fdr"),
               across(where(is.numeric), ~ round(.x, 3)))) %>% 
  bind_rows(.id = "guild") %>% 
  left_join(rdf, by = join_by(guild)) %>% 
  mutate(`pseudo_F_(df)` = paste0(statistic, " (", df, ", ", rdf, ")"),
         term = str_remove(term, "\\+ ")) %>% 
  select(term, `pseudo_F_(df)`, p.value, p.adj) %>% 
  kable(format = "pandoc")
```

| term     | pseudo_F\_(df) | p.value | p.adj |
|:---------|:---------------|--------:|------:|
| gf_index | 2.077 (1, 8)   |   0.006 | 0.013 |
| pH       | 2.049 (1, 8)   |   0.009 | 0.013 |
| pl_rich  | 1.748 (1, 8)   |   0.028 | 0.028 |
| gf_index | 4.028 (1, 9)   |   0.002 | 0.004 |
| pH       | 3.133 (1, 9)   |   0.009 | 0.009 |
| OM       | 1.793 (1, 7)   |   0.011 | 0.015 |
| gf_index | 1.865 (1, 7)   |   0.008 | 0.015 |
| pl_rich  | 1.826 (1, 7)   |   0.009 | 0.015 |
| NO3      | 1.56 (1, 7)    |   0.041 | 0.041 |

### Global tests

``` r
list(
  all_fungi   = mod_glax,
  amf         = amf_mod_glax,
  saprotrophs = sapro_mod_glax
) %>% map(\(df) df %>% 
            tidy() %>% 
            filter(term != "Residual") %>% 
            mutate(p.adj = p.adjust(p.value, "fdr"),
                   across(where(is.numeric), ~ round(.x, 3)))) %>% 
  bind_rows(.id = "guild") %>% 
  left_join(rdf, by = join_by(guild)) %>% 
  mutate(`pseudo_F_(df)` = paste0(statistic, " (", df, ", ", rdf, ")")) %>% 
  select(guild, term, `pseudo_F_(df)`, p.value, p.adj) %>% 
  kable(format = "pandoc")
```

| guild       | term  | pseudo_F\_(df) | p.value | p.adj |
|:------------|:------|:---------------|--------:|------:|
| all_fungi   | Model | 2.151 (3, 8)   |       0 |     0 |
| amf         | Model | 4.01 (2, 9)    |       0 |     0 |
| saprotrophs | Model | 1.992 (4, 7)   |       0 |     0 |

### Component axes

``` r
list(
  all_fungi   = mod_inax,
  amf         = amf_mod_inax,
  saprotrophs = sapro_mod_inax
) %>% map(\(df) df %>% 
            tidy() %>% 
            filter(term != "Residual") %>% 
            mutate(p.adj = p.adjust(p.value, "fdr"),
                   across(where(is.numeric), ~ round(.x, 3)))) %>% 
  bind_rows(.id = "guild") %>% 
  left_join(rdf, by = join_by(guild)) %>% 
  mutate(`pseudo_F_(df)` = paste0(statistic, " (", df, ", ", rdf, ")"),
         term = str_remove(term, "\\+ ")) %>% 
  select(guild, term, `pseudo_F_(df)`, p.value, p.adj) %>% 
  kable(format = "pandoc")
```

| guild       | term   | pseudo_F\_(df) | p.value | p.adj |
|:------------|:-------|:---------------|--------:|------:|
| all_fungi   | dbRDA1 | 2.991 (1, 8)   |   0.001 | 0.003 |
| all_fungi   | dbRDA2 | 2.13 (1, 8)    |   0.004 | 0.004 |
| all_fungi   | dbRDA3 | 1.96 (1, 8)    |   0.004 | 0.004 |
| amf         | dbRDA1 | 5.546 (1, 9)   |   0.000 | 0.001 |
| amf         | dbRDA2 | 2.749 (1, 9)   |   0.002 | 0.002 |
| saprotrophs | dbRDA1 | 2.844 (1, 7)   |   0.000 | 0.002 |
| saprotrophs | dbRDA2 | 2.396 (1, 7)   |   0.002 | 0.002 |
| saprotrophs | dbRDA3 | 2.297 (1, 7)   |   0.002 | 0.002 |
| saprotrophs | dbRDA4 | 1.774 (1, 7)   |   0.018 | 0.018 |

### Biplot panels

All soil fungi

``` r
fig6a <- 
  ggplot(mod_scor_site, aes(x = dbRDA1, y = dbRDA2)) +
  geom_segment(data = mod_scor_bp, 
               aes(x = origin, xend = dbRDA1, y = origin, yend = dbRDA2), 
               arrow = arrow(length = unit(2, "mm"), type = "closed"),
               color = c("darkblue", "darkblue", "gray20", "gray20")) +
  geom_text(data = mod_scor_bp, 
            aes(x = labx, y = laby, label = envlabs), 
            # nudge_x = c(-0.1, 0.1, 0), nudge_y = c(0.06, -0.06, 0),
            size = 3, color = "black", fontface = 2) +
  geom_point(aes(fill = field_type), size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  labs(
    x = paste0("db-RDA 1 (", mod_axpct[1], "%)"),
    y = paste0("db-RDA 2 (", mod_axpct[2], "%)")) +
  lims(x = c(-0.95,1.6)) +
  scale_fill_manual(values = ft_pal[2:3]) +
  scale_y_continuous(breaks = c(-1, 0, 1)) +
  theme_ord +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
```

AMF

``` r
fig6b <-
  ggplot(amf_mod_scor_site, aes(x = -1 * (dbRDA1), y = dbRDA2)) +
  geom_segment(data = amf_mod_scor_bp,
               aes(x = origin, xend = -1 * (dbRDA1), y = origin, yend = dbRDA2),
               arrow = arrow(length = unit(2, "mm"), type = "closed"),
               color = c("darkblue", "darkblue", "gray20")) +
  geom_text(data = amf_mod_scor_bp,
            aes(x = -1 * (labx), y = laby, label = envlabs),
            # nudge_x = (c(0.05, 0.2, -0.2)), nudge_y = c(0.1, 0.04, -0.04),
            size = 3, color = "gray20", fontface = 2) +
  geom_point(aes(fill = field_type), size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  labs(
    x = paste0("db-RDA 1 (", amf_mod_axpct[1], "%)"),
    y = paste0("db-RDA 2 (", amf_mod_axpct[2], "%)")) +
  lims(x = c(-1.2,1.2)) +
  scale_fill_manual(values = ft_pal[2:3]) +
  theme_ord +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
```

Pathogens, PCoA fig

``` r
fig6c <-
  ggplot(patho_mod_scor_site, aes(x = (MDS1), y = MDS2)) +
  geom_point(aes(fill = field_type), size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  labs(
    x = paste0("PCoA 1 (", patho_mod_eig[1], "%)"),
    y = paste0("PCoA 2 (", patho_mod_eig[2], "%)")) +
  # lims(x = c(-1.1,1.05), y = c(-1.6,0.9)) +
  scale_fill_manual(values = ft_pal[2:3]) +
  scale_y_continuous(breaks = c(-1, 0, 1)) +
  theme_ord +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
```

Saprotrophs

``` r
fig6d <-
  ggplot(sapro_mod_scor_site, aes(x = (dbRDA1), y = dbRDA2)) +
  geom_segment(data = sapro_mod_scor_bp,
               aes(x = origin, xend = (dbRDA1), y = origin, yend = dbRDA2),
               arrow = arrow(length = unit(2, "mm"), type = "closed"),
               color = c("gray20", "gray20", "darkblue", "darkblue", "gray20")) +
  geom_text(data = sapro_mod_scor_bp,
            aes(x = (labx), y = laby, label = paste0("bold(", envlabs, ")")), parse = TRUE,
            # nudge_x = (c(0.05, 0.2, -0.2)), nudge_y = c(0.1, 0.04, -0.04),
            size = 3, color = "gray20") +
  geom_point(aes(fill = field_type), size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  labs(
    x = paste0("db-RDA 1 (", sapro_mod_axpct[1], "%)"),
    y = paste0("db-RDA 2 (", sapro_mod_axpct[2], "%)")) +
  lims(x = c(-0.9,1.7)) +
  scale_fill_manual(name = "Field type", values = ft_pal[2:3]) +
  theme_ord +
  theme(legend.position = c(0.98, 0.03),
        legend.justification = c(1, 0),
        legend.title = element_text(size = 9, face = 1),
        legend.text = element_text(size = 8, face = 1),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2),
        legend.key = element_rect(fill = "white"),
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
```

#### Unified figure

Display results of constrained analyses

``` r
fig6up <- (fig6a | plot_spacer() | fig6b) +
  plot_layout(widths = c(0.50, 0.01, 0.50))
fig6dn <- (fig6c | plot_spacer() | fig6d) +
  plot_layout(widths = c(0.50, 0.01, 0.50))
fig6 <- (fig6up / fig6dn) +
  plot_layout(widths = c(0.50, 0.50)) +
  plot_annotation(tag_levels = 'A')
```

``` r
fig6
```

![](resources/fungal_ecology_files/figure-gfm/fig6-1.png)<!-- -->

Fungal community ordinations which are constrained or unconstrained by
explanatory variables. Panels show results for all soil fungi **a**, amf
**b**, pathogens **c**, and saprotrophs **d**. Percent of constrained
(db-RDA) and unconstrained (PCoA) variation explained is shown with axis
labels. For explanatory variables with significant community
correlations, blue arrows show the grass-forb index with labels
indicating the direction of relative increase in C4 grasses or forbs,
respectively, along the index. The black arrows show other significant
constraining variables. Points show locations of restored fields (green)
and remnant fields (blue) in Wisconsin.
