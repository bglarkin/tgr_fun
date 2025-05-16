#' ---
#' title: "Correlations: guilds and environment"
#' author: "Beau Larkin\n"
#' date: "Last updated: `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   github_document:
#'     toc: true
#'     toc_depth: 2
#'     df_print: paged
#' ---
#'
#' # Description
#' Site‑average guild abundances were correlated with plant functional groups. Sequence 
#' data are compositional, so abundances were centered via the additive‑log‑ratio (ALR) 
#' method [Gloor 2017](https://www.frontiersin.org/articles/10.3389/fmicb.2017.02224) and
#' [Greenacre 2021](https://doi.org/10.3389/fmicb.2021.727398) before differential tests.
#' Plant cover (10×1 m quadrats, WI only; M. Healy 2016) was analysed without log‑transformation
#' because values are not instrument‑bounded and may exceed 100% [Gloor 2017].
#' Trait data was obtained from the [TRY database](https://www.try-db.org)  
#' ([Kattge 2010](https://doi.org/10.1111/j.2041-210X.2010.00067.x),
#' [Kattge 2011](https://doi.org/10.1111/j.1365-2486.2011.02451.x)).

#' # Packages and libraries
packages_needed <- c(
  "tidyverse", "colorspace", "vegan", "knitr", "conflicted",
  "geosphere", "ape", "performance", "patchwork"
)
to_install <- setdiff(packages_needed, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)
invisible(lapply(packages_needed, library, character.only = TRUE))

#' ## Root path function
root_path <- function(...) rprojroot::find_rstudio_root_file(...)

#+ conflicts,message=FALSE
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

#+ styles
source(root_path("resources", "styles.R"))

#' # Data
# Data ———————— ####
#' 
#' ## Site metadata and design
sites <- read_csv(paste0(getwd(), "/clean_data/sites.csv"), show_col_types = FALSE) %>% 
  mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant")))

#' ## Sites-species tables
#' List *spe* holds average sequence abundances from repeated samples in each site. 
#' CSV files were produced in `sequence_data.R`. 
spe <- list(
  its_avg     = read_csv(root_path("clean_data/spe_ITS_avg.csv"), show_col_types = FALSE),
  amf_avg     = read_csv(root_path("clean_data/spe_18S_avg.csv"), show_col_types = FALSE),
  amf_avg_uni = read_delim(root_path("otu_tables/18S/18S_avg_4unifrac.tsv"), show_col_types = FALSE)
)
#' 
#' ## Microbial species metadata
spe_meta <- list(
  its = read_csv(root_path("clean_data/spe_ITS_metadata.csv"), show_col_types = FALSE),
  amf = read_csv(root_path("clean_data/spe_18S_metadata.csv"), show_col_types = FALSE)
) %>% 
  map(. %>% mutate(across(everything(), ~ replace_na(., "unidentified"))))

#' ## Plant data
#' Abundance in functional groups and by species are only available from Wisconsin sites. 
#' Only C4_grass and forbs are used. Others: C3_grass, legume, and shrubTree were found 
#' previously to have high VIF in models or were not chosen in forward selection. 
pfg <- read_csv(root_path("clean_data", "plant_traits.csv"), show_col_types = FALSE) %>% 
  select(field_name, C4_grass, forb)

#' # Data wrangling
# Data wrangling ———————— ####
#' 
#' - C4 grass and forb cover are transformed into a single index using PCA in restored sites only.
#' - The OTU abundance tables must be wrangled to perform a log-ratio transformation, which reduces
#' data skewness and compositionality bias (Aitchison 1986, Gloor et al. 2017). 
#' The transformation will be applied across guilds for whole soil fungi and families for AMF. 
#' Raw abundances are kept for plotting.
#' Abundance data are also joined with site and env paramaters to facilitate downstream analyses.
#' 
#' ## Grass-forb index
#' C4 grass and forb cover are highly correlated (*r* = `r round(cor(pfg$forb, pfg$C4_grass), 2)`) 
#' in restored prairies. In models or constrained ordinations, they are collinear and cannot be
#' used simultaneously. An index of grass-forb cover is created to solve this problem. 
pfg_pca <- 
  pfg %>% 
  left_join(sites %>% select(field_name, field_type), by = join_by(field_name)) %>% 
  filter(field_type == "restored") %>% 
  select(-field_type) %>% 
  column_to_rownames(var = "field_name") %>% 
  decostand(method = "standardize") %>% 
  rda()
pfg_pca %>% summary() # 92% variation on first axis
gf_index = scores(pfg_pca, choices = 1, display = "sites") %>% 
  data.frame() %>% 
  rename(gf_index = PC1) %>% 
  rownames_to_column(var = "field_name")
#' 
#' ## Whole soil fungi
#' ### Raw abundances
its_guab <- 
  spe$its_avg %>% 
  pivot_longer(starts_with("otu"), names_to = "otu_num", values_to = "abund") %>% 
  left_join(spe_meta$its %>% select(otu_num, primary_lifestyle), by = join_by(otu_num)) %>% 
  group_by(field_name, primary_lifestyle) %>% summarize(abund = sum(abund), .groups = "drop") %>% 
  arrange(field_name, -abund) %>% 
  pivot_wider(names_from = "primary_lifestyle", values_from = "abund") %>% 
  mutate(saprotroph = rowSums(across(ends_with("_saprotroph")))) %>% 
  select(field_name, unidentified, saprotroph, plant_pathogen, everything()) %>% 
  select(-ends_with("_saprotroph"))
its_guab_pfg <- 
  its_guab %>% 
  left_join(pfg, by = join_by(field_name)) %>% 
  left_join(gf_index, by = join_by(field_name)) %>% 
  left_join(sites %>% select(field_name, field_type, region, yr_since), by = join_by(field_name)) %>% 
  select(field_name, field_type, yr_since, region, everything())
#' 
#' ### Log ratio transformed abundances
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
#' 
#' ## AMF
#' ### Raw abundances
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
#' 
#' ### Log ratio transformed abundances
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

#' # AMF abundance in families
# AMF abundance in families ———————— ####
#' Display raw abundances in a table but separate means with log ratio transformed data
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
#' 
#' Test RCLR transformed abundances across field types for each family
glom_lm <- lm(Glomeraceae ~ field_type, data = amf_fmlr_pfg)
summary(glom_lm)
#' NS
clar_lm <- lm(Claroideoglomeraceae ~ field_type, data = amf_fmlr_pfg)
summary(clar_lm)
performance::check_distribution(clar_lm) 
leveneTest(Claroideoglomeraceae ~ field_type, data = amf_fmlr_pfg)
TukeyHSD(aov(Claroideoglomeraceae ~ field_type, data = amf_fmlr_pfg))
#' Model R2_adj 0.24, p<0.02
ggplot(amf_fmlr_pfg, aes(x = field_type, y = Paraglomeraceae)) + geom_boxplot()
para_lm <- lm(Paraglomeraceae ~ field_type, data = amf_fmlr_pfg)
summary(para_lm)
#' NS
dive_lm <- lm(Diversisporaceae ~ field_type, data = amf_fmlr_pfg)
summary(dive_lm)
#' NS
giga_lm <- lm(Gigasporaceae ~ field_type, data = amf_fmlr_pfg)
summary(giga_lm)
#' NS





# Pathogens and plants
# Formally incorporate into this script later...
ggplot(its_guab_pfg %>% filter(field_type == "restored", region != "FL"), aes(x = gf_index, y = plant_pathogen)) +
  geom_smooth(method = "lm") +
  geom_point()


gf_patho_lm <- lm(plant_pathogen ~ gf_index, data = its_gulr_pfg %>% filter(field_type == "restored", region != "FL"))
#' Diagnostics
par(mfrow = c(2,2))
plot(gf_patho_lm) 
#' Some residual structure, but a smooth qq fit. 
#' Minor leverage with point 4 pulls the slope to more level, risking a type II error rather than type I. 
#' Model is still significant with point 4 removed. 
performance::check_distribution(gf_patho_lm) 
#' Response and residuals normal, no transformations warranted and linear model appropriate.
summary(gf_patho_lm)


