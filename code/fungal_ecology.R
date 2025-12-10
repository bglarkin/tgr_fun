#' ---
#' title: "Results: Soil Fungal Communities"
#' author: "Beau Larkin\n"
#' date: "Last updated: `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   github_document:
#'     toc: true
#'     toc_depth: 2
#' ---
#'
#' # Description
#' **Scope** – Biomass (PLFA/NLFA), OTU richness and diversity, and β‑diversity of soil fungi across corn, 
#' restored, and remnant prairie fields.
#' 
#' **Alpha diversity** – 97 %-OTUs (ITS & 18S); site means are replicates; means separation model selection
#' based on response and residuals distributions; sequencing depth used as covariate 
#' per [Bálint 2015](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.13018) when warranted; 
#' pairwise LSMs via *emmeans*.
#' 
#' **Beta diversity** – Workflow after [Song 2015](https://doi.org/10.1371/journal.pone.0127234):  
#'  1. PCoA of Bray (ITS) or UNIFRAC (18S) distances
#'  1. homogeneity test diagnostics 
#'  1. PERMANOVA (+ pairwise) 
#' 
#' Cartesian inter‑site distance enters models as a covariate per [Redondo 2020](https://doi.org/10.1093/femsec/fiaa082).
#' 
#' # Packages and libraries
# Libraries ———————— ####
#+ packages,message=FALSE
packages_needed <- c(
  "colorspace", "emmeans", "gridExtra", "knitr", "tidyverse", "vegan",
  "rprojroot", "phyloseq", "ape", "phangorn", "geosphere", "conflicted",
  "ggpubr", "patchwork", "car", "performance", "boot", "indicspecies",
  "MASS", "DHARMa", "broom", "rlang", "rsq", "purrr", "sandwich", "lmtest"
)

to_install <- setdiff(packages_needed, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)
invisible(lapply(packages_needed, library, character.only = TRUE))

#' ## Root path function
root_path <- function(...) rprojroot::find_rstudio_root_file(...)

#+ conflicts,message=FALSE
conflicts_prefer(
  dplyr::filter(),
  dplyr::select(),
  dplyr::where(),
  vegan::diversity(),
  purrr::map()
)
#' 
#+ graphics_styles
source(root_path("resources", "styles.R"))
#' 

#' # Data
# Data ———————— ####
#' 
#' ## Site metadata and design
sites <- read_csv(root_path("clean_data/sites.csv"), show_col_types = FALSE) %>% 
  mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant")))
#' 
#' ### Wrangle site metadata
#' Intersite geographic distance will be used as a covariate in clustering. 
#' Raw coordinates in data file aren't distances; convert to distance matrix and summarize with PCoA
field_dist <- as.dist(distm(sites[, c("long", "lat")], fun = distHaversine))
field_dist_pcoa <- pcoa(field_dist)
field_dist_pcoa$values[c(1,2), c(1,2)] %>% 
  kable(format = "pandoc")
#' First axis of geographic distance PCoA explains 91% of the variation among sites. 
sites$dist_axis_1 <- field_dist_pcoa$vectors[, 1]
#' 
#' ## Fatty Acids: Biomass
#' Use only 18.2 for soil fungi
fa <- read_csv(root_path("clean_data/plfa.csv"), show_col_types = FALSE) %>% 
  rename(fungi_18.2 = fa_18.2) %>% 
  select(field_name, fungi_18.2, amf) %>%
  left_join(
    sites %>% select(field_name, field_type),
    by = join_by(field_name)
  )
#' 
#' ## Sites-species tables
#' CSV files were produced in `sequence_data.R`. Average sequence abundance at sites included here.
#' Amf_avg_uni table is in species-samples format to enable use of `Unifrac()` later.
its_avg     = read_csv(root_path("clean_data/spe_ITS_avg.csv"), show_col_types = FALSE)
amf_avg     = read_csv(root_path("clean_data/spe_18S_avg.csv"), show_col_types = FALSE)
#' 
#' ## Microbial species metadata
its_meta = read_csv(root_path("clean_data/spe_ITS_metadata.csv"), show_col_types = FALSE) %>% 
  mutate(primary_lifestyle = case_when(str_detect(primary_lifestyle, "_saprotroph$") ~ "saprotroph",
                                       TRUE ~ primary_lifestyle),
         across(everything(), ~ replace_na(., "unidentified")))
amf_meta = read_csv(root_path("clean_data/spe_18S_metadata.csv"), show_col_types = FALSE) %>% 
  mutate(across(everything(), ~ replace_na(., "unidentified")))
#' 
#' ### Wrangle additional species and metadata objects
# Species data wrangling ———————— ####
#' 
#' 1. Proportional species abundance corrected for site biomass
#' 1. Unifrac products for AMF
#' 1. Phyloseq products to process the Unifrac distance cleanly
#' 
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
#' 
#' ## Plant data
#' Abundance in functional groups and by species are only available from Wisconsin sites. 
#' Only C4_grass and forbs are used. Others: C3_grass, legume, and shrubTree were found 
#' previously to have high VIF in models or were not chosen in forward selection. 
pfg <- read_csv(root_path("clean_data", "plant_traits.csv"), show_col_types = FALSE) 
#' 
#' ## Soil properties
soil <- read_csv(root_path("clean_data/soil.csv"), show_col_types = FALSE)[-c(26:27), ]
#' 
#' # Species, environment, and metadata
# Metadata wrangling ———————— ####
#' 
#' ## Plant communities
#' ### Plant species richness
plant <- read_csv(root_path("clean_data/plant_avg.csv"), show_col_types = FALSE)
prich <- plant %>% 
  select(-BARESOIL, -LITTER, -ROSA, -SALIX) %>% # remove non-species entries
  rowwise() %>% 
  mutate(pl_rich = sum(c_across(where(is.numeric)) > 0)) %>% 
  select(field_name = SITE, pl_rich) %>% 
  left_join(sites, by = join_by(field_name)) %>% 
  filter(region != "FL", field_type == "restored") %>% 
  ungroup()
#' Visualize how richness varies in fields...restored in Wisconsin only as used later
with(prich, cor.test(yr_since, pl_rich))
#' Years since restoration isn't obviously related to plant species richness. The trend is fewer
#' species with older fields, but the correlation isn't significant, again KORP is a high 
#' leverage point.
#' 
#' ## Plant functional groups
#' C4 grass and forb cover are transformed into a single index using PCA in restored sites only.
#' Abundance data are also joined with site and env paramaters to facilitate downstream analyses.
#' 
#' ### Grass-forb index
#' C4 grass and forb cover are highly correlated (*r* = `r round(cor(pfg$forb, pfg$C4_grass), 2)`) 
#' in restored prairies. In models or constrained ordinations, they are collinear and cannot be
#' used simultaneously. An index of grass-forb cover is created to solve this problem. 
pfg_pca <- 
  pfg %>%
  select(field_name, C4_grass, forb) %>% 
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
#' Are field age and gf_index correlated?
gfi_yrs <- gf_index %>% 
  left_join(sites %>% select(field_name, yr_since), by = join_by(field_name)) %>% 
  arrange(-gf_index)
with(gfi_yrs, cor.test(yr_since, gf_index))
#' The relatively strong correlation suggests that different restoration methods over time
#' are still reflected in plant composition. Years since restoration is highly related to 
#' plant community change. 
#' 
#' Visualize grass forb gradient compared with plant composition, grass and forb cover,
#' and years since restoration. 
pfg_comp <- 
  pfg %>% 
  select(field_name, C3_grass:shrubTree) %>% 
  rowwise() %>% 
  mutate((across(where(is.numeric), ~ .x / sum(c_across(where(is.numeric))))) * 100) %>% 
  ungroup() %>% 
  pivot_longer(C3_grass:shrubTree, names_to = "pfg", values_to = "pct_comp") %>% 
  left_join(sites, by = join_by(field_name)) %>% 
  left_join(gf_index, by = join_by(field_name)) %>% 
  filter(field_type == "restored", region != "FL") %>% 
  select(field_name, yr_since, gf_index, pfg, pct_comp) %>% 
  mutate(pfg = factor(pfg, levels = rev(c("forb", "C4_grass", "C3_grass", "legume", "shrubTree")),
                      labels = rev(c("forb", "grass (C4)", "grass (C3)", "legume", "shrub, tree"))))
pfg_comp_fig <- 
  ggplot(pfg_comp, aes(x = fct_reorder(field_name, -gf_index), y = pct_comp, group = pfg)) +
  geom_col(aes(fill = pfg)) +
  labs(x = NULL, y = "Percent composition") +
  scale_fill_discrete_qualitative(name = "Functional group", palette = "pfg-col", rev = TRUE) +
  theme_cor+
  theme(plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
pfg_pct <- 
  pfg %>% 
  select(field_name, C4_grass, forb) %>%
  pivot_longer(C4_grass:forb, names_to = "pfg", values_to = "pct_cvr") %>% 
  left_join(sites, by = join_by(field_name)) %>% 
  left_join(gf_index, by = join_by(field_name)) %>% 
  filter(field_type == "restored", region != "FL") %>% 
  select(field_name, yr_since, gf_index, pfg, pct_cvr)  %>% 
  mutate(pfg = factor(pfg, levels = rev(c("C4_grass", "forb")),
                      labels = c("forb", "grass (C4)")))
gf_pct_fig <- 
  ggplot(pfg_pct, aes(x = fct_reorder(field_name, -gf_index), y = pct_cvr, group = pfg)) +
  geom_step(aes(color = pfg)) +
  geom_point(aes(color = pfg)) +
  scale_color_manual(name = "Functional group", values = gfi_cols) +
  labs(x = NULL, y = "Percent cover") +
  theme_cor+
  theme(plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
yrs_fig <- 
  ggplot(gfi_yrs %>% mutate(grp = "a"), aes(x = fct_reorder(field_name, -gf_index), y = yr_since, group = grp)) +
  geom_step() +
  geom_point() +
  labs(x = NULL, y = "Yrs. since resto.") +
  theme_cor +
  theme(plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1.1))
# Wrangle gf_index and label positions to align panel D with previous x axes
loc_space <- (max(gfi_yrs$gf_index) - min(gfi_yrs$gf_index)) / length(gfi_yrs$gf_index)
gfi_loc_data <- 
  bind_rows(
    gfi_yrs %>% mutate(grp = "Value"),
    gfi_yrs %>% 
      mutate(temp = c(max(gfi_yrs$gf_index) - (loc_space * 0.15), # Alignment not exact, may need post-production work
                      max(gfi_yrs$gf_index) - (loc_space * (c(1:9) * 1.15))),
             grp = "Site") %>% 
      select(-gf_index, field_name, gf_index = temp, yr_since)
  ) %>% 
  mutate(grp = factor(grp, levels = c("Value", "Site")))
gfi_seg_data <- 
  gfi_yrs %>% 
  mutate(y0 = 1) %>% 
  select(field_name, x0 = gf_index, y0) %>% 
  left_join(
    gfi_loc_data %>% 
      filter(grp == "Site") %>% 
      mutate(y1 = 2) %>% 
      select(field_name, x1 = gf_index, y1),
    by = join_by(field_name)
    ) %>% 
  arrange(x0)
gfi_loc_fig <- 
  ggplot(gfi_loc_data, aes(x = -gf_index, y = grp)) + 
  geom_point(shape = NA) +
  geom_segment(data = gfi_seg_data, aes(x = -x0, y = y0, xend = -x1, yend = y1)) +
  geom_point() +
  labs(y = NULL, x = "Grass-forb index") +
  theme_cor +
  theme(plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1.1))
#' 
#' #### Unified figure
#+ pfg_fig_patchwork,warning=FALSE
pfg_pct_fig <- (pfg_comp_fig / plot_spacer() / gf_pct_fig / plot_spacer() / yrs_fig / plot_spacer() / gfi_loc_fig) +
  plot_layout(heights = c(1,0.01,1,0.01,0.5,0.01,0.3))  +
  plot_annotation(tag_levels = 'A') 
#+ pfg_fig,warning=FALSE,fig.height=7,fig.width=7
pfg_pct_fig
#+ pfg_fig_save,warning=FALSE,echo=FALSE
ggsave(root_path("figs", "figS8.png"),
       plot = pfg_pct_fig,
       width = 7.5,
       height = 7.5,
       units = "in",
       dpi = 600)
#' 
#' ## Whole soil fungi
#' Wrangle data to produce proportional biomass in guilds for its and families for amf
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
#' Wrangle a second set to compare raw sequence abundances and proportion of biomass values together
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
#' 
#' #### AMF
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
#' 
#' # Functions
#' Executed from a separate script to save lines here; to view the function navigate to 
#' `functions.R` in the code folder, accessible from the root dir of the repo.
# Functions ———————— ####
source(root_path("code", "functions.R"))


#' 
#' # Whole Soil Fungi
# Whole soil fungi ———————— ####
#' 
#' ## Diversity Indices
#+ its_diversity
its_div <- calc_div(its_avg, sites)
#' 
#' ### Richness
#' Account for sequencing depth as a covariate. Test covar transformations for best model performance.
its_rich_covar <- covar_shape_test(
  data  = its_div,
  y     = "richness",       
  covar = "depth",   
  group = "field_type"           
)
its_rich_covar$compare
#' Best model does not use transformation of covariate.
its_rich_lm <- lm(richness ~ depth + field_type, data = its_div) # Interaction NS (not shown)
#' Diagnostics
#+ its_rich_covar_diagnostics,warning=FALSE,fig.width=7,fig.height=9
its_rich_covar$diagnostics
#' Long tails, some midrange structure, no leverage points
distribution_prob(its_rich_lm)
#' residuals distribution normal or long-tailed, response log
leveneTest(richness ~ field_type, data = its_div) %>% as.data.frame() %>% kable(format = "pandoc")
leveneTest(residuals(its_rich_lm) ~ its_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc")
#' Residuals/response distributions do not suggest the need for transformation.
#' Levene's p > 0.05 → fail to reject = variances can be considered equal across groups.
#' 
#' Model results, group means, and post-hoc. Use Type II SS for test of variables due to unbalanced design.
its_rich_covar$anova_t2
#' Sequence depth is significant, less so than field type. Check relationship of depth and field type. 
its_div %>% 
    group_by(field_type) %>% 
    summarize(across(c(depth, richness), ~ round(mean(.x), 0))) %>% 
    kable(format = "pandoc")
#' Sequence depth isn't obviously related to field type, but they're weakly inversely related. 
#' The possibility for an interaction between depth and field_type was tested, it did not improve 
#' model fit (tested with anova(), not shown).
#' Proceed with means separation by obtaining estimated marginal means for field type.
#' Arithmetic means calculated in this case.
its_rich_em <- emmeans(its_rich_lm, ~ field_type, type = "response")
#' Results tables below show the emmeans summary of group means and confidence intervals, 
#' with sequencing depth as a covariate, and the post hoc contrast of richness among field types. 
#' Main effect and covariate in model significant after p value adjustment (see summary section): pairwise
#' contrast warranted.
#+ its_rich_em_summary,echo=FALSE
kable(summary(its_rich_em), 
      format = "pandoc", 
      caption = "Confidence level used: 0.95")
#+ its_rich_em_posthoc,echo=FALSE
kable(pairs(its_rich_em), 
      format = "pandoc", 
      caption = "P value adjustment: tukey method for comparing a family of 3 estimates")
#' OTU richness in cornfields is significantly less than in restored or remnant fields (p<0.001), which 
#' don't differ. 
#+ its_richness_fig,fig.width=4,fig.height=4
its_rich_fig <- 
  ggplot(summary(its_rich_em), aes(x = field_type, y = emmean)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = emmean, ymax = upper.CL), width = 0, linewidth = lw) +
  geom_text(aes(y = upper.CL, label = c("a", "b", "b")), vjust = -1.5, family = "sans", size = 4) +
  labs(x = NULL, y = "Richness") +
  lims(y = c(0, 760)) +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1))
#' 
#' ### Shannon's diversity
#' Account for sequencing depth as a covariate. Test covariate transformations.
its_shan_covar <- covar_shape_test(
  data  = its_div,
  y     = "shannon",       
  covar = "depth",   
  group = "field_type"           
)
its_shan_covar$compare
#' Log transformation of depth selected; difference between models is slight. Produce model 
#' with centered, log transformed depth. 
its_shan_lm <- lm(shannon ~ depth_clg + field_type, 
                  data = its_div %>% mutate(depth_clg = log(depth) - mean(log(depth))))
#' Diagnostics
#+ its_shan_covar_diagnostics,warning=FALSE,fig.width=7,fig.height=9
its_shan_covar$diagnostics
#' Lots of residual structure, no leverage points, no evidence for increasing mean/var relationship.
distribution_prob(its_shan_lm)
#' residuals distribution most likely cauchy/normal; symmetric but long tails, response gamma
leveneTest(shannon ~ field_type, data = its_div) %>% as.data.frame() %>% kable(format = "pandoc")
leveneTest(residuals(its_shan_lm) ~ its_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc")
#' Residuals distribution does not suggest the need for transformation.
#' Levene's p > 0.05 → fail to reject = variances can be considered equal.
#' Response more suspicious. Examine CV in groups to assess changes in variance. 
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
#' CV constant to declining.
#' Relatively low p value likely due to unequal variance in restored and remnant despite similar means.
#' Unbalanced data. 
#' 
#' Model results, group means, and post-hoc. Type II SS used due to unbalanced design.
its_shan_covar$anova_t2
#' Sequence depth is not a significant predictor of Shannon diversity.
#' Proceed with means separation by obtaining estimated marginal means for field type.
#' Arithmetic means calculated in this case.
its_shan_em <- emmeans(its_shan_lm, ~ field_type, type = "response")
#' Results tables below show the emmeans summary of group means and confidence intervals,
#' with sequencing depth as a covariate, and the post hoc contrast of richness among field types. 
#+ its_shan_em_summary,echo=FALSE
kable(summary(its_shan_em), 
      format = "pandoc", 
      caption = "Confidence level used: 0.95")
#+ its_shan_em_posthoc,echo=FALSE
kable(pairs(its_shan_em), 
      format = "pandoc", 
      caption = "P value adjustment: tukey method for comparing a family of 3 estimates")
#' Shannon's diversity in cornfields is significantly less than in restored or remnant fields, which 
#' don't differ.
#+ its_shan_fig,fig.width=4,fig.height=4
its_shan_fig <- 
  ggplot(summary(its_shan_em), aes(x = field_type, y = emmean)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = emmean, ymax = upper.CL), width = 0, linewidth = lw) +
  geom_text(aes(y = upper.CL, label = c("a", "b", "b")), vjust = -1.5, family = "sans", size = 4) +
  labs(x = NULL, y = "Shannon diversity") +
  lims(y = c(0, 160)) +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
#' 
#' ## Abundance (PLFA biomass)
plfa_lm <- lm(fungi_18.2 ~ field_type, data = fa)
par(mfrow = c(2,2))
plot(plfa_lm) 
#' variance differs slightly in groups. Tails on qq plot diverge, lots of groups structure visible.
distribution_prob(plfa_lm)
#' Residuals distribution fits normal, response normal-ish
leveneTest(residuals(plfa_lm) ~ fa$field_type) %>% as.data.frame() %>% kable(format = "pandoc") # No covariate, response and residuals tests equivalent
#' Residuals distribution does not suggest the need for transformation.
#' Levene's p > 0.05 → fail to reject = variances can be considered equal.
#' 
#' Model results, group means, and post-hoc, with arithmetic means from emmeans
anova(plfa_lm)
plfa_em <- emmeans(plfa_lm, ~ field_type, type = "response")
#+ plfa_em_summary,echo=FALSE
kable(summary(plfa_em), 
      format = "pandoc", 
      caption = "Confidence level used: 0.95")
#+ plfa_em_posthoc,echo=FALSE
kable(pairs(plfa_em), 
      format = "pandoc", 
      caption = "P value adjustment: tukey method for comparing a family of 3 estimates")
#+ plfa_fig,fig.width=4,fig.height=4,fig.align='center'
plfa_fig <- 
  ggplot(summary(plfa_em), aes(x = field_type, y = emmean)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = emmean, ymax = upper.CL), width = 0, linewidth = lw) +
  labs(x = NULL, y = expression(Biomass~(nmol[PLFA]%*%g[soil]^-1))) +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1.1))
#' 
#' ## Beta Diversity
#' 
#' 1. Using biomass-weighted relative abundance [Waller et al. 2023](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/1365-2745.14392)
#' 1. Using sequence-based relative abundance [McKnight et al. 2019](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13115)
#' 1. Contrast the two with procrustes
#' 
#' ### ITS, biomass-weighted relative abundance
#+ its_ord_ma
d_its_ma <- its_avg_ma %>% 
  data.frame(row.names = 1) %>% 
  vegdist("bray")
mva_its_ma <- mva(d = d_its_ma, env = sites)
#+ its_ord_ma_results
mva_its_ma$dispersion_test
mva_its_ma$permanova
mva_its_ma$pairwise_contrasts[c(1,3,2), c(1,2,4,3,8)] %>% 
  kable(format = "pandoc", caption = "Pairwise permanova contrasts")
#' No eignevalue correction was needed. Two relative eigenvalues exceeded broken stick model. 
#' Based on the homogeneity of variance test, the null hypothesis of equal variance among groups is 
#' accepted across all clusters and in pairwise comparison of clusters (both p>0.05), supporting the application of 
#' a PERMANOVA test. 
#' 
#' Clustering revealed that community variation was related to geographic distance, the covariate in 
#' the model. With geographic distance accounted for, the test variable 'field type' significantly explained 
#' variation in fungal communities, with a post-hoc test revealing that communities in corn fields differed from
#' communities in restored and remnant fields. 
#' 
#' Plotting results: 
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
    x = paste0("Axis 1 (", mva_its_ma$axis_pct[1], "%)"),
    y = paste0("Axis 2 (", mva_its_ma$axis_pct[2], "%)")) +
  theme_ord +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
#' 
#' #### Unified figure
#+ fig2_patchwork,warning=FALSE
fig2_ls <- (its_rich_fig / plot_spacer() / plfa_fig) +
  plot_layout(heights = c(1,0.01,1)) 
fig2 <- (fig2_ls | plot_spacer() | its_ma_ord) +
  plot_layout(widths = c(0.35, 0.01, 0.64)) +
  plot_annotation(tag_levels = 'A') 
#+ fig2,warning=FALSE,fig.height=4,fig.width=6.5
fig2
#' **Fig 2.** Whole-soil fungal communities in **corn**, **restored**, and **remnant** prairie fields.
#' **a** OTU richness and **b** fungal biomass (PLFA) are shown as columns with 95 % CIs; lowercase
#' letters mark significant pairwise differences.
#' **c** Principal-coordinate (PCoA) ordination of ITS-based (97 % OTU) community
#' distances: small points = sites, large circles = field-type centroids (error bars =
#' 95 % CI). Cornfields cluster apart from restored or remnant prairies.
#' Numbers in black circles give years since restoration. Axis labels show the
#' percent variation explained. 

#+ fig2_save,warning=FALSE,echo=FALSE
ggsave(root_path("figs", "fig2.png"),
       plot = fig2,
       width = 6.5,
       height = 4,
       units = "in",
       dpi = 600)
#' 
#' ### ITS, sequence-based relative abundance
#' Comparison figure and stats for supplemental
#+ its_ord
d_its <- its_avg %>% 
  data.frame(row.names = 1) %>%
  decostand("total") %>%
  vegdist("bray")
mva_its <- mva(d = d_its, env = sites)
#+ its_ord_results
mva_its$dispersion_test
mva_its$permanova
mva_its$pairwise_contrasts[c(1,3,2), c(1,2,4,3,8)] %>% 
    kable(format = "pandoc", caption = "Pairwise permanova contrasts")
#' No eignevalue correction was needed. Two relative eigenvalues exceeded broken stick model. 
#' Based on the homogeneity of variance test, the null hypothesis of equal variance among groups is 
#' accepted across all clusters and in pairwise comparison of clusters (both p>0.05), supporting the application of 
#' a PERMANOVA test. 
#' 
#' Clustering revealed that community variation was related to geographic distance, the covariate in 
#' the model. With geographic distance accounted for, the test variable 'field type' significantly explained 
#' variation in fungal communities, with a post-hoc test revealing that communities in corn fields differed from
#' communities in restored and remnant fields. 
#' 
#' Plotting results: 
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
        x = paste0("Axis 1 (", mva_its$axis_pct[1], "%)"),
        y = paste0("Axis 2 (", mva_its$axis_pct[2], "%)")) +
  scale_fill_manual(values = ft_pal) +
  theme_ord +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
#' 
#' #### Supplemental figure
#+ its_shan_ord_sup_patchwork,warning=FALSE
its_shan_ord_sup <- (its_shan_fig | plot_spacer() | its_ord) +
  plot_layout(widths = c(0.45, 0.01, 0.55)) +
  plot_annotation(tag_levels = 'A') 
#+ its_shan_ord_sup,warning=FALSE,fig.height=4,fig.width=6.5
its_shan_ord_sup
#+ its_shan_ord_sup_save,warning=FALSE,echo=FALSE
ggsave(root_path("figs", "figS3.png"),
       plot = its_shan_ord_sup,
       width = 7.5,
       height = 4,
       units = "in",
       dpi = 600)
#' 
#' #### Contrast community metrics
#' Procrustes test on PCoA values using axes with eigenvalues exceeding a broken stick model
#+ its_protest
set.seed(20251111)
its_protest <- protest(
  pcoa(d_its)$vectors[, 1:2],
  pcoa(d_its_ma)$vectors[, 1:2],
  permutations = 1999
)
its_protest
#' Including biomass changes little. The spatial configurations of both ordinations are highly correlated.
#' $R^{2}=$ `r round(its_protest$scale^2, 2)`, p<0.001. 
#' 
#' ## Constrained Analysis
#' Test explanatory variables for correlation with site ordination. Using plant data, 
#' so the analysis is restricted to Wisconsin sites. Edaphic variables are too numerous 
#' to include individually, so transform micro nutrients using PCA. Forb and grass 
#' cover is highly collinear; use the grass-forb index produced previously with PCA. 
soil_micro_pca <- 
  soil %>% 
  left_join(sites %>% select(field_name, field_type, region), by = join_by(field_name)) %>% 
  filter(field_type == "restored", !(region %in% "FL")) %>% 
  select(field_name, SO4, Zn, Fe, Mn, Cu, Ca, Mg, Na, -field_key, -field_type, -region) %>% 
  column_to_rownames(var = "field_name") %>% 
  decostand(method = "standardize") %>% 
  rda()
summary(soil_micro_pca) # 70% on first two axes
soil_micro_index <- scores(soil_micro_pca, choices = c(1, 2), display = "sites") %>% 
  data.frame() %>% 
  rename(soil_micro_1 = PC1, soil_micro_2 = PC2) %>% 
  rownames_to_column(var = "field_name")
soil_macro <- 
  soil %>% 
  left_join(sites %>% select(field_name, field_type, region), by = join_by(field_name)) %>% 
  filter(field_type == "restored", !(region %in% "FL")) %>% 
  select(-c(field_key, field_type, region, SO4, Zn, Fe, Mn, Cu, Ca, Mg, Na))
#' 
#' Assemble explanatory variables and begin iterative selection process. 
#' Plant functional groups and traits not included here were eliminated in previous forward selection
#' procedures (not shown). 
#' Check the VIF for each explanatory variable to test for collinearity if model overfitting is 
#' detected. Then run forward selection in `dbrda()`. 
#' 
env_vars <- sites %>% 
  filter(field_type == "restored", !(region %in% "FL")) %>% 
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
#' Check VIF
env_expl %>% 
  cor() %>% 
  solve() %>% 
  diag() %>% 
  sort()
#' OM, K, and soil_micro_1 with high VIF in initial VIF check. 
#' Removed soil_micro_1 and K to maintain OM in the model.
#' No overfitting detected in full model; proceed with forward selection. 
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
#' 
#' ### Constrained Analysis Results
mod_step
(mod_glax <- anova(mod_step, permutations = 1999))
(mod_inax <- anova(mod_step, by = "axis", permutations = 1999))
(mod_r2   <- RsquareAdj(mod_step, permutations = 1999))
anova(mod_step, by = "margin", permutations = 1999) %>% 
  as.data.frame() %>% 
  mutate(p.adj = p.adjust(`Pr(>F)`, "fdr")) %>% 
  kable(, format = "pandoc")
#' Based on permutation tests with n=1999 permutations, the model shows a significant 
#' correlation between the site ordination on fungal communities
#' and the selected explanatory variables (p<0.001). The first two constrained axes are 
#' also significant (p<0.001, P<0.01). The selected variables explain $R^{2}_{\text{Adj}}$=21.3% of the community 
#' variation. Selected explanatory variables are pH and the grass-forb index; see table for 
#' individual p values and statistics. 
#' 
#' Create the figure, combine with pathogen-plant correlation figure in patchwork later:
mod_pars <- 
  dbrda(
    spe_its_wi_resto ~ gf_index + pH + Condition(env_cov),
    data = env_expl,
    distance = "bray"
  )
mod_pars_eig <- round(mod_pars$CCA$eig * 100, 1)
mod_scor <- scores(
  mod_pars,
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
    mutate(envlabs = c(">forb", "pH")),
  data.frame(
    envvar = "gf_index",
    dbRDA1 = -0.847304624873555,
    dbRDA2 = -0.448178789544163,
    envlabs = ">grass")
  ) %>% 
  arrange(envvar, envlabs) %>% 
  mutate(
    origin = 0,
    m = dbRDA2 / dbRDA1, 
    d = sqrt(dbRDA1^2 + dbRDA2^2), 
    dadd = sqrt((max(dbRDA1)-min(dbRDA2))^2 + (max(dbRDA2)-min(dbRDA2))^2)*0.05,
    labx = ((d+dadd)*cos(atan(m)))*(dbRDA1/abs(dbRDA1)), 
    laby = ((d+dadd)*sin(atan(m)))*(dbRDA1/abs(dbRDA1)))
#+ fig6a,warning=FALSE,fig.height=4,fig.width=4
fig6a <- 
  ggplot(mod_scor_site, aes(x = dbRDA1, y = dbRDA2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  geom_segment(data = mod_scor_bp, 
               aes(x = origin, xend = dbRDA1, y = origin, yend = dbRDA2), 
               arrow = arrow(length = unit(2, "mm"), type = "closed"),
               color = c("darkblue", "darkblue", "gray20")) +
  geom_text(data = mod_scor_bp, 
            aes(x = labx, y = laby, label = envlabs), 
            nudge_x = c(-0.1, 0.1, 0), nudge_y = c(0.06, -0.06, 0),
            size = 3, color = "black") +
  geom_point(fill = ft_pal[2], size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  labs(
    x = paste0("Constr. Axis 1 (", mod_pars_eig[1], "%)"),
    y = paste0("Constr. Axis 2 (", mod_pars_eig[2], "%)")) +
  theme_ord +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
#' Will combine with 6b (AMF)
#' 
#' ## Fungi‒pfg correlations
fungi_resto <- its_div %>% 
  left_join(fa %>% select(field_name, fungi_mass = fungi_18.2), by = join_by(field_name)) %>% 
  left_join(sites, by = join_by(field_name, field_type)) %>% 
  left_join(gf_index, by = join_by(field_name)) %>% 
  filter(field_type == "restored", region != "FL") %>% 
  select(field_name, fungi_ab = depth, fungi_mass, gf_index)

fuma_rest_m <- lm(fungi_mass ~ gf_index, data = fungi_resto)
#+ cm1,warning=FALSE,fig.width=7,fig.height=9
check_model(fuma_rest_m)
summary(fuma_rest_m)
#' PFG doesn't strongly predict fungal biomass at sites. How are mass and sequence 
#' abundance related?
furest_m_raw  <- lm(fungi_ab ~ fungi_mass + gf_index, data = fungi_resto)
furest_m_logy <- lm(log(fungi_ab) ~ fungi_mass + gf_index, data = fungi_resto)
furest_m_logx <- lm(fungi_ab ~ log(fungi_mass) + gf_index, data = fungi_resto)
furest_m_both <- lm(log(fungi_ab) ~ log(fungi_mass) + gf_index, data = fungi_resto)
#' 
compare_performance(furest_m_raw, furest_m_logy, furest_m_logx, furest_m_both,
                    metrics = c("AIC", "RMSE","R2"), rank = TRUE)
#+ cm2,warning=FALSE,fig.width=7,fig.height=9
check_model(furest_m_logy)
summary(furest_m_logy)
ggplot(fungi_resto, aes(x = gf_index, y = fungi_mass)) +
  geom_text(label = rownames(fungi_resto))
#' No relationshp detected with gf_index. High leverage point is again KORP1.
#' ITS based community turnover is correlated with grass/forb change, but abundances
#' of grasses/forbs aren't correlated with these fungi.

#' 
#' # Arbuscular mycorrhizal fungi
# AMF ———————— ####
#' ## Diversity Indices
#+ amf_diversity
amf_div <- calc_div(amf_avg, sites)
#' 
#' ### Richness
#' Account for sequencing depth as a covariate. Test covar transformations 
#' for best model performance.
amf_rich_covar <- covar_shape_test(
  data  = amf_div,
  y     = "richness",       
  covar = "depth",   
  group = "field_type"           
)
amf_rich_covar$compare
#' Log transformation performs best, but very little difference among models
amf_rich_lm <- lm(richness ~ depth_clg + field_type, 
                  data = amf_div %>% mutate(depth_clg = log(depth) - mean(log(depth))))
#' Diagnostics
#+ amf_rich_covar_diagnostics,warning=FALSE,fig.width=7,fig.height=9
amf_rich_covar$diagnostics
#' Long tails, one outlier without significant leverage...mean/variance relationship shows no trend...
distribution_prob(amf_rich_lm)
#' Residuals distribution most likely normal, response bimodal (ignore)
leveneTest(richness ~ field_type, data = amf_div) %>% as.data.frame() %>% kable(format = "pandoc")
leveneTest(residuals(amf_rich_lm) ~ amf_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc")
#' Residuals/response distributions do not suggest the need for transformation.
#' Levene's p > 0.05 → fail to reject = variances can be considered equal.
#' 
#' Model results, group means, and post-hoc
amf_rich_covar$anova_t2
#' Sequencing depth not a significant predictor of amf richness
amf_rich_em <- emmeans(amf_rich_lm, ~ field_type, type = "response")
#' Results tables below show the emmeans summary of group arithmetic means and confidence intervals, 
#' and the post hoc contrast of richness among field types. 
#' Main effect in model significant after p value adjustment (see summary section): pairwise
#' contrast warranted.
#+ amf_rich_em_summary,echo=FALSE
kable(summary(amf_rich_em), 
      format = "pandoc", 
      caption = "Confidence level used: 0.95")
#+ amf_rich_em_posthoc,echo=FALSE
kable(pairs(amf_rich_em), 
      format = "pandoc", 
      caption = "P value adjustment: tukey method for comparing a family of 3 estimates")
#' OTU richness in cornfields is significantly less than in restored or remnant fields, which 
#' don't differ. Plot the results
#+ amf_richness_fig,fig.width=4,fig.height=4,fig.align='center'
amf_rich_fig <- 
  ggplot(summary(amf_rich_em), aes(x = field_type, y = emmean)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = emmean, ymax = upper.CL), width = 0, linewidth = lw) +
  geom_text(aes(y = upper.CL, label = c("a", "b", "b")), vjust = -1.5, family = "sans", size = 4) +
  labs(x = NULL, y = "Richness") +
  lims(y = c(0, 75)) +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1))
#' 
#' ### Shannon diversity
#' Account for sequencing depth as a covariate. Test covar transformations for best model performance.
amf_shan_covar <- covar_shape_test(
  data  = amf_div,
  y     = "shannon",       
  covar = "depth",   
  group = "field_type"           
)
amf_shan_covar$compare
#' Models are equivalent; no transformation selected on parsimony grounds
amf_shan_lm <- lm(shannon ~ depth + field_type, data = amf_div)
#' Diagnostics
#+ amf_shan_covar_diagnostics,warning=FALSE,fig.width=7,fig.height=9
amf_shan_covar$diagnostics
par(mfrow = c(2,2))
plot(amf_shan_lm)
#' Variance somewhat non-constant in groups, qqplot fit is poor, 
#' one leverage point (Cook's > 0.5), a cornfield with high richness. Mean
#' richness in corn fields is lowest; this outlier would make the pairwise contrast
#' less significant, possible Type II error which is more acceptable.
distribution_prob(amf_shan_lm)
#' Residuals/response distributions most likely normal. 
leveneTest(shannon ~ field_type, data = amf_div) %>% as.data.frame() %>% kable(format = "pandoc")
leveneTest(residuals(amf_shan_lm) ~ amf_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc")
#' Residuals/response distributions do not suggest the need for transformation.
#' Covariate adds little added explanatory value.
#' Levene's p > 0.05 → fail to reject = variances can be considered equal.
#' 
#' Model results, group means, and post-hoc
Anova(amf_shan_lm, type = 2)
#' Sequencing depth not a significant predictor of Shannon diversity. Produce arithmetic means
#' in groups and post hoc contrasts
amf_shan_em <- emmeans(amf_shan_lm, ~ field_type, type = "response")
#' Results tables below show the emmeans summary of group means and confidence intervals,
#' with sequencing depth as a covariate, and the post hoc contrast of richness among field types. 
#+ amf_shan_em_summary,echo=FALSE
kable(summary(amf_shan_em), 
      format = "pandoc", 
      caption = "Confidence level used: 0.95")
#+ amf_shan_em_posthoc,echo=FALSE
kable(pairs(amf_shan_em), 
      format = "pandoc", 
      caption = "P value adjustment: tukey method for comparing a family of 3 estimates")
#' Shannon's diversity in cornfields is significantly less than in restored or remnant fields, which 
#' don't differ.
#+ amf_shannons_fig,fig.width=4,fig.height=4,fig.align='center'
amf_shan_fig <- 
  ggplot(summary(amf_shan_em), aes(x = field_type, y = emmean)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = emmean, ymax = upper.CL), width = 0, linewidth = lw) +
  geom_text(aes(y = upper.CL, label = c("a", "b", "b")), vjust = -1.5, family = "sans", size = 4) +
  labs(x = NULL, y = "Shannon diversity") +
  lims(y = c(0, 32)) +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
#' 
#' ## Abundance (NLFA biomass)
nlfa_lm <- lm(amf ~ field_type, data = fa)
#' Diagnostics
par(mfrow = c(2,2))
plot(nlfa_lm) # variance obviously not constant in groups
distribution_prob(nlfa_lm)
# response distribution gamma; resids likely normal
leveneTest(residuals(nlfa_lm) ~ fa$field_type) # No covariate, response and residuals tests equivalent
#' Residuals distribution variance may not be equal in groups.
#' Levene's p = 0.054, close to rejecting the null of equal variance. Check CV in groups.
fa %>%
  mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant"))) %>%
  group_by(field_type) %>%
  summarize(mean = mean(amf),
            cv = sd(amf) / mean) %>%
  mutate(across(mean:cv, ~ round(.x, 2))) %>%
  kable(format = "pandoc", caption = "Mean and CV relationship in groups")
#' CV increases with mean, suggesting > proportional mean/variance relationship. 
#' Determine best model choice of log-transformed response or gamma glm.
nlfa_lm_log   <- lm(log(amf) ~ field_type, data = fa)
par(mfrow = c(2,2))
plot(nlfa_lm_log) # qqplot ok, one high leverage point in remnants
ncvTest(nlfa_lm_log) # p=0.16, null of constant variance not rejected

nlfa_glm  <- glm(amf ~ field_type, family = Gamma(link = "log"), data = fa)
nlfa_glm_diag <- glm.diag(nlfa_glm)
glm.diag.plots(nlfa_glm, nlfa_glm_diag) # qqplot shows strong fit; no leverage >0.5
performance::check_overdispersion(nlfa_glm) # not detected
#' Gamma glm is the best choice; no high-leverage point
#' 
#' Model results, group means, and post-hoc
summary(nlfa_glm)
anova(nlfa_glm) # Decline in residual deviance worth the cost in df
nlfa_em <- emmeans(nlfa_glm, ~ field_type, type = "response")
#+ nlfa_em_summary,echo=FALSE
kable(summary(nlfa_em),
      format = "pandoc",
      caption = "Confidence level used: 0.95.\nIntervals are back-transformed from the log scale")
#+ nlfa_em_posthoc,echo=FALSE
kable(pairs(nlfa_em),
      format = "pandoc",
      caption = "P value adjustment: tukey method for comparing a family of 3 estimates.\nTests are performed on the log scale")
#+ nlfa_fig,fig.width=4,fig.height=4,
nlfa_fig <-
  ggplot(summary(nlfa_em), aes(x = field_type, y = response)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = response, ymax = upper.CL), width = 0, linewidth = lw) +
  geom_text(aes(y = upper.CL, label = c("a", "b", "b")), vjust = -1.5, family = "sans", size = 4) +
  labs(x = NULL, y = Biomass~(nmol[NLFA]%*%g[soil]^-1)) +
  scale_fill_manual(values = ft_pal) +
  lims(y = c(0, 75)) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1.1))
#' 
#' ## Beta Diversity
#' 
#' 1. Using biomass-weighted relative abundance and b-c distance (unifrac is scale invariant; it's based on the proportion 
#' of reads on each descending branch, multiplying rows by any constant doesn't change this). 
#' 1. Using sequence-based relative abundance and unifrac distance to display phylogenetically aware information.
#' 1. Contrast the two with procrustes.
#' 
#' ### AMF, biomass-weighted relative abundance
#+ amf_ord_ma
d_amf_ma <- amf_avg_ma %>% 
  data.frame(row.names = 1) %>% 
  vegdist("bray")
mva_amf_ma <- mva(d = d_amf_ma, env = sites, corr = "lingoes")
#+ amf_ord_ma_results
mva_amf_ma$dispersion_test
mva_amf_ma$permanova
mva_amf_ma$pairwise_contrasts[c(1,3,2), c(1,2,4,3,8)] %>% 
  kable(format = "pandoc", caption = "Pairwise permanova contrasts")
#' Lingoes correction was applied to negative eignevalues. Three relative eigenvalues exceeded broken stick model. 
#' Based on the homogeneity of variance test, the null hypothesis of equal variance among groups is 
#' accepted across all clusters and in pairwise comparison of clusters (both p>0.05), supporting the application of 
#' a PERMANOVA test. 
#' 
#' Clustering revealed that community variation was not related to geographic distance, the covariate in 
#' the model. With geographic distance accounted for, the test variable 'field type' significantly explained 
#' variation in fungal communities, with a post-hoc test revealing that communities in corn fields differed from
#' communities in restored and remnant fields. 
#' 
#' Plotting results: 
amf_ma_ord_data <- mva_amf_ma$ordination_scores %>% mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant")))
p_amf_ma_centers <- amf_ma_ord_data %>% 
  group_by(field_type) %>% 
  summarize(across(starts_with("Axis"), list(mean = mean, ci_l = ci_l, ci_u = ci_u), .names = "{.fn}_{.col}"), .groups = "drop") %>% 
  mutate(across(c(ci_l_Axis.1, ci_u_Axis.1), ~ mean_Axis.1 + .x),
         across(c(ci_l_Axis.2, ci_u_Axis.2), ~ mean_Axis.2 + .x),
         across(ends_with("Axis.1"), ~ .x * -1)) # reverse axis values to be consistent with other plots
amf_ma_ord <- 
  ggplot(amf_ma_ord_data, aes(x = Axis.1 * -1, y = Axis.2)) +  # reverse axis values to be consistent with other plots
  geom_linerange(data = p_amf_ma_centers, aes(x = mean_Axis.1, y = mean_Axis.2, xmin = ci_l_Axis.1, xmax = ci_u_Axis.1), linewidth = lw) +
  geom_linerange(data = p_amf_ma_centers, aes(x = mean_Axis.1, y = mean_Axis.2, ymin = ci_l_Axis.2, ymax = ci_u_Axis.2), linewidth = lw) +
  geom_point(data = p_amf_ma_centers, aes(x = mean_Axis.1, y = mean_Axis.2, fill = field_type), size = lg_size, stroke = lw, shape = 21) +
  geom_point(aes(fill = field_type), size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  scale_fill_manual(values = ft_pal) +
  labs(
    x = paste0("Axis 1 (", mva_amf_ma$axis_pct[1], "%)"),
    y = paste0("Axis 2 (", mva_amf_ma$axis_pct[2], "%)")) +
  theme_ord +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
#' 
#' #### Unified figure
#+ fig3_patchwork,warning=FALSE
fig3_ls <- (amf_rich_fig / plot_spacer() / nlfa_fig) +
  plot_layout(heights = c(1,0.01,1)) 
fig3 <- (fig3_ls | plot_spacer() | amf_ma_ord) +
  plot_layout(widths = c(0.35, 0.01, 0.64)) +
  plot_annotation(tag_levels = 'A') 
#+ fig3,warning=FALSE,fig.height=4,fig.width=6.5
fig3
#' **Fig 3.** AMF communities in corn, restored, and remnant prairie fields. OTU richness **a**;
#' NLFA biomass **b** (95 % CI, letters = Tukey groups). PCoA of BC distances on proportion of biomass abundance data
#' (18S, 97 % OTUs) **c**: small points = sites, large rings = field‑type centroids ±95 % CI. 
#' Numbers in points give years since restoration. Axes show % variance. Corn clusters apart from both 
#' prairie types. Shading: corn grey, restored black, remnant white.

#+ fig3_save,warning=FALSE,fig.height=5,fig.width=7,echo=FALSE
ggsave(root_path("figs", "fig3.png"),
       plot = fig3,
       width = 6.5,
       height = 4,
       units = "in",
       dpi = 600)
#' 
#' ### AMF, sequence-based relative abundance, unifrac distance
#' 
#+ amf_ord
d_amf <- UniFrac(amf_ps, weighted = TRUE, normalized = TRUE)
mva_amf <- mva(d = d_amf, env = sites, corr = "lingoes")
#+ amf_ord_results
mva_amf$dispersion_test
mva_amf$permanova
mva_amf$pairwise_contrasts[c(1,3,2), c(1,2,4,3,8)] %>% 
  kable(format = "pandoc", caption = "Pairwise permanova contrasts")
#' Lingoes eigenvalue correction was used. The first three relative eigenvalues exceeded broken stick model. 
#' Based on the homogeneity of variance test, the null hypothesis of equal variance among groups is 
#' accepted across all clusters and in pairwise comparison of clusters (both p>0.05), supporting the application of 
#' a PERMANOVA test. 
#' 
#' Clustering revealed that geographic distance among sites did not significantly explain AMF community variation. 
#' With geographic distance accounted for, the test variable field type significantly explained 
#' variation in AMF communities, with a post-hoc test revealing that communities in corn fields differed from
#' communities in restored and remnant fields. 
#' 
#' Plotting the result:
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
    x = paste0("Axis 1 (", mva_amf$axis_pct[1], "%)"),
    y = paste0("Axis 2 (", mva_amf$axis_pct[2], "%)")) +
  theme_ord +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
#' 
#' #### Supplemental figure
#+ amf_shan_ord_sup_patchwork,warning=FALSE
amf_shan_ord_sup <- (amf_shan_fig | plot_spacer() | amf_ord) +
  plot_layout(widths = c(0.45, 0.01, 0.55)) +
  plot_annotation(tag_levels = 'A') 
#+ amf_shan_ord_sup,warning=FALSE,fig.height=4,fig.width=6.5
amf_shan_ord_sup
#+ amf_shan_ord_sup_save,warning=FALSE,echo=FALSE
ggsave(root_path("figs", "figS5.png"),
       plot = amf_shan_ord_sup,
       width = 7.5,
       height = 4,
       units = "in",
       dpi = 600)
#' 
#' #### Contrast community metrics
#' Procrustes test on PCoA values using axes with eigenvalues exceeding a broken stick model
#+ amf_protest
set.seed(20251111)
amf_protest <- protest(
  pcoa(d_amf, correction = "lingoes")$vectors[, 1:3],
  pcoa(d_amf_ma, correction = "lingoes")$vectors[, 1:3],
  permutations = 1999
)
amf_protest
#' The ordinations differ in spatial arrangement somewhat, with a correlation of
#' $R^{2}=$ `r round(amf_protest$scale^2, 2)`, however, the null that these solutions are unrelated
#' is still rejected at p<0.001. Clearly, the low biomass in cornfields is a driving difference in 
#' the biomass-aware ordination, which, as a result, should possibly be preferred in this case.
#' 
#' ## AMF abundance in families
#' Test proportion of biomass across field types for each family. A similar contrast with 
#' sequence abundance only returns largely NS results, highlighting the incorrect inference 
#' when using compositional data to conduct site-differential analysis.
#' 
#' For each AMF family, fit a Gamma GLM with log link predicting biomass-weighted relative abundance by 
#' field type. We controlled the false discovery rate across the four family-level tests using a
#' Benjamini–Hochberg (BH) false discovery rate correction. For pairwise contrasts, within models 
#' showing a significant field-type effect, compare estimated marginal means with Tukey adjustment 
#' for all pairwise contrasts (emmeans). 
#' 
#' ### Glomeraceae
glom_lm <- lm(Glmrc_mass ~ field_type, data = amf_fam_ma)
par(mfrow = c(2,2))
plot(glom_lm) # variance non-constant among groups
distribution_prob(glom_lm)
#' residuals normal, response chi/gamma...mean variance relationship suggested
leveneTest(residuals(glom_lm) ~ amf_fam_ma$field_type)
#' Levene's p<0.05, null of equal variance rejected
amf_fam_ma %>% 
  mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant"))) %>%
  group_by(field_type) %>%
  summarize(mean = mean(Glmrc_mass),
            cv = sd(Glmrc_mass) / mean) %>%
  mutate(across(mean:cv, ~ round(.x, 2))) %>%
  kable(format = "pandoc", caption = "Mean and CV relationship in groups")
#' Corroborates the mean/variance relationship. Proceed with log transformation or gamma glm
glom_lm_log   <- lm(log(Glmrc_mass) ~ field_type, data = amf_fam_ma)
par(mfrow = c(2,2))
plot(glom_lm_log) # some improvement, leverage point still exists
ncvTest(glom_lm_log) # p=0.29, null of constant variance not rejected, model fit is improved
#' Gamma glm to reduce leverage point...
glom_glm <- glm(Glmrc_mass ~ field_type, family = Gamma(link = "log"), data = amf_fam_ma)
#' Diagnostics
glom_glm_diag <- glm.diag(glom_glm)
glm.diag.plots(glom_glm, glom_glm_diag) # qqplot shows strong fit; no leverage >0.5
#+ cm3,warning=FALSE,fig.width=7,fig.height=9
check_model(glom_glm) # corroborates
performance::check_overdispersion(glom_glm) # not detected
#' Gamma glm is the best choice; no high-leverage point
#' 
#' Model results, group means, and post-hoc
anova(glom_glm) # Decline in residual deviance worth the cost in df
glom_em <- emmeans(glom_glm, ~ field_type, type = "response")
#' 
#' ### Claroideoglomeraceae
clar_lm <- lm(Clrdg_mass ~ field_type, data = amf_fam_ma)
par(mfrow = c(2,2))
plot(clar_lm) # variance non-constant among groups
distribution_prob(clar_lm)
#' residuals and response very long tailed
leveneTest(residuals(clar_lm) ~ amf_fam_ma$field_type)
#' Levene's p>0.05, null of equal variance not rejected
#' Look to gamma to reduce outlier and remain consistent with this model set
clar_glm <- glm(Clrdg_mass ~ field_type, family = Gamma(link = "log"), data = amf_fam_ma)
#' Diagnostics
clar_glm_diag <- glm.diag(clar_glm)
glm.diag.plots(clar_glm, clar_glm_diag) # qqplot shows strong fit; no outlier point
#+ cm4,warning=FALSE,fig.width=7,fig.height=9
check_model(clar_glm) # corroborates
performance::check_overdispersion(clar_glm) # not detected
#' Gamma glm is the best choice; no high-leverage point
#' 
#' Model results, group means, and post-hoc
anova(clar_glm) # Decline in residual deviance worth the cost in df
clar_em <- emmeans(clar_glm, ~ field_type, type = "response")
#' 
#' ### Paraglomeraceae
para_lm <- lm(Clrdg_mass ~ field_type, data = amf_fam_ma)
par(mfrow = c(2,2))
plot(para_lm) # variance non-constant among groups
distribution_prob(para_lm)
#' residuals and response very long tailed
leveneTest(residuals(para_lm) ~ amf_fam_ma$field_type)
#' Levene's p>0.05, null of equal variance not rejected
#' Look to gamma to reduce outlier and remain consistent with this model set
para_glm <- glm(Prglm_mass ~ field_type, family = Gamma(link = "log"), data = amf_fam_ma)
#' Diagnostics
para_glm_diag <- glm.diag(para_glm)
glm.diag.plots(para_glm, para_glm_diag) # qqplot shows strong fit; no outlier point
#+ cm5,warning=FALSE,fig.width=7,fig.height=9
check_model(para_glm) # corroborates
performance::check_overdispersion(para_glm) # not detected
#' Gamma glm is the best choice; no high-leverage point
#' 
#' Model results, group means, and post-hoc
anova(para_glm) # field_type not significant
para_em <- emmeans(para_glm, ~ field_type, type = "response")
#' Pairwise significance not appropriate due to overall variable NS
#' 
#' ### Diversisporaceae
diver_glm <- glm(Dvrss_mass ~ field_type, family = Gamma(link = "log"), data = amf_fam_ma)
#' Diagnostics
diver_glm_diag <- glm.diag(diver_glm)
glm.diag.plots(diver_glm, diver_glm_diag) # qqplot shows strong fit; no outlier point
#+ cm6,warning=FALSE,fig.width=7,fig.height=9
check_model(diver_glm) # corroborates
performance::check_overdispersion(diver_glm) # not detected
#' Gamma glm is the best choice; no high-leverage point
#' Model results, group means, and post-hoc
anova(diver_glm) # Decline in residual deviance worth the cost in df
diver_em <- emmeans(diver_glm, ~ field_type, type = "response")
#' 
#' ### Gigasporaceae and others
#' Zeroes in data; comparison not warranted
#' 
#' ### Results table
#' Combine later with confidence intervals and corrected p values
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
#' Pairwise contrasts
#' Paraglomeraceae model was NS; no pairwise indicated
list(glom = glom_em, clar = clar_em, diver = diver_em) %>% 
  map(\(df) pairs(df, adjust = "tukey"))
#' Model results
list(glom = glom_glm, clar = clar_glm, para = para_glm, diver = diver_glm) %>% 
  map(\(df) df %>% anova() %>% as.data.frame()) %>% 
  bind_rows(.id = "family") %>% 
  mutate(p.adj = p.adjust(`Pr(>F)`, "fdr") %>% round(., 5))
#' 
#' ## Indicator species analysis
amf_ind <- inspan(spe=amf_avg_ma, meta=amf_meta, guild=NULL, site_dat=sites)
amf_ind %>% 
  select(field_type, family, taxon, starts_with("corn"), starts_with("restor"), starts_with("rem")) %>% 
  filter(taxon != "unidentified") %>% 
  arrange(field_type) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% 
  kable(format = "pandoc", caption = "Indicator species analysis results with avg biomass")
#' None of the identified species seem relevant for further discussion...
#' 
#' ## Constrained Analysis
#' Env covars processed in the ITS section (see above)
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
#' 
#' ### Constrained Analysis Results
amf_mod_step
(amf_mod_glax <- anova(amf_mod_step, permutations = 1999))
(amf_mod_inax <- anova(amf_mod_step, by = "axis", permutations = 1999))
(amf_mod_r2   <- RsquareAdj(amf_mod_step, permutations = 1999))
amf_mod_step$anova %>% kable(, format = "pandoc")
#' Based on permutation tests with n=1999 permutations, 
#' after accounting for inter-site pairwise distance as a covariate, the model shows a significant
#' correlation between the site ordination on fungal communities
#' and the selected explanatory variables (p<0.001). The first two constrained axes are
#' also significant (p<0.001, p<0.02). The selected variables explain $R^{2}_{\text{Adj}}$=35.1% of the community
#' variation. Selected explanatory variables are pH and the grass-forb index; see table for
#' individual p values and statistics.
#' 
#' ### AMF constrained figure
amf_mod_pars <-
  dbrda(
    spe_amf_wi_resto ~ gf_index + OM + Condition(env_cov),
    data = env_expl,
    distance = "bray"
  )
amf_mod_pars_eig <- round(amf_mod_pars$CCA$eig * 100, 1)

amf_mod_scor <- scores(
  amf_mod_pars,
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
    mutate(envlabs = c(">forb", "OM")),
  data.frame(
    envvar = "gf_index",
    dbRDA1 = 0.8450446,
    dbRDA2 = -0.4524256,
    envlabs = ">grass")
) %>% 
  arrange(envvar, envlabs) %>% 
  mutate(
    origin = 0,
    m = dbRDA2 / dbRDA1,
    d = sqrt(dbRDA1^2 + dbRDA2^2),
    dadd = sqrt((max(dbRDA1)-min(dbRDA2))^2 + (max(dbRDA2)-min(dbRDA2))^2)*0.1,
    labx = ((d+dadd)*cos(atan(m)))*(dbRDA1/abs(dbRDA1)),
    laby = ((d+dadd)*sin(atan(m)))*(dbRDA1/abs(dbRDA1)))
#+ fig6b,warning=FALSE,fig.height=4,fig.width=4
fig6b <-
  ggplot(amf_mod_scor_site, aes(x = (dbRDA1 * -1), y = dbRDA2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  geom_segment(data = amf_mod_scor_bp,
               aes(x = origin, xend = (dbRDA1 * -1), y = origin, yend = dbRDA2),
               arrow = arrow(length = unit(2, "mm"), type = "closed"),
               color = c("gray20", "darkblue", "darkblue")) +
  geom_text(data = amf_mod_scor_bp,
            aes(x = (labx * -1), y = laby, label = envlabs),
            nudge_x = (c(0.05, 0.2, -0.2) * -1), nudge_y = c(0.1, 0.04, -0.04),
            size = 3, color = "gray20") +
  geom_point(fill = ft_pal[2], size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  labs(
    x = paste0("Constr. Axis 1 (", amf_mod_pars_eig[1], "%)"),
    y = paste0("Constr. Axis 2 (", amf_mod_pars_eig[2], "%)")) +
  lims(x = c(-1.1,1.05), y = c(-1.6,0.9)) +
  theme_ord +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
#' 
#' #### Unified figure
#' Display results of constrained analysis for ITS and AMF
#+ fig6_patchwork,warning=FALSE
fig6 <- (fig6a | plot_spacer() | fig6b) +
  plot_layout(widths = c(0.50, 0.01, 0.50)) +
  plot_annotation(tag_levels = 'A') 
#+ fig6,warning=FALSE,fig.height=4,fig.width=6.5
fig6
#' Results of constrained analysis on **a** whole-soil fungi and **b** arbuscular mycorrhizal fungi.
#' Percent variation explained by each constrained axis is shown with axis labels. Points show locations
#' of restored fields in Wisconsin based on fungal community distance. Blue arrows show the grass-forb 
#' index with labels indicating the direction of relative increase in grasses or forbs, respectively,
#' along the index. The black arrows show significant constraints from edaphic properties. 
#+ fig6_save,warning=FALSE,echo=FALSE
ggsave(
  root_path("figs", "fig6.png"),
  plot = fig6,
  width = 6.5,
  height = 3.5,
  units = "in",
  dpi = 600
)
#' 
#' ## AMF correlations with pfg
amf_resto <- amf_fam %>% 
  rowwise() %>% 
  mutate(amf_ab = sum(c_across(Glmrc_ab:Ambsp_ab))) %>% 
  select(-c(Glmrc_ab:Ambsp_ab)) %>% 
  filter(field_type == "restored", region != "FL") %>% 
  left_join(fa %>% select(field_name, amf_mass = amf), by = join_by(field_name))
amma_rest_m <- lm(amf_mass ~ gf_index, data = amf_resto)
#+ cm7,warning=FALSE,fig.width=7,fig.height=9
check_model(amma_rest_m)
summary(amma_rest_m)
#' No simple linear relationship. Check contributions of biomass and sequence abundance.
amrest_m_raw  <- lm(amf_ab ~ amf_mass + gf_index, data = amf_resto)
amrest_m_logy <- lm(log(amf_ab) ~ amf_mass + gf_index, data = amf_resto)
amrest_m_logx <- lm(amf_ab ~ log(amf_mass) + gf_index, data = amf_resto)
amrest_m_both <- lm(log(amf_ab) ~ log(amf_mass) + gf_index, data = amf_resto)
compare_performance(amrest_m_raw, amrest_m_logy, amrest_m_logx, amrest_m_both,
                    metrics = c("AIC", "RMSE","R2"), rank = TRUE)
#+ cm8,warning=FALSE,fig.width=7,fig.height=9
check_model(amrest_m_logx)
summary(amrest_m_logx)
#+ amf_resto_fig,warning=FALSE,fig.height=5,fig.width=5
ggplot(amf_resto, aes(x = gf_index, y = amf_mass)) +
  geom_text(label = rownames(amf_resto))
#' No relationshp detected with gf_index. High leverage point is again KORP1.

#' 
#' # Putative plant pathogens
# Putative plant pathogens ———————— ####
#' 
#' Retrieve pathogen sequence abundance
patho <- guildseq(its_avg, its_meta, "plant_pathogen")
#' 
#' ## Diversity Indices
#+ patho_diversity
patho_div <- calc_div(patho, sites)
#' 
#' ### Richness
#' Account for sequencing depth as a covariate. Compare models with raw, sqrt, and log transformed depth.
patho_rich_covar <- covar_shape_test(
  data  = patho_div,
  y     = "richness",       
  covar = "depth",   
  group = "field_type"           
)
patho_rich_covar$compare
#' Sqrt transform selected based on model performance. Model with centered, sqrt depth covar
patho_rich_lm <- lm(richness ~ depth_csq + field_type, # Interaction NS (not shown)
                    data = patho_div %>% mutate(depth_csq = sqrt(depth) - mean(sqrt(depth))))
#' Diagnostics
#+ patho_rich_covar_diagnostics,warning=FALSE,fig.width=7,fig.height=9
patho_rich_covar$diagnostics
distribution_prob(patho_rich_lm)
#' residuals distribution normal or close, response showing group divisions
leveneTest(richness ~ field_type, data = patho_div) %>% as.data.frame() %>% 
  kable(format = "pandoc", caption = "Response var in groups")
leveneTest(residuals(patho_rich_lm) ~ patho_div$field_type) %>% as.data.frame() %>% 
  kable(format = "pandoc", caption = "Residuals var in groups")
#' Residuals/response distributions do not suggest the need for transformation.
#' Levene's p > 0.05 → fail to reject = variances can be considered equal.
#' 
#' Model results, group means, and post-hoc
patho_rich_covar$anova_t2
#' Sequence depth is highly significant (also after p value adjustment); 
#' richness doesn't vary in groups. 
#' Calculate confidence intervals for figure.
#' Arithmetic means calculated in this case.
patho_rich_em <- emmeans(patho_rich_lm, ~ field_type, type = "response")
#+ patho_rich_em_summary,echo=FALSE
kable(summary(patho_rich_em), 
      format = "pandoc", 
      caption = "Confidence level used: 0.95")
#+ patho_rich_em_posthoc,echo=FALSE
kable(pairs(patho_rich_em), 
      format = "pandoc", 
      caption = "P value adjustment: tukey method for comparing a family of 3 estimates")
#+ patho_richness_fig,fig.width=4,fig.height=4
patho_rich_fig <- 
  ggplot(summary(patho_rich_em), aes(x = field_type, y = emmean)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = emmean, ymax = upper.CL), width = 0, linewidth = lw) +
  labs(x = NULL, y = "Richness") +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1))
#' 
#' ### Shannon's diversity
#' Account for sequencing depth as a covariate. Test transformations of covariate.
patho_shan_covar <- covar_shape_test(
  data  = patho_div,
  y     = "shannon",       
  covar = "depth",   
  group = "field_type"           
)
patho_shan_covar$compare
#' Log transform selected, center the covariate
patho_shan_lm <- lm(shannon ~ depth_clg + field_type, 
                    data = patho_div %>% mutate(depth_clg = log(depth) - mean(log(depth))))
#' Diagnostics
#+ patho_shan_covar_diagnostics,warning=FALSE,fig.width=7,fig.height=9
patho_shan_covar$diagnostics # variance similar in groups 
distribution_prob(patho_shan_lm)
#' residuals distribution most likely cauchy/normal; symmetric but long tails
#' response normal
leveneTest(shannon ~ field_type, data = patho_div) %>% as.data.frame() %>% kable(format = "pandoc")
leveneTest(residuals(patho_shan_lm) ~ patho_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc")
#' Residuals distribution does not suggest the need for further model selection.
#' Levene's p > 0.05 → fail to reject = variances can be considered equal.
#' 
#' Model results, group means, and post-hoc
patho_shan_covar$anova_t2
#' Sequence depth is a significant predictor of Shannon diversity, field type is not.
patho_shan_em <- emmeans(patho_shan_lm, ~ field_type, type = "response")
#' Results tables below show the emmeans summary of group means and confidence intervals,
#' with sequencing depth as a covariate, and the post hoc contrast of richness among field types.
#+ patho_shan_em_summary,echo=FALSE
kable(summary(patho_shan_em), 
      format = "pandoc", 
      caption = "Confidence level used: 0.95")
#+ patho_shan_em_posthoc,echo=FALSE
kable(pairs(patho_shan_em), 
      format = "pandoc", 
      caption = "P value adjustment: tukey method for comparing a family of 3 estimates")
#+ patho_shan_fig,fig.width=4,fig.height=4
patho_shan_fig <- 
  ggplot(summary(patho_shan_em), aes(x = field_type, y = emmean)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = emmean, ymax = upper.CL), width = 0, linewidth = lw) +
  labs(x = NULL, y = "Shannon diversity") +
  # lims(y = c(0, 160)) +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
#' 
#' ## Abundance
patho_ma_lm <- lm(patho_mass ~ field_type, data = its_guild_ma)
par(mfrow = c(2,2))
plot(patho_ma_lm) 
#' no serious violations observed
distribution_prob(patho_ma_lm)
#' Residuals distribution fits normal, response gamma?
leveneTest(residuals(patho_ma_lm) ~ fa$field_type) %>% as.data.frame() %>% kable(format = "pandoc") 
#' No covariate, response and residuals tests equivalent.
#' Residuals distribution does not suggest the need for transformation.
#' Levene's p > 0.05 → fail to reject = variances can be considered equal.
#' 
#' Model results, group means, and post-hoc, with arithmetic means from emmeans
anova(patho_ma_lm)
patho_ma_em <- emmeans(patho_ma_lm, ~ field_type, type = "response")
#+ patho_plfa_em_summary,echo=FALSE
kable(summary(patho_ma_em), 
      format = "pandoc", 
      caption = "Confidence level used: 0.95")
#+ patho_plfa_em_posthoc,echo=FALSE
kable(pairs(patho_ma_em), 
      format = "pandoc", 
      caption = "P value adjustment: tukey method for comparing a family of 3 estimates")
#+ patho_plfa_fig,fig.width=4,fig.height=4,fig.align='center'
patho_ma_fig <- 
  ggplot(summary(patho_ma_em), aes(x = field_type, y = emmean)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = emmean, ymax = upper.CL), width = 0, linewidth = lw) +
  labs(x = NULL, y = expression(Biomass~(nmol[PLFA]%*%g[soil]^-1))) +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1.1))
#' 
#' ## Beta Diversity
#' Community distances handled similarly to previous
#' ### Biomass-weighted relative abundance
patho_ma <- guildseq(its_avg_ma, its_meta, "plant_pathogen")
#+ patho_ma_ord
d_patho_ma <- patho_ma %>% 
  data.frame(row.names = 1) %>% 
  vegdist("bray")
mva_patho_ma <- mva(d = d_patho_ma, env = sites, corr = "lingoes")
#+ patho_ord_results
mva_patho_ma$dispersion_test
mva_patho_ma$permanova
mva_patho_ma$pairwise_contrasts[c(1,3,2), c(1,2,4,3,8)] %>% 
  kable(format = "pandoc", caption = "Pairwise permanova contrasts")
#' Lingoes correction was needed. Three axes were significant based on broken stick test. 
#' Based on the homogeneity of variance test, the null hypothesis of equal variance among groups is 
#' accepted across all clusters and in pairwise comparison of clusters (both p>0.05), supporting the application of 
#' a PERMANOVA test. 
#' 
#' An effect of geographic distance (covariate) on pathogen communities was not supported. 
#' With geographic distance accounted for, the test variable 'field type' significantly explained 
#' variation in fungal communities, with a post-hoc test revealing that communities in corn fields differed from
#' communities in restored and remnant fields. 
#' 
#' Plotting results: 
patho_ma_ord_data <- mva_patho_ma$ordination_scores %>% mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant")))
p_patho_ma_centers <- patho_ma_ord_data %>% 
  group_by(field_type) %>% 
  summarize(across(starts_with("Axis"), list(mean = mean, ci_l = ci_l, ci_u = ci_u), .names = "{.fn}_{.col}"), .groups = "drop") %>% 
  mutate(across(c(ci_l_Axis.1, ci_u_Axis.1), ~ mean_Axis.1 + .x),
         across(c(ci_l_Axis.2, ci_u_Axis.2), ~ mean_Axis.2 + .x))
patho_ma_ord <- 
  ggplot(patho_ma_ord_data, aes(x = Axis.1, y = Axis.2)) +
  geom_linerange(data = p_patho_ma_centers, aes(x = mean_Axis.1, y = mean_Axis.2, xmin = ci_l_Axis.1, xmax = ci_u_Axis.1), linewidth = lw) +
  geom_linerange(data = p_patho_ma_centers, aes(x = mean_Axis.1, y = mean_Axis.2, ymin = ci_l_Axis.2, ymax = ci_u_Axis.2), linewidth = lw) +
  geom_point(data = p_patho_ma_centers, aes(x = mean_Axis.1, y = mean_Axis.2, fill = field_type), size = lg_size, stroke = lw, shape = 21) +
  geom_point(aes(fill = field_type), size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  scale_fill_manual(values = ft_pal) +
  labs(
    x = paste0("Axis 1 (", mva_patho_ma$axis_pct[1], "%)"),
    y = paste0("Axis 2 (", mva_patho_ma$axis_pct[2], "%)")) +
  theme_ord +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1))
#' 
#' ## Unified figure
#+ fig4_patchwork,warning=FALSE
fig4_ls <- (patho_rich_fig / plot_spacer() / patho_ma_fig) +
  plot_layout(heights = c(1,0.01,1)) 
fig4 <- (fig4_ls | plot_spacer() | patho_ma_ord) +
  plot_layout(widths = c(0.35, 0.01, 0.64)) +
  plot_annotation(tag_levels = 'A') 
#+ fig4,warning=FALSE,fig.height=4,fig.width=6.5
fig4
#' **Fig 4.** Putative plant pathogen communities in **corn**, **restored**, and **remnant** prairie fields.
#' **a** OTU richness and **b** sequence abundance are shown as columns with 95 % CIs.
#' **c** Principal-coordinate (PCoA) ordination of ITS-based (97 % OTU) community
#' distances: small points = sites, large circles = field-type centroids (error bars =
#' 95 % CI). Cornfields cluster apart from restored or remnant prairies (P < 0.01).
#' Numbers in black circles give years since restoration. Axis labels show the
#' percent variation explained. 
#' 
#+ fig4_save,warning=FALSE,fig.height=5,fig.width=7,echo=FALSE
ggsave(root_path("figs", "fig4.png"),
       plot = fig4,
       width = 6.5,
       height = 4,
       units = "in",
       dpi = 600)
#' 
#' ### Sequence-based relative abundance
d_patho <- patho %>% 
  data.frame(row.names = 1) %>% 
  decostand("total") %>%
  vegdist("bray")
mva_patho <- mva(d = d_patho, env = sites, corr = "lingoes")
#' Diagnostics/results
mva_patho$dispersion_test
mva_patho$permanova
mva_patho$pairwise_contrasts[c(1,3,2), c(1,2,4,3,8)] %>% 
  kable(format = "pandoc", caption = "Pairwise permanova contrasts")
#' 
#' Lingoes correction was needed. Three axes were significant based on a broken stick test. 
#' Based on the homogeneity of variance test, the null hypothesis 
#' of equal variance among groups is accepted across all clusters and in pairwise comparison of 
#' clusters (both p>0.05), supporting the application of a PERMANOVA test.
#' An effect of geographic distance (covariate) on pathogen communities was not supported. 
#' With geographic distance accounted for, the test variable ‘field type’ significantly explained 
#' variation in fungal communities, with a post-hoc test revealing that communities in corn 
#' fields differed from communities in restored and remnant fields.
#' 
#' Plot results
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
    x = paste0("Axis 1 (", mva_patho$axis_pct[1], "%)"),
    y = paste0("Axis 2 (", mva_patho$axis_pct[2], "%)")) +
  theme_ord +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1))
#' 
#' #### Supplemental figure
#+ patho_shan_ord_sup_patchwork,warning=FALSE
patho_shan_ord_sup <- (patho_shan_fig | plot_spacer() | patho_ord) +
  plot_layout(widths = c(0.45, 0.01, 0.55)) +
  plot_annotation(tag_levels = 'A') 
#+ patho_shan_ord_sup,warning=FALSE,fig.height=4,fig.width=6.5
patho_shan_ord_sup
#+ patho_shan_ord_sup_save,warning=FALSE,echo=FALSE
ggsave(root_path("figs", "figS6.png"),
       plot = patho_shan_ord_sup,
       width = 7.5,
       height = 4,
       units = "in",
       dpi = 600)
#' 
#' #### Contrast community metrics
#' Procrustes test on PCoA values using axes with eigenvalues exceeding a broken stick model
#+ patho_protest
set.seed(20251112)
patho_protest <- protest(
  pcoa(d_patho)$vectors[, 1:3],
  pcoa(d_patho_ma)$vectors[, 1:3],
  permutations = 1999
)
patho_protest
#' Including biomass changes little. The spatial configurations of both ordinations are highly correlated.
#' $R^{2}=$ `r round(patho_protest$scale^2, 2)`, p<0.001. 
#' 
#' ## Pathogen Indicator Species
#' Use as a tool to find species for discussion. Unbalanced design and bias to agricultural soil
#' research may make the indicator stats less appropriate for other use.
patho_ind <- inspan(its_avg, its_meta, "plant_pathogen", sites)
patho_ind %>% 
  select(A, B, stat, p_val_adj, field_type, species, starts_with("corn"), starts_with("restor"), starts_with("rem")) %>% 
  filter(species != "unidentified") %>% 
  arrange(field_type, p_val_adj) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>% 
  kable(format = "pandoc", caption = "Indicator species analysis results with abundances")
#' 
#' ## Pathogen—pfg correlations
#' Pathogen abundance correlated with grass and forbs. 
patho_resto <- its_guild %>% 
  filter(field_type == "restored", region != "FL") %>% 
  left_join(its_guild_ma %>% select(field_name, patho_mass), by = join_by(field_name)) %>% 
  mutate(
    patho_prop = patho_abund / fungi_abund, # no zeroes present...
    patho_logit = qlogis(patho_prop)
  ) %>% 
  select(-sapro_abund)
#' Inspect simple linear relationship.
pama_rest_m <- lm(patho_mass ~ gf_index, data = patho_resto)
#+ cm9,warning=FALSE,fig.width=7,fig.height=9
check_model(pama_rest_m)
summary(pama_rest_m)
#' GF index and biomass weighted relative abundance of pathogens aren't strongly 
#' related. We need to find out how total fungal biomass relates to gf_index in defining the 
#' relationship with pathogen absolute sequence abundance. Essentially, using total biomass
#' as a covariate removes the compositional element of total sequence abundance. 
#' Since sequence abundance * biomass, the composite variable, is a product, the diagnostic
#' model should probably be log-log. Check to make sure:
parest_m_raw  <- lm(patho_abund ~ fungi_mass + gf_index, data = patho_resto)
parest_m_logy <- lm(log(patho_abund) ~ fungi_mass + gf_index, data = patho_resto)
parest_m_logx <- lm(patho_abund ~ log(fungi_mass) + gf_index, data = patho_resto)
parest_m_both <- lm(log(patho_abund) ~ log(fungi_mass) + gf_index, data = patho_resto)
compare_performance(parest_m_raw, parest_m_logy, parest_m_logx, parest_m_both, 
                    metrics = c("AIC", "RMSE","R2"), rank = TRUE)
#' Log-log is best. For further investigation, Compute a pathogen proportion (row proportions 
#' of guilds).
#' 
#' ### Does plant composition shift the relative pathogen proportion?
#+ patho_logit_gfi_fig,warning=FALSE,fig.width=5,fig.height=5
ggplot(patho_resto, aes(x = gf_index, y = patho_logit)) +
  geom_text(label = rownames(patho_resto))

parest_m_rel <- lm(patho_logit ~ gf_index, data = patho_resto)
distribution_prob(parest_m_rel)
shapiro.test(parest_m_rel$residuals)
#+ cm10,warning=FALSE,fig.width=7,fig.height=9
check_model(parest_m_rel)
summary(parest_m_rel)

slopes <- map_dbl(seq_len(nrow(patho_resto)), function(i){
  coef(lm(patho_logit ~ gf_index, data = patho_resto[-i, ]))["gf_index"]
})
summary(slopes); slopes[which.min(slopes)]; slopes[which.max(slopes)]

coeftest(parest_m_rel, vcov. = vcovHC(parest_m_rel, type = "HC3"))
coefci(parest_m_rel, vcov. = vcovHC(parest_m_rel, type = "HC3"))
par(mfrow = c(1,1))
crPlots(parest_m_rel, terms = ~ gf_index)
ncvTest(parest_m_rel)
#' A linear model indicated a positive association between grass–forb index 
#' and the logit share of pathogen sequences (β = 0.67, 95% CI = 0.26-1.08, R2adj = 0.63, n = 10).
#' Results were robust to HC3 standard errors and to robust regression (bisquare), 
#' with similar slope estimates in a leave-one-out test. Model residuals did not deviate from normal 
#' based on a shapiro test. 
#' 
#' ### Does plant composition affect total fungal biomass?
#+ patho_resto_gfi_fig,warning=FALSE,fig.width=5,fig.height=5
ggplot(patho_resto, aes(x = gf_index, y = log(fungi_mass))) +
  geom_text(label = rownames(patho_resto))
parest_m_biom <- lm(log(fungi_mass) ~ gf_index, data = patho_resto)
distribution_prob(parest_m_biom)
shapiro.test(parest_m_biom$residuals) # residuals distribution non-normal
par(mfrow = c(2,2))
plot(parest_m_biom)
summary(parest_m_biom)

coeftest(parest_m_biom, vcov = vcovHC(parest_m_biom, type = "HC3")) # robust slope and SEs also NS
coefci(parest_m_biom, vcov. = vcovHC(parest_m_biom, type = "HC3"))
#' Non-normal residuals distribution suggests the need for corroboration. Try a nonparametric test.
with(patho_resto, cor.test(gf_index, log(fungi_mass), method = "spearman", exact = FALSE)) # nonparametric correlation NS
#' OLS and robust estimates agree on direction/magnitude; effect not significant
par(mfrow = c(1,1))
crPlots(parest_m_biom, terms = ~ gf_index)
ncvTest(parest_m_biom)
#' The model likely works well enough to show that there's not a strong relationship between GF index
#' and total fungal biomass, which agrees with what we found before with biomass as an untransformed
#' response variable. 
#' 
#' ### Does gf_index still matter for absolute pathogens once biomass is in the model?
#+ patho_resto_massmass_fig,warning=FALSE,fig.width=5,fig.height=5
ggplot(patho_resto, aes(x = log(fungi_mass), y = log(patho_mass))) +
  geom_text(label = rownames(patho_resto))
parest_m_abs <- lm(log(patho_mass) ~ log(fungi_mass) + gf_index, data = patho_resto)

augment(parest_m_abs)
distribution_prob(parest_m_abs)
shapiro.test(parest_m_abs$residuals)
#' Residuals distribution difference from normal rejected by shapiro test and machine learning approach
par(mfrow = c(2,2))
plot(parest_m_abs)
#' Sturcture, leverage point detected
crPlots(parest_m_abs)
avPlots(parest_m_abs)
#' Relationships appear monotonic (in log-log space). 
coeftest(parest_m_abs, vcov. = vcovHC(parest_m_abs, type = "HC3")) # Robust Wald t test
coefci(parest_m_abs, vcov. = vcovHC(parest_m_abs, type = "HC3"))
par(mfrow = c(1,1))
crPlots(parest_m_abs, terms = ~ gf_index)
ncvTest(parest_m_abs)
#' Diagnostics (Shapiro p=0.64; NCV p=0.54; HC3 and LOOCV stable) support the additive log–log model.
#' 
#' Results:
Anova(parest_m_abs, type = 2)
#+ parest_m_abs_rsq,warning=FALSE,message=FALSE
rsq.partial(parest_m_abs, adj=TRUE)$partial.rsq
#' Elasticity of pathogen mass to fungal biomass: a 1% increase in biomass ≈ 1.54% increase in pathogen mass.
#' Gf_index is multiplicative on the original scale. One-unit increase in gf_index ≈ exp(0.698) ≈ 2.01 in pathogen mass.
#' Including HC3 SE: 0.698 ± 1.96·0.2106 = exp(0.284, 1.112) = (1.33, 3.04).
#' 
#' View results:
#' Median fungal biomass on the original scale and model prediction for figure.
med_fungi <- median(patho_resto$fungi_mass, na.rm = TRUE)
newdat <- tibble(
  gf_index   = seq(min(patho_resto$gf_index, na.rm = TRUE),
                   max(patho_resto$gf_index, na.rm = TRUE),
                   length.out = 200),
  fungi_mass = med_fungi
)
# Predict on the log scale, then back-transform to the original scale
pred <- augment(parest_m_abs, newdata = newdat, se_fit = TRUE) %>%
  mutate(
    fit_med = exp(.fitted),                       # median on raw scale
    lwr_med = exp(.fitted - 1.96 * .se.fit),
    upr_med = exp(.fitted + 1.96 * .se.fit)
  )
#+ fig7_data,warning=FALSE
fig7 <- 
  ggplot(pred, aes(x = gf_index, y = fit_med)) +
  geom_ribbon(aes(ymin = lwr_med, ymax = upr_med), fill = "gray90") +
  geom_line(color = "black", linewidth = lw) +
  geom_point(data = patho_resto, aes(x = gf_index, y = patho_mass),
             fill = ft_pal[2], size = sm_size, stroke = lw, shape = 21) +
  geom_text(data = patho_resto, aes(x = gf_index, y = patho_mass, label = yr_since), 
            size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  labs(
    x = "Grass–forb index",
    y = "Pathogen mass",
    caption = paste0("Curve shown at median fungal biomass")
  ) +
  theme_cor
#+ fig7,warning=FALSE,fig.height=4,fig.width=4
fig7
#+ fig7_save,warning=FALSE,echo=FALSE
ggsave(
  root_path("figs", "fig7.png"),
  plot = fig7,
  width = 4,
  height = 4,
  units = "in",
  dpi = 600
)
#' ### Alternative explanations
#' Is plant richness related to pathogens?
#+ patho_resto_prich_fig,warning=FALSE,fig.width=5,fig.height=5
patho_resto %>% left_join(prich, by = join_by(field_name)) %>% 
  ggplot(aes(x = pl_rich, y = patho_mass)) + 
  geom_point()
#' not at all...and it wasn't selected in whole soil fungi either
with(gf_index %>% left_join(prich, by = join_by(field_name)), cor.test(gf_index, pl_rich))
#+ patho_resto_gfi_prich_fig,warning=FALSE,fig.width=5,fig.height=5
gf_index %>% left_join(prich, by = join_by(field_name)) %>% 
  ggplot(aes(x = gf_index, y = pl_rich)) +
  geom_point()
#' Plant richness relationship with GF index is weak and driven by a single outlier site. 
#' The relative proportion of grasses and forbs is related to pathogen abundance 
#' irrespective of the number of plant species. 

#' 
#' # Putative saprotrophs
# Putative saprotrophs ———————— ####
#' 
#' Retrieve pathogen sequence abundance
sapro <- guildseq(its_avg, its_meta, "saprotroph")
#' 
#' ## Diversity Indices
#+ sapro_diversity
sapro_div <- calc_div(sapro, sites)
#' 
#' ### Richness
#' Account for sequencing depth as a covariate
sapro_rich_covar <- covar_shape_test(
  data  = sapro_div,
  y     = "richness",       
  covar = "depth",   
  group = "field_type"           
)
sapro_rich_covar$compare
#' Log transform selected based on model performance. Model with centered, log depth covar
sapro_rich_lm <- lm(richness ~ depth_clg + field_type,
                    data = sapro_div %>% mutate(depth_clg = log(depth) - mean(log(depth))))
#' Diagnostics
#+ sapro_rich_covar_diagnostics,warning=FALSE,fig.width=7,fig.height=9
sapro_rich_covar$diagnostics # heavy residual structure, poor qq alignment
distribution_prob(sapro_rich_lm)
#' residuals distribution normal or close, response showing group divisions and count 
#' overdispersion (binomial family)
leveneTest(richness ~ field_type, data = sapro_div) %>% as.data.frame() %>% 
  kable(format = "pandoc", caption = "Response var in groups")
leveneTest(residuals(sapro_rich_lm) ~ sapro_div$field_type) %>% as.data.frame() %>% 
  kable(format = "pandoc", caption = "Residuals var in groups")
#' Levene's p > 0.05 → fail to reject = variances can be considered equal.
#' Despite lack of unequal variance, does gamma glm improve other diagnostics?
sapro_rich_glm  <- glm(richness ~ depth_clg + field_type, family = Gamma(link = "log"), 
                       data = sapro_div %>% mutate(depth_clg = log(depth) - mean(log(depth))))
sapro_glm_diag <- glm.diag(sapro_rich_glm)
glm.diag.plots(sapro_rich_glm, sapro_glm_diag) 
#' Slight improvement to qq plot shape, slightly reduced hatvalue (not shown)
performance::check_overdispersion(sapro_rich_glm) # not detected
sapro_glm_sim <- simulateResiduals(sapro_rich_glm)
plot(sapro_glm_sim) # DHARMa passes all tests
#' Gamma glm is the best choice; no high-leverage point
#' 
#' Model results, group means, and post-hoc
summary(sapro_rich_glm)
Anova(sapro_rich_glm, type = 2)
#' Sequence depth is significant; richness doesn't vary in groups. View trend:
sapro_div %>% 
  group_by(field_type) %>% 
  summarize(avg_richness = mean(richness),
            avg_depth = mean(depth), .groups = "drop") %>% 
  mutate(across(starts_with("avg"), ~ round(.x, 1))) %>% 
  kable(format = "pandoc")
#' Inverse relationship between depth and field_type. There was a mildly significant
#' interaction, but how to work out what that would even mean?
#' 
#' Calculate confidence intervals for figure.
#' Arithmetic means calculated in this case, back-transformed.
sapro_rich_em <- emmeans(sapro_rich_glm, ~ field_type, type = "response")
#+ sapro_rich_em_summary,echo=FALSE
kable(summary(sapro_rich_em), 
      format = "pandoc", 
      caption = "Confidence level used: 0.95")
#' Model NS; no post hoc comparison...
#+ sapro_richness_fig,fig.width=4,fig.height=4
sapro_rich_fig <- 
  ggplot(summary(sapro_rich_em), aes(x = field_type, y = response)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = response, ymax = upper.CL), width = 0, linewidth = lw) +
  labs(x = NULL, y = "Richness") +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1))
#' 
#' ### Shannon's diversity
#' Account for sequencing depth as a covariate
sapro_shan_covar <- covar_shape_test(
  data = sapro_div,
  y = "shannon",
  covar = "depth",
  group = "field_type"
)
sapro_shan_covar$compare
#' Log transform selected based on model performance. Model with centered, log depth covar
sapro_shan_lm <- lm(shannon ~ depth_clg + field_type,
                    data = sapro_div %>% mutate(depth_clg = log(depth) - mean(log(depth))))
#' Diagnostics
#+ sapro_shan_covar_diagnostics,warning=FALSE,fig.width=7,fig.height=9
sapro_shan_covar$diagnostics # variance similar in groups 
distribution_prob(sapro_shan_lm)
#' residuals distribution most likely normal, qq fit good, no evidence of mean/variance increase
#' response non-normal, check variance in groups though
leveneTest(shannon ~ field_type, data = sapro_div) %>% as.data.frame() %>% kable(format = "pandoc")
leveneTest(residuals(sapro_shan_lm) ~ sapro_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc")
#' Residuals distribution does not suggest the need for transformation.
#' Levene's p > 0.05 → fail to reject = variances can be considered equal.

#' Model results, group means, and post-hoc
sapro_shan_covar$anova_t2
#' Sequence depth is not a significant predictor of Shannon diversity, nor field type
sapro_shan_em <- emmeans(sapro_shan_lm, ~ field_type, type = "response")
#' Results tables below show the emmeans summary of group means and confidence intervals,
#' with sequencing depth as a covariate, and the post hoc contrast of richness among field types.
#+ sapro_shan_fig,fig.width=4,fig.height=4
sapro_shan_fig <- 
  ggplot(summary(sapro_shan_em), aes(x = field_type, y = emmean)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = emmean, ymax = upper.CL), width = 0, linewidth = lw) +
  labs(x = NULL, y = "Shannon diversity") +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
#' 
#' ## Abundance
#' Proportional biomass
sapro_ma_lm <- lm(sapro_mass ~ field_type, data = its_guild_ma)
par(mfrow = c(2,2))
plot(sapro_ma_lm) 
#' Variance looks consistent, no leverage points, poor qq fit
distribution_prob(sapro_ma_lm)
#' Residuals distribution fits normal, so do residuals
leveneTest(residuals(sapro_ma_lm) ~ sapro_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc") 
#' No covariate; response and residuals tests equivalent
#' Residuals distribution does not suggest the need for transformation.
#' Levene's p > 0.05 → fail to reject = variances can be considered equal (aka homoscedastic 
#' by group).
#' 
#' Produce model results, group means, and post-hoc, with arithmetic means from emmeans
Anova(sapro_ma_lm)
sapro_ma_em <- emmeans(sapro_ma_lm, ~ field_type, type = "response")
#+ sapro_ab_em_summary,echo=FALSE
kable(summary(sapro_ma_em),
      format = "pandoc",
      caption = "Confidence level used: 0.95")
#+ sapro_ab_fig,fig.width=4,fig.height=4
sapro_ma_fig <- 
  ggplot(summary(sapro_ma_em), aes(x = field_type, y = emmean)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = emmean, ymax = upper.CL), width = 0, linewidth = lw) +
  labs(x = NULL, y = expression(Biomass~(nmol[PLFA]%*%g[soil]^-1))) +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1.1))
#' 
#' ## Beta Diversity
#' Community distance handled similarly to previous
#' ### Biomass-weighted relative abundance
sapro_ma <- guildseq(its_avg_ma, its_meta, "saprotroph")
#+ sapro_ma_ord
d_sapro_ma <- sapro_ma %>% 
  data.frame(row.names = 1) %>% 
  vegdist("bray")
mva_sapro_ma <- mva(d = d_sapro_ma, env = sites)
#+ sapro_ma_ord_results
mva_sapro_ma$dispersion_test
mva_sapro_ma$permanova
mva_sapro_ma$pairwise_contrasts[c(1,3,2), c(1,2,4,3,8)] %>% 
  kable(format = "pandoc", caption = "Pairwise permanova contrasts")
#' Lingoes correction was not necessary. Based on the homogeneity of variance test, the null hypothesis of equal variance among groups is 
#' accepted across all clusters and in pairwise comparison of clusters (both p>0.05), supporting the application of 
#' a PERMANOVA test. 
#' 
#' An effect of geographic distance (covariate) on pathogen communities was detected. 
#' With geographic distance accounted for, the test variable 'field type' significantly explained 
#' variation in fungal communities, with a post-hoc test revealing that communities in corn fields differed from
#' communities in restored and remnant fields. 
#' 
#' Plotting results: 
sapro_ma_ord_data <- mva_sapro_ma$ordination_scores %>% mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant")))
p_sapro_ma_centers <- sapro_ma_ord_data %>% 
  group_by(field_type) %>% 
  summarize(across(starts_with("Axis"), list(mean = mean, ci_l = ci_l, ci_u = ci_u), .names = "{.fn}_{.col}"), .groups = "drop") %>% 
  mutate(across(c(ci_l_Axis.1, ci_u_Axis.1), ~ mean_Axis.1 + .x),
         across(c(ci_l_Axis.2, ci_u_Axis.2), ~ mean_Axis.2 + .x))
sapro_ma_ord <- 
  ggplot(sapro_ma_ord_data, aes(x = Axis.1, y = Axis.2)) +
  geom_linerange(data = p_sapro_ma_centers, aes(x = mean_Axis.1, y = mean_Axis.2, xmin = ci_l_Axis.1, xmax = ci_u_Axis.1), linewidth = lw) +
  geom_linerange(data = p_sapro_ma_centers, aes(x = mean_Axis.1, y = mean_Axis.2, ymin = ci_l_Axis.2, ymax = ci_u_Axis.2), linewidth = lw) +
  geom_point(data = p_sapro_ma_centers, aes(x = mean_Axis.1, y = mean_Axis.2, fill = field_type), size = lg_size, stroke = lw, shape = 21) +
  geom_point(aes(fill = field_type), size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  scale_fill_manual(values = ft_pal) +
  labs(
    x = paste0("Axis 1 (", mva_sapro_ma$axis_pct[1], "%)"),
    y = paste0("Axis 2 (", mva_sapro_ma$axis_pct[2], "%)")) +
  theme_ord +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1))
#' 
#' ## Unified figure
#+ fig5_patchwork,warning=FALSE
fig5_ls <- (sapro_rich_fig / plot_spacer() / sapro_ma_fig) +
  plot_layout(heights = c(1,0.01,1)) 
fig5 <- (fig5_ls | plot_spacer() | sapro_ma_ord) +
  plot_layout(widths = c(0.35, 0.01, 0.64)) +
  plot_annotation(tag_levels = 'A') 
#+ fig5,warning=FALSE,fig.height=4,fig.width=6.5
fig5
#' **Fig 5.** Putative plant pathogen communities in **corn**, **restored**, and **remnant** prairie fields.
#' **a** OTU richness and **b** sequence abundance are shown as columns with 95 % CIs.
#' **c** Principal-coordinate (PCoA) ordination of ITS-based (97 % OTU) community
#' distances: small points = sites, large circles = field-type centroids (error bars =
#' 95 % CI). Cornfields cluster apart from restored or remnant prairies (P < 0.01).
#' Numbers in black circles give years since restoration. Axis labels show the
#' percent variation explained. Colours/shading: corn = grey, restored = black,
#' remnant = white.
#+ fig5_save,warning=FALSE,fig.height=5,fig.width=7,echo=FALSE
ggsave(root_path("figs", "fig5.png"),
       plot = fig5,
       width = 6.5,
       height = 4,
       units = "in",
       dpi = 600)
#' 
#' ### Sequence-based relative abundance
#+ sapro_ord
d_sapro <- sapro %>%
  data.frame(row.names = 1) %>%
  decostand("total") %>%
  vegdist("bray")
mva_sapro <- mva(d = d_sapro, env = sites)
#+ sapro_ord_results
mva_sapro$dispersion_test
mva_sapro$permanova
mva_sapro$pairwise_contrasts[c(1,3,2), c(1,2,4,3,8)] %>%
  kable(format = "pandoc", caption = "Pairwise permanova contrasts")
#' Lingoes correction was not necessary. Based on the homogeneity of variance test, the null hypothesis of equal variance among groups is
#' accepted across all clusters and in pairwise comparison of clusters (both p>0.05), supporting the application of
#' a PERMANOVA test.
#'
#' An effect of geographic distance (covariate) on pathogen communities was detected
#' With geographic distance accounted for, the test variable 'field type' significantly explained
#' variation in fungal communities, with a post-hoc test revealing that communities in corn fields differed from
#' communities in restored and remnant fields.
#'
#' Plotting results:
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
    x = paste0("Axis 1 (", mva_sapro$axis_pct[1], "%)"),
    y = paste0("Axis 2 (", mva_sapro$axis_pct[2], "%)")) +
  theme_ord +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1))
#' 
#' #### Supplemental figure
#+ sapro_shan_ord_sup_patchwork,warning=FALSE
sapro_shan_ord_sup <- (sapro_shan_fig | plot_spacer() | sapro_ord) +
  plot_layout(widths = c(0.45, 0.01, 0.55)) +
  plot_annotation(tag_levels = 'A') 
#+ sapro_shan_ord_sup,warning=FALSE,fig.height=4,fig.width=6.5
sapro_shan_ord_sup
#+ sapro_shan_ord_sup_save,warning=FALSE,echo=FALSE
ggsave(root_path("figs", "figS7.png"),
       plot = sapro_shan_ord_sup,
       width = 7.5,
       height = 4,
       units = "in",
       dpi = 600)
#' 
#' #### Contrast community metrics
#' Procrustes test on PCoA values using axes with eigenvalues exceeding a broken stick model
#+ sapro_protest
set.seed(20251119)
sapro_protest <- protest(
  pcoa(d_sapro)$vectors[, 2],
  pcoa(d_sapro_ma)$vectors[, 2],
  permutations = 1999
)
sapro_protest
#' Including biomass changes little. The spatial configurations of both ordinations are highly correlated.
#' $R^{2}=$ `r round(sapro_protest$scale^2, 2)`, p<0.001. 
#' 
#' ## Saprotroph Indicator Species
sapro_ind <- inspan(its_avg, its_meta, "saprotroph", sites)
sapro_ind %>% 
  select(A, B, stat, p_val_adj, field_type, species, starts_with("corn"), starts_with("restor"), starts_with("rem")) %>% 
  filter(species != "unidentified") %>% 
  arrange(field_type, p_val_adj) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>% 
  kable(format = "pandoc", caption = "Indicator species analysis results with abundances")
#' 
#' ## Saprotroph correlations with pfg
sapro_resto <- its_guild %>% 
  filter(field_type == "restored", region != "FL") %>% 
  mutate(
    sapro_prop = sapro_abund / fungi_abund, # no zeroes present...
    sapro_logit = qlogis(sapro_prop)
  ) %>% 
  left_join(its_guild_ma %>% select(field_name, sapro_mass), by = join_by(field_name)) %>% 
  select(-patho_abund)

sarest_m_raw  <- lm(sapro_abund ~ fungi_mass + gf_index, data = sapro_resto)
sarest_m_logy <- lm(log(sapro_abund) ~ fungi_mass + gf_index, data = sapro_resto)
sarest_m_logx <- lm(sapro_abund ~ log(fungi_mass) + gf_index, data = sapro_resto)
sarest_m_both <- lm(log(sapro_abund) ~ log(fungi_mass) + gf_index, data = sapro_resto)

compare_performance(sarest_m_raw, sarest_m_logy, sarest_m_logx, sarest_m_both, 
                    metrics = c("AIC", "RMSE","R2"), rank = TRUE)
par(mfrow = c(2,2))
plot(sarest_m_raw)
summary(sarest_m_raw)
#+ sapro_resto_gfi_fig,warning=FALSE,fig.width=5,fig.height=5
ggplot(sapro_resto, aes(x = gf_index, sapro_prop)) +
  geom_text(label = rownames(sapro_resto))
#' No relationshp detected. High leverage point is KORP1 but it's difficult to see a 
#' significant inference based on changing it's position.