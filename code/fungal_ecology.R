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
  mutate(pl_rich = sum(c_across(where(is.numeric)) > 0),
         pl_shan = exp(diversity(c_across(where(is.numeric))))
         ) %>% 
  select(field_name = SITE, pl_rich, pl_shan) %>% 
  left_join(sites, by = join_by(field_name)) %>% 
  filter(field_type != "corn") %>% 
  ungroup()
with(prich[-c(3,6,9), ], cor.test(yr_since, pl_rich)) # Remnant fields don't have an age
#' Years since restoration isn't obviously related to plant species richness.
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
#' Define the grass_forb index
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
                                  labels = c(expression("Richness"), expression(paste("Shannon (", italic(e)^italic(H), ")")))) +
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
  geom_step(aes(color = pfg)) +
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
#' 
#' #### Unified figure
#+ pfg_fig_patchwork,warning=FALSE
pfg_pct_fig <- (plt_div / plot_spacer() / pfg_comp_fig / plot_spacer() / gf_pct_fig / plot_spacer() / gfi_loc_fig) +
  plot_layout(heights = c(1,0.01,1,0.01,0.7,0.01,0.2)) +
  plot_annotation(tag_levels = 'A') 
#+ pfg_fig,warning=FALSE,fig.height=7,fig.width=7
pfg_pct_fig
#+ pfg_fig_save_svg,warning=FALSE,echo=FALSE
ggsave(root_path("figs", "figS8.svg"), plot = pfg_pct_fig, device = "svg",
       width = 7.5, height = 7, units = "in")
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
its_div <- calc_div(its_avg, sites) %>% 
  mutate(depth_csq = sqrt(depth) - mean(sqrt(depth)))
#' 
#' ### Richness
#' Sequence depth square root transformed and centered. Negative binomial model used to handle count 
#' data. Poisson model was overdispersed (not shown). 
its_rich_glm_i <- glm.nb(richness ~ depth_csq * field_type, data = its_div)
Anova(its_rich_glm_i, type = 3, test.statistic = "LR") # no interaction detected
#' Fit additive model
its_rich_glm <- glm.nb(richness ~ depth_csq + field_type, data = its_div)
#' Diagnostics
#+ its_rich_covar_diagnostics,warning=FALSE,fig.width=7,fig.height=9
check_model(its_rich_glm)
check_overdispersion(its_rich_glm)
check_collinearity(its_rich_glm)
#' Long tails, some midrange structure, no leverage points
distribution_prob(its_rich_glm)
#' residuals distribution normal or long-tailed, response log
leveneTest(richness ~ field_type, data = its_div) %>% as.data.frame() %>% kable(format = "pandoc")
leveneTest(residuals(its_rich_glm) ~ its_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc")
#' Residuals/response distributions do not suggest the need for transformation.
#' Levene's p > 0.05 → fail to reject = variances can be considered equal across groups.
#' 
#' Model results, group means, and post-hoc. Use Type II LR test of variables due to unbalanced design.
Anova(its_rich_glm, type = 2, test.statistic = "LR")
#' Sequence depth is significant, less so than field type. 
#' Proceed with means separation by obtaining estimated marginal means for field type.
its_rich_em <- emmeans(its_rich_glm, ~ field_type, type = "response")
#' Results tables below show the emmeans summary of group means and confidence intervals,
#' with sequencing depth as a covariate, and the post hoc contrast of richness among field types. 
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
#' 
#' ### Shannon's diversity
#' Sequence depth square root transformed and centered  
its_shan_lm <- lm(shannon ~ depth_csq + field_type, data = its_div)
#' Diagnostics
#+ its_shan_covar_diagnostics,warning=FALSE,fig.width=7,fig.height=9
check_model(its_shan_lm)
#' Some residual structure, no leverage points, no evidence for increasing mean/var relationship.
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
#' Residuals' CV constant to declining.
#' Relatively low Levene's p value likely due to unequal variance in restored and remnant despite similar means.
#' Unbalanced data and possible biological reality likely causing this. No need for further transformation.
#' 
#' Model results, group means, and post-hoc. Type II SS used due to unbalanced design.
Anova(its_shan_lm, type = 2)
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
  geom_text(aes(y = upper.CL, label = c("a", "b", "b")), vjust = -1, family = "sans", size = 3.5) +
  labs(x = "Field type", y = expression(Shannon~diversity~paste("(", italic(e)^italic(H), ")"))) +
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
  labs(x = "Field type", y = expression(atop("Biomass", paste("(", nmol[PLFA], " × ", g[soil]^{-1}, ")")))) +
  # labs(x = "Field Type", y = "Biomass") +
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
#' ### Biomass-weighted relative abundance
#+ its_ord_ma
d_its_ma <- its_avg_ma %>% 
  data.frame(row.names = 1) %>% 
  vegdist("bray")
mva_its_ma <- mva(d = d_its_ma, env = sites)
#+ its_ord_ma_results
mva_its_ma$dispersion_test
mva_its_ma$permanova
mva_its_ma$pairwise_contrasts[c(1,3,2), c(1,2,4,3,8)] %>% 
  arrange(group1, desc(group2)) %>% 
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
    x = paste0("PCoA 1 (", mva_its_ma$axis_pct[1], "%)"),
    y = paste0("PCoA 2 (", mva_its_ma$axis_pct[2], "%)")) +
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
#' **a** OTU richness and **b** fungal biomass (nmol PLFA g soil^-1)
#' are shown as columns with 95 % CIs; lowercase
#' letters mark significant pairwise differences.
#' **c** Principal-coordinate (PCoA) ordination of ITS-based (97 % OTU) community
#' distances: small points = sites, large circles = field-type centroids (error bars =
#' 95 % CI). Cornfields cluster apart from restored or remnant prairies.
#' Numbers in circles give years since restoration. Axis labels show the
#' percent variation explained. 
#+ fig2_save,warning=FALSE,echo=FALSE
ggsave(root_path("figs", "fig2.svg"), plot = fig2, device = "svg",
       width = 18, height = 10.5, units = "cm")
#' 
#' ### Sequence-based relative abundance
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
  arrange(group1, desc(group2)) %>% 
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
    x = paste0("PCoA 1 (", mva_its$axis_pct[1], "%)"),
    y = paste0("PCoA 2 (", mva_its$axis_pct[2], "%)")) +
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
ggsave(root_path("figs", "figS3.svg"), plot = its_shan_ord_sup, device = "svg",
       width = 7.5, height = 4, units = "in")
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
  filter(field_type != "corn", region != "FL") %>% 
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
  filter(field_type != "corn", region != "FL") %>% 
  select(-c(field_key, field_type, region, SO4, Zn, Fe, Mn, Cu, Ca, Mg, Na))
#' 
#' Assemble explanatory variables and begin iterative selection process. 
#' Plant functional groups and traits not included here were eliminated in previous forward selection
#' procedures (not shown). 
#' Check the VIF for each explanatory variable to test for collinearity if model overfitting is 
#' detected. Then run forward selection in `dbrda()`. 
#' 
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
(mod_r2   <- RsquareAdj(mod_step, permutations = 1999))
(mod_glax <- anova(mod_step, permutations = 1999))
(mod_inax <- anova(mod_step, by = "axis", permutations = 1999))
(mod_axpct <- round(100 * mod_step$CCA$eig / sum(mod_step$CCA$eig), 1))
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
#' Create the figure objects. Figure 6 will be produced with panels from other groups, see code
#' in the saprotrophs section. 
#+ fig6_objects
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
#' 
#' ## Soil fungi and plant community correlations
#' Data for these tests
fungi_resto <- its_div %>% 
  left_join(fa %>% select(field_name, fungi_mass = fungi_18.2), by = join_by(field_name)) %>% 
  left_join(sites, by = join_by(field_name, field_type)) %>% 
  left_join(gf_index, by = join_by(field_name)) %>% 
  left_join(prich %>% select(field_name, pl_rich, pl_shan), by = join_by(field_name)) %>% 
  filter(field_type != "corn", region != "FL") %>% 
  select(field_name, fungi_ab = depth, fungi_mass, gf_index, pl_rich, pl_shan)
#' 
#' ### Plant alpha diversity and fungal biomass
#' Is plant richness related to pathogen mass?
fa_prich_lm <- lm(fungi_mass ~ pl_rich, data = fungi_resto)
summary(fa_prich_lm)
#' Fungal mass and plant richness are weakly correlated but driven by a high-leverage point (not shown).
#' When seq abund is a response and log(mass) included as a covariate, no relationship 
#' is detected (not shown).  
#' 
#' Is plant diversity related to fungal mass?
fa_pshan_lm <- lm(fungi_mass ~ pl_shan, data = fungi_resto)
summary(fa_pshan_lm)
#' Fungal biomass and plant diversity are negatively related but the correlation 
#' is not significant. It's driven almost entirely by KORP (not shown) and wouldn't be close to 
#' significant otherwise, no further testing warranted. 
#' 
#' ### Fungal biomass and grass/forb composition
#' Inspect simple linear relationship. Naïve model.
fuma_rest_m <- lm(fungi_mass ~ gf_index, data = fungi_resto)
#+ cm1,warning=FALSE,fig.width=7,fig.height=9
check_model(fuma_rest_m)
summary(fuma_rest_m)
#' The relationship is poor and needs no further analysis

#' 
#' # Arbuscular mycorrhizal fungi
# AMF ———————— ####
#' ## Diversity Indices
#+ amf_diversity
amf_div <- calc_div(amf_avg, sites) %>% 
  mutate(depth_csq = sqrt(depth) - mean(sqrt(depth)))
#' 
#' ### Richness
#' Sequence depth square root transformed and centered. Negative binomial model was underdispersed 
#' and failed to converge at default iterations; use poisson glm instead. 
amf_rich_glm_i <- glm(richness ~ depth_csq * field_type, data = amf_div, family = poisson(link = "log")) 
Anova(amf_rich_glm_i, type = 3, test.statistic = "LR") # interaction detected
check_overdispersion(amf_rich_glm_i) # not overdispersed
augment(amf_rich_glm_i) # corn site has cooks >0.9
check_collinearity(amf_rich_glm_i) # depth and field_type VIF > 26
#' An interaction was detected, but including it in the model leads to very poor diagnostics.
#' It's driven by one site in corn with high leverage, and it introduces high 
#' multicollinearity. Further, the outlier point would tend to lead to a Type II
#' error of inference, making it a conservative choice to stick with the additive model. 
#' 
#' Fit additive model
amf_rich_glm <- glm(richness ~ depth_csq + field_type, data = amf_div, family = poisson(link = "log")) 
#' Diagnostics
#+ amf_rich_covar_diagnostics,warning=FALSE,fig.width=7,fig.height=9
check_model(amf_rich_glm)
check_overdispersion(amf_rich_glm)
check_collinearity(amf_rich_glm)
#' Long tails, some midrange structure, no leverage points, overdispersion, or multicollinearity
distribution_prob(amf_rich_glm)
#' residuals distribution normal or long-tailed, response count-distributed
leveneTest(richness ~ field_type, data = amf_div) %>% as.data.frame() %>% kable(format = "pandoc")
leveneTest(residuals(amf_rich_glm) ~ amf_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc")
#' Residuals/response distributions do not suggest the need for transformation.
#' Levene's p > 0.05 → fail to reject = variances can be considered equal across groups.
#' 
#' Model results, group means, and post-hoc. Use Type II LR test of variables due to unbalanced design.
Anova(amf_rich_glm, type = 2, test.statistic = "LR")
#' Sequencing depth not a significant predictor of amf richness
amf_rich_em <- emmeans(amf_rich_glm, ~ field_type, type = "response")
#' Results tables below show the emmeans summary of estimated marginal means and confidence intervals, 
#' and the post hoc contrast of richness among field types. 
#' Main effect in model significant; pairwise contrast warranted.
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
#' 
#' ### Shannon diversity
#' Sequence depth square root transformed and centered 
amf_shan_lm <- lm(shannon ~ depth_csq + field_type, data = amf_div)
#' Diagnostics
#+ amf_shan_covar_diagnostics,warning=FALSE,fig.width=7,fig.height=9
check_model(amf_shan_lm)
#' Variance appears somewhat non-constant in groups, qqplot fit is poor, 
#' one leverage point (Cook's > 0.5), a cornfield with high richness. Mean
#' richness in corn fields is lowest; this outlier would make the pairwise contrast
#' less significant, possible Type II error which is more acceptable.
distribution_prob(amf_shan_lm)
#' Residuals/response distributions most likely normal. 
leveneTest(shannon ~ field_type, data = amf_div) %>% as.data.frame() %>% kable(format = "pandoc")
leveneTest(residuals(amf_shan_lm) ~ amf_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc")
#' Residuals/response distributions do not suggest the need for transformation.
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
  geom_text(aes(y = upper.CL, label = c("a", "b", "b")),  vjust = -1, family = "sans", size = 3.5) +
  labs(x = "Field type", y = expression(Shannon~diversity~paste("(", italic(e)^italic(H), ")"))) +
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
#' Log:
nlfa_lm_log <- lm(log(amf) ~ field_type, data = fa)
par(mfrow = c(2,2))
plot(nlfa_lm_log) # qqplot ok, one high leverage point in remnants
ncvTest(nlfa_lm_log) # p=0.16, null of constant variance not rejected
#' Gamma glm:
nlfa_glm  <- glm(amf ~ field_type, family = Gamma(link = "log"), data = fa)
nlfa_glm_diag <- glm.diag(nlfa_glm)
glm.diag.plots(nlfa_glm, nlfa_glm_diag) # qqplot shows strong fit; no leverage >0.5
performance::check_overdispersion(nlfa_glm) # not detected
#' Gamma glm is the best choice; no high-leverage point
#' 
#' Model results, group means, and post-hoc
Anova(nlfa_glm, test.statistic = "LR") 
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
  geom_text(aes(y = upper.CL, label = c("a", "b", "b")),  vjust = -1, family = "sans", size = 3.5) +
  labs(x = "Field type", y = expression(atop("Biomass", paste("(", nmol[NLFA], " × ", g[soil]^{-1}, ")")))) +
  # labs(x = "Field Type", y = "Biomass") +
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
#' ### Biomass-weighted relative abundance
#+ amf_ord_ma
d_amf_ma <- amf_avg_ma %>% 
  data.frame(row.names = 1) %>% 
  vegdist("bray")
mva_amf_ma <- mva(d = d_amf_ma, env = sites, corr = "lingoes")
#+ amf_ord_ma_results
mva_amf_ma$dispersion_test
mva_amf_ma$permanova
mva_amf_ma$pairwise_contrasts[c(1,3,2), c(1,2,4,3,8)] %>% 
  arrange(group1, desc(group2)) %>% 
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
#' biomass **b** (nmol NLFA g soil^-1), 95 % CI, letters = Tukey groups). 
#' PCoA of BC distances on proportion of biomass abundance data
#' (18S, 97 % OTUs) **c**: small points = sites, large rings = field‑type centroids ±95 % CI. 
#' Numbers in points give years since restoration. Axes show % variance. Corn clusters apart from both 
#' prairie types. 
#+ fig3_save,warning=FALSE,fig.height=5,fig.width=7,echo=FALSE
ggsave(root_path("figs", "fig3.svg"), plot = fig3, device = "svg",
       width = 18, height = 10.5, units = "cm")
#' 
#' ### Sequence-based relative abundance, unifrac distance
#' 
#+ amf_ord
d_amf <- UniFrac(amf_ps, weighted = TRUE, normalized = TRUE)
mva_amf <- mva(d = d_amf, env = sites, corr = "lingoes")
#+ amf_ord_results
mva_amf$dispersion_test
mva_amf$permanova
mva_amf$pairwise_contrasts[c(1,3,2), c(1,2,4,3,8)] %>% 
  arrange(group1, desc(group2)) %>% 
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
    x = paste0("PCoA 1 (", mva_amf$axis_pct[1], "%)"),
    y = paste0("PCoA 2 (", mva_amf$axis_pct[2], "%)")) +
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
ggsave(root_path("figs", "figS5.svg"), plot = amf_shan_ord_sup, device = "svg",
       width = 7.5, height = 4, units = "in")
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
#' The null that these solutions are unrelated
#' is rejected at p<0.001. However, the alignment isn't perfect. 
#' Clearly, the low biomass in cornfields is a driving difference in 
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
(amf_mod_r2   <- RsquareAdj(amf_mod_step, permutations = 1999))
(amf_mod_glax <- anova(amf_mod_step, permutations = 1999))
(amf_mod_inax <- anova(amf_mod_step, by = "axis", permutations = 1999))
(amf_mod_axpct <- round(100 * amf_mod_step$CCA$eig / sum(amf_mod_step$CCA$eig), 1))
amf_mod_step$anova %>% 
  as.data.frame() %>% 
  mutate(p.adj = p.adjust(`Pr(>F)`, "fdr")) %>% 
  kable(, format = "pandoc")
#' Based on permutation tests with n=1999 permutations, 
#' after accounting for inter-site pairwise distance as a covariate, the model shows a significant
#' correlation between the site ordination on fungal communities
#' and the selected explanatory variables (p<0.001). The first two constrained axes are
#' also significant (p<0.001, p<0.02). The selected variables explain $R^{2}_{\text{Adj}}$=`r round(amf_mod_r2$adj.r.squared, 3) * 100` of the community
#' variation. Selected explanatory variables are pH and the grass-forb index; see table for
#' individual p values and statistics.
#' 
#' #### AMF constrained figure
#' Produce figure objects. Code for multipanel fig 6 is shown in the saprotroph section.
#+ amf_fig6_objects
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
#' 
#' ## AM fungi and plant community correlations
#' Data for these tests
amf_resto <- amf_div %>% 
  left_join(fa %>% select(field_name, amf_mass = amf), by = join_by(field_name)) %>% 
  left_join(sites, by = join_by(field_name, field_type)) %>% 
  left_join(gf_index, by = join_by(field_name)) %>% 
  left_join(prich %>% select(field_name, pl_rich, pl_shan), by = join_by(field_name)) %>% 
  filter(field_type != "corn", region != "FL") %>% 
  select(field_name, amf_ab = depth, amf_mass, gf_index, pl_rich, pl_shan) 
#' 
#' ### Plant richness and fungal biomass
#' Is plant richness related to am fungal mass?
amfa_prich_lm <- lm(amf_mass ~ pl_rich, data = amf_resto)
summary(amfa_prich_lm)
#' AM fungal mass and plant richness are nearly perfectly unrelated. 
#' 
#' Is plant diversity related to am fungal mass?
amfa_pshan_lm <- lm(amf_mass ~ pl_shan, data = amf_resto)
summary(amfa_pshan_lm)
#' AM fungal biomass and plant diversity are positively related but only weakly so,
#' no further testing warranted. 
#' 
#' ### AM fungal biomass and grass/forb composition
#' Inspect simple linear relationship. Naïve model.
amma_rest_m <- lm(amf_mass ~ gf_index, data = amf_resto)
#+ cm7,warning=FALSE,fig.width=7,fig.height=9
check_model(amma_rest_m)
summary(amma_rest_m)
#' AM fungal mass increases with grass-forb index slightly with p value just 
#' above 0.95 alpha cutoff. Account for mass and gf_index as explanatory factors for 
#' sequence abundance.
summary(lm(amf_ab ~ amf_mass + gf_index, data = amf_resto))
#' AM fungal mass predicts increase in amf sequences, but after accounting for that, 
#' gf_index doesn't explain any new variation. 
#' Logistic regression unavailable as AM fungal sequence abundance aren't proportional of a total. 

#' 
#' # Putative plant pathogens
# Putative plant pathogens ———————— ####
#' 
#' Retrieve pathogen sequence abundance
patho <- guildseq(its_avg, its_meta, "plant_pathogen")
#' 
#' ## Diversity Indices
#+ patho_diversity
patho_div <- calc_div(patho, sites) %>% 
  mutate(depth_csq = sqrt(depth) - mean(sqrt(depth)))
#' 
#' ### Richness
#' Sequence depth square root transformed and centered. Negative binomial model was underdispersed 
#' and failed to converge at default iterations; use poisson glm instead. 
patho_rich_glm_i <- glm(richness ~ depth_csq * field_type, data = patho_div, family = poisson(link = "log")) 
Anova(patho_rich_glm_i, type = 3, test.statistic = "LR") # no interaction detected
#' Fit additive model
patho_rich_glm <- glm(richness ~ depth_csq + field_type, data = patho_div, family = poisson(link = "log")) 
#' Diagnostics
#+ patho_rich_covar_diagnostics,warning=FALSE,fig.width=7,fig.height=9
check_model(patho_rich_glm)
check_overdispersion(patho_rich_glm)
check_collinearity(patho_rich_glm)
#' Some midrange structure, no leverage points, overdispersion, or multicollinearity
distribution_prob(patho_rich_glm)
#' residuals distribution normal or long-tailed, response count-distributed
leveneTest(richness ~ field_type, data = patho_div) %>% as.data.frame() %>% kable(format = "pandoc")
leveneTest(residuals(patho_rich_glm) ~ patho_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc")
#' Residuals/response distributions do not suggest the need for transformation.
#' Levene's p > 0.05 → fail to reject = variances can be considered equal across groups.
#' 
#' Model results, group means, and post-hoc. Use Type II LR test of variables due to unbalanced design.
Anova(patho_rich_glm, type = 2, test.statistic = "LR")
#' Sequence depth is highly significant; richness doesn't vary in groups. 
#+ patho_depth_ft_cor
patho_div %>% 
  group_by(field_type) %>% 
  summarize(across(c(depth, richness), ~ round(mean(.x), 0))) %>% 
  kable(format = "pandoc", caption = "Average sequence depth and pathogen richness in field types")
#' Depth is correlated with richness in field types. Differences in richness are small
#' and with depth variance removed first, this explains why richness isn't significantly 
#' different. 
#' 
#' Calculate confidence intervals for figure.
#' Arithmetic means calculated in this case.
patho_rich_em <- emmeans(patho_rich_glm, ~ field_type, type = "response")
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
  ggplot(summary(patho_rich_em), aes(x = field_type, y = rate)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = rate, ymax = asymp.UCL), width = 0, linewidth = lw) +
  labs(x = NULL, y = expression(atop("Richness", paste("(", italic(n), " OTUs)")))) +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1))
#' 
#' ### Shannon's diversity
#' Sequence depth square root transformed and centered 
patho_shan_lm <- lm(shannon ~ depth_csq + field_type, data = patho_div)
#' Diagnostics
#+ patho_shan_covar_diagnostics,warning=FALSE,fig.width=7,fig.height=9
check_model(patho_shan_lm)
distribution_prob(patho_shan_lm)
#' residuals distribution most likely cauchy/normal; symmetric but long tails
#' response normal
leveneTest(shannon ~ field_type, data = patho_div) %>% as.data.frame() %>% kable(format = "pandoc")
leveneTest(residuals(patho_shan_lm) ~ patho_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc")
#' Residuals distribution does not suggest the need for further model selection.
#' Levene's p > 0.05 → fail to reject = variances can be considered equal.
#' 
#' Model results, group means, and post-hoc
Anova(patho_shan_lm, type = 2)
#' Neither predictor is significant
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
  labs(x = "Field type", y = expression(Shannon~diversity~paste("(", italic(e)^italic(H), ")"))) +
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
  labs(x = "Field type", y = expression(atop("Biomass (scaled)", paste(bold(`(`), "(", nmol[PLFA], " × ", g[soil]^{-1}, ")", " × ", paste("(rel. abund)", bold(`)`)))))) +
  # labs(x = "Field Type", y = "Biomass (scaled)") +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1.1))
#' 
#' ## Beta Diversity
#' Community distances handled similarly to previous
#' 
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
  arrange(group1, desc(group2)) %>% 
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
#' 
#' #### Unified figure
#+ fig4_patchwork,warning=FALSE
fig4_ls <- (patho_rich_fig / plot_spacer() / patho_ma_fig) +
  plot_layout(heights = c(1,0.01,1)) 
fig4 <- (fig4_ls | plot_spacer() | patho_ma_ord) +
  plot_layout(widths = c(0.35, 0.01, 0.64)) +
  plot_annotation(tag_levels = 'A') 
#+ fig4,warning=FALSE,fig.height=4,fig.width=6.5
fig4
#' **Fig 4.** Putative plant pathogen communities in **corn**, **restored**, and **remnant** prairie fields.
#' **a** OTU richness and **b** biomass (nmol PLFA g soil^-1 * (proportional sequence abundance)) 
#' are shown as columns with 95 % CIs.
#' **c** Principal-coordinate (PCoA) ordination of ITS-based (97 % OTU) community
#' distances: small points = sites, large circles = field-type centroids (error bars =
#' 95 % CI). Cornfields cluster apart from restored or remnant prairies (P < 0.01).
#' Numbers in black circles give years since restoration. Axis labels show the
#' percent variation explained. 
#' 
#+ fig4_save,warning=FALSE,fig.height=5,fig.width=7,echo=FALSE
ggsave(root_path("figs", "fig4.svg"), plot = fig4, device = "svg",
       width = 18, height = 10.5, units = "cm")
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
  arrange(group1, desc(group2)) %>% 
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
    x = paste0("PCoA 1 (", mva_patho$axis_pct[1], "%)"),
    y = paste0("PCoA 2 (", mva_patho$axis_pct[2], "%)")) +
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
ggsave(root_path("figs", "figS6.svg"), plot = patho_shan_ord_sup, device = "svg",
       width = 7.5, height = 4, units = "in")
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
#' research may make the indicator stats less appropriate for other use. Using biomass-weighted relative
#' abundance for this differential analysis. 
patho_ind <- inspan(its_avg_ma, its_meta, "plant_pathogen", sites)
patho_ind %>% 
  select(A, B, stat, p_val_adj, field_type, species, starts_with("corn"), starts_with("restor"), starts_with("rem")) %>% 
  filter(species != "unidentified") %>% 
  arrange(field_type, p_val_adj) %>% 
  mutate(across(A:p_val_adj, ~ round(.x, 3)),
         across(corn_avg:remnant_ci, ~ num(.x, notation = "sci"))) %>% 
  kable(format = "pandoc", caption = "Indicator species analysis results with biomass-aware relative abundances in field types")
#' 
#' ## Constrained Analysis
#' Env covars processed in the ITS section (see above)
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
#' 
#' ### Constrained Analysis Results
patho_mod_step
(patho_mod_r2   <- RsquareAdj(patho_mod_step, permutations = 1999))
(patho_mod_glax <- anova(patho_mod_step, permutations = 1999))
(patho_mod_inax <- anova(patho_mod_step, by = "axis", permutations = 1999))
(patho_mod_axpct <- round(100 * patho_mod_step$CCA$eig / sum(patho_mod_step$CCA$eig), 1))
#' Based on permutation tests with n=1999 permutations, 
#' after accounting for inter-site pairwise distance as a covariate, the model shows 
#' no significant correlation between pathogen community turnover and explanatory variables.
#' 
#' #### Pathogen (un)constrained figure
#' Produce figure objects. In this case the panel will show a PCoA ordination of sites based on 
#' pathogen communities. Code for multipanel fig 6 is shown in the saprotroph section.
#+ patho_fig6_objects
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
#' 
#' ## Pathogen and plant community correlations
#' Data for these tests
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
#' 
#' ### Plant richness and pathogen biomass
#' Is plant richness related to pathogen mass?
pathofa_prich_lm = lm(patho_mass ~ pl_rich, data = patho_resto)
summary(pathofa_prich_lm)
#' Pathogen mass and plant richness aren't correlated, though the direction is negative. 
#' Relationship is weak enough that no further tests are warranted.
#' 
#' Is plant diversity related to pathogen mass?
pathofa_pshan_lm = lm(patho_mass ~ pl_shan, data = patho_resto)
summary(pathofa_pshan_lm)
#' NS  
#' 
#' ### Pathogen biomass and grass/forb composition
#' Inspect simple linear relationship. Naïve model.
patho_gf_lm <- lm(patho_mass ~ gf_index, data = patho_resto)
#+ cm9,warning=FALSE,fig.width=7,fig.height=9
check_model(patho_gf_lm)
summary(patho_gf_lm)
#' The naïve model shows a positive relationship that isn't significant. 
#' Diagnostic reveals noisy fit and lots of structure. Decompose biomass and 
#' sequence abundance in a model to test changes in each given the other. 
#' 
#' Multiple, weighted logistic glm
patho_gf_glm <- glm(patho_prop ~ fungi_mass_lc + gf_index,
                    data = patho_resto, family = quasibinomial(link = "logit"),
                    weights = fungi_abund)
summary(patho_gf_glm) # dispersion parameter >117 justifies quasibinomial
#' Diagnostics
check_model(patho_gf_glm)
check_collinearity(patho_gf_glm)
augment(patho_gf_glm)
#' Long tails and low n showing structure. Moderate leverage at LPRP1: high 
#' pathogens, high gf_index but very low biomass...this is evidence of the noise that 
#' caused the naïve model to fail. 
distribution_prob(patho_gf_glm)
#' Residuals distribution normal
loocv_paglm_gfi <- map_dbl(seq_len(nrow(patho_resto)), function(i){
  exp(coef(glm(patho_prop ~ fungi_mass_lc + gf_index, 
           data = patho_resto[-i, ], 
           family = quasibinomial(link = "logit"),
           weights = fungi_abund))["gf_index"])
})
summary(loocv_paglm_gfi)
(cv_paglm <- (sd(loocv_paglm_gfi) / mean(loocv_paglm_gfi) * 100) %>% round(., 1))
#' Grass-forb index LOOCV variation of 12.6% on the back-transformed scale shows that the influential 
#' points (LPRP1) and three other potential outliers from the qq plot do not 
#' significantly affect fit. Sign and magnitude of LOO slopes wouldn't change inference. 
loocv_paglm_fma <- map_dbl(seq_len(nrow(patho_resto)), function(i){
  exp(coef(glm(patho_prop ~ fungi_mass_lc + gf_index, 
               data = patho_resto[-i, ], 
               family = quasibinomial(link = "logit"),
               weights = fungi_abund))["fungi_mass_lc"])
})
summary(loocv_paglm_fma)
(cv_paglm <- (sd(loocv_paglm_fma) / mean(loocv_paglm_fma) * 100) %>% round(., 1))
#' Similarly, fungal mass LOOCV variation of 12.9% on the back-transformed scale shows that the influential 
#' points (LPRP1) and three other potential outliers from the qq plot do not 
#' significantly affect fit. Sign and magnitude of LOO slopes wouldn't change inference. 
#' The higher variability here shows that the noise of PLFA variation is substantial.
#' View partial regression plots for consistency and save data for plotting.
paglm_crpldata <- as.data.frame(crPlots(patho_gf_glm))
#' Noise in fungal mass data is obvious here. Fit of partial gf_index is clean.
#+ parest_m_abs_rsq,warning=FALSE,message=FALSE
rsq.partial(patho_gf_glm, adj = TRUE)$partial.rsq
#' Create objects for plotting
paglm_med_fungi <- median(patho_resto$fungi_mass_lc, na.rm = TRUE)
paglm_med_abund <- median(patho_resto$fungi_abund, na.rm = TRUE) # Needed for weight context
paglm_newdat <- tibble(
  gf_index = seq(min(patho_resto$gf_index, na.rm = TRUE),
                 max(patho_resto$gf_index, na.rm = TRUE),
                 length.out = 200),
  fungi_mass_lc = paglm_med_fungi,
  fungi_abund = paglm_med_abund 
)
#' Predict on link scale, back-transform with plogis
paglm_pred <- predict(patho_gf_glm, newdata = paglm_newdat, type = "link", se.fit = TRUE) %>%
  as_tibble() %>%
  bind_cols(paglm_newdat) %>%
  mutate(
    fit_prob = plogis(fit),
    lwr_prob = plogis(fit - 1.96 * se.fit),
    upr_prob = plogis(fit + 1.96 * se.fit)
  )
#+ fig7a,warning=FALSE
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
#+ fig7b,warning=FALSE
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
#+ fig7c,warning=FALSE
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
#+ fig7_patchwork
fig7rh <- (fig7b / plot_spacer() / fig7c) +
  plot_layout(heights = c(0.5, 0.01,0.5))
fig7 <- (fig7a | plot_spacer() | (wrap_elements(full = fig7yax_grob) & theme(plot.tag = element_blank())) | fig7rh) +
  plot_layout(widths = c(0.62, 0.004, 0.02, 0.38))
#+ fig7,warning=FALSE,message=FALSE,fig.height=5,fig.width=7
fig7
#+ fig7_save,warning=FALSE,message=FALSE,echo=FALSE
ggsave(root_path("figs", "fig7.svg"), plot = fig7, device = "svg",
       width = 18, height = 11, units = "cm")

#' 
#' # Putative saprotrophs
# Putative saprotrophs ———————— ####
#' 
#' Retrieve pathogen sequence abundance
sapro <- guildseq(its_avg, its_meta, "saprotroph")
#' 
#' ## Diversity Indices
#+ sapro_diversity
sapro_div <- calc_div(sapro, sites) %>% 
  mutate(depth_csq = sqrt(depth) - mean(sqrt(depth)))
#' 
#' ### Richness
#' Sequence depth square root transformed and centered. Poisson model was overdispersed (not shown), 
#' use negative binomial instead.  
sapro_rich_glm_i <- glm.nb(richness ~ depth_csq * field_type, data = sapro_div) 
Anova(sapro_rich_glm_i, type = 3, test.statistic = "LR") # interaction detected
check_model(sapro_rich_glm_i)
check_overdispersion(sapro_rich_glm_i) # not overdispersed
augment(sapro_rich_glm_i) # corn site has cooks >0.9
check_collinearity(sapro_rich_glm_i) # depth and interaction VIF > 6
#' An interaction was detected, but including it in the model leads to very poor diagnostics.
#' It's driven by one site in corn with high leverage, and it introduces high 
#' multicollinearity.  
#' 
#' Fit additive model
sapro_rich_glm <- glm.nb(richness ~ depth_csq + field_type, data = sapro_div) 
#' Diagnostics
#+ sapro_rich_covar_diagnostics,warning=FALSE,fig.width=7,fig.height=9
check_model(sapro_rich_glm)
check_overdispersion(sapro_rich_glm)
check_collinearity(sapro_rich_glm)
#' Long tails, some structure throughout, no leverage points, overdispersion, or multicollinearity
distribution_prob(sapro_rich_glm)
#' residuals distribution normal or long-tailed, response count-distributed
leveneTest(richness ~ field_type, data = sapro_div) %>% as.data.frame() %>% kable(format = "pandoc")
leveneTest(residuals(sapro_rich_glm) ~ sapro_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc")
#' Residuals/response distributions do not suggest the need for transformation.
#' Levene's p > 0.05 → fail to reject = variances can be considered equal across groups.
#' 
#' Model results, group means, and post-hoc. Use Type II LR test of variables due to unbalanced design.
Anova(sapro_rich_glm, type = 2, test.statistic = "LR")
#' Differences in richness are very close to significance. Calculate confidence intervals for figure.
#' Estimated marginal means calculated in this case
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
  geom_errorbar(aes(ymin = response, ymax = asymp.UCL), width = 0, linewidth = lw) +
  labs(x = NULL, y = expression(atop("Richness", paste("(", italic(n), " OTUs)")))) +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1))
#' 
#' ### Shannon's diversity
#' Sequence depth square root transformed and centered 
sapro_shan_lm <- lm(shannon ~ depth_csq + field_type, data = sapro_div)
#' Diagnostics
#+ sapro_shan_covar_diagnostics,warning=FALSE,fig.width=7,fig.height=9
check_model(sapro_shan_lm)
distribution_prob(sapro_shan_lm)
#' residuals distribution most likely normal, qq fit good, no evidence of mean/variance increase
#' response non-normal, check variance in groups though
leveneTest(shannon ~ field_type, data = sapro_div) %>% as.data.frame() %>% kable(format = "pandoc")
leveneTest(residuals(sapro_shan_lm) ~ sapro_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc")
#' Residuals distribution does not suggest the need for transformation.
#' Levene's p > 0.05 → fail to reject = variances can be considered equal.
#' 
#' Model results, group means, and post-hoc
Anova(sapro_shan_lm, type = 2)
#' Sequence depth is not a significant predictor of Shannon diversity, nor field type
sapro_shan_em <- emmeans(sapro_shan_lm, ~ field_type, type = "response")
#' Results tables below show the emmeans summary of group means and confidence intervals,
#' with sequencing depth as a covariate, and the post hoc contrast of richness among field types.
#+ sapro_shan_fig,fig.width=4,fig.height=4
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
anova(sapro_ma_lm)
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
  labs(x = "Field type", y = expression(atop("Biomass (scaled)", paste(bold(`(`), "(", nmol[PLFA], " × ", g[soil]^{-1}, ")", " × ", paste("(rel. abund)", bold(`)`)))))) +
  # labs(x = "Field Type", y = "Biomass (scaled)") +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1.1))
#' 
#' ## Beta Diversity
#' Community distance handled similarly to previous
#' 
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
  arrange(group1, desc(group2)) %>% 
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
#' 
#' #### Unified figure
#+ fig5_patchwork,warning=FALSE
fig5_ls <- (sapro_rich_fig / plot_spacer() / sapro_ma_fig) +
  plot_layout(heights = c(1,0.01,1)) 
fig5 <- (fig5_ls | plot_spacer() | sapro_ma_ord) +
  plot_layout(widths = c(0.35, 0.01, 0.64)) +
  plot_annotation(tag_levels = 'A') 
#+ fig5,warning=FALSE,fig.height=4,fig.width=6.5
fig5
#' **Fig 5.** Putative plant pathogen communities in **corn**, **restored**, and **remnant** prairie fields.
#' **a** OTU richness and **b** biomass (nmol PLFA g soil^-1 * (proportional sequence abundance)) 
#' are shown as columns with 95 % CIs.
#' **c** Principal-coordinate (PCoA) ordination of ITS-based (97 % OTU) community
#' distances: small points = sites, large circles = field-type centroids (error bars =
#' 95 % CI). Cornfields cluster apart from restored or remnant prairies (P < 0.01).
#' Numbers in black circles give years since restoration. Axis labels show the
#' percent variation explained. Colours/shading: corn = grey, restored = black,
#' remnant = white.
#+ fig5_save,warning=FALSE,fig.height=5,fig.width=7,echo=FALSE
ggsave(root_path("figs", "fig5.svg"), plot = fig5, device = "svg",
       width = 18, height = 10.5, units = "cm")
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
  arrange(group1, desc(group2)) %>% 
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
    x = paste0("PCoA 1 (", mva_sapro$axis_pct[1], "%)"),
    y = paste0("PCoA 2 (", mva_sapro$axis_pct[2], "%)")) +
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
ggsave(root_path("figs", "figS7.svg"), plot = sapro_shan_ord_sup, device = "svg",
       width = 7.5, height = 4, units = "in")
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
sapro_ind <- inspan(its_avg_ma, its_meta, "saprotroph", sites)
sapro_ind %>% 
  select(A, B, stat, p_val_adj, field_type, species, starts_with("corn"), starts_with("restor"), starts_with("rem")) %>% 
  filter(species != "unidentified") %>% 
  arrange(field_type, p_val_adj) %>% 
  mutate(across(A:p_val_adj, ~ round(.x, 3)),
         across(corn_avg:remnant_ci, ~ num(.x, notation = "sci"))) %>% 
  kable(format = "pandoc", caption = "Indicator species analysis results with biomass-aware relative abundances in field types")
#' 
#' ## Constrained Analysis
#' Env covars processed in the ITS section (see above)
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
#' 
#' ### Constrained Analysis Results
sapro_mod_step
(sapro_mod_r2   <- RsquareAdj(sapro_mod_step, permutations = 1999))
(sapro_mod_glax <- anova(sapro_mod_step, permutations = 1999))
(sapro_mod_inax <- anova(sapro_mod_step, by = "axis", permutations = 1999))
(sapro_mod_axpct <- round(100 * sapro_mod_step$CCA$eig / sum(sapro_mod_step$CCA$eig), 1))
sapro_mod_step$anova %>% 
  as.data.frame() %>% 
  mutate(p.adj = p.adjust(`Pr(>F)`, "fdr")) %>% 
  kable(, format = "pandoc")
#' Based on permutation tests with n=1999 permutations, 
#' after accounting for inter-site pairwise distance as a covariate, the model shows 
#' correlations between the site ordination on saprotroph communities
#' and the selected explanatory variables (p<0.001). The first four constrained axes are
#' also significant (p<0.05). The selected variables explain $R^{2}_{\text{Adj}}$=
#' `r round(sapro_mod_r2$adj.r.squared * 100, 1)`% of the community
#' variation. Selected explanatory variables are SOM, grass-forb index, plant richness,
#' and nitrate; see table for individual p values and statistics.
#' 
#' #### Saprotroph constrained figure
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
#' 
#' ## Saprotroph and plant community correlations
#' Data for these tests
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
#' 
#' ### Plant richness and saprotroph biomass
#' Is plant richness related to saprotroph mass?
saprofa_prich_lm <- lm(sapro_mass ~ pl_rich, data = sapro_resto)
summary(saprofa_prich_lm)
#' Strong negative relationship worthy of closer examination. It isn't driven by 
#' grass-forb index (see below), so this isn't a confounding with C4 grass abundance...
#' KORP site appears to have some leverage (not shown), conduct logistic model 
#' and diagnostics as before.
#' 
#' Multiple, weighted logistic glm
sapro_prich_glm <- glm(sapro_prop ~ fungi_mass_lc + pl_rich,
                       data = sapro_resto, family = quasibinomial(link = "logit"),
                       weights = fungi_abund)
summary(sapro_prich_glm) # dispersion parameter >62 justifies quasibinomial
#' Diagnostics
check_model(sapro_prich_glm)
check_collinearity(sapro_prich_glm)
augment(sapro_prich_glm)
#' Long tails and low n showing moderate structure. 
distribution_prob(sapro_prich_glm)
#' Residuals distribution long-tailed / normal
loocv_saglm_gfi <- map_dbl(seq_len(nrow(sapro_resto)), function(i){
  exp(coef(glm(sapro_prop ~ fungi_mass_lc + pl_rich,
               data = sapro_resto[-i, ],
               family = quasibinomial(link = "logit"),
               weights = fungi_abund))["pl_rich"])
})
summary(loocv_saglm_gfi)
(cv_saglm <- (sd(loocv_saglm_gfi) / mean(loocv_saglm_gfi) * 100) %>% round(., 1))
#' Grass-forb index LOOCV variation of 0.1% on the back-transformed scale shows 
#' remarkably consistent prediction. 
loocv_saglm_fma <- map_dbl(seq_len(nrow(sapro_resto)), function(i){
  exp(coef(glm(sapro_prop ~ fungi_mass_lc + pl_rich,
               data = sapro_resto[-i, ],
               family = quasibinomial(link = "logit"),
               weights = fungi_abund))["fungi_mass_lc"])
})
summary(loocv_saglm_fma)
(cv_saglm <- (sd(loocv_saglm_fma) / mean(loocv_saglm_fma) * 100) %>% round(., 1))
#' Slope CV on fungal mass of 4.6% is low but shows greater noise with the covariate
#' than the test variable. 
saglm_crpldata <- as.data.frame(crPlots(sapro_prich_glm))
#' Noise in fungal mass data is obvious here. Fit of partial gf_index is clean.
#+ sarest_m_abs_rsq,warning=FALSE,message=FALSE
rsq.partial(sapro_prich_glm, adj = TRUE)$partial.rsq
#' Create objects for plotting
saglm_med_fungi <- median(sapro_resto$fungi_mass_lc, na.rm = TRUE)
saglm_med_abund <- median(sapro_resto$fungi_abund, na.rm = TRUE) # Needed for weight context
saglm_newdat <- tibble(
  pl_rich = seq(min(sapro_resto$pl_rich, na.rm = TRUE),
                max(sapro_resto$pl_rich, na.rm = TRUE),
                length.out = 200),
  fungi_mass_lc = saglm_med_fungi,
  fungi_abund = saglm_med_abund
)
#' Predict on link scale, back-transform with plogis
saglm_pred <- predict(sapro_prich_glm, newdata = saglm_newdat, type = "link", se.fit = TRUE) %>%
  as_tibble() %>%
  bind_cols(saglm_newdat) %>%
  mutate(
    fit_prob = plogis(fit),
    lwr_prob = plogis(fit - 1.96 * se.fit),
    upr_prob = plogis(fit + 1.96 * se.fit)
  )
#+ fig8a,warning=FALSE
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
#+ fig8b,warning=FALSE
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
#+ fig8c,warning=FALSE
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
#+ fig8_patchwork
fig8rh <- (fig8b / plot_spacer() / fig8c) +
  plot_layout(heights = c(0.5, 0.01,0.5))
fig8 <- (fig8a | plot_spacer() | (wrap_elements(full = fig8yax_grob) & theme(plot.tag = element_blank())) | fig8rh) +
  plot_layout(widths = c(0.62, 0.004, 0.02, 0.38))
#+ fig8,warning=FALSE,message=FALSE,fig.height=5,fig.width=7
fig8
#+ fig8_save,warning=FALSE,message=FALSE,echo=FALSE
ggsave(root_path("figs", "fig8.svg"), plot = fig8, device = "svg",
       width = 18, height = 11, units = "cm")
#' 
#' ### Plant diversity and saprotroph biomass
#' Is plant diversity related to saprotroph mass?
saprofa_pshan_lm <- lm(sapro_mass ~ pl_shan, data = sapro_resto)
summary(saprofa_pshan_lm)
#' Saprotroph mass and plant diversity may also deserve further examination.
#' 
sapro_pshan_glm <- glm(sapro_prop ~ fungi_mass_lc + pl_shan,
                       data = sapro_resto, family = quasibinomial(link = "logit"),
                       weights = fungi_abund)
summary(sapro_pshan_glm)
#' NS
#' 
#' ### Saprotroph biomass and grass/forb composition
sama_rest_m <- lm(sapro_mass ~ gf_index, data = sapro_resto)
summary(sama_rest_m)
#' Likely not enough of a relationship to warrant further attention
#' Examine mass+richness models, similar to pathogen models
sapro_gf_glm <- glm(sapro_prop ~ fungi_mass_lc + gf_index,
                    data = sapro_resto, family = quasibinomial(link = "logit"),
                    weights = fungi_abund)
summary(sapro_gf_glm)
#' NS

#' 
#' # Summary results
#' The following outputs display unified summary reports where needed for tables, etc.
#' 
#' ## Richness summary
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
#' 
#' ## Shannon diversity summary
#' Fungal OTU Shannon Diversity differences across field types accounting for sequencing depth.
#' Field type effects were evaluated using Type II Analysis of Deviance.
#' P-values for field type were adjusted for multiple comparisons
#' across fungal groups using the Benjamini-Hochberg procedure.
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
#' 
#' ## Biomass (abundance) summary
#' Fungal biomass differences differences across field types.
#' Field type effects were evaluated using ANOVA or Analysis of Deviance.
#' P-values for field type were adjusted for multiple comparisons
#' across fungal groups using the Benjamini-Hochberg procedure."
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
#' 
#' ## Beta diversity summary
#' ### Using abundance-scaled biomass
#' Fungal community differences differences among field types.
#' Field type effects were evaluated using Permanova.
#' P-values for field type were adjusted for multiple comparisons
#' across fungal groups using the Benjamini-Hochberg procedure."
#+ permanova_summary,warning=FALSE,message=FALSE
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
  select(guild, term, `Pseudo_F_(df)`, R2, p.value, p.adj) %>% 
  kable(format = "pandoc")
#' 
#' #' ### Using relative sequence abundance
#' Fungal community differences differences among field types.
#' Field type effects were evaluated using Permanova.
#' P-values for field type were adjusted for multiple comparisons
#' across fungal groups using the Benjamini-Hochberg procedure."
#+ permanova_summary,warning=FALSE,message=FALSE
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
  select(guild, term, `Pseudo_F_(df)`, R2, p.value, p.adj) %>% 
  kable(format = "pandoc")
#' 
#' ## Constrained analysis summary
#' Environmental drivers were identified via partial distance-based Redundancy Analysis (db-RDA) 
#' using forward selection. 
#' Geographic distance (PCoA Axis 1) was included as a conditional term to partial out spatial effects. 
#' Radj2 represents the cumulative variance explained by the final selected model. 
#' P-values are based on 1,999 permutations; Padj reflects FDR correction within the guild.
#= dbrda_r2
list(
  all_fungi = mod_r2$adj.r.squared,
  amf = amf_mod_r2$adj.r.squared,
  saprotrophs = sapro_mod_r2$adj.r.squared
) %>% map(\(val) tibble(r2 = round(val, 3))) %>% 
  bind_rows(.id = "guild") %>% 
  kable(format = "pandoc")
#+ dbrda_var_summary
list(
  all_fungi = mod_step$anova %>% 
    as.data.frame() %>% 
    mutate(p.adj = p.adjust(`Pr(>F)`, "fdr"),
           across(where(is.numeric), ~ round(.x, 3)), 
           `F_(df)` = paste0(F, " (", Df, ", ", mod_inax["Residual", "Df"], ")")),
  amf = amf_mod_step$anova %>% 
    as.data.frame() %>% 
    mutate(p.adj = p.adjust(`Pr(>F)`, "fdr"),
           across(where(is.numeric), ~ round(.x, 3)), 
           `F_(df)` = paste0(F, " (", Df, ", ", amf_mod_inax["Residual", "Df"], ")")),
  saprotrophs = sapro_mod_step$anova %>% 
    as.data.frame() %>% 
    mutate(p.adj = p.adjust(`Pr(>F)`, "fdr"),
           across(where(is.numeric), ~ round(.x, 3)), 
           `F_(df)` = paste0(F, " (", Df, ", ", sapro_mod_inax["Residual", "Df"], ")"))
) %>% 
  map(\(df) df %>% select(`F_(df)`, P=`Pr(>F)`, P_adj = p.adj)) %>% 
  kable(format = "pandoc")
  
#+ dbrda_axis_summary,message=FALSE,warning=FALSE
tidy(mod_inax)


#' All soil fungi
#+ fig6a
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
#' AMF
#+ fig6b
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
#' Pathogens, PCoA fig
#+ fig6c
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
#' Saprotrophs
#+ fig6d
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
#' 
#' #### Unified figure
#' Display results of constrained analyses
#+ fig6_patchwork,warning=FALSE
fig6up <- (fig6a | plot_spacer() | fig6b) +
  plot_layout(widths = c(0.50, 0.01, 0.50))
fig6dn <- (fig6c | plot_spacer() | fig6d) +
  plot_layout(widths = c(0.50, 0.01, 0.50))
fig6 <- (fig6up / fig6dn) +
  plot_layout(widths = c(0.50, 0.50)) +
  plot_annotation(tag_levels = 'A')
#+ fig6,warning=FALSE,fig.height=7,fig.width=7
fig6
#' Fungal community ordinations which are constrained or unconstrained by explanatory variables. 
#' Panels show results for all soil fungi **a**, amf **b**, pathogens **c**, and saprotrophs **d**.
#' Percent of constrained (db-RDA) and unconstrained (PCoA) variation explained is shown with axis labels.
#' For explanatory variables with significant community correlations, blue arrows show the grass-forb index 
#' with labels indicating the direction of relative increase in 
#' C4 grasses or forbs, respectively, along the index. The black arrows show other significant constraining
#' variables. Points show locations of restored fields (green) and remnant fields (blue) in Wisconsin. 
#+ fig6_save,warning=FALSE,echo=FALSE
ggsave(root_path("figs", "fig6.svg"), plot = fig6, device = "svg",
       width = 18, height = 18, units = "cm")


