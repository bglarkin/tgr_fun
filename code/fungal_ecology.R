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
#' based on response and residuals distributions; √‑transformation of sequencing depth used as covariate 
#' per [Bálint 2015](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.13018); 
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
  "ggpubr", "patchwork", "car", "performance", "broom", "boot", "indicspecies",
  "MASS", "DHARMa"
)

to_install <- setdiff(packages_needed, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)
invisible(lapply(packages_needed, library, character.only = TRUE))

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
#' CSV files were produced in `sequence_data.R`. Amf_avg_uni table is in species-samples format
#' to enable use of `Unifrac()` later.
spe <- list(
    its_avg     = read_csv(root_path("clean_data/spe_ITS_avg.csv"), show_col_types = FALSE),
    amf_avg     = read_csv(root_path("clean_data/spe_18S_avg.csv"), show_col_types = FALSE),
    amf_avg_uni = read_delim(root_path("otu_tables/18S/18S_avg_4unifrac.tsv"), show_col_types = FALSE)
)
#' 
#' ## Microbial species metadata
spe_meta <- list(
    its = read_csv(root_path("clean_data/spe_ITS_metadata.csv"), show_col_types = FALSE) %>% 
      mutate(primary_lifestyle = case_when(str_detect(primary_lifestyle, "_saprotroph$") ~ "saprotroph", 
                                           TRUE ~ primary_lifestyle)),
    amf = read_csv(root_path("clean_data/spe_18S_metadata.csv"), show_col_types = FALSE),
    amf_avg_uni = read_delim(root_path("otu_tables/18S/18S_avg_4unifrac.tsv"), show_col_types = FALSE)
) %>% 
  map(. %>% mutate(across(everything(), ~ replace_na(., "unidentified"))))

#' ## Plant data
#' Abundance in functional groups and by species are only available from Wisconsin sites. 
#' Only C4_grass and forbs are used. Others: C3_grass, legume, and shrubTree were found 
#' previously to have high VIF in models or were not chosen in forward selection. 
pfg <- read_csv(root_path("clean_data", "plant_traits.csv"), show_col_types = FALSE) %>% 
  select(field_name, C4_grass, forb)

#' ## Soil properties
soil <- read_csv(root_path("clean_data/soil.csv"), show_col_types = FALSE)[-c(26:27), ]

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
#' ### Wrangle species and metadata
#' Raw and log ratio transformed abundances
#' 
#' #### Whole soil fungi
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
its_gulr_pfg <- 
  its_guab %>% 
  column_to_rownames(var = "field_name") %>% 
  decostand("rclr", MARGIN = 1) %>%
  rownames_to_column(var = "field_name") %>% 
  as_tibble() %>% 
  left_join(pfg, by = join_by(field_name)) %>% 
  left_join(gf_index, by = join_by(field_name)) %>% 
  left_join(sites %>% select(field_name, field_type, region, yr_since), by = join_by(field_name)) %>% 
  select(field_name, field_type, yr_since, region, everything())
its_guab_fa <- 
  its_guab %>% 
  rowwise() %>% 
  mutate(
    total = sum(c_across(where(is.numeric)), na.rm = TRUE),
    across(unidentified:unspecified_pathotroph, ~ if_else(total > 0, .x / total, 0))
  ) %>% 
  select(-total) %>%
  left_join(fa %>% select(-amf), by = join_by(field_name)) %>% 
  mutate(across(unidentified:unspecified_pathotroph, ~ .x * fungi_18.2)) %>% 
  select(field_name:unspecified_pathotroph)
its_gulr_fa <- 
  its_guab_fa %>% 
  column_to_rownames(var = "field_name") %>% 
  decostand("rclr", MARGIN = 1) %>%
  rownames_to_column(var = "field_name") %>% 
  as_tibble() %>% 
  left_join(sites %>% select(field_name, field_type, region, yr_since), by = join_by(field_name)) %>% 
  select(field_name, field_type, yr_since, region, everything())
#' 
#' #### AMF
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
amf_fmlr_pfg <- 
  amf_fmab %>% 
  column_to_rownames(var = "field_name") %>% 
  decostand("rclr", MARGIN = 2) %>% 
  rownames_to_column(var = "field_name") %>% 
  as_tibble() %>% 
  left_join(pfg, by = join_by(field_name)) %>% 
  left_join(gf_index, by = join_by(field_name)) %>% 
  left_join(sites %>% select(field_name, field_type, region, yr_since), by = join_by(field_name)) %>% 
  select(field_name, field_type, yr_since, region, everything())

#' ### Phyloseq databases
#' Only AMF needed here for UNIFRAC distance in PCoA
amf_ps <- phyloseq(
    otu_table(data.frame(spe$amf_avg_uni, row.names = 1) %>%
                decostand(method = "total", MARGIN = 2),
              taxa_are_rows = TRUE),
    tax_table(as.matrix(data.frame(spe_meta$amf, row.names = 2))),
    read.dna(root_path("otu_tables/18S/18S_sequences.fasta"), format = "fasta") %>%
        phyDat(type = "DNA") %>% dist.hamming() %>% NJ(),
    sample_data(sites %>% column_to_rownames(var = "field_name"))
)

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
its_div <- calc_div(spe$its_avg, sites)
#' 
#' ### Richness
#' Account for sequencing depth as a covariate
its_rich_lm <- lm(richness ~ sqrt(depth) + field_type, data = its_div)
#' Diagnostics
par(mfrow = c(2,2))
plot(its_rich_lm) # variance similar in groups
distribution_prob(its_rich_lm)
#' residuals distribution normal or close, response log
leveneTest(richness ~ field_type, data = its_div) %>% as.data.frame() %>% kable(format = "pandoc")
leveneTest(residuals(its_rich_lm) ~ its_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc")
#' Residuals/response distributions do not suggest the need for transformation.
#' Levene's p > 0.05 → fail to reject = variances can be considered equal.

#' Model results, group means, and post-hoc
summary(its_rich_lm)
#' Sequence depth is significant, less so than field type. Check relationship of depth and field type. 
its_div %>% 
    group_by(field_type) %>% 
    summarize(across(c(depth, richness), ~ round(mean(.x), 0))) %>% 
    kable(format = "pandoc")
#' Sequence depth isn't obviously related to field type. 
#' Proceed with means separation by obtaining estimated marginal means for field type.
#' Arithmetic means calculated in this case.
its_rich_em <- emmeans(its_rich_lm, ~ field_type, type = "response")
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
  ggplot(summary(its_rich_em), aes(x = field_type, y = emmean)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = emmean, ymax = upper.CL), width = 0, linewidth = lw) +
  geom_text(aes(y = upper.CL, label = c("a", "b", "b")), vjust = -1.5, family = "serif", size = 4) +
  labs(x = NULL, y = "Richness") +
  lims(y = c(0, 760)) +
  scale_fill_manual(values = c("gray", "black", "white")) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1))

#' 
#' ### Shannon's diversity
#' Account for sequencing depth as a covariate
its_shan_lm <- lm(shannon ~ sqrt(depth) + field_type, data = its_div)
#' Diagnostics
par(mfrow = c(2,2))
plot(its_shan_lm) # variance similar in groups 
distribution_prob(its_shan_lm)
#' residuals distribution most likely cauchy/normal; symmetric but long tails
#' residuals distribution normal or close, response gamma
leveneTest(shannon ~ field_type, data = its_div) %>% as.data.frame() %>% kable(format = "pandoc")
leveneTest(residuals(its_shan_lm) ~ its_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc")
#' Residuals distribution does not suggest the need for transformation.
#' Levene's p > 0.05 → fail to reject = variances can be considered equal.
#' Residuals distribution does not suggest the need for transformation. 
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

#' Model results, group means, and post-hoc
summary(its_shan_lm)
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
  geom_text(aes(y = upper.CL, label = c("a", "b", "b")), vjust = -1.5, family = "serif", size = 4) +
  labs(x = NULL, y = "Shannon diversity") +
  lims(y = c(0, 160)) +
  scale_fill_manual(values = c("gray", "black", "white")) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))

#' 
#' ## PLFA
plfa_lm <- lm(fungi_18.2 ~ field_type, data = fa)
par(mfrow = c(2,2))
plot(plfa_lm) # variance differs slightly in groups. Tails on qq plot diverge, lots of groups structure
distribution_prob(plfa_lm)
#' Residuals distribution fits normal, response normal-ish
leveneTest(residuals(plfa_lm) ~ fa$field_type) %>% as.data.frame() %>% kable(format = "pandoc") # No covariate, response and residuals tests equivalent
#' Residuals distribution does not suggest the need for transformation.
#' Levene's p > 0.05 → fail to reject = variances can be considered equal.

#' Model results, group means, and post-hoc, with arithmetic means from emmeans
summary(plfa_lm)
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
  labs(x = NULL, y = expression(PLFA~(nmol%*%g[soil]^-1))) +
  scale_fill_manual(values = c("gray", "black", "white")) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1.02))

#' 
#' ## Beta Diversity
#' Abundances were transformed by row proportions in sites before producing a distance matrix per
#' [McKnight et al.](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13115)
#+ its_ord
d_its <- spe$its_avg %>% 
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
    geom_point(aes(fill = field_type), size = sm_size, stroke = lw, shape = 21) +
    geom_text(aes(label = yr_since), size = yrtx_size, family = "serif", fontface = 2, color = "white") +
    geom_linerange(data = p_its_centers, aes(x = mean_Axis.1, y = mean_Axis.2, xmin = ci_l_Axis.1, xmax = ci_u_Axis.1), linewidth = lw) +
    geom_linerange(data = p_its_centers, aes(x = mean_Axis.1, y = mean_Axis.2, ymin = ci_l_Axis.2, ymax = ci_u_Axis.2), linewidth = lw) +
    geom_point(data = p_its_centers, aes(x = mean_Axis.1, y = mean_Axis.2, fill = field_type), size = lg_size, stroke = lw, shape = 21) +
    scale_fill_manual(values = c("gray", "black", "white")) +
    labs(
        x = paste0("Axis 1 (", mva_its$axis_pct[1], "%)"),
        y = paste0("Axis 2 (", mva_its$axis_pct[2], "%)")) +
    theme_ord +
    theme(legend.position = "none",
          plot.tag = element_text(size = 14, face = 1),
          plot.tag.position = c(0, 1.01))
    # guides(fill = guide_legend(position = "inside")) +
    # theme(legend.justification = c(0.03, 0.98))

#' 
#' ## Unified figure
#+ fig2_patchwork,warning=FALSE
fig2_ls <- (its_rich_fig / plot_spacer() / plfa_fig) +
    plot_layout(heights = c(1,0.01,1)) 
fig2 <- (fig2_ls | plot_spacer() | its_ord) +
    plot_layout(widths = c(0.35, 0.01, 0.64)) +
    plot_annotation(tag_levels = 'a') 
#+ fig2,warning=FALSE,fig.height=4,fig.width=6.5
fig2
#' **Fig 2.** Whole-soil fungal communities in **corn**, **restored**, and **remnant** prairie fields.
#' **a** OTU richness and **b** fungal biomass (PLFA) are shown as columns with 95 % CIs; lowercase
#' letters mark significant pairwise differences (P < 0.001).
#' **c** Principal-coordinate (PCoA) ordination of ITS-based (97 % OTU) community
#' distances: small points = sites, large circles = field-type centroids (error bars =
#' 95 % CI). Cornfields cluster apart from restored or remnant prairies (P < 0.01).
#' Numbers in black circles give years since restoration. Axis labels show the
#' percent variation explained. Colours/shading: corn = grey, restored = black,
#' remnant = white.

#+ fig2_save,warning=FALSE,echo=FALSE
ggsave(root_path("figs", "fig2.png"),
       plot = fig2,
       width = 6.5,
       height = 4,
       units = "in",
       dpi = 600)

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

#' Assemble explanatory variables and begin iterative selection process. 
#' Check the VIF for each explanatory variable to test for collinearity if model overfitting is 
#' detected. Then run forward selection in `dbrda()`. 
env_vars <- sites %>% 
  filter(field_type == "restored", !(region %in% "FL")) %>% 
  select(field_name, dist_axis_1) %>% # 90% on axis 1
  left_join(soil_micro_index, by = join_by(field_name)) %>% # 70% on first two axes
  left_join(soil_macro, by = join_by(field_name)) %>% 
  left_join(gf_index, by = join_by(field_name)) %>% # 92% on axis 1
  select(-starts_with("field_key"), -soil_micro_1, -K) %>% # soil_micro_1 removed based on initial VIF check
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
spe_its_wi_resto <- spe$its_avg %>% 
  filter(field_name %in% rownames(env_expl)) %>% 
  column_to_rownames(var = "field_name") %>% 
  decostand("total")

mod_null <- dbrda(spe_its_wi_resto ~ 1 + Condition(env_cov), data = env_expl, distance = "bray")
mod_full <- dbrda(spe_its_wi_resto ~ . + Condition(env_cov), data = env_expl, distance = "bray")
mod_step <- ordistep(mod_null, 
                     scope = formula(mod_full), 
                     direction = "forward", 
                     permutations = 1999, 
                     trace = FALSE)

#' ### Constrained Analysis Results
mod_step
(mod_glax <- anova(mod_step, permutations = 1999))
(mod_inax <- anova(mod_step, by = "axis", permutations = 1999))
(mod_r2   <- RsquareAdj(mod_step, permutations = 1999))
mod_step$anova %>% kable(, format = "pandoc")
#' Based on permutation tests with n=1999 permutations, the model shows a significant 
#' correlation between the site ordination on fungal communities
#' and the selected explanatory variables (p=0.001). The first two constrained axes are 
#' also significant (p<0.05). The selected variables explain $R^{2}_{\text{Adj}}$=21.3% of the community 
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
mod_scor_bp <- mod_scor$biplot %>% 
  data.frame() %>% 
  rownames_to_column(var = "envvar") %>% 
  mutate(
    origin = 0,
    m = dbRDA2 / dbRDA1, 
    d = sqrt(dbRDA1^2 + dbRDA2^2), 
    dadd = sqrt((max(dbRDA1)-min(dbRDA2))^2 + (max(dbRDA2)-min(dbRDA2))^2)*0.1,
    labx = ((d+dadd)*cos(atan(m)))*(dbRDA1/abs(dbRDA1)), 
    laby = ((d+dadd)*sin(atan(m)))*(dbRDA1/abs(dbRDA1)))

#+ fig6,warning=FALSE,fig.height=4,fig.width=4
fig6 <- 
  ggplot(mod_scor_site, aes(x = dbRDA1, y = dbRDA2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  geom_point(fill = "black", size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "serif", fontface = 2, color = "white") +
  geom_segment(data = mod_scor_bp, 
               aes(x = origin, xend = dbRDA1, y = origin, yend = dbRDA2), 
               arrow = arrow(length = unit(2, "mm"), type = "closed"),
               color = "gray20") +
  geom_text(data = mod_scor_bp, 
            aes(x = labx, y = laby, label = c("grass—forb\nindex", "pH")), 
            nudge_x = c(-0.2,-0.06), nudge_y = c(0.08,0.1),
            size = 3, color = "gray20") +
  labs(
    x = paste0("Constr. Axis 1 (", mod_pars_eig[1], "%)"),
    y = paste0("Constr. Axis 2 (", mod_pars_eig[2], "%)")) +
  theme_ord

#+ fig6_save,warning=FALSE,echo=FALSE
ggsave(
  root_path("figs", "fig6.png"),
  plot = fig6,
  width = 4,
  height = 4,
  units = "in",
  dpi = 600
)

#' 
#' # Arbuscular mycorrhizal fungi
# AMF ———————— ####
#' ## Diversity Indices
#+ amf_diversity
amf_div <- calc_div(spe$amf_avg, sites)
#' 
#' ### Richness
amf_rich_lm <- lm(richness ~ sqrt(depth) + field_type, data = amf_div)
#' Diagnostics
par(mfrow = c(2,2))
plot(amf_rich_lm) # variance similar in groups with an outlier
distribution_prob(amf_rich_lm)
#' Residuals distribution most likely normal, response bimodal (ignore)
leveneTest(richness ~ field_type, data = amf_div) %>% as.data.frame() %>% kable(format = "pandoc")
leveneTest(residuals(amf_rich_lm) ~ amf_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc")
#' Residuals/response distributions do not suggest the need for transformation.
#' Levene's p > 0.05 → fail to reject = variances can be considered equal.

#' Model results, group means, and post-hoc
summary(amf_rich_lm)
#' Sequencing depth not a significant predictor of amf richness
amf_rich_em <- emmeans(amf_rich_lm, ~ field_type, type = "response")
#' Results tables below show the emmeans summary of group arithmetic means and confidence intervals, 
#' and the post hoc contrast of richness among field types. 
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
  geom_text(aes(y = upper.CL, label = c("a", "b", "b")), vjust = -1.5, family = "serif", size = 4) +
  labs(x = NULL, y = "Richness") +
  lims(y = c(0, 75)) +
  scale_fill_manual(values = c("gray", "black", "white")) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1.02))

#' 
#' ### Shannon diversity
amf_shan_lm <- lm(shannon ~ sqrt(depth) + field_type, data = amf_div)
#' Diagnostics
par(mfrow = c(2,2))
plot(amf_shan_lm) 
#' Variance somewhat non-constant in groups, qqplot fit is poor, 
#' one leverage point (Cook's > 0.5)
distribution_prob(amf_shan_lm)
#' Residuals/response distributions most likely normal. 
leveneTest(shannon ~ field_type, data = amf_div) %>% as.data.frame() %>% kable(format = "pandoc")
leveneTest(residuals(amf_shan_lm) ~ amf_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc")
#' Residuals/response distributions do not suggest the need for transformation.
#' Covariate adds little added explanatory value.
#' Levene's p > 0.05 → fail to reject = variances can be considered equal.

#' Model results, group means, and post-hoc
summary(amf_shan_lm)
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
  geom_text(aes(y = upper.CL, label = c("a", "b", "b")), vjust = -1.5, family = "serif", size = 4) +
  labs(x = NULL, y = NULL) +
  lims(y = c(0, 32)) +
  scale_fill_manual(values = c("gray", "black", "white")) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(-0.05, 1))

#' 
#' ### Shannon's diversity figure
#' Includes both ITS and AMF for supplemental figure
#+ its_amf_shan_fig_pw,warning=FALSE
FigS1 <- (its_shan_fig | plot_spacer() | amf_shan_fig) +
    plot_layout(widths = c(1, 0.01, 1)) +
    plot_annotation(tag_levels = 'a')
#+ its_amf_shan_fig,warning=FALSE,fig.height=4,fig.width=7
FigS1
#+ amf_div_fig_save
ggsave(
    root_path("figs", "figS1.png"),
    plot = FigS1,
    height = 3,
    width = 6,
    units = "in",
    dpi = 600
)

#' 
#' ## NLFA
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
plot(nlfa_lm_log) # qqplot ok, one high leverage point in remnants
ncvTest(nlfa_lm_log) # p=0.16, null of constant variance not rejected

nlfa_glm  <- glm(amf ~ field_type, family = Gamma(link = "log"), data = fa)
nlfa_glm_diag <- glm.diag(nlfa_glm)
glm.diag.plots(nlfa_glm, nlfa_glm_diag) # qqplot shows strong fit; no leverage >0.5
performance::check_overdispersion(nlfa_glm) # not detected
#' Gamma glm is the best choice; no high-leverage point

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
  geom_text(aes(y = upper.CL, label = c("a", "b", "b")), vjust = -1.5, family = "serif", size = 4) +
  labs(x = NULL, y = expression(NLFA~(nmol%*%g[soil]^-1))) +
  scale_fill_manual(values = c("gray", "black", "white")) +
  lims(y = c(0, 75)) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1.02))

#' 
#' ## Beta Diversity
#' AMF (18S sequences)
#' 
#' Abundances were transformed by row proportions in sites before producing a distance matrix per
#' [McKnight et al.](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13115). 
#' Row proportions were calculated on the raw abundance data before creating the phyloseq data (see 
#' above). UNIFRAC distance matrix is created on the standardized abundance data.  
#+ amf_ord
d_amf <- UniFrac(amf_ps, weighted = TRUE)
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
         across(c(ci_l_Axis.2, ci_u_Axis.2), ~ mean_Axis.2 + .x))
amf_ord <- 
  ggplot(amf_ord_data, aes(x = Axis.1, y = Axis.2)) +
  geom_point(aes(fill = field_type), size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "serif", fontface = 2, color = "white") +
  geom_linerange(data = p_amf_centers, aes(x = mean_Axis.1, y = mean_Axis.2, xmin = ci_l_Axis.1, xmax = ci_u_Axis.1), linewidth = lw) +
  geom_linerange(data = p_amf_centers, aes(x = mean_Axis.1, y = mean_Axis.2, ymin = ci_l_Axis.2, ymax = ci_u_Axis.2), linewidth = lw) +
  geom_point(data = p_amf_centers, aes(x = mean_Axis.1, y = mean_Axis.2, fill = field_type), size = lg_size, stroke = lw, shape = 21) +
  scale_fill_manual(values = c("gray", "black", "white")) +
  labs(
    x = paste0("Axis 1 (", mva_amf$axis_pct[1], "%)"),
    y = paste0("Axis 2 (", mva_amf$axis_pct[2], "%)")) +
  theme_ord +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1.01))

#' ## Unified figure
#+ fig3_patchwork,warning=FALSE
fig3_ls <- (amf_rich_fig / plot_spacer() / nlfa_fig) +
  plot_layout(heights = c(1,0.01,1)) 
fig3 <- (fig3_ls | plot_spacer() | amf_ord) +
  plot_layout(widths = c(0.35, 0.01, 0.64)) +
  plot_annotation(tag_levels = 'a') 
#+ fig3,warning=FALSE,fig.height=4,fig.width=6.5
fig3
#' **Fig 3.** AMF communities in corn, restored, and remnant prairie fields. OTU richness **a**;
#' NLFA biomass **b** (95 % CI, letters = Tukey groups; P < 0.05 / 0.0001). PCoA of weighted‑UniFrac distances 
#' (18S, 97 % OTUs) **c**: small points = sites, large rings = field‑type centroids ±95 % CI. 
#' Numbers in points give years since restoration. Axes show % variance. Corn clusters apart from both 
#' prairie types (P < 0.01). Shading: corn grey, restored black, remnant white.

#+ fig3_save,warning=FALSE,fig.height=5,fig.width=7,echo=FALSE
ggsave(root_path("figs", "fig3.png"),
       plot = fig3,
       width = 6.5,
       height = 4,
       units = "in",
       dpi = 600)

#' ## AMF abundance in families
#' Display raw abundances in a table but separate means with log ratio transformed data
amf_fmab_ft <- 
  amf_fmab_pfg %>% 
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
distribution_prob(clar_lm)
leveneTest(Claroideoglomeraceae ~ field_type, data = amf_fmlr_pfg) %>% as.data.frame() %>% kable(format = "pandoc")
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

#' 
#' # Putative plant pathogens
# Putative plant pathogens ———————— ####
#' 
#' Retrieve pathogen sequence abundance
patho <- guildseq(spe$its_avg, spe_meta$its, "plant_pathogen")
#' 
#' ## Diversity Indices
#+ patho_diversity
patho_div <- calc_div(patho, sites)
#' 
#' ### Richness
#' Account for sequencing depth as a covariate
patho_rich_lm <- lm(richness ~ sqrt(depth) + field_type, data = patho_div)
#' Diagnostics
par(mfrow = c(2,2))
plot(patho_rich_lm) # variance similar in groups
distribution_prob(patho_rich_lm)
#' residuals distribution normal or close, response showing group divisions
leveneTest(richness ~ field_type, data = patho_div) %>% as.data.frame() %>% 
  kable(format = "pandoc", caption = "Response var in groups")
leveneTest(residuals(patho_rich_lm) ~ patho_div$field_type) %>% as.data.frame() %>% 
  kable(format = "pandoc", caption = "Residuals var in groups")
#' Residuals/response distributions do not suggest the need for transformation.
#' Levene's p > 0.05 → fail to reject = variances can be considered equal.

#' Model results, group means, and post-hoc
summary(patho_rich_lm)
#' Sequence depth is highly significant; richness doesn't vary in groups. 
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
  # lims(y = c(0, 760)) +
  scale_fill_manual(values = c("gray", "black", "white")) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1))

#' 
#' ### Shannon's diversity
#' Account for sequencing depth as a covariate
patho_shan_lm <- lm(shannon ~ sqrt(depth) + field_type, data = patho_div)
#' Diagnostics
par(mfrow = c(2,2))
plot(its_shan_lm) # variance similar in groups 
distribution_prob(patho_shan_lm)
#' residuals distribution most likely cauchy/normal; symmetric but long tails
#' response normal
leveneTest(shannon ~ field_type, data = patho_div) %>% as.data.frame() %>% kable(format = "pandoc")
leveneTest(residuals(patho_shan_lm) ~ patho_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc")
#' Residuals distribution does not suggest the need for transformation.
#' Levene's p > 0.05 → fail to reject = variances can be considered equal.

#' Model results, group means, and post-hoc
summary(patho_shan_lm)
#' Sequence depth is not a significant predictor of Shannon diversity, nor field type
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
  scale_fill_manual(values = c("gray", "black", "white")) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))

#' 
#' ## Abundance
#' Proportional sequence abundance as a proportion of biomass on a per-site basis.
#' Use log ratio transformed data, which handles composition bias so depth not 
#' needed as a covariate.
#' Note: previously used sequence abundance alone (data = its_gulr_pfg). Correcting for 
#' biomass differences across fields seems more appropriate.
patho_ab_lm <- lm(plant_pathogen ~ field_type, data = its_gulr_fa)
par(mfrow = c(2,2))
plot(patho_ab_lm) 



# NEED TO EXAMINE GLM MODEL or RLM, SEE SAPROTROPH SHANNONS FOR WORKFLOW TO HANDLE THE BIG OUTLIER
# Rememver that the saprotrophs abundance plot will likely not have negative values, might replace
# robust lm with gamma. 
# If you can go to gamma with all of them, do it, so you don't have one more model type to talk about. 






#' Variance differs slightly in groups. Tails on qq plot diverge. Row 16 an outlier.
#' This is the Lake Petite cornfield, very high. Results from pathogens being proportionally 
#' high at this site but PLFA being relatively low.
distribution_prob(patho_ab_lm)
#' Residuals distribution fits normal
leveneTest(residuals(patho_ab_lm) ~ patho_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc") # No covariate, response and residuals tests equivalent
#' Residuals distribution does not suggest the need for transformation.
#' Levene's p > 0.05 → fail to reject = variances can be considered equal.
#' Outlier will lean toward a Type II error, which is potentially conservative.

#' Model results, post hoc
summary(patho_ab_lm)
emmeans(patho_ab_lm, ~ field_type, type = "response")

#' Produce results with abundance corrected biomass (not log transformed) for a more intuitive 
#' figure. 
#' 



# FIG needs to use the emmeans data I think, particularly if I use gamma or rlm


#+ patho_fig,fig.width=4,fig.height=4
patho_ab_fig <- 
  its_gulr_fa %>% 
  # left_join(sites, by = join_by(field_name)) %>%
  group_by(field_type) %>% 
  summarize(
    mean = mean(plant_pathogen),
    ci = qnorm(0.975) * sd(plant_pathogen) / sqrt(n())
  ) %>% 
  ggplot(., aes(x = field_type, y = mean)) +
  geom_errorbar(aes(ymin = mean-ci, ymax = mean+ci), width = 0, linewidth = lw) +
  geom_point(aes(fill = field_type), shape = 21, size = sm_size) +
  labs(x = NULL, y = "Seq. prop. mass (LRT)") +
  scale_fill_manual(values = c("gray", "black", "white")) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1.02))

#' 
#' ## Beta Diversity
#' Abundances were transformed by row proportions in sites before producing a distance matrix per
#' [McKnight et al.](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13115)
#+ patho_ord
d_patho <- patho %>% 
  data.frame(row.names = 1) %>% 
  decostand("total") %>%
  vegdist("bray")
mva_patho <- mva(d = d_patho, env = sites, corr = "lingoes")
#+ patho_ord_results
mva_patho$dispersion_test
mva_patho$permanova
mva_patho$pairwise_contrasts[c(1,3,2), c(1,2,4,3,8)] %>% 
  kable(format = "pandoc", caption = "Pairwise permanova contrasts")
#' Lingoes correction was needed. Based on the homogeneity of variance test, the null hypothesis of equal variance among groups is 
#' accepted across all clusters and in pairwise comparison of clusters (both p>0.05), supporting the application of 
#' a PERMANOVA test. 
#' 
#' An effect of geographic distance (covariate) on pathogen communities was not supported. 
#' With geographic distance accounted for, the test variable 'field type' significantly explained 
#' variation in fungal communities, with a post-hoc test revealing that communities in corn fields differed from
#' communities in restored and remnant fields. 
#' 
#' Plotting results: 
patho_ord_data <- mva_patho$ordination_scores %>% mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant")))
p_patho_centers <- patho_ord_data %>% 
  group_by(field_type) %>% 
  summarize(across(starts_with("Axis"), list(mean = mean, ci_l = ci_l, ci_u = ci_u), .names = "{.fn}_{.col}"), .groups = "drop") %>% 
  mutate(across(c(ci_l_Axis.1, ci_u_Axis.1), ~ mean_Axis.1 + .x),
         across(c(ci_l_Axis.2, ci_u_Axis.2), ~ mean_Axis.2 + .x))
patho_ord <- 
  ggplot(patho_ord_data, aes(x = Axis.1, y = Axis.2)) +
  geom_point(aes(fill = field_type), size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "serif", fontface = 2, color = "white") +
  geom_linerange(data = p_patho_centers, aes(x = mean_Axis.1, y = mean_Axis.2, xmin = ci_l_Axis.1, xmax = ci_u_Axis.1), linewidth = lw) +
  geom_linerange(data = p_patho_centers, aes(x = mean_Axis.1, y = mean_Axis.2, ymin = ci_l_Axis.2, ymax = ci_u_Axis.2), linewidth = lw) +
  geom_point(data = p_patho_centers, aes(x = mean_Axis.1, y = mean_Axis.2, fill = field_type), size = lg_size, stroke = lw, shape = 21) +
  scale_fill_manual(values = c("gray", "black", "white")) +
  labs(
    x = paste0("Axis 1 (", mva_patho$axis_pct[1], "%)"),
    y = paste0("Axis 2 (", mva_patho$axis_pct[2], "%)")) +
  theme_ord +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1.01))
# guides(fill = guide_legend(position = "inside")) +
# theme(legend.justification = c(0.03, 0.98))

#' 
#' ## Unified figure
#+ fig4_patchwork,warning=FALSE
fig4_ls <- (patho_rich_fig / plot_spacer() / patho_ab_fig) +
  plot_layout(heights = c(1,0.01,1)) 
fig4 <- (fig4_ls | plot_spacer() | patho_ord) +
  plot_layout(widths = c(0.35, 0.01, 0.64)) +
  plot_annotation(tag_levels = 'a') 
#+ fig4,warning=FALSE,fig.height=4,fig.width=6.5
fig4
#' **Fig 4.** Putative plant pathogen communities in **corn**, **restored**, and **remnant** prairie fields.
#' **a** OTU richness and **b** sequence abundance are shown as columns with 95 % CIs.
#' **c** Principal-coordinate (PCoA) ordination of ITS-based (97 % OTU) community
#' distances: small points = sites, large circles = field-type centroids (error bars =
#' 95 % CI). Cornfields cluster apart from restored or remnant prairies (P < 0.01).
#' Numbers in black circles give years since restoration. Axis labels show the
#' percent variation explained. Colours/shading: corn = grey, restored = black,
#' remnant = white.

#+ fig4_save,warning=FALSE,fig.height=5,fig.width=7,echo=FALSE
ggsave(root_path("figs", "fig4.png"),
       plot = fig4,
       width = 6.5,
       height = 4,
       units = "in",
       dpi = 600)

#' ## Pathogen Indicator Species
#' Use as a tool to find species for discussion. Unbalanced design and bias to agricultural soil
#' research may make the indicator stats less appropriate for other use.
patho_ind <- inspan(spe$its_avg, spe_meta$its, "plant_pathogen", sites)
patho_ind %>% 
  select(-otu_num, -primary_lifestyle, -p.value) %>% 
  arrange(field_type, p_val_adj) %>% 
  kable(format = "pandoc", caption = paste("Indicator species analysis, plant pathogens"))

patho_abund_ft <- 
  guildseq(spe$its_avg, spe_meta$its, "plant_pathogen") %>% 
  pivot_longer(starts_with("otu"), names_to = "otu_num", values_to = "abund") %>% 
  left_join(spe_meta$its, by = join_by(otu_num)) %>% 
  left_join(sites %>% select(field_name, field_type), by = join_by(field_name)) %>% 
  select(-otu_ID, -otu_num, -primary_lifestyle, -field_name) %>%
  filter(species != "unidentified", abund > 0) %>% 
  pivot_wider(names_from = "field_type", values_from = "abund", values_fn = ~ round(mean(.x), 1), values_fill = 0) %>% 
  select(phylum:species, corn, restored, remnant) %>% 
  rowwise() %>% 
  mutate(avg = mean(c_across(corn:remnant)) %>% round(., 1)) %>% 
  arrange(-corn)
patho_abund_ft %>% 
  filter(avg >= 1) %>%
  kable(format = "pandoc", caption = "Named pathogen species and abundances in field types\n(Mean abundance >= 1 shown)")

#' ## Pathogen—Plant Correlations
#' Whole-soil fungi correlated with grass and forbs; investigate whether pathogens do specifically.
ggplot(its_guab_pfg %>% filter(field_type == "restored", region != "FL"), aes(x = gf_index, y = plant_pathogen)) +
  geom_smooth(method = "lm") +
  geom_point(shape = 1, size = 7) +
  geom_text(aes(label = yr_since)) +
  labs(x = "Index: C4_grass <—> Forb abundance", y = "Sequence abundance, plant pathogens") +
  theme_classic()

#' Model the relationship
gf_patho_lm <- lm(plant_pathogen ~ gf_index, data = its_gulr_pfg %>% filter(field_type == "restored", region != "FL"))
#' Diagnostics
par(mfrow = c(2,2))
plot(gf_patho_lm) 
#' Some residual structure, but a smooth qq fit. 
#' Minor leverage with point 4 pulls the slope to more level, risking a type II error rather than type I. 
#' Model is still significant with point 4 removed. 
distribution_prob(gf_patho_lm)
#' Response and residuals normal, no transformations warranted and linear model appropriate.
summary(gf_patho_lm)

#' 
#' # Putative saprotrophs
# Putative saprotrophs ———————— ####
#' 
#' Retrieve pathogen sequence abundance
sapro <- guildseq(spe$its_avg, spe_meta$its, "saprotroph")
#' 
#' ## Diversity Indices
#+ sapro_diversity
sapro_div <- calc_div(sapro, sites)
#' 
#' ### Richness
#' Account for sequencing depth as a covariate
sapro_rich_lm <- lm(richness ~ sqrt(depth) + field_type, data = sapro_div)
#' Diagnostics
par(mfrow = c(2,2))
plot(sapro_rich_lm) # heavy residual structure, poor qq alignment, borderline leverage
distribution_prob(sapro_rich_lm)
#' residuals distribution normal or close, response showing group divisions and count 
#' overdispersion (binomial family)
leveneTest(richness ~ field_type, data = sapro_div) %>% as.data.frame() %>% 
  kable(format = "pandoc", caption = "Response var in groups")
leveneTest(residuals(sapro_rich_lm) ~ sapro_div$field_type) %>% as.data.frame() %>% 
  kable(format = "pandoc", caption = "Residuals var in groups")
#' Levene's p > 0.05 → fail to reject = variances can be considered equal.
#' LM-log transformation did not improve fit (not shown). Despite lack of unequal variance,
#' does gamma glm improve other diagnostics?
sapro_rich_glm  <- glm(richness ~ sqrt(depth) + field_type, family = Gamma(link = "log"), data = sapro_div)
sapro_glm_diag <- glm.diag(sapro_rich_glm)
glm.diag.plots(sapro_rich_glm, sapro_glm_diag) 
#' scale-location much more uniform but still a little heavy tails,
#' qqplot shows much better fit; leverage point now ~0.3
performance::check_overdispersion(sapro_rich_glm) # not detected
sapro_glm_sim <- simulateResiduals(sapro_rich_glm)
plot(sapro_glm_sim) # DHARMa passes all tests
#' Gamma glm is the best choice; no high-leverage point

#' Model results, group means, and post-hoc
summary(sapro_rich_glm)
#' Sequence depth is significant; richness doesn't vary in groups. View trend:
sapro_div %>% 
  group_by(field_type) %>% 
  summarize(avg_richness = mean(richness),
            avg_depth = mean(depth), .groups = "drop") %>% 
  mutate(across(starts_with("avg"), ~ round(.x, 1))) %>% 
  kable(format = "pandoc")
#' depth and richness in a negative relationship; depth not driving richness...
 
#' Calculate confidence intervals for figure.
#' Arithmetic means calculated in this case, back-transformed.
sapro_rich_em <- emmeans(sapro_rich_glm, ~ field_type, type = "response")
#+ sapro_rich_em_summary,echo=FALSE
kable(summary(sapro_rich_em), 
      format = "pandoc", 
      caption = "Confidence level used: 0.95")
#+ sapro_rich_em_posthoc,echo=FALSE
kable(pairs(sapro_rich_em), 
      format = "pandoc", 
      caption = "P value adjustment: tukey method for comparing a family of 3 estimates")
#+ sapro_richness_fig,fig.width=4,fig.height=4
sapro_rich_fig <- 
  ggplot(summary(sapro_rich_em), aes(x = field_type, y = response)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = response, ymax = upper.CL), width = 0, linewidth = lw) +
  labs(x = NULL, y = "Richness") +
  # lims(y = c(0, 760)) +
  scale_fill_manual(values = c("gray", "black", "white")) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1))

#' 
#' ### Shannon's diversity
#' Account for sequencing depth as a covariate
sapro_shan_lm <- lm(shannon ~ sqrt(depth) + field_type, data = sapro_div)
#' Diagnostics
par(mfrow = c(2,2))
plot(sapro_shan_lm) # variance similar in groups 
distribution_prob(sapro_shan_lm)
#' residuals distribution most likely cauchy/normal; symmetric but long tails
#' response non-normal, check variance in groups though
leveneTest(shannon ~ field_type, data = sapro_div) %>% as.data.frame() %>% kable(format = "pandoc")
leveneTest(residuals(sapro_shan_lm) ~ sapro_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc")
#' Residuals distribution does not suggest the need for transformation.
#' Levene's p > 0.05 → fail to reject = variances can be considered equal.

#' Model results, group means, and post-hoc
summary(sapro_shan_lm)
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
  # lims(y = c(0, 160)) +
  scale_fill_manual(values = c("gray", "black", "white")) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))

#' 
#' ## Sequence abundance
#' Use log ratio transformed data to detrend compositional data before differential analysis. 
#' Transformation also corrects for different sample depths; no covariate needed. 
sapro_ab_lm <- lm(saprotroph ~ field_type, data = its_gulr_pfg)
par(mfrow = c(2,2))
plot(sapro_ab_lm) 
#' variance differs slightly in groups. Tails on qq plot diverge. Variance unequal but not bad, 
#' two leverage points.
distribution_prob(sapro_ab_lm)
#' Residuals distribution fits normal, so do residuals
leveneTest(residuals(sapro_ab_lm) ~ sapro_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc") # No covariate, response and residuals tests equivalent
#' Residuals distribution does not suggest the need for transformation.
#' Levene's p > 0.05 → fail to reject = variances can be considered equal (aka homoscedastic 
#' by group).
#'  
#' Try a model that 
#' handles the heavy tails on the residual distribution better? Transformed data has values < 0, 
#' so gamma glm is not available. Try robust regression (robust M-estimator with bisquare ψ).
sapro_ab_rlm <- rlm(saprotroph ~ field_type, data = its_gulr_pfg)
par(mfrow = c(2,2))
plot(sapro_ab_rlm) 
#' Minor visual improvements to fit 
distribution_prob(sapro_ab_rlm)
#' Residuals distribution fits normal, so do residuals
leveneTest(residuals(sapro_ab_rlm) ~ sapro_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc")
broom::tidy(sapro_ab_lm, conf.int = TRUE)
broom::tidy(sapro_ab_rlm, conf.int = TRUE) # rlm produces a better fit
#' Points with a high Cook's statistic remain high, but have less leverage. 
#' Iteratively re-weighted least squares (rlm) lowers leverage for the two extreme sites, 
#' but their Cook’s distances remain ≈ 0.5 because the robust fit reduces the 
#' residual scale estimate, inflating studentised residuals. Coefficient estimates 
#' changed little, confirming results are not driven by those points.
#' 
#' Produce model results, group means, and post-hoc, with arithmetic means from emmeans
summary(sapro_ab_rlm)
sapro_ab_em <- emmeans(sapro_ab_rlm, ~ field_type, type = "response", vcov. = vcov(sapro_ab_rlm))
#+ sapro_ab_em_summary,echo=FALSE
kable(summary(sapro_ab_em),
      format = "pandoc",
      caption = "Confidence level used: 0.95")
#+ sapro_ab_em_posthoc,echo=FALSE
kable(pairs(sapro_ab_em),
      format = "pandoc",
      caption = "P value adjustment: tukey method for comparing a family of 3 estimates")
#' Aligning statistical result with visual appearance requires a dot line plot
#+ sapro_ab_fig,fig.width=4,fig.height=4
sapro_ab_fig <- 
  # sapro_div %>% 
  # group_by(field_type) %>% 
  # summarize(seq_ab = mean(depth), upper.CL = seq_ab + ci_u(depth), .groups = "drop") %>% 
  ggplot(summary(sapro_ab_em), aes(x = field_type, y = emmean)) +
  # geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0, linewidth = lw) +
  geom_point(aes(fill = field_type), shape = 21, size = sm_size) +
  geom_text(aes(y = asymp.UCL, label = c("a", "b", "c")), vjust = -1.5, family = "serif", size = 4) + 
  labs(x = NULL, y = "Seq. abund. (LRT)") +
  scale_fill_manual(values = c("gray", "black", "white")) +
  lims(y = c(-0.6, 0.7)) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1.02))

#' 
#' ## Beta Diversity
#' Abundances were transformed by row proportions in sites before producing a distance matrix per
#' [McKnight et al.](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13115)
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
#' An effect of geographic distance (covariate) on pathogen communities was not supported. 
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
  geom_point(aes(fill = field_type), size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "serif", fontface = 2, color = "white") +
  geom_linerange(data = p_sapro_centers, aes(x = mean_Axis.1, y = mean_Axis.2, xmin = ci_l_Axis.1, xmax = ci_u_Axis.1), linewidth = lw) +
  geom_linerange(data = p_sapro_centers, aes(x = mean_Axis.1, y = mean_Axis.2, ymin = ci_l_Axis.2, ymax = ci_u_Axis.2), linewidth = lw) +
  geom_point(data = p_sapro_centers, aes(x = mean_Axis.1, y = mean_Axis.2, fill = field_type), size = lg_size, stroke = lw, shape = 21) +
  scale_fill_manual(values = c("gray", "black", "white")) +
  labs(
    x = paste0("Axis 1 (", mva_sapro$axis_pct[1], "%)"),
    y = paste0("Axis 2 (", mva_sapro$axis_pct[2], "%)")) +
  theme_ord +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1.01))
# guides(fill = guide_legend(position = "inside")) +
# theme(legend.justification = c(0.03, 0.98))

#' 
#' ## Unified figure
#+ fig5_patchwork,warning=FALSE
fig5_ls <- (sapro_rich_fig / plot_spacer() / sapro_ab_fig) +
  plot_layout(heights = c(1,0.01,1)) 
fig5 <- (fig5_ls | plot_spacer() | sapro_ord) +
  plot_layout(widths = c(0.35, 0.01, 0.64)) +
  plot_annotation(tag_levels = 'a') 
#+ fig5,warning=FALSE,fig.height=4,fig.width=6.5
fig5
#' **Fig 4.** Putative plant pathogen communities in **corn**, **restored**, and **remnant** prairie fields.
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

#' ## Saprotroph Indicator Species
sapro_ind <- inspan(spe$its_avg, spe_meta$its, "saprotroph", sites)
sapro_ind_table <- 
  sapro_ind %>% 
  select(-otu_num, -primary_lifestyle, -p.value) %>%
  arrange(field_type, p_val_adj)
sapro_ind_table[1:20, ] %>% 
  kable(format = "pandoc", caption = paste("Indicator species analysis, saprotrophs\n(First 20 rows shown)"))
sapro_abund_ft <- 
  guildseq(spe$its_avg, spe_meta$its, "saprotroph") %>% 
  pivot_longer(starts_with("otu"), names_to = "otu_num", values_to = "abund") %>% 
  left_join(spe_meta$its, by = join_by(otu_num)) %>% 
  left_join(sites %>% select(field_name, field_type), by = join_by(field_name)) %>% 
  select(-otu_ID, -otu_num, -primary_lifestyle, -field_name) %>%
  filter(species != "unidentified", abund > 0) %>% 
  pivot_wider(names_from = "field_type", values_from = "abund", values_fn = ~ round(mean(.x), 1), values_fill = 0) %>% 
  select(phylum:species, corn, restored, remnant) %>% 
  rowwise() %>% 
  mutate(avg = mean(c_across(corn:remnant)) %>% round(., 1)) %>% 
  arrange(-corn)
sapro_abund_ft %>% 
  filter(avg >= 10) %>% 
  kable(format = "pandoc", caption = "Named saprotroph species and abundances in field types\n(Mean abundance >= 10 shown)")

