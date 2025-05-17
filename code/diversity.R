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
#'  1. PCoA of Bray (ITS) or UNIFRAC (18S) distances,  
#'  2. homogeneity test diagnostics 
#'  3. PERMANOVA (+ pairwise) 
#' Cartesian inter‑site distance enters models as a covariate per [Redondo 2020](https://doi.org/10.1093/femsec/fiaa082).

#' # Packages and libraries
# Libraries ———————— ####
#+ packages,message=FALSE
packages_needed <- c(
  "colorspace", "emmeans", "gridExtra", "knitr", "tidyverse", "vegan",
  "rprojroot", "phyloseq", "ape", "phangorn", "geosphere", "conflicted",
  "ggpubr", "patchwork", "car", "performance", "broom", "boot"
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
sites <- read_csv(paste0(getwd(), "/clean_data/sites.csv"), show_col_types = FALSE) %>% 
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
sites$dist_axis_1 = field_dist_pcoa$vectors[, 1]

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

#' 
#' ### Wrangle species and metadata
#' Raw and log ratio transformed abundances
its_guab <- 
  spe$its_avg %>% 
  pivot_longer(starts_with("otu"), names_to = "otu_num", values_to = "abund") %>% 
  left_join(spe_meta$its %>% select(otu_num, primary_lifestyle), by = join_by(otu_num)) %>% 
  group_by(field_name, primary_lifestyle) %>% summarize(abund = sum(abund), .groups = "drop") %>% 
  arrange(field_name, -abund) %>% 
  pivot_wider(names_from = "primary_lifestyle", values_from = "abund") %>% 
  select(field_name, unidentified, saprotroph, plant_pathogen, everything())
its_gulr <- 
  its_guab %>% 
  column_to_rownames(var = "field_name") %>% 
  decostand("rclr", MARGIN = 2) %>%
  rownames_to_column(var = "field_name") %>% 
  as_tibble() %>% 
  left_join(sites %>% select(field_name, field_type), by = join_by(field_name)) %>% 
  select(field_name, field_type, everything())

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
#' ## Fatty Acids: Biomass
#' Use only 18.2 for soil fungi
fa <- read_csv(file.path(getwd(), "clean_data/plfa.csv"), show_col_types = FALSE) %>% 
    rename(fungi_18.2 = fa_18.2) %>% 
    select(field_name, fungi_18.2, amf) %>%
    left_join(
        sites %>% select(field_name, field_type),
        by = join_by(field_name)
    )

#' 
#' # Functions
# Functions ———————— ####
#' 
#' ## Alpha diversity calculations
#' Returns a dataframe of alpha diversity (richness, Shannon's) for analysis and plotting.
#+ calc_diversity_function
calc_div <- function(spe) {
    div_data <- 
        spe %>% 
        rowwise() %>% 
        mutate(
            depth = sum(c_across(starts_with("otu"))),
            richness = sum(c_across(starts_with("otu")) > 0),
            shannon = exp(diversity(c_across(starts_with("otu"))))
        ) %>% 
        select(-starts_with("otu")) %>% 
        as_tibble() %>% 
        left_join(sites %>% select(field_type, field_name), by = join_by(field_name)) %>% 
        mutate(across(starts_with("field"), ~ factor(.x, ordered = FALSE)))
    
    return(div_data)
    
}

#' ## Confidence intervals
#' Calculate upper and lower confidence intervals with alpha=0.05
#+ ci_function
ci_u <- function(x) {(sd(x) / sqrt(length(x))) * qnorm(0.975)}
ci_l <- function(x) {(sd(x) / sqrt(length(x))) * qnorm(0.025)}

#' 
#' ## Multivariate analysis
#' Ordination → dispersion check → global & pairwise PERMANOVA → envfit.
#' Args: *d* dist, *env* metadata, *corr* PCoA correction, *nperm* permutations.
#' 
#+ mva_function
mva <- function(d, env=sites, corr="none", nperm=1999) {
    # Ordination
    p <- pcoa(d, correction = corr)
    p_vals <- data.frame(p$values) %>% 
        rownames_to_column(var = "Dim") %>% 
        mutate(Dim = as.integer(Dim))
    p_eig <- p_vals[1:2, grep("Rel", colnames(p_vals))] %>% round(., 3) * 100
    p_vec <- data.frame(p$vectors)
    p_sco <-
        p_vec[, 1:2] %>%
        rownames_to_column(var = "field_name") %>%
        left_join(env, by = join_by(field_name))
    p_fit <- envfit(
        p_vec ~ yr_since, 
        data.frame(env, row.names = 1), 
        choices = c(1,2),
        na.rm = TRUE,
        permutations = nperm)
    p_fit_sco <- scores(p_fit, display = "bp")
    # Test multivariate dispersions
    disper <- betadisper(d, env$field_type, bias.adjust = TRUE)
    mvdisper <- permutest(disper, pairwise = TRUE, permutations = nperm)
    # Global PERMANOVA
    gl_permtest <- adonis2(
        d ~ dist_axis_1 + field_type,
        data = env,
        permutations = nperm,
        by = "terms")
    # Pairwise PERMANOVA
    group_var <- as.character(env$field_type)
    groups <- as.data.frame(t(combn(unique(group_var), m = 2)))
    contrasts <- data.frame(
        group1 = groups$V1,
        group2 = groups$V2,
        R2 = NA,
        F_value = NA,
        df1 = NA,
        df2 = NA,
        p_value = NA
    )
    for (i in seq(nrow(contrasts))) {
        group_subset <-
            group_var == contrasts$group1[i] |
            group_var == contrasts$group2[i]
        # Contrast matrices for Unifrac and Bray distance aren't compatible
        contrast_matrix <- as.matrix(d)[group_subset, group_subset]
        fit <- adonis2(
            contrast_matrix ~ env$dist_axis_1[group_subset] + group_var[group_subset],
            permutations = nperm,
            by = "terms")
        # Prepare contrasts table
        contrasts$R2[i] <- round(fit[grep("group_var", rownames(fit)), "R2"], digits = 3)
        contrasts$F_value[i] <- round(fit[grep("group_var", rownames(fit)), "F"], digits = 3)
        contrasts$df1[i] <- fit[grep("group_var", rownames(fit)), "Df"]
        contrasts$df2[i] <- fit[grep("Residual", rownames(fit)), "Df"]
        contrasts$p_value[i] <- fit[grep("group_var", rownames(fit)), 5]
    }
    contrasts$p_value_adj <- p.adjust(contrasts$p_value, method = "fdr") %>% round(., 4)
    
    # Results
    out <- list(
        correction_note    = p$note,
        ordination_values  = p_vals[1:10, ],
        axis_pct           = p_eig,
        ordination_scores  = p_sco,
        dispersion_test    = mvdisper,
        permanova          = gl_permtest,
        pairwise_contrasts = contrasts,
        vector_fit_result  = p_fit,
        vector_fit_scores  = p_fit_sco
    )
    
    return(out)
    
}

#' ## Model distribution probabilities
#' Probable distributions of response and residuals. Package performance prints
#' javascript which doesn't render on github documents.
#+ distribution_prob_function
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

#' ## Filter spe to a guild
#' Create samp-spe matrix of sequence abundance in a guild
#+ guildseq_function
guildseq <- function(spe, guild) {
  guab <- 
    spe %>% 
    pivot_longer(starts_with("otu"), names_to = "otu_num", values_to = "abund") %>% 
    left_join(spe_meta$its %>% select(otu_num, primary_lifestyle), by = join_by(otu_num)) %>% 
    filter(primary_lifestyle == guild) %>% 
    select(-primary_lifestyle) %>% 
    pivot_wider(names_from = otu_num, values_from = abund)
  return(guab)
}

#' 
#' # Whole Soil Fungi
# Whole soil fungi ———————— ####
#' 
#' ## Diversity Indices

#+ its_diversity
its_div <- calc_div(spe$its_avg)
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
#' Sequence depth is not a significant predictor of Shannon diversity
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
mva_its <- mva(d = d_its)
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

#+ fig2_save,warning=FALSE,fig.height=5,fig.width=7,echo=FALSE
ggsave(root_path("figs", "fig2.png"),
       plot = fig2,
       width = 6.5,
       height = 4,
       units = "in",
       dpi = 600)

#' 
#' # Arbuscular mycorrhizal fungi
# AMF ———————— ####
#' ## Diversity Indices
#+ amf_diversity
amf_div <- calc_div(spe$amf_avg)
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
#' Residuals distribution does not suggest the need for transformation.
#' Levene's p ~ 0.05 → reject the null of equal variance. Check CV in groups.
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
mva_amf <- mva(d = d_amf, corr = "lingoes")
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

#' 
#' # Putative plant pathogens
# Putative plant pathogens ———————— ####
#' 
#' Retrieve pathogen sequence abundance
patho <- guildseq(spe$its_avg, "plant_pathogen")
#' 
#' ## Diversity Indices
#+ patho_diversity
patho_div <- calc_div(patho)
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
#' Sequence depth is highly significant; richness doesn't vary in groups 
#' Calculate confidence intervals for figure.
#' Arithmetic means calculated in this case.
patho_rich_em <- emmeans(patho_rich_lm, ~ field_type, type = "response")
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
#' ## Sequence abundance
#' Use log ratio transformed data, which handles composition bias so depth not 
#' needed as a covariate.
patho_ab_lm <- lm(plant_pathogen ~ field_type, data = its_gulr)
par(mfrow = c(2,2))
plot(patho_ab_lm) # variance differs slightly in groups. Tails on qq plot diverge
distribution_prob(patho_ab_lm)
#' Residuals distribution fits normal
leveneTest(residuals(patho_ab_lm) ~ patho_div$field_type) %>% as.data.frame() %>% kable(format = "pandoc") # No covariate, response and residuals tests equivalent
#' Residuals distribution does not suggest the need for transformation.
#' Levene's p > 0.05 → fail to reject = variances can be considered equal.

#' Model results, group means, and post-hoc, with arithmetic means from emmeans
summary(patho_ab_lm)
patho_ab_em <- emmeans(patho_ab_lm, ~ field_type, type = "response")
#+ patho_fig,fig.width=4,fig.height=4
patho_ab_fig <- 
  patho_div %>% 
  group_by(field_type) %>% 
  summarize(seq_ab = mean(depth), upper.CL = seq_ab + ci_u(depth), .groups = "drop") %>% 
  ggplot(aes(x = field_type, y = seq_ab)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = seq_ab, ymax = upper.CL), width = 0, linewidth = lw) +
  labs(x = NULL, y = "Sequence abundance") +
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
mva_patho <- mva(d = d_patho, corr = "lingoes")
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
