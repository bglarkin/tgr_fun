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
#' Main results analysis and figures
#' - Biomass: PLFA/NLFA
#' - OTU richness and diversity
#' - Beta diversity
#' 
#' ## Alpha diversity
#' Sequence abundance and fungal biomass are used to evaluate alpha diversity. Sequence abundances were determined
#' in 97% similar OTUs for ITS and 18S genes. Site averages were used because sites are the replicate in this study. 
#' 
#' This report presents visualizations and results of tests contrasting OTU richness and diversity, and 
#' fungal biomass, among field types. 
#' Model selection consisted of visual inspection of diagnostic plots 
#' to choose the model that best minimized heteroscedasticity and the residuals' departure from a 
#' normal distribution. Potential outlier points were examined for leverage. 
#' 
#' Since OTU counts considered here lie far from zero and produce relatively normal or log-normal distributions, 
#' the central limit theorem justifies that richness be interpreted as a continuous variable. Further model selection 
#' used linear models and then generalized linear models, which provided the best fits with either Gaussian 
#' error distributions and log-link functions. Sequence depth was square-root transformed to improve fit 
#' as in [Bálint et al. 2015](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.13018).
#' 
#' Means separation among field types was accomplished using least-squares means, implemented in package emmeans. 
#' 
#' 
#' ## Beta diversity
#' Beta diversity is tested with multivariate analysis, conducted using species abundances averaged
#' in replicate sites. Pooling abundances from multiple subsamples has been shown to increase 
#' recovery of microbial richness [Song et al. 2015](https://doi.org/10.1371/journal.pone.0127234) The analysis flow is:
#' 
#' 1. Ordination of sites (PCoA)
#'  1. ITS-based abundances are standardized to proportions within sites and converted into a bray-curtis distance matrix 
#'  [McKnight et al. 2019](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13115).
#'  1. 18S-based abundances are standardized to proportions within sites and converted into a UNIFRAC distance matrix
#' 1. Clustering by field type (corn, restored, remnant)
#'  1. Multivariate homogeneity test
#'  1. Global PERMANOVA
#'  1. Pairwise PERMANOVA
#' 1. Post-hoc fitting of restored site age (Envfit)
#' 
#' Beta diverstiy could depend on intersite distance, which may limit propagule dispersal [Redondo et al. 2020](https://doi.org/10.1093/femsec/fiaa082).
#' We use cartesian intersite distance as a covariate in statistical tests to account for this. 
#' 
#' 




#' # Packages and libraries
# Libraries ———————— ####
packages_needed = c(
    "colorspace",
    "emmeans",
    "gridExtra",
    "knitr",
    "tidyverse",
    "vegan",
    "rprojroot",
    "phyloseq", 
    "ape", 
    "phangorn", 
    "geosphere", 
    "conflicted",
    "ggpubr",
    "patchwork",
    "car",
    "performance",
    "broom"
)
packages_installed = packages_needed %in% rownames(installed.packages())
#+ packages,message=FALSE
if (any(!packages_installed)) {
    install.packages(packages_needed[!packages_installed])
}
#+ libraries,message=FALSE
for (i in 1:length(packages_needed)) {
    library(packages_needed[i], character.only = T)
}

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
sites <- read_csv(paste0(getwd(), "/clean_data/sites.csv"), show_col_types = FALSE) 

#' ### Wrangle site metadata
#' Intersite geographic distance will be used as a covariate in clustering. 
#' Raw coordinates in data file aren't distances; convert to distance matrix and summarize with PCoA
field_dist <- as.dist(distm(sites[, c("long", "lat")], fun = distHaversine))
field_dist_pcoa <- pcoa(field_dist)
field_dist_pcoa$values[c(1,2), c(1,2)] %>% 
    kable(format = "pandoc")
#' First axis of geographic distance PCoA explains 91% of the variation among sites. 
sites$dist_axis_1 = field_dist_pcoa$vectors[, 1]

#' ## Sites-species tables
#' List *spe* holds average sequence abundances for the top 6 samples per field. 
#' CSV files were produced in `process_data.R`
spe <- list(
    its_avg   = read_csv(root_path("clean_data/spe_ITS_avg.csv"), 
                         show_col_types = FALSE),
    amf_avg   = read_csv(root_path("clean_data/spe_18S_avg.csv"), 
                         show_col_types = FALSE),
    amf_avg_uni = read_delim(root_path("otu_tables/18S/18S_avg_4unifrac.tsv"),
                         show_col_types = FALSE)
)
#' 
#' ## Microbial species metadata
spe_meta <- list(
    its = read_csv(root_path("clean_data/spe_ITS_metadata.csv"),
                   show_col_types = FALSE),
    amf = read_csv(root_path("clean_data/spe_18S_metadata.csv"),
                   show_col_types = FALSE)
) %>% map(. %>% mutate(across(everything(), ~ replace_na(., "unidentified"))))

#' ### Phyloseq databases
#' Sequence abundance data aren't standardized here; 
#' they will be after the data are subsetted by guilds or taxonomic ranks.
#' #### Whole soil fungi
its_ps <- phyloseq(
    otu_table(
        data.frame(spe$its_avg, row.names = 1),
        taxa_are_rows = FALSE
    ),
    sample_data(
        sites %>% column_to_rownames(var = "field_name")
    ),
    tax_table(
        data.frame(spe_meta$its, row.names = 1) %>% 
            as.matrix()
    )
)
#' 
#' #### AMF
amf_avg_ps <- phyloseq(
    otu_table(
        data.frame(spe$amf_avg_uni, row.names = 1) %>%
            decostand(method = "total", MARGIN = 2),
        taxa_are_rows = TRUE
    ),
    read.dna(
        root_path("otu_tables/18S/18S_sequences.fasta"),
        format = "fasta"
    ) %>%
        phyDat(type = "DNA") %>%
        dist.hamming() %>%
        NJ(),
    sample_data(
        sites %>% column_to_rownames(var = "field_name")
    )
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
#' Perform ordination, diagnostics, PERMANOVA, and post-hoc
#' 
#' 1. d = object of class dist
#' 1. s = species abundance matrix
#' 1. env = site metadata
#' 1. corr = any correction used in `pcoa()`
#' 1. distype = distance method, unifrac or other
#' 1. nperm = number of permutations for tests
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
#' 








#' # Whole Soil Fungi
# Whole soil fungi ———————— ####
#' ## Diversity Indices

#+ its_diversity
its_div <- calc_div(spe$its_avg)
#' 
#' ### Richness


#' Account for sequencing depth as a covariate
its_rich_lm <- lm(richness ~ sqrt(depth) + field_type, data = its_div)
par(mfrow = c(2,2))
plot(its_rich_lm) # variance similar in groups 
shapiro.test(residuals(its_rich_lm)) # p=0.17
performance::check_distribution(its_rich_lm) # residuals distribution most likely cauchy; symmetric but long tails
#' Based on visual inspection and normality of residuals, log transformation or link function may not be warranted. 
#' Examine the CV in groups; if relatively constant, variance is proportional to the mean and 
#' no transformations are needed. Calculate CV based on model fitted values and residuals due to 
#' presence of covariate.

augment(its_rich_lm) %>%
  mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant"))) %>% 
  group_by(field_type) %>%
  summarise(
    mean_fitted = mean(.fitted),
    sd_resid    = sd(.resid),
    cv_resid    = sd_resid / mean_fitted
  ) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>% 
  kable(format = "pandoc", caption = "CV of residuals and fitted means in groups")
#' CV is relatively constant

summary(its_rich_lm)
#' Sequence depth is a significant predictor of richness, but less so than field type. Does
#' depth confound interpretation of OTU richness in field types?
its_div %>% 
    group_by(field_type) %>% 
    summarize(across(c(depth, richness), ~ round(mean(.x), 0))) %>% 
    kable(format = "pandoc")
#' Sequence depth is similar in corn and restored fields, but richness in restored fields 
#' is approx 28% higher in restored fields. Sequence depth was lowest in remnant fields, 
#' where richness was the highest. This suggests that sequence depth isn't confounding the 
#' signal of field type, but it was appropriate to treat it as a covariate in the model
#' because it has it's own relationship to field type. 
#' 
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
#+ its_richness_fig,fig.width=4,fig.height=4,fig.align='center'
its_rich_fig <- 
  data.frame(summary(its_rich_em)) %>% 
  mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant"))) %>% 
    ggplot(aes(x = field_type, y = emmean)) +
    geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
    geom_errorbar(aes(ymin = emmean, ymax = upper.CL), width = 0, linewidth = lw) +
    geom_text(aes(y = upper.CL, label = c("a", "b", "b")), vjust = -1.5, family = "serif", size = 4) +
    labs(x = NULL, y = "Richness") +
    lims(y = c(0, 760)) +
    scale_fill_manual(values = c("gray", "black", "white")) +
    theme_cor +
    theme(legend.position = "none",
          plot.tag = element_text(size = 14, face = 1),
          plot.tag.position = c(2, 0.98))



#' ### Shannon's diversity


#' Account for sequencing depth as a covariate
its_shan_lm <- lm(shannon ~ sqrt(depth) + field_type, data = its_div)
par(mfrow = c(2,2))
plot(its_shan_lm) # variance similar in groups 
shapiro.test(residuals(its_shan_lm)) # p=0.81
performance::check_distribution(its_shan_lm) # residuals distribution most likely cauchy/normal; symmetric but long tails
#' Based on visual inspection and normality of residuals, log transformation or link function may not be warranted. 
#' Examine the CV in groups; if relatively constant, variance is proportional to the mean and 
#' no transformations are needed. Calculate CV based on model fitted values and residuals due to 
#' presence of covariate.

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
#' CV is not related to the mean

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
#+ its_richness_fig,fig.width=4,fig.height=4
its_shan_fig <- 
  data.frame(summary(its_shan_em)) %>% 
  mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant"))) %>% 
    ggplot(aes(x = field_type, y = emmean)) +
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




#' ## PLFA

plfa_lm <- lm(fungi_18.2 ~ field_type, data = fa)
par(mfrow = c(2,2))
plot(plfa_lm) # variance differs slightly in groups. Tails on qq plot diverge, lots of groups structure
shapiro.test(residuals(plfa_lm)) # passes, but close to alpha at 0.10
performance::check_distribution(plfa_lm) # response distribution fits normal best
#' Based on visual inspection and normality of residuals, log transformation or link function may not be warranted. 
#' Examine the CV in groups to look for any other evidence of non-constant variance.

fa %>%
  mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant"))) %>% 
  group_by(field_type) %>%
  summarize(mean = mean(fungi_18.2),
            cv = sd(fungi_18.2) / mean) %>% 
  mutate(across(mean:cv, ~ round(.x, 2))) %>% 
  kable(format = "pandoc", caption = "Mean and CV relationship in groups")
#' CV declines with the mean

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
  data.frame(summary(plfa_em)) %>% 
  mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant"))) %>% 
    ggplot(aes(x = field_type, y = emmean)) +
    geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
    geom_errorbar(aes(ymin = emmean, ymax = upper.CL), width = 0, linewidth = lw) +
    labs(x = NULL, y = expression(PLFA~(nmol%*%g[soil]^-1))) +
    scale_fill_manual(values = c("gray", "black", "white")) +
    theme_cor +
    theme(legend.position = "none",
          plot.tag = element_text(size = 14, face = 1),
          plot.tag.position = c(0, 1.02))




#' ## Beta Diversity

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





#' ## Unified figure
#+ fig2_patchwork,warning=FALSE
fig2_ls <- (its_rich_fig / plot_spacer() / plfa_fig) +
    plot_layout(heights = c(1,0.01,1)) 
fig2 <- (fig2_ls | plot_spacer() | its_ord) +
    plot_layout(widths = c(0.35, 0.01, 0.64)) +
    plot_annotation(tag_levels = 'a') 
#+ fig2,warning=FALSE,fig.height=4,fig.width=6.5
fig2
#' **Fig 2** Whole soil fungal communities from cornfields, restored, or remnant prairies,
#' with column charts showing **a** OTU richness and **b** fungal
#' biomass (PLFA). Error bars show 95% confidence intervals and lowercase letters show 
#' significant pairwise contrasts (P < 0.001). An ordination of community data **c** shows
#' a  principal coordinate analysis of soil fungal composition as inferred by clustering 
#' ITS sequences into 97% similar OTUs. Small circles depict individual sites and large circles show
#' centroids of clusters based on field type. Horizontal and vertical error bars around centroids encompass 95% confidence intervals around the 
#' mean location of sites in the cluster. In pairwise contrasts, cornfields clustered separately from 
#' restored or remnant prairies (P < 0.01). 
#' Text within the black circles indicates the number of years between restoration and 
#' collection of field samples. Percentages included in the axis titles indicate the percent of community variation explained on each axis 
#' out of the entire ordination. Across all plots, shading represents field type, with corn shaded gray, restored shaded black, 
#' and remnant shaded white. 

#+ fig2_save,warning=FALSE,fig.height=5,fig.width=7,echo=FALSE
ggsave(root_path("figs", "fig2.png"),
       plot = fig2,
       width = 6.5,
       height = 4,
       units = "in",
       dpi = 600)







#' # Arbuscular mycorrhizal fungi
# AMF ———————— ####
#' ## Diversity Indices
#+ amf_diversity
amf_div <- calc_div(spe$amf_avg)
#' 
#' ### Richness
#' Account for sequencing depth as a covariate
amf_rich_lm <- lm(richness ~ sqrt(depth) + field_type, data = amf_div)
par(mfrow = c(2,2))
plot(amf_rich_lm) # variance similar in groups with an outlier
shapiro.test(residuals(amf_rich_lm)) # p=0.20
performance::check_distribution(amf_rich_lm) # residuals distribution most likely normal. 
#' Based on visual inspection and normality of residuals, log transformation or link function may not be warranted. 
#' Examine the CV in groups; if relatively constant, variance is proportional to the mean and 
#' no transformations are needed. Calculate CV based on model fitted values and residuals due to 
#' presence of covariate.

augment(amf_rich_lm) %>%
  mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant"))) %>% 
  group_by(field_type) %>%
  summarise(
    mean_fitted = mean(.fitted),
    sd_resid    = sd(.resid),
    cv_resid    = sd_resid / mean_fitted
  ) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>% 
  kable(format = "pandoc", caption = "CV of residuals and fitted means in groups")
#' CV is most influenced by n in groups but shows no dependence on the mean

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
  data.frame(summary(amf_rich_em)) %>% 
  mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant"))) %>% 
    ggplot(aes(x = field_type, y = emmean)) +
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
#' ### Shannon's diversity

amf_shan_lm <- lm(shannon ~ sqrt(depth) + field_type, data = amf_div)
par(mfrow = c(2,2))
plot(amf_shan_lm) # variance somewhat non-constant in groups, qqplot fit is poor, one leverage point (Cook's >0.5)
shapiro.test(residuals(amf_shan_lm)) # p=0.38
performance::check_distribution(amf_shan_lm) # residuals distribution most likely normal. 
#' Based on visual inspection and normality of residuals, log transformation or link function may not be warranted. 
#' Examine the CV in groups; if relatively constant, variance is proportional to the mean and 
#' no transformations are needed. 

augment(amf_shan_lm) %>%
  mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant"))) %>% 
  group_by(field_type) %>%
  summarise(
    mean_fitted = mean(.fitted),
    sd_resid    = sd(.resid),
    cv_resid    = sd_resid / mean_fitted
  ) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>% 
  kable(format = "pandoc", caption = "CV of residuals and fitted means in groups")
#' CV is very stable

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
  data.frame(amf_shan_em) %>% 
  mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant"))) %>% 
    ggplot(aes(x = field_type, y = emmean)) +
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
#+ its_amf_shan_fig,warning=FALSE,fig.height=4,fig.width=7
FigS1 <- (its_shan_fig | plot_spacer() | amf_shan_fig) +
    plot_layout(widths = c(1, 0.01, 1)) +
    plot_annotation(tag_levels = 'a')

#+ amf_div_fig_save
ggsave(
    root_path("figs", "figS1.png"),
    plot = FigS1,
    height = 3,
    width = 6,
    units = "in",
    dpi = 600
)


#' ## NLFA

nlfa_lm <- lm(amf ~ field_type, data = fa)
par(mfrow = c(2,2))
plot(nlfa_lm) # variance obviously not constant in groups
shapiro.test(residuals(nlfa_lm)) # passes, but close to alpha 0.08
performance::check_distribution(nlfa_lm) # response distribution non-normal; resids likely normal
#' Based on visual inspection and normality of residuals, log transformation or link function warranted
#' based on non-constant variance. 
#' Examine the CV in groups; if relatively constant, variance is proportional to the mean and 
#' log-transformation is possibly acceptable. If the cv increases strongly with the mean, choose
#' choose gamma glm to account for increasing variance 

fa %>%
  mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant"))) %>% 
  group_by(field_type) %>%
  summarize(mean = mean(amf),
            cv = sd(amf) / mean) %>% 
  mutate(across(mean:cv, ~ round(.x, 2))) %>% 
  kable(format = "pandoc", caption = "Mean and CV relationship in groups")
#' CV increases strongly with the mean

nlfa_lm_log   <- lm(log(amf) ~ field_type, data = fa)
plot(nlfa_lm_log) # qqplot ok, one high leverage point in remnants
ncvTest(nlfa_lm_log) # p=0.16, null of constant variance not rejected

nlfa_glm  <- glm(amf ~ field_type, family = Gamma(link = "log"), data = fa)
nlfa_glm_diag <- glm.diag(nlfa_glm)
glm.diag.plots(nlfa_glm, nlfa_glm_diag) # qqplot shows strong fit; no leverage >0.5
performance::check_overdispersion(nlfa_glm) # not detected
#' Gamma glm is the best choice

#' Model results, group means, and post-hoc
#' - Arithmetic means will be returned from this model by emmeans()
summary(nlfa_glm)
anova(nlfa_glm) # Decline in residual deviance worth the cost in df
nlfa_em <- emmeans(nlfa_glm, ~ field_type, type = "response")
#+ plfa_em_summary,echo=FALSE
kable(summary(nlfa_em), 
      format = "pandoc", 
      caption = "Confidence level used: 0.95.\nIntervals are back-transformed from the log scale")
#+ plfa_em_posthoc,echo=FALSE
kable(pairs(nlfa_em), 
      format = "pandoc", 
      caption = "P value adjustment: tukey method for comparing a family of 3 estimates.\nTests are performed on the log scale")



#+ plfa_fig,fig.width=4,fig.height=4,
nlfa_fig <-
  data.frame(summary(nlfa_em)) %>% mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant"))) %>% 
  ggplot(aes(x = field_type, y = response)) +
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




#' ## Beta Diversity

#' ## AMF (18S sequences)
#' UNIFRAC distance matrix must be created. 
#+ amf_ord
d_amf <- UniFrac(amf_avg_ps, weighted = TRUE)
mva_amf <- mva(d = d_amf, corr = "lingoes")
#+ amf_ord_results
mva_amf$dispersion_test
mva_amf$permanova
mva_amf$pairwise_contrasts
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
#' **Fig 3** AMF communities from cornfields, restored, or remnant prairies,
#' with column charts showing **a** OTU richness and **b** fungal
#' biomass (NLFA). Error bars show 95% confidence intervals and lowercase letters show 
#' significant pairwise contrasts (**a** P < 0.05, **b** P < 0.0001). 
#' Ordination of sites **c** from principal coordinate analysis of AM fungal composition was based on UNIFRAC distance
#' of 18S sequences clustered into 97% similar OTUs. Small circles depict locations of individual sites and large circles show
#' centroids of clusters based on field type. Shading represents field type, with corn shaded gray, restored shaded black, 
#' and remnant shaded white. Horizontal and vertical error bars around centroids encompass 95% confidence intervals around the 
#' mean location of sites in the cluster. Text within the black circles indicates the number of years between restoration and 
#' collection of field samples. Percentages included in the axis titles indicate the percent of community variation explained on each axis 
#' out of the entire ordination. In pairwise contrasts, cornfields clustered separately from 
#' restored or remnant prairies (P < 0.01). 
#' Text within the black circles indicates the number of years between restoration and 
#' collection of field samples. Percentages included in the axis titles indicate the percent of community variation explained on each axis 
#' out of the entire ordination. Across all plots, shading represents field type, with corn shaded gray, restored shaded black, 
#' and remnant shaded white. 

#+ fig3_save,warning=FALSE,fig.height=5,fig.width=7,echo=FALSE
ggsave(root_path("figs", "fig3.png"),
       plot = fig3,
       width = 6.5,
       height = 4,
       units = "in",
       dpi = 600)



