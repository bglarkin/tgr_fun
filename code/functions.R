#' ---
#' title: "Supplement: Functions"
#' author: "Beau Larkin\n"
#' date: "Last updated: `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   github_document:
#'     toc: true
#'     toc_depth: 2
#' ---
#' 
#' # Description
#' Functions that accompany the repository, sourced from this file to save space elsewhere
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

#' ## Permanova on soil data
#' Simplified version of `mva()` for use with the soil properties data
soilperm <- function(clust_vec, clust_var) {
  #' ### Permanova
  soil_d <- dist(site_sco)
  
  # Test multivariate dispersions
  soil_disper <- betadisper(soil_d, clust_vec, bias.adjust = TRUE)
  soil_mvdisper <- permutest(soil_disper, pairwise = TRUE, permutations = 1999)
  # Global PERMANOVA
  soil_gl_permtest <- adonis2(
    as.formula(paste("soil_d ~", clust_var)),
    data = soil_ord_scores,
    permutations = 1999,
    by = "terms")
  # Pairwise PERMANOVA
  group_var <- as.character(clust_vec)
  groups <- as.data.frame(t(combn(unique(group_var), m = 2)))
  soil_contrasts <- data.frame(
    group1 = groups$V1,
    group2 = groups$V2,
    R2 = NA,
    F_value = NA,
    df1 = NA,
    df2 = NA,
    p_value = NA
  )
  for (i in seq(nrow(soil_contrasts))) {
    group_subset <-
      group_var == soil_contrasts$group1[i] |
      group_var == soil_contrasts$group2[i]
    contrast_d <- as.matrix(soil_d)[group_subset, group_subset]
    fit <- adonis2(
      contrast_d ~ group_var[group_subset],
      permutations = 1999,
      by = "terms")
    # Prepare contrasts table
    soil_contrasts$R2[i] <- round(fit[grep("group_var", rownames(fit)), "R2"], digits = 3)
    soil_contrasts$F_value[i] <- round(fit[grep("group_var", rownames(fit)), "F"], digits = 3)
    soil_contrasts$df1[i] <- fit[grep("group_var", rownames(fit)), "Df"]
    soil_contrasts$df2[i] <- fit[grep("Residual", rownames(fit)), "Df"]
    soil_contrasts$p_value[i] <- fit[grep("group_var", rownames(fit)), 5]
  }
  soil_contrasts$p_value_adj <- p.adjust(soil_contrasts$p_value, method = "fdr") %>% round(., 4)
  
  out = list(
    mvdisper = soil_mvdisper,
    gl_permtest = soil_gl_permtest,
    contrasts = soil_contrasts
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

#' ## Perform Indicator Species Analysis
#' Function `inspan()` takes a combined species and sites data frame and 
#' filters OTUs for indicators of field types. 
#+ inspan_function
inspan <- function(guild, nperm=1999) {
  data <- guildseq(spe$its_avg, guild) %>% 
    left_join(sites, by = join_by(field_name))
  spe <- data.frame(
    data %>% select(field_name, starts_with("otu")),
    row.names = 1
  )
  grp = data$field_type
  mp <- multipatt(
    spe, grp, max.order = 1, 
    control = how(nperm = nperm))
  si <- mp$sign %>% 
    select(index, stat, p.value) %>% 
    mutate(field_type = case_when(index == 1 ~ "corn", 
                                  index == 2 ~ "restored", 
                                  index == 3 ~ "remnant"),
           p_val_adj = p.adjust(p.value, "fdr")) %>% 
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
    left_join(spe_meta$its %>% select(-otu_ID), by = join_by(otu_num)) %>% 
    select(otu_num, A, B, stat, p.value, p_val_adj, 
           field_type, primary_lifestyle, everything()) %>% 
    arrange(field_type, -stat)
  
  return(out)
  
}
