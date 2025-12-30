Supplement: Functions
================
Beau Larkin

Last updated: 30 December, 2025

- [Description](#description)
  - [Sequence data processing
    functions](#sequence-data-processing-functions)
  - [Confidence interval helper](#confidence-interval-helper)
  - [Alpha diversity calculations](#alpha-diversity-calculations)
  - [RMSE](#rmse)
  - [Confidence intervals](#confidence-intervals)
  - [Multivariate analysis](#multivariate-analysis)
  - [Permanova on soil data](#permanova-on-soil-data)
  - [Model distribution
    probabilities](#model-distribution-probabilities)
  - [Filter spe to a guild](#filter-spe-to-a-guild)
  - [Perform Indicator Species
    Analysis](#perform-indicator-species-analysis)
  - [Calculate pairwise distances among sites and present summary
    statistics](#calculate-pairwise-distances-among-sites-and-present-summary-statistics)
  - [Test raw and transformed covariates in linear
    models](#test-raw-and-transformed-covariates-in-linear-models)
  - [Site mapping functions](#site-mapping-functions)

# Description

Functions that accompany the repository, sourced from this file to save
space elsewhere

## Sequence data processing functions

### ETL: clean OTU data and return formatted objects

``` r
etl <- function(spe, env = sites, taxa, traits = NULL, varname, gene, cluster_type = "otu",
                colname_prefix, folder) {
  varname <- enquo(varname)
  data <- spe %>% left_join(taxa, by = join_by(`#OTU ID`))
  
  meta <- if (gene == "ITS") {
    data %>%
      mutate(!!varname := paste0(cluster_type, "_", row_number())) %>%
      select(!starts_with(gene)) %>%
      rename(otu_ID = `#OTU ID`) %>%
      select(!!varname, everything()) %>%
      separate_wider_delim(taxonomy, delim = ";",
                           names = c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
                           too_few = "align_start") %>%
      mutate(across(kingdom:species, ~ str_sub(.x, 4))) %>%
      left_join(traits, by = join_by(phylum, class, order, family, genus)) %>%
      select(-kingdom, -Confidence)
  } else {
    data %>%
      mutate(!!varname := paste0(cluster_type, "_", row_number())) %>%
      select(!starts_with(gene)) %>%
      rename(otu_ID = `#OTU ID`) %>%
      select(!!varname, everything()) %>%
      separate(taxonomy,
               into = c("class", "order", "family", "genus", "taxon", "accession"),
               sep = ";", fill = "right") %>%
      select(-Confidence)
  }
  
  spe_samps <- data %>%
    mutate(!!varname := paste0(cluster_type, "_", row_number())) %>%
    select(!!varname, starts_with(gene)) %>%
    column_to_rownames(var = as_name(varname)) %>%
    t() %>% as.data.frame() %>% rownames_to_column("rowname") %>%
    mutate(rowname = str_remove(rowname, colname_prefix)) %>%
    separate_wider_delim("rowname", delim = "_", names = c("field_key", "sample")) %>%
    mutate(across(c(field_key, sample), as.numeric)) %>%
    left_join(env %>% select(field_key, field_name), by = join_by(field_key)) %>%
    select(field_name, sample, everything(), -field_key) %>%
    arrange(field_name, sample)
  
  spe_avg <- spe_samps %>%
    group_by(field_name) %>%
    summarize(across(starts_with(cluster_type), mean), .groups = "drop") %>%
    arrange(field_name)
  
  samples_fields <- spe_samps %>%
    count(field_name, name = "n") %>%
    left_join(sites, by = join_by(field_name)) %>%
    select(field_name, region, n) %>%
    arrange(n, field_name) %>%
    kable(format = "pandoc", caption = paste("Number of samples per field", gene, sep = ",\n"))
  
  write_csv(meta, root_path(folder, paste("spe", gene, "metadata.csv", sep = "_")))
  write_csv(spe_samps, root_path(folder, paste("spe", gene, "samples.csv", sep = "_")))
  write_csv(spe_avg, root_path(folder, paste("spe", gene, "avg.csv", sep = "_")))
  
  list(samples_fields = samples_fields, spe_meta = meta, spe_samps = spe_samps, spe_avg = spe_avg)
}
```

### Species accumulation

``` r
spe_accum <- function(data) {
  df <- data.frame(
    samples = specaccum(data[, -c(1, 2)], method = "exact")$site,
    richness = specaccum(data[, -c(1, 2)], method = "exact")$richness,
    sd = specaccum(data[, -c(1, 2)], method = "exact")$sd
  )
  df
}
```

## Confidence interval helper

``` r
ci <- function(x) std.error(x) * qnorm(0.975)
```

## Alpha diversity calculations

Returns a dataframe of alpha diversity (richness, Shannon’s) for
analysis and plotting.

``` r
calc_div <- function(spe, site_dat) {
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
    ungroup() %>% 
    left_join(site_dat %>% select(field_type, field_name), by = join_by(field_name)) %>% 
    select(field_name, field_type, depth, richness, shannon)
  
  return(div_data)
  
}
```

## RMSE

Compare model results

``` r
rmse <- function(obs, pred) {
  sqrt(mean((obs - pred)^2, na.rm = TRUE))
}
```

## Confidence intervals

Calculate upper and lower confidence intervals with alpha=0.05

``` r
ci_u <- function(x) {(sd(x) / sqrt(length(x))) * qnorm(0.975)}
ci_l <- function(x) {(sd(x) / sqrt(length(x))) * qnorm(0.025)}
```

## Multivariate analysis

Ordination → dispersion check → global & pairwise PERMANOVA → envfit.
Args: *d* dist, *env* metadata, *corr* PCoA correction, *nperm*
permutations.

``` r
mva <- function(d, env, corr = "none", nperm = 1999, seed = 20251101) {
  stopifnot(is.data.frame(env))
  if (!all(c("field_type", "dist_axis_1") %in% names(env))) {
    stop("`env` must contain columns `field_type` and `dist_axis_1`.")
  }
  
  # Distance labels and env alignment
  if (inherits(d, "dist")) {
    lab <- attr(d, "Labels")
  } else if (is.matrix(d)) {
    if (is.null(rownames(d))) stop("Distance matrix `d` must have row names.")
    lab <- rownames(d)
    d <- as.dist(d)  # coerce for betadisper/pcoa convenience
  } else {
    stop("`d` must be a 'dist' or a symmetric distance matrix.")
  }
  
  # Align feature names across data sources
  env <- as.data.frame(env)
  if ("field_name" %in% names(env)) rownames(env) <- env$field_name
  if (!setequal(rownames(env), lab)) {
    miss_env <- setdiff(lab, rownames(env))
    miss_d   <- setdiff(rownames(env), lab)
    stop("Sample mismatch between `d` and `env`.\n",
         "In d not in env: ", paste(miss_env, collapse = ", "),
         "\nIn env not in d: ", paste(miss_d, collapse = ", "))
  }
  env <- env[lab, , drop = FALSE]  
  
  # Set grouping variables
  g_chr     <- as.character(env$field_type)
  g_levels  <- sort(unique(g_chr))         # deterministic order
  clust_vec <- factor(g_chr, levels = g_levels)
  
  # Ordination (PCoA)
  p <- pcoa(d, correction = corr)
  p_vals <- data.frame(p$values) %>% 
    rownames_to_column(var = "Dim") %>% 
    mutate(Dim = as.integer(Dim))
  p_eig <- p_vals[1:2, grep("Rel", colnames(p_vals))] %>% round(., 3) * 100
  
  p_vec <- data.frame(p$vectors, check.names = FALSE)
  p_sco <- p_vec[, 1:2, drop = FALSE] %>% 
    rownames_to_column(var = "field_name") %>% 
    left_join(env, by = join_by(field_name))
  
  # Homogeneity of multivariate dispersion
  disper <- betadisper(d, clust_vec, bias.adjust = TRUE)
  if (!is.null(seed)) set.seed(seed + 2L)
  mvdisper <- permutest(disper, pairwise = TRUE, permutations = nperm)
  
  # Global PERMANOVA
  if (!is.null(seed)) set.seed(seed + 3L)
  gl_permtest <- adonis2(
    d ~ dist_axis_1 + field_type,
    data = env,
    permutations = nperm,
    by = "terms"
  )
  
  # Pairwise PERMANOVA
  groups <- combn(g_levels, m = 2) %>%  t() %>%  as.data.frame()
  names(groups) <- c("V1", "V2")
  
  contrasts <- data.frame(
    group1   = groups$V1,
    group2   = groups$V2,
    R2       = NA_real_,
    F_value  = NA_real_,
    df1      = NA_integer_,
    df2      = NA_integer_,
    p_value  = NA_real_
  )
  
  d_mat <- as.matrix(d)
  for (i in seq_len(nrow(contrasts))) {
    g1 <- contrasts$group1[i]
    g2 <- contrasts$group2[i]
    keep <- clust_vec %in% c(g1, g2)
    
    contrast_mat <- d_mat[keep, keep, drop = FALSE]
    
    if (!is.null(seed)) set.seed(seed + 100L + i)
    
    fit <- adonis2(
      contrast_mat ~ env$dist_axis_1[keep] + factor(g_chr[keep]),
      permutations = nperm,
      by = "terms"
    )
    
    rn <- rownames(fit)
    term_row <- grep("factor\\(g_chr\\[keep\\]\\)", rn, fixed = FALSE)
    if (length(term_row) != 1L) term_row <- grep("g_chr", rn, fixed = TRUE)
    
    contrasts$R2[i]      <- round(fit[term_row, "R2"], 3)
    contrasts$F_value[i] <- round(fit[term_row, "F"], 3)
    contrasts$df1[i]     <- fit[term_row, "Df"]
    contrasts$df2[i]     <- fit[grep("Residual", rn), "Df"]
    contrasts$p_value[i] <- fit[term_row, 5]
  }
  contrasts$p_value_adj <- p.adjust(contrasts$p_value, method = "fdr") %>% round(4)
  
  # Results
  list(
    correction_note    = p$note,
    ordination_values  = p_vals[1:min(10, nrow(p_vals)), ],
    axis_pct           = p_eig,
    ordination_scores  = p_sco,
    dispersion_test    = mvdisper,
    permanova          = gl_permtest,
    pairwise_contrasts = contrasts
  )
}
```

## Permanova on soil data

Simplified version of `mva()` for use with the soil properties data

``` r
soilperm <- function(ord_scores, clust_var, nperm = 1999, seed = 20251103) {
  # Distance object
  soil_d <- ord_scores %>%
    select(field_name, PC1, PC2) %>%
    column_to_rownames("field_name") %>%
    dist()
  
  # Clusters
  clust_vec <- ord_scores %>% pull({{clust_var}})
  g_chr <- as.character(clust_vec)
  
  # Control and reproducibility
  set.seed(seed)
  ctrl <- how(nperm = nperm)  
  
  # Homogeneity of multivariate dispersion
  soil_disper <- betadisper(soil_d, clust_vec, bias.adjust = TRUE)
  soil_mvdisper <- permutest(soil_disper, pairwise = TRUE, permutations = ctrl)
  
  # Global PERMANOVA
  perm_glob <- shuffleSet(n = attr(soil_d, "Size"), control = ctrl)
  soil_gl_permtest <- adonis2(soil_d ~ clust_vec,
                              permutations = perm_glob,
                              by = "terms",
                              parallel = 1)  # force serial for determinism
  
  # Pairwise PERMANOVA
  lev <- sort(unique(g_chr))
  pairs <- t(combn(lev, 2))
  res <- vector("list", nrow(pairs))
  
  for (i in seq_len(nrow(pairs))) {
    keep <- g_chr %in% pairs[i, ]
    d_sub <- as.dist(as.matrix(soil_d)[keep, keep])
    g_sub <- g_chr[keep]
    
    perm_pw <- shuffleSet(n = attr(d_sub, "Size"), control = ctrl)
    
    fit <- adonis2(d_sub ~ g_sub,
                   permutations = perm_pw,
                   by = "terms",
                   parallel = 1)
    
    res[[i]] <- data.frame(
      group1  = pairs[i, 1],
      group2  = pairs[i, 2],
      R2      = unname(fit$R2[1]),
      F_value = unname(fit$F[1]),
      df1     = unname(fit$Df[1]),
      df2     = unname(fit$Df[nrow(fit)]),
      p_value = unname(fit$`Pr(>F)`[1])
    )
  }
  
  soil_contrasts <- bind_rows(res) %>%
    mutate(p_value_adj = p.adjust(p_value, method = "fdr"))
  
  # Results
  list(
    mvdisper    = soil_mvdisper,
    gl_permtest = soil_gl_permtest,
    contrasts   = soil_contrasts
  )
}
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

## Filter spe to a guild

Create samp-spe matrix of sequence abundance in a guild

``` r
guildseq <- function(spe, meta, guild) {
  guab <- 
    spe %>% 
    pivot_longer(starts_with("otu"), names_to = "otu_num", values_to = "abund") %>% 
    left_join(meta %>% select(otu_num, primary_lifestyle), by = join_by(otu_num)) %>% 
    filter(primary_lifestyle == guild) %>% 
    select(-primary_lifestyle) %>% 
    pivot_wider(names_from = otu_num, values_from = abund)
  return(guab)
}
```

## Perform Indicator Species Analysis

Function `inspan()` takes a combined species and sites data frame and
filters OTUs for indicators of field types.

``` r
inspan <- function(spe, meta, guild, site_dat, nperm=1999) {
  if(is.null(guild)) {
    data <- spe %>% left_join(site_dat, by = join_by(field_name))
    spe_g <- data.frame(
      data %>% select(field_name, starts_with("otu")),
      row.names = 1)
    grp <- data$field_type
  } else {
    data <- guildseq(spe, meta, guild) %>% 
      left_join(site_dat, by = join_by(field_name))
    spe_g <- data.frame(
      data %>% select(field_name, starts_with("otu")),
      row.names = 1)
    grp <- data$field_type
  }
  
  # Indicator species analysis
  mp <- multipatt(
    x = spe_g, cluster = grp, func = "IndVal.g", max.order = 1, 
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
  
  # Join to abundance in field types
  seq_abund <- 
    spe_g %>% 
    rownames_to_column(var = "field_name") %>% 
    pivot_longer(starts_with("otu"), names_to = "otu_num", values_to = "abund") %>% 
    left_join(site_dat %>% select(field_name, field_type), by = join_by(field_name)) %>% 
    group_by(field_type, otu_num) %>% 
    summarize(avg = mean(abund),
              ci = qnorm(0.975) * sd(abund) / sqrt(n()),
              .groups = "drop") %>% 
    pivot_wider(names_from = "field_type", values_from = c("avg", "ci"), names_glue = "{field_type}_{.value}") %>% 
    select(otu_num, starts_with("corn"), starts_with("restor"), starts_with("rem"))
  
  out <- 
    si %>% 
    left_join(A, by = join_by(otu_num, field_type)) %>% 
    left_join(B, by = join_by(otu_num, field_type)) %>% 
    left_join(meta %>% select(-otu_ID), by = join_by(otu_num)) %>% 
    left_join(seq_abund, by = join_by(otu_num)) %>% 
    select(otu_num, A, B, stat, p.value, p_val_adj, 
           field_type, everything()) %>% 
    arrange(field_type, -stat)
  
  return(out)
  
}
```

## Calculate pairwise distances among sites and present summary statistics

Function `reg_dist_stats()` requires a haversine distance matrix, site
metadata, and is filtered by regions to produce the desired output.

``` r
reg_dist_stats <- function(dist_mat,
                           sites_df,
                           filt_rg) {
  dm  <- round(as.matrix(dist_mat) / 1000, 1)
  idx <- which(upper.tri(dm), arr.ind = TRUE)
  
  pairs <- tibble(site1 = rownames(dm)[idx[, 1]],
                  site2 = colnames(dm)[idx[, 2]],
                  dist  = dm[idx])
  
  meta <- sites_df %>%
    select(site = field_key, ft = field_type, rg = region) %>% 
    mutate(site = as.character(site), ft = as.character(ft)) %>% 
    filter(rg == filt_rg)
  
  pairs %>%
    filter(site1 %in% meta$site, site2 %in% meta$site) %>% 
    left_join(meta %>% select(site, ft), by = c("site1" = "site")) %>% rename(ft1 = ft) %>%
    left_join(meta %>% select(site, ft), by = c("site2" = "site")) %>% rename(ft2 = ft) %>%
    mutate(group_pair = paste(pmin(ft1, ft2), pmax(ft1, ft2), sep = "-")) %>%
    group_by(group_pair) %>%
    summarize(
      min_dist = min(dist, na.rm = TRUE),
      median_dist = median(dist, na.rm = TRUE),
      max_dist = max(dist, na.rm = TRUE),
      .groups = "drop"
    )
}
```

## Test raw and transformed covariates in linear models

Function `covar_shape_test()` compares a series of linear models with
raw, square root, or log transforms of the covariate. It selects the
best model based on various criteria.

``` r
covar_shape_test <- function(data, y, covar, group = field_type) {
  
  resp <- ensym(y)
  x    <- ensym(covar)
  g    <- ensym(group)
  
  df <- data %>% 
    select(!!resp, !!g, x_raw = !!x) %>% 
    mutate(x_lin  = as.numeric(scale(x_raw, center = TRUE, scale = FALSE)),
           x_sqrt = as.numeric(scale(sqrt(pmax(x_raw, 0)), center = TRUE, scale = FALSE)),
           x_log  = as.numeric(scale(log1p(pmax(x_raw, 0)), center = TRUE, scale = FALSE))
    ) %>% drop_na()
  
  f_lin  <- new_formula(expr(!!resp), expr(x_lin  + !!g))
  f_sqrt <- new_formula(expr(!!resp), expr(x_sqrt + !!g))
  f_log  <- new_formula(expr(!!resp), expr(x_log  + !!g))
  
  # Fit candidate models
  m_lin  <- lm(f_lin,  data = df)
  m_sqrt <- lm(f_sqrt, data = df)
  m_log  <- lm(f_log,  data = df)
  
  # Compare candidate fits
  cmp <- compare_performance(
    lin = m_lin, sqrt = m_sqrt, log = m_log,
    metrics = c("AIC","RMSE","R2"), rank = TRUE
  )
  best_name <- cmp$Name[1]
  best_mod  <- switch(best_name, m_lin = m_lin, m_sqrt = m_sqrt, m_log = m_log)
  
  # Type-II tests (unbalanced design; additive model)
  anova_t2 <- Anova(best_mod, type = 2)
  
  # LS-means for group at centered covariate
  emm <- emmeans(best_mod, specs = as_name(g))
  
  chmd <- check_model(best_mod)
  
  # return everything
  list(
    data      = df,
    compare   = cmp,
    best_name = best_name,
    best_model= best_mod,
    anova_t2  = anova_t2,
    emmeans   = emm,
    diagnostics = chmd
  )
}
```

## Site mapping functions

Functions facilitate the creation of mapping objects

### Create buffered bounding boxes

``` r
bbox_buffer_km <- function(pts_sf, buffer_km = 20) {
  albers <- 5070 # NAD83 / Conus Albers
  pts_sf %>%
    st_transform(albers) %>%
    st_bbox() %>%
    st_as_sfc(crs = albers) %>%
    st_buffer(dist = buffer_km * 1000) %>%
    st_transform(4326) %>%
    st_bbox()
}
```

### Create custom point nudges

Add field names manually: - nudge_e_m = meters to move EAST (positive =
east, negative = west) - nudge_n_m = meters to move NORTH (positive =
north, negative = south)

``` r
nudge_coords <- function(sites, lon_col = "long", lat_col = "lat",
                         east_col = "nudge_e_m", north_col = "nudge_n_m") {
  stopifnot(all(c(lon_col, lat_col) %in% names(sites)))
  # If nudge columns are missing, treat as zeros
  if (!(east_col  %in% names(sites))) sites[[east_col]]  <- 0
  if (!(north_col %in% names(sites))) sites[[north_col]] <- 0
  
  # Vectorized meters → degrees (WGS84 approximation)
  m_per_deg_lat <- 111320                      # ~ meters per degree latitude
  m_per_deg_lon <- m_per_deg_lat * cospi(sites[[lat_col]] / 180)
  
  dx_deg <- sites[[east_col]]  / m_per_deg_lon
  dy_deg <- sites[[north_col]] / m_per_deg_lat
  
  sites %>%
    mutate(
      long_plot = .data[[lon_col]] + ifelse(is.finite(dx_deg), dx_deg, 0),
      lat_plot  = .data[[lat_col]] + ifelse(is.finite(dy_deg), dy_deg, 0)
    )
}
```

### Retrieve road layer data from open street maps.

``` r
get_osm_roads <- function(bb, density = 8) {
  
  # validate density
  types <- c("motorway","trunk","primary","secondary",
             "tertiary","unclassified","residential","service")
  if (!is.numeric(density) || length(density) != 1L || !is.finite(density)) {
    stop("`density` must be a single finite number in 1:8")
  }
  density <- as.integer(max(1L, min(length(types), density)))
  vals <- types[seq_len(density)]
  
  # build layer
  res <- 
    opq(bbox = bb) %>% 
    add_osm_feature(key = "highway", value = vals) %>% 
    osmdata_sf()
  
  return(st_crop(res$osm_lines, st_as_sfc(bb)))
  
}
```

### Build site-level map panels

``` r
make_zoom_map <- function(bb, panel_tag = NULL, pos = c(0,1), show_counties = FALSE, road_data = NULL) {
  
  crop_states   <- st_crop(cont, bb)
  crop_counties <- if (show_counties) st_crop(st_transform(counties, 4326), bb) else NULL
  roads         <- road_data
  
  pts <- sites_sf %>%
    filter(long >= bb["xmin"], long <= bb["xmax"],
           lat  >= bb["ymin"], lat  <= bb["ymax"]) %>%
    st_drop_geometry()
  
  # Labels only where yr_since is available (restored sites)
  pts_lab <- sites_plot %>%
    filter(!is.na(yr_since)) %>%
    mutate(lbl = as.character(round(yr_since, 0)))
  
  g <- ggplot() +
    geom_sf(data = crop_states, fill = "ivory", color = "black", linewidth = 0.5) +
    { if (!is.null(crop_counties)) geom_sf(data = crop_counties, fill = NA, color = "gray85", linewidth = 0.3) } +
    { if (!is.null(roads)) geom_sf(data = roads, color = "grey70", linewidth = 0.3) } +
    geom_point(
      data = sites_plot,
      aes(x = long_plot, y = lat_plot, fill = field_type),
      shape = 21, size = sm_size, stroke = lw, color = "black"
    ) +
    geom_text(
      data = pts_lab,
      aes(x = long_plot, y = lat_plot, label = lbl),
      size = yrtx_size, family = "sans", fontface = 2, color = "black"
    ) +
    scale_fill_manual(values = ft_pal) +
    annotation_scale(location = "bl", width_hint = 0.35, height = grid::unit(0.15, "cm")) +
    # geom_label_npc(
    #   aes(npcx = 0.02, npcy = 0.98, label = panel_tag %||% ""),
    #   color = "indianred",
    #   fill = "snow",
    #   hjust = "left", vjust = "top",
    #   size = 3, fontface = "bold",
    #   label.r = unit(0.3, "mm"),
    #   label.size = 0.4,
    #   label.padding = unit(c(0.4, 0.3, 0.15, 0.3), "lines")
    # ) +
    coord_sf(
      xlim = c(bb["xmin"], bb["xmax"]),
      ylim = c(bb["ymin"], bb["ymax"]),
      expand = FALSE
    ) +
    labs(tag = panel_tag) +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "aliceblue", color = "black", linewidth = 0.5),
      legend.position = "none",
      plot.tag = element_text(size = 14, face = 1, hjust = 0),
      plot.tag.position = pos
    )
  
  g
  
}
```
