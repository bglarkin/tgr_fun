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
#' Inter‑site distance enters models as a covariate per [Redondo 2020](https://doi.org/10.1093/femsec/fiaa082); 
#' Moran's eigenvalues tested and used where significant.
#' 
#' # Packages and libraries
# Libraries ———————— ####
#+ packages,message=FALSE
packages_needed <- c(
  # Analysis
  "emmeans", "vegan", "phyloseq", "ape", "phangorn", "geosphere", 
  "car", "rlang", "rsq", "sandwich", "lmtest", "performance", "boot",
  "MASS", "DHARMa", "broom", "ALDEx2", "adespatial",
  # Scripting
  "rprojroot", "conflicted", "purrr", "knitr", "tidyverse", 
  # Graphics
  "colorspace", "grid", "gridExtra", "ggpubr", "patchwork", "ggpattern"
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
#' # Functions
#' Executed from a separate script to save lines here; to view the function navigate to 
#' `functions.R` in the code folder, accessible from the root dir of the repo.
# Functions ———————— ####
source(root_path("code", "functions.R"))

#' 
#' # Data
# Data ———————— ####
#' Loading order reflects downstream dependencies
#' 
#' ## Site metadata
sites <- read_csv(root_path("clean_data/sites.csv"), show_col_types = FALSE) %>% 
  mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant")))
sites_wi <- sites %>% 
  filter(region != "FL", field_type != "corn")
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
its_avg = read_csv(root_path("clean_data/spe_ITS_avg.csv"), show_col_types = FALSE)
its_wi  = its_avg %>% 
  filter(field_name %in% sites_wi$field_name) %>% 
  select(field_name, where(~ is.numeric(.x) && sum(.x) > 0))
amf_avg = read_csv(root_path("clean_data/spe_18S_avg.csv"), show_col_types = FALSE)
amf_wi  = amf_avg %>% 
  filter(field_name %in% sites_wi$field_name) %>% 
  select(field_name, where(~ is.numeric(.x) && sum(.x) > 0))
#' 
#' ### Microbial species metadata
its_meta = read_csv(root_path("clean_data/spe_ITS_metadata.csv"), show_col_types = FALSE) %>% 
  mutate(primary_lifestyle = case_when(str_detect(primary_lifestyle, "_saprotroph$") ~ "saprotroph",
                                       str_detect(primary_lifestyle, "unspecified_path") ~ "unidentified",
                                       TRUE ~ primary_lifestyle),
         across(everything(), ~ replace_na(., "unidentified")))
amf_meta = read_csv(root_path("clean_data/spe_18S_metadata.csv"), show_col_types = FALSE) %>% 
  mutate(across(everything(), ~ replace_na(., "unidentified")))
#' 
#' ### Spe subsets for guilds, regions
patho <- guildseq(its_avg, its_meta, "plant_pathogen")
patho_wi <- guildseq(its_wi, its_meta, "plant_pathogen")
sapro <- guildseq(its_avg, its_meta, "saprotroph")
sapro_wi <- guildseq(its_wi, its_meta, "saprotroph")
#' 
#' ### Additional species and metadata objects
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
#+ amf_ps_object,message=FALSE,warning=FALSE
amf_ps <- phyloseq(
  otu_table(data.frame(amf_avg_uni, row.names = 1), taxa_are_rows = TRUE),
  tax_table(as.matrix(data.frame(amf_meta, row.names = 2))),
  read.dna(root_path("otu_tables/18S/18S_sequences.fasta"), format = "fasta") %>%
    phyDat(type = "DNA") %>% dist.hamming() %>% NJ(),
  sample_data(sites %>% column_to_rownames(var = "field_name"))
)
#' 
#' ### Species distance matrices
#' #### All sites
d_all <- list(
  d_its = its_avg,
  d_patho = patho,
  d_sapro = sapro
) %>% map(\(df) df %>% 
            data.frame(row.names = 1) %>%
            decostand("total") %>%
            vegdist("bray"))
d_all$d_amf <- UniFrac(amf_ps, weighted = TRUE, normalized = TRUE)
d_all$d_amf_ma <- amf_avg_ma %>% 
  data.frame(row.names = 1) %>% 
  vegdist("bray")
#' 
#' #### Wisconsin sites
d_wi <- 
  list(
    d_its_wi   = its_wi,
    d_patho_wi = patho_wi,
    d_sapro_wi = sapro_wi
  ) %>% 
  map(\(df) df %>% 
        data.frame(row.names = 1) %>% 
        decostand("total") %>% 
        vegdist("bray"))
#+ amf_uni_prune,warning=FALSE,message=FALSE
amf_ps_wi <- prune_samples(
  sites %>% filter(region != "FL", field_type != "corn") %>% pull(field_name), 
  amf_ps
) %>% prune_taxa(taxa_sums(.) > 0, .)
d_wi$d_amf_wi <- UniFrac(amf_ps_wi, weighted = TRUE, normalized = TRUE)
#' 
#' ## Inter-site distance
#' Inter-site geographic distance will be considered as a covariate in clustering and 
#' regression analyses. Compute and use Moran's eigenvalues separately for all sites
#' and restored+remnant in Wisconsin. 
#' 
#' ### All sites
#' db-MEM
coord_tbl <- sites %>% select(long, lat) %>% as.matrix()
rownames(coord_tbl) <- sites$field_name
mem <- dbmem(coord_tbl) %>% as.data.frame()
setequal(rownames(d_all$d_its), rownames(mem))
setequal(rownames(d_all$d_amf), rownames(mem))
setequal(rownames(d_all$d_patho), rownames(mem))
setequal(rownames(d_all$d_sapro), rownames(mem))
#' 
#' #### ITS fungi
mem_null_its <- dbrda(d_all$d_its ~ 1, data = mem)
mem_full_its <- dbrda(d_all$d_its ~ ., data = mem)
mem_step_its <- ordistep(mem_null_its, scope = formula(mem_full_its), direction = "forward", 
                         permutations = 1999, trace = FALSE)
RsquareAdj(mem_step_its, permutations = 1999)$adj.r.squared
anova(mem_step_its, by = "margin", permutations = 1999) %>% 
  as.data.frame() %>% 
  mutate(p.adj = p.adjust(`Pr(>F)`, "fdr")) %>% 
  kable(, format = "pandoc")
#' MEM2
#' 
#' #### AMF
#' Unifrac distance
mem_null_amf <- dbrda(d_all$d_amf ~ 1, data = mem)
mem_full_amf <- dbrda(d_all$d_amf ~ ., data = mem)
mem_step_amf <- ordistep(mem_null_amf, scope = formula(mem_full_amf), direction = "forward", 
                         permutations = 1999, trace = FALSE)
RsquareAdj(mem_step_amf, permutations = 1999)$adj.r.squared
anova(mem_step_amf, by = "margin", permutations = 1999) %>% 
  as.data.frame() %>% 
  mutate(p.adj = p.adjust(`Pr(>F)`, "fdr")) %>% 
  kable(, format = "pandoc")
#' None
#' 
#' #### Pathogens
mem_null_patho <- dbrda(d_all$d_patho ~ 1, data = mem)
mem_full_patho <- dbrda(d_all$d_patho ~ ., data = mem)
mem_step_patho <- ordistep(mem_null_patho, scope = formula(mem_full_patho), direction = "forward", 
                         permutations = 1999, trace = FALSE)
RsquareAdj(mem_step_patho, permutations = 1999)$adj.r.squared
anova(mem_step_patho, by = "margin", permutations = 1999) %>% 
  as.data.frame() %>% 
  mutate(p.adj = p.adjust(`Pr(>F)`, "fdr")) %>% 
  kable(, format = "pandoc")
#' None
#' 
#' #### Saprotrophs
mem_null_sapro <- dbrda(d_all$d_sapro ~ 1, data = mem)
mem_full_sapro <- dbrda(d_all$d_sapro ~ ., data = mem)
mem_step_sapro <- ordistep(mem_null_sapro, scope = formula(mem_full_sapro), direction = "forward", 
                           permutations = 1999, trace = FALSE)
RsquareAdj(mem_step_sapro, permutations = 1999)$adj.r.squared
anova(mem_step_sapro, by = "margin", permutations = 1999) %>% 
  as.data.frame() %>% 
  mutate(p.adj = p.adjust(`Pr(>F)`, "fdr")) %>% 
  kable(, format = "pandoc")
#' MEM1, MEM2, MEM3
#' Join eigenvectors to sites
sites <- sites %>% left_join(mem %>% rownames_to_column(var = "field_name"), by = join_by(field_name))
#' 
#' ### Wisconsin sites
#' db-MEM
coord_tbl_wi <- sites_wi %>% select(long, lat) %>% as.matrix()
rownames(coord_tbl_wi) <- sites_wi$field_name
mem_wi <- dbmem(coord_tbl_wi) %>% as.data.frame()
setequal(rownames(d_wi$d_its_wi), rownames(mem_wi))
setequal(rownames(d_wi$d_amf_wi), rownames(mem_wi))
setequal(rownames(d_wi$d_patho_wi), rownames(mem_wi))
setequal(rownames(d_wi$d_sapro_wi), rownames(mem_wi))
#' 
#' #### ITS fungi
mem_null_its_wi <- dbrda(d_wi$d_its_wi ~ 1, data = mem_wi)
mem_full_its_wi <- dbrda(d_wi$d_its_wi ~ ., data = mem_wi)
mem_step_its_wi <- ordistep(mem_null_its_wi, scope = formula(mem_full_its_wi), direction = "forward", 
                         permutations = 1999, trace = FALSE)
RsquareAdj(mem_step_its_wi, permutations = 1999)$adj.r.squared
anova(mem_step_its_wi, by = "margin", permutations = 1999) %>% 
  as.data.frame() %>% 
  mutate(p.adj = p.adjust(`Pr(>F)`, "fdr")) %>% 
  kable(, format = "pandoc")
#' MEM2
#' 
#' #### AMF
#' Unifrac distance
mem_null_amf_wi <- dbrda(d_wi$d_amf_wi ~ 1, data = mem_wi)
mem_full_amf_wi <- dbrda(d_wi$d_amf_wi ~ ., data = mem_wi)
mem_step_amf_wi <- ordistep(mem_null_amf_wi, scope = formula(mem_full_amf_wi), direction = "forward", 
                         permutations = 1999, trace = FALSE)
RsquareAdj(mem_step_amf_wi, permutations = 1999)$adj.r.squared
anova(mem_step_amf_wi, by = "margin", permutations = 1999) %>% 
  as.data.frame() %>% 
  mutate(p.adj = p.adjust(`Pr(>F)`, "fdr")) %>% 
  kable(, format = "pandoc")
#' None
#' 
#' #### Pathogens
mem_null_patho_wi <- dbrda(d_wi$d_patho_wi ~ 1, data = mem_wi)
mem_full_patho_wi <- dbrda(d_wi$d_patho_wi ~ ., data = mem_wi)
mem_step_patho_wi <- ordistep(mem_null_patho_wi, scope = formula(mem_full_patho_wi), direction = "forward", 
                           permutations = 1999, trace = FALSE)
RsquareAdj(mem_step_patho_wi, permutations = 1999)$adj.r.squared
anova(mem_step_patho_wi, by = "margin", permutations = 1999) %>% 
  as.data.frame() %>% 
  mutate(p.adj = p.adjust(`Pr(>F)`, "fdr")) %>% 
  kable(, format = "pandoc")
#' MEM2
#' 
#' #### Saprotrophs
mem_null_sapro_wi <- dbrda(d_wi$d_sapro_wi ~ 1, data = mem_wi)
mem_full_sapro_wi <- dbrda(d_wi$d_sapro_wi ~ ., data = mem_wi)
mem_step_sapro_wi <- ordistep(mem_null_sapro_wi, scope = formula(mem_full_sapro_wi), direction = "forward", 
                           permutations = 1999, trace = FALSE)
RsquareAdj(mem_step_sapro_wi, permutations = 1999)$adj.r.squared
anova(mem_step_sapro_wi, by = "margin", permutations = 1999) %>% 
  as.data.frame() %>% 
  mutate(p.adj = p.adjust(`Pr(>F)`, "fdr")) %>% 
  kable(, format = "pandoc")
#' MEM2, MEM1
#' Join eigenvectors to sites
sites_wi <- sites_wi %>% left_join(mem_wi %>% rownames_to_column(var = "field_name"), by = join_by(field_name))
#' #' ## Environmental data
## Env data ———————— ####
#' 
#' ### Plant communities
#' #### Plant functional groups
#' Abundance in functional groups and by species are only available from Wisconsin sites. 
#' Only C4_grass and forbs are used. Others: C3_grass, legume, and shrubTree were found 
#' previously to have high VIF in models or were not chosen in forward selection. 
pfg <- read_csv(root_path("clean_data", "plant_traits.csv"), show_col_types = FALSE) 
#' 
#' #### Plant species and richness
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
#' #### Grass-forb axis
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
gf_axis = scores(pfg_pca, choices = 1, display = "sites") %>% 
  data.frame() %>% 
  rename(gf_axis = PC1) %>% 
  rownames_to_column(var = "field_name")
#' 
#' Are field age and gf_axis correlated?
gfi_yrs <- gf_axis %>% 
  left_join(sites %>% select(field_name, yr_since), by = join_by(field_name)) %>% 
  arrange(-gf_axis)
with(gfi_yrs, cor.test(yr_since, gf_axis))
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
  left_join(gf_axis, by = join_by(field_name)) %>% 
  left_join(fs8_xlab, by = join_by(field_name)) %>% 
  select(field_name, xlab, gf_axis, pl_rich, pl_shan) %>%
  pivot_longer(pl_rich:pl_shan, names_to = "var", values_to = "value") %>% 
  ggplot(aes(x = fct_reorder(xlab, gf_axis), y = value, group = var)) +
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
  left_join(gf_axis, by = join_by(field_name)) %>% 
  filter(field_type != "corn") %>% 
  select(field_name, yr_since, gf_axis, pfg, pct_comp) %>% 
  mutate(pfg = factor(pfg, levels = c("shrubTree", "legume", "C3_grass", "C4_grass", "forb"),
                      labels = c("shrub, tree", "legume", "grass (C3)", "grass (C4)", "forb"))) %>% 
  left_join(fs8_xlab, by = join_by(field_name))
pfg_comp_fig <- 
  ggplot(pfg_comp, aes(x = fct_reorder(xlab, gf_axis), y = pct_comp, group = pfg)) +
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
  left_join(gf_axis, by = join_by(field_name)) %>% 
  filter(field_type != "corn") %>% 
  select(field_name, yr_since, gf_axis, pfg, pct_cvr)  %>% 
  mutate(pfg = factor(pfg, levels = c("C4_grass", "forb"),
                      labels = c("grass (C4)", "forb"))) %>% 
  left_join(fs8_xlab, by = join_by(field_name))
gf_pct_fig <- 
  ggplot(pfg_pct, aes(x = fct_reorder(xlab, gf_axis), y = pct_cvr, group = pfg)) +
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
  ggplot(aes(x = gf_axis, y = rep("PCA 1", nrow(gfi_yrs)))) + 
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.3, color = "gray20") +
  geom_point(aes(color = field_type), shape = 21, size = 1.8, fill = "white", stroke = 0.9) +
  labs(y = NULL, x = "Grass-forb axis") +
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
ggsave(root_path("figs", "figS7.svg"), plot = pfg_pct_fig, device = svglite::svglite,
       width = 7.5, height = 7, units = "in")
#' 
#' ### Soil properties
soil <- read_csv(root_path("clean_data/soil.csv"), show_col_types = FALSE)[-c(26:27), ]
#' 
#' ### Omnibus spe, env, and metadata files
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
  left_join(gf_axis, by = join_by(field_name)) %>%
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
  left_join(gf_axis, by = join_by(field_name)) %>% 
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
  left_join(gf_axis, by = join_by(field_name)) %>% 
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
  left_join(gf_axis, by = join_by(field_name)) %>% 
  left_join(sites %>% select(field_name, field_type, region, yr_since), by = join_by(field_name)) %>% 
  select(field_name, field_type, yr_since, region, everything())

#'  
#' # Composition in guilds
# Composition in guilds ———————— ####
#' 
#' ## Fungi
its_meta %>% 
  count(primary_lifestyle) %>% 
  mutate(composition = round(n / sum(n) * 100, 1)) %>% 
  arrange(-composition) %>% 
  kable(format = "pandoc", caption = "ITS-detectable fungi: composition in guilds")
#' 
#' ## AM fungi
#' Composition in families
amf_meta %>% 
  count(family) %>% 
  mutate(composition = round(n / sum(n) * 100, 1)) %>% 
  arrange(-composition) %>% 
  kable(format = "pandoc", caption = "AM fungi: composition in families")

#' 
#' # Alpha diversity
# Alpha diversity ———————— ####
#' Preprocess data for diversity indices
#+ its_diversity
its_div <- calc_div(its_avg, sites) %>% 
  mutate(depth_csq = sqrt(depth) - mean(sqrt(depth)))
#+ amf_diversity
amf_div <- calc_div(amf_avg, sites) %>% 
  mutate(depth_csq = sqrt(depth) - mean(sqrt(depth)))
#+ patho_diversity
patho_div <- calc_div(patho, sites) %>% 
  mutate(depth_csq = sqrt(depth) - mean(sqrt(depth)))
#+ sapro_diversity
sapro_div <- calc_div(sapro, sites) %>% 
  mutate(depth_csq = sqrt(depth) - mean(sqrt(depth)))

#' 
#' ## Richness
## Richness ———————— ####
#' ### ITS fungi
#' Sequence depth square root transformed and centered. Negative binomial model used to handle count 
#' data. Poisson model was overdispersed (not shown). 
#' 
#' Test interaction
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
#' 
#' ### AM fungi
#' Sequence depth square root transformed and centered. Negative binomial model was underdispersed 
#' and failed to converge at default iterations; use poisson glm instead. 
#' 
#' Test interaction
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
#' don't differ.
#' 
#' ### Pathogens
#' Sequence depth square root transformed and centered. Negative binomial model was underdispersed 
#' and failed to converge at default iterations; use poisson glm instead. 
#' 
#' Test interaction
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
#' 
#' ### Saprotrophs
#' Sequence depth square root transformed and centered. Poisson model was overdispersed (not shown), 
#' use negative binomial instead.  
#' 
#' Test interaction
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

#' 
#' ## Shannon diversity
## Shannon diversity ———————— ####
#' ### ITS fungi
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
#' Shannon diversity in cornfields is significantly less than in restored or remnant fields, which 
#' don't differ.
#' 
#' ### AM fungi
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
#' 
#' ### Pathogens
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
#' 
#' ### Saprotrophs
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
#' 
#' ## Unified results
## Unified results ———————— ####
#' Summary statistics for richness models
#' Fungal OTU richness differences across field types accounting for sequencing depth.
#' Field type effects were evaluated using Type II Analysis of Deviance.
#' P-values for field type were adjusted for multiple comparisons
#' across fungal groups using the Benjamini-Hochberg procedure.
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
  kable(format = "pandoc")
#' 
#' Summary statistics for Shannon models
#' Fungal OTU Shannon diversity differences across field types accounting for sequencing depth.
#' Field type effects were evaluated using Type II Analysis of Variance
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
#' Results summary and figures
div_tagpos <- c(0.18, 1.1)
#+ its_div_fig
its_div_fig <- 
  bind_rows(
    rich = summary(its_rich_em) %>% 
      select(field_type, mean = response, lcl = asymp.LCL, ucl = asymp.UCL),
    shan = summary(its_shan_em) %>% 
      select(field_type, mean = emmean, lcl = lower.CL, ucl = upper.CL),
    .id = "index"
  ) %>% 
  ggplot(aes(x = field_type, y = mean)) +
  geom_col_pattern(
    aes(fill = field_type, pattern = index),
    position = position_dodge(width = div_dodw), width = div_colw, color = "black", linewidth = lw,
    pattern_fill = div_patfil, pattern_colour = div_patcol, pattern_density = div_patden, pattern_spacing = div_patspa
  ) +
  geom_errorbar(aes(ymin = mean, ymax = ucl, group = index), 
                position = position_dodge(width = div_dodw), width = 0, linewidth = lw) +
  geom_text(aes(y = ucl, label = rep(c("a", "b", "b"), 2), group = index), 
            position = position_dodge(width = div_dodw), vjust = -1, family = "sans", size = 3.5) +
  labs(x = NULL) +
  scale_y_continuous(name = expression(atop("Soil fungal", paste("Richness (", italic(n), " OTUs)"))), limits = c(0, 700), 
                     sec.axis = sec_axis(~ . , name = expression(Shannon~diversity~paste("(", italic(e)^italic(H), ")")), breaks = c(0, 100, 200))) +
  scale_pattern_manual(values = c("none", "stripe")) +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = div_tagpos)
#' 
#+ amf_div_fig
amf_div_fig <- 
  bind_rows(
    rich = summary(amf_rich_em) %>% 
      select(field_type, mean = rate, lcl = asymp.LCL, ucl = asymp.UCL),
    shan = summary(amf_shan_em) %>% 
      select(field_type, mean = emmean, lcl = lower.CL, ucl = upper.CL),
    .id = "index"
  ) %>% 
  ggplot(aes(x = field_type, y = mean)) +
  geom_col_pattern(
    aes(fill = field_type, pattern = index),
    position = position_dodge(width = div_dodw), width = div_colw, color = "black", linewidth = lw,
    pattern_fill = div_patfil, pattern_colour = div_patcol, pattern_density = div_patden, pattern_spacing = div_patspa
  ) +
  geom_errorbar(aes(ymin = mean, ymax = ucl, group = index), 
                position = position_dodge(width = div_dodw), width = 0, linewidth = lw) +
  geom_text(aes(y = ucl, label = rep(c("a", "b", "b"), 2), group = index), 
            position = position_dodge(width = div_dodw), vjust = -1, family = "sans", size = 3.5) +
  labs(x = NULL) +
  scale_y_continuous(name = expression(atop("AM fungal", paste("Richness (", italic(n), " OTUs)"))), limits = c(0, 80), 
                     sec.axis = sec_axis(~ . , name = expression(Shannon~diversity~paste("(", italic(e)^italic(H), ")")), breaks = c(0, 15, 30))) +
  scale_pattern_manual(values = c("none", "stripe")) +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = div_tagpos)
#' 
#+ patho_div_fig
patho_div_fig <- 
  bind_rows(
    rich = summary(patho_rich_em) %>% 
      select(field_type, mean = rate, lcl = asymp.LCL, ucl = asymp.UCL),
    shan = summary(patho_shan_em) %>% 
      select(field_type, mean = emmean, lcl = lower.CL, ucl = upper.CL),
    .id = "index"
  ) %>% 
  ggplot(aes(x = field_type, y = mean)) +
  geom_col_pattern(
    aes(fill = field_type, pattern = index),
    position = position_dodge(width = div_dodw), width = div_colw, color = "black", linewidth = lw,
    pattern_fill = div_patfil, pattern_colour = div_patcol, pattern_density = div_patden, pattern_spacing = div_patspa
  ) +
  geom_errorbar(aes(ymin = mean, ymax = ucl, group = index),
                position = position_dodge(width = div_dodw), width = 0, linewidth = lw) +
  labs(x = NULL) +
  scale_y_continuous(name = expression(atop("Pathogen fungal", paste("Richness (", italic(n), " OTUs)"))),  
                     sec.axis = sec_axis(~ . , name = expression(Shannon~diversity~paste("(", italic(e)^italic(H), ")")), breaks = c(0, 5, 10, 15))) +
  scale_pattern_manual(values = c("none", "stripe")) +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = div_tagpos)
#' 
#+ sapro_div_fig
sapro_div_fig <- 
  bind_rows(
    rich = summary(sapro_rich_em) %>% 
      select(field_type, mean = response, lcl = asymp.LCL, ucl = asymp.UCL),
    shan = summary(sapro_shan_em) %>% 
      select(field_type, mean = emmean, lcl = lower.CL, ucl = upper.CL),
    .id = "index"
  ) %>% 
  ggplot(aes(x = field_type, y = mean)) +
  geom_col_pattern(
    aes(fill = field_type, pattern = index),
    position = position_dodge(width = div_dodw), width = div_colw, color = "black", linewidth = lw,
    pattern_fill = div_patfil, pattern_colour = div_patcol, pattern_density = div_patden, pattern_spacing = div_patspa
  ) +
  geom_errorbar(aes(ymin = mean, ymax = ucl, group = index),
                position = position_dodge(width = div_dodw), width = 0, linewidth = lw) +
  labs(x = NULL) +
  scale_y_continuous(name = expression(atop("Soil fungal", paste("Richness (", italic(n), " OTUs)"))),  
                     sec.axis = sec_axis(~ . , name = expression(Shannon~diversity~paste("(", italic(e)^italic(H), ")")), breaks = c(0, 20, 40))) +
  scale_pattern_manual(values = c("none", "stripe")) +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = div_tagpos)
#' 
#' ### Unified figure
#' Display diversity index results
#+ div_fig_patchwork,warning=FALSE
fig2 <- (its_div_fig / plot_spacer() / amf_div_fig / plot_spacer() / patho_div_fig / plot_spacer() / sapro_div_fig) +
  plot_layout(heights = c(rep(c(1, 0.1), 3), 1), axis_titles = "collect") +
  plot_annotation(tag_levels = 'A')
#+ div_fig,warning=FALSE,fig.height=9,fig.width=4.5
fig2
#+ div_fig_save,warning=FALSE,echo=FALSE
ggsave(root_path("figs", "fig2.svg"), plot = fig2, device = svglite::svglite,
       width = 8.5, height = 17, units = "cm")

#' 
#' # Abundance
# Abundance ———————— ####
#' Biomass and abundance-scaled biomass
#' 
#' ## ITS fungi (PLFA)
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
#' 
#' ## AM fungi (NLFA)
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
#' 
#' ## Pathogens
#' Abundance-scaled biomass
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
#' 
#' ## Saprotrophs
#' Abundance-scaled biomass
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
#+ sapro_ma_em_summary,echo=FALSE
kable(summary(sapro_ma_em),
      format = "pandoc",
      caption = "Confidence level used: 0.95")
#+ sapro_ma_fig,fig.width=4,fig.height=4
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
#' ## Unified results
## Unified results ———————— ####
#' Fungal biomass differences differences across field types.
#' Field type effects were evaluated using ANOVA or Analysis of Deviance.
#' P-values for field type were adjusted for multiple comparisons
#' across fungal groups using the Benjamini-Hochberg procedure.
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
#' Figures
#+ plfa_fig,fig.width=4,fig.height=4,fig.align='center'
plfa_fig <- 
  ggplot(summary(plfa_em), aes(x = field_type, y = emmean)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = emmean, ymax = upper.CL), width = 0, linewidth = lw) +
  labs(x = NULL, y = expression(atop("Biomass", paste("(", nmol[PLFA], " × ", g[soil]^{-1}, ")")))) +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1))
#+ nlfa_fig,fig.width=4,fig.height=4,
nlfa_fig <-
  ggplot(summary(nlfa_em), aes(x = field_type, y = response)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = response, ymax = upper.CL), width = 0, linewidth = lw) +
  geom_text(aes(y = upper.CL, label = c("a", "b", "b")),  vjust = -1, family = "sans", size = 3.5) +
  labs(x = NULL, y = expression(atop("Biomass", paste("(", nmol[NLFA], " × ", g[soil]^{-1}, ")")))) +
  scale_fill_manual(values = ft_pal) +
  lims(y = c(0, 75)) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1))
#+ patho_plfa_fig,fig.width=4,fig.height=4,fig.align='center'
patho_ma_fig <- 
  ggplot(summary(patho_ma_em), aes(x = field_type, y = emmean)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = emmean, ymax = upper.CL), width = 0, linewidth = lw) +
  labs(x = NULL, y = expression(atop("Biomass (scaled)", paste(bold(`(`), "(", nmol[PLFA], " × ", g[soil]^{-1}, ")", " × ", paste("(rel. abund)", bold(`)`)))))) +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1))
#+ sapro_ab_fig,fig.width=4,fig.height=4
sapro_ma_fig <- 
  ggplot(summary(sapro_ma_em), aes(x = field_type, y = emmean)) +
  geom_col(aes(fill = field_type), color = "black", width = 0.5, linewidth = lw) +
  geom_errorbar(aes(ymin = emmean, ymax = upper.CL), width = 0, linewidth = lw) +
  labs(x = NULL, y = expression(atop("Biomass (scaled)", paste(bold(`(`), "(", nmol[PLFA], " × ", g[soil]^{-1}, ")", " × ", paste("(rel. abund)", bold(`)`)))))) +
  scale_fill_manual(values = ft_pal) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1))
#' ### Patchwork and export figure
#' Unified figure for supplemental
#+ figS3_patchwork,warning=FALSE
figS3up <- (plfa_fig | plot_spacer() | nlfa_fig) +
  plot_layout(widths = c(0.50, 0.01, 0.50))
figS3dn <- (patho_ma_fig | plot_spacer() | sapro_ma_fig) +
  plot_layout(widths = c(0.50, 0.01, 0.50))
figS3 <- (figS3up / plot_spacer() / figS3dn) +
  plot_layout(heights = c(0.50, 0.01, 0.50)) +
  plot_annotation(tag_levels = 'A')
#+ figS3_fig,warning=FALSE,fig.height=7,fig.width=7
figS3
#+ figS3_save,warning=FALSE,echo=FALSE
ggsave(root_path("figs", "figS3.svg"), plot = figS3, device = svglite::svglite,
       width = 7.5, height = 4.5, units = "in")

#' 
#' # Beta diversity
# Beta diversity ———————— ####
#' PCoA of B-C dissimilarity matrix of relative sequence abundance (row proportion) used for ITS2 
#' OTU sets where biomass did not differ among field types. For AM fungi, the ordination is based on 
#' UNIFRAC distance, and because biomass did differ among field types, the UNIFRAC results are 
#' constrasted with B-C dissimilarity of abundance-scaled biomass. 
#' 
#' Inter-site distance covariate needed for ITS fungi and saprotrophs
#' 
#' ## ITS fungi
#+ its_ord
mva_its <- mva(d = d_all$d_its, env = sites, covar = "MEM1")
#+ its_ord_results
mva_its$dispersion_test
mva_its$permanova
mva_its$pairwise_contrasts[c(1,3,2), c(1,2,4,3,7,8)] %>% 
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
#' ## AM fungi
#' ### Standard ordination
#' Using sequence-based relative abundance, unifrac distance. No inter-site distance covariate.
#+ amf_ord
mva_amf <- mva(d = d_all$d_amf, env = sites, corr = "lingoes")
#+ amf_ord_results
mva_amf$dispersion_test
mva_amf$permanova
mva_amf$pairwise_contrasts[c(1,3,2), c(1,2,4,3,7,8)] %>% 
  arrange(group1, desc(group2)) %>% 
  kable(format = "pandoc", caption = "Pairwise permanova contrasts")
#' Lingoes eigenvalue correction was used. The first three relative eigenvalues exceeded broken stick model. 
#' Based on the homogeneity of variance test, the null hypothesis of equal variance among groups is 
#' accepted across all clusters and in pairwise comparison of clusters (both p>0.05), supporting the application of 
#' a PERMANOVA test. 
#' 
#' Clustering revealed that geographic distance among sites did not significantly explain AMF community variation.
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
  geom_point(data = p_amf_centers, 
             aes(x = mean_Axis.1, y = mean_Axis.2, fill = field_type), 
             size = lg_size, stroke = lw, shape = 21, show.legend = c(fill = FALSE)) +
  geom_point(aes(fill = field_type), size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  scale_fill_manual(name = "Field type", values = ft_pal) +
  labs(
    x = paste0("PCoA 1 (", mva_amf$axis_pct[1], "%)"),
    y = paste0("PCoA 2 (", mva_amf$axis_pct[2], "%)")) +
  theme_ord +
  theme(legend.position = c(0.98, 0.5),
        legend.justification = c(1, 0),
        legend.title = element_text(size = 9, face = 1),
        legend.text = element_text(size = 8, face = 1),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2),
        legend.key = element_rect(fill = "white"),
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
#' 
#' ### Biomass-aware ordination
#' Using abundance-scaled biomass, B-C distance
#+ amf_ord_ma
mva_amf_ma <- mva(d = d_all$d_amf_ma, env = sites, corr = "lingoes")
#+ amf_ord_ma_results
mva_amf_ma$dispersion_test
mva_amf_ma$permanova
mva_amf_ma$pairwise_contrasts[c(1,3,2), c(1,2,4,3,7,8)] %>% 
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
  geom_point(data = p_amf_ma_centers, 
             aes(x = mean_Axis.1, y = mean_Axis.2, fill = field_type), 
             size = lg_size, stroke = lw, shape = 21, show.legend = c(fill = FALSE)) +
  geom_point(aes(fill = field_type), size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25)) +
  scale_fill_manual(name = "Field Type", values = ft_pal) +
  labs(
    x = paste0("PCoA 1 (", mva_amf_ma$axis_pct[1], "%)"),
    y = paste0("PCoA 2 (", mva_amf_ma$axis_pct[2], "%)")) +
  theme_ord +
  theme(legend.title = element_text(size = 9, face = 1),
        legend.text = element_text(size = 8, face = 1))
#' 
#' ### Supplemental figure
#+ figS4,warning=FALSE
amf_ma_ord
#+ figS4_save,warning=FALSE,fig.height=5,fig.width=7,echo=FALSE
ggsave(root_path("figs", "figS6.svg"), plot = amf_ma_ord, device = svglite::svglite,
       width = 5.25, height = 4.25, units = "in")
#' 
#' ### Contrast AMF ordinations
#' Procrustes test on PCoA values using axes with eigenvalues exceeding a broken stick model
#+ amf_protest
set.seed(20251111)
amf_protest <- protest(
  pcoa(d_all$d_amf, correction = "lingoes")$vectors[, 1:3],
  pcoa(d_all$d_amf_ma, correction = "lingoes")$vectors[, 1:3],
  permutations = 1999
)
amf_protest
#' The null that these solutions are unrelated
#' is rejected at p<0.001. However, the alignment isn't perfect. 
#' Clearly, the low biomass in cornfields is a driving difference in 
#' the biomass-aware ordination. Inference would be nearly identical in both cases, 
#' all diagnostics also the same.
#' 
#' ## Pathogens
mva_patho <- mva(d = d_all$d_patho, env = sites, corr = "lingoes")
#' Diagnostics/results
mva_patho$dispersion_test
mva_patho$permanova
mva_patho$pairwise_contrasts[c(1,3,2), c(1,2,4,3,7,8)] %>% 
  arrange(group1, desc(group2)) %>% 
  kable(format = "pandoc", caption = "Pairwise permanova contrasts")
#' 
#' Lingoes correction was needed. Three axes were significant based on a broken stick test. 
#' Based on the homogeneity of variance test, the null hypothesis 
#' of equal variance among groups is accepted across all clusters and in pairwise comparison of 
#' clusters (both p>0.05), supporting the application of a PERMANOVA test.
#' An effect of geographic distance (covariate) on pathogen communities was not supported. 
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
#' ## Saprotrophs
#' Account for spatial effects
#+ sapro_ord
mva_sapro <- mva(d = d_all$d_sapro, env = sites, covar = c("MEM1", "MEM2", "MEM3"))
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
#' ## Beta diversity summary
## Unified results ———————— ####
#' ### Model summary statistics
#' Relative sequence abundance results
#' 
#' Fungal community differences differences among field types.
#' Field type effects were evaluated using Permanova.
#' P-values for field type were adjusted for multiple comparisons
#' across fungal groups using the Benjamini-Hochberg procedure.
#+ unified_permanova_summary,warning=FALSE,message=FALSE
gl_perms <- list(
  its   = mva_its$permanova,
  amf   = mva_amf$permanova,
  patho = mva_patho$permanova,
  sapro = mva_sapro$permanova
) %>% map(\(df) tidy(df) %>% select(term, pseudo_F = statistic, df, R2, p.value))
gl_perms_rdf <- gl_perms %>% 
  map(\(df) df %>% filter(term == "Residual") %>% select(rdf = df)) %>% 
  bind_rows(.id = "guild")
gl_perms %>% 
  bind_rows(.id = "guild") %>% 
  left_join(gl_perms_rdf, by = join_by(guild)) %>% 
  mutate(p.adj = if_else(term == "field_type", p.adjust(p.value, "fdr"), NA_real_),
         across(where(is.numeric), ~ round(.x, 4)),
         `Pseudo_F_(df)` = paste0(pseudo_F, " (", df, " ", rdf, ")")) %>% 
  filter(term %in% c("MEM1", "MEM2", "MEM3", "field_type")) %>% 
  select(guild, term, `Pseudo_F_(df)`, R2, p.value, p.adj) %>% 
  kable(format = "pandoc", caption = "PERMANOVA summary")
#' 
#' Model summary for biomass-aware AM fungi results
#+ permanova_summary,warning=FALSE,message=FALSE
list(amf_ma = mva_amf_ma$permanova) %>% 
  map(\(df) tidy(df) %>% select(term, pseudo_F = statistic, df, R2, p.value)) %>% 
  bind_rows(.id = "guild") %>% 
  mutate(p.adj = if_else(term == "field_type", p.adjust(p.value, "fdr"), NA_real_),
         across(where(is.numeric), ~ round(.x, 4)),
         `Pseudo_F_(df)` = paste0(pseudo_F, " (", df, ", 22)")) %>% 
  filter(term == "field_type") %>% 
  select(guild, term, `Pseudo_F_(df)`, R2, p.value, p.adj) %>% 
  kable(format = "pandoc")
#' 
#' ### Unified figure
#' Display community ordinations
#+ betadiv_patchwork,warning=FALSE
fig3up <- (its_ord | plot_spacer() | amf_ord) +
  plot_layout(widths = c(0.50, 0.01, 0.50))
fig3dn <- (patho_ord | plot_spacer() | sapro_ord) +
  plot_layout(widths = c(0.50, 0.01, 0.50))
fig3 <- (fig3up / plot_spacer() / fig3dn) +
  plot_layout(heights = c(0.50, 0.01, 0.50)) +
  plot_annotation(tag_levels = 'A')
#+ betadiv_fig,warning=FALSE,fig.height=7,fig.width=7
fig3
#+ betadiv_save,warning=FALSE,echo=FALSE
ggsave(root_path("figs", "fig3.svg"), plot = fig3, device = svglite::svglite,
       width = 18, height = 18, units = "cm")

#' 
#' # Fungal communities and the environment
# FungComm-env corr ———————— ####
#' 
#' Do soil properties and plant communities explain variation in fungal communities? What is the 
#' relative explanatory power of each, and which particular variable correlate with fungal 
#' communities?
#' 
#' Restored and remnant prairies in Wisconsin are used to explore these questions. Spatial covariate
#' needed only with saprotrophs.
#' 
#' ## Wrangle explanatory vars
soil_micro_pca <- 
  soil %>% 
  filter(field_name %in% sites_wi$field_name) %>% 
  select(field_name, SO4, Zn, Fe, Mn, Cu, Ca, Mg, Na) %>% 
  column_to_rownames(var = "field_name") %>% 
  decostand(method = "standardize") %>% 
  rda()
summary(soil_micro_pca) # 63% on first two axes
soil_micro_index <- scores(soil_micro_pca, choices = c(1, 2), display = "sites") %>% 
  data.frame() %>% 
  rename(soil_micro_1 = PC1, soil_micro_2 = PC2) %>% 
  rownames_to_column(var = "field_name")
soil_macro <- 
  soil %>% 
  filter(field_name %in% sites_wi$field_name) %>% 
  select(field_name, pH, OM, NO3, P, K)
#' 
#' Assemble explanatory variables and begin iterative selection process. 
#' Plant functional groups and traits not included here were eliminated in previous forward selection
#' procedures (not shown). 
#' Check the VIF for each explanatory variable to test for collinearity if model overfitting is 
#' detected. Then run forward selection in `dbrda()`. 
#' 
env_vars <- sites_wi %>% 
  select(field_name, MEM1, MEM2) %>% 
  left_join(soil_micro_index, by = join_by(field_name)) %>% 
  left_join(soil_macro, by = join_by(field_name)) %>% 
  left_join(gf_axis, by = join_by(field_name)) %>% # 92% on axis 1
  left_join(prich %>% select(field_name, pl_rich), by = join_by(field_name)) %>% # plant richness
  left_join(pfg %>% select(field_name, C3_grass, legume, shrubTree), by = join_by(field_name)) %>% 
  select(-soil_micro_1, -shrubTree, -legume, -C3_grass) %>% # variables removed after VIF check
  column_to_rownames(var = "field_name") %>% 
  as.data.frame()
env_cov <- env_vars[,c("MEM1", "MEM2"), drop = TRUE]
env_expl <- env_vars[, setdiff(colnames(env_vars), c("MEM1", "MEM2")), drop = FALSE] %>% 
  decostand("standardize")
#' Check VIF
env_expl %>% scale() %>% cor() %>% solve() %>% diag() %>% sort() %>% round(2)
#' High VIF or less informative vars iteratively removed with VIF > 10
#' 
#' ## Constrained analyses
## db-RDA ———————— ####
#' 
#' Test explanatory variables for correlation with site ordination. Using plant data, 
#' so the analysis is restricted to Wisconsin sites. Edaphic variables are too numerous 
#' to include individually, so transform micro nutrients using PCA. Forb and grass 
#' cover is highly collinear; use the grass-forb index produced previously with PCA.
#' 
#' Geographic distance covariate was significant with ITS (MEM2), pathogens (MEM2), and
#' saprotrophs (MEM1 & MEM2)
#' 
#' ### ITS fungi
#' Condition MEM2
mod_null <- dbrda(d_wi$d_its_wi ~ 1 + Condition(env_cov[, "MEM2"]), data = env_expl)
mod_full <- dbrda(d_wi$d_its_wi ~ . + Condition(env_cov[, "MEM2"]), data = env_expl)
mod_step <- ordistep(mod_null, 
                     scope = formula(mod_full), 
                     direction = "forward", 
                     permutations = 1999, 
                     trace = FALSE)
#' Results
mod_step
(mod_r2   <- RsquareAdj(mod_step, permutations = 1999))
(mod_glax <- anova(mod_step, permutations = 1999))
(mod_inax <- anova(mod_step, by = "axis", permutations = 1999))
(mod_axpct <- round(100 * mod_step$CCA$eig / sum(mod_step$CCA$eig), 1))
anova(mod_step, by = "margin", permutations = 1999) %>% 
  as.data.frame() %>% 
  mutate(p.adj = p.adjust(`Pr(>F)`, "fdr")) %>% 
  kable(, format = "pandoc")
#' 
#' Create the figure objects. Figure will be produced with panels from other groups. 
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
    mutate(envlabs = c(">forb", "plant spp.")),
  data.frame(
    envvar = "gf_axis",
    dbRDA1 = -mod_scor$biplot["gf_axis", 1],
    dbRDA2 = -mod_scor$biplot["gf_axis", 2],
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
#' ### AM fungi
#' Relative sequence abundance
#' Env covars processed in the ITS section (see above)
amf_mod_null <- dbrda(d_wi$d_amf_wi ~ 1, data = env_expl)
amf_mod_full <- dbrda(d_wi$d_amf_wi ~ ., data = env_expl)
amf_mod_step <- ordistep(amf_mod_null,
                         scope = formula(amf_mod_full),
                         direction = "forward",
                         permutations = 1999,
                         trace = FALSE)
#' Results
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
#' and the selected explanatory variables.
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
    envvar = "gf_axis",
    dbRDA1 = -amf_mod_scor$biplot["gf_axis", 1],
    dbRDA2 = -amf_mod_scor$biplot["gf_axis", 2],
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
#' ### Pathogens
#' Env covars processed in the ITS section (see above)
patho_mod_null <- dbrda(d_wi$d_patho_wi ~ 1 + Condition(env_cov[, "MEM2"]), data = env_expl)
patho_mod_full <- dbrda(d_wi$d_patho_wi ~ . + Condition(env_cov[, "MEM2"]), data = env_expl)
patho_mod_step <- ordistep(patho_mod_null,
                           scope = formula(patho_mod_full),
                           direction = "forward",
                           permutations = 1999,
                           trace = FALSE)
#' Results
patho_mod_step
(patho_mod_r2   <- RsquareAdj(patho_mod_step, permutations = 1999))
(patho_mod_glax <- anova(patho_mod_step, permutations = 1999))
(patho_mod_inax <- anova(patho_mod_step, by = "axis", permutations = 1999))
(patho_mod_axpct <- round(100 * patho_mod_step$CCA$eig / sum(patho_mod_step$CCA$eig), 1))
patho_mod_step$anova %>% 
  as.data.frame() %>% 
  mutate(p.adj = p.adjust(`Pr(>F)`, "fdr")) %>% 
  kable(, format = "pandoc")
#' Based on permutation tests with n=1999 permutations, 
#' after accounting for inter-site pairwise distance as a covariate, the model shows 
#' no significant correlation between pathogen community turnover and explanatory variables.
#' 
#' #### Pathogen constrained figure
patho_mod_step_eig <- c(round(patho_mod_step$CCA$eig * 100, 1), round(patho_mod_step$CA$eig * 100, 1)[1])
patho_mod_scor <- scores(
  patho_mod_step,
  choices = c(1, 2),
  display = c("bp", "sites"),
  tidy = FALSE
)
patho_mod_scor_site <- patho_mod_scor$sites %>%
  data.frame() %>%
  rownames_to_column(var = "field_name") %>%
  left_join(sites, by = join_by(field_name))
patho_mod_scor_bp <- bind_rows(
  patho_mod_scor$biplot %>%
    data.frame() %>%
    rownames_to_column(var = "envvar") %>%
    mutate(envlabs = c("K", ">forb")),
  data.frame(
    envvar = "gf_axis",
    dbRDA1 = -patho_mod_scor$biplot["gf_axis", 1],
    dbRDA2 = -patho_mod_scor$biplot["gf_axis", 2],
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
#' ### Saprotrophs
#' Env covars processed in the ITS section (see above)
#' Two significant spatial vars
sapro_mod_null <- dbrda(d_wi$d_sapro_wi ~ 1 + Condition(MEM1 + MEM2), data = cbind(env_expl, env_cov))
sapro_mod_full <- dbrda(d_wi$d_sapro_wi ~ soil_micro_2 + pH + OM + NO3 + P + K + gf_axis + pl_rich + Condition(MEM1 + MEM2), data = cbind(env_expl, env_cov))
sapro_mod_step <- ordistep(sapro_mod_null,
                           scope = formula(sapro_mod_full),
                           direction = "forward",
                           permutations = 1999,
                           trace = FALSE)
#' Results
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
#' and the selected explanatory variables.
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
    mutate(envlabs = c(">forb", "OM", "plant spp.")),
  data.frame(
    envvar = "gf_axis",
    dbRDA1 = -sapro_mod_scor$biplot["gf_axis", 1],
    dbRDA2 = -sapro_mod_scor$biplot["gf_axis", 2],
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
#' ### Constrained analysis unified summary
#' Environmental drivers were identified via partial distance-based Redundancy Analysis (db-RDA) 
#' using forward selection. 
#' Geographic distance (PCoA Axis 1) was included as a conditional term to partial out spatial effects. 
#' Radj2 represents the cumulative variance explained by the final selected model. 
#' P-values are based on 1,999 permutations; Padj reflects FDR correction within the guild.
#' 
#' Produce objects with explanatory power and degrees of freedom for reporting
#+ dbrda_r2
dbrda_r2 <- data.frame(
  guild = c("all_fungi", "amf", "pathogens", "saprotrophs"),
  r2adj    = round(c(mod_r2$adj.r.squared, amf_mod_r2$adj.r.squared, patho_mod_r2$adj.r.squared, sapro_mod_r2$adj.r.squared), 3)
)
#+ rdf
dbrda_rdf <- data.frame(
  guild = c("all_fungi", "amf", "pathogens", "saprotrophs"),
  rdf   = c(mod_inax["Residual", "Df"], amf_mod_inax["Residual", "Df"], patho_mod_inax["Residual", "Df"], sapro_mod_inax["Residual", "Df"])
)
#' 
#' #### Global tests
#+ dbrda_global_summary,message=FALSE,warning=FALSE
list(
  all_fungi   = mod_glax,
  amf         = amf_mod_glax,
  pathogens   = patho_mod_glax,
  saprotrophs = sapro_mod_glax
) %>% map(\(df) df %>% tidy() %>% filter(term != "Residual")) %>% 
  bind_rows(.id = "guild") %>% 
  left_join(dbrda_rdf, by = join_by(guild)) %>% 
  left_join(dbrda_r2, by = join_by(guild)) %>% 
  mutate(`pseudo_F_(df)` = paste0(round(statistic, 2), " (", df, ", ", rdf, ")"),
         p.adj = p.adjust(p.value, "fdr"),
         across(where(is.numeric), ~ round(.x, 4))) %>% 
  select(guild, term, `pseudo_F_(df)`, r2adj, p.value, p.adj) %>% 
  kable(format = "pandoc")
#' 
#' #### Component axes
#+ dbrda_axis_summary,message=FALSE,warning=FALSE
list(
  all_fungi   = mod_inax,
  amf         = amf_mod_inax,
  pathogens   = patho_mod_inax,
  saprotrophs = sapro_mod_inax
) %>% map(\(df) df %>% tidy() %>% filter(term != "Residual")) %>% 
  bind_rows(.id = "guild") %>% 
  left_join(dbrda_rdf, by = join_by(guild)) %>% 
  mutate(`pseudo_F_(df)` = paste0(round(statistic, 2), " (", df, ", ", rdf, ")"),
         term = str_remove(term, "\\+ "),
         p.adj = p.adjust(p.value, "fdr"),
         across(where(is.numeric), ~ round(.x, 4))) %>% 
  select(guild, term, `pseudo_F_(df)`, p.value, p.adj) %>% 
  kable(format = "pandoc")
#' 
#' #### Selected constraining variables
#+ dbrda_var_summary
list(
  all_fungi   = mod_step$anova,
  amf         = amf_mod_step$anova,
  pathogens   = patho_mod_step$anova,
  saprotrophs = sapro_mod_step$anova
) %>% 
  map(\(df) df %>% tidy()) %>% 
  bind_rows(.id = "guild") %>% 
  left_join(dbrda_rdf, by = join_by(guild)) %>% 
  mutate(statistic = round(statistic, 3),
         `pseudo_F_(df)` = paste0(statistic, " (", df, ", ", rdf, ")"),
         term = str_remove(term, "\\+ "),
         p.adj = p.adjust(p.value, "fdr"),
         across(where(is.numeric), ~ round(.x, 4))) %>% 
  select(guild, term, `pseudo_F_(df)`, p.value, p.adj) %>% 
  arrange(guild, p.value) %>% 
  kable(format = "pandoc")
#' 
#' #### Biplot panels
#' All soil fungi
#+ fig4a
fig4a <- 
  ggplot(mod_scor_site, aes(x = dbRDA1, y = dbRDA2)) +
  geom_segment(data = mod_scor_bp, 
               aes(x = origin, xend = dbRDA1, y = origin, yend = dbRDA2), 
               arrow = arrow(length = unit(2, "mm"), type = "closed"),
               color = c(pfg_col[5], pfg_col[4], "gray20")) +
  geom_text(data = mod_scor_bp, 
            aes(x = labx, y = laby, label = envlabs), 
            size = 3, color = "gray20", fontface = 2) +
  geom_point(aes(fill = field_type), size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  labs(
    x = paste0("db-RDA 1 (", mod_axpct[1], "%)"),
    y = paste0("db-RDA 2 (", mod_axpct[2], "%)")) +
  scale_x_continuous(limits = c(-1.5,1.4), breaks = c(-1, 0, 1)) +
  scale_y_continuous(breaks = c(-1, 0, 1)) +
  scale_fill_manual(values = ft_pal[2:3]) +
  theme_ord +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
#' AMF
#+ fig4b
fig4b <-
  ggplot(amf_mod_scor_site, aes(x = dbRDA1, y = dbRDA2)) +
  geom_segment(data = amf_mod_scor_bp,
               aes(x = origin, xend = dbRDA1, y = origin, yend = dbRDA2),
               arrow = arrow(length = unit(2, "mm"), type = "closed"),
               color = c(pfg_col[5], pfg_col[4], "gray20")) +
  geom_text(data = amf_mod_scor_bp,
            aes(x = labx, y = laby, label = envlabs),
            size = 3, color = "gray20", fontface = 2) +
  geom_point(aes(fill = field_type), size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  labs(
    x = paste0("db-RDA 1 (", amf_mod_axpct[1], "%)"),
    y = paste0("db-RDA 2 (", amf_mod_axpct[2], "%)")) +
  scale_x_continuous(limits = c(-1.3,1.3), breaks = c(-1, 0, 1)) +
  scale_y_continuous(limits = c(-0.95, 1.2), breaks = c(-1, 0, 1)) +
  scale_fill_manual(values = ft_pal[2:3]) +
  theme_ord +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
#' Pathogens, PCoA fig
#+ fig4c
fig4c <-
  ggplot(patho_mod_scor_site, aes(x = -1 * dbRDA1, y = dbRDA2)) +
  geom_segment(data = patho_mod_scor_bp,
               aes(x = origin, xend = -1 * dbRDA1, y = origin, yend = dbRDA2),
               arrow = arrow(length = unit(2, "mm"), type = "closed"),
               color = c("gray20", pfg_col[5], pfg_col[4])) +
  geom_text(data = patho_mod_scor_bp,
            aes(x = -1 * labx, y = laby, label = envlabs),
            size = 3, color = "gray20", fontface = 2) +
  geom_point(aes(fill = field_type), size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  labs(
    x = paste0("db-RDA 1 (", patho_mod_step_eig[1], "%)"),
    y = paste0("db-RDA 2 (", patho_mod_step_eig[2], "%)")) +
  scale_x_continuous(limits = c(-1.3,1.3), breaks = c(-1, 0, 1)) +
  scale_y_continuous(breaks = c(-1, 0, 1)) +
  scale_fill_manual(values = ft_pal[2:3]) +
  theme_ord +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1, hjust = 0),
        plot.tag.position = c(0, 1))
#' Saprotrophs
#+ fig4d
fig4d <-
  ggplot(sapro_mod_scor_site, aes(x = -1 * dbRDA1, y = dbRDA2)) +
  geom_segment(data = sapro_mod_scor_bp,
               aes(x = origin, xend = -1 * dbRDA1, y = origin, yend = dbRDA2),
               arrow = arrow(length = unit(2, "mm"), type = "closed"),
               color = c("gray20", pfg_col[5], pfg_col[4], "gray20")) +
  geom_text(data = sapro_mod_scor_bp,
            aes(x = -1 * labx, y = laby, label = envlabs),
            size = 3, color = "gray20", fontface = 2) +
  geom_point(aes(fill = field_type), size = sm_size, stroke = lw, shape = 21) +
  geom_text(aes(label = yr_since), size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  labs(
    x = paste0("db-RDA 1 (", sapro_mod_axpct[1], "%)"),
    y = paste0("db-RDA 2 (", sapro_mod_axpct[2], "%)")) +
  lims(x = c(-1.5,1.5)) +
  scale_y_continuous(breaks = c(-1, 0, 1)) +
  scale_fill_manual(name = "Field type", values = ft_pal[2:3]) +
  theme_ord +
  theme(legend.position = c(0.98, 0.65),
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
#+ fig4_patchwork,warning=FALSE
fig4up <- (fig4a | plot_spacer() | fig4b) +
  plot_layout(widths = c(0.50, 0.01, 0.50))
fig4dn <- (fig4c | plot_spacer() | fig4d) +
  plot_layout(widths = c(0.50, 0.01, 0.50))
fig4 <- (fig4up / plot_spacer() / fig4dn) +
  plot_layout(heights = c(0.50, 0.01, 0.50)) +
  plot_annotation(tag_levels = 'A')
#+ fig4,warning=FALSE,fig.height=7,fig.width=7
fig4
#' Fungal community ordinations which are constrained or unconstrained by explanatory variables. 
#' Panels show results for all soil fungi **a**, amf **b**, pathogens **c**, and saprotrophs **d**.
#' Percent of constrained (db-RDA) and unconstrained (PCoA) variation explained is shown with axis labels.
#' For explanatory variables with significant community correlations, blue arrows show the grass-forb index 
#' with labels indicating the direction of relative increase in 
#' C4 grasses or forbs, respectively, along the index. The black arrows show other significant constraining
#' variables. Points show locations of restored fields (green) and remnant fields (blue) in Wisconsin. 
#+ fig6_save,warning=FALSE,echo=FALSE
ggsave(root_path("figs", "fig4.svg"), plot = fig4, device = svglite::svglite,
       width = 18, height = 18, units = "cm")

#' 
#' # Fungal abundance and the environment
# FungAbund-env corr ———————— ####
#' 
#' Plant community establishment has varied over time. How do plant communities relate
#' to fungal abundance/proportion in restored and remnant fields?
#' 
#' ## ITS fungi
#' 
#' How variable is biomass across sites?
(its_ma_cv <- 
   sd(fa %>% filter(field_name %in% sites_wi$field_name) %>% pull(fungi_18.2)) / 
   mean(fa %>% filter(field_name %in% sites_wi$field_name) %>% pull(fungi_18.2)) * 100)
#' 
#' Data for tests
fungi_resto <- its_div %>% 
  left_join(fa %>% select(field_name, fungi_mass = fungi_18.2), by = join_by(field_name)) %>% 
  left_join(sites, by = join_by(field_name, field_type)) %>% 
  left_join(gf_axis, by = join_by(field_name)) %>% 
  left_join(prich %>% select(field_name, pl_rich, pl_shan), by = join_by(field_name)) %>% 
  filter(field_type != "corn", region != "FL") %>% 
  select(field_name, fungi_ab = depth, fungi_mass, gf_axis, pl_rich, pl_shan)
#' 
#' ### Plant alpha diversity and fungal biomass
#' Is plant richness related to pathogen mass?
fa_prich_lm <- lm(fungi_mass ~ pl_rich, data = fungi_resto)
summary(fa_prich_lm)
#' Fungal mass and plant richness are weakly correlated but driven by a high-leverage point (not shown).
#' When seq proportion is a response and log(mass) included as a covariate, no relationship 
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
#' Inspect simple linear relationship.
fuma_rest_m <- lm(fungi_mass ~ gf_axis, data = fungi_resto)
summary(fuma_rest_m)
#' The relationship is poor and needs no further analysis
#' 
#' ## AM fungi
#' 
#' How variable is biomass across sites?
(amf_ma_cv <- 
    sd(fa %>% filter(field_name %in% sites_wi$field_name) %>% pull(amf)) / 
    mean(fa %>% filter(field_name %in% sites_wi$field_name) %>% pull(amf)) * 100)
#' 
#' Data for these tests
amf_resto <- amf_div %>% 
  left_join(fa %>% select(field_name, amf_mass = amf), by = join_by(field_name)) %>% 
  left_join(sites, by = join_by(field_name, field_type)) %>% 
  left_join(gf_axis, by = join_by(field_name)) %>% 
  left_join(prich %>% select(field_name, pl_rich, pl_shan), by = join_by(field_name)) %>% 
  filter(field_type != "corn", region != "FL") %>% 
  select(field_name, amf_ab = depth, amf_mass, gf_axis, pl_rich, pl_shan) 
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
amma_rest_m <- lm(amf_mass ~ gf_axis, data = amf_resto)
summary(amma_rest_m)
#' AM fungal mass increases with grass-forb index slightly with p value just 
#' above 0.95 alpha cutoff. 
#' 
#' ## Pathogens
## Pathogens ———————— ####
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
#' How variable is biomass-scaled abundance across sites?
(patho_ma_cv <- sd(patho_resto$patho_mass) / mean(patho_resto$patho_mass) * 100)
#' 
#' ### Plant richness and pathogens
#' Is plant richness related to pathogen mass or proportion?
pathofa_prich_lm = lm(patho_mass ~ pl_rich, data = patho_resto)
summary(pathofa_prich_lm)
#' Pathogen mass and plant richness aren't correlated, though the direction is negative. 
#' Relationship is weak enough that no further tests are warranted.
#' 
patho_prich_glm <- glm(patho_prop ~ fungi_mass_lc + pl_rich,
                    data = patho_resto, family = quasibinomial(link = "logit"),
                    weights = fungi_abund)
summary(patho_prich_glm) 
#' No relationship detected.
#'  
#' Is plant diversity related to pathogen mass or porportion?
pathofa_pshan_lm = lm(patho_mass ~ pl_shan, data = patho_resto)
summary(pathofa_pshan_lm)
#' NS  
#' 
patho_pshan_glm <- glm(patho_prop ~ fungi_mass_lc + pl_shan,
                       data = patho_resto, family = quasibinomial(link = "logit"),
                       weights = fungi_abund)
summary(patho_pshan_glm) 
#' NS
#'  
#' ### Plant functional groups and pathogens
#' #### PFG and pathogen mass
patho_gf_lm <- lm(patho_mass ~ gf_axis, data = patho_resto)
#+ cm9,warning=FALSE,fig.width=7,fig.height=9
check_model(patho_gf_lm) 
#' No obvious issues
summary(patho_gf_lm)
#' The model shows a positive relationship that isn't significant. 
#' Diagnostic reveals noisy fit and lots of structure. 
#' 
#' #### PFG and pathogen proportion
#' Note on interpretation: exponentiated coefficients are interpreted as odds ratios 
#' for pathogen dominance within the fungal community. Sequence abundances are used as 
#' analytic weights so that sites with higher sequencing depth contributed proportionally 
#' more information to the likelihood.
patho_gf_glm <- glm(patho_prop ~ fungi_mass_lc + gf_axis,
                    data = patho_resto, family = quasibinomial(link = "logit"),
                    weights = fungi_abund)
summary(patho_gf_glm) # dispersion parameter >117 justifies quasibinomial
#' Diagnostics
check_model(patho_gf_glm)
check_collinearity(patho_gf_glm)
augment(patho_gf_glm)
#' Long tails and low n showing structure. Moderate leverage at LPRP1: high 
#' pathogens, high gf_axis but very low biomass...this is evidence of the noise that 
#' caused the naïve model to fail. 
distribution_prob(patho_gf_glm)
#' Residuals distribution normal
loocv_paglm_gfi <- map_dbl(seq_len(nrow(patho_resto)), function(i){
  exp(coef(glm(patho_prop ~ fungi_mass_lc + gf_axis, 
           data = patho_resto[-i, ], 
           family = quasibinomial(link = "logit"),
           weights = fungi_abund))["gf_axis"])
})
summary(loocv_paglm_gfi)
(cv_paglm <- (sd(loocv_paglm_gfi) / mean(loocv_paglm_gfi) * 100) %>% round(., 1))
#' Grass-forb index LOOCV variation of 12.6% on the back-transformed scale shows that the influential 
#' points (LPRP1) and three other potential outliers from the qq plot do not 
#' significantly affect fit. Sign and magnitude of LOO slopes wouldn't change inference. 
loocv_paglm_fma <- map_dbl(seq_len(nrow(patho_resto)), function(i){
  exp(coef(glm(patho_prop ~ fungi_mass_lc + gf_axis, 
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
#' View partial regression plots for consistency.
avPlots(patho_gf_glm)
#' Noise in fungal mass data is obvious here. Fit of partial gf_axis is clean. No 
#' non-linear structure is obvious. Both variables seem valuable.
#' 
#' Partial R2 values
#+ parest_m_abs_rsq,warning=FALSE,message=FALSE
data.frame(
  term = c("fungal_mass", "gf_axis"),
  partial_R2 = rsq.partial(patho_gf_glm, adj = TRUE)$partial.rsq
) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% 
  kable(format = "pandoc", caption = "Partial R2 from weighted logistic regression")
#' Model summary
#+ parest_m_abs_summary
patho_null_glm <- glm(patho_prop ~ 1,
                  data = patho_resto, family = quasibinomial(link = "logit"),
                  weights = fungi_abund)
anova(patho_null_glm, patho_gf_glm, test = "F")
#' Deviance explained
(patho_gf_glm_pr2 <- round(1-(summary(patho_gf_glm)$deviance / summary(patho_gf_glm)$null.deviance), 3))
#' Summary of terms
#+ parest_m_abs_terms
tidy(patho_gf_glm) %>% 
  mutate(odds_ratio = exp(estimate), exp_std.error = exp(std.error),
         across(where(is.numeric), ~ round(.x, 3))) %>% 
  select(term, estimate, odds_ratio, std.error, exp_std.error, statistic, p.value) %>% 
  kable(format = "pandoc", caption = "Summary of terms from weighted logistic regression")
(patho_or_pct <- round((exp(coef(patho_gf_glm)[3])^0.1)-1, 3))
#' Confidence intervals on the prediction scale, results on the increment of 0.1 increase
#' in grass-forb index desired due to scale of that variable. Note: in the following output,
#' percent predicted changes are calculated *in excess* of 100% (e.g., 1.124 = 12.4%).
#+ parest_m_abs_confint,warning=FALSE,message=FALSE
((exp(confint(patho_gf_glm))^0.1)-1) %>%
  as.data.frame() %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% 
  kable(format = "pandoc", caption = "95% confidence intervals, back transformed from the log-scale")
#' Create objects for plotting
paglm_med_fungi <- median(patho_resto$fungi_mass_lc, na.rm = TRUE)
paglm_med_abund <- median(patho_resto$fungi_abund, na.rm = TRUE) # Needed for weight context
paglm_newdat <- tibble(
  gf_axis = seq(min(patho_resto$gf_axis, na.rm = TRUE),
                 max(patho_resto$gf_axis, na.rm = TRUE),
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
#' 
#' #### PFG and pathogen species
#' Test which species co-vary with grass-forb axis across sites using a compositionality-aware
#' robust test.
# patho_wi <- guildseq(its_avg, its_meta, "plant_pathogen") %>% # spe matrix
#   left_join(sites %>% select(field_name, field_type, region), by = join_by(field_name)) %>% 
#   filter(field_type != "corn", region != "FL") %>% 
#   select(field_name, where(~ is.numeric(.x) && sum(.x) > 0))
#' Uses function `aldex_gradient()`
#+ patho_aldex_fun,message=FALSE,warning=FALSE
patho_gf_specor <- aldex_gradient(
  spe_tbl = patho_wi,
  covar_tbl = gf_axis,
  covar_col = "gf_axis",
  replicate_multiplier = 10,
  mc.samples = 256,
  denom = "all",
  seed = 20260129
)
#+ patho_aldex_results
patho_gf_spe <- 
  patho_gf_specor$ranked %>% 
  left_join(its_meta %>% 
              select(-otu_ID, -phylum, -primary_lifestyle), 
            by = join_by(otu == otu_num)) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% 
  arrange(rho_p) %>% 
  select(cov_est, rho:species) %>% 
  as_tibble()
#+ patho_aldex20
kable(
  patho_gf_spe %>% filter(abs(rho) >= 0.4) %>% arrange(-rho),
  format = "pandoc", caption = "Pathogen species correlates with grass-forb axis"
)
patho_gf_spe %>% 
  filter(abs(rho) >= 0.4) %>% 
  summarise(
    n_total = n(),
    n_positive = sum(rho > 0),
    n_negative = sum(rho < 0)
  )
#' 
#' ## Saprotrophs
## Saprotrophs ———————— ####
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
#' How variable is biomass-scaled abundance across sites?
(sapro_ma_cv <- sd(sapro_resto$sapro_mass) / mean(sapro_resto$sapro_mass) * 100)
#' 
#' ### Plant richness and saprotrophs
#' Is plant richness related to saprotroph mass?
saprofa_prich_lm <- lm(sapro_mass ~ pl_rich, data = sapro_resto)
distribution_prob(saprofa_prich_lm)
check_model(saprofa_prich_lm)
#' Passes visual diagnostics
summary(saprofa_prich_lm)
#' Strong negative relationship worthy of closer examination. It isn't driven by 
#' grass-forb index (see below), so this isn't a confounding with C4 grass abundance...
#' Residuals distribution normal. 
#' KORP site appears to have some leverage (not shown). With possible leverage, conduct a 
#' LOOCV check.
loocv_sapro_prich_lm <- map_dbl(seq_len(nrow(sapro_resto)), function(i){
  coef(lm(sapro_mass ~ pl_rich, data = sapro_resto[-i, ]))["pl_rich"]
})
summary(loocv_sapro_prich_lm)
(cv_saprich_lm <- (sd(loocv_sapro_prich_lm) / mean(loocv_sapro_prich_lm) * 100) %>% round(., 1) %>% abs(.))
#' All slopes negative and similar in magnitude, with CV = `r cv_saprich_lm`.
#' 
#' Is plant richness related to saprotroph proportion?
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
avPlots(sapro_prich_glm)
#' Noise in fungal mass data is obvious here. Fit of partial gf_axis is clean. No non-linear
#' behavior is obvious, increasing spread with fungal mass expected for model type, 
#' no overdispersion was detected earlier though. 
#' 
#' Partial R2 values
#+ sarest_m_abs_rsq,warning=FALSE,message=FALSE
data.frame(
  term = c("fungal_mass", "pl_rich"),
  partial_R2 = rsq.partial(sapro_prich_glm, adj = TRUE)$partial.rsq
) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% 
  kable(format = "pandoc", caption = "Partial R2 from weighted logistic regression")
#' Model summary
#+ sarest_m_abs_summary
sapro_null_glm <- glm(sapro_prop ~ 1,
                      data = sapro_resto, family = quasibinomial(link = "logit"),
                      weights = fungi_abund)
anova(sapro_null_glm, sapro_prich_glm, test = "F")
#' Deviance explained
(sapro_prich_glm_pr2 <- round(1-(summary(sapro_prich_glm)$deviance / summary(sapro_prich_glm)$null.deviance), 3))
#' Summary of terms
#' Odds ratio prediction and confidence intervals on the prediction scale, results on the increment of an increase
#' of 10 plant species desired due to scale of that variable. Note: in the following output,
#' percent predicted changes are calculated *in excess* of 100% (e.g., 0.085 = -15%).
#+ sarest_m_abs_terms
tidy(sapro_prich_glm) %>% 
  mutate(odds_ratio = exp(estimate), exp_std.error = exp(std.error),
         across(where(is.numeric), ~ round(.x, 3))) %>% 
  select(term, estimate, odds_ratio, std.error, exp_std.error, statistic, p.value) %>% 
  kable(format = "pandoc", caption = "Summary of terms from weighted logistic regression")
(sapro_or_pct <- round(1-(exp(coef(sapro_prich_glm)[3])^10), 3))
#+ sarest_m_abs_confint,warning=FALSE,message=FALSE
(1-(exp(confint(sapro_prich_glm))^10)) %>%
  as.data.frame() %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% 
  kable(format = "pandoc", caption = "95% confidence intervals, back transformed from the log-scale")
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
#' 
#' #### Plant richness and saprotroph species
#' Identify saprotroph species that co-vary with richness across sites.
sapro_wi <- guildseq(its_avg, its_meta, "saprotroph") %>% # spe matrix
  left_join(sites %>% select(field_name, field_type, region), by = join_by(field_name)) %>% 
  filter(field_type != "corn", region != "FL") %>% 
  select(field_name, where(~ is.numeric(.x) && sum(.x) > 0)) # Back-transform to field sums
#' Using function `aldex_gradient`.
#+ sapro_aldex_fun,message=FALSE,warning=FALSE
sapro_rich_specor <- aldex_gradient(
  spe_tbl = sapro_wi,
  covar_tbl = prich %>% select(field_name, pl_rich),
  covar_col = "pl_rich",
  replicate_multiplier = 10,
  mc.samples = 256,
  denom = "all",
  seed = 20260129
)
#+ sapro_aldex_results
sapro_rich_spe <- 
  sapro_rich_specor$ranked %>% 
  left_join(its_meta %>% 
              select(-otu_ID, -phylum, -primary_lifestyle), 
            by = join_by(otu == otu_num)) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% 
  arrange(rho_p) %>% 
  select(cov_est, rho:species) %>% 
  as_tibble()
#+ sapro_aldex20
kable(
  sapro_rich_spe %>% filter(abs(rho) >= 0.4) %>% arrange(rho),
  format = "pandoc", caption = "Saprotroph species correlates with plant richness"
)
sapro_rich_spe %>% 
  filter(abs(rho) >= 0.4) %>% 
  summarise(
    n_total = n(),
    n_positive = sum(rho > 0),
    n_negative = sum(rho < 0)
  )
sapro_rich_spe %>% filter(abs(rho) >= 0.4) %>% 
  mutate(sign = case_when(rho < 0 ~ "negative", rho > 0 ~ "positive")) %>% 
  group_by(sign, class) %>% 
  count() %>% 
  pivot_wider(names_from = "sign", values_from = "n") %>% 
  kable(format = "pandoc", caption = "Saprotroph classes and correlations with plant richness")
#' 
#' ### Plant diversity and saprotrophs
#' Is plant diversity related to saprotroph mass?
saprofa_pshan_lm <- lm(sapro_mass ~ pl_shan, data = sapro_resto)
distribution_prob(saprofa_pshan_lm)
check_model(saprofa_pshan_lm)
summary(saprofa_pshan_lm)
#' Saprotroph mass has a weak, significant relationship with saprotroph mass. 
#' Two extreme values exist, making a bit of a dumbbell shape. Residuals distribution
#' normal. Two points with borderline leverage, conduct LOOCV. 
loocv_sapshlm_fma <- map_dbl(seq_len(nrow(sapro_resto)), function(i){
  coef(lm(sapro_mass ~ pl_shan, data = sapro_resto[-i, ]))["pl_shan"]
})
summary(loocv_sapshlm_fma)
(cv_sapshlm_fma <- (sd(loocv_sapshlm_fma) / mean(loocv_sapshlm_fma) * 100) %>% round(., 1) %>% abs(.))
#' Model relies heavily on two extreme values, with a LOOCV variation of `r cv_sapshlm_fma`%. 
#' Suggest that makes the model unimportant. All slopes negative, but many very slight.
#' Is plant diversity related to saprotroph proportion?
sapro_pshan_glm <- glm(sapro_prop ~ fungi_mass_lc + pl_shan,
                       data = sapro_resto, family = quasibinomial(link = "logit"),
                       weights = fungi_abund)
summary(sapro_pshan_glm)
#' NS
#' 
#' ### Plant functional groups and saprotrophs
sama_rest_m <- lm(sapro_mass ~ gf_axis, data = sapro_resto)
summary(sama_rest_m)
#' Likely not enough of a relationship to warrant further attention
#' Examine mass+richness models, similar to pathogen models
sapro_gf_glm <- glm(sapro_prop ~ fungi_mass_lc + gf_axis,
                    data = sapro_resto, family = quasibinomial(link = "logit"),
                    weights = fungi_abund)
summary(sapro_gf_glm)
#' NS
#' 
#' ### Guild-plant relationships
## Unified results ———————— ####
#' Create multipanel figure, post-production in editing software will be necessary.
#+ fig7a,warning=FALSE
fig5a <-
  ggplot(paglm_pred, aes(x = gf_axis, y = fit_prob)) +
  geom_line(color = "black", linewidth = lw) +
  geom_point(data = patho_resto, aes(x = gf_axis, y = patho_prop, fill = field_type),
             size = sm_size, stroke = lw, shape = 21) +
  geom_text(data = patho_resto, aes(x = gf_axis, y = patho_prop, label = yr_since),
            size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  labs(
    x = "Grass–forb axis",
    y = "Pathogen proportion",
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
#+ fig7a_rug_data
gfa_fgc <- # grass-forb axis, forb-grass composition
  pfg_pct %>% 
  group_by(field_name) %>% 
  mutate(pct_comp = pct_cvr / sum(pct_cvr)) %>% 
  filter(pfg == "forb") %>% 
  select(field_name, gf_axis, forb_comp = pct_comp) %>% 
  arrange(gf_axis)
#+ fig7a_rug,warning=FALSE
fig5a_rug <- add_fig7_rug(
  fig5a,
  comp_df = gfa_fgc,
  y0 = 0.05,
  h  = 0.010,
  forb_fill  = pfg_col[4],
  grass_fill = pfg_col[5]
) +
  expand_limits(y = 0.05) +
  geom_text(data = data.frame(x = c(-0.45, 0.45), y = c(0.04, 0.04),
                              lab = c(paste0("bold(grass~(C[4]))"), paste0("bold(forb)"))),
            aes(x = x, y = y, label = lab), parse = TRUE, size = 2.8, family = "Helvetica")
#+ fig7b,warning=FALSE
fig5b <-
  ggplot(saglm_pred, aes(x = pl_rich, y = fit_prob)) +
  geom_line(color = "black", linewidth = lw) +
  geom_point(data = sapro_resto, aes(x = pl_rich, y = sapro_prop, fill = field_type),
             size = sm_size, stroke = lw, shape = 21) +
  geom_text(data = sapro_resto, aes(x = pl_rich, y = sapro_prop, label = yr_since),
            size = yrtx_size, family = "sans", fontface = 2, color = "black") +
  labs(
    x = expression(paste("Plant richness (", italic(n), " species)")),
    y = "Saprotroph proportion",
    tag = "A"
  ) +
  scale_fill_manual(name = "Field type", values = ft_pal[2:3]) +
  theme_cor +
  theme(legend.position = "none",
        plot.tag = element_text(size = 14, face = 1),
        plot.tag.position = c(0, 1))
#+ fig7_patchwork
fig5 <- (fig5a_rug | plot_spacer() | fig5b) +
  plot_layout(widths = c(0.50, 0.01, 0.50), axis_titles = "collect_y") +
  plot_annotation(tag_levels = 'A')
#+ fig7_display,warning=FALSE
fig5
#+ fig7_save,warning=FALSE,echo=FALSE
ggsave(root_path("figs", "fig5.svg"), plot = fig5, device = svglite::svglite,
       width = 18, height = 9, units = "cm")
#' 
