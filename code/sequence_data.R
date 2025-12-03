#' ---
#' title: "Species Data: ETL and Diagnostics"
#' author: "Beau Larkin\n"
#' date: "Last updated: `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   github_document:
#'     toc: true
#'     toc_depth: 3
#'     df_print: paged
#' ---
#' 
#' # Description
#' 
#' - Load and clean QIIME2 sequence data
#' - Apply fungal traits
#' - Create sample and site OTU tables
#' - Export UNIFRAC tables for AMF
#' - Evaluate sampling effort with rarefaction and accumulation
#' 
#' **Note:** individual sample based data needed to produce collection curves.
#'
#' # Resources

#' ## Packages
packages_needed <- c("tidyverse", "vegan", "knitr", "colorspace", "plotrix", "rprojroot", 
                     "rlang", "patchwork", "cowplot")

to_install <- setdiff(packages_needed, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)
invisible(lapply(packages_needed, library, character.only = TRUE))
#' 
#' # Functions
#' ## Root path function
root_path <- function(...) rprojroot::find_rstudio_root_file(...)
#' 
#' Repo functions loaded from a separate script to save lines here; to view the function navigate to 
#' `functions.R` in the code folder, accessible from the root dir of the repo.
source(root_path("code", "functions.R"))
#' 
#' ## Styles
#+ graphics_styles
source(root_path("resources", "styles.R"))


#' 
#' # Load and process data

#' ## Import data
#+ species_taxa_traits,message=FALSE
its_otu <- read_delim(root_path("otu_tables/ITS/ITS_otu_raw.txt"), show_col_types = FALSE)
its_taxa <- read_delim(root_path("otu_tables/ITS/ITS_otu_taxonomy.txt"), show_col_types = FALSE)
amf_otu <- read_delim(root_path("otu_tables/18S/18S_otu_raw.txt"), show_col_types = FALSE) %>% select(-last_col())
amf_taxa <- read_delim(root_path("otu_tables/18S/18S_otu_taxonomy.txt"), show_col_types = FALSE)
traits <- read_csv(root_path("otu_tables/2023-02-23_fungal_traits.csv"), show_col_types = FALSE) %>%
    select(phylum:primary_lifestyle)

#+ import_sites,message=FALSE
sites <- read_csv(root_path("clean_data/sites.csv"), show_col_types = FALSE) %>%
    mutate(field_type = factor(field_type, levels = c("corn", "restored", "remnant"))) %>%
    select(-lat, -long, -yr_restore)

#' ## ETL processing
#+ process_its,message=FALSE
its <- etl(spe = its_otu, taxa = its_taxa, traits = traits, varname = "otu_num", gene = "ITS",
           colname_prefix = "ITS_TGP_", folder = "clean_data")

#+ process_18S,message=FALSE,warning=FALSE
amf <- etl(spe = amf_otu, taxa = amf_taxa, varname = "otu_num", gene = "18S",
           colname_prefix = "18S_TGP_", folder = "clean_data")

#' # Sampling depth and coverage
#' Script running `rarecurve()` is commented out because it takes so long to execute.
#' Data were saved to the wd and are used for making figures. These files are too large
#' to upload to GitHub and are ignored. Please run the calls to `rarecurve()` to create
#' your own rarefaction and species accumulation data files. 

#' ## Rarefaction: ITS sample
# its_rc <- rarecurve(
#     its$spe_samps %>%
#         mutate(field_sample = paste(field_name, sample, sep = "_")) %>%
#         column_to_rownames("field_sample") %>%
#         select(-field_name, -sample),
#     step = 1, tidy = TRUE)
# write_csv(its_rc, root_path("clean_data", "its_rare_samp.csv"))

#' Read in the data already produced with `rarecurve()`.
its_rc <- read_csv(root_path("clean_data", "its_rare_samp.csv"), show_col_types = FALSE)

#+ its_rarefaction,fig.width=4,fig.height=7
its_rc %>%
    separate_wider_delim(Site, delim = "_", names = c("field_name", "sample_key"), cols_remove = FALSE) %>% 
    rename(seq_abund = Sample, otus = Species, field_sample = Site) %>%
    left_join(sites, by = "field_name") %>%
    ggplot(aes(x = seq_abund, y = otus, group = field_sample)) +
    facet_wrap(vars(field_type), ncol = 1, scales = "free") +
    geom_line(aes(color = field_type), linewidth = 0.4) +
    scale_color_discrete_qualitative(palette = "Harmonic") +
    labs(x = "Sequence abundance", y = "OTUs", title = "Rarefaction of ITS samples") +
    theme_corf +
    theme(legend.position = "none")

#' ## Rarefaction: ITS site-averaged
# its_rc_site <- rarecurve(
#     its$spe_samps %>%
#         group_by(field_name) %>%
#         summarize(across(starts_with("otu"), sum), .groups = "drop") %>%
#         column_to_rownames("field_name"),
#     step = 1, tidy = TRUE)
# write_csv(its_rc_site, root_path("clean_data", "its_rare_site.csv"))

#' Read in the data already produced by `rarecurve()`.
its_rc_site <- read_csv(root_path("clean_data", "its_rare_site.csv"), show_col_types = FALSE)

#+ its_rarefaction_site_avg,fig.width=4,fig.height=7
its_rc_site %>%
    rename(seq_abund = Sample, otus = Species, field_name = Site) %>%
    left_join(sites, by = "field_name") %>%
    ggplot(aes(x = seq_abund, y = otus, group = field_name)) +
    facet_wrap(vars(field_type), ncol = 1, scales = "free_y") +
    geom_line(aes(color = field_type), linewidth = 0.4) +
  scale_color_discrete_qualitative(palette = "Harmonic") +
    labs(x = "Sequence abundance", y = "OTUs", title = "Rarefaction of ITS (site-averaged)") +
    theme_corf +
    theme(legend.position = "none")

#' ## Rarefaction: AMF sample
#+ amf_rc,message=FALSE,warning=FALSE
# amf_rc <- rarecurve(
#     amf$spe_samps %>%
#         mutate(field_sample = paste(field_name, sample, sep = "_")) %>%
#         column_to_rownames("field_sample") %>%
#         select(-field_name, -sample),
#     step = 1, tidy = TRUE)
# write_csv(amf_rc, root_path("clean_data", "amf_rare_samp.csv"))

#' Read in data produced by `rarecurve()`.
amf_rc <- read_csv(root_path("clean_data", "amf_rare_samp.csv"), show_col_types = FALSE)

#+ amf_rarefaction,fig.width=4,fig.height=7
amf_rc %>%
    separate_wider_delim(Site, delim = "_", names = c("field_name", "sample_key"), cols_remove = FALSE) %>% 
    rename(seq_abund = Sample, otus = Species, field_sample = Site) %>%
    left_join(sites, by = "field_name") %>%
    ggplot(aes(x = seq_abund, y = otus, group = field_sample)) +
    facet_wrap(vars(field_type), ncol = 1, scales = "free") +
    geom_line(aes(color = field_type), linewidth = 0.4) +
  scale_color_discrete_qualitative(palette = "Harmonic") +
    labs(x = "Sequence abundance", y = "OTUs", title = "Rarefaction of 18S samples") +
    theme_corf +
    theme(legend.position = "none")

#' ## Rarefaction: AMF site-averaged
#+ amf_rc_site,message=FALSE,warning=FALSE
# amf_rc_site <- rarecurve(
#     amf$spe_samps %>%
#         group_by(field_name) %>%
#         summarize(across(starts_with("otu"), sum), .groups = "drop") %>%
#         column_to_rownames("field_name"),
#     step = 1, tidy = TRUE)
# write_csv(amf_rc_site, root_path("clean_data", "amf_rare_site.csv"))

#' Read in data already produced by `rarecurve()`.
amf_rc_site <- read_csv(root_path("clean_data", "amf_rare_site.csv"), show_col_types = FALSE)

#+ amf_rarefaction_site_avg,fig.width=4,fig.height=7
amf_rc_site %>%
    rename(seq_abund = Sample, otus = Species, field_name = Site) %>%
    left_join(sites, by = "field_name") %>%
    ggplot(aes(x = seq_abund, y = otus, group = field_name)) +
    facet_wrap(vars(field_type), ncol = 1, scales = "free_y") +
    geom_line(aes(color = field_type), linewidth = 0.4) +
  scale_color_discrete_qualitative(palette = "Harmonic") +
    labs(x = "Sequence abundance", y = "OTUs", title = "Rarefaction of 18S (site-averaged)") +
    theme_corf +
    theme(legend.position = "none")

#' # Figure data
#+ species_accumulation,message=FALSE,warning=FALSE
accum <- bind_rows(
    list(
        ITS = bind_rows(split(its$spe_samps, ~ field_name) %>% map(spe_accum), .id = "field_name"),
        AMF = bind_rows(split(amf$spe_samps, ~ field_name) %>% map(spe_accum), .id = "field_name")
    ),
    .id = "dataset"
) %>%
    mutate(dataset = factor(dataset, levels = c("ITS", "AMF"), ordered = TRUE)) %>%
    left_join(sites, by = "field_name")
#+ rarefaction_unified,message=FALSE,warning=FALSE
rarefac <- bind_rows(
  list(
    ITS = its_rc_site,
    AMF = amf_rc_site
  ),
  .id = "dataset"
) %>%
  select(field_name = Site, everything()) %>% 
  mutate(dataset = factor(dataset, levels = c("ITS", "AMF"), ordered = TRUE)) %>%
  left_join(sites, by = "field_name")

#' # Output figures
its_rare_fig <- ggplot(rarefac %>% filter(dataset == "ITS"), aes(x = Sample, y = Species, group = field_name)) +
  facet_grid(rows = vars(dataset), cols = vars(field_type), scales = "fixed") +
  geom_line(aes(color = field_type)) +
  scale_color_manual(values = ft_pal) +
  labs(
    x = NULL,
    y = NULL) +
  theme_corf +
  scale_x_continuous(breaks = seq(0, 90000, 30000)) +
  theme(legend.position = "none",
        plot.margin = unit(rep(0,4), "mm"))
amf_rare_fig <- ggplot(rarefac %>% filter(dataset == "AMF"), aes(x = Sample, y = Species, group = field_name)) +
  facet_grid(rows = vars(dataset), cols = vars(field_type), scales = "fixed") +
  geom_line(aes(color = field_type)) +
  scale_color_manual(values = ft_pal) +
  labs(
    x = NULL,
    y = NULL) +
  theme_corf +
  scale_x_continuous(breaks = seq(0, 45000, 15000)) +
  theme(legend.position = "none",
        strip.text.x = element_blank(),   
        strip.background.x = element_blank(),
        plot.margin = unit(rep(0,4), "mm"))
#' Create figure panels
rare_panels <- (its_rare_fig / plot_spacer() / amf_rare_fig) +
  plot_layout(heights = c(1,0.01,1))
x_lab <- ggdraw() + 
  draw_label("Sequence abundance", hjust = 0.5, vjust = 0.5, size = 10) 
y_lab <- ggdraw() + 
  draw_label("OTUs", angle = 90, hjust = 0.5, vjust = 0.5, size = 10)
rare_fig_h <- (y_lab | rare_panels) + plot_layout(widths = c(0.03, 1))
rare_fig <- rare_fig_h / x_lab + plot_layout(heights = c(1, 0.10))
rare_fig <- rare_fig + plot_annotation(
  theme = theme(plot.margin = unit(rep(0,4), "mm"))
)
#+ rarefaction_fig,fig.width=7,fig.height=4
rare_fig
#+ rarefaction_fig_save
ggsave(
  root_path("figs", "figS1.png"),
  plot = rare_fig,
  height = 4.24,
  width = 7.5,
  units = "in",
  dpi = 600
)

#+ species_accumulation_fig,fig.width=7,fig.height=4
accum_fig <- 
  ggplot(accum, aes(x = samples, y = richness, group = field_name)) +
  facet_grid(rows = vars(dataset), cols = vars(field_type), scales = "free_y") +
  geom_line(aes(color = field_type)) +
  geom_segment(aes(xend = samples, y = richness - sd, yend = richness + sd, color = field_type)) +
  scale_color_manual(values = ft_pal) +
  labs(x = "Samples", y = "OTUs") +
    scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
    theme_corf +
    theme(legend.position = "none", 
          plot.margin = unit(c(0,2,4,2), "mm"))
#+ accum_fig_save
ggsave(
  root_path("figs", "figS2.png"),
  plot = accum_fig,
  height = 4.25,
  width = 7.5,
  units = "in",
  dpi = 600
)





