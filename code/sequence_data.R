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

#' # Resources

#' ## Packages
packages_needed <- c("tidyverse", "vegan", "knitr", "colorspace", "plotrix", "rprojroot", "rlang")
if (any(!packages_needed %in% rownames(installed.packages())))
    install.packages(packages_needed[!packages_needed %in% rownames(installed.packages())])
invisible(lapply(packages_needed, library, character.only = TRUE))

#' ## Styles
#+ graphics_styles
source(root_path("resources", "styles.txt"))

#' # Functions

#' ## Root path function
root_path <- function(...) rprojroot::find_rstudio_root_file(...)

#' ## ETL: clean OTU data and return formatted objects
#+ function_etl
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

#' ## Species accumulation
#+ function_spe_accum
spe_accum <- function(data) {
    df <- data.frame(
        samples = specaccum(data[, -c(1, 2)], method = "exact")$site,
        richness = specaccum(data[, -c(1, 2)], method = "exact")$richness,
        sd = specaccum(data[, -c(1, 2)], method = "exact")$sd
    )
    df
}

#' ## Confidence interval helper
#+ function_confint
ci <- function(x) std.error(x) * qnorm(0.975)

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
    mutate(field_type = factor(field_type, ordered = TRUE, levels = c("corn", "restored", "remnant"))) %>%
    select(-lat, -long, -yr_restore)

#' ## ETL processing
#+ process_its,message=FALSE
its <- etl(spe = its_otu, taxa = its_taxa, traits = traits, varname = "otu_num", gene = "ITS",
           colname_prefix = "ITS_TGP_", folder = "clean_data")

#+ process_18S,message=FALSE,warning=FALSE
amf <- etl(spe = amf_otu, taxa = amf_taxa, varname = "otu_num", gene = "18S",
           colname_prefix = "18S_TGP_", folder = "clean_data")

#' ## UNIFRAC tables for AMF
#+ wrangle_amf_samps
amf$spe_samps %>%
    mutate(field_sample = paste(field_name, sample, sep = "_")) %>%
    select(field_sample, everything(), -field_name, -sample) %>%
    column_to_rownames("field_sample") %>%
    t() %>% as.data.frame() %>% rownames_to_column("otu_num") %>%
    left_join(amf$spe_meta %>% select(otu_num, otu_ID), by = "otu_num") %>%
    select(otu_ID, everything(), -otu_num) %>%
    write_tsv(root_path("otu_tables/18S/18S_samps_4unifrac.tsv"))

#+ wrangle_amf_avg
amf$spe_avg %>%
    column_to_rownames("field_name") %>%
    t() %>% as.data.frame() %>% rownames_to_column("otu_num") %>%
    left_join(amf$spe_meta %>% select(otu_num, otu_ID), by = "otu_num") %>%
    select(otu_ID, everything(), -otu_num) %>%
    write_tsv(root_path("otu_tables/18S/18S_avg_4unifrac.tsv"))

#' # Sampling depth and coverage
#' Script running `rarecurve()` is commented out because it takes so long to execute.
#' Data were saved to the wd and are used for making figures. These files are too large
#' to upload to GitHub and are ignored. Please run these calls to `rarecurve()` to create
#' your own rarefaction and species accumulation data files. 

#' ## Rarefaction: ITS sample
# its_rc <- rarecurve(
#     its$spe_samps %>%
#         mutate(field_sample = paste(field_name, sample, sep = "_")) %>%
#         column_to_rownames("field_sample") %>%
#         select(-field_name, -sample),
#     step = 1, tidy = TRUE)
# write_csv(its_rc, root_path("clean_data", "its_rare_samp.csv"))

its_rc <- read_csv(root_path("clean_data", "its_rare_samp.csv"), show_col_types = FALSE)

#+ its_rarefaction,fig.width=6,fig.height=8
its_rc %>%
    separate_wider_delim(Site, delim = "_", names = c("field_name", "sample_key"), cols_remove = FALSE) %>% 
    rename(seq_abund = Sample, otus = Species, field_sample = Site) %>%
    left_join(sites, by = "field_name") %>%
    ggplot(aes(x = seq_abund, y = otus, group = field_sample)) +
    facet_wrap(vars(field_type), ncol = 1, scales = "free") +
    geom_line(aes(color = field_type), linewidth = 0.4) +
    scale_color_discrete_qualitative(palette = "Harmonic") +
    labs(x = "Sequence abundance", y = "OTUs", title = "Rarefaction of ITS samples") +
    theme_bw()

#' ## Rarefaction: ITS site-averaged
# its_rc_site <- rarecurve(
#     its$spe_samps %>%
#         group_by(field_name) %>%
#         summarize(across(starts_with("otu"), sum), .groups = "drop") %>%
#         column_to_rownames("field_name"),
#     step = 1, tidy = TRUE)
# write_csv(its_rc_site, root_path("clean_data", "its_rare_site.csv"))

its_rc_site <- read_csv(root_path("clean_data", "its_rare_site.csv"), show_col_types = FALSE)

#+ its_rarefaction_site_avg,fig.width=6,fig.height=8
its_rc_site %>%
    rename(seq_abund = Sample, otus = Species, field_name = Site) %>%
    left_join(sites, by = "field_name") %>%
    ggplot(aes(x = seq_abund, y = otus, group = field_name)) +
    facet_wrap(vars(field_type), ncol = 1, scales = "free_y") +
    geom_line(aes(color = field_type), linewidth = 0.4) +
    scale_color_discrete_qualitative(palette = "Harmonic") +
    labs(x = "Sequence abundance", y = "OTUs", title = "Rarefaction of ITS (site-averaged)") +
    theme_bw()

#' ## Rarefaction: AMF sample
#+ amf_rc,message=FALSE,warning=FALSE
# amf_rc <- rarecurve(
#     amf$spe_samps %>%
#         mutate(field_sample = paste(field_name, sample, sep = "_")) %>%
#         column_to_rownames("field_sample") %>%
#         select(-field_name, -sample),
#     step = 1, tidy = TRUE)
# write_csv(amf_rc, root_path("clean_data", "amf_rare_samp.csv"))

amf_rc <- read_csv(root_path("clean_data", "amf_rare_samp.csv"), show_col_types = FALSE)

#+ amf_rarefaction,fig.width=6,fig.height=8
amf_rc %>%
    separate_wider_delim(Site, delim = "_", names = c("field_name", "sample_key"), cols_remove = FALSE) %>% 
    rename(seq_abund = Sample, otus = Species, field_sample = Site) %>%
    left_join(sites, by = "field_name") %>%
    ggplot(aes(x = seq_abund, y = otus, group = field_sample)) +
    facet_wrap(vars(field_type), ncol = 1, scales = "free") +
    geom_line(aes(color = field_type), linewidth = 0.4) +
    scale_color_discrete_qualitative(palette = "Harmonic") +
    labs(x = "Sequence abundance", y = "OTUs", title = "Rarefaction of 18S samples") +
    theme_bw()

#' ## Rarefaction: AMF site-averaged
#+ amf_rc_site,message=FALSE,warning=FALSE
# amf_rc_site <- rarecurve(
#     amf$spe_samps %>%
#         group_by(field_name) %>%
#         summarize(across(starts_with("otu"), sum), .groups = "drop") %>%
#         column_to_rownames("field_name"),
#     step = 1, tidy = TRUE)
# write_csv(amf_rc_site, root_path("clean_data", "amf_rare_site.csv"))

amf_rc_site <- read_csv(root_path("clean_data", "amf_rare_site.csv"), show_col_types = FALSE)

#+ amf_rarefaction_site_avg,fig.width=6,fig.height=8
amf_rc_site %>%
    rename(seq_abund = Sample, otus = Species, field_name = Site) %>%
    left_join(sites, by = "field_name") %>%
    ggplot(aes(x = seq_abund, y = otus, group = field_name)) +
    facet_wrap(vars(field_type), ncol = 1, scales = "free_y") +
    geom_line(aes(color = field_type), linewidth = 0.4) +
    scale_color_discrete_qualitative(palette = "Harmonic") +
    labs(x = "Sequence abundance", y = "OTUs", title = "Rarefaction of 18S (site-averaged)") +
    theme_bw()

#' # Species accumulation curves
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

#+ species_accumulation_fig,fig.width=8,fig.height=5
ggplot(accum, aes(x = samples, y = richness, group = field_name)) +
    facet_wrap(vars(dataset), scales = "free_y") +
    geom_line(aes(color = field_type)) +
    geom_segment(aes(xend = samples, y = richness - sd, yend = richness + sd, color = field_type)) +
    scale_color_discrete_qualitative(palette = "Harmonic") +
    labs(
        x = "Samples",
        y = expression(N[0]),
        caption = "Species accumulation using exact method; error = moment-based SD"
    ) +
    scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
    theme_bw()











