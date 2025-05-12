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
#' REVISE into more of a TOC than explanations
#' 
#' This is the sequence data only
#' 
#' Microbial sequence abundances were produced by Lorinda Bullington in QIIME2. Files output from
#' QIIME2 need to be extracted, transformed, and loaded (ETL) for downstream analysis. Once ETL
#' is complete, the resulting samples-species matrices must be examined for sampling depth and 
#' community coverage. 
#' 
#' ## Source data
#' Sequences for fungi (ITS sequence based identification) and AMF (18S sequence based 
#' identification) were grouped into in 97% similar OTUs. Taxonomy and metadata are included as
#' separate files. Singletons were removed in QIIME2. 
#' A few samples failed to amplify. Most sites are characterized by 10 samples, but a few
#' comprise 9. These will be identified here.
#' 
#' ## Desired outcome
#' For each raw abundance table, species OTU codes must be aligned with short, unique keys, and then species
#' tables must be transposed into samples-species tables. Once these are done, a second set 
#' of files will be produced with sample abundances averaged into sites (site is the replicate in this study),
#' resulting in site-species tables. For the 18S data, transposed species-samples and species-sites tables are needed to produce 
#' UNIFRAC distance matrices later; these tables will be produced and stored in this script.
#' 
#' For each taxonomy table, taxonomic strings must be parsed and unnecessary characters
#' removed. A function is used to streamline the pipeline and to reduce errors. 
#' [Fungal traits](https://link.springer.com/article/10.1007/s13225-020-00466-2) data will be joined
#' with the ITS taxonomy.  
#' 
#' ## Diagnostics
#' 
#' - Individual-based rarefaction is done at the sample level to examine the effect of sequencing depth on OTU recovery: are samples adequate?
#' - Species accumulation is done at the site level to examine the effect of OTU recovery on site community characterization: do samples characterize sites adequately?
#' 
#' # Resources
#' ## Packages and libraries
packages_needed = c("tidyverse", "vegan", "knitr", "colorspace", "plotrix")
packages_installed = packages_needed %in% rownames(installed.packages())
#+ packages,message=FALSE
if (any(!packages_installed)) {
    install.packages(packages_needed[!packages_installed])
}
#+ libraries,message=FALSE
for (i in 1:length(packages_needed)) {
    library(packages_needed[i], character.only = T)
}
#' 
#' ## Functions
#' ### Extract, transform, and load
#+ function_etl
etl <- function(spe, env=sites, taxa, traits=NULL, varname, gene, cluster_type="otu", colname_prefix, folder) {
    
    # Variable definitions
    # spe             = Dataframe or tibble with QIIME2 sequence abundances output, 
    #                   OTUs in rows and samples in columns.
    # env             = Site metadata
    # taxa            = Dataframe or tibble with QIIME2 taxonomy outputs; OTUs in
    #                   rows and metadata in columns. 
    # traits          = Additional dataframe of traits or guilds.
    # varname         = An unique key will be created to replace the cumbersome cluster 
    #                   hash which is produced by QIIME2. Varname, a string, begins the
    #                   column name for this new, short key. Unquoted. Example: otu_num
    # gene            = Gene region targeted in sequencing, e.g.: "ITS", "18S", "16S". 
    #                   Quoted. Used to select() column names, so must match text in 
    #                   column names. Also used to create distinct file names.
    # cluster_type    = Clustering algorithm output, e.g.: "otu", "sv". Quoted.
    #                   Used to create simple cluster IDs. Only OTUs used here.
    # colname_prefix  = Existing, verbose prefix to text of OTU column names. The function removes
    #                   this prefix to make OTU names more concise. 
    # folder          = The function creates output files in the working directory by
    #                   default. To use a subfolder, use this variable. Quoted. 
    #                   Include the "/" before the folder name. If no folder 
    #                   name is desired, use "".
    
    varname <- enquo(varname)
    
    data <- spe %>% left_join(taxa, by = join_by(`#OTU ID`))
    
    # Produce metadata file which is joinable later
    if(gene == "ITS") {
        meta <-
            data %>%
            mutate(!!varname := paste0(cluster_type, "_", row_number())) %>%
            select(!starts_with(gene)) %>%
            rename(otu_ID = `#OTU ID`) %>%
            select(!!varname, everything()) %>%
            separate_wider_delim(
                taxonomy,
                delim = ";",
                names = c(
                    "kingdom",
                    "phylum",
                    "class",
                    "order",
                    "family",
                    "genus",
                    "species"
                ),
                cols_remove = TRUE,
                too_few = "align_start"
            ) %>%
            mutate(kingdom = str_sub(kingdom, 4, nchar(kingdom)),
                   phylum  = str_sub(phylum,  4, nchar(phylum)),
                   class   = str_sub(class,   4, nchar(class)),
                   order   = str_sub(order,   4, nchar(order)),
                   family  = str_sub(family,  4, nchar(family)),
                   genus   = str_sub(genus,   4, nchar(genus)),
                   species = str_sub(species, 4, nchar(species))) %>%
            left_join(traits, by = join_by(phylum, class, order, family, genus)) %>%
            select(-kingdom, -Confidence)
    } else {
        meta <-
            data %>%
            mutate(!!varname := paste0(cluster_type, "_", row_number())) %>%
            select(!starts_with(gene)) %>%
            rename(otu_ID = `#OTU ID`) %>%
            select(!!varname, everything()) %>%
            separate(taxonomy,
                     c("class", "order", "family", "genus", "taxon", "accession"),
                     sep = ";", remove = TRUE, fill = "right") %>%
            select(-Confidence)
    }
    
    # Produce the samples-species table
    spe_samps <-
        data.frame(
            data %>%
                mutate(!!varname := paste0(cluster_type, "_", row_number())) %>%
                select(!!varname, starts_with(gene)),
            row.names = 1
        ) %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        mutate(rowname = str_remove(rowname, colname_prefix)) %>%
        separate_wider_delim(cols = rowname, delim = "_", names = c("field_key", "sample")) %>% 
        mutate(across(c("field_key", "sample"), as.numeric)) %>% 
        left_join(env %>% select(field_key, field_name), by = join_by(field_key)) %>% 
        select(field_name, sample, everything(), -field_key) %>% 
        arrange(field_name, sample)
    # Produce the sites-species table
    spe_avg <-
        spe_samps %>%
        group_by(field_name) %>%
        summarize(across(starts_with(cluster_type), ~ mean(.x)), .groups = "drop") %>%
        select(field_name, everything()) %>% 
        arrange(field_name)
    # Display number of samples in each field to consider survey balance and effort
    samples_fields <-
        spe_samps %>%
        group_by(field_name) %>%
        summarize(n = n(), .groups = "drop") %>%
        left_join(sites, by = join_by(field_name)) %>%
        select(field_name, region, n) %>% 
        arrange(n, field_name) %>% 
        kable(format = "pandoc", caption = paste("Number of samples available in each field", gene, sep = ",\n"))
    # Clean files needed in working directory for use elsewhere
    write_csv(meta, paste0(getwd(), "/", folder, "/spe_", gene, "_metadata.csv"))
    write_csv(spe_samps, paste0(getwd(), "/", folder, "/spe_", gene, "_samples.csv"))
    write_csv(spe_avg, paste0(getwd(), "/", folder, "/spe_", gene, "_avg.csv"))
    # Assign objects to variables for use in the current environment
    out <- list(
        samples_fields      = samples_fields,
        spe_meta            = meta,
        spe_samps           = spe_samps,
        spe_avg             = spe_avg
    )
    return(out)
    
}
#' 
#' ### Species accumulation
#+ function_spe_accum
spe_accum <- function(data) {
    df <- data.frame(
        samples = specaccum(data[, -c(1,2)], method = "exact", conditioned = FALSE)$site,
        richness = specaccum(data[, -c(1,2)], method = "exact", conditioned = FALSE)$richness,
        sd = specaccum(data[, -c(1,2)], method = "exact", conditioned = FALSE)$sd
    )
    return(df)
}
#'
#' ### Confidence interval calculation
#+ function_confint
ci = function(x){std.error(x) * qnorm(0.975)}
#'
#' # Load and process data
#' ## Import data files
#+ species_taxa_traits,message=FALSE
# Species and metadata
its_otu  <- read_delim(file.path(getwd(), "otu_tables/ITS/ITS_otu_raw.txt"), 
                       show_col_types = FALSE)
its_taxa <- read_delim(file.path(getwd(), "otu_tables/ITS/ITS_otu_taxonomy.txt"), 
                       show_col_types = FALSE)
# The 18S OTU file contains an unknown site label in the last column; remove it
amf_otu  <- read_delim(file.path(getwd(), "otu_tables/18S/18S_otu_raw.txt"), 
                       show_col_types = FALSE) %>% select(-last_col())
amf_taxa <- read_delim(file.path(getwd(), "otu_tables/18S/18S_otu_taxonomy.txt"), 
                       show_col_types = FALSE)
traits   <- read_csv(file.path(getwd(),   "otu_tables/2023-02-23_fungal_traits.csv"), 
                     show_col_types = FALSE) %>% select(phylum:primary_lifestyle)
#+ import_sites,message=FALSE
# Site metadata
sites    <- read_csv(file.path(getwd(), "clean_data/sites.csv"), show_col_types = FALSE) %>% 
    mutate(field_type = factor(field_type, ordered = TRUE, levels = c("corn", "restored", "remnant"))) %>% 
    select(-lat, -long, -yr_restore)
#' 
#' ## ETL
#' With the ITS data, two sites are characterized by nine samples, the rest produced ten. 
#+ process_its,message=FALSE
# ITS
(its <-
    etl(
        spe = its_otu,
        taxa = its_taxa,
        traits = traits,
        varname = "otu_num",
        gene = "ITS",
        colname_prefix = "ITS_TGP_",
        folder = "clean_data"
    ))
#' With the 18S data, four sites are characterized by nine samples, the rest produced ten.
#+ process_18S,message=FALSE,warning=FALSE
# 18S
(amf <-
    etl(
        spe = amf_otu,
        taxa = amf_taxa,
        varname = "otu_num",
        gene = "18S",
        colname_prefix = "X18S_TGP_",
        folder = "clean_data"
    ))
#' 
#' ## Post-processing 18S data
#' The raw 18S OTU abundance table must be transposed into species-samples or species-sites
#' tables so that UNIFRAC distances my be computed later. 
#+ wrangle_amf_samps
# Species-samples table of 18S OTUs
data.frame(
    amf$spe_samps %>%
        mutate(field_sample = paste(field_name, sample, sep = "_")) %>% 
        select(field_sample, everything(), -field_name, -sample),
    row.names = 1
) %>%
t() %>%
as.data.frame() %>%
rownames_to_column(var = "otu_num") %>%
left_join(amf$spe_meta %>% select(otu_num, otu_ID), by = join_by(otu_num)) %>%
select(otu_ID, everything(), -otu_num) %>% 
write_tsv(file.path(getwd(), "otu_tables/18S/18S_samps_4unifrac.tsv"))
#+ wrangle_amf_avg
# Species-sites table of 18S OTUs
data.frame(
    amf$spe_avg %>% select(field_name, everything()),
    row.names = 1
) %>%
t() %>%
as.data.frame() %>%
rownames_to_column(var = "otu_num") %>%
left_join(amf$spe_meta %>% select(otu_num, otu_ID), by = join_by(otu_num)) %>%
select(otu_ID, everything(), -otu_num) %>% 
write_tsv(file.path(getwd(), "otu_tables/18S/18S_avg_4unifrac.tsv"))
#' 
#' # Sampling depth and coverage
#' 1. Did sequencing exhaust OTU accumulation at varying sequencing depths?
#' 1. Were sites adequately characterized by samples?
#' 
#' ## Rarefaction
#' ### ITS: sample based
its_rc <- rarecurve(its$spe_samps %>% 
                        mutate(field_sample = paste(field_name, sample, sep = "_")) %>% 
                        column_to_rownames(var = "field_sample") %>% 
                        select(-field_name, -sample), 
                    step = 1,
                    tidy = TRUE)
#+ its_rarefaction,warning=FALSE,message=FALSE,fig.width=6,fig.height=8,fig.align='center'
its_rc %>% 
    separate_wider_delim(Site, delim = "_", names = c("field_name", "sample_key"), cols_remove = FALSE) %>% 
    rename(seq_abund = Sample, otus = Species, field_sample = Site) %>% 
    left_join(sites, by = join_by(field_name)) %>% 
    ggplot(aes(x = seq_abund, y = otus, group = field_sample)) +
    facet_wrap(vars(field_type), ncol = 1, scales = "free") +
    geom_line(aes(color = field_type), linewidth = 0.4) +
    scale_color_discrete_qualitative(palette = "Harmonic") +
    labs(x = "Number of individuals (sequence abundance)",
         y = "OTUs",
         title = "Rarefaction of ITS samples") +
    theme_bw()
#' 
#' ### ITS: site based
#' Sequences summed in sites
its_rc_site <- rarecurve(
    its$spe_samps %>% 
        group_by(field_name) %>% 
        summarize(across(starts_with("otu"), sum)) %>% 
        column_to_rownames(var = "field_name"),
    step = 1,
    tidy = TRUE)
#+ its_rarefaction,warning=FALSE,message=FALSE,fig.width=6,fig.height=8,fig.align='center'
its_rc_site %>% 
    rename(seq_abund = Sample, otus = Species, field_name = Site) %>% 
    left_join(sites, by = join_by(field_name)) %>% 
    ggplot(aes(x = seq_abund, y = otus, group = field_name)) +
    facet_wrap(vars(field_type), ncol = 1, scales = "free_y") +
    geom_line(aes(color = field_type), linewidth = 0.4) +
    scale_color_discrete_qualitative(palette = "Harmonic") +
    labs(x = "Number of individuals (sequence abundance)",
         y = "OTUs",
         title = "Rarefaction of ITS samples") +
    theme_bw()
#' 
#' ### AMF: sample based
amf_rc <- rarecurve(amf$spe_samps %>% 
                        mutate(field_sample = paste(field_name, sample, sep = "_")) %>% 
                        column_to_rownames(var = "field_sample") %>% 
                        select(-field_name, -sample), 
                    step = 1,
                    tidy = TRUE)
#+ amf_rarefaction,warning=FALSE,message=FALSE,fig.width=6,fig.height=8,fig.align='center'
amf_rc %>% 
    separate_wider_delim(Site, delim = "_", names = c("field_name", "sample_key"), cols_remove = FALSE) %>% 
    rename(seq_abund = Sample, otus = Species, field_sample = Site) %>% 
    left_join(sites, by = join_by(field_name)) %>% 
    ggplot(aes(x = seq_abund, y = otus, group = field_sample)) +
    facet_wrap(vars(field_type), ncol = 1, scales = "free") +
    geom_line(aes(color = field_type), linewidth = 0.4) +
    scale_color_discrete_qualitative(palette = "Harmonic") +
    labs(x = "Number of individuals (sequence abundance)",
         y = "OTUs",
         title = "Rarefaction of 18S samples") +
    theme_bw()
#' 
#' ### AMF: site based
#' Sequences summed in sites
amf_rc_site <- rarecurve(
    amf$spe_samps %>% 
        group_by(field_name) %>% 
        summarize(across(starts_with("otu"), sum)) %>% 
        column_to_rownames(var = "field_name"),
    step = 1,
    tidy = TRUE)
#+ its_rarefaction,warning=FALSE,message=FALSE,fig.width=6,fig.height=8,fig.align='center'
amf_rc_site %>% 
    rename(seq_abund = Sample, otus = Species, field_name = Site) %>% 
    left_join(sites, by = join_by(field_name)) %>% 
    ggplot(aes(x = seq_abund, y = otus, group = field_name)) +
    facet_wrap(vars(field_type), ncol = 1, scales = "free_y") +
    geom_line(aes(color = field_type), linewidth = 0.4) +
    scale_color_discrete_qualitative(palette = "Harmonic") +
    labs(x = "Number of individuals (sequence abundance)",
         y = "OTUs",
         title = "Rarefaction of 18S samples") +
    theme_bw()
#' 
#' Sequence abundance varied widely among samples and sites, but flat curves for each suggest that OTU recovery was 
#' exhausted for each of them.
#' 
#' ## Species accumulation
#' #+ species_accumulation,warning=FALSE,message=FALSE
accum <- bind_rows(
    list(
        ITS = bind_rows(
            split(its$spe_samps, ~ field_name) %>% 
                map(spe_accum),
            .id = "field_name"
        ),
        AMF = bind_rows(
            split(amf$spe_samps, ~ field_name) %>% 
                map(spe_accum),
            .id = "field_name"
        )
    ),
    .id = "dataset"
) %>% 
    mutate(dataset = factor(dataset, ordered = TRUE, levels = c("ITS", "AMF"))) %>% 
    left_join(sites, by = join_by(field_name))
#+ species_accumulation_fig,warning=FALSE,message=FALSE,fig.width=8,fig.height=5,fig.align='center'
ggplot(accum, aes(x = samples, y = richness, group = field_name)) +
    facet_wrap(vars(dataset), scales = "free_y") +
    geom_line(aes(color = field_type)) +
    geom_segment(aes(x = samples, y = richness-sd, xend = samples, yend = richness+sd, color = field_type)) +
    scale_color_discrete_qualitative(palette = "Harmonic") +
    labs(x = "Samples", y = expression(N[0]), 
         caption = "Species accumulation by the \"exact\" method with unconditioned, moment-based standard deviation") +
    scale_x_continuous(breaks = c(0,2,4,6,8,10)) +
    theme_bw()
#' Steep curves suggest that sampling was inadequate to recover all richness in sites. With ITS data, 
#' OTU recovery curves were similar across sites. For AMF, sites vary in curve steepness, suggesting that 
#' some were better characterized than others. 

