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
#' 
#' 
#' 
#' # Packages and libraries
packages_needed = c(
    "conflicted",
    "colorspace",
    "emmeans",
    "gridExtra",
    "knitr",
    "tidyverse",
    "vegan"
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
#+ conflicts,message=FALSE
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("diversity", "vegan")
#' 
#+ graphics_styles
source(file.path(getwd(), "resources", "styles.txt"))
#' 
#' # Data
#' ## Sites-species tables
#' List *spe* holds average sequence abundances for the top 6 samples per field. 
#' CSV files were produced in `process_data.R`
spe <- list(
    its_avg   = read_csv(paste0(getwd(), "/clean_data/spe_ITS_avg.csv"), 
                         show_col_types = FALSE),
    amf_avg   = read_csv(paste0(getwd(), "/clean_data/spe_18S_avg.csv"), 
                         show_col_types = FALSE)
)
#' 
#' ## Site metadata and design
sites <- read_csv(paste0(getwd(), "/clean_data/sites.csv"), show_col_types = FALSE) %>% 
    mutate(field_type = factor(field_type, ordered = TRUE, levels = c("corn", "restored", "remnant"))) %>% 
    select(-lat, -long, -yr_restore)
#' 
#' ## PLFA
#' Use only 18.2 for soil fungi
fa <- read_csv(file.path(getwd(), "clean_data/plfa.csv"), show_col_types = FALSE) %>% 
    rename(fungi_18.2 = fa_18.2) %>% 
    select(field_name, fungi_18.2, starts_with("gram"), bacteria, actinomycetes, amf)
#' 
#' 
#' # Functions
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
#' 