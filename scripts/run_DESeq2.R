#!/usr/bin/R
# Custom parameters for the report
suppressMessages(library('tidyverse'))
suppressMessages(library('yaml'))
suppressMessages(library('DESeq2'))
suppressMessages(library('gtools'))
suppressMessages(library('biomaRt'))
suppressMessages(library('BiocParallel'))

# assume this is being run from within the R project
projectdir <- here::here()

source(here::here("scripts","setup_functions.R"))
source(here::here("scripts","data_functions.R"))
source(here::here("scripts","file_functions.R"))
source(here::here("scripts","DESeq_functions.R"))

config <- yaml::read_yaml(here::here("config","config.yaml"), eval.expr = T)
params <- config$params

paths <- set_up_paths(params)
get_analysis_id <- get_analysis_id(params)

# Identify where metadata can be found
SampleKeyFile <- file.path(paths$metadata, "metadata.QC_applied.txt")
ContrastsFile <- file.path(paths$metadata, "contrasts.txt")

# Read in metadata
DESeqDesign <- read.delim(SampleKeyFile,
                          stringsAsFactors = FALSE,
                          sep = "\t",
                          header = TRUE,
                          quote = "\"",
                          row.names = 1) # Column must have unique IDs!!
DESeqDesign$original_names <- rownames(DESeqDesign)
DESeqDesignAsRead <- DESeqDesign

# read in contrasts
contrasts <- read.delim(ContrastsFile, stringsAsFactors = FALSE, sep = "\t", header = FALSE,  quote = "\"")

# set interesting groups
intgroup <- params$intgroup # "Interesting groups" - experimental group/covariates
# TODO: maybe here's where the factor conversion should happen?
# TODO: what about nuisance variables?
# TODO: faceting might make certain experimental designs illegal. Should we check for those or let DESeq do it?


species_data <- load_species(params$species)
ensembl <- useMart("ensembl",
                   dataset = species_data$ensembl_species,
                   host = "useast.ensembl.org")


if (params$platform == "RNA-Seq") {
  SampleDataFile <- file.path(paths$processed, "genes.data.tsv")
  params$sampledata_sep = "\t"
  
  params$threshold <- 1000000 # Number of aligned reads per sample required
  params$MinCount <- 1
  params$alpha <- pAdjValue <- 0.05 # Relaxed from 0.01
  params$linear_fc_filter <- 1.5
  params$biomart_filter <- "ensembl_gene_id"
} else if (params$platform == "TempO-Seq") {
  SampleDataFile <- file.path(paths$processed, "count_table.csv")
  params$sampledata_sep = ","

  params$threshold = 100000 # Number of aligned reads per sample required
  params$MinCount <- 0.5
  params$alpha <- pAdjValue <- 0.05 
  params$linear_fc_filter <- 1.5
  
  bs <- load_biospyder(biospyder_dbs, temposeq_manifest)
  params$biospyder_ID <- bs$biospyder_ID
  params$biomart_filter <- bs$biomart_filter
  params$biospyder_filter <- bs$biospyder_filter
  params$biospyder <- bs$biospyder
} else { 
  stop("Platform/technology not recognized") 
}

# Se this variable to be TRUE if you want to have separate plots of top genes as defined in the R-ODAF template
params$R_ODAF_plots = FALSE


#################9
# TODO: set aside saving/loading data for now
# it might be less important anyway, since the script will save the DESeq output
# just focus on the DEseq calculation, including faceting and saving the output
# maybe if no facets, can format the data the same with facet name "all" or something.
#################
# if(params$use_cached_RData){
#     if(is.na(params$group_facet)){
#         dds <- load_cached_data(paths$RData, sampleData, params)
#     } else {
#         ddsList <- load_cached_data(paths$RData, sampleData, params)
#     }
# } else {
#     sampleData <- load_count_data(SampleDataFile, params$sampledata_sep)
#     processed <- process_data(sampledata, DESeqDesign, intgroup, params)
#     sampleData <- processed$sampleData
#     DESeqDesign <- processed$DESeqDesign

#     if(is.na(params$group_facet)){
#         dds <- learn_deseq_model(sampledata, DESeqDesign, intgroup, params)
#         save_cached_data(dds, paths$RData, params)
#     } else {
#         for (current_filter in facets) {
            
#         }
#     }

# }

sampleData <- load_count_data(SampleDataFile, params$sampledata_sep)

processed <- process_data_and_metadata(sampledata, DESeqDesign, contrasts, intgroup, params)
sampleData <- processed$sampleData
DESeqDesign <- processed$DESeqDesign
contrasts <- processed$contrasts

# set up facets if necessary
# facets will be all facets if group_filter is not set, and the filter otherwise
if(!is.na(params$group_facet)){
    if(!is.na(params$group_filter)){
        facets <- params$group_filter
    }else {
        facets <- DESeqDesign %>%
            filter(!(params$group_facet) %in% c(params$exclude_groups, skip_extra)) %>%
            filter(!solvent_control) %>%
            pull(params$group_facet) %>% 
            unique()
    }
}

if(is.na(params$group_facet)){
    dds <- learn_deseq_model(sampledata, DESeqDesign, intgroup, params)
    #save_cached_data(dds, paths$RData, params)
} else {
    ddsList <- list()
    metadata_subset <- subset_metadata(DESeqDesign, params, contrasts)
    DESeqDesign_subset <- metadata_subset$DESeqDesign
    contrasts_subset <- metadata_subset$contrasts
    sampleData_subset <- subset_data(sampleData, DESeqDesign_subset)

    for (current_filter in facets) {
        ddsList[[current_filter]] <- learn_deseq_model(sampleData_subset, DESeqDesign_subset, intgroup, params)
    }
}



if (is.na(params$group_facet)) { # all data in one facet
    message("Learning a single model for the whole experiment.")
    res <- get_DESeq_results(dds, contrasts, params, NA, paths$DEG_output)
} else {
    resList <- list()
    if (any(!is.na(params$group_filter))) { # filter facets
        message(paste0("The group(s) of interest is (are) ",
                 paste(params$group_filter, collapse = " and "),".\n",
                 "Learning a single model for that (those) groups."))
    } else { # do all facets separately
        message(paste0("Learning multiple models based on ", params$group_facet, "..."))
    }
    for (current_filter in facets) {
        res <- get_DESeq_results(ddsList[[current_filter]], contrasts, params, current_filter, paths$DEG_output)
        resList[[current_filter]] <- res
    }
}

