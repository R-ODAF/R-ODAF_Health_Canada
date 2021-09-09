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
print(projectdir)

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

DESeqDesign <- filter_metadata(DESeqDesign, params)
if(!is.na(params$sortcol)){
  sorted <- sort_metadata(DESeqDesign, contrasts, params)
  DESeqDesign <- sorted$design
  contrasts <- sorted$contrasts
}

species_data <- load_species(params$species)


# Identify where count data can be found
if (params$platform == "TempO-Seq") {
  SampleDataFile <- file.path(paths$processed, "count_table.csv")
  sampledata_sep = ","
} else {
  SampleDataFile <- file.path(paths$processed, "genes.data.tsv")
  sampledata_sep = "\t"
}

ensembl <- useMart("ensembl",
                   dataset = species_data$ensembl_species,
                   host = "useast.ensembl.org")


if (params$platform == "RNA-Seq") {
  params$threshold <- 1000000 # Number of aligned reads per sample required
  params$MinCount <- 1
  params$alpha <- pAdjValue <- 0.05 # Relaxed from 0.01
  params$linear_fc_filter <- 1.5
  params$biomart_filter <- "ensembl_gene_id"
} else if (params$platform == "TempO-Seq") {
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

if (is.na(params$group_facet)) { # all data in one facet
  message("Writing a single report for whole experiment.")
  if(params$use_cached_RData){
    load_cached_data(paths$RData, params)
  } else {
    sampleData <- load_count_data(SampleDataFile, sampledata_sep)
    processed <- process_data(sampledata, DESeqDesign, intgroup, params)
    sampleData <- processed$sampleData
    DESeqDesign <- processed$DESeqDesign

    dds <- learn_deseq_model(sampledata, DESeqDesign, intgroup, params)

    save_cached_data(dds, paths$RData, params)
  }
} else {
  stop("not implemented yet")
    facets <- DESeqDesign %>%
        filter(!(params$group_facet) %in% c(params$exclude_groups, skip_extra)) %>%
        filter(!solvent_control) %>%
        pull(params$group_facet) %>% 
        unique()

    if (any(!is.na(params$group_filter))) { # filter facets
        message(paste0("The group(s) of interest is (are) ",
                 paste(params$group_filter, collapse = " and "),".\n",
                 "Writing a single report for that (those) groups."))

    } else { # do all facets separately
      message(paste0("Making multiple reports based on ", params$group_facet, "..."))

    }
}