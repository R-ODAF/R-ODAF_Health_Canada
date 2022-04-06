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


##############################################################################################
# SETUP
##############################################################################################

config <- yaml::read_yaml(here::here("config","config.yaml"), eval.expr = T)
params <- config$deseq

projectdir <- params$projectdir
if(is.null(projectdir)){
  projectdir <- here::here()
  params$projectdir <- projectdir
}

paths <- set_up_paths(params)
get_analysis_id <- get_analysis_id(params)

species_data <- load_species(params$species)
ensembl <- useMart("ensembl",
                   dataset = species_data$ensembl_species,
                   host = "useast.ensembl.org")

# set some additional parameters based on platform
if (params$platform == "RNA-Seq") {
  SampleDataFile <- file.path(paths$processed, "genes.data.tsv")
  params$sampledata_sep = "\t"
  
  params$threshold <- 1000000 # Number of aligned reads per sample required
  params$MinCount <- 1
  params$alpha <- pAdjValue <- 0.05 # Relaxed from 0.01
  params$linear_fc_filter <- 1.5
  params$biomart_filter <- "ensembl_gene_id"
} else if (params$platform == "TempO-Seq") {
  SampleDataFile <- file.path(paths$processed, "count_table.tsv")
  params$sampledata_sep = "\t"

  params$threshold = 100000 # Number of aligned reads per sample required
  params$MinCount <- 0.5
  params$alpha <- pAdjValue <- 0.05 
  params$linear_fc_filter <- 1.5
  
  bs <- load_biospyder(params$biospyder_dbs, species_data$temposeq_manifest)
  params$biospyder_ID <- bs$biospyder_ID
  params$biomart_filter <- bs$biomart_filter
  params$biospyder_filter <- bs$biospyder_filter
  params$biospyder <- bs$biospyder
} else { 
  stop("Platform/technology not recognized") 
}

# Set this variable to be TRUE if you want to have separate plots of top genes as defined in the R-ODAF template
params$R_ODAF_plots = FALSE


##############################################################################################
# DATA LOADING AND PROCESSING
##############################################################################################

# Identify where metadata can be found
MetadataFile <- file.path(paths$metadata, "metadata.QC_applied.txt")
ContrastsFile <- file.path(paths$metadata, "contrasts.txt")

# Read in metadata
DESeqDesign <- read.delim(MetadataFile,
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
design_to_use <- params$design

# if multiple intgroups, combine into new group variable
if (length(intgroup) > 1){
    new_group_name <- paste(intgroup,collapse="_")
    if(new_group_name %in% colnames(DESeqDesign)){
        stop(paste0("Metadata cannot contain the column name ",new_group_name))
    }
    DESeqDesign <- DESeqDesign %>% unite((!!sym(new_group_name)), intgroup, remove = FALSE)
    # now redo the contrasts
    contrasts <- contrasts %>%
        left_join(DESeqDesign, by=c("V1"=params$design)) %>%
        dplyr::select(V1, V1.new = new_group_name, V2) %>%
        unique() %>%
        left_join(DESeqDesign, by=c("V2"=params$design)) %>%
        dplyr::select(V1, V1.new, V2, V2.new=new_group_name) %>%
        unique() %>%
        dplyr::select(V1=V1.new, V2=V2.new)
    design_to_use <- new_group_name
    intgroup <- new_group_name
}

# load count data
sampleData <- load_count_data(SampleDataFile, params$sampledata_sep)
print(typeof(sampleData))

processed <- process_data_and_metadata(sampledata, DESeqDesign, contrasts, intgroup, design_to_use, params)
sampleData <- processed$sampleData
DESeqDesign <- processed$DESeqDesign
contrasts <- processed$contrasts

# set up facets if necessary
# the facets array will be all facets if group_filter is not set, and the filter otherwise
if(!is.na(params$group_facet)){
    if(!is.na(params$group_filter)){
        facets <- params$group_filter
    }else {
        facets <- DESeqDesign %>%
            filter(!(params$group_facet) %in% c(params$exclude_groups)) %>%
            filter(!solvent_control) %>%
            pull(params$group_facet) %>% 
            unique()
    }
} else {
    facets <- NA
}

stopifnot((is.na(params$group_facet) || length(facets) > 0))


ddsList <- list()
designList <- list()
overallResList <- list()
rldList <- list()

if(is.na(params$group_facet)){
    message("### Learning a single model for the whole experiment. ###")
    dds <- learn_deseq_model(sampleData, DESeqDesign, intgroup, design_to_use, params)
    # TODO: do this. need nuisance params
    # rld <- regularize_data(dds, covariates, nuisance)
    DESeq_results <- get_DESeq_results(dds, DESeqDesign, contrasts, design_to_use, params, NA, paths$DEG_output)
    resList <- DESeq_results$resList
    ddsList[['all']] <- dds
    overallResList[['all']] <- DESeq_results$resList
    designList[['all']] <- DESeqDesign
    # rldList[['all']] <- rld
} else {
    for (current_filter in facets) {
        message(paste0("### Learning model for ", current_filter, ". ###"))
        metadata_subset <- subset_metadata(DESeqDesign, design_to_use, contrasts, current_filter)
        DESeqDesign_subset <- metadata_subset$DESeqDesign
        contrasts_subset <- metadata_subset$contrasts
        sampleData_subset <- subset_data(sampleData, DESeqDesign_subset)

        check_data(sampleData_subset, DESeqDesign_subset, contrasts_subset)

        ddsList[[current_filter]] <- learn_deseq_model(sampleData_subset, DESeqDesign_subset, intgroup, design_to_use, params)
        designList[[current_filter]] <- DESeqDesign_subset
        # TODO: do this. need nuisance params
        # rldList[[current_filter]] <- regularize_data(dds, covariates, nuisance)
        DESeq_results <- get_DESeq_results(ddsList[[current_filter]], designList[[current_filter]], contrasts_subset, design_to_use, params, current_filter, paths$DEG_output)
        overallResList[[current_filter]] <- DESeq_results$resList
    }
}


summary_counts <- data.frame()
if(is.na(params$group_facet)){
    comparisons <- names(resList)
    for(comp in comparisons){ # by comparison
        res <- resList[[comp]]
        counts <- nrow(res)
        row <- data.frame(comparison=comp, DEG=counts)
        summary_counts <- rbind(summary_counts, row)
    }
} else {
    for (current_filter in facets) {
        resList <- overallResList[[current_filter]]
        comparisons <- names(resList)
        for(comp in comparisons){ # by comparison
            res <- resList[[comp]]
            counts <- nrow(res)
            row <- data.frame(facet=current_filter, comparison=comp, DEG=counts)
            summary_counts <- rbind(summary_counts, row)
        }
    }
}

message(paste0(sum(summary_counts$DEG), " total DEG counts found. Missing rows indicate 0 DEGs passed filters"))
message(paste(capture.output(summary_counts), collapse="\n"))

# write the table of DEG counts to a file
write.table(summary_counts,
            file = file.path(paths$DEG_output, paste0(params$project_name, "_DEG_counts_summary.txt")),
            sep = "\t",
            quote = FALSE)

# save DESeq results to a file

save(ddsList, designList, overallResList, rldList, DESeqDesign, facets, contrasts, intgroup, design_to_use, paths, file=file.path(paths$DEG_output, paste0(params$project_name, "_DEG_data.RData")))