#!/usr/bin/R
# Custom parameters for the report
suppressMessages(library('tidyverse'))
suppressMessages(library('yaml'))
suppressMessages(library('DESeq2'))
suppressMessages(library('gtools'))
suppressMessages(library('biomaRt'))
suppressMessages(library('BiocParallel'))

source(here::here("scripts","setup_functions.R"))
source(here::here("scripts","data_functions.R"))
source(here::here("scripts","file_functions.R"))
source(here::here("scripts","DESeq_functions.R"))


##############################################################################################
# SETUP
##############################################################################################

config <- yaml::read_yaml(here::here("config","config.yaml"), eval.expr = T)

# Combine required params from config
params <- c(config$common, config$DESeq2)
# replace nulls in params with NA
params <- replace_nulls_in_config(params)
# If projectdir is not set, figure out current project root directory
projectdir <- params$projectdir
if (is.na(projectdir)) {
  projectdir <- here::here()
  params$projectdir <- projectdir
}




paths <- set_up_paths(params)

get_analysis_id <- get_analysis_id(params)

species_data <- load_species(params$species, params$wikipathways_filename, params$biospyder_manifest_file)
params$species_data <- species_data
# ensembl <- useMart("ensembl",
#                    dataset = species_data$ensembl_species,
#                    host = "useast.ensembl.org")

params <- set_up_platform_params(params)

# Set this variable to be TRUE if you want to have separate plots of top genes as defined in the R-ODAF template
params$R_ODAF_plots <- FALSE

check_required_params(params)

##############################################################################################
# DATA LOADING AND PROCESSING
##############################################################################################

# Identify where metadata can be found
MetadataFile <- file.path(paths$metadata, "metadata.QC_applied.txt")
ContrastsFile <- file.path(paths$metadata, "contrasts.txt")

sample_count_metadata <- list()

# Read in metadata
exp_metadata <- read.delim(MetadataFile,
                          stringsAsFactors = FALSE,
                          sep = "\t",
                          header = TRUE,
                          quote = "\"",
                          colClasses = "character",
                          row.names = 1) # Column must have unique IDs!!
exp_metadata$original_names <- rownames(exp_metadata)
exp_metadataAsRead <- exp_metadata

sample_count_metadata$samples_postQC <- nrow(exp_metadata)

# read in contrasts
contrasts <- read.delim(ContrastsFile, stringsAsFactors = FALSE, sep = "\t", header = FALSE,  quote = "\"", 
                        colClasses = "character")
# set interesting groups
intgroup <- params$intgroup # "Interesting groups" - experimental group/covariates
design_to_use <- params$design

# if multiple intgroups, combine into new group variable
if (length(intgroup) > 1){
    new_group_name <- paste(intgroup,collapse="_")
    if(new_group_name %in% colnames(exp_metadata)){
        stop(paste0("Metadata cannot contain the column name ",new_group_name))
    }
    exp_metadata <- exp_metadata %>% unite((!!sym(new_group_name)), intgroup, remove = FALSE)
    # now redo the contrasts
      # contrasts <- contrasts %>%
      #     left_join(exp_metadata, by=c("V1"=params$design)) %>%
      #     dplyr::select(V1, V1.new = new_group_name, V2) %>%
      #     unique() %>%
      #     left_join(exp_metadata, by=c("V2"=params$design)) %>%
      #     dplyr::select(V1, V1.new, V2, V2.new=new_group_name) %>%
      #     unique() %>%
      #     dplyr::select(V1=V1.new, V2=V2.new)
    design_to_use <- new_group_name
    covariates <- intgroup[intgroup != params$design]
    intgroup <- new_group_name
} else {
  covariates <- NA
}
original_design <- params$design

# load count data
count_data <- load_count_data(params$count_data_file, params$sampledata_sep)


processed <- process_data_and_metadata(count_data, exp_metadata, contrasts, intgroup, design_to_use, params)
count_data <- processed$count_data
exp_metadata <- processed$exp_metadata
contrasts <- processed$contrasts 

sample_count_metadata$samples_filtered <- nrow(exp_metadata)

# set up facets if necessary
# the facets array will be all facets if group_filter is not set, and the filter otherwise
if(!is.na(params$group_facet)){
  if(!is.na(params$group_filter)){
    facets <- params$group_filter
  }else {
    facets <- exp_metadata %>%
      filter(!(params$group_facet) %in% c(params$exclude_groups)) %>%
      filter(!(solvent_control==TRUE)) %>%
      pull(params$group_facet) %>% 
      unique()
  }
} else {
  facets <- NA
}

stopifnot((is.na(params$group_facet) || length(facets) > 0))


# set up the rest of the output paths (requires facets)
paths <- set_up_paths_2(paths,params,facets)

ddsList <- list()
designList <- list()
contrastsList <- list()
overallResListAll <- list()
overallResListFiltered <- list()
overallResListDEGs <- list()
rldList <- list()
mergedDEGsList <- list()
filtered_table <- data.frame()

if(is.na(params$group_facet)){
    message("### Learning a single model for the whole experiment. ###")
    if(params$write_additional_output){
      write_additional_output(count_data, exp_metadata, design_to_use, params)
    }
    dds <- learn_deseq_model(count_data, exp_metadata, design_to_use, params)
    rld <- regularize_data(dds, original_design, covariates, params$batch_var)
    DESeq_results <- get_DESeq_results(dds, exp_metadata, contrasts, design_to_use, params, NA, paths$DEG_output)
    ddsList[['all']] <- dds
    overallAllGenes <- DESeq_results$dfGenes
    overallResListAll[['all']] <- DESeq_results$resListAll
    overallResListFiltered[['all']] <- DESeq_results$resListFiltered
    overallResListDEGs[['all']] <- DESeq_results$resListDEGs
    designList[['all']] <- exp_metadata
    contrastsList[['all']] <- contrasts
    rldList[['all']] <- rld
    mergedDEGsList[['all']] <- DESeq_results$mergedDEGs
    filtered_table <- rbind(filtered_table, DESeq_results$filtered_table)
} else {
    for (current_filter in facets) {
        message(paste0("### Learning model for ", current_filter, ". ###"))
        metadata_subset <- subset_metadata(exp_metadata, design_to_use, contrasts, params$group_facet, current_filter)
        exp_metadata_subset <- metadata_subset$exp_metadata
        contrasts_subset <- metadata_subset$contrasts
        count_data_subset <- subset_data(count_data, exp_metadata_subset)
        
        check_data(count_data_subset, exp_metadata_subset, contrasts_subset)
        
        if(params$write_additional_output){
          write_additional_output(count_data_subset, exp_metadata_subset, design_to_use, params)
        }
        ddsList[[current_filter]] <- learn_deseq_model(count_data_subset, exp_metadata_subset, design_to_use, params)
        designList[[current_filter]] <- exp_metadata_subset
        contrastsList[[current_filter]] <- contrasts_subset
        rldList[[current_filter]] <- regularize_data(ddsList[[current_filter]], original_design, covariates, params$batch_var)
        DESeq_results <- get_DESeq_results(ddsList[[current_filter]], designList[[current_filter]], contrasts_subset, design_to_use, params, current_filter, paths$DEG_output)
        overallResListAll[[current_filter]] <- DESeq_results$resListAll
        overallResListFiltered[[current_filter]] <- DESeq_results$resListFiltered
        overallResListDEGs[[current_filter]] <- DESeq_results$resListDEGs
        mergedDEGsList[[current_filter]] <- DESeq_results$mergedDEGs
        filtered_table <- rbind(filtered_table, DESeq_results$filtered_table)
    }
}

summary_counts <- data.frame()
if(is.na(params$group_facet)){
    resList <- overallResListDEGs[['all']]
    comparisons <- names(resList)
    for(comp in comparisons){ # by comparison
        res <- resList[[comp]]
        counts <- nrow(res)
        row <- data.frame(comparison=comp, DEG=counts)
        summary_counts <- rbind(summary_counts, row)
    }
} else {
    for (current_filter in facets) {
        resList <- overallResListDEGs[[current_filter]]
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

# save DESeq results to a file
save(ddsList, designList, contrastsList, overallResListAll, overallResListFiltered, overallResListDEGs, rldList, mergedDEGsList, exp_metadata, facets, count_data, contrasts, intgroup, design_to_use, paths, filtered_table, sample_count_metadata, file=file.path(paths$RData, paste0(params$project_title, "_DEG_data.RData")))

source(here::here("scripts","write_output_tables.R"))
