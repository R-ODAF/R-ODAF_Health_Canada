#!/usr/bin/R
# Custom parameters for the report
suppressMessages(library(dplyr))
suppressMessages(library(DESeq2))
suppressMessages(library(BiocParallel))
suppressMessages(library(R.ODAF.utils))

# Parse command line arguments
run_args <- commandArgs(trailingOnly = TRUE)

# Check if at least one argument is provided
if (length(run_args) > 0) {
  # Assume the first argument is the new location
  results_location_arg <- run_args[1]
  # Source functions and pass along analysis directory argument
} else {
  message("Error: Missing argument. Provide the analysis directory name as an argument.")
}

# Load project parameters
params <- R.ODAF.utils::get_params(context = "analysis")

# Set up initial project paths
paths <- R.ODAF.utils::set_up_filepaths(params, results_location_arg)

##############################################################################################
# DATA LOADING AND PROCESSING
##############################################################################################

# Read in metadata
exp_metadata <- R.ODAF.utils::get_metadata(file.path(paths$metadata, "metadata.QC_applied.txt"), paths)
sample_count_metadata <- list()
sample_count_metadata$samples_postQC <- nrow(exp_metadata)

# read in contrasts
exp_contrasts <- R.ODAF.utils::get_contrasts(file.path(paths$contrasts, "contrasts.txt"), paths)
original_design <- params$design

# load count data
count_data <- load_count_data(file.path(paths$processed, "count_table.tsv"), params$sampledata_sep)
# TODO - let user define count table file path in config?

processed <- process_data_and_metadata(count_data, exp_metadata, exp_contrasts, params)
count_data <- processed$count_data
exp_metadata <- processed$exp_metadata
exp_contrasts <- processed$contrasts

sample_count_metadata$samples_filtered <- nrow(exp_metadata)

# Set up facets if necessary
# the facets array will be all facets if deseq_filter is not set, and the filter otherwise
facets <- get_facets(exp_metadata, params)
stopifnot((is.na(params$deseq_facet) || length(facets) > 0))

# set up the rest of the output paths (requires facets)
paths <- R.ODAF.utils::set_up_filepaths(params,
                                        results_location_arg,
                                        exp_metadata,
                                        make_deseq2_dirs = TRUE)

ddsList <- list()
designList <- list()
contrastsList <- list()
overallResListAll <- list()
overallResListFiltered <- list()
overallResListDEGs <- list()
rldList <- list()
mergedDEGsList <- list()
filtered_table <- data.frame()

if (is.na(params$deseq_facet)){
  message("### Learning a single model for the whole experiment. ###")
  if(params$write_additional_output){
    write_additional_output(count_data, exp_metadata, params[["design"]], params)
  }
  dds <- learn_deseq_model(count_data, exp_metadata, params[["design"]], params)
  rld <- regularize_data(dds, original_design, covariates, params$batch_var)
  DESeq_results <- get_DESeq_results(dds, exp_metadata, exp_contrasts, params[["design"]], params, NA, paths$DEG_output)
  ddsList[['all']] <- dds
  overallAllGenes <- DESeq_results$dfGenes
  overallResListAll[['all']] <- DESeq_results$resListAll
  overallResListFiltered[['all']] <- DESeq_results$resListFiltered
  overallResListDEGs[['all']] <- DESeq_results$resListDEGs
  designList[['all']] <- exp_metadata
  contrastsList[['all']] <- exp_contrasts
  rldList[['all']] <- rld
  mergedDEGsList[['all']] <- DESeq_results$mergedDEGs
  filtered_table <- rbind(filtered_table, DESeq_results$filtered_table)
} else {
  for (current_filter in facets) {
    message(paste0("### Learning model for ", current_filter, ". ###"))
    metadata_subset <- subset_metadata(exp_metadata, params[["design"]], exp_contrasts, params$deseq_facet, current_filter)
    exp_metadata_subset <- metadata_subset$exp_metadata
    contrasts_subset <- metadata_subset$contrasts
    count_data_subset <- subset_data(count_data, exp_metadata_subset)

    check_data(count_data_subset, exp_metadata_subset, contrasts_subset)

    if (params$write_additional_output) {
      write_additional_output(count_data_subset, exp_metadata_subset, params[["design"]], params)
    }
    ddsList[[current_filter]] <- learn_deseq_model(count_data_subset, exp_metadata_subset, params[["design"]], params)
    designList[[current_filter]] <- exp_metadata_subset
    contrastsList[[current_filter]] <- contrasts_subset
    rldList[[current_filter]] <- regularize_data(ddsList[[current_filter]], original_design, covariates = NA, params$batch_var)
    DESeq_results <- get_DESeq_results(ddsList[[current_filter]], designList[[current_filter]], contrasts_subset, params[["design"]], params, current_filter)
    overallResListAll[[current_filter]] <- DESeq_results$resListAll
    overallResListFiltered[[current_filter]] <- DESeq_results$resListFiltered
    overallResListDEGs[[current_filter]] <- DESeq_results$resListDEGs
    mergedDEGsList[[current_filter]] <- DESeq_results$mergedDEGs
    filtered_table <- rbind(filtered_table, DESeq_results$filtered_table)
  }
}

summary_counts <- data.frame()
if (is.na(params$deseq_facet)) {
  resList <- overallResListDEGs[['all']]
  comparisons <- names(resList)
  for (comp in comparisons) { # by comparison
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
message(paste(capture.output(summary_counts), collapse = "\n"))

# save DESeq results to a file
save(ddsList,
     designList,
     contrastsList,
     overallResListAll,
     overallResListFiltered,
     overallResListDEGs,
     rldList,
     mergedDEGsList,
     exp_metadata,
     facets,
     count_data,
     exp_contrasts,
     params,
     paths,
     filtered_table,
     sample_count_metadata,
     file = file.path(paths$RData, paste0(params$project_title, "_DEG_data.RData")))

if (is.na(params$deseq_facet)) {
  facets <- c("all")
}

if (params$parallel) {
  parallel::mcmapply(
    FUN = R.ODAF.utils::write_tables,
    facet = facets,
    MoreArgs = list(params = params), # Additional static argument
    mc.cores = round(params$cpus * 0.6)
  )
} else {
  mapply(
    FUN = R.ODAF.utils::write_tables,
    facet = facets,
    MoreArgs = list(params = params) # The list of additional static arguments
  )
}
