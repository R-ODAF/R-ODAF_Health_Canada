run_analysis_for_facet <- function(count_data, exp_metadata, exp_contrasts, design, params, current_filter) {
  metadata_subset <- subset_metadata(exp_metadata, design, exp_contrasts, params$deseq_facet, current_filter)
  exp_metadata_subset <- metadata_subset$exp_metadata
  contrasts_subset <- metadata_subset$contrasts
  count_data_subset <- subset_data(count_data, exp_metadata_subset)
  
  check_data(count_data_subset, exp_metadata_subset, contrasts_subset)
  
  if(params$write_additional_output) {
    write_additional_output(count_data_subset, exp_metadata_subset, design, params)
  }
  dds <- learn_deseq_model(count_data_subset, exp_metadata_subset, design, params)
  rld <- regularize_data(dds, original_design, covariates = NA, params$batch_var)
  DESeq_results <- get_DESeq_results(dds, exp_metadata_subset, contrasts_subset, design, params, current_filter, paths$DEG_output)

  list(
    dds = dds,
    rld = rld,
    resListAll = DESeq_results$resListAll,
    resListFiltered = DESeq_results$resListFiltered,
    resListDEGs = DESeq_results$resListDEGs,
    mergedDEGs = DESeq_results$mergedDEGs,
    filtered_table = DESeq_results$filtered_table
  )
}

# # Handle the no facet case
# if (is.na(params$deseq_facet)){
#   message("### Learning a single model for the whole experiment. ###")
#   results <- run_analysis_for_facet(count_data, exp_metadata, exp_contrasts, params[["design"]], params, NA)
#   # Now assign each of the results to the overall lists
#   ddsList[['all']] <- results$dds
#   overallResListAll[['all']] <- results$resListAll
#   # ... and so on for each result
# } else {
#   # Handle multiple facet case by iterating through facets
#   for (current_filter in facets) {
#     message(paste0("### Learning model for ", current_filter, ". ###"))
#     results <- run_analysis_for_facet(count_data, exp_metadata, exp_contrasts, params[["design"]], params, current_filter)
#     # Assign to overall collections
#     ddsList[[current_filter]] <- results$dds
#     overallResListAll[[current_filter]] <- results$resListAll
#     # ... and the rest of the assignments
#   }
# }
