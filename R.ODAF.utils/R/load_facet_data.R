#' Prepare Data For Report
#'
#' This function selects the appropriate data preparation function based on the
#' facet and display facet parameters and applies it to load necessary data objects
#' for the analysis report.
#'
#' @param paths A list of paths including the RData directory.
#' @param params A list of parameters used for the analysis.
#' @param facet_override A character string to override the facet parameter. Mostly used for troubleshooting.
#' @return An environment containing prepared data objects for reporting.
#' @importFrom matrixStats rowVars
#' @export
load_facet_data <- function(paths, params, facet_override = NA, filter_override = NA) {
  # Load data file into a new environment
  data_env <- new.env()
  data_file <- file.path(paths$RData, paste0(params$project_title, "_DEG_data.RData"))
  load(data_file, envir = data_env)
   if (!is.na(facet_override)) {
    if (is.na(filter_override)) {
      stop("You must also set the filter parameter when manually setting the facet parameter.")
    }
    message("You have manually set the facet parameter (i.e., column name from metadata) to: ", facet_override)
    message("Filtering the data based on the following filter: ", filter_override)
    params$reports_facet <- facet_override
    params$reports_filter <- filter_override
    if (length(data_env$ddsList) > 1) {
      params$deseq_facet <- facet_override
      params$deseq_filter <- filter_override
    }
   }
  # Create a list to store the result of the data preparation
  result <- list()

  # Case 1: DESeq2 on all samples; make reports for all samples.
  # Both "deseq_facet" and "reports_facet" are NA (unset).
  # (i.e., all groups and data in a single report)
  if (is.na(params$deseq_facet) && is.na(params$reports_facet)) {
    result <- prepare_data_case1(
      data_env$ddsList,
      data_env$overallResListAll,
      data_env$overallResListDEGs,
      data_env$rldList,
      data_env$mergedDEGsList,
      data_env$exp_metadata,
      data_env$contrastsList,
      data_env$allBiomarkers
    )
  # Case 2: DESeq2 on all samples; but, make faceted reports.
  # "deseq_facet" is NA, but "reports_facet" is set.
  } else if (is.na(params$deseq_facet) && !is.na(params$reports_facet)) {
    result <- prepare_data_case2(
      params,
      data_env$ddsList,
      data_env$overallResListAll,
      data_env$overallResListDEGs,
      data_env$rldList,
      data_env$mergedDEGsList,
      data_env$exp_metadata,
      data_env$designList, # Not in case 1
      data_env$contrastsList,
      data_env$allBiomarkers
    )
  # Case 3: DESeq2 is faceted; reports are faceted.
  # The two facets must match.
  } else if (!is.na(params$deseq_facet) && !is.na(params$reports_facet)) {
    if (params$deseq_facet != params$reports_facet) {
      stop("Error: reports_facet must match deseq_facet, otherwise DESeq2 results get mixed and matched.")
    }
    result <- prepare_data_case3(
      params,
      data_env$ddsList,
      data_env$overallResListAll,
      data_env$overallResListDEGs,
      data_env$rldList,
      data_env$mergedDEGsList,
      data_env$designList,
      data_env$contrastsList,
      data_env$allBiomarkers
    )
  # Case 4: DESeq2 is faceted, reports are not: this one doesn't make sense, since it could mislead end-users.
  } else {
    stop("Making a single report for faceted data not supported. Did you forget to set reports_facet?")
  }
  # Populate the environment with the result
  list2env(result, envir = .GlobalEnv)

  # filter the regularized data a couple ways for different displays
  result$rld_DEGs <- rld[row.names(assay(rld)) %in% mergedDEGs]

  rv <-  matrixStats::rowVars(assay(rld), useNames = FALSE)
  select <- order(rv, decreasing = TRUE)[1:params$nBest]
  result$rld_top <- rld[select, ]
  select_heatmap <- order(rv, decreasing = TRUE)[1:params$nHeatmapDEGs]
  result$rld_top_heatmap <- rld[select_heatmap, ]

  result$allResults <- annotate_deseq_table(resultsListAll, params, filter_results = FALSE)
  result$significantResults <- annotate_deseq_table(resultsListDEGs, params, filter_results = FALSE)

  result$ordered_contrast_strings <- contrasts_subset %>% mutate(contrast_string = paste(V1, 'vs', V2, sep = " ")) %>% pull(contrast_string)

  result$allResults$contrast <- factor(result$allResults$contrast, levels = result$ordered_contrast_strings)
  result$significantResults$contrast <- factor(result$significantResults$contrast, levels = result$ordered_contrast_strings)
  result$sample_count_metadata <- data_env$sample_count_metadata
  result$filtered_table <- data_env$filtered_table
  result$exp_metadata <- data_env$exp_metadata
  return(list2env(result, .GlobalEnv))
}
