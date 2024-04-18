#' Get Facets from Metadata
#'
#' This function extracts unique values from a specified column in the metadata,
#' excluding any groups specified by the user and an additional skip value.
#'
#' @param metadata A data frame containing the experimental metadata.
#' @param skip_extra A character string specifying an additional group name to be
#'   excluded from the facets.
#' @param single_facet_constant A character string specifying the value to return
#'  if no facets are found.
#' @param params A list of parameters for the analysis.
#' @return A character vector of unique facet values after exclusion.
#' @export
get_facets <- function(metadata,
                       params,
                       skip_extra = c("DMSO"),
                       single_facet_constant = "single_facet_constant_12345") {
  if (is.null(metadata)) {
    message("No metadata provided. Skipping facet extraction.")
    return(NA)
  }
  # Case 1: DESeq2 on all samples; make reports for all samples.
  # Both "deseq_facet" and "reports_facet" are NA (unset).
  # (i.e., all groups and data in a single report)
  if (is.na(params[["deseq_facet"]]) && is.na(params[["reports_facet"]])) {
    report_facets <- single_facet_constant
    # Case 2: DESeq2 on all samples; but, make faceted reports.
    # "deseq_facet" is NA, but "reports_facet" is set.
  } else if (is.na(params[["deseq_facet"]]) && !is.na(params[["reports_facet"]])) {
    report_facets <- parse_facets(metadata, params, skip_extra)
    # Case 3: DESeq2 is faceted; reports are faceted.
    # The two facets must match.
  } else if (!is.na(params[["deseq_facet"]]) && !is.na(params[["reports_facet"]])) {
    if (params[["deseq_facet"]] != params[["reports_facet"]]) {
      stop("Error: reports_facet must match deseq_facet, otherwise DESeq2 results get mixed and matched.")
    }
    report_facets <- parse_facets(metadata, params, skip_extra)
    # Case 4: DESeq2 is faceted, reports are not: this one doesn't make sense, since it could mislead end-users.
  } else {
    stop("Making a single report for faceted data not supported. Did you forget to set reports_facet?")
  }
  return(report_facets)
}


#' Parse facet names for generating reports
#'
#' This function takes in metadata, parameters, and a flag to skip extra groups,
#' and returns the facets based on which multiple reports will be generated.
#'
#' @param metadata The metadata containing information about the groups.
#' @param params The parameters specifying the reports facet and exclude groups.
#' @param skip_extra A vector indicating extra groups to skip.
#' @importFrom dplyr filter pull
#'
#' @return A vector of unique facets based on which multiple reports will be generated.
#' If the reports facet is not specified, it returns NA indicating a single report for all groups.
#'
#' @export
parse_facets <- function(metadata, params, skip_extra) {
  if (!is.na(params[["reports_facet"]])) {
    exclude_groups <- c(params[["exclude_groups"]], params[["solvent_control"]], skip_extra)

	facets <- metadata %>%
	  dplyr::filter(
		!(params[["reports_facet"]] %in% exclude_groups),
		solvent_control == FALSE
	  ) %>%
	  dplyr::pull(params[["reports_facet"]]) %>%
	  unique()

    message(paste0("Facets will be based on ", params[["reports_facet"]], "."))
    return(facets)
  } else {
    message("Making a single report for all groups...")
    return(NA)
  }
}
