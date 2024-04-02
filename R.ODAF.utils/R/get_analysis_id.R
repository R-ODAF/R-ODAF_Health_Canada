#' Generate an Analysis ID
#'
#' This function generates a unique ID for the analysis based on the current time and the project title.
#' If a `deseq_filter` is provided, it is appended to the ID string.
#'
#' @param params A list of parameters from which the ID is constructed.
#' @return A character string representing the unique analysis ID.
#' @export
get_analysis_id <- function(params) {
  # Set analysis ID. This ID will be used as prefix for the output files
  # Normally, as follows: year - project_title - deseq_filter
  if (is.na(params$deseq_filter) || is.null(params$deseq_filter)) {
    analysisID <- paste(format(Sys.time(), '%Y'), params$project_title, sep = "_")
  } else {
    analysisID <- paste(format(Sys.time(), '%Y'),
      params$project_title,
      paste(params$deseq_filter, collapse = "_"),
      sep = "_")
  }
  return(analysisID)
}
