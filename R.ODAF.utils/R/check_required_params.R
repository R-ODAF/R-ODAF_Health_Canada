#' Check for Required Parameters in a Configuration List
#'
#' This function ensures that all required parameters are set in the configuration list.
#' Stops and throws an error if any required parameter is missing.
#'
#' @param params Configuration list to check.
#' @export
check_required_params <- function(params) {
  required_params <- c("projectdir", "dose", "platform", "nmr_threshold")
  for (p in required_params) {
    if (is.na(params[p])) {
      stop(paste0("Required param ", p, " was NA. You should set that param to something in config.yaml."))
    } else if (is.null(params[p])) {
      stop(paste0("Required param ", p, " was null. This shouldn't happen, as you should be running 'replace_nulls_in_config"))
    } else if (!is.null(params$deseq_facet) && !is.na(params$deseq_facet)) {
      if (params$deseq_facet %in% params$intgroup_to_plot) {
        stop(paste0("The column for faceting (", params$deseq_facet, ") should not be an element in your intgroup_to_plot list."))
      }
    }
  }
}
