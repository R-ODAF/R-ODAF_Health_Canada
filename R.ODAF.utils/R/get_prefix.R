#' Generate a Filename Prefix
#'
#' Constructs a filename prefix based on provided parameters and an optional facet.
#' The prefix includes the platform, project title, and current time, with optional
#' inclusion of a facet or group filter. Spaces are replaced with underscores and
#' the resulting string is sanitized to ensure it is a valid filesystem path.
#'
#' @param prefix_pars A list containing parameters used to construct the prefix,
#'   including 'platform', 'project_title', and optionally 'reports_facet'
#'   and 'reports_filter'.
#' @param prefix_facet An optional character string representing the facet to include
#'   in the prefix. If 'reports_facet' is NA, the facet is ignored.
#' @return A character string representing the sanitized filename prefix.
#' @importFrom fs path_sanitize
#' @export
get_prefix <- function(prefix_pars, prefix_facet) {
   pars <- prefix_pars
   facet <- prefix_facet
   # Are there ways this logic might break?
   if (is.na(pars$reports_facet)) {
      message("Writing a single report for whole experiment.")
      prefix <- paste0(pars$platform, "_",
         pars$project_title, "_",
         format(Sys.time(), '%d-%m-%Y.%H.%M'))
   } else if (any(!is.na(pars$reports_filter))) {
      message(paste0("The group(s) of interest is (are) ",
         paste(pars$reports_filter, collapse = " and "), ".\n",
         "Writing a single report for that (those) groups."))
      prefix <- paste0(pars$platform, "_",
         pars$project_title, "_",
         paste(pars$reports_filter, collapse = "_"), "_",
         format(Sys.time(), '%d-%m-%Y.%H.%M'))
   } else {
      message(paste0("Building report for ", facet, "..."))
      prefix <- paste0(pars$platform, "_",
         pars$project_title, "_",
         facet, "_",
         format(Sys.time(), '%d-%m-%Y.%H.%M'))
   }
   prefix <- gsub(" ", "_", prefix)
   prefix <- fs::path_sanitize(prefix)
   return(prefix)
}
