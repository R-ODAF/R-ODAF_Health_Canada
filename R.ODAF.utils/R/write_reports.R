
#' Generate Main Report
#'
#' Renders the main report HTML file for a given facet using the specified parameters.
#' The function updates the `reports_filter` in `pars` based on the facet provided.
#'
#' @param pars A list of parameters used for report generation.
#' @param facet A character string specifying the facet for which to generate the report.
#' @param paths A list of paths used throughout the analysis.
#' @return Invisible NULL. The function is called for its side effect of rendering an HTML report.
#' @export
make_main_reports <- function(pars, paths, facet) {
  message("Making main report")
  if (is.na(pars$deseq_facet) && is.na(pars$reports_facet)) {
    pars$reports_filter <- NULL
  } else {
    pars$reports_filter <- facet
  }
  if (is.null(paths)) {
    stop("Paths not provided")
  }
  prefix <- get_prefix(prefix_pars = pars, prefix_facet = facet)
  if (pars$generate_main_report) {
    main_report <- file.path(pars$projectdir, "Rmd", "DESeq2_report_new.Rmd")
    main_file <- file.path(paths$reports_dir, paste0(prefix, ".html"))
    render_report(main_report, main_file, pars)
  }
}

#' Generate Statistics Report
#'
#' Renders the statistics report HTML file for a given facet using the specified parameters.
#' The function updates the `reports_filter` in `pars` based on the facet provided.
#'
#' @param pars A list of parameters used for report generation.
#' @param facet A character string specifying the facet for which to generate the report.
#' @param paths A list of paths used throughout the analysis.
#' @return Invisible NULL. The function is called for its side effect of rendering an HTML report.
#' @export
make_stats_reports <- function(pars, paths, facet) {
  if (is.na(pars$deseq_facet) && is.na(pars$reports_facet)) {
    pars$reports_filter <- NULL
  } else {
    pars$reports_filter <- facet
  }
  prefix <- get_prefix(prefix_pars = pars, prefix_facet = facet)
  if (pars$generate_stats_report) {
    message("Generating stats report")
    stats_report <- file.path(paths$projectdir, "Rmd", "stats_report.Rmd")
    stats_file <- file.path(paths$reports_dir, paste0("stats_", prefix, ".html"))
    options(pandoc.stack.size = "128m")
    render_report(stats_report, stats_file, pars)
  }
}

#' Generate Data Explorer Report
#'
#' Renders the data explorer report HTML file for a given facet using the specified parameters.
#' The function updates the `reports_filter` in `pars` based on the facet provided.
#'
#' @param pars A list of parameters used for report generation.
#' @param facet A character string specifying the facet for which to generate the report.
#' @param paths A list of paths used throughout the analysis.
#' @return Invisible NULL. The function is called for its side effect of rendering an HTML report.
#' @export
make_data_reports <- function(pars, paths, facet) {
  if (is.na(pars$deseq_facet) && is.na(pars$reports_facet)) {
    pars$reports_filter <- NULL
  } else {
    pars$reports_filter <- facet
  }
  prefix <-get_prefix(prefix_pars = pars, prefix_facet = facet)
  if (pars$generate_data_explorer_report) {
    message("Generating data explorer report")
    data_explorer_report <- file.path(paths$projectdir, "Rmd", "data_explorer_report.Rmd")
    data_explorer_file <- file.path(paths$reports_dir, paste0("data_explorer_", prefix, ".html"))  
    render_report(data_explorer_report, data_explorer_file, pars)
  }
}

#' Generate GO and Pathway Analysis Report
#'
#' Renders the GO and pathway analysis report HTML file for a given facet using the specified parameters.
#' The function updates the `reports_filter` in `pars` based on the facet provided.
#'
#' @param pars A list of parameters used for report generation.
#' @param facet A character string specifying the facet for which to generate the report.
#' @param paths A list of paths used throughout the analysis.
#' @return Invisible NULL. The function is called for its side effect of rendering an HTML report.
#' @export
make_pathway_reports <- function(pars, paths, facet)  {
  if (is.na(pars$deseq_facet) && is.na(pars$reports_facet)) {
    pars$reports_filter <- NULL
  } else {
    pars$reports_filter <- facet
  }
  prefix <- get_prefix(prefix_pars = pars, prefix_facet = facet)
  if (pars$generate_go_pathway_report) {
    message("Generating GO and pathway analysis report")
    go_pathway_report <- file.path(paths$projectdir, "Rmd", "go_pathway_report.Rmd")
    go_pathway_file <- file.path(paths$reports_dir, paste0("go_pathway_", prefix, ".html"))
    render_report(go_pathway_report, go_pathway_file, pars)
  }
}

#' Generate TGx-DDI Report
#'
#' Renders the TGx-DDI report HTML file for a given facet using the specified parameters.
#' The function updates the `reports_filter` in `pars` based on the facet provided.
#'
#' @param pars A list of parameters used for report generation.
#' @param facet A character string specifying the facet for which to generate the report.
#' @param paths A list of paths used throughout the analysis.
#' @return Invisible NULL. The function is called for its side effect of rendering an HTML report.
#' @export
make_tgxddi_reports <- function(pars, paths, facet) {
  if (is.na(pars$deseq_facet) && is.na(pars$reports_facet)) {
    pars$reports_filter <- NULL
  } else {
    pars$reports_filter <- facet
  }
  prefix <-get_prefix(prefix_pars = pars, prefix_facet = facet)
  if (pars$generate_tgxddi_report) {
    message("Generating TGX-DDI report")
    tgxddi_report <- file.path(paths$projectdir, "Rmd", "tgx_ddi.Rmd")
    tgxddi_file <- file.path(paths$reports_dir, paste0("tgxddi_", prefix, ".html"))  
    render_report(tgxddi_report, tgxddi_file, pars)
  }
}

#' Generate TGx-HDACi Report
#'
#' Renders the TGx-HDACi report HTML file for a given facet using the specified parameters.
#' The function updates the `reports_filter` in `pars` based on the facet provided.
#'
#' @param pars A list of parameters used for report generation.
#' @param facet A character string specifying the facet for which to generate the report.
#' @param paths A list of paths used throughout the analysis.
#' @return Invisible NULL. The function is called for its side effect of rendering an HTML report.
#' @export
make_hdaci_reports <- function(pars, paths, facet) {
  if (is.na(pars$deseq_facet) && is.na(pars$reports_facet)) {
    pars$reports_filter <- NULL
  } else {
    pars$reports_filter <- facet
  }
  prefix <-get_prefix(prefix_pars = pars, prefix_facet = facet)
  if (pars$generate_tgxhdaci_report) {
    message("Generating TGX-HDACi report")
    tgxhdaci_report <- file.path(paths$projectdir, "Rmd", "tgx-hdaci.Rmd")
    tgxhdaci_file <- file.path(paths$reports_dir, paste0("tgxhdaci_", prefix, ".html"))  
    render_report(tgxhdaci_report, tgxhdaci_file, pars)
  }
}