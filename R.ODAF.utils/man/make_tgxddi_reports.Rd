% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_reports.R
\name{make_tgxddi_reports}
\alias{make_tgxddi_reports}
\title{Generate TGx-DDI Report}
\usage{
make_tgxddi_reports(pars, paths, facet)
}
\arguments{
\item{pars}{A list of parameters used for report generation.}

\item{paths}{A list of paths used throughout the analysis.}

\item{facet}{A character string specifying the facet for which to generate the report.}
}
\value{
Invisible NULL. The function is called for its side effect of rendering an HTML report.
}
\description{
Renders the TGx-DDI report HTML file for a given facet using the specified parameters.
The function updates the \code{reports_filter} in \code{pars} based on the facet provided.
}
