#' Set Up Platform Specific Parameters
#'
#' Configures additional analysis parameters such as file paths and thresholds based on
#' the specified analysis platform.
#'
#' @param params A list containing initial parameters, including 'platform'
#'        and 'species_data' which are used to determine the specific settings.
#' @return A list of parameters updated with platform-specific settings.
#' @export
set_up_platform_params <- function(params) {
  species_data <- params$species_data
  params$sampledata_sep <- "\t"
  # set some additional parameters based on platform
  if (params$platform == "RNA-Seq") {
    params$MinCount <- 1
    params$alpha <- 0.05 # Relaxed from 0.01
    params$feature_id <- "Ensembl_Gene_ID"
  } else if (params$platform == "TempO-Seq") {
    params$MinCount <- 0.5
    params$alpha <- 0.05
    bs <- load_biospyder_new(params$biospyder_dbs, species_data$temposeq_manifest)
    params$feature_id <- bs$feature_id # Probe_Name
    params$biospyder <- bs$biospyder # manifest
  } else {
    stop("Platform/technology not recognized")
  }
  return(params)
}
