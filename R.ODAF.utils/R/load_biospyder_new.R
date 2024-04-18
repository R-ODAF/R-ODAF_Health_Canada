#' Load BioSpyder Data
#'
#' This function loads the BioSpyder assay data from a specified file, typically a CSV format.
#'
#' @param biospyder_dbs Path to the BioSpyder databases.
#' @param temposeq_manifest Name of the TempO-Seq assay manifest file.
#' @return A list containing the BioSpyder data frame and related information.
#' @export
#' @importFrom utils read.delim
load_biospyder_new <- function(biospyder_dbs, temposeq_manifest) {
  return_data = list()
  biospyder <- read.delim(file.path(biospyder_dbs, temposeq_manifest), # Assay manifest...
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    quote = "\"")

  feature_id <- "Probe_Name"
  biospyder[biospyder == "NULL"] <- NA
  return_data$biospyder <- biospyder
  return_data$feature_id <- feature_id
  return(return_data)
}
