#' Read Contrasts File
#'
#' Reads the contrasts from a specified tab-delimited file.
#'
#' @param file_path A string indicating the file path of the contrasts file to read.
#' @param paths A list of paths used throughout the analysis.
#' @return A data frame containing the contrast definitions with columns as character data.
#' @export
get_contrasts <- function(file_path = NULL, paths) {
  if (is.null(file_path)) {
    message("No file path provided for contrasts. Using default path.")
    file_path <- file.path(paths[["contrasts"]], "contrasts.txt")
  }
  # Make sure file exists before attempting to read
  if (!file.exists(file_path)) {
    stop("Contrasts file does not exist at the specified path: ", file_path)
  }

  contrast_data <- read.delim(file_path,
                              stringsAsFactors = FALSE,
                              sep = "\t",
                              header = FALSE,
                              quote = "\"",
                              colClasses = "character")

  return(contrast_data)
}
