#' Read Experiment Metadata
#'
#' Reads the experimental metadata from a specified file and adds a
#' new column with the original row names.
#'
#' @param file_path A string indicating the file path of the metadata file to read.
#' @param paths A list of paths used throughout the analysis.
#' @return A data frame containing the experiment metadata with an additional
#' column, 'original_names', which stores the original row names.
#' @export
get_metadata <- function(file_path = NULL, paths) {
  if (is.null(file_path)) {
    message("No file path provided for metadata. Using default path.")
    file_path <- file.path(paths[["qc"]], "metadata.QC_applied.txt")
  }
  # Make sure file exists before attempting to read
  if (!file.exists(file_path)) {
    stop("Metadata file does not exist at the specified path: ", file_path)
  }
  exp_metadata <- read.delim(file_path,
                             stringsAsFactors = FALSE,
                             sep = "\t",
                             header = TRUE,
                             quote = "\"",
                             colClasses = "character",
                             row.names = 1) # Ensure the column has unique IDs

  exp_metadata[["original_names"]] <- rownames(exp_metadata)

  return(exp_metadata)
}
