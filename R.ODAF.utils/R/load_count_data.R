#' Load Count Data from a File
#'
#' Reads count data from a specified file, typically a delimited text file.
#'
#' @param count_data_file A string giving the path to the count data file.
#' @param sampledata_sep The character string to be used for separating columns in the file. Default is "\\t".
#' @return A data frame with genes as rows and samples as columns.
#' @export
load_count_data <- function(count_data_file, sampledata_sep) {
  count_data <- read.delim(count_data_file,
    sep = sampledata_sep,
    stringsAsFactors = FALSE,
    header = TRUE,
    quote = "\"",
    row.names = 1,
    check.names = FALSE)
  return(count_data)
}
