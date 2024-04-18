#' Subset Count Data to Match Samples in Metadata
#'
#' Aligns and subsets count data to correspond with samples present in the provided experimental metadata.
#'
#' @param count_data A data frame (or DESeq2 model object) containing count data.
#' @param exp_metadata A data frame containing the experimental metadata to match against.
#' @return A subsetted count data frame with the columns ordered to
#'         correspond to the samples specified in the metadata.
#' @export
subset_data <- function(count_data, exp_metadata) {
  # Reorder the metadata table to correspond to the order of columns in the count data
  exp_metadata_sorted <- exp_metadata[exp_metadata$original_names %in% colnames(count_data), ]
  count_data_subset <- count_data[, exp_metadata_sorted$original_names]
  return(count_data_subset)
}
