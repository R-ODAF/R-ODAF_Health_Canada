#' Process Count Data and Experimental Metadata
#'
#' This function applies a series of processing steps to count data and experimental metadata,
#' which includes filtering, formatting, sorting and performing sanity checks.
#'
#' @param count_data A matrix or data frame of raw count data from transcriptomic analysis.
#' @param exp_metadata A data frame containing the experimental metadata.
#' @param exp_contrasts A data frame of contrasts used for differential expression.
#' @param intgroup A character vector with names of columns used as interest groups in the analysis.
#' @param design The name of the column used for the design in the experimental metadata.
#' @param params A list of parameters defining specific aspects of processing.
#' @return A list containing processed count data, experimental metadata, and contrasts.
#' @importFrom gtools mixedsort
#' @export
process_data_and_metadata <- function(count_data, exp_metadata, exp_contrasts, params) {
  message("Filtering metadata...")
  exp_metadata <- filter_metadata(exp_metadata, params, exp_contrasts)
  message("Formatting and sorting metadata...")
  exp_metadata <- format_and_sort_metadata(exp_metadata, params$design, params$sortcol)
  message("Filtering count data...")
  count_data <- filter_data(count_data, exp_metadata, params$nmr_threshold)
  if (!is.na(params$sortcol)) {
    message("Sorting contrasts...")
    exp_contrasts <- sort_contrasts(exp_metadata, exp_contrasts, params, contrast_col = 1)
  }
  message("Sanity checks for data...")
  check_data(count_data, exp_metadata, exp_contrasts)
  return(list(count_data = count_data, exp_metadata = exp_metadata, contrasts = exp_contrasts))
}

#' Filter experimental metadata based on inclusion/exclusion criteria
#'
#' This function processes the experimental metadata for an
#' experiment by excluding or including certain samples or groups based
#' on user-defined criteria.
#'
#' @param exp_metadata A data frame containing the experiment metadata.
#' @param params The list of parameters describing the study.
#' @param design The name of the column used for the experimental design grouping.
#' @import rlang
#' @import dplyr
#' @return A filtered metadata data frame.
#' @export
filter_metadata <- function(exp_metadata, params, exp_contrasts) {
  # exclude samples
  if (any(!is.na(params$exclude_samples))) {
    exp_metadata <- exp_metadata %>%
      dplyr::filter(!original_names %in% params$exclude_samples)
  }
  # exclude groups
  if (any(!is.na(params$exclude_groups))) {
    exp_metadata <- exp_metadata %>%
      dplyr::filter(!.data[[params$design]] %in% params$exclude_groups)
    contrasts_to_filter <- exp_metadata %>%
      dplyr::filter(!.data[[params$design]] %in% params$exclude_groups) %>%
      pull(.data[[params$design]]) %>%
      unique()
    exp_contrasts <- exp_contrasts %>%
      dplyr::filter(V1 %in% contrasts_to_filter)
    if (params$strict_contrasts == TRUE) {
      exp_contrasts <- exp_contrasts %>%
        dplyr::filter(V2 %in% contrasts_to_filter)
    }
  }
  if (!is.na(params$include_only_column) && !is.na(params$include_only_group)) {
    exp_metadata <- exp_metadata %>%
      dplyr::filter((!!sym(params$include_only_column)) %in% params$include_only_group)
    limit_contrasts <- exp_metadata %>%
      pull(!!sym(params$design)) %>%
      unique() %>%
      as.character()
    exp_contrasts <- exp_contrasts %>% dplyr::filter(V1 %in% limit_contrasts)
  }
  return(exp_metadata)
}


#' Format and sort experiment metadata for downstream analysis
#'
#' This function formats and sorts the experimental metadata by converting
#' grouping variables (`intgroup`) to factors and sorting by another
#' column (`sortcol`) if specified.
#'
#' @param exp_metadata A data frame of the experiment metadata.
#' @param intgroup A character vector of column names to be converted to factors for grouping.
#' @param design The name of the column used for the experimental design grouping.
#' @param sortcol An optional column used for sorting the `design` grouping.
#' @return The formatted and sorted metadata data frame.
#' @importFrom gtools mixedorder
#' @export
format_and_sort_metadata <- function(exp_metadata, design, sortcol) {
  # Design needs to be a factor for DESeq2
  exp_metadata[design] <- factor(exp_metadata[[design]], levels = gtools::mixedsort(unique(exp_metadata[[design]])))

  # if sortcol is defined, sort the design variable based on that
  if (!is.na(sortcol)) {
    design_factor_reordered <- factor(exp_metadata[[design]],
      levels = unique(exp_metadata[[design]][gtools::mixedorder(exp_metadata[[sortcol]])]),
      ordered = FALSE)
    exp_metadata[[design]] <- design_factor_reordered

    # also sort the other interesting groups we want to plot
    for (ig in params$intgroup_to_plot) {
      intgroup_reordered <- factor(exp_metadata[[ig]],
        levels = unique(exp_metadata[[ig]][mixedorder(exp_metadata[[sortcol]])]),
        ordered = FALSE)
      exp_metadata[[ig]] <- intgroup_reordered
    }
  }
  return(exp_metadata)
}

#' Filter Count Data Based on Minimum Counts
#'
#' Removes samples from count data that have total read counts below a specified threshold.
#' NA values within the count data are replaced with zeros.
#'
#' @param count_data A data frame containing count data.
#' @param exp_metadata A data frame containing the corresponding experimental metadata.
#' @param threshold A numeric threshold for minimum allowed count sum across samples.
#' @return A filtered count data matrix or data frame.
#' @export
filter_data <- function(count_data, exp_metadata, threshold) {
  # First data clean-up: replace NA & remove samples with total readcount < threshold
  count_data[is.na(count_data)] <- 0
  count_data <- count_data[, (colSums(count_data) > threshold)] # reads required per sample
  # count_data <- count_data[(rowSums(count_data) > 1),] # reads required per gene
  exp_metadata <- exp_metadata[exp_metadata$original_names %in% colnames(count_data), ]
  count_data <- count_data[, exp_metadata$original_names]
  return(count_data)
}

#' Sort contrasts based on experimental metadata and a sorting column
#'
#' Sorts the contrasts based on the order of entries in a specified sorting column
#' within the experimental metadata.
#'
#' @param exp_metadata A data frame containing the experiment metadata.
#' @param exp_contrasts A data frame containing contrast pairs to sort.
#' @param design The name of the column in `exp_metadata` that corresponds to the design variable.
#' @param sortcol The name of the column to sort by.
#' @return A data frame of sorted contrast pairs.
#' @importFrom gtools mixedorder
#' @importFrom rlang .data
#' @export
sort_contrasts <- function(exp_metadata, exp_contrasts, params, contrast_col = 1) {
  ordered_design <- exp_metadata[gtools::mixedorder(exp_metadata[, params$sortcol]), ] %>%
    dplyr::select(params$design) %>%
    dplyr::pull()
  message("Ordered experimental groups: ", ordered_design)
  message("Sorting contrasts based on:", names(exp_contrasts[[contrast_col]]))
  ordered_contrasts <- exp_contrasts %>%
    dplyr::slice(match(ordered_design, exp_contrasts[[contrast_col]])) %>%
    unique()
  if(nrow(exp_contrasts) == nrow(ordered_contrasts)) {
    return(ordered_contrasts)
  } else {
    return(exp_contrasts)
    warning("NOTE: Your contrasts will NOT be sorted, because the sortcol does not appear to have a one-to-one relationship with contrasts.")
  }
}

#' Perform Sanity Checks on Data
#'
#' Runs a series of sanity checks on count data and corresponding metadata
#' to ensure that they are consistent with each other for further analysis.
#'
#' @param sd A data frame of count data to check.
#' @param exp_metadata A data frame of the experimental metadata to check against count data.
#' @param exp_contrasts A data frame of contrasts to validate.
#' @return None; this function is called for its side effects.
#' @export
check_data <- function(sd, exp_metadata, exp_contrasts) {
  message("Sanity checks for data")
  # make sure they're not empty
  stopifnot(exprs = {
    ncol(sd) > 0
    nrow(exp_metadata) > 0
    nrow(exp_contrasts) > 0
  })
  # Sanity check: each sample (row) in the metadata should have a corresponding column in the count data
  stopifnot(all(exp_metadata$original_names %in% colnames(sd)))
  # Sanity check: each column in the count data should have a corresponding sample (row) in the metadata
  stopifnot(all(colnames(sd) %in% exp_metadata$original_names))
  message("All OK \U0001f44d")
}
