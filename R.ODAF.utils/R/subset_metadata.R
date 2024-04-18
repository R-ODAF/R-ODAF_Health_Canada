#' Subset Experiment Metadata Based on Facet and Filter Criteria
#'
#' This function subsets the experimental metadata based on facet and filter criteria.
#' It also subsets the contrasts data frame based on the subsets defined by this criteria.
#'
#' @param exp_metadata A data frame containing the experimental metadata.
#' @param design The name of the column used for the experimental design grouping.
#' @param contrasts A data frame of contrasts used for differential expression analysis.
#' @param current_facet A character string indicating the facet column to filter on.
#' @param current_filter A vector of values defining the filter criteria for the facet.
#' @return A list containing the subsetted experimental metadata and contrasts (each as a data frame).
#' @export
subset_metadata <- function(exp_metadata, design, contrasts, current_facet, current_filter) {
  contrasts_to_filter <- exp_metadata %>%
    dplyr::filter(!!sym(current_facet) %in% current_filter) %>% # NOTE: Not sure if %in% or == is better here.
    pull(design) %>%
    unique()
  contrasts_subset <- contrasts %>% dplyr::filter(V1 %in% contrasts_to_filter)
  if (params$strict_contrasts == TRUE) {
    contrasts_subset <- contrasts_subset %>% dplyr::filter(V2 %in% contrasts_to_filter)
  }
  exp_metadata_subset <- exp_metadata %>%
    dplyr::filter(!!sym(design) %in% (unlist(contrasts_subset) %>% unique())) # %>%
  # dplyr::filter(!!sym(current_facet) %in% current_filter)
  # The line above was added to deal with an edge case where samples were not properly filtered because current_facet wasn't in the contrast names
  # That edge case is uncommon and there are probably easier ways to deal with it
  # Including the commented dplyr code above breaks the analysis if controls aren't matched by chemical name
  # Ex. works for BPA_10 vs BPA_0 (as in test data), but not for BPA_10 vs DMSO

  # relevel the design and interesting groups
  exp_metadata_subset[[design]] <- factor(exp_metadata_subset[[design]],
    levels = unique(unlist(contrasts_subset)),
    ordered = FALSE)
  if (!is.na(params$sortcol)) {
    design_factor_reordered <- factor(exp_metadata_subset[[design]],
      levels = unique(exp_metadata_subset[[design]][mixedorder(exp_metadata_subset[[params$sortcol]])]),
      ordered = FALSE)
    exp_metadata_subset[[design]] <- design_factor_reordered

    for (ig in params$intgroup_to_plot) {
      intgroup_reordered <- factor(exp_metadata_subset[[ig]],
        levels = unique(exp_metadata_subset[[ig]][mixedorder(exp_metadata_subset[[params$sortcol]])]),
        ordered = FALSE)
      exp_metadata_subset[[ig]] <- intgroup_reordered
    }
  }
  return(list(exp_metadata = exp_metadata_subset, contrasts = contrasts_subset))
}
