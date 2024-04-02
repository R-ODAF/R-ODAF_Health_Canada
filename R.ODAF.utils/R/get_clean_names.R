#' Obtain Clean Names for Visualization
#'
#' This function joins input names with annotations to create a clean set of identifiers, such as gene symbols,
#' for use in visualizations like heatmaps. It handles cases where multiple genes are mapped from a single feature.
#'
#' @param res Data frame containing the annotation information.
#' @param vector_of_names Character vector of feature identifiers to map to gene symbols.
#' @param feature_id The feature identifier column name in res used for mapping to vector_of_names.
#' @param species_gene_symbol The gene symbol column name in res used for mapping to gene names.
#' @return A character vector of clean names suitable for heatmaps or other visualizations.
#' @importFrom dplyr left_join distinct group_by mutate pull filter coalesce
#' @export
get_clean_names <- function(res, vector_of_names, feature_id, species_gene_symbol) {
  # Rename this to map_features_to_gene_symbols?
  hm_names <- data.frame(genes = vector_of_names)
  names(hm_names)[names(hm_names) == "genes"] <- feature_id
  heatmap_annotation_table <- res %>%
    dplyr::select(!!ensym(feature_id),
      !!ensym(species_gene_symbol))
  hm_names <- dplyr::left_join(hm_names,
    heatmap_annotation_table,
    by = feature_id) %>%
    dplyr::distinct() %>%
    dplyr::group_by(!!ensym(feature_id)) %>%
    # May be tricky to determine if more than one gene should be combined in a row.
    # Instead, currently using first non-NA value, since it's possible to have NAs sneak in
    # particularly at this point... high variance genes might not be picked up as DEGs.
    # dplyr::mutate(combined = paste0(!!ensym(species_gene_symbol), collapse = ", ")) %>%
    dplyr::mutate(combined = coalesce(!!ensym(species_gene_symbol), !!ensym(feature_id))) %>%
    dplyr::select(-!!ensym(species_gene_symbol)) %>%
    slice(1) %>% # This is bit of a hack - ran into problems with TempO-Seq probes that map to multiple genes.
    dplyr::distinct()
  hm_names <- hm_names[match(vector_of_names, hm_names[[feature_id]]), ]
  hm_names <- hm_names %>% dplyr::pull(combined)
  return(hm_names)
}
