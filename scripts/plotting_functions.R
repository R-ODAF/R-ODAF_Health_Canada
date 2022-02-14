get_clean_names <- function(res, vector_of_names, biomart_filter, species_gene_symbol) {
  hm_names <- data.frame(genes = vector_of_names)
  names(hm_names)[names(hm_names) == "genes"] <- biomart_filter
  heatmap_annotation_table <- res %>%
    dplyr::select(!!ensym(biomart_filter),
                  !!ensym(species_gene_symbol))
  hm_names <- dplyr::left_join(hm_names,
                               heatmap_annotation_table,
                               by = biomart_filter) %>%
    dplyr::distinct() %>%
    dplyr::group_by(!!ensym(biomart_filter)) %>%
    # May be tricky to determine if more than one gene should be combined in a row.
    # Instead, currently using first non-NA value, since it's possible to have NAs sneak in
    # particularly at this point... high variance genes might not be picked up as DEGs.
    #dplyr::mutate(combined = paste0(!!ensym(species_gene_symbol), collapse = ", ")) %>%
    dplyr::mutate(combined = coalesce(!!ensym(species_gene_symbol), !!ensym(biomart_filter))) %>%
    dplyr::select(-!!ensym(species_gene_symbol)) %>%
    slice(1) %>% # This is bit of a hack - ran into problems with TempO-Seq probes that map to multiple genes.
    dplyr::distinct()
  hm_names <- hm_names[match(vector_of_names, hm_names[[biomart_filter]]), ]
  hm_names <- hm_names %>% dplyr::pull(combined)
  return(hm_names)
}
