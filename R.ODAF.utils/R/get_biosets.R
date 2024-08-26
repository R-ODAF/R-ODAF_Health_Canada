#' Generate biosets
#'
#' This function takes DEGs (after R-ODAF filters are applied) and applies additional filters
#' It writes biosets files and returns a biosets objects
#'
#' @param params A list of parameters from which the ID is constructed.
#' @param facet A character string indicating the current facet (subset of data) being processed.
#' @importFrom AnnotationDbi loadDb dbfile
#' @importFrom dplyr left_join distinct mutate arrange filter group_by ungroup
#' @importFrom stringr str_trunc str_replace_all str_split
#' @return A character string representing the unique analysis ID.
#' @export
#' 
 write_biosets <- function(params, facet) {
  current_filter <- facet
  message(paste0("Generating biosets for ", current_filter))

  biosets_unfiltered <- bioset_input[[current_filter]]

  filteredResults <- annotate_deseq_table(biosets_unfiltered, params, filter_results = TRUE, biosets_filter = TRUE)
  resultsContrasts <- stringr::str_split(filteredResults$contrast," vs ",2,TRUE)[,1]
  
  biosetsFilteredResults <- filteredResults %>%
      bind_cols(dose=resultsContrasts) %>%
      dplyr::select(Gene_Symbol, padj, linearFoldChange, dose)

  # Set up output paths
  if(is.na(params$deseq_facet)){
    biosets_folder <- paths$biosets_output
  } else{
    biosets_folder <- paths$biosets_output[[current_filter]]
  }
  
  biosets <- list() 

  all_doses <- unique(biosetsFilteredResults$dose)
  for(d in all_doses){
    bfr <- biosetsFilteredResults %>%
      filter(dose==d) %>%
      dplyr::select(Gene_Symbol, padj, linearFoldChange) %>%
      arrange(Gene_Symbol,-abs(linearFoldChange)) %>%
      distinct(Gene_Symbol, .keep_all=TRUE) 
    
    colnames(bfr) <- c("Gene","pval","fc")
    bfr <- bfr %>% arrange(desc(abs(fc)))
    # assume that current_filter includes the chemical + timepoint + dose
    biosets_name <- paste(d,params$units,params$celltype, sep='_')
    biosets_name <- paste0(str_replace_all(biosets_name, " ", "_"))
    biosets_fname <- paste0(biosets_name,".txt")
  
    write.table(bfr %>% mutate(across(where(is.numeric), ~ round(., digits = params$output_digits))),
                file = file.path(biosets_folder,biosets_fname),
                quote = FALSE, sep = "\t", col.names = NA)
    
    biosets[[biosets_name]] <- bfr %>% mutate(across(where(is.numeric), ~ round(., digits = params$output_digits)))
  } 
  return(biosets) # Need to figure out what to return (how to format)
 }
