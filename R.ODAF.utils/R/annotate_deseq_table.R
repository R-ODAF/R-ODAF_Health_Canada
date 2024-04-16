#' Annotate and Transform DESeq2 Results Table
#'
#' This function takes a list of DESeq2 results and annotates them with gene symbols
#' and descriptions. It performs transformation of log2 fold changes into linear
#' fold changes and applies optional filters.
#'
#' @param deseq_results_list A list of data frames, each containing DESeq2 results.
#' @param params List containing parameters such as platform, feature_id, species_data, etc.
#' @param filter_results Logical; if TRUE, apply additional filters to the table.
#' @param biosets_filter Logical; if TRUE, apply filters specifically for 'biosets'.
#' @return An annotated and transformed data frame of DESeq2 results.
#' @importFrom dplyr left_join mutate distinct
#' @importFrom data.table rbindlist
#' @export
annotate_deseq_table <- function(deseq_results_list,
                                 params,
                                 filter_results = FALSE,
                                 biosets_filter = FALSE) {
  annotated_results <- lapply(deseq_results_list, function(deg_table) {
    # Skip if the table is NULL or empty
    if (is.null(deg_table) || nrow(deg_table) == 0) {
      return(NULL)
    }

    # Proceed with the annotation and transformation
    feature_ids <- rownames(deg_table)
    contrast <- gsub(pattern = paste0("log2.*", params$design, "\ "), replacement = "", deg_table@elementMetadata[[2]][2])
    deg_table <- cbind(Feature_ID = feature_ids, as.data.frame(deg_table), contrast = contrast)

    if (params$platform == "TempO-Seq") {
      deg_table <- dplyr::left_join(deg_table, params$biospyder, by = c("Feature_ID" = params$feature_id))
    } else {
      # Annotation using orgdb package
      tryCatch({
        annotations <- AnnotationDbi::select(
          AnnotationDbi::loadDb(params$species_data$orgdb),
          columns = c("ENSEMBL", "SYMBOL", "GENENAME"),
          keys = feature_ids,
          keytype = "ENSEMBL")
        # Ensuring unique rows for annotations
        annotations <- dplyr::distinct(annotations, ENSEMBL, .keep_all = TRUE)
        colnames(annotations) <- c("Feature_ID", "Gene_Symbol", "description")
        deg_table <- dplyr::left_join(deg_table, annotations, by = "Feature_ID")
      }, error = function(e) {
        message("Error during annotation: ", e$message)
      })
    }

    # Transform log2 fold changes into linear fold changes
    deg_table <- deg_table %>%
      dplyr::mutate(linearFoldChange = ifelse(log2FoldChange > 0, 2^log2FoldChange, -1 / (2^log2FoldChange)),
        Gene_Symbol_2 = dplyr::coalesce(Gene_Symbol, Feature_ID)) %>%
      dplyr::select(Feature_ID,
#        Ensembl_Gene_ID,
        Gene_Symbol = Gene_Symbol_2,
        baseMean,
        log2FoldChange,
        linearFoldChange,
        lfcSE,
        pvalue,
        padj,
        contrast)

    # Apply biosets-specific filtering
    if (biosets_filter) {
      deg_table <- deg_table[
        !is.na(deg_table$pvalue) &
          deg_table$pvalue < params$alpha &
          abs(deg_table$linearFoldChange) > params$linear_fc_filter_biosets, ]
    }
    return(deg_table %>% dplyr::distinct())
  })

  # Combine the individual tables into one large table
  combined_results <- data.table::rbindlist(annotated_results, fill = TRUE)
  return(combined_results)
}





#   x <- deseq_results_list
#   annotated_results <- list()
#   for (i in 1:length(x)) {
#     deg_table <- x[[i]]
#     # Add taxonomy
#     if (is.null(deg_table)) {
#       next
#     } else if (nrow(deg_table) == 0) {
#       next
#     } else {
#       deg_table <- cbind(Feature_ID = row.names(x[[i]]),
#                          as(deg_table, "data.frame"),
#                          contrast = gsub(pattern = paste0("log2.*", params$design, "\ "),
#                                          replacement =  "",
#                                          x = x[[i]]@elementMetadata[[2]][2]))
#       if(params$platform == "TempO-Seq"){
#         deg_table <- dplyr::left_join(deg_table,
#                                       params$biospyder,
#                                       by = c(Feature_ID = params$feature_id))
#       } else{
#         # need to catch a testForValidKeys error in the case where the only resulting genes have ensembl IDs that aren't in the AnnotationDBI database
#         result = tryCatch({
#           descriptions <- AnnotationDbi::select(get(params$species_data$orgdb),
#                                                 columns = c("ENSEMBL", "SYMBOL", "GENENAME"),
#                                                 keys = deg_table$Feature_ID,
#                                                 keytype="ENSEMBL") %>%
#             distinct(ENSEMBL, .keep_all=TRUE)
#           colnames(descriptions) <- c("Ensembl_Gene_ID","Gene_Symbol","description")
#           descriptions$Feature_ID <- descriptions$Ensembl_Gene_ID
#           deg_table <- dplyr::left_join(deg_table, descriptions, by="Feature_ID")
#         }, error = function(e) {
#           message("Annotation error: ", e$message)
#         })
#       }
#       if(!("Gene_Symbol" %in% colnames(deg_table))){
#         deg_table$Ensembl_Gene_ID <- deg_table$Feature_ID
#         deg_table$Gene_Symbol <- NA
#       }
#       deg_table <- deg_table %>%
#         mutate(linearFoldChange = ifelse(log2FoldChange > 0, 2 ^ log2FoldChange, -1 / (2 ^ log2FoldChange))) %>%
#         mutate(Gene_Symbol_2 = coalesce(Gene_Symbol, Ensembl_Gene_ID)) %>%
#         dplyr::select(Feature_ID, Ensembl_Gene_ID, Gene_Symbol = Gene_Symbol_2, baseMean, log2FoldChange, linearFoldChange, lfcSE, pvalue, padj, contrast)
#
#       ## FILTERS ##
#       if (biosets_filter == TRUE) {
#         # for biosets, filter on unadjusted p-value
#         deg_table <- deg_table[!is.na(deg_table$pval) & deg_table$pval < params$alpha & abs(deg_table$linearFoldChange) > params$linear_fc_filter_biosets, ]
#       }
#       annotated_results[[i]] <- deg_table %>% dplyr::distinct()
#     }
#   }
#   y <- data.table::rbindlist(annotated_results)
#   return(y)
# }
