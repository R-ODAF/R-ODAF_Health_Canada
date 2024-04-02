#' Write Additional Output for BMDexpress and Biomarker Analysis
#'
#' This function processes count data and experimental metadata to write out
#' files formatted for use with BMDexpress and biomarker analysis. It performs
#' the CPM normalization, log2 transformation, and aggregation as needed.
#'
#' @param count_data A matrix or data frame containing count data.
#' @param exp_metadata A data frame containing the experimental metadata.
#' @param design_to_use The column in `exp_metadata` that specifies the design.
#' @param params A list containing parameters for the analysis.
#' @importFrom edgeR cpm
#' @importFrom dplyr group_by summarize across
#' @importFrom utils write.table
#' @export
write_additional_output <- function(count_data,
                                    exp_metadata,
                                    design_to_use,
                                    params) {
  if (!is.na(params[["dose"]])) {
    cpm_data <- edgeR::cpm(count_data)
    bmdexpress <- as.data.frame(log2(cpm_data + 1))
    bmdexpress <- cbind(SampleID = c(row.names(bmdexpress)),
      bmdexpress,
      stringsAsFactors = FALSE)

    if (params$platform == "TempO-Seq") {
      biomarkers <- count_data
      biomarkers$gene <- rownames(count_data)
      biomarkers$gene <- gsub("_.*", "", biomarkers$gene)
      biomarkers <- biomarkers %>%
        dplyr::group_by(gene) %>%
        dplyr::summarize(across(where(is.numeric), sum))
      colnames(biomarkers)[[1]] <- "SampleID" # For reasons
      biomarkers[-1] <- cpm(biomarkers[-1])
      biomarkers[-1] <- as.data.frame(log2(biomarkers[-1] + 1))
    } else { biomarkers <- bmdexpress } # Still includes all genes

    bmdexpress <- bmdexpress[rowSums(count_data) > 5, ]
    # add a dose header line to both files
    bmdexpress <- rbind(c("Dose", as.character(exp_metadata[colnames(bmdexpress)[-1], ][[params$dose]])),
      bmdexpress,
      stringsAsFactors = FALSE)
    biomarkers <- rbind(c("Dose", as.character(exp_metadata[colnames(biomarkers)[-1], ][[params$dose]])),
      biomarkers,
      stringsAsFactors = FALSE)

    # Determine names of dose groups in which n per group > 1
    groups_for_bmdexpress <- which(table(t(bmdexpress[1, ])) > 1) %>% names()
    # Rewrite bmdexpress table
    # Manually include the "Dose" and gene name column
    bmdexpress <- bmdexpress[, (bmdexpress[1, ]) %in% c("Dose", groups_for_bmdexpress)]

    if (!is.na(params$deseq_facet)) {
      fname <- paste0("bmdexpress_input_",
        paste(current_filter,
          collapse = "_"),
        ".txt")
      fname2 <- paste0("biomarker_input_",
        paste(current_filter,
          collapse = "_"),
        ".txt")
      write.table(bmdexpress,
        file = file.path(paths$BMD_output,
          fname),
        quote = FALSE,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE)
      write.table(biomarkers,
        file = file.path(paths$BMD_output,
          fname2),
        quote = FALSE,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE)
    } else {
      write.table(bmdexpress,
        file = file.path(paths$BMD_output, "bmdexpress_input.txt"),
        quote = FALSE,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE)
      write.table(biomarkers,
        file = file.path(paths$BMD_output, "biomarkers_input.txt"),
        quote = FALSE,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE)
    }
  }

}
