#' Write DESeq2 Results and Annotations to Files
#'
#' This function takes the results from `DESeq2` analysis, annotates them, and writes out several files
#' including those necessary for downstream analysis and visualization.
#'
#' @param facet A character string indicating the current facet (subset of data) being processed.
#' @importFrom dplyr left_join distinct mutate arrange filter group_by ungroup
#' @importFrom data.table setDT setnames
#' @importFrom openxlsx createWorkbook addWorksheet writeDataTable createStyle freezePane saveWorkbook modifyBaseFont mergeCells writeData conditionalFormatting setColWidths
#' @importFrom AnnotationDbi loadDb dbfile
#' @importFrom edgeR cpm
#' @importFrom stats quantile median
#' @importFrom tidyr pivot_wider
#' @importFrom stringr str_trunc str_replace_all str_split
#' @importFrom fs path_sanitize
#' @importFrom DESeq2 counts
#' @export

write_tables <- function(facet, params) {
  db <- AnnotationDbi::loadDb(params$species_data$orgdb)
  current_filter <- facet
  message(paste0("Writing tables for ", current_filter))
  resultsListAll <- overallResListAll[[current_filter]] 
  resultsListDEGs <- overallResListDEGs[[current_filter]]
  if (length(resultsListDEGs) < 1) { return(message("No output for this facet.")) }
  resultsListFiltered <- overallResListFiltered[[current_filter]] # For BMDExpress
  dds <- ddsList[[current_filter]] 

  allResults <- annotate_deseq_table(resultsListAll, params, filter_results = FALSE)
  significantResults <- annotate_deseq_table(resultsListDEGs, params, filter_results = FALSE)


  Counts <- counts(dds, normalized = TRUE)
  CPMdds <- cpm(counts(dds, normalized = TRUE))

  if(params$platform == "TempO-Seq"){
    descriptions <- AnnotationDbi::select(db, columns = c("ENSEMBL", "GENENAME"), keys = allResults$Ensembl_Gene_ID, keytype="ENSEMBL") %>% distinct()
    colnames(descriptions) <- c("Ensembl_Gene_ID","description")
    id_table <- params$biospyder %>% left_join(descriptions) %>% dplyr::select(Feature_ID=Probe_Name, Gene_Symbol, Ensembl_Gene_ID, description) # this is annoying: could select columns using contains("Gene_Symbol", ignore.case =TRUE)
  } else {
    id_table <- AnnotationDbi::select(db, columns = c("ENSEMBL", "SYMBOL", "GENENAME"), keys = overallAllGenes$Ensembl_Gene_ID, keytype="ENSEMBL") %>% distinct()
    colnames(id_table) <- c("Ensembl_Gene_ID","Gene_Symbol","description")
    id_table$Feature_ID <- id_table$Ensembl_Gene_ID
  }
  summaryTable <- allResults %>%
    dplyr::select(Feature_ID, baseMean) %>%
    distinct()

  contrastsInSummary <- c()

  prefix <- paste0(params$platform, "_",
                   params$project_title, "_",
                   current_filter, "_",
                   format(Sys.time(),'%d-%m-%Y.%H.%M'))  
  prefix <- str_replace_all(prefix, " ", "_")
  prefix <- fs::path_sanitize(prefix)

  for (i in seq_along(resultsListAll)) {
    message(resultsListAll[[i]]@elementMetadata[[2]][2])
    q <- gsub(pattern = paste0("log2\ fold\ change\ \\(MMSE\\):\ ", params$design, "\ "),
              replacement =  "",
              x = resultsListAll[[i]]@elementMetadata[[2]][2])
    toJoin <- as.data.frame(resultsListAll[[i]])
    setDT(toJoin, keep.rownames = TRUE)[]
    setnames(toJoin, 1, "Feature_ID")
    toJoin <- mutate(toJoin, linearFoldChange = ifelse(log2FoldChange > 0,
                                                       2 ^ log2FoldChange,
                                                       -1 / (2 ^ log2FoldChange)))
    toJoin <- toJoin[, c(1:3, 7, 4:6)]
    summaryTable <- dplyr::left_join(summaryTable, dplyr::select(toJoin, !c(baseMean, pvalue, lfcSE)), by = "Feature_ID")

    names(summaryTable)[[ncol(summaryTable) - 2]] <- paste0("log2FoldChange_", i)
    names(summaryTable)[[ncol(summaryTable) - 1]] <- paste0("linearFoldChange_", i)
    names(summaryTable)[[ncol(summaryTable)]] <- paste0("FDR_", i)
    contrastsInSummary[i] <- q
    message(summaryTable %>% nrow())
  }

  message("getting final output tables")
  maxFCs <- allResults %>%
    dplyr::group_by(Feature_ID) %>%
    dplyr::filter(abs(linearFoldChange) == max(abs(linearFoldChange))) %>%
    dplyr::ungroup() %>%
    dplyr::select(Feature_ID, linearFoldChange)

  minPvals <- allResults %>%
    dplyr::group_by(Feature_ID) %>%
    dplyr::filter(padj == min(padj)) %>%
    dplyr::ungroup() %>%
    dplyr::select(Feature_ID, padj)

  summaryTable <- summaryTable %>%
    left_join(id_table, by = "Feature_ID") %>%
    left_join(maxFCs, by = "Feature_ID") %>%
    left_join(minPvals, by = "Feature_ID") %>%
    dplyr::rename(maxFoldChange = linearFoldChange,
                  minFDR_pval = padj) %>%
    dplyr::distinct() %>%
    mutate(maxFoldChange = abs(maxFoldChange)) # This eliminates the direction of change: this way it's easy to sort.

  numColsToPrepend <- ncol(summaryTable) - 3*length(resultsListAll) - 2 # Number of columns per contrast = 3. Subtract two for the baseMean and genes columns.
  colPositionsToPrependSTART <- ncol(summaryTable) - numColsToPrepend + 1
  colPositionsOfData <- ncol(summaryTable) - numColsToPrepend
  summaryTable <- as.data.frame(summaryTable)
  summaryTable <- summaryTable[, c(1,
                                   colPositionsToPrependSTART:ncol(summaryTable),
                                   2:colPositionsOfData)]
  summaryTable <- summaryTable %>% dplyr::distinct() # Just in case duplicates snuck by

  CPMddsDF <- data.frame(genes = row.names(CPMdds), CPMdds, check.names = FALSE)
  CPMddsDF <- dplyr::left_join(CPMddsDF, id_table, by = c("genes" = "Feature_ID"))
  numColsToPrepend <- ncol(CPMddsDF) - ncol(CPMdds) - 1
  colPositionsToPrependSTART <- ncol(CPMddsDF) - numColsToPrepend + 1
  colPositionsOfData <- ncol(CPMddsDF) - numColsToPrepend
  CPMddsDF <- CPMddsDF[, c(1, colPositionsToPrependSTART:ncol(CPMddsDF), 2:colPositionsOfData)]

  if(is.na(params$deseq_facet)){
    output_folder <- paths$DEG_output
  } else{
    output_folder <- paths$DEG_output[[current_filter]]
  }

  #######################################
  ### Write results table from DESeq2
  #######################################
  message("write results tables to txt")
  message(paste0("Rounding numeric data to ",params$output_digits," digits."))
  write.table(allResults %>% mutate(across(where(is.numeric), ~ round(., digits = params$output_digits))),
              file = file.path(output_folder,
                               paste0(prefix,"-DESeq_output_ALL.txt")),
              quote = FALSE, sep = "\t", col.names = NA)
  write.table(significantResults %>% mutate(across(where(is.numeric), ~ round(., digits = params$output_digits))),
              file = file.path(output_folder,
                               paste0(prefix, "-DESeq_output_significant.txt")),
              quote = FALSE, sep = "\t", col.names = NA)
  write.table(summaryTable %>% mutate(across(where(is.numeric), ~ round(., digits = params$output_digits))),
              file = file.path(output_folder,
                               paste0(prefix, "-DESeq_output_all_genes.txt")),
              quote = FALSE, sep = "\t", col.names = NA)
  write.table(CPMddsDF %>% mutate(across(where(is.numeric), ~ round(., digits = params$output_digits))),
              file = file.path(output_folder,
                               paste0(prefix, "-Per_sample_CPM.txt")),
              quote = FALSE, sep = "\t", col.names = NA)
  write.table(Counts %>% as.data.frame() %>% mutate(across(where(is.numeric), ~ round(., digits = params$output_digits))),
              file = file.path(output_folder,
                               paste0(prefix, "-Per_sample_normalized_counts.txt")),
              quote = FALSE, sep = "\t", col.names = NA)



  ##########################
  ### Write results in Excel
  ##########################
  message("write results tables to xls")
  ### Global options
  options("openxlsx.borderColour" = "#4F80BD")
  options("openxlsx.borderStyle" = "thin")
  options("openxlsx.maxWidth" = 50)
  hs1 <- createStyle(textDecoration = "Bold",
                     border = "Bottom",
                     fontColour = "black")
  hs2 <- createStyle(textDecoration = "Bold",
                     border = c("top", "bottom", "left", "right"),
                     fontColour = "black",
                     fgFill = "#C5D9F1")

  ### Summary results - one gene per line, columns are contrasts
  wb1 <- createWorkbook()
  modifyBaseFont(wb1, fontSize = 10, fontName = "Arial Narrow")
  addWorksheet(wb1, "DESeq_results_per_gene")
  for (j in 1:length(contrastsInSummary)) {
    myStartcol = 8 + ((j - 1) * 3)
    myEndcol = 10 + ((j - 1) * 3)
    openxlsx::mergeCells(wb1,
               sheet = 1,
               cols = myStartcol:myEndcol,
               rows = 1)
    openxlsx::writeData(
      wb1,
      sheet = 1,
      x = contrastsInSummary[j],
      startCol = myStartcol,
      startRow = 1)
  }
  openxlsx::conditionalFormatting(wb1,
                        sheet = 1,
                        rows = 1,
                        cols = 1:ncol(summaryTable),
                        type = "contains",
                        rule = "",
                        style = hs2)
  freezePane(wb1, sheet = 1, firstActiveRow = 3, firstActiveCol = 4)
  writeDataTable(wb1,
                 sheet = 1,
                 startRow = 2,
                 x = summaryTable,
                 colNames = TRUE,
                 rowNames = FALSE,
                 tableStyle = "none",
                 headerStyle = hs1,
                 keepNA = TRUE,
                 na.string = "NA")
  setColWidths(wb1, sheet = 1, cols = 1:6, widths = "auto") # This is hard-coded, so prone to error; will only impact auto adjustment of col widths.
  setColWidths(wb1, sheet = 1, cols = 7:ncol(summaryTable), widths = 13) # This is hard-coded, so prone to error; will only impact auto adjustment of col widths.
  fname1 <- file.path(output_folder, paste0("1.", prefix, "-DESeq_by_gene.xlsx"))
  saveWorkbook(wb1, fname1, overwrite = TRUE)
  
  ### All results in one table
  wb2 <- createWorkbook()
  modifyBaseFont(wb2, fontSize = 10, fontName = "Arial Narrow")
  addWorksheet(wb2, paste0("FDR", params$alpha, ".Linear.FC.", params$linear_fc_filter))
  freezePane(wb2, sheet = 1, firstRow = TRUE, firstActiveCol = 4)
  writeDataTable(wb2,
                 sheet = 1,
                 x = significantResults,
                 colNames = TRUE,
                 rowNames = FALSE,
                 tableStyle = "none",
                 headerStyle = hs1,
                 keepNA = TRUE,
                 na.string = "NA")
  setColWidths(wb2, sheet = 1, cols = 1:ncol(significantResults), widths = "auto")
  addWorksheet(wb2, "DESeq_all_results")
  freezePane(wb2, sheet = 2, firstRow = TRUE, firstActiveCol = 4)
  writeDataTable(wb2,
                 sheet = 2,
                 x = allResults,
                 colNames = TRUE,
                 rowNames = FALSE,
                 tableStyle = "none",
                 headerStyle = hs1,
                 keepNA = TRUE,
                 na.string = "NA")
  setColWidths(wb2, sheet = 2, cols = 1:ncol(allResults), widths = "auto")
  fname2 <- file.path(output_folder, paste0("2.", prefix, "-DESeq_all.xlsx"))
  saveWorkbook(wb2, fname2, overwrite = TRUE)

  ### All results with different tabs for each contrast
  wb3 <- createWorkbook()
  modifyBaseFont(wb3, fontSize = 10, fontName = "Arial Narrow")

  short_contrast_names <- paste(exp_contrasts$V1, "v.", exp_contrasts$V2)
  short_contrast_names <- stringr::str_trunc(short_contrast_names,
                                             31,
                                             side = "right",
                                             ellipsis = "")
  # Get rid of illegal characters... I'm sure there will be more some day
  short_contrast_names <- gsub(pattern = ":",
                               replacement = ".",
                               x = short_contrast_names)

  for (i in 1:length(levels(factor(allResults$contrast)))) {
    print(i)
    dataToWrite <- allResults[allResults$contrast == levels(factor(allResults$contrast))[i],]
    addWorksheet(wb3, short_contrast_names[i])
    freezePane(wb3, sheet = i, firstRow = TRUE, firstActiveCol = 4)
    writeDataTable(wb3,
                   sheet = i,
                   x = dataToWrite,
                   colNames = TRUE,
                   rowNames = FALSE,
                   tableStyle = "none",
                   headerStyle = hs1,
                   keepNA = TRUE,
                   na.string = "NA")
    setColWidths(wb3, sheet = i, cols = 1:ncol(dataToWrite), widths = "auto")
  }
  fname3 <- file.path(output_folder, paste0("3.", prefix, "-DESeq_by_contrast.xlsx"))
  saveWorkbook(wb3, fname3, overwrite = TRUE)

  ### CPM
  wb4 <- createWorkbook()
  modifyBaseFont(wb4, fontSize = 10, fontName = "Arial Narrow")
  addWorksheet(wb4, "Counts per million")
  freezePane(wb4, sheet = 1, firstRow = TRUE, firstActiveCol = 4)
  writeDataTable(wb4,
                 sheet = 1,
                 x = as.data.frame(CPMddsDF),
                 colNames = TRUE,
                 rowNames = FALSE,
                 tableStyle = "none",
                 headerStyle = hs1,
                 keepNA = TRUE,
                 na.string = "NA")
  setColWidths(wb4, sheet = 1, cols = 1:ncol(CPMddsDF), widths = "auto")
  fname4 <- file.path(output_folder, paste0("4.", prefix, "-CPM.xlsx"))
  saveWorkbook(wb4, fname4, overwrite = TRUE)

  ### IPA
  IPA <- allResults %>% 
    distinct() %>% 
    pivot_wider(names_from = contrast, values_from = c(log2FoldChange, linearFoldChange, lfcSE, pvalue, padj))
  wb5 <- createWorkbook()
  modifyBaseFont(wb5, fontSize = 10, fontName = "Arial Narrow")
  addWorksheet(wb5, "For IPA upload")
  freezePane(wb5, sheet = 1, firstRow = TRUE, firstActiveCol = 5)
  writeDataTable(wb5,
                 sheet = 1,
                 x = as.data.frame(IPA),
                 colNames = TRUE,
                 rowNames = FALSE,
                 tableStyle = "none",
                 headerStyle = hs1,
                 keepNA = TRUE,
                 na.string = "NA")
  setColWidths(wb5, sheet = 1, cols = 1:ncol(IPA), widths = "auto")
  fname5 <- file.path(output_folder, paste0("5.", prefix, "-IPA.xlsx"))
  saveWorkbook(wb5, fname5, overwrite = TRUE)
  DBI::dbDisconnect(dbconn(db))
}
