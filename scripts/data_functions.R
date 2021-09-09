
filter_metadata <- function(DESeqDesign, params){
    # exclude samples
    if (any(!is.na(params$exclude_samples))) {
        DESeqDesign <- DESeqDesign %>% 
            dplyr::filter(!original_names %in% params$exclude_samples)
    }
    # exclude groups
    if (any(!is.na(params$exclude_groups))) {
        DESeqDesign <- DESeqDesign %>%
            dplyr::filter(!(!!sym(params$design)) %in% params$exclude_groups)
        contrasts_to_filter <- DESeqDesign %>% 
            dplyr::filter(!(!!sym(params$design)) %in% params$exclude_groups) %>%
            pull(params$design) %>% 
            unique()
        contrasts <- contrasts %>%
            dplyr::filter(V1 %in% contrasts_to_filter)
        if (params$strict_contrasts == T) {
            contrasts <- contrasts %>%
                dplyr::filter(V2 %in% contrasts_to_filter)
        }
    }
    if (!is.na(params$include_only_column) & !is.na(params$include_only_group)) {
        DESeqDesign <- DESeqDesign %>%
            dplyr::filter((!!sym(params$include_only_column)) %in% params$include_only_group)
        limit_contrasts <- DESeqDesign %>%
            pull(!!sym(params$design)) %>%
            unique() %>%
            as.character()
        contrasts <- contrasts %>% dplyr::filter(V1 %in% limit_contrasts)
    }
    return(DESeqDesign)
}

format_and_sort_metadata <- function(DESeqDesign, intgroup){
    # Intgroups need to be factors for DESeq2
    # make sure the levels are sorted for plotting later
    for(int in intgroup){
        DESeqDesign[int] <- factor(DESeqDesign[[int]], levels=mixedsort(unique(DESeqDesign[[int]])))
    }
    # if sortcol is defined, sort the design variable based on that
    if (!is.na(params$sortcol)){
        design_factor_reordered <- factor(DESeqDesign[[params$design]],
                                        levels = unique(DESeqDesign[[params$design]][mixedorder(DESeqDesign[[params$sortcol]])]),
                                        ordered = FALSE)
        DESeqDesign[[params$design]] <- design_factor_reordered
    }
    return(DESeqDesign)
}

process_data <- function(sampledata, DESeqDesign, intgroup, params){
  sampleData <- filter_data(sampleData, DESeqDesign, params$threshold)
  DESeqDesign <- format_and_sort_metadata(DESeqDesign, intgroup)
  # need to fix this still
  #check_data(sampleData, DESeqDesign)
  return(list(sampleData=sampleData, DESeqDesign=DESeqDesign))
}

filter_data <- function(sampleData, DESeqDesign, threshold){
    # First data clean-up: replace NA & remove samples with total readcount < threshold
    sampleData[ is.na(sampleData) ] <- 0 
    sampleData <- sampleData[,(colSums(sampleData) > threshold)] # reads required per sample
    DESeqDesign <- DESeqDesign[DESeqDesign$original_names %in% colnames(sampleData),]
    sampleData <- sampleData[,DESeqDesign$original_names]
    return(sampleData)
}

#TODO: make this more better. Should run actual tests on this to make sure that stuff matches
check_data <- function(sampleData, DESeqDesign){
  # Sanity check: each sample (row) in the metadata should have a corresponding column in the count data
  metadata_in_sampledata <- all(DESeqDesign$original_names %in% colnames(sampleData))
  # Sanity check: each column in the count data should have a corresponding sample (row) in the metadata
  sampledata_in_metadata <- all(colnames(sampleData) %in% DESeqDesign$original_names)
  # Find samples that were removed because they weren't in metadata
  removed <- colnames(sampleData[which(!colnames(sampleData) %in% DESeqDesign$original_names)])
  # Reorder the metadata table to correspond to the order of columns in the count data
  DESeqDesign <- DESeqDesign[DESeqDesign$original_names %in% colnames(sampleData),]
  # DESeqDesign <- na.omit(DESeqDesign) # This can cause issues since it removes any lines with missing data. Should instead check for NAs in required columns.
  sampleData <- sampleData[,DESeqDesign$original_names]
  samples_after <- nrow(DESeqDesign)

  head(DESeqDesign$original_names)
  head(colnames(sampleData)) # Output should match
  return(sampleData)
}


load_count_data <- function(SampleDataFile, sampledata_sep){
  sampleData <- read.delim(SampleDataFile,
                         sep = sampledata_sep,
                         stringsAsFactors = FALSE,
                         header = TRUE, 
                         quote = "\"",
                         row.names = 1,
                         check.names = FALSE)
  return(sampleData)
}