# filter metadata
# this only applies to specifically included/excluded data, not the facet filtering
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

sort_contrasts <- function(DESeqDesign, contrasts, params){
    ordered_design <- DESeqDesign[mixedorder(DESeqDesign[,params$sortcol]),] %>%
        dplyr::select(params$design) %>%
        dplyr::pull()
    ordered_contrasts <- contrasts %>%
        dplyr::slice(match(ordered_design, V1)) %>%
        unique()
    return(ordered_contrasts)
}

process_data_and_metadata <- function(sampledata, DESeqDesign, contrasts, intgroup, params){
    sampleData <- filter_data(sampleData, DESeqDesign, params$threshold)
    DESeqDesign <- filter_metadata(DESeqDesign, params)
    DESeqDesign <- format_and_sort_metadata(DESeqDesign, intgroup)
    if(!is.na(params$sortcol)){
        contrasts <- sort_contrasts(DESeqDesign, contrasts, params)
    }
    check_data(sampleData, DESeqDesign, contrasts)
    return(list(sampleData=sampleData, DESeqDesign=DESeqDesign, contrasts=contrasts))
}

# filter count data based on minimum counts. Does not do facet filtering
filter_data <- function(sampleData, DESeqDesign, threshold){
    # First data clean-up: replace NA & remove samples with total readcount < threshold
    sampleData[ is.na(sampleData) ] <- 0 
    sampleData <- sampleData[,(colSums(sampleData) > threshold)] # reads required per sample
    DESeqDesign <- DESeqDesign[DESeqDesign$original_names %in% colnames(sampleData),]
    sampleData <- sampleData[,DESeqDesign$original_names]
    return(sampleData)
}

# sanity checks
check_data <- function(sampleData, DESeqDesign, contrasts){
    message("Sanity checks for data")
    # make sure they're not empty
    stopifnot(exprs = {
        nrow(sampleData) > 0
        nrow(DESeqDesign) > 0
        nrow(contrasts) > 0
    })
  # Sanity check: each sample (row) in the metadata should have a corresponding column in the count data
  stopifnot(all(DESeqDesign$original_names %in% colnames(sampleData)))
  # Sanity check: each column in the count data should have a corresponding sample (row) in the metadata
  stopifnot(all(colnames(sampleData) %in% DESeqDesign$original_names))
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

# subset metadata based on facet + filter
subset_metadata <- function(DESeqDesign, params, contrasts){
    contrasts_to_filter <- DESeqDesign %>%
        dplyr::filter(!!sym(params$group_facet) %in% params$group_filter) %>% # NOTE: Not sure if %in% or == is better here.
        pull(params$design) %>% 
        unique()
    contrasts_subset <- contrasts %>% dplyr::filter(V1 %in% contrasts_to_filter)
    if (params$strict_contrasts == T) {
        contrasts_subset <- contrasts_subset %>% dplyr::filter(V2 %in% contrasts_to_filter)
    }
    DESeqDesign_subset <- DESeqDesign %>% dplyr::filter(!!sym(params$design) %in% (unlist(contrasts_subset) %>% unique()) )
    return(list(DESeqDesign=DESeqDesign_subset, contrasts=contrasts_subset))
}

# subset count data to samples in metadata
# TODO: add some tests here. Maybe test for zero samples left or something
# TODO: also should be a check earlier on, separate from this, that checks that all sample names are in metadata and vice versa to start with
# TODO: finally, there should be a check after all the processing is done that the actual count data hasn't changed
subset_data <- function(sampleData, DESeqDesign){
    # Reorder the metadata table to correspond to the order of columns in the count data
    DESeqDesign_sorted <- DESeqDesign[DESeqDesign$original_names %in% colnames(sampleData),]
    sampleData_subset <- sampleData[,DESeqDesign_sorted$original_names]
    return(sampleData_subset)
}