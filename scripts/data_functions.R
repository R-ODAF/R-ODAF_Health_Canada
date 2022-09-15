# filter metadata
# this only applies to specifically included/excluded data, not the facet filtering
filter_metadata <- function(DESeqDesign, params, design){
    # exclude samples
    if (any(!is.na(params$exclude_samples))) {
        DESeqDesign <- DESeqDesign %>% 
            dplyr::filter(!original_names %in% params$exclude_samples)
    }
    # exclude groups
    if (any(!is.na(params$exclude_groups))) {
        DESeqDesign <- DESeqDesign %>%
            dplyr::filter(!(!!sym(design)) %in% params$exclude_groups)
        contrasts_to_filter <- DESeqDesign %>% 
            dplyr::filter(!(!!sym(design)) %in% params$exclude_groups) %>%
            pull(design) %>% 
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
            pull(!!sym(design)) %>%
            unique() %>%
            as.character()
        contrasts <- contrasts %>% dplyr::filter(V1 %in% limit_contrasts)
    }
    return(DESeqDesign)
}

format_and_sort_metadata <- function(DESeqDesign, intgroup, design, sortcol){
    # Intgroups need to be factors for DESeq2
    # make sure the levels are sorted for plotting later
    for(int in intgroup){
        DESeqDesign[int] <- factor(DESeqDesign[[int]], levels=mixedsort(unique(DESeqDesign[[int]])))
    }
    # if sortcol is defined, sort the design variable based on that
    if (!is.na(sortcol)){
        design_factor_reordered <- factor(DESeqDesign[[design]],
                                        levels = unique(DESeqDesign[[design]][mixedorder(DESeqDesign[[sortcol]])]),
                                        ordered = FALSE)
        DESeqDesign[[design]] <- design_factor_reordered
    }
    return(DESeqDesign)
}

sort_contrasts <- function(DESeqDesign, contrasts, design, sortcol){
    ordered_design <- DESeqDesign[mixedorder(DESeqDesign[,params$sortcol]),] %>%
        dplyr::select(design) %>%
        dplyr::pull()
    ordered_contrasts <- contrasts %>%
        dplyr::slice(match(ordered_design, V1)) %>%
        unique()
    return(ordered_contrasts)
}

process_data_and_metadata <- function(sampleData, DESeqDesign, contrasts, intgroup, design, params){
    sampleData <- filter_data(sampleData, DESeqDesign, params$nmr_threshold)
    DESeqDesign <- filter_metadata(DESeqDesign, params, design)
    DESeqDesign <- format_and_sort_metadata(DESeqDesign, intgroup, design, params$sortcol)
    if(!is.na(params$sortcol)){
        contrasts <- sort_contrasts(DESeqDesign, contrasts, design, params$sortcol)
    }
    check_data(sampleData, DESeqDesign, contrasts)
    return(list(sampleData=sampleData, DESeqDesign=DESeqDesign, contrasts=contrasts))
}

# filter count data based on minimum counts. Does not do facet filtering
filter_data <- function(sampleData, DESeqDesign, threshold){
    # First data clean-up: replace NA & remove samples with total readcount < threshold 
    sampleData[ is.na(sampleData) ] <- 0 
    sampleData <- sampleData[,(colSums(sampleData) > threshold)] # reads required per sample
    #sampleData <- sampleData[(rowSums(sampleData) > 1),] # reads required per gene
    DESeqDesign <- DESeqDesign[DESeqDesign$original_names %in% colnames(sampleData),]
    sampleData <- sampleData[,DESeqDesign$original_names]
    return(sampleData)
}

# sanity checks
check_data <- function(sd, des, con){
    message("Sanity checks for data")
    # make sure they're not empty
    stopifnot(exprs = {
        ncol(sd) > 0
        nrow(des) > 0
        nrow(con) > 0
    })
    # Sanity check: each sample (row) in the metadata should have a corresponding column in the count data
    stopifnot(all(des$original_names %in% colnames(sd)))
    # Sanity check: each column in the count data should have a corresponding sample (row) in the metadata
    stopifnot(all(colnames(sd) %in% des$original_names))
    message("All OK 👍")
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
subset_metadata <- function(DESeqDesign, design, contrasts, current_facet, current_filter){
    contrasts_to_filter <- DESeqDesign %>%
        dplyr::filter(!!sym(current_facet) %in% current_filter) %>% # NOTE: Not sure if %in% or == is better here.
        pull(design) %>% 
        unique()
    contrasts_subset <- contrasts %>% dplyr::filter(V1 %in% contrasts_to_filter)
    if (params$strict_contrasts == T) {
        contrasts_subset <- contrasts_subset %>% dplyr::filter(V2 %in% contrasts_to_filter)
    }
    DESeqDesign_subset <- DESeqDesign %>% dplyr::filter(!!sym(design) %in% (unlist(contrasts_subset) %>% unique()) )
    #levels(DESeqDesign_subset[design,]) <- unlist(contrasts_subset) %>% unique()
    return(list(DESeqDesign=DESeqDesign_subset, contrasts=contrasts_subset))
}

# subset count data to samples in metadata
subset_data <- function(sampleData, DESeqDesign){
    # Reorder the metadata table to correspond to the order of columns in the count data
    DESeqDesign_sorted <- DESeqDesign[DESeqDesign$original_names %in% colnames(sampleData),]
    sampleData_subset <- sampleData[,DESeqDesign_sorted$original_names]
    return(sampleData_subset)
}

subset_results <- function(res, DESeqDesign){
  # Reorder the metadata table to correspond to the order of columns in the count data
  DESeqDesign_sorted <- DESeqDesign[DESeqDesign$original_names %in% colnames(sampleData),]
  res_subset <- res[,DESeqDesign_sorted$original_names]
  return(res_subset)
  
}