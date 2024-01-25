library(gtools)
library(data.table)

# filter metadata
# this only applies to specifically included/excluded data, not the facet filtering
filter_metadata <- function(exp_metadata, params, design){
    # exclude samples
    if (any(!is.na(params$exclude_samples))) {
        exp_metadata <- exp_metadata %>% 
            dplyr::filter(!original_names %in% params$exclude_samples)
    }
    # exclude groups
    if (any(!is.na(params$exclude_groups))) {
        exp_metadata <- exp_metadata %>%
            dplyr::filter(!(!!sym(design)) %in% params$exclude_groups)
        contrasts_to_filter <- exp_metadata %>% 
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
        exp_metadata <- exp_metadata %>%
            dplyr::filter((!!sym(params$include_only_column)) %in% params$include_only_group)
        limit_contrasts <- exp_metadata %>%
            pull(!!sym(design)) %>%
            unique() %>%
            as.character()
        contrasts <- contrasts %>% dplyr::filter(V1 %in% limit_contrasts)
    }
    return(exp_metadata)
}

format_and_sort_metadata <- function(exp_metadata, intgroup, design, sortcol){
    # Intgroups need to be factors for DESeq2
    # make sure the levels are sorted for plotting later
    for(int in intgroup){
        exp_metadata[int] <- factor(exp_metadata[[int]], levels=mixedsort(unique(exp_metadata[[int]])))
    }
    # if sortcol is defined, sort the design variable based on that
    if (!is.na(sortcol)){
        design_factor_reordered <- factor(exp_metadata[[design]],
                                        levels = unique(exp_metadata[[design]][mixedorder(exp_metadata[[sortcol]])]),
                                        ordered = FALSE)
        exp_metadata[[design]] <- design_factor_reordered
        
        # also sort the other interesting groups we want to plot
        for (ig in params$intgroup_to_plot){
          intgroup_reordered <- factor(exp_metadata[[ig]],
                                            levels = unique(exp_metadata[[ig]][mixedorder(exp_metadata[[sortcol]])]),
                                            ordered = FALSE)
          exp_metadata[[ig]] <- intgroup_reordered
        }
    }
    return(exp_metadata)
}

sort_contrasts <- function(exp_metadata, contrasts, design, sortcol){
    ordered_design <- exp_metadata[mixedorder(exp_metadata[,params$sortcol]),] %>%
        dplyr::select(design) %>%
        dplyr::pull()
    ordered_contrasts <- contrasts %>%
        dplyr::slice(match(ordered_design, V1)) %>%
        unique()
    return(ordered_contrasts)
}

process_data_and_metadata <- function(count_data, exp_metadata, contrasts, intgroup, design, params){
    exp_metadata <- filter_metadata(exp_metadata, params, design)
    exp_metadata <- format_and_sort_metadata(exp_metadata, intgroup, design, params$sortcol)
    count_data <- filter_data(count_data, exp_metadata, params$nmr_threshold)
    if(!is.na(params$sortcol)){
        contrasts <- sort_contrasts(exp_metadata, contrasts, design, params$sortcol)
    }
    check_data(count_data, exp_metadata, contrasts)
    return(list(count_data=count_data, exp_metadata=exp_metadata, contrasts=contrasts))
}

# filter count data based on minimum counts. Does not do facet filtering
filter_data <- function(count_data, exp_metadata, threshold){
    # First data clean-up: replace NA & remove samples with total readcount < threshold 
    count_data[ is.na(count_data) ] <- 0 
    count_data <- count_data[,(colSums(count_data) > threshold)] # reads required per sample
    #count_data <- count_data[(rowSums(count_data) > 1),] # reads required per gene
    exp_metadata <- exp_metadata[exp_metadata$original_names %in% colnames(count_data),]
    count_data <- count_data[,exp_metadata$original_names]
    return(count_data)
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


load_count_data <- function(count_data_file, sampledata_sep, collapse_probes = FALSE){
  count_data <- read.delim(count_data_file,
                         sep = sampledata_sep,
                         stringsAsFactors = FALSE,
                         header = TRUE, 
                         quote = "\"",
                         row.names = 1,
                         check.names = FALSE)
  if (collapse_probes == TRUE) {
    count_data$gene <- rownames(count_data)
    count_data$gene <- gsub("_.*", "", count_data$gene)
    data.table::setDT(count_data)
    count_data[, lapply(.SD, sum),
               by = gene, .SDcols = which(sapply(count_data, is.numeric))]
    
    # count_data <- count_data %>%
    #   dplyr::group_by(gene) %>%
    #   dplyr::summarize(across(where(is.numeric), sum))
    #count_data <- as.data.frame(count_data)
    
    rownames(count_data) <- count_data$gene
    count_data <- count_data[,-1]
  }
  return(count_data)
}

# subset metadata based on facet + filter
subset_metadata <- function(exp_metadata, design, contrasts, current_facet, current_filter){
    contrasts_to_filter <- exp_metadata %>%
        dplyr::filter(!!sym(current_facet) %in% current_filter) %>% # NOTE: Not sure if %in% or == is better here.
        pull(design) %>% 
        unique()
    contrasts_subset <- contrasts %>% dplyr::filter(V1 %in% contrasts_to_filter)
    if (params$strict_contrasts == T) {
        contrasts_subset <- contrasts_subset %>% dplyr::filter(V2 %in% contrasts_to_filter)
    }
    exp_metadata_subset <- exp_metadata %>% dplyr::filter(!!sym(design) %in% (unlist(contrasts_subset) %>% unique()) )
    # relevel the design and interesting groups
    exp_metadata_subset[[design]] <- factor(exp_metadata_subset[[design]],
                                            levels = unique(unlist(contrasts_subset)),
                                            ordered = FALSE)
    if (!is.na(params$sortcol)){
        design_factor_reordered <- factor(exp_metadata_subset[[design]],
                                          levels = unique(exp_metadata_subset[[design]][mixedorder(exp_metadata_subset[[params$sortcol]])]),
                                          ordered = FALSE)
        exp_metadata_subset[[design]] <- design_factor_reordered
        
      for (ig in params$intgroup_to_plot){
        intgroup_reordered <- factor(exp_metadata_subset[[ig]],
                                     levels = unique(exp_metadata_subset[[ig]][mixedorder(exp_metadata_subset[[params$sortcol]])]),
                                     ordered = FALSE)
        exp_metadata_subset[[ig]] <- intgroup_reordered
      }
    }
    return(list(exp_metadata=exp_metadata_subset, contrasts=contrasts_subset))
}

# subset count data to samples in metadata
subset_data <- function(count_data, exp_metadata){
    # Reorder the metadata table to correspond to the order of columns in the count data
    exp_metadata_sorted <- exp_metadata[exp_metadata$original_names %in% colnames(count_data),]
    count_data_subset <- count_data[,exp_metadata_sorted$original_names]
    return(count_data_subset)
}

subset_results <- function(res, exp_metadata){
  # Reorder the metadata table to correspond to the order of columns in the count data
  exp_metadata_sorted <- exp_metadata[exp_metadata$original_names %in% colnames(count_data),]
  res_subset <- res[,exp_metadata_sorted$original_names]
  return(res_subset)
  
}