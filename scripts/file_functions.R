set_up_paths <- function(params) {
    paths <- list()
    # Other important system paths to specify in config
    paths$wikipathways <- params$wikipathways_directory
    # For project structure
    # Should probably update this to use the file.path() function.
    paths$root <- params$projectdir
    paths$data <- file.path(paths$root, "data")
    paths$raw <- file.path(paths$data, "raw")
    paths$processed <- file.path(paths$data, "processed")
    paths$metadata <- file.path(paths$data, "metadata")
    paths$results <- file.path(paths$root, "analysis")
    paths$reports <- file.path(paths$results, "reports")
    paths$BMD_output <- file.path(paths$results, "BMD_and_biomarker_files")
    paths$RData <- file.path(paths$results, "DEG_RData")
    paths$pathway_analysis <- file.path(paths$results, "pathway_analysis")
    lapply(paths, function(x) if(!dir.exists(x)) dir.create(x, recursive = TRUE))
    return(paths)
}

set_up_paths_2 <- function(paths, params, facets){
  if (is.na(params$group_facet) || is.null(params$group_facet)) {
    paths$DEG_output <- file.path(paths$results, "DEG_lists")
    paths$pathway_analysis <- file.path(paths$results, "pathway_analysis")
  } else {
    # make multiple outputs for different facets
    for(f in facets){
      paths$DEG_output[[f]] <- file.path(paths$results, "DEG_lists", paste0(f))
      paths$pathway_analysis[[f]] <- file.path(paths$results, "pathway_analysis", paste0(f))
    }
  }
  lapply(paths$DEG_output, function(x) if(!dir.exists(x)) dir.create(x, recursive = TRUE))
  lapply(paths$pathway_analysis, function(x) if(!dir.exists(x)) dir.create(x, recursive = TRUE))
  return(paths)
}
  
  

load_cached_data <- function(RDataPath, params, sampleData, facets=NULL){
    if(!is.na(params$group_facet)){
        ddsList = list()
        for (current_filter in facets) {
            dds <- readRDS(file = file.path(RDataPath, paste0("dds_", paste(current_filter, collapse = "_"), ".RData")))
            if (!identical(as.data.frame(round(counts(dds))), round(sampleData), 0)) {
                stop("Attempted to load a cached file that contained non-identical count data, exiting")
            }
            ddsList[[current_filter]] <- dds
        }
        return(ddsList)
  } else {
    if (file.exists(file.path(RDataPath, "dds.RData")) ) {
        print(paste("Already found DESeq2 object from previous run; loading from disk."))
        dds <- readRDS(file.path(RDataPath, "dds.RData"))
        if (!identical(as.data.frame(round(counts(dds))), round(sampleData), 0)) {
            stop("Attempted to load a cached file that contained non-identical count data, exiting")
        }
    }
    return(dds)
  }
}

save_cached_data <- function(dds, RDataPath, current_filter=NULL){
    if (is.na(current_filter)) {
        saveRDS(dds, file = file.path(RDataPath, "dds.RData"))
    } else {
        saveRDS(dds, file = file.path(RDataPath, paste0("dds_", paste(current_filter, collapse = "_"), ".RData")))
    }
}
