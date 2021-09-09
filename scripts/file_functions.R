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
    paths$reports <- file.path(paths$root, "reports")
    paths$results <- file.path(paths$root, "results")
    if (is.na(params$group_facet)) {
      paths$DEG_output <- file.path(paths$results, "DEG_output")
    } else {
      paths$DEG_output <- file.path(paths$results, "DEG_output", paste0("group_", paste(params$group_filter, collapse = "_")))
    }
    paths$pathway_analysis <- file.path(paths$DEG_output, "/pathway_analysis")
    paths$RData <- file.path(paths$DEG_output, "/RData")
    paths$BMD_output <- file.path(paths$results, "/DEG_output/BMD_and_biomarker_files")
    lapply(paths, function(x) if(!dir.exists(x)) dir.create(x))
    return(paths)
}

load_cached_data <- function(RDataPath, params){
  if(!is.na(params$group_facet)){
    stop("TODO: implement cached data loading for faceted data")
  }
  if (file.exists(file.path(RDataPath, "dds.RData")) ) {
    print(paste("Already found DESeq2 object from previous run; loading from disk."))
    load(file.path(RDataPath, "dds.RData"))
  }
  if (!identical(as.data.frame(round(counts(dds))), round(sampleData), 0)) {
    stop("Attempted to load a cached file that contained non-identical count data, exiting")
  }
}

save_cached_data <- function(dds, RDataPath, params){
  if (is.na(params$group_facet)) {
    save(dds, file = file.path(RDataPath, "dds.RData"))
  } else {
    save(dds, file = file.path(RDataPath, paste0("dds_", paste(params$group_filter, collapse = "_"), ".RData")))
  }
}
