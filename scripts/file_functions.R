library(edgeR)


set_up_paths <- function(params) {
    paths <- list()
    # Other important system paths to specify in config
    paths$wikipathways <- params$wikipathways_directory
    # For project structure
    paths$root <- params$projectdir
    paths$inputs <- file.path(paths$root, "inputs")
    paths$output <- file.path(paths$root, "output")
    paths$raw <- file.path(paths$inputs, "raw")
    paths$processed <- file.path(paths$output, "processed")
    paths$metadata <- file.path(paths$inputs, "metadata")
    paths$contrasts <- file.path(paths$inputs, "contrasts")
    paths$results <- file.path(paths$output, "analysis", results_location_arg)
    paths$BMD_output <- file.path(paths$results, "BMD_and_biomarker_files")
    paths$RData <- file.path(paths$results, "RData")
    paths$pathway_analysis <- file.path(paths$results, "pathway_analysis")
    lapply(paths, function(x) if(!dir.exists(x)) dir.create(x, recursive = TRUE))
    return(paths)
}

set_up_paths_2 <- function(paths, params, facets, display_facets){
  if (is.na(params$group_facet) || is.null(params$group_facet)) {
    paths$DEG_output <- file.path(paths$results, "DEG_lists")
    paths$biosets_output <- paths$BMD_output
  } else {
    # make multiple outputs for different facets
    for(f in facets){
      paths$DEG_output[[f]] <- file.path(paths$results, "DEG_lists", paste0(f))
      paths$biosets_output[[f]] <- file.path(paths$BMD_output, paste0(f))
    }
  }
  lapply(paths$DEG_output, function(x) if(!dir.exists(x)) dir.create(x, recursive = TRUE))
  lapply(paths$biosets_output, function(x) if(!dir.exists(x)) dir.create(x, recursive = TRUE))
  return(paths)
}

set_up_paths_3 <- function(paths,params,display_facets){
  if (is.na(params$display_group_facet) || is.null(params$display_group_facet)) {
    paths$pathway_analysis <- file.path(paths$results, "pathway_analysis")
  } else {
    # make multiple outputs for different facets
    paths$pathway_analysis <- c()
    for(f in display_facets){
      paths$pathway_analysis[[f]] <- file.path(paths$results, "pathway_analysis", paste0(f))
    }
  }
  lapply(paths$pathway_analysis, function(x) if(!dir.exists(x)) dir.create(x, recursive = TRUE))
  return(paths)
}

  
  

load_cached_data <- function(RDataPath, params, count_data, facets=NULL){
    if(!is.na(params$group_facet)){
        ddsList = list()
        for (current_filter in facets) {
            dds <- readRDS(file = file.path(RDataPath, paste0("dds_", paste(current_filter, collapse = "_"), ".RData")))
            if (!identical(as.data.frame(round(counts(dds))), round(count_data), 0)) {
                stop("Attempted to load a cached file that contained non-identical count data, exiting")
            }
            ddsList[[current_filter]] <- dds
        }
        return(ddsList)
  } else {
    if (file.exists(file.path(RDataPath, "dds.RData")) ) {
        print(paste("Already found DESeq2 object from previous run; loading from disk."))
        dds <- readRDS(file.path(RDataPath, "dds.RData"))
        if (!identical(as.data.frame(round(counts(dds))), round(count_data), 0)) {
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


write_additional_output <- function(count_data, exp_metadata, design_to_use, params){
  if (!is.na(params$dose)) {
    cpm_data <- cpm(count_data)
    bmdexpress <- as.data.frame(log2(cpm_data + 1))
    bmdexpress <- cbind(SampleID = c(row.names(bmdexpress)),
                        bmdexpress,
                        stringsAsFactors = F)
    
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
    
    bmdexpress <- bmdexpress[rowSums(count_data) > 5,]
    # add a dose header line to both files
    bmdexpress <- rbind(c("Dose", as.character(exp_metadata[colnames(bmdexpress)[-1],][[params$dose]])),
                        bmdexpress,
                        stringsAsFactors = F)
    biomarkers <- rbind(c("Dose", as.character(exp_metadata[colnames(biomarkers)[-1],][[params$dose]])),
                        biomarkers,
                        stringsAsFactors = F)
    
    # Determine names of dose groups in which n per group > 1
    groups_for_bmdexpress <- which(table(t(bmdexpress[1,])) > 1) %>% names()
    # Rewrite bmdexpress table
    # Manually include the "Dose" and gene name column
    bmdexpress <- bmdexpress[,(bmdexpress[1,]) %in% c("Dose",groups_for_bmdexpress)]
    
    if (!is.na(params$group_facet)) {
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
                  quote = F,
                  sep = "\t",
                  row.names = F,
                  col.names = T)
      write.table(biomarkers,
                  file = file.path(paths$BMD_output,
                                   fname2),
                  quote = F,
                  sep = "\t",
                  row.names = F,
                  col.names = T)
    } else {
      write.table(bmdexpress,
                  file = file.path(paths$BMD_output, "bmdexpress_input.txt"),
                  quote = F,
                  sep = "\t",
                  row.names = F,
                  col.names = T)
      write.table(biomarkers,
                  file = file.path(paths$BMD_output, "biomarkers_input.txt"),
                  quote = F,
                  sep = "\t",
                  row.names = F,
                  col.names = T)
    }
  }
  
}
