filter_data <- function(sampleData, DESeqDesign, threshold){
    # First data clean-up: replace NA & remove samples with total readcount < threshold
    sampleData[ is.na(sampleData) ] <- 0 
    sampleData <- sampleData[,(colSums(sampleData) > threshold)] # reads required per sample

    DESeqDesign <- DESeqDesign[DESeqDesign$original_names %in% colnames(sampleData),]
    sampleData <- sampleData[,DESeqDesign$original_names]
}

#TODO: make this more better
check_data <- function(sampleData){
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


get_design <- function(intgroup){
  return(formula(paste0("~", paste0(intgroup, collapse = " + ")))
}


learn_deseq_model <- function(sampledata, DESeqDesign, intgroup, params){
  current_design <- get_design(intgroup)
  dds <- DESeqDataSetFromMatrix(countData = round(sampleData),
                                colData   = as.data.frame(DESeqDesign),
                                design    = current_design)
  # TODO: what is this filtering for?
  #dds <- dds[, DESeqDesign$original_names]
  #dds <- dds[rowSums(counts(dds)) > 1]
  bpparam <- MulticoreParam(params$cpus)
  dds <- DESeq(dds, parallel = TRUE, BPPARAM = bpparam)
  return(dds)
}

# covariates are used to calculate within-group variability. Nuisance parameters are removed (e.g. for visualization) See:
# http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-after-vst-are-there-still-batches-in-the-pca-plot
regularize_data <- function(dds, covariates, nuisance, blind=FALSE){
  if (!is.na(nuisance)) {
    rld <- vst(dds, blind)
    mat <- assay(rld)
    condition <- formula(paste0("~",
                                paste0(covariates[!covariates %in% nuisance],collapse = " + ")))
    mm <- model.matrix(condition, colData(rld))
    mat <- limma::removeBatchEffect(mat, batch = rld[[nuisance]], design = mm)
    assay(rld) <- mat
  } else {
    rld <- vst(dds, blind)
  }
  return(rld)
}