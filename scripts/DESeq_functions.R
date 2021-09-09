get_design <- function(intgroup){
  return(formula(paste0("~", paste0(intgroup, collapse = " + "))))
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