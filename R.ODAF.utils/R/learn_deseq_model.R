#' Learn a DESeq2 Model from Count Data
#'
#' This function creates a DESeqDataSet and runs the DESeq2 pipeline to learn a model based on count data.
#'
#' @param sd A data frame of count data.
#' @param metadata A data frame with the same number of rows as there are columns in sd, containing the experimental metadata.
#' @param design Character vector specifying design factor (a column in the metadata representing experimental groups).
#' @param params The list of configuration parameters for the experiment.
#' @return A DESeq2DataSet object after running the DESeq pipeline.
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq counts
#' @importFrom BiocParallel MulticoreParam
#' @export
learn_deseq_model <- function(sd, metadata, design, params) {
  current_design <- get_design(design, params)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(sd),
    colData = as.data.frame(metadata),
    design = current_design)
  # if(params$filter_gene_counts){ # filter gene counts to decrease runtime. Not recommended for biomarker input!
  #   dds <- dds[rowSums(DESeq2::counts(dds)) > 1]
  # }
  bpparam <- BiocParallel::MulticoreParam(params$cpus)
  dds <- dds[rowSums(DESeq2::counts(dds)) > 1]
  dds <- DESeq2::DESeq(dds, parallel = TRUE, BPPARAM = bpparam)
  return(dds)
}

#' Generate Design Formula for DESeq2 Analysis
#'
#' This function generates a design formula for use in DESeq2 analysis.
#' If a formula override is provided in the params, that formula is used instead.
#'
#' @param design Character vector specifying design factor(s).
#' @param params A list of parameters containing formula_override.
#' @return A design formula for DESeq2.
#' @export
get_design <- function(design, params) {
  # return(formula(paste0("~", paste0(c(design), collapse = " + "))))
  if (!is.na(params$formula_override)) {
    des <- formula(params$formula_override)
  } else { des <- formula(paste0("~", design)) }
  return(des)
}
