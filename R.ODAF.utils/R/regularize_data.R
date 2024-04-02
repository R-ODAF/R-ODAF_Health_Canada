#' Regularize Data for Visualization and Further Analysis
#'
#' Applies variance stabilizing transformation (VST) or removes batch effects from the data.
#'
#' @param dds A DESeq2DataSet object.
#' @param design A character string of the design variable.
#' @param covariates A character vector of covariate factors. These will be excluded from batch effect removal.
#' @param batch_param A character string specifying the batch parameter: for example, "batch" or "plate".
#' @param blind Logical, if TRUE variance stabilizing transformations are blind to experimental design.
#' @return A DESeqTransform object with regularized data.
#' @importFrom DESeq2 vst varianceStabilizingTransformation
#' @importFrom SummarizedExperiment assay
#' @importFrom stats model.matrix
#' @importFrom limma removeBatchEffect
#' @export
regularize_data <- function(dds,
                            design,
                            covariates = NA,
                            batch_param,
                            blind = FALSE) {
  if (!is.na(batch_param)) {
    rld <- DESeq2::vst(dds, blind)
    mat <- SummarizedExperiment::assay(rld)
    if (!is.na(covariates)) {
      condition <- formula(
        paste0("~",
               design,
               paste0(covariates[!covariates %in% batch_param], collapse = " + "))
      )

    } else {
      condition <- formula(paste0("~", design))
    }
    mm <- stats::model.matrix(condition, SummarizedExperiment::colData(rld))
    if (length(unique(rld[[batch_param]])) > 1) {
      mat <- limma::removeBatchEffect(
        mat,
        batch = rld[[batch_param]],
        design = mm
      )
      SummarizedExperiment::assay(rld) <- mat
    }
  } else {
    if (nrow(dds) < 1000) {
      rld <- varianceStabilizingTransformation(dds, blind)
    } else {
      rld <- vst(dds, blind)
    }
  }
  return(rld)
}
