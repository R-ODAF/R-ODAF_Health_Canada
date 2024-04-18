#' Calculate Row Variances
#'
#' This function computes the variance for each row in a numeric matrix or data frame,
#' with support for `NA` values. Rows with one or fewer non-`NA` values will have their
#' variance set to `NA`.
#'
#' @param x A numeric matrix or data frame with rows representing variables and columns representing observations.
#' @param ... Additional arguments to be passed to `rowSums` and `rowMeans`.
#'
#' @return A numeric vector containing the variance for each row.
#' @export
row_vars <- function(x, ...) {
  sqr     <- function(x)  x * x
  n       <- rowSums(!is.na(x))
  n[n <= 1] <- NA
  return(rowSums(sqr(x - rowMeans(x, ...)), ...) / (n - 1))
}
