#' Render an R Markdown Report
#'
#' This function renders an R Markdown report to the specified output file,
#' using the provided parameters. It generates a random temporary directory
#' for storing intermediate files during the rendering process and cleans up
#' afterwards.
#'
#' @param report_in The path to the R Markdown source file.
#' @param report_out The path where the rendered report should be saved.
#' @param render_pars A list of parameters to pass to the R Markdown document.
#' @return Invisible NULL. The function is called for its side effect of rendering a report.
#' @importFrom stringi stri_rand_strings
#' @export
render_report <- function(report_in, report_out, render_pars) {
  random_tmp <- file.path("/tmp", paste0("intermediates_", create_random_string()))
  rmarkdown::render(input = report_in,
    encoding = "UTF-8",
    output_file = report_out,
    # output_dir = random_tmp,
    params = render_pars,
    envir = new.env(),
    clean = TRUE,
    run_pandoc = TRUE,
    intermediates_dir = random_tmp)
  system(paste0("rm -rf ", random_tmp))
}


#' Create a random string.
#'
#' This function generates a random string by sampling characters from the
#' combination of uppercase letters and digits.
#'
#' @param digits A vector of digits to include in the random string.
#'               Default is 0:9.
#' @param letters_set A vector of letters to include in the random string.
#'
#' @return A random string of length 10.
#'
#' @export
create_random_string <- function(len = 20, digits = 0:9, letters_set = c(LETTERS, letters)) {
  v <- sample(c(letters_set, digits), len, replace = TRUE)
  return(paste0(v, collapse = ""))
}
