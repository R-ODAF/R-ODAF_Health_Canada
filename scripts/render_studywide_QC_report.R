#!/usr/bin/R
# Custom parameters for the report
library(tidyverse)
require(yaml)

#Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio/bin/pandoc")

config <- yaml::read_yaml(file.path(here::here(),
                                    "config/config.yaml"),
                          eval.expr = T)

# Combine required params from config
params <-c(config$common, config$QC)
projectdir <- params$projectdir
# If projectdir is not set, figure out current project root directory
if (is.null(projectdir)) {
  projectdir <- here::here()
  params$projectdir <- projectdir
}

# Replace any other NULL in params with NA
replace_nulls <- function(x) {ifelse(is.null(x), NA, x)}
params <- lapply(params, replace_nulls)

# Input file - Rmd
inputFile <- file.path(projectdir, "Rmd", "Sample_QC.Rmd")

message("Writing QC report for all samples in the experiment.")
# Output file - HTML
filename <- paste0("Study-wide_Sample_QC_",
                   params$platform, "_",
                   params$project_name, "_",
                   format(Sys.time(),'%d-%m-%Y.%H.%M'),
                   ".html")

outFile <- file.path(projectdir,
                     "analysis", "QC",
                     filename)

rmarkdown::render(input = inputFile,
                  encoding = "UTF-8",
                  output_file = outFile,
                  params = params,
                  envir = new.env())
