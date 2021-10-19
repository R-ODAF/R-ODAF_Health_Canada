#!/usr/bin/R
# Custom parameters for the report
library(tidyverse)
require(yaml)

config <- yaml::read_yaml(file.path(here::here(),
                                    "config/config.yml"),
                          eval.expr = T)

projectdir <- config$QC$projectdir
if (is.null(projectdir)) {
  projectdir <- here::here()
  config$QC$projectdir <- projectdir
}

# Input file - Rmd
inputFile <- file.path(projectdir, "Rmd", "Sample_QC.Rmd")

message("Writing QC report for all samples in the experiment.")
# Output file - HTML
filename <- paste0("Study-wide_Sample_QC_",
                   config$QC$platform, "_",
                   config$QC$project_name, "_",
                   format(Sys.time(),'%d-%m-%Y.%H.%M'),
                   ".html")

outFile <- file.path(projectdir,
                     "reports",
                     filename)

dir.create(file.path(projectdir,"reports"))

rmarkdown::render(input = inputFile,
                  encoding = "UTF-8",
                  output_file = outFile,
                  params = config$QC,
                  envir = new.env())
