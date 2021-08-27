#!/usr/bin/R
# Custom parameters for the report
library(tidyverse)
require(yaml)

config <- yaml::read_yaml(file.path(here::here(),
                                    "Rmd/config.qc.yml"),
                          eval.expr = T)

# Input file - Rmd
inputFile <- file.path(config$params$projectdir, "Rmd", "Sample_QC.Rmd")

  message("Writing QC report for all samples in the experiment.")
  # Output file - HTML
  filename <- paste0("Study-wide_Sample_QC",
                     config$params$platform, "_",
                     config$params$project_name, "_",
                     format(Sys.time(),'%d-%m-%Y.%H.%M'),
                     ".html")
  
  outFile <- file.path(config$params$projectdir,
                       "reports",
                       filename)
  
  dir.create(file.path(config$params$projectdir,"reports"))
  
  rmarkdown::render(input = inputFile,
                    encoding = "UTF-8",
                    output_file = outFile,
                    params = config$params,
                    envir = new.env())
  