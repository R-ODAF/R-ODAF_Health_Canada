#!/usr/bin/R
# Custom parameters for the report
library(tidyverse)
require(yaml)

source(here::here("scripts","setup_functions.R"))

#Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio/bin/pandoc")

config <- yaml::read_yaml(here::here("inputs","config","config.yaml"), eval.expr = T)

# Combine required params from config
params <- c(config$common, config$QC)
# replace nulls in params with NA
params <- replace_nulls_in_config(params)
# If projectdir is not set, figure out current project root directory
projectdir <- params$projectdir
if (is.na(projectdir)) {
  projectdir <- here::here()
  params$projectdir <- projectdir
}

# Input file - Rmd
input_file <- file.path(projectdir, "Rmd", "Sample_QC.Rmd")

message("Writing QC report for all samples in the experiment.")
# Output file - HTML
filename <- paste0("Study-wide_Sample_QC_",
                   params$platform, "_",
                   params$project_title, "_",
                   format(Sys.time(),'%d-%m-%Y.%H.%M'),
                   ".html")
out_dir <- file.path(projectdir, "output", "QC")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

out_file <- file.path(out_dir, filename)

rmarkdown::render(input = input_file,
                  encoding = "UTF-8",
                  output_file = out_file,
                  params = params,
                  envir = new.env())
