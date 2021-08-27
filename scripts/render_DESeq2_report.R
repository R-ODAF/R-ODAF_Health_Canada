#!/usr/bin/R
# Custom parameters for the report
library(tidyverse)
require(yaml)

config <- yaml::read_yaml(file.path(here::here(),
                                    "Rmd/config.yml"),
                          eval.expr = T)

skip_extra <- c("DMSO", "DMSO 24 h", "DMSO 96 h") # Remove DMSO controls as a facet

# Input file - Rmd
inputFile <- file.path(config$params$projectdir, "Rmd", "DESeq2_report.rnaseq.Rmd")

# Identify where metadata can be found
SampleKeyFile <- file.path(config$params$projectdir,
                           "metadata/metadata.QC_applied.txt")

# Read in metadata
DESeqDesign <- read.delim(SampleKeyFile,
                          stringsAsFactors = FALSE,
                          sep = "\t",
                          header = TRUE,
                          quote = "\"",
                          row.names = 1) # Column must have unique IDs!!
DESeqDesign$original_names <- rownames(DESeqDesign)

# Run DESeq2 and make reports
if (is.na(config$params$group_facet)) {
  message("Writing a single report for whole experiment.")
  # Output file - HTML
  filename <- paste0(config$params$platform, "_",
                     config$params$project_name, "_",
                     format(Sys.time(),'%d-%m-%Y.%H.%M'),
                     ".html")
  outFile <- file.path(config$params$projectdir,
                       "reports",
                       filename)
  rmarkdown::render(input = inputFile,
                    encoding = "UTF-8",
                    output_file = outFile,
                    params = config$params,
                    envir = new.env())
} else if (any(!is.na(config$params$group_filter))) {
  message(paste0("The group(s) of interest is (are) ",
                 paste(config$params$group_filter, collapse = " and "),".\n",
                 "Writing a single report for that (those) groups."))
  # Output file - HTML
  filename <- paste0(config$params$platform, "_",
                     config$params$project_name, "_",
                     paste(config$params$group_filter, collapse = "_"), "_",
                     format(Sys.time(),'%d-%m-%Y.%H.%M'),
                     ".html")
  outFile <- file.path(config$params$projectdir,
                       "reports",
                       filename)
  rmarkdown::render(input = inputFile,
                    encoding = "UTF-8",
                    output_file = outFile,
                    params = config$params,
                    envir = new.env())
} else {
  # Remove config$params$exclude_groups
  facets <- DESeqDesign %>%
    filter(!(!!sym(config$params$group_facet)) %in%
             c(config$params$exclude_groups, skip_extra)) %>%
    pull(config$params$group_facet) %>% 
    unique()
  message(paste0("Making multiple reports based on ",
                 config$params$group_facet ,"..."))
  for (i in facets) {
    message(paste0("Building report for ", i, "..."))
    config$params$group_filter <- i
    filename <- paste0(config$params$platform, "_",
                       config$params$project_name, "_",
                       i, "_",
                       format(Sys.time(),'%d-%m-%Y.%H.%M'),
                       ".html")
    outFile <- file.path(config$params$projectdir,
                         "reports",
                         filename)
    rmarkdown::render(input = inputFile,
                      encoding = "UTF-8",
                      output_file = outFile,
                      params = config$params,
                      envir = new.env())
  }
}
