#!/usr/bin/R
# Custom parameters for the report
library(tidyverse)
require(yaml)

Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio/bin/pandoc")


config <- yaml::read_yaml(file.path(here::here(),
                                    "config/config.yaml"),
                          eval.expr = T)

projectdir <- config$QC$projectdir
if (is.null(projectdir)) {
  projectdir <- here::here()
  config$DESeq2$projectdir <- projectdir
}

skip_extra <- c("DMSO Pool 1", "DMSO Pool 2", "DMSO Pool 3", "DMSO Pool 4") # Remove DMSO controls as a facet

# Input file - Rmd
inputFile <- file.path(config$DESeq2$projectdir, "Rmd", "DESeq2_report.rnaseq.Rmd")

# Identify where metadata can be found
SampleKeyFile <- file.path(config$DESeq2$projectdir,
                           "data/metadata/metadata.QC_applied.txt")

# replace nulls in params with NA
params = list()
for (name in names(config$DESeq2)) {
  param <- config$DESeq2[[name]]
  if(is.null(param)){
    params[[name]] <- NA
  } else {
    params[[name]] <- param
  }
}

# Read in metadata
DESeqDesign <- read.delim(SampleKeyFile,
                          stringsAsFactors = FALSE,
                          sep = "\t",
                          header = TRUE,
                          quote = "\"",
                          row.names = 1) # Column must have unique IDs!!
DESeqDesign$original_names <- rownames(DESeqDesign)

# Run DESeq2 and make reports
if (is.na(params$group_facet)) {
  message("Writing a single report for whole experiment.")
  # Output file - HTML
  filename <- paste0(params$platform, "_",
                     params$project_name, "_",
                     format(Sys.time(),'%d-%m-%Y.%H.%M'),
                     ".html")
  outFile <- file.path(params$projectdir,
                       "reports",
                       filename)
  rmarkdown::render(input = inputFile,
                    encoding = "UTF-8",
                    output_file = outFile,
                    params = params,
                    envir = new.env())
} else if (any(!is.na(params$group_filter))) {
  message(paste0("The group(s) of interest is (are) ",
                 paste(params$group_filter, collapse = " and "),".\n",
                 "Writing a single report for that (those) groups."))
  # Output file - HTML
  filename <- paste0(params$platform, "_",
                     params$project_name, "_",
                     paste(params$group_filter, collapse = "_"), "_",
                     format(Sys.time(),'%d-%m-%Y.%H.%M'),
                     ".html")
  outFile <- file.path(params$projectdir,
                       "reports",
                       filename)
  rmarkdown::render(input = inputFile,
                    encoding = "UTF-8",
                    output_file = outFile,
                    params = params,
                    envir = new.env())
} else {
  # Remove params$exclude_groups
  facets <- DESeqDesign %>%
    filter(!(!!sym(params$group_facet)) %in%
             c(params$exclude_groups, skip_extra)) %>%
    pull(params$group_facet) %>% 
    unique()
  message(paste0("Making multiple reports based on ",
                 params$group_facet ,"..."))
  for (i in facets[1:length(facets)]) {
    message(paste0("Building report for ", i, "..."))
    params$group_filter <- i
    filename <- paste0(params$platform, "_",
                       params$project_name, "_",
                       i, "_",
                       format(Sys.time(),'%d-%m-%Y.%H.%M'),
                       ".html")
    outFile <- file.path(config$DESeq2$projectdir,
                         "reports",
                         filename)
    rmarkdown::render(input = inputFile,
                      encoding = "UTF-8",
                      output_file = outFile,
                      params = params,
                      envir = new.env())
  }
}
