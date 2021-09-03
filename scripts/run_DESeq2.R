#!/usr/bin/R
# Custom parameters for the report
library(tidyverse)
library(yaml)


projectdir <- here::here()
print(projectdir)
config <- yaml::read_yaml(file.path(projectdir, "Rmd/config.yml"), eval.expr = T)

# Identify where metadata can be found
SampleKeyFile <- file.path(config$params$projectdir, "metadata/metadata.QC_applied.txt")

# Read in metadata
DESeqDesign <- read.delim(SampleKeyFile,
                          stringsAsFactors = FALSE,
                          sep = "\t",
                          header = TRUE,
                          quote = "\"",
                          row.names = 1) # Column must have unique IDs!!
DESeqDesign$original_names <- rownames(DESeqDesign)

paths <- set_up_paths(params)

# Identify where count data can be found
if (Platform == "TempO-Seq") {
  SampleDataFile <- file.path(paths$processed, "count_table.csv")
  sampledata_sep = ","
} else {
  SampleDataFile <- file.path(paths$processed, "genes.data.tsv")
  sampledata_sep = "\t"
}


if (is.na(config$params$group_facet)) { # all data in one facet
  message(paste0("Making multiple reports based on ",
                 config$params$group_facet ,"..."))
    # load data
    # run DEseq2

} else {
    facets <- DESeqDesign %>%
        filter(!(config$params$group_facet) %in% c(config$params$exclude_groups, skip_extra)) %>%
        filter(!solvent_control) %>%
        pull(config$params$group_facet) %>% 
        unique()

    if (any(!is.na(config$params$group_filter))) { # filter facets
        message(paste0("The group(s) of interest is (are) ",
                 paste(config$params$group_filter, collapse = " and "),".\n",
                 "Writing a single report for that (those) groups."))

    } else { # do all facets separately
        message("Writing a single report for whole experiment.")

    }
}


set_up_paths <- function(params) {
    paths <- list()
    # Other important system paths to specify in config
    paths$wikipathways <- params$wikipathways_directory
    # For project structure
    # Should probably update this to use the file.path() function.
    paths$root <- params$projectdir
    paths$data <- file.path(paths$root, "data")
      paths$raw <- file.path(paths$data, "raw")
      paths$processed <- file.path(paths$data, "processed")
      paths$metadata <- file.path(paths$data, "metadata")
    paths$reports <- file.path(paths$root, "reports")
    paths$results <- file.path(paths$root, "results")
    if (is.na(params$group_facet)) {
      paths$DEG_output <- file.path(paths$results, "DEG_output")
    } else {
      paths$DEG_output <- file.path(paths$results, "DEG_output", paste0("group_", paste(params$group_filter, collapse = "_")))
    }
    paths$pathway_analysis <- file.path(paths$DEG_output, "/pathway_analysis")
    paths$RData <- file.path(paths$DEG_output, "/RData")
    paths$BMD_output <- file.path(paths$results, "/DEG_output/BMD_and_biomarker_files")
}