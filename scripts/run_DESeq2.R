#!/usr/bin/R
# Custom parameters for the report
suppressMessages(library('tidyverse'))
suppressMessages(library('yaml'))
suppressMessages(library('DESeq2'))

source("setup_functions.R")

projectdir <- here::here()
print(projectdir)
config <- yaml::read_yaml(file.path(projectdir, "Rmd/config.yml"), eval.expr = T)

params <- config$params

paths <- set_up_paths(params)
get_analysis_id <- get_analysis_id(params)

# Identify where metadata can be found
SampleKeyFile <- file.path(params$projectdir, "metadata/metadata.QC_applied.txt")

# Read in metadata
DESeqDesign <- read.delim(SampleKeyFile,
                          stringsAsFactors = FALSE,
                          sep = "\t",
                          header = TRUE,
                          quote = "\"",
                          row.names = 1) # Column must have unique IDs!!
DESeqDesign$original_names <- rownames(DESeqDesign)


# Identify where count data can be found
if (Platform == "TempO-Seq") {
  SampleDataFile <- file.path(paths$processed, "count_table.csv")
  sampledata_sep = ","
} else {
  SampleDataFile <- file.path(paths$processed, "genes.data.tsv")
  sampledata_sep = "\t"
}



if (is.na(params$group_facet)) { # all data in one facet
  message(paste0("Making multiple reports based on ",
                 params$group_facet ,"..."))
    # load data
    # run DEseq2

} else {
    facets <- DESeqDesign %>%
        filter(!(params$group_facet) %in% c(params$exclude_groups, skip_extra)) %>%
        filter(!solvent_control) %>%
        pull(params$group_facet) %>% 
        unique()

    if (any(!is.na(params$group_filter))) { # filter facets
        message(paste0("The group(s) of interest is (are) ",
                 paste(params$group_filter, collapse = " and "),".\n",
                 "Writing a single report for that (those) groups."))

    } else { # do all facets separately
        message("Writing a single report for whole experiment.")

    }
}
