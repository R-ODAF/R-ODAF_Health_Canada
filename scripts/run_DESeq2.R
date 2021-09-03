#!/usr/bin/R
# Custom parameters for the report
suppressMessages(library('tidyverse'))
suppressMessages(library('yaml'))
suppressMessages(library('DESeq2'))
suppressMessages(library('gtools'))

# assume this is being run from within the R project
projectdir <- here::here()
print(projectdir)

source(here::here("scripts","setup_functions.R"))

config <- yaml::read_yaml(here::here("config","config.yaml"), eval.expr = T)
params <- config$params

paths <- set_up_paths(params)
get_analysis_id <- get_analysis_id(params)

# Identify where metadata can be found
SampleKeyFile <- file.path(paths$metadata, "metadata.QC_applied.txt")
ContrastsFile <- file.path(paths$metadata, "contrasts.txt")

# Read in metadata
DESeqDesign <- read.delim(SampleKeyFile,
                          stringsAsFactors = FALSE,
                          sep = "\t",
                          header = TRUE,
                          quote = "\"",
                          row.names = 1) # Column must have unique IDs!!
DESeqDesign$original_names <- rownames(DESeqDesign)
DESeqDesignAsRead <- DESeqDesign
print(DESeqDesign)

# read in contrasts
contrasts <- read.delim(ContrastsFile, stringsAsFactors = FALSE, sep = "\t", header = FALSE,  quote = "\"")

# set interesting groups
intgroup <- params$intgroup # "Interesting groups" - experimental group/covariates
# TODO: maybe here's where the factor conversion should happen?
# TODO: what about nuisance variables?

DESeqDesign <- filter_metadata(DESeqDesign, params)
if(!is.na(params$sortcol)){
  sorted <- sort_metadata(DESeqDesign, contrasts, params)
  DESeqDesign <- sorted$design
  contrasts <- sorted$contrasts
}

print(DESeqDesign)

# Identify where count data can be found
if (params$platform == "TempO-Seq") {
  SampleDataFile <- file.path(paths$processed, "count_table.csv")
  sampledata_sep = ","
} else {
  SampleDataFile <- file.path(paths$processed, "genes.data.tsv")
  sampledata_sep = "\t"
}


if (is.na(params$group_facet)) { # all data in one facet
  message("Writing a single report for whole experiment.")
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
      message(paste0("Making multiple reports based on ", params$group_facet, "..."))

    }
}
