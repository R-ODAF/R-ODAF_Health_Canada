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




if (is.na(config$params$group_facet)) { # all data in one facet
  message(paste0("Making multiple reports based on ",
                 config$params$group_facet ,"..."))

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
