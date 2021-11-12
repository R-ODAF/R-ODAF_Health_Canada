#!/usr/bin/R
# Custom parameters for the report

# If you want to run this outside RStudio, you need to tell R where Pandoc is:
# Find it in RStudio:
# Sys.getenv("RSTUDIO_PANDOC")
# Set it in the terminal:
# Sys.setenv(RSTUDIO_PANDOC="OUTPUT FROM ABOVE COMMAND")
# Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio-server/bin/pandoc")

library(tidyverse)
require(yaml)

source(here::here("scripts","file_functions.R"))

config <- yaml::read_yaml(file.path(here::here(),
                                    "config/config.yaml"),
                          eval.expr = T)

# Combine required params from config
params <- c(config$common, config$DESeq2)
# If projectdir is not set, figure out current project root directory
projectdir <- params$projectdir
if (is.null(projectdir)) {
  projectdir <- here::here()
  params$projectdir <- projectdir
}

paths <- set_up_paths(config$DESeq2)

# replace nulls in params with NA
replace_nulls_in_config(config$DESeq2)

skip_extra <- c("DMSO") # Remove DMSO controls as a facet

# Input file - Rmd
inputFile <- file.path(projectdir, "Rmd", "DESeq2_report_new.Rmd")

# Identify where metadata can be found
SampleKeyFile <- file.path(projectdir, "data/metadata/metadata.QC_applied.txt")

# Read in metadata
DESeqDesign <- read.delim(SampleKeyFile,
                          stringsAsFactors = FALSE,
                          sep = "\t",
                          header = TRUE,
                          quote = "\"",
                          row.names = 1) # Column must have unique IDs!!
DESeqDesign$original_names <- rownames(DESeqDesign)

# Make directory for DESeq2 Reports
report_dir <- file.path(projectdir, "analysis", "DEG_reports")
deglist_dir <- file.path(projectdir, "analysis", "DEG_lists")
if (!dir.exists(report_dir)) {dir.create(report_dir, recursive = TRUE)}
if (!dir.exists(deglist_dir)) {dir.create(deglist_dir, recursive = TRUE)}

# Run DESeq2 and make reports
if (is.na(params$display_group_facet)) {
  message("Writing a single report for whole experiment.")
  # Output file - HTML
  filename <- paste0(params$platform, "_",
                     params$project_name, "_",
                     format(Sys.time(),'%d-%m-%Y.%H.%M'),
                     ".html")
  outFile <- file.path(report_dir, filename)
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
  outFile <- file.path(report_dir, filename)
  rmarkdown::render(input = inputFile,
                    encoding = "UTF-8",
                    output_file = outFile,
                    params = params,
                    envir = new.env())
} else {
  # Remove params$exclude_groups
  facets <- DESeqDesign %>%
    filter(!(!!sym(params$display_group_facet)) %in%
             c(params$exclude_groups, skip_extra)) %>%
    pull(params$display_group_facet) %>% 
    unique()
  facets <- facets[grep(pattern = "DMSO", x = facets, invert = T)]
  message(paste0("Making multiple reports based on ",
                 params$display_group_facet ,"..."))
  print(facets)
  for (i in facets[1:length(facets)]) {
    message(paste0("Building report for ", i, "..."))
    params$group_filter <- i
    filename <- paste0(params$platform, "_",
                       params$project_name, "_",
                       i, "_",
                       format(Sys.time(),'%d-%m-%Y.%H.%M'),
                       ".html")
    outFile <- file.path(report_dir, filename)
    rmarkdown::render(input = inputFile,
                      encoding = "UTF-8",
                      output_file = outFile,
                      params = params,
                      envir = new.env())
  }
  deg_files <- fs::dir_ls(deglist_dir, regexp = "\\-DEG_summary.txt$", recurse = T)
  # This depends on 'cat' being available on the command line (i.e., linux-specific)
  # Also some insane quoting going on here, but I don't see an easier way
  system(paste0('cat "', paste(deg_files, collapse='"  "'),
                '"  >  ', file.path(deglist_dir, "DEG_summary.txt")))
  
  # This would probably fail in cases where different numbers of contrasts exists across facets.
  # But could otherwise be useful?
  # results <- deg_files %>% map_dfr(read_tsv, col_names=T, .id="source") 
}
