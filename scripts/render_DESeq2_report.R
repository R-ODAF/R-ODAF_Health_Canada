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

#Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio/bin/pandoc")


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
                          row.names = NULL) # Column must have unique IDs!!
# 
# DESeqDesignPreQC <- read.delim(file.path(config$DESeq2$projectdir,
#                                     "data/metadata/metadata.txt"),
#                           stringsAsFactors = FALSE,
#                           sep = "\t",
#                           header = TRUE,
#                           quote = "\"",
#                           row.names = NULL) # Column must have unique IDs!!
# 
# # How many groups in total are being input?
# DESeqDesign %>% group_by(group) %>% count()
# DESeqDesignPreQC  %>% group_by(group) %>% count()
# # 386
# # 420
# 
# # Manual removal of some groups (e.g., cytotoxicity?)
# remove <- read.table(file.path(config$DESeq2$projectdir,"./data/metadata/remove.txt"),
#                      sep="\t", header = F) %>% pull()
# length(remove)
# # 48
# # Expect 420 - 48 = 372 to remain
# 
# test <- DESeqDesignPreQC %>% dplyr::filter(!group %in% remove)
# cytotoxic <- DESeqDesignPreQC %>% dplyr::filter(group %in% remove)
# test %>% group_by(group) %>% count()
# # 373??? Mix6 33 10d - not in original metadata...?
# 
# # Once above is correct, go ahead and make final filtered metadata
# DESeqDesign <- DESeqDesign %>% dplyr::filter(!group %in% remove)
# # If you've already done this, don't overwrite the backup!
# if (!file.exists(paste0(SampleKeyFile,".bak"))) {
#   system(paste0("cp ",SampleKeyFile," ",paste0(SampleKeyFile,".bak")))
# }
# write.table(DESeqDesign,
#             SampleKeyFile,
#             sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

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
  facets <- facets[grep(pattern = "DMSO", x = facets, invert = T)]
  message(paste0("Making multiple reports based on ",
                 params$group_facet ,"..."))
  print(facets)
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
  deg_files <- fs::dir_ls(file.path(config$DESeq2$projectdir, "DEG_output"),
                          regexp = "\\-DEG_summary.txt$", recurse = T)
  # This depends on 'cat' being available on the command line (i.e., linux-specific)
  # Also some insane quoting going on here, but I don't see an easier way
  system(paste0('cat "', paste(deg_files, collapse='"  "'),
                '"  >  ', file.path(config$DESeq2$projectdir,
                                 "DEG_output/DEG_summary.txt")))
  
  # This would probably fail in cases where different numbers of contrasts exists across facets.
  # But could otherwise be useful?
  # results <- deg_files %>% map_dfr(read_tsv, col_names=T, .id="source") 
}
