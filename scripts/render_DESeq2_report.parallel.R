#!/usr/bin/R
# Custom parameters for the report

# If you want to run this outside RStudio, you need to tell R where Pandoc is:
# Find it in RStudio:
# Sys.getenv("RSTUDIO_PANDOC")
# Set it in the terminal:
# Sys.setenv(RSTUDIO_PANDOC="OUTPUT FROM ABOVE COMMAND")
# Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio-server/bin/pandoc")

library('tidyverse')
library('devtools')
library('DESeq2')
library('data.table')
library('yaml')
library('BiocParallel')
library('rmarkdown')
clean_tmpfiles_mod <- function() {
  message("Calling clean_tmpfiles_mod()")
}

single_facet_constant = "single_facet_constant_12345" # Used internally; don't have a facet named this.

assignInNamespace("clean_tmpfiles", clean_tmpfiles_mod, ns = "rmarkdown")

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if at least one argument is provided
if (length(args) > 0) {
  # Assume the first argument is the new location
  results_location_arg <- args[1]

  # Source functions and pass along analysis directory argument
  source(here::here("scripts","file_functions.R"))
  
} else {
  message("Error: Missing argument. Provide the analysis directory name as an argument.\n")
}

source(here::here("scripts","setup_functions.R"))

config <- yaml::read_yaml(here::here("inputs","config","config.yaml"), eval.expr = T)

# Combine required params from config
params <- c(config$common, config$DESeq2)
# replace nulls in params with NA
params <- replace_nulls_in_config(params)
# If projectdir is not set, figure out current project root directory
projectdir <- params$projectdir
if (is.na(projectdir)) {
  projectdir <- here::here()
  params$projectdir <- projectdir
}

paths <- set_up_paths(params)
species_data <- load_species(params$species, params$wikipathways_filename, params$biospyder_manifest_file)
params$species_data <- species_data
params <- set_up_platform_params(params)
check_required_params(params)

# Load DEGs list
dataFile <- file.path(paths$RData, paste0(params$project_title, "_DEG_data.RData"))
#load(dataFile) # metadata, contrasts, counts, resultsList
attach(dataFile)
mergedDEGsList <- mergedDEGsList
detach()

# Identify where metadata can be found
exp_metadata_file <- file.path(paths$metadata, "metadata.QC_applied.txt")

# Read in metadata
exp_metadata <- read.delim(exp_metadata_file,
                          stringsAsFactors = FALSE,
                          sep = "\t",
                          header = TRUE,
                          quote = "\"",
                          row.names = 1) # Column must have unique IDs!!
exp_metadata$original_names <- rownames(exp_metadata)

# Determine whether facets are needed
# And if so, what they should be
# Case 1: no facet, no display facet
if(is.na(params$group_facet) && is.na(params$reports_facet)){
  display_facets <- single_facet_constant
  # Case 2: no facet, yes display facet
} else if(is.na(params$group_facet) && !is.na(params$reports_facet)){
  display_facets <- get_facets()
  # Case 3: yes facet, yes display facet
} else if(!is.na(params$group_facet) && !is.na(params$reports_facet)){
  if(params$group_facet != params$reports_facet) {
    stop("Error: reports_facet must match group_facet, otherwise DESeq2 results get mixed and matched.")
  }
  display_facets <- get_facets()
  # Which facets have DEGs?
  hasDEGs <- names(which(sapply(X = mergedDEGsList,
                                FUN = function(i) length(i)>=1),
                         arr.ind = T))
  display_facets <- display_facets[display_facets %in% hasDEGs]
  # Case 4: yes facet, no display facet, this one doesn't make sense
} else {
  stop("Making a single report for faceted data not supported. Did you forget to set reports_facet?")
}

paths <- set_up_paths_3(paths,params,display_facets)

# Make directory for DESeq2 Reports
report_dir <- file.path(paths$results, "DEG_reports")
deglist_dir <- file.path(paths$results, "DEG_lists")
if (!dir.exists(report_dir)) {dir.create(report_dir, recursive = TRUE)}
if (!dir.exists(deglist_dir)) {dir.create(deglist_dir, recursive = TRUE)}

if (params$parallel){
  biocluster <- BiocParallel::MulticoreParam(workers = round(params$cpus*0.9))
  BiocParallel::bpmapply(FUN = make_main_reports, facet = display_facets, MoreArgs = list(pars = params), BPPARAM = biocluster)
  if (params$generate_main_report == T) {Sys.sleep(60)}
  BiocParallel::bpmapply(FUN = make_data_reports, facet = display_facets, MoreArgs = list(pars = params), BPPARAM = biocluster)
  if (params$generate_data_explorer_report  == T) {Sys.sleep(60)}
  BiocParallel::bpmapply(FUN = make_pathway_reports, facet = display_facets, MoreArgs = list(pars = params), BPPARAM = biocluster)
  if (params$generate_go_pathway_report  == T) {Sys.sleep(60)}
  # Why does this use so much memory!?
  # It is knitr::kable that is the problem.
  # See the extra stats Rmd for details, the chunk name metadata-report
  #reduced_cpus <- round(params$cpus*0.5)
  #if (reduced_cpus < 1) { reduced_cpus = 1 }
  #biocluster <- BiocParallel::MulticoreParam(workers = reduced_cpus )
  BiocParallel::bpmapply(FUN = make_stats_reports, facet = display_facets, MoreArgs = list(pars = params), BPPARAM = biocluster)
} else {
  base::mapply(FUN = make_main_reports, facet = display_facets, MoreArgs = list(pars = params))
  base::mapply(FUN = make_data_reports, facet = display_facets, MoreArgs = list(pars = params))
  base::mapply(FUN = make_pathway_reports, facet = display_facets, MoreArgs = list(pars = params))
  base::mapply(FUN = make_stats_reports, facet = display_facets, MoreArgs = list(pars = params))
}

# Add back after troubleshooting above code...
if (!is.na(params$group_facet)) {
  source(here::here(file.path("scripts","summarize_across_facets.R")))
}

# NOTE Manually clean up temporary files
# This is required because of the clean_tmpfiles_mod() workaround!
system("rm -rf /tmp/intermediates_*")

# save config and contrasts file too
file.copy(file.path(paths$inputs, "config", "config.yaml"), paths$record)
file.copy(file.path(paths$contrasts,"contrasts.txt"), paths$record)
