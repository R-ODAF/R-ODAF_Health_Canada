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

assignInNamespace("clean_tmpfiles", clean_tmpfiles_mod, ns = "rmarkdown")

source(here::here("scripts","file_functions.R"))
source(here::here("scripts","setup_functions.R"))

config <- yaml::read_yaml(file.path(here::here(),
                                    "config/config.yaml"),
                          eval.expr = T)

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
species_data <- load_species(params$species)
params$species_data <- species_data
params <- set_up_platform_params(params)
check_required_params(params)

# Load DEGs list
dataFile <- file.path(paths$RData, paste0(params$project_title, "_DEG_data.RData"))
#load(dataFile) # metadata, contrasts, counts, resultsList
attach(dataFile)
mergedDEGsList <- mergedDEGsList
detach()

digits = 2 # For rounding numbers

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

render_report <- function(report_in, report_out, pars) {
  message("Generating report...")
  random_tmp <- file.path("/tmp",paste0("intermediates_", stringi::stri_rand_strings(1, 10)))
  rmarkdown::render(input = report_in,
                    encoding = "UTF-8",
                    output_file = report_out,
                    #output_dir = random_tmp,
                    params = pars,
                    envir = new.env(),
                    clean = TRUE,
                    intermediates_dir = random_tmp)
  system(paste0("rm -rf ", random_tmp))
}

# Determine filename prefix based on existing parameters
get_prefix <- function(params = pars, facet) {
    # Are there ways this logic might break?
  # params <- as.list(params)
  if (is.na(params$display_group_facet)) {
    message("Writing a single report for whole experiment.")
    prefix <- paste0(params$platform, "_",
                     params$project_title, "_",
                     format(Sys.time(),'%d-%m-%Y.%H.%M'))
  } else if (any(!is.na(params$display_group_filter))) {
    message(paste0("The group(s) of interest is (are) ",
                   paste(params$display_group_filter, collapse = " and "),".\n",
                   "Writing a single report for that (those) groups."))
    prefix <- paste0(params$platform, "_",
                     params$project_title, "_",
                     paste(params$display_group_filter, collapse = "_"), "_",
                     format(Sys.time(),'%d-%m-%Y.%H.%M'))
  } else {
    message(paste0("Building report for ", facet, "..."))
    prefix <- paste0(params$platform, "_",
                     params$project_title, "_",
                     facet, "_",
                     format(Sys.time(),'%d-%m-%Y.%H.%M'))
  }
  prefix <- gsub(" ", "_", prefix)
  prefix <- fs::path_sanitize(prefix)
  return(prefix)
}

# Write the reports
# These are the functions that run in parallel
# Necessary to run them per-facet and get the file prefix within each instance
make_main_reports <- function(params = pars, facet) {
  params$display_group_filter <- facet
  prefix <- get_prefix(facet = facet)
  if(params$generate_main_report){
    main_report <- file.path(projectdir, "Rmd", "DESeq2_report_new.Rmd")
    main_file <- file.path(report_dir, paste0(prefix,".html"))
    render_report(main_report, main_file, params)
  }
}

make_stats_reports <- function(params = pars, facet) {
  params$display_group_filter <- facet
  prefix <- get_prefix(facet = facet)
  if(params$generate_extra_stats_report){
    message("Generating extra stats report")
    extra_stats_report <- file.path(projectdir, "Rmd", "extra_stats_report.Rmd")
    extra_stats_file <- file.path(report_dir, paste0("extra_stats_",prefix,".html"))
    render_report(extra_stats_report, extra_stats_file, params)
  }
}

make_data_reports <- function(params = pars, facet) {
  params$display_group_filter <- facet
  prefix <- get_prefix(facet = facet)
  if(params$generate_data_explorer_report){
    message("Generating data explorer report")
    data_explorer_report <- file.path(projectdir, "Rmd", "data_explorer_report.Rmd")
    data_explorer_file <- file.path(report_dir, paste0("data_explorer_",prefix,".html"))  
    render_report(data_explorer_report, data_explorer_file, params)
  }
}

make_pathway_reports <- function(params = pars, facet)  {
  params$display_group_filter <- facet
  prefix <- get_prefix(facet = facet)
  if(params$generate_go_pathway_report){
    message("Generating GO and pathway analysis report")
    go_pathway_report <- file.path(projectdir, "Rmd", "go_pathway_report.Rmd")
    go_pathway_file <- file.path(report_dir, paste0("go_pathway_",prefix,".html"))
    render_report(go_pathway_report, go_pathway_file, params)
  }
}

# Generic function to build and filter facets
get_facets <- function(metadata = DESeqDesign,
                       exclude = params$exclude_groups,
                       display_facet = params$display_group_facet,
                       skip_extra = "DMSO") {
  facets <- metadata %>%
    filter(!(!!sym(display_facet)) %in%
             c(exclude, skip_extra)) %>%
    pull(display_facet) %>% 
    unique()
  message(paste0("Making multiple reports based on ",
                 display_facet ,"..."))
  return(facets)
}

# Determine whether facets are needed
# And if so, what they should be
  # Case 1: no facet, no display facet
if(is.na(params$group_facet) && is.na(params$display_group_facet)){
  facets <- NULL
  # Case 2: no facet, yes display facet
} else if(is.na(params$group_facet) && !is.na(params$display_group_facet)){
  facets <- get_facets()
  # Case 3: yes facet, yes display facet
} else if(!is.na(params$group_facet) && !is.na(params$display_group_facet)){
  if(params$group_facet != params$display_group_facet) {
    stop("Error: display_group_facet must match group_facet, otherwise DESeq2 results get mixed and matched.")
  }
  facets <- get_facets()
  # Which facets have DEGs?
  hasDEGs <- names(which(sapply(X = mergedDEGsList,
                                FUN = function(i) length(i)>=1),
                         arr.ind = T))
  facets <- facets[facets %in% hasDEGs]
  # Case 4: yes facet, no display facet, this one doesn't make sense
} else {
  stop("Making a single report for faceted data not supported. Did you forget to set display_group_facet?")
}

pars <- params

BPPARAM <- BiocParallel::MulticoreParam(workers = round(params$cpus*0.9))

BiocParallel::bpmapply(FUN = make_main_reports, facet = facets, BPPARAM = BPPARAM)
BiocParallel::bpmapply(FUN = make_data_reports, facet = facets, BPPARAM = BPPARAM)
BiocParallel::bpmapply(FUN = make_pathway_reports, facet = facets, BPPARAM = BPPARAM)

# Why does this use so much memory!?
# It is knitr::kable that is the problem.
# See the extra stats Rmd for details, the chunk named metadata-report
reduced_cpus <- round(params$cpus*0.2)
if (reduced_cpus < 1) { reduced_cpus = 1 }
BPPARAM <- BiocParallel::MulticoreParam(workers = reduced_cpus )
BiocParallel::bpmapply(FUN = make_stats_reports, facet = facets, BPPARAM = BPPARAM)

# Add back after troubleshooting above code...
source(here::here(file.path("scripts","summarize_across_facets.R")))

# NOTE Manually clean up temporary files
# This is required because of the clean_tmpfiles_mod() workaround!
system("rm -rf /tmp/intermediates_*")
