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

skip_extra <- c("DMSO") # Remove DMSO controls as a facet
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
  rmarkdown::render(input = report_in,
                    encoding = "UTF-8",
                    output_file = report_out,
                    params = pars,
                    envir = new.env())
}

# convenience function
make_reports <- function(params, facet) {
  # Determine filename prefix based on existing parameters
  # Are there ways this logic might break?
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
    params$display_group_filter <- facet
    prefix <- paste0(params$platform, "_",
                     params$project_title, "_",
                     facet, "_",
                     format(Sys.time(),'%d-%m-%Y.%H.%M'))
  }
  prefix <- gsub(" ", "_", prefix)
  prefix <- fs::path_sanitize(prefix)

  # Write the reports
  if(params$generate_main_report){
    main_report <- file.path(projectdir, "Rmd", "DESeq2_report_new.Rmd")
    main_file <- file.path(report_dir, paste0(file_prefix,".html"))
    render_report(main_report, main_file, params)
  }
  if(params$generate_extra_stats_report){
    message("Generating extra stats report")
    extra_stats_report <- file.path(projectdir, "Rmd", "extra_stats_report.Rmd")
    extra_stats_file <- file.path(report_dir, paste0("extra_stats_",file_prefix,".html"))
    render_report(extra_stats_report, extra_stats_file, params)
  }
  if(params$generate_data_explorer_report){
    message("Generating data explorer report")
    data_explorer_report <- file.path(projectdir, "Rmd", "data_explorer_report.Rmd")
    data_explorer_file <- file.path(report_dir, paste0("data_explorer_",file_prefix,".html"))  
    render_report(data_explorer_report, data_explorer_file, params)
  }
  if(params$generate_go_pathway_report){
    message("Generating GO and pathway analysis report")
    go_pathway_report <- file.path(projectdir, "Rmd", "go_pathway_report.Rmd")
    go_pathway_file <- file.path(report_dir, paste0("go_pathway_",file_prefix,".html"))
    render_report(go_pathway_report, go_pathway_file, params)
  }
}

# Determine whether facets are needed
if (!is.na(params$display_group_facet) && is.na(params$display_group_filter)) {
  facets <- DESeqDesign %>%
    filter(!(!!sym(params$display_group_facet)) %in%
             c(params$exclude_groups, skip_extra)) %>%
    pull(params$display_group_facet) %>% 
    unique()
  facets <- facets[grep(pattern = "DMSO", x = facets, invert = T)]
  message(paste0("Making multiple reports based on ",
                 params$display_group_facet ,"..."))
} else {
  facets <- NULL
}

#### make_reports(params, facets)
  
  # parallel::mcmapply(FUN = render_reports_parallel, facets, mc.cores = params$cpus/2)
  
  #library(doParallel)
  n_cores <- parallel::detectCores()
  cluster <- parallel::makeCluster(n_cores-1)
  doParallel::registerDoParallel(cluster)
  # 
  parallel::clusterMap(cl = cluster, make_reports, params = params, facet = facets)
  #   foreach(i=seq_along(facets), .combine='c', .export = ls(globalenv())) %dopar% { # Changing to %dopar% fails.
  #     print(facets[i])
  #     render_reports_parallel(facets[i])
  #   }
  source(here::here(file.path("scripts","summarize_across_facets.R")))
  # TODO: reproduce these files
  # deg_files <- fs::dir_ls(deglist_dir, regexp = "\\-DEG_summary.txt$", recurse = T)
  # # This depends on 'cat' being available on the command line (i.e., linux-specific)
  # # Also some insane quoting going on here, but I don't see an easier way
  # system(paste0('cat "', paste(deg_files, collapse='"  "'),
  #               '"  >  ', file.path(deglist_dir, "DEG_summary.txt")))
  # 
  # This would probably fail in cases where different numbers of contrasts exists across facets.
  # But could otherwise be useful?
  # results <- deg_files %>% map_dfr(read_tsv, col_names=T, .id="source") 
}
