#!/usr/bin/R
# Custom parameters for the report

# If you want to run this outside RStudio, you need to tell R where Pandoc is:
# Find it in RStudio:
# Sys.getenv("RSTUDIO_PANDOC")
# Set it in the terminal:
# Sys.setenv(RSTUDIO_PANDOC="OUTPUT FROM ABOVE COMMAND")
# Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio-server/bin/pandoc")
# Sys.setenv(RSTUDIO_PANDOC = "/usr/lib/rstudio-server/bin/quarto/bin/tools/pandoc")

library("BiocParallel")
library("rmarkdown")
devtools::load_all("./R.ODAF.utils")

clean_tmpfiles_mod <- function() {
  message("Calling clean_tmpfiles_mod()")
}

assignInNamespace("clean_tmpfiles", clean_tmpfiles_mod, ns = "rmarkdown")

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if at least one argument is provided
if (length(args) > 0) {
  # Assume the first argument is the new location
  results_location_arg <- args[1]
} else {
  message("Error: Missing argument. Provide the analysis directory name as an argument.\n")
}
# Load project parameters
params <- R.ODAF.utils::get_params(context = "analysis")

# Set up initial project paths
paths <- R.ODAF.utils::set_up_filepaths(params, results_location_arg)

# Load DEGs list
data_file <- file.path(paths$RData, paste0(params$project_title, "_DEG_data.RData"))
data_env <- new.env()
#load(dataFile) # metadata, contrasts, counts, resultsList
load(data_file, envir = data_env)
mergedDEGsList <- data_env$mergedDEGsList
overallResListAll <- data_env$overallResListAll
overallResListDEGs <- data_env$overallResListDEGs
exp_contrasts <- data_env$exp_contrasts
filtered_table <-  data_env$filtered_table
allBiomarkers <- data_env$allBiomarkers
rm(data_env)
gc()

# Read in metadata
exp_metadata <- R.ODAF.utils::get_metadata(file.path(paths$metadata, "metadata.QC_applied.txt"), paths)

facets <- R.ODAF.utils::get_facets(exp_metadata, params)

# Only make reports for facets with DEGS
hasDEGs <- names(which(sapply(X = mergedDEGsList,
                                FUN = function(i) length(i)>=1),
                         arr.ind = T))
display_facets <- facets[facets %in% hasDEGs]

# set up the rest of the output paths (requires facets)
paths <- R.ODAF.utils::set_up_filepaths(params,
                                        results_location_arg,
                                        exp_metadata,
                                        make_report_dirs = TRUE)

if (params$parallel) {
  biocluster <- BiocParallel::MulticoreParam(workers = round(params$cpus * 0.9))
  BiocParallel::bpmapply(FUN = make_main_reports, facet = display_facets, MoreArgs = list(pars = params, paths = paths), BPPARAM = biocluster)
  if (params$generate_main_report == TRUE) {
    Sys.sleep(60)
  }
  BiocParallel::bpmapply(FUN = make_data_reports, facet = display_facets, MoreArgs = list(pars = params, paths = paths), BPPARAM = biocluster)
  if (params$generate_data_explorer_report  == TRUE) {
    Sys.sleep(60)
  }
  BiocParallel::bpmapply(FUN = make_pathway_reports, facet = display_facets, MoreArgs = list(pars = params, paths = paths), BPPARAM = biocluster)
  if (params$generate_go_pathway_report  == TRUE) {
    Sys.sleep(60)
  }
  # Why does this use so much memory!?
  # It is knitr::kable that is the problem.
  # See the extra stats Rmd for details, the chunk name metadata-report
  #reduced_cpus <- round(params$cpus*0.5)
  #if (reduced_cpus < 1) { reduced_cpus = 1 }
  #biocluster <- BiocParallel::MulticoreParam(workers = reduced_cpus )
  BiocParallel::bpmapply(FUN = make_stats_reports, facet = display_facets, MoreArgs = list(pars = params, paths = paths), BPPARAM = biocluster)
} else {
  base::mapply(FUN = make_main_reports, facet = display_facets, MoreArgs = list(pars = params, paths = paths))
  base::mapply(FUN = make_data_reports, facet = display_facets, MoreArgs = list(pars = params, paths = paths))
  base::mapply(FUN = make_pathway_reports, facet = display_facets, MoreArgs = list(pars = params, paths = paths))
  base::mapply(FUN = make_stats_reports, facet = display_facets, MoreArgs = list(pars = params, paths = paths))
}

if (params$generate_tgxddi_report) {
  base::mapply(FUN = make_tgxddi_reports, facet = display_facets, MoreArgs = list(pars = params, paths = paths))
}


# Add back after troubleshooting above code...
# if (!is.na(params$reports_facet)) {
#   summarize_across_facets(overallResListAll, overallResListDEGs, filtered_table, facets, params)
# }

# NOTE Manually clean up temporary files
# This is required because of the clean_tmpfiles_mod() workaround!
system("rm -rf /tmp/intermediates_*")

# save config and contrasts file too
file.copy(file.path(paths$inputs, "config", "config.yaml"), paths$record)
file.copy(file.path(paths$contrasts, "contrasts.txt"), paths$record)
