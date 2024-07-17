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
exp_metadata <- R.ODAF.utils::get_metadata(file.path(paths$qc, "metadata.QC_applied.txt"), paths)

facets <- R.ODAF.utils::get_facets(exp_metadata, params)

# Only make reports for facets with DEGS
hasDEGs <- names(which(sapply(X = mergedDEGsList,
                                FUN = function(i) length(i)>=1),
                         arr.ind = T))
if(!is.na(params$deseq_facet)) {
  display_facets <- facets[facets %in% hasDEGs]
} else {
  display_facets <- facets
}

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
  if (params$species == "human") {
    base::mapply(FUN = make_tgxddi_reports, facet = facets, MoreArgs = list(pars = params, paths = paths))
    # Concatenate all the TGxDDI output csv files into one
    tgxddi_files <- list.files(paths$reports_dir, pattern = "_tgx-ddi_results.csv", full.names = TRUE)
    tgxddi_df <- readr::read_csv(tgxddi_files[1], show_col_types = FALSE)
    for (i in 2:length(tgxddi_files)) {
      tgxddi_df <- dplyr::bind_rows(tgxddi_df, readr::read_csv(tgxddi_files[i], show_col_types = FALSE))
    }
    # Write out the concatenated file
    readr::write_csv(tgxddi_df, file.path(paths$reports_dir, paste0("tgx-ddi_results_summary.csv")))
    # Delete the individual files
    file.remove(tgxddi_files)
  }
  else {
    message(paste("TGx-DDI analysis is only available for human datasets. Your parameters indicate that the data is from", params$species, ". Skipping TGx-DDI analysis."))
  }
}

if (params$generate_tgxhdaci_report) {
  if (params$platform == "TempO-Seq") {
    base::mapply(FUN = make_hdaci_reports, facet = facets, MoreArgs = list(pars = params, paths = paths))
    # Concatenate all the TGxHDACi output csv files into one
    tgxhdaci_files <- list.files(paths$reports_dir, pattern = "_tgx-HDACi_results.csv", full.names = TRUE)
    tgxhdaci_df <- readr::read_csv(tgxhdaci_files[1], show_col_types = FALSE)
    for (i in 2:length(tgxhdaci_files)) {
      tgxhdaci_df <- dplyr::bind_rows(tgxhdaci_df, readr::read_csv(tgxhdaci_files[i], show_col_types = FALSE))
    }
    # Write out the concatenated file
    readr::write_csv(tgxhdaci_df, file.path(paths$reports_dir, paste0("tgx-hdaci_results_summary.csv")))
    # Delete the individual files
    file.remove(tgxhdaci_files)
  }
  else {
    message("TGxhdaci report generation is currently only supported for TempO-Seq data.")
  }
}

# Add back after troubleshooting above code...
if (!is.na(params$reports_facet)) {
  summarize_across_facets(overallResListAll, overallResListDEGs, filtered_table, facets, params)
}

# NOTE Manually clean up temporary files
# This is required because of the clean_tmpfiles_mod() workaround!
system("rm -rf /tmp/intermediates_*")

# save config and contrasts file too
file.copy(file.path(paths$inputs, "config", "config.yaml"), paths$record)
file.copy(file.path(paths$contrasts, params$contrasts_file), paths$record)
