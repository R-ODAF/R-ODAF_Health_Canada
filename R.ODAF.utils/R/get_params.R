#' Retrieve and Set Up Parameters for a Given Context
#'
#' This function reads configuration settings from a YAML file and sets up the parameters for the
#' specified context, either 'analysis' or 'QC'. It prepares the parameters list by combining common
#' settings with context-specific settings and replaces any NULLs with NA. It also determines the
#' project directory and loads species-specific data.
#'
#' @param context Character string specifying the context for parameter setup; accepts 'analysis' or 'QC'.
#' @importFrom yaml read_yaml
#' @importFrom here here
#' @return A list of parameters set up for the specified context.
#' @export
get_params <- function(context = "analysis") {
  message("Reading config file...")
  config <- yaml::read_yaml(here::here("inputs", "config", "config.yaml"), eval.expr = TRUE)
  message("Config file read successfully.")
  message("Setting up parameters...")
  if (context == "analysis") {
    params <- c(config$common, config$DESeq2)
    message("Setting up species data...")
    params$species_data <- load_species(params$species,
                                        params$wikipathways_filename,
                                        params$biospyder_manifest_file)
    message("Setting up platform specific parameters...")
    params <- set_up_platform_params(params)
    check_required_params(params) # TODO - run a different set of checks for QC params?
  } else if (context == "QC") {
    params <- c(config$common, config$QC)
  } else {
    stop("Invalid context: pick one of 'analysis' or 'QC'")
  }
  message("Replacing nulls in params with NA...")
  # replace nulls in params with NA
  params <- replace_nulls_in_config(params)

  message("Setting up project directory...")
  # If projectdir is not set, figure out current project root directory
  projectdir <- params$projectdir
  if (is.na(projectdir)) {
    projectdir <- here::here()
    params$projectdir <- projectdir
  }
  return(params)
}

#' Replace NULL Values in Configuration List with NAs
#'
#' This function replaces `NULL` values in a configuration list with `NA`
#' to avoid issues with unassigned list elements during the configuration setup.
#'
#' @param params A list of configuration parameters.
#' @return The configuration list with `NULL` values replaced by `NA`.
#' @export
replace_nulls_in_config <- function(params) {
  params_new <- list()
  for (name in names(params)) {
    param <- params[[name]]
    if (is.null(param)) {
      params_new[[name]] <- NA
    } else {
      params_new[[name]] <- param
    }
  }
  return(params_new)
}
