#!/usr/bin/R
# Custom parameters for the report
library(R.ODAF.utils)

# Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio/bin/pandoc")

params <- R.ODAF.utils::get_params(context = "QC")

# Input file - Rmd
input_file <- file.path(params$projectdir, "Rmd", "Sample_QC.Rmd")

message("Writing QC report for all samples in the experiment.")
# Output file - HTML
filename <- paste0("Study-wide_Sample_QC_",
                   params$platform, "_",
                   params$project_title, "_",
                   format(Sys.time(),'%d-%m-%Y.%H.%M'),
                   ".html")
out_dir <- file.path(params$projectdir, "output", "QC")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

out_file <- file.path(out_dir, filename)

rmarkdown::render(input = input_file,
                  encoding = "UTF-8",
                  output_file = out_file,
                  params = params,
                  envir = new.env())
