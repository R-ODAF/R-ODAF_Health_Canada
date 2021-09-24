#!/usr/bin/R

suppressMessages(library('tidyverse'))
suppressMessages(library('yaml'))

projectdir <- here::here()

config <- yaml::read_yaml(here::here("config","config.new.yaml"), eval.expr = T)
params <- config$params

RMarkdown_file <- ""

# load data
# should load ddsList, designList, overallResList, rldList, DESeqDesign, facets, contrasts, intgroup, paths
# TODO: take params out of save on other script
load(file.path(paths$DEG_output, paste0(params$project_name, "_DEG_data.RData")))

# figure out if we're doing one or multiple reports
if(is.na(params$display_group_facet)){
    # one report for all data

    # list of DESeq results (indexed by comparison)
    resList <- overallResList[['all']]
} else {
    for (current_filter in facets) {
        # list of DESeq results (indexed by comparison)
        resList <- overallResList[[current_filter]]

        output_report_file = "asdf"

        params$current_filter = current_filter

        rmarkdown::render(input = RMarkdown_file,
                    encoding = "UTF-8",
                    output_file = output_report_file,
                    params = params,
                    envir = new.env())
    }
}