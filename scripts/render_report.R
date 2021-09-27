#!/usr/bin/R

suppressMessages(library('tidyverse'))
suppressMessages(library('yaml'))

projectdir <- here::here()

config <- yaml::read_yaml(here::here("config","config.new.yaml"), eval.expr = T)
params <- config$params


RMarkdown_file <- here::here("Rmd","DESeq2_report_new.Rmd")

source(here::here("scripts","file_functions.R"))
source(here::here("scripts","data_functions.R"))

paths <- set_up_paths(params)


# load data
# should load ddsList, designList, overallResList, rldList, DESeqDesign, facets, contrasts, intgroup, paths
load(file.path(paths$DEG_output, paste0(params$project_name, "_DEG_data.RData")))
# TODO: take params out of save on other script
#config <- yaml::read_yaml(here::here("config","config.new.yaml"), eval.expr = T)
#params <- config$params
# TODO ok it makes sense for the params to come from the yaml file, not the save. But there needs to be some way for the DESeq script to save stuff that was changed if a new group was generated
# basically the DESeq script shouldn't modify the parameters, that should be read-only.
# I can make some edits to the other branch to facilitate this


# figure out if we're doing one or multiple reports
if(is.na(params$display_group_facet)){
    message("Generating a single report for all data")
    # one report for all data

    # list of DESeq results (indexed by comparison)
    resultsList <- overallResList[['all']]

} else {
    message(paste0("Faceting report by ",params$display_group_facet))
    if(!is.na(params$display_group_filter)){
        facets <- params$display_group_filter
    }else {
        facets <- DESeqDesign %>%
            filter(!(params$display_group_facet) %in% c(params$exclude_groups)) %>%
            pull(params$display_group_facet) %>% 
            unique()
    }
    message(facets)
    for (current_filter in facets) {
        message(paste0("Generating report for ",current_filter))
        # list of DESeq results (indexed by comparison)
        resultsList <- overallResList[[current_filter]]

        output_report_file = "asdf"

        params$current_filter = current_filter


        # subset data
        metadata_subset <- subset_metadata(DESeqDesign, params, contrasts, current_filter)
        metadata <- metadata_subset$DESeqDesign
        contrasts <- metadata_subset$contrasts
        counts <- data.frame() # subset_data(sampleData, metadata)

        print("data subsetted!!!")
        print(contrasts)

        dataFile <- file.path(paths$DEG_output, paste0(params$project_name, "_report_data.RData"))
        save(metadata, contrasts, counts, resultsList, file=dataFile)
        params$dataFile <- dataFile


        rmarkdown::render(input = RMarkdown_file,
                    encoding = "UTF-8",
                    output_file = output_report_file,
                    params = params,
                    envir = new.env())
    }
}