#!/usr/bin/R
# Custom parameters for the report
library('tidyverse')
library('yaml')
library('DESeq2')                                                                                                                                                            

set_up_paths <- function(params) {
    paths <- list()
    # Other important system paths to specify in config
    paths$wikipathways <- params$wikipathways_directory
    # For project structure
    # Should probably update this to use the file.path() function.
    paths$root <- params$projectdir
    paths$data <- file.path(paths$root, "data")
      paths$raw <- file.path(paths$data, "raw")
      paths$processed <- file.path(paths$data, "processed")
      paths$metadata <- file.path(paths$data, "metadata")
    paths$reports <- file.path(paths$root, "reports")
    paths$results <- file.path(paths$root, "results")
    if (is.na(params$group_facet)) {
      paths$DEG_output <- file.path(paths$results, "DEG_output")
    } else {
      paths$DEG_output <- file.path(paths$results, "DEG_output", paste0("group_", paste(params$group_filter, collapse = "_")))
    }
    paths$pathway_analysis <- file.path(paths$DEG_output, "/pathway_analysis")
    paths$RData <- file.path(paths$DEG_output, "/RData")
    paths$BMD_output <- file.path(paths$results, "/DEG_output/BMD_and_biomarker_files")
}


load_species <- function(species){
  species_data = list()
  species_data$loaded <- FALSE
  if (species == "human") {
    # Human:
    library('org.Hs.eg.db')
    species_data$orgdb <- "org.Hs.eg.db"
    species_data$species_sci <- "Homo sapiens"
    species_data$wiki <- "wikipathways-20210810-gmt-Homo_sapiens.gmt"
    species_data$ensembl_species <- "hsapiens_gene_ensembl"
    species_data$species_gene_symbol <- "external_gene_name"
    species_data$kegg_organism <- "hsa"
    species_data$temposeq_manifest <- "181019_Human_S1500_Surrogate_1.2_Manifest.txt"
    species_data$loaded <- TRUE
    # options:
    # "191113_Human_S1500_Surrogate_2.0_Manifest.csv"
    # "181019_Human_S1500_Surrogate_1.2_Manifest.txt"
    # "191004_Human_Whole_Transcriptome_2.0_Manifest.txt"
  } else if (species == "mouse") {
    # Mouse:
    library('org.Mm.eg.db')
    species_data$orgdb <- "org.Mm.eg.db"
    species_data$species_sci <- "Mus musculus"
    species_data$wiki <- "wikipathways-20210810-gmt-Mus_musculus.gmt"
    species_data$ensembl_species <- "mmusculus_gene_ensembl"
    species_data$species_gene_symbol <- "mgi_symbol"
    species_data$kegg_organism <- "mmu"
    species_data$temposeq_manifest <- "181130 Mouse S1500+ Surrogate 1.2 Manifest.xlsx" # or 190603 Mouse Whole Transcriptome 1.1 Manifest.xlsx
    species_data$loaded <- TRUE
  } else if (species == "rat") {
    # Rat: 
    library('org.Rn.eg.db')
    species_data$orgdb <- "org.Rn.eg.db"
    species_data$species_sci <- "Rattus norvegicus"
    species_data$wiki <- "wikipathways-20210810-gmt-Rattus_norvegicus.gmt"
    species_data$ensembl_species <- "rnorvegicus_gene_ensembl"
    species_data$species_gene_symbol <- "rgd_symbol"
    species_data$kegg_organism <- "rno"
    species_data$temposeq_manifest <- "190809 Rat Whole Transcriptome 1.0 Manifest.xlsx"
    species_data$loaded <- TRUE
  } else if (species == "hamster") {
    # Golden hamster:
    species_data$OrgDb.Ma <- query(AnnotationHub(), c("OrgDb", "Mesocricetus auratus"))[[1]]
    species_data$orgdb <- "OrgDb.Ma"
    species_data$species_sci <- "Mesocricetus auratus"
    species_data$ensembl_species <- "mauratus_gene_ensembl"
    species_data$species_gene_symbol <- "external_gene_name"
    species_data$loaded <- TRUE
  } else {
    stop("No species picked!")
  }
  return(species_data)
}


projectdir <- here::here()
print(projectdir)
config <- yaml::read_yaml(file.path(projectdir, "Rmd/config.yml"), eval.expr = T)

# Identify where metadata can be found
SampleKeyFile <- file.path(config$params$projectdir, "metadata/metadata.QC_applied.txt")

# Read in metadata
DESeqDesign <- read.delim(SampleKeyFile,
                          stringsAsFactors = FALSE,
                          sep = "\t",
                          header = TRUE,
                          quote = "\"",
                          row.names = 1) # Column must have unique IDs!!
DESeqDesign$original_names <- rownames(DESeqDesign)

paths <- set_up_paths(params)

# Identify where count data can be found
if (Platform == "TempO-Seq") {
  SampleDataFile <- file.path(paths$processed, "count_table.csv")
  sampledata_sep = ","
} else {
  SampleDataFile <- file.path(paths$processed, "genes.data.tsv")
  sampledata_sep = "\t"
}



if (is.na(config$params$group_facet)) { # all data in one facet
  message(paste0("Making multiple reports based on ",
                 config$params$group_facet ,"..."))
    # load data
    # run DEseq2

} else {
    facets <- DESeqDesign %>%
        filter(!(config$params$group_facet) %in% c(config$params$exclude_groups, skip_extra)) %>%
        filter(!solvent_control) %>%
        pull(config$params$group_facet) %>% 
        unique()

    if (any(!is.na(config$params$group_filter))) { # filter facets
        message(paste0("The group(s) of interest is (are) ",
                 paste(config$params$group_filter, collapse = " and "),".\n",
                 "Writing a single report for that (those) groups."))

    } else { # do all facets separately
        message("Writing a single report for whole experiment.")

    }
}
