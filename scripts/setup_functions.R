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
    lapply(paths, function(x) if(!dir.exists(x)) dir.create(x))
    return(paths)
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

get_analysis_id <- function(params){
    # Set analysis ID. This ID will be used as prefix for the output files
    # Normally, as follows: year - project_name - group_filter
    if (is.na(params$group_filter)) {
    analysisID <- paste(format(Sys.time(), '%Y'), params$project_name, sep = "_")
    } else {
    analysisID <- paste(format(Sys.time(), '%Y'),
                        params$project_name,
                        paste(params$group_filter, collapse = "_"),
                        sep = "_")
    }
    return(analysisID)
}

filter_metadata <- function(DESeqDesign, params){
    # exclude samples
    if (any(!is.na(params$exclude_samples))) {
        DESeqDesign <- DESeqDesign %>% 
            dplyr::filter(!original_names %in% params$exclude_samples)
    }
    # exclude groups
    if (any(!is.na(params$exclude_groups))) {
        DESeqDesign <- DESeqDesign %>%
            dplyr::filter(!(!!sym(params$design)) %in% params$exclude_groups)
        contrasts_to_filter <- DESeqDesign %>% 
            dplyr::filter(!(!!sym(params$design)) %in% params$exclude_groups) %>%
            pull(params$design) %>% 
            unique()
        contrasts <- contrasts %>%
            dplyr::filter(V1 %in% contrasts_to_filter)
        if (params$strict_contrasts == T) {
            contrasts <- contrasts %>%
                dplyr::filter(V2 %in% contrasts_to_filter)
        }
    }
    if (!is.na(params$include_only_column) & !is.na(params$include_only_group)) {
        DESeqDesign <- DESeqDesign %>%
            dplyr::filter((!!sym(params$include_only_column)) %in% params$include_only_group)
        limit_contrasts <- DESeqDesign %>%
            pull(!!sym(params$design)) %>%
            unique() %>%
            as.character()
        contrasts <- contrasts %>% dplyr::filter(V1 %in% limit_contrasts)
    }
    return(DESeqDesign)
}

# TODO: this might need work. Conversion to factors might require sorting?
sort_metadata <- function(DESeqDesign, contrasts, params){
    # reorder contrast list by the specified column
    ordered_metadata <- DESeqDesign[mixedorder(DESeqDesign[,params$sortcol]),]
    ordered_design <- ordered_metadata %>%
        dplyr::select(params$design) %>%
        dplyr::pull()
    contrasts <- contrasts %>%
        dplyr::slice(match(ordered_design, V1)) %>%
        unique()
    return(list(design=ordered_metadata,contrasts=contrasts))
}

load_count_data <- function(SampleDataFile, sampledata_sep){
  sampleData <- read.delim(SampleDataFile,
                         sep = sampledata_sep,
                         stringsAsFactors = FALSE,
                         header = TRUE, 
                         quote = "\"",
                         row.names = 1,
                         check.names = FALSE)
  return(sampleData)
}

load_biospyder <- function(biospyder_dbs, temposeq_manifest){
  return_data = list()
  biospyder <- read.delim(file.path(biospyder_dbs, temposeq_manifest), # Assay manifest...
                        stringsAsFactors = FALSE,
                        sep = "\t",
                        header = TRUE,
                        quote = "\"")
  # Annoyingly, the manifests are different depending on platform and version.
  if (colnames(biospyder)[1] == "PROBE_NAME") {
    biospyder_ID = "ENSEMBL_GENE_ID"
    biomart_filter <- "PROBE_NAME"
    biospyder_filter = "ensembl_gene_id"
  } else if (colnames(biospyder)[1] == "Probe.name") {
    biospyder_ID = "Reference.Transcript"
    biomart_filter = "Probe.name"
    biospyder_filter = "refseq_mrna"
  }
  # Fill in the blanks for TempO-Seq Manifest...
  # Set NULL values to NA
  biospyder[ biospyder == "NULL" ] <- NA
  return_data$biospyder <- biospyder
  return_data$biospyder_ID <- biospyder_ID
  return_data$biomart_filter <- biomart_filter
  return_data$biospyder_filter <- biospyder_filter
  return(return_data)
}