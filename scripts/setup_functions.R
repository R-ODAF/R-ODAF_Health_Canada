replace_nulls_in_config <- function(params) {
    params_new = list()
    for (name in names(params)) {
    param <- params[[name]]
    if(is.null(param)){
        params_new[[name]] <- NA
    } else {
        params_new[[name]] <- param
    }
    }
    return(params_new)
}

load_species <- function(species){
  species_data = list()
  species_data$loaded <- FALSE
  if (species == "human") {
    # Human:
    library('org.Hs.eg.db')
    species_data$orgdb <- "org.Hs.eg.db"
    species_data$species_sci <- "Homo sapiens"
    species_data$wiki <- params$wikipathways_filename
    species_data$ensembl_species <- "hsapiens_gene_ensembl"
    species_data$species_gene_symbol <- "external_gene_name"
    species_data$kegg_organism <- "hsa"
    species_data$temposeq_manifest <- params$biospyder_manifest_file
    species_data$loaded <- TRUE
  } else if (species == "mouse") {
    # Mouse:
    library('org.Mm.eg.db')
    species_data$orgdb <- "org.Mm.eg.db"
    species_data$species_sci <- "Mus musculus"
    species_data$wiki <- params$wikipathways_filename
    species_data$ensembl_species <- "mmusculus_gene_ensembl"
    species_data$species_gene_symbol <- "mgi_symbol"
    species_data$kegg_organism <- "mmu"
    species_data$temposeq_manifest <- params$biospyder_manifest_file
    species_data$loaded <- TRUE
  } else if (species == "rat") {
    # Rat: 
    library('org.Rn.eg.db')
    species_data$orgdb <- "org.Rn.eg.db"
    species_data$species_sci <- "Rattus norvegicus"
    species_data$wiki <- params$wikipathways_filename
    species_data$ensembl_species <- "rnorvegicus_gene_ensembl"
    species_data$species_gene_symbol <- "rgd_symbol"
    species_data$kegg_organism <- "rno"
    species_data$temposeq_manifest <- params$biospyder_manifest_file
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
    # Normally, as follows: year - project_title - group_filter
    if (is.na(params$group_filter) || is.null(params$group_filter)) {
    analysisID <- paste(format(Sys.time(), '%Y'), params$project_title, sep = "_")
    } else {
    analysisID <- paste(format(Sys.time(), '%Y'),
                        params$project_title,
                        paste(params$group_filter, collapse = "_"),
                        sep = "_")
    }
    return(analysisID)
}

load_biospyder_new <- function(biospyder_dbs, temposeq_manifest){
  return_data = list()
  biospyder <- read.delim(file.path(biospyder_dbs, temposeq_manifest), # Assay manifest...
                          stringsAsFactors = FALSE,
                          sep = ",",
                          header = TRUE,
                          quote = "\"")
  
  feature_id <- "Probe_Name"
  biospyder[ biospyder == "NULL" ] <- NA
  return_data$biospyder <- biospyder
  return_data$feature_id <- feature_id
  return(return_data)
}


set_up_platform_params <-function(params){
    species_data <- params$species_data
    SampleDataFile <- file.path(paths$processed, "count_table.tsv")
    params$SampleDataFile <- SampleDataFile
    params$sampledata_sep = "\t"
  # set some additional parameters based on platform
  if (params$platform == "RNA-Seq") {
    params$threshold <- 1000000
    params$MinCount <- 1
    params$alpha <- pAdjValue <- 0.05 # Relaxed from 0.01
    params$linear_fc_filter <- 1.5
    params$feature_id <- "Ensembl_Gene_ID"
  } else if (params$platform == "TempO-Seq") {
    params$threshold = 100000
    params$MinCount <- 0.5
    params$alpha <- pAdjValue <- 0.05 
    params$linear_fc_filter <- 1.5
    
    bs <- load_biospyder_new(params$biospyder_dbs, species_data$temposeq_manifest)
    params$feature_id <- bs$feature_id # Probe_Name
    params$biospyder <- bs$biospyder # manifest
  } else { 
    stop("Platform/technology not recognized") 
  }
  return(params)
}