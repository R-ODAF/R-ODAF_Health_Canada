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
    species_data$wiki <- params$wikipathways_filename
    species_data$ensembl_species <- "mmusculus_gene_ensembl"
    species_data$species_gene_symbol <- "mgi_symbol"
    species_data$kegg_organism <- "mmu"
    species_data$temposeq_manifest <- "181130 Mouse S1500+ Surrogate 1.2 Manifest.txt" # or 190603 Mouse Whole Transcriptome 1.1 Manifest.xlsx
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
  } else if (colnames(biospyder)[1] == "PROBE_ID") {
    # This should be an "or" statement instead... these manifests are so inconsistent.
    # Same general idea as the first case, but slight differences in the Excel file
    # Had to alter PROBE_NAME manually to include the number as XYZ_1
    biospyder_ID = "ENSEMBL_GENE_ID"
    biomart_filter <- "PROBE_NAME"
    biospyder_filter = "ensembl_gene_id"
  } else if (colnames(biospyder)[1] == "Probe.Name" && colnames(biospyder)[2] == "Gene.Symbol") {
    biospyder_ID = "Gene.Symbol"
    biomart_filter <- "Gene.Symbol"
    biospyder_filter = "ensembl_gene_id"
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