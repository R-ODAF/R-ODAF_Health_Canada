#' Load Species Specific Data
#'
#' This function loads information for a specified organism used in the analysis.
#' It requires Bioconductor packages corresponding to the species being analyzed.
#'
#' @param species Target species for the analysis.
#' @param wiki Path or URL to the WikiPathways data.
#' @param manifest Path to the manifest file required for TempO-Seq analysis.
#' @return A list containing various objects pertaining to the species data.
#' @export
load_species <- function(species, wiki, manifest) {
  species_data <- list()
  species_data$loaded <- FALSE
  if (species == "human") {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      stop("Package 'org.Hs.eg.db' is required for human species data.")
    }
    library(org.Hs.eg.db)
    species_data$orgdb <- AnnotationDbi::dbfile(get("org.Hs.eg.db"))
    species_data$species_sci <- "Homo sapiens"
    species_data$wiki <- wiki
    species_data$ensembl_species <- "hsapiens_gene_ensembl"
    species_data$species_gene_symbol <- "external_gene_name"
    species_data$kegg_organism <- "hsa"
    species_data$temposeq_manifest <- manifest
    species_data$loaded <- TRUE
  } else if (species == "mouse") {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      stop("Package 'org.Mm.eg.db' is required for mouse species data.")
    }
    library(org.Mm.eg.db)
    species_data$orgdb <- AnnotationDbi::dbfile(get("org.Mm.eg.db"))
    species_data$species_sci <- "Mus musculus"
    species_data$wiki <- wiki
    species_data$ensembl_species <- "mmusculus_gene_ensembl"
    species_data$species_gene_symbol <- "mgi_symbol"
    species_data$kegg_organism <- "mmu"
    species_data$temposeq_manifest <- manifest
    species_data$loaded <- TRUE
  } else if (species == "rat") {
    if (!requireNamespace("org.Rn.eg.db", quietly = TRUE)) {
      stop("Package 'org.Rn.eg.db' is required for rat species data.")
    }
    library(org.Rn.eg.db)
    species_data$orgdb <-  AnnotationDbi::dbfile(get("org.Rn.eg.db"))
    species_data$species_sci <- "Rattus norvegicus"
    species_data$wiki <- wiki
    species_data$ensembl_species <- "rnorvegicus_gene_ensembl"
    species_data$species_gene_symbol <- "rgd_symbol"
    species_data$kegg_organism <- "rno"
    species_data$temposeq_manifest <- manifest
    species_data$loaded <- TRUE
  } else if (species == "hamster") {
    if (!requireNamespace("AnnotationHub", quietly = TRUE)) {
      stop("Package 'AnnotationHub' is required for golden hamster species data.")
    }
    species_data$OrgDb.Ma <- AnnotationHub::query(AnnotationHub::AnnotationHub(), c("OrgDb", "Mesocricetus auratus"))[[1]]
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
