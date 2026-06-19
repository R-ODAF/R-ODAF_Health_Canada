library(data.table)
library(Matrix)
library(tools)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Arguments must be supplied in this order: (1) full path to mtx file (2) Library name (3) path to table of barcodes and sampleIDs", call. = FALSE)
}

infile <- args[1]
library <- args[2]
barcodefile <- args[3]
matrix_dir <- dirname(infile)
fname <- file_path_sans_ext(basename(infile))
outfile <- paste0("output/", library, "/", library, "_", fname, ".tsv") 

print(paste0("Importing ", infile))

f <- file(infile, "r")
mat <- as.data.frame(as.matrix(readMM(f)))
close(f)
feature.names = fread(paste0(matrix_dir, "/features.tsv"), header = FALSE, stringsAsFactors = FALSE, data.table = F)
barcode.names = fread(paste0(matrix_dir, "/barcodes.tsv"), header = FALSE, stringsAsFactors = FALSE, data.table = F)
colnames(mat) <- barcode.names$V1
rownames(mat) <- feature.names$V1

# Replace barcodes with sampleIDs in column names
md <- fread(barcodefile, header = TRUE, stringsAsFactors = FALSE, data.table = F)
barcode_to_sample <- setNames(md$sample_id, md$barcode)
colnames(mat) <- barcode_to_sample[colnames(mat)]

# Filter only for samples in metadata
# To get rid of columns for barcodes that weren't used in this library
samples <- md$sample_id
mat <- mat[, colnames(mat) %in% samples]

print(paste0("Creating ", outfile))
fwrite(mat, file = outfile, sep = "\t", quote = F, row.names = T, col.names = T)