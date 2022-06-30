#!/usr/bin/env Rscript
library(parallel)
library(QuasR)
library(rtracklayer)
library(GenomicFeatures)

args = commandArgs(trailingOnly=TRUE)

samplefile <- args[1]
genome <- args[2]
annotfile <- args[3]
count_table_file <- args[4]
num_threads   <- as.numeric(args[5])

cl <- makeCluster(num_threads)

# set up cache dir in tempdir
cache_dir <- paste0(tempdir(), "/proba")
dir.create(cache_dir)

message("Running qAlign function...")
proj2 <- qAlign(samplefile,
                genome,
                paired = "no",
                cacheDir = cache_dir,
                clObj = cl)

txStart <- import.gff(annotfile, format = "gtf")
names(txStart) <- txStart@seqnames

cnt <- qCount(proj2, txStart, clObj = cl)

# Make dataframe as count table
cnmat           <- as.data.frame(cnt, header = TRUE)
cnmat$width     <- NULL
write.table(cnmat, count_table_file, sep = '\t', quote = F, row.names = T)
message("Counts Table Completed")
