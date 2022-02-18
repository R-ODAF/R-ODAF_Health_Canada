#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

samplefile <- args[1]
genome <- args[2]
annotfile <- args[3]
count_table_file <- args[4]
mapped_unmapped_file <- args[5]

# set up cache dir in tempdir
cache_dir <- paste0(tempdir(),"/proba")
dir.create(cache_dir)


# Import the QuasR library
library(QuasR)

message("Running qAlign function...")
proj2 <- qAlign(samplefile,
                genome,
                paired="no",
                cacheDir=cache_dir)

message("Loading rtracklayer and GenomicFeatures...")
library(rtracklayer)
library(GenomicFeatures)

txStart <- import.gff(annotfile, format="gtf")
names(txStart) <- txStart@seqnames

cnt <- qCount(proj2, txStart)

# Make dataframe as count table
cnmat           <- as.data.frame(cnt, header=TRUE)
#colnames(cnmat) <- gsub("fastqs/", "", colnames(cnmat))
cnmat$width     <- NULL
write.csv(cnmat, count_table_file)
cat("Counts Table Completed\n")

# Get Mapped/Unmapped Counts and Output
a_frame     <- data.frame(alignmentStats(proj2), header=TRUE, check.names=FALSE) 
unmapped    <- a_frame$unmapped 
mapped      <- colSums(cnmat) 
two         <- rbind(mapped) 
three       <- data.frame(two, check.names=FALSE)     
three$width <- NULL 
ll          <- rbind(unmapped=unmapped, mapped=three)

write.csv(as.matrix(ll), mapped_unmapped_file)
message("Mapped/Unmapped Table Completed\n")