# This script takes in a deduplicated count table (arg1), a matching non-deduplicated count table (arg2), and an output directory
# And calculates the ratio of deduplication for each gene in each sample
# File names must be {library}_umiDedup-{method}.tsv
# Ex. Rscript dedup_stats.R L1306H03_umiDedup-1MM_Directional.tsv L1306H03_umiDedup-No-dedup.tsv output/QC

library(tidyverse)
library(matrixStats)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Argument must be supplied: deduplicated count table, non-deduplicated count table.", call. = FALSE)
}

outdir <- args[3]


# Get info from file name
fname <- basename(args[1])
library <- str_extract(fname, "^[^_]+")
umimethod <- str_extract(fname, "(?<=umiDedup-)[^.]+")

outfile <- paste0(outdir, "/", library, "_", umimethod, "_dedupratios.txt")

print(paste0("Library: ", library))
print(paste0("UMI Method: ", umimethod))

# Import count tables
dedup <- read.delim(args[1], sep = "\t", header = TRUE, row.names = 1)
nodedup <- read.delim(args[2], sep = "\t", header = TRUE, row.names = 1)

# Sanity checks
stopifnot(nrow(dedup) == nrow(nodedup))
stopifnot(ncol(dedup) == ncol(nodedup))
stopifnot(all(rownames(dedup) == rownames(nodedup)))
stopifnot(all(colnames(dedup) == colnames(nodedup)))

# Convert to long format
dedup_long <- dedup %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "dedup_counts")

nodedup_long <- nodedup %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "nodedup_counts")

# Join and calculate ratio
ratio_long <- dedup_long %>%
  inner_join(nodedup_long, by = c("gene", "sample")) %>%
  mutate(
    ratio = nodedup_counts / dedup_counts,
    ratio = ifelse(is.infinite(ratio) | is.nan(ratio), 0, ratio)
  )

# Convert back to wide format 
ratio <- ratio_long %>%
  select(gene, sample, ratio) %>%
  pivot_wider(names_from = sample, values_from = ratio, id_cols = gene) %>%
  column_to_rownames("gene")

# Max ratio of nodedup to dedup per sample
maxratio_stats <- tibble(
  sample = colnames(ratio),
  max_ratio = colMaxs(as.matrix(ratio), na.rm = TRUE)
) %>%
  arrange(sample) %>%
  as.data.frame()

write.table(maxratio_stats, outfile, row.names = FALSE, quote = FALSE, sep = "\t")
