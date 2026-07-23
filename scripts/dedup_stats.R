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

# Import count tables
dedup <- read.delim(args[1], sep = "\t", header = TRUE, row.names = 1)
nodedup <- read.delim(args[2], sep = "\t", header = TRUE, row.names = 1)

# Get outdir from arguments
outdir <- args[3]

# Get library and umi dedup method from arguments
library <- args[4]
umimethod <- args[5]

print(paste0("Library: ", library))
print(paste0("UMI Method: ", umimethod))

# Set up output files

outfile_max <- paste0(outdir, "/", library, "_umi", umimethod, "_dedup_max_ratios.txt")
outfile_hist <- paste0(outdir, "/", library, "_umi", umimethod, "_dedup_ratios_hist.pdf")
outfile_sums <- paste0(outdir, "/", library, "_umi", umimethod, "_dedup_countsums.txt")


# Sanity checks
stopifnot(nrow(dedup) == nrow(nodedup))
stopifnot(ncol(dedup) == ncol(nodedup))
stopifnot(all(rownames(dedup) == rownames(nodedup)))
stopifnot(all(colnames(dedup) == colnames(nodedup)))

################################################################################
# Compare total read counts
################################################################################
# Sum counts for each sample
df_sums <- data.frame(colSums(dedup), colSums(nodedup)) %>%
  dplyr::rename(sumreads_dedup = colSums.dedup.,
                sumreads_nodedup = colSums.nodedup.) %>%
  dplyr::mutate(ratio = sumreads_dedup/sumreads_nodedup,
                diff = sumreads_nodedup - sumreads_dedup)

write.table(df_sums, outfile_sums, row.names = TRUE, quote = FALSE, sep = "\t")
print(paste0("Created file: ", outfile_sums))

################################################################################
# Calculate dedup ratios
################################################################################
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
  ) %>%
  dplyr::filter(ratio > 0)


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

write.table(maxratio_stats, outfile_max, row.names = TRUE, quote = FALSE, sep = "\t")
print(paste0("Created file: ", outfile_max))

# Plot deduplication ratios

p <- ggplot(ratio_long, aes(x = ratio)) +
  geom_histogram(binwidth = 1) +
  scale_y_log10() +
  theme_bw() +
  labs (title = paste0("Dedup Ratio Histogram. Library: ", library, ". UMI method: ", umimethod),
        y = "N mapped genes (log10 scale)",
        x = "Deduplication ratio (non-dedup count / dedup count)")

ggsave(outfile_hist, p)
print(paste0("Created file: ", outfile_hist))
