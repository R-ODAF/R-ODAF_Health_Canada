library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Argument must be supplied: full path to input count table, full path to desired output file.", call. = FALSE)
}

infile <- args[1]
outfile <- args[2]

df <- read.delim(infile, sep = "\t", row.names="X", header = TRUE)


df_ercc <- df %>%
  rownames_to_column("gene") %>%
  filter(grepl("ERCC", gene)) %>%
  column_to_rownames("gene")


# Calculate sums of reads per sample
ercc_sums <- colSums(df_ercc)
total_sums <- colSums(df)

# Create a dataframe with the results
ercc_stats <- tibble(
  sample_id = names(total_sums),
  ercc_reads = as.numeric(ercc_sums),
  total_reads = as.numeric(total_sums),
  ercc_proportion = ercc_reads / total_reads
) %>% as.data.frame() %>%
  dplyr::arrange(sample_id)

write.table(ercc_stats, outfile, quote = FALSE, sep = "\t", row.names = FALSE)
