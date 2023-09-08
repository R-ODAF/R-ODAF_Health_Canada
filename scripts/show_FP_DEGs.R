# Load count data
count_data <- read.table("data/processed/count_table.tsv", sep = "\t", header = T)

#spheroid_DEGs <- readRDS("bootstrap_DEGs/output/1000_iterations/nobatch_2023-05-10/spheroids.RDS")
#tk6_DEGs <- readRDS("bootstrap_DEGs/output/1000_iterations/nobatch_2023-05-10/tk6")
#mcf7_DEGs <- readRDS("bootstrap_DEGs/output/1000_iterations/nobatch_2023-05-10/mcf7.RDS")
#heparg_DEGs <- readRDS("bootstrap_DEGs/output/1000_iterations/nobatch_2023-05-10/heparg.RDS")



# Find interesting examples of false positive DEGs.
poor_experiments <- all_cell_lines %>%
  dplyr::filter(nDEGs>10, n_per_group == 4, Filter == "n_final_DEGs")

# examples...
# TK6_10_10_27 TK6_10_10_529

# Really bad ones at n = 4
# TK6_4_4_136
# MCF-7_4_4_59
# Liver spheroids_4_4_712
# HepaRG_4_4_582

bad_experiments <- list(spheroid_DEGs$`Liver spheroids_4_4_712`,
                        tk6_DEGs$`TK6_4_4_136`,
                        mcf7_DEGs$`MCF-7_4_4_59`,
                        heparg_DEGs$`HepaRG_4_4_582`)


# Get DEGs from the bad experiments
# Return a list that can be used as function input
get_bad_DEGs <- function(bad_exp) {
  bad_DEG_table <- as.data.frame(
    bad_exp$DESeqResFiltered
  )
  bad_DEGs <- bad_DEG_table |>
    dplyr::arrange(desc(abs(bad_DEG_table$log2FoldChange))) %>%
    dplyr::filter(padj<=0.05) %>% head() %>% row.names()
  return(bad_DEGs)
}


# Search for a gene name by pattern...
rownames(count_data)[grep(pattern = "CYP", x = rownames(count_data))]


plot_FDR_gene_counts_by_group <- function(gene_to_plot = "CYP1A1_1698",
                                          cell_line_to_plot = "HepaRG",
                                          mdata_subset = heparg_DEGs$HepaRG_3_3_1$exp_metadata,
                                          chemical = "DMSO",
                                          mdata = metadata
                                          ) {
  sample_ids_to_plot <- mdata %>%
    dplyr::filter(cell_line == cell_line_to_plot & chemical == chemical) %>%
    pull(original_names)
  count_data_cells <- count_data %>% dplyr::select(sample_ids_to_plot)
  plot_data <- count_data_cells[gene_to_plot,] %>%
    pivot_longer(cols = everything(),
                 names_to = "sample",
                 values_to = "count")
  plot_data$gene <- gene_to_plot
  
  samples_A <- mdata_subset %>%
    dplyr::filter(group == "A") %>% dplyr::pull(sample_id)
  samples_B <- mdata_subset %>%
    dplyr::filter(group == "B") %>% dplyr::pull(sample_id)
  
  plot_data$group <- ifelse(plot_data$sample %in% samples_A, "A",
                            ifelse(plot_data$sample %in% samples_B, "B",
                                   "All samples"))
  
  # Identify rows with group equal to a group name
  rows_to_duplicate <- plot_data$group != "All samples"
  
  # Create a new data frame with duplicated rows and group set to NA
  df_duplicate <- plot_data[rows_to_duplicate, ]
  df_duplicate$group <- "All samples"
  
  # Bind the duplicated rows with the original data frame
  result_df <- rbind(plot_data, df_duplicate)
  result_df$group <- factor(result_df$group, levels = c("All samples","B","A"))
  
  p <- ggplot(result_df, aes(x = gene, y = count, color = group)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.2)) +
    scale_y_log10() +
    facet_grid(~group) +
    ggtitle(paste0(gene_to_plot,", ",cell_line_to_plot)) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank()) 
  return(p)
}

# False positive TK6
plot_FDR_gene_counts_by_group("WBP1_12591", "TK6", tk6_DEGs$TK6_10_10_27$exp_metadata)
# Same gene in TK6 but different dataset
plot_FDR_gene_counts_by_group("WBP1_12591", "TK6", tk6_DEGs$TK6_10_10_28$exp_metadata)

plot_FDR_gene_counts_by_group("WBP1_12591", "HepaRG", heparg_DEGs$HepaRG_10_10_2$exp_metadata)

# The worst experiments for each cell line:
# TK6
FP_deg_plots <- list()
for (i in seq_along(bad_experiments)) {
  exp_name <- bad_experiments[[i]][[11]]
  current_cell_line <- strsplit(exp_name, "_")[[1]][[1]]
  current_bad_degs <- get_bad_DEGs(bad_experiments[[i]])
  current_plots <- list()
  for (j in seq_along(current_bad_degs)) {
    current_plot <- plot_FDR_gene_counts_by_group(current_bad_degs[[j]], current_cell_line, bad_experiments[[i]]$exp_metadata)
    current_plots[[j]] <- current_plot
  }
  FP_deg_plots[[i]] <- current_plots
}


# What if we look at genes identified in FPs in one experiment - look across cell lines and experiments

# Do experiments with many numerous DEGs have something different compared to experiments with few?

# Is there something weird (in common) about the samples included in the bad experiments??
 # Do the same genes pop up within a cell model as FP DEGs? Across models? What about other correlations - GC content?
