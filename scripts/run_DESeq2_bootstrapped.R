# Load libraries #############
library(readxl)
library(tidyverse)
library(DESeq2)
library(edgeR)
library(ashr)
library(gtools)

res.design <- "group"
js_predictor_variable <- "group"
js_model_formula <- formula(paste0("~",js_predictor_variable))
js_all_variables <- paste0(js_predictor_variable)

nboot <- 1000
cell_lines <- c("Liver spheroids", "TK6", "HepaRG", "MCF-7")
n_vector <- 3:10
boot_vector <- 1:nboot

covariates <- NA
params <- list()
params$cpus <- 40
params$alpha <- 0.05
params$MinCount <- 0.5
params$linear_fc_filter <- 1.5

# Set R-ODAF parameters ############
params$nmr_threshold <- 100000
params$cooks <- FALSE
params$filter_gene_counts <- FALSE
params$parallel <- FALSE
params$R_ODAF_plots <- FALSE

# Import/create files ##############
metadata <- read.delim("data/metadata/metadata.QC_applied.txt")
metadata <- metadata %>%
  dplyr::filter(chemical == "DMSO") %>%
  dplyr::select(-group)
count_data <- read.delim("data/processed/count_table.tsv")
#count_data[is.na(count_data) ] <- 0
#count_data <- count_data[,(colSums(count_data) > params$nmr_threshold)]

DESeq_outputs <- tibble(rownames(count_data)) %>%
  rename("rownames(count_data)"="probe_id") # Odd name, but there is a reason for it
contrasts <- as_tibble(matrix(c("A", "B"), nrow=1, ncol=2), .name_repair="minimal")


# Set up sampling #############
# Create lists
ls_run_name <- list()
ls_run_number <- list()
#ls_metadata <- list()
ls_metadata <- readRDS("bootstrap_DEGs/output/1000_iterations/nobatch_2023-05-10/ls_metadata.RDS")

# Select a subset of samples for each iteration, and create the metadata file
for(cell_line_num in seq_along(cell_lines)){
  cell_line_name <- cell_lines[cell_line_num]
  for (j in seq_along(n_vector)) {
    group1_n <- n_vector[j]
    group2_n <- group1_n
    i <- 1
    while (i < nboot+1) {
      # Create integer for this iteration of the loop
      run_number <- i+((j-1)*length(boot_vector))+((cell_line_num-1)*length(boot_vector)*length(n_vector))
      print(run_number)
      # Sample
      # fastq_id <- metadata %>% 
      #   dplyr::filter(cell_line == cell_line_name) %>% 
      #   dplyr::select(Sample_ID) %>% 
      #   as.vector() %>% 
      #   unlist()
      #sample_id <- unname(sample(fastq_id, size=(group1_n+group2_n), replace=FALSE))
      #group <- c(rep("A",group1_n), rep("B",(group2_n)))
      #metadata_suba <- data.frame(sample_id, group)
      #metadata_sub <- metadata_suba %>% 
      #  left_join(metadata, by = c("sample_id" = "Sample_ID")) %>% 
      #  dplyr::select(sample_id, group)
      #metadata_sub$group <- as.factor(metadata_sub$group)
      ls_run_number[[run_number]] <- run_number
      ls_run_name[[run_number]] <- paste0(cell_lines[cell_line_num], "_", group1_n, "_", group2_n, "_", boot_vector[i])
      print(ls_run_name[[run_number]])
      #ls_metadata[[run_number]] <- metadata_sub
      #names(ls_metadata) <- ls_run_name
      i <- i+1}
    }
}

js_fun_everything <- function(i,
                              ls_metadata = ls_metadata_use, 
                              count_data = countdata_use, 
                              params = params_use,  
                              intgroup = js_predictor_variable_use,
                              contrasts = contrasts_use,
                              js_model_formula = js_model_formula_use, 
                              ls_run_name = ls_run_name_use) {
  exp_metadata <- ls_metadata[[i]]
  exp_metadata$original_names <- exp_metadata$sample_id

  filtered_table <- data.frame()
  
  count_data <- count_data[,(colSums(count_data) > params$nmr_threshold)] # reads required per sample
  exp_metadata <- exp_metadata[exp_metadata$original_names %in% colnames(count_data),]
  count_data <- count_data[,exp_metadata$original_names]
  
  # Intgroups need to be factors for DESeq2
  # make sure the levels are sorted for plotting later
  for(int in intgroup){
    exp_metadata[int] <- factor(exp_metadata[[int]], levels=mixedsort(unique(exp_metadata[[int]])))
  }
  
  check_data <- function(sd, des, con){
    message("Sanity checks for data")
    # make sure they're not empty
    stopifnot(exprs = {
      ncol(sd) > 0
      nrow(des) > 0
      nrow(con) > 0
    })
    # Sanity check: each sample (row) in the metadata should have a corresponding column in the count data
    stopifnot(all(des$original_names %in% colnames(sd)))
    # Sanity check: each column in the count data should have a corresponding sample (row) in the metadata
    stopifnot(all(colnames(sd) %in% des$original_names))
    #message("All OK 👍")
  }
  
  check_data(count_data, exp_metadata, contrasts)
  processed <- list(count_data=count_data, exp_metadata=exp_metadata, contrasts=contrasts)
  
  count_data <- processed$count_data
  exp_metadata <- processed$exp_metadata
  contrasts <- processed$contrasts 
  
  num_rows_metadata <- nrow(exp_metadata)
  
  dds <- DESeqDataSetFromMatrix(countData = round(count_data),
                                colData   = as.data.frame(exp_metadata),
                                design    = js_model_formula)
  
  if(params$filter_gene_counts){ # filter gene counts to decrease runtime. Not recommended for biomarker input!
    dds <- dds[rowSums(counts(dds)) > 1]
  }
  dds <- dds[rowSums(counts(dds)) > 1]
  
  dds <- DESeq(dds,
               parallel = params$parallel)
               #BPPARAM = bpparam)

  filtered_table <- data.frame()
  Counts  <- counts(dds, normalized = TRUE)
  CPMdds  <- edgeR::cpm(counts(dds, normalized = TRUE))
  current_group_filter <- NA
  
  mergedDEGs <- c()
  allCounts_all_filters <- data.frame()
  DECounts_real <- data.frame()
  js_filtered_out <- list()
  
  
  for (x in 1:nrow(contrasts)) {  # For all comparisons to be done  
    condition1 <- unname(unlist(contrasts[1, 2])) # Control
    condition2 <- unname(unlist(contrasts[1, 1])) # Experimental
    contrast_string <- paste(condition2, "vs", condition1, sep = "_")
    message(contrast_string)
    exp_metadata_subset <- as.matrix(exp_metadata)
    
    # sanity checks
    stopifnot(exprs = {
      nrow(exp_metadata_subset) > 0
    })
    
    # Filter results using R-ODAF filters
    Filter <- matrix(data = NA, ncol = 3, nrow = nrow(Counts))
    rownames(Filter) <- rownames(Counts) # genes
    colnames(Filter) <- c("Low", "quantile", "spike")
    
    # Apply the "Relevance" condition
    message(paste0("Filtering genes: 75% of at least 1 group need to be above ", params$MinCount, " CPM"))
    SampPerGroup <- table(exp_metadata_subset[, js_predictor_variable])
    if (!SampPerGroup[condition2] > 1) { next }
    
    for (gene in 1:nrow(dds)) {
      CountsPass <- NULL
      for (group in 1:length(SampPerGroup)) {
        sampleCols <- grep(names(SampPerGroup)[group], exp_metadata_subset[, js_predictor_variable], fixed = T)
        sampleNamesGroup <- exp_metadata_subset[sampleCols, "original_names"]
        Check <- sum(CPMdds[gene,sampleNamesGroup] >= params$MinCount) >= 0.75 * SampPerGroup[group]
        CountsPass <- c(CountsPass, Check)
      }
      if (sum(CountsPass) > 0) {Filter[gene, 1] <- 1 }  else {Filter[gene,1] <- 0 }
    }
    
    compte <- Counts[Filter[,1] == 1,]
    rm_MinCount <- Counts[Filter[,1] != 1,]
    Filter <- Filter[rownames(Filter) %in% rownames(compte), , drop = F]
    
    intitial_count <- nrow(dds)
    num_relevance_filtered <- nrow(dds) - nrow(Filter)
    message(paste0("Relevance filtering removed ", num_relevance_filtered,
                   " genes from the ", nrow(dds)," assessed. ",
                   nrow(Filter), " genes remaining"))
    # run the DEseq p-value calculation
    message("Obtaining the DESeq2 results")
    currentContrast <- c(res.design, condition2, condition1)
    
    
    res <- DESeq2::results(dds[rownames(compte),],
                           parallel = params$parallel,
                           #BPPARAM = bpparam,
                           contrast = currentContrast,
                           alpha = params$alpha,
                           pAdjustMethod = 'fdr',
                           cooksCutoff = params$cooks) # If Cooks cutoff disabled - manually inspect.
    
    res <- lfcShrink(dds,
                     contrast = currentContrast,
                     parallel = params$parallel,
                     #BPPARAM = bpparam,
                     res = res,
                     type = "ashr")
    dds <- NULL
    gc()
    
    DEsamples <<- subset(res, res$padj < params$alpha)
    if (nrow(DEsamples) == 0) {
      print("No significant results were found for this contrast. Moving on...")
      next
    }
    DECounts <- compte[rownames(compte) %in% rownames(DEsamples), , drop = F]
    Filter <- Filter[rownames(Filter) %in% rownames(DECounts), , drop = F]
    # message("Check median against third quantile" )
    # message("AND")
    # message("Check the presence of a spike" )
    
    for (gene in 1:nrow(DECounts)) {
      # Check the median against third quantile
      quantilePass <- NULL
      sampleColsg1 <- grep(dimnames(SampPerGroup)[[1]][1],exp_metadata_subset[,js_predictor_variable], fixed = T)
      sampleColsg2 <- grep(dimnames(SampPerGroup)[[1]][2],exp_metadata_subset[,js_predictor_variable], fixed = T)
      
      # Same problem as above, use names explictly
      sampleNames_g1 <- exp_metadata_subset[sampleColsg1, "original_names"]
      sampleNames_g2 <- exp_metadata_subset[sampleColsg2, "original_names"]
      
      Check <- median(as.numeric(DECounts[gene, sampleNames_g1])) > quantile(DECounts[gene, sampleNames_g2], 0.75)[[1]]
      quantilePass <- c(quantilePass, Check)
      
      Check <- median(as.numeric(DECounts[gene, sampleNames_g2])) > quantile(DECounts[gene, sampleNames_g1], 0.75)[[1]]
      quantilePass <- c(quantilePass, Check)
      
      if (sum(quantilePass) > 0) {
        Filter[gene, 2] <- 1
      }  else {
        Filter[gene, 2] <- 0
      }
      
      # Check for spikes 
      spikePass <- NULL
      for (group in 1:length(SampPerGroup)) {
        sampleColsSpike <- grep(dimnames(SampPerGroup)[[1]][group], exp_metadata_subset[ ,js_predictor_variable], fixed = T)
        sampleNamesSpike <- exp_metadata_subset[sampleColsSpike, "original_names"]
        if (max(DECounts[gene,sampleColsSpike]) == 0) {Check <- FALSE} else {
          Check <- (max(DECounts[gene, sampleNamesSpike]) / sum(DECounts[gene, sampleNamesSpike])) >=
            1.4 * (SampPerGroup[group])^(-0.66)
          spikePass <- c(spikePass, Check)
        }
      }
      if (sum(spikePass) > 1) {
        Filter[gene, 3] <- 0
      }  else {
        Filter[gene, 3] <- 1
      }
    }
    
    #message(paste0("Filtering by linear fold-change: linear FC needs to be above ", params$linear_fc_filter))
    
    allCounts_all_filters <- res[rowSums(Filter) == 3 ,]
    DECounts_real <- DEsamples[rowSums(Filter) == 3 & !is.na(DEsamples$padj) &  abs(DEsamples$log2FoldChange) > log2(params$linear_fc_filter) ,]
  
    DECounts_no_quant <- DEsamples[Filter[, 2] == 0 ,] # save these to output later 
    DECounts_spike <- DEsamples[Filter[, 3] == 0 ,] # save these to output later
    
    js_filtered_out <- list(rm_MinCount, DECounts_no_quant, DECounts_spike)
    
    #TODO: output quantile rule failing and spike failing genes
    
    # message(paste0("A total of ", nrow(DECounts_real),
    #                " DEGs were selected (out of ", nrow(DECounts) ,"), after ", nrow(DECounts_no_quant),
    #                " genes(s) removed by the quantile rule, ", nrow(DECounts_spike),
    #                " gene(s) with a spike, and linear fold-change filtering was applied"))
    # message("DESeq2 Done")
  
    mergedDEGs <- c(mergedDEGs, rownames(DECounts_real))
    
    filtered_table <- rbind(filtered_table, data.frame(facet = current_group_filter,
                                                       contrast = contrast_string,
                                                       initial = intitial_count,
                                                       relevance_filtered = num_relevance_filtered,
                                                       quantile_filtered = nrow(DECounts_no_quant),
                                                       spike_filtered = nrow(DECounts_spike),
                                                       passed_all_filters = nrow(DECounts_real)))
    
  }
  
  mergedDEGs <- unique(mergedDEGs)

  output <- list(#dds = dds,
                 DESeqRes = res,
                 DESeqResFiltered = allCounts_all_filters,
                 DESeqResDEGs = DECounts_real,
                 DESeqmergedDEGs = mergedDEGs,
                 DESeqfiltered_table = filtered_table,
                 DESeqjs_filtered_out = js_filtered_out,
                 exp_metadata = exp_metadata,
                 contrasts = contrasts,
                 intgroup = intgroup,
                 num_rows_metadata = num_rows_metadata,
                 ls_run_name[[i]])
  return(output)
  gc()
}

ls_metadata_use <- ls_metadata
countdata_use <- count_data
params_use <- params
js_predictor_variable_use <- js_predictor_variable
contrasts_use <- contrasts
js_model_formula_use <- js_model_formula
res.design_use <- res.design
ls_run_name_use <- ls_run_name

path <- paste0("bootstrap_DEGs/output/",nboot,"_iterations/nobatch_",Sys.Date())
dir.create(path=path,recursive = T)

saveRDS(ls_metadata, file=paste0(path, "/ls_metadata.RDS"))

spheroids <- parallel::mclapply(X = ls_run_number[1:8000],
                                FUN = js_fun_everything,
                                mc.cores = 40)
names(spheroids) <- ls_run_name[1:8000]
saveRDS(spheroids, file=paste0(path, "/spheroids.RDS"))

tk6 <- parallel::mclapply(X = ls_run_number[8001:16000],
                          FUN = js_fun_everything,
                          mc.cores = 40)
names(tk6) <- ls_run_name[8001:16000]
saveRDS(tk6, file=paste0(path, "/tk6.RDS"))

heparg <- parallel::mclapply(X = ls_run_number[16001:24000],
                          FUN = js_fun_everything,
                          mc.cores = 40)
names(heparg) <- ls_run_name[16001:24000]
saveRDS(heparg, file=paste0(path, "/heparg.RDS"))

mcf7 <- parallel::mclapply(X = ls_run_number[24001:32000],
                             FUN = js_fun_everything,
                             mc.cores = 40)
names(mcf7) <- ls_run_name[24001:32000]
saveRDS(heparg, file=paste0(path, "/mcf7.RDS"))
