get_design <- function(design){
  #return(formula(paste0("~", paste0(c(design), collapse = " + "))))
  return(formula(paste0("~",design)))
}


learn_deseq_model <- function(sd, des, intgroup, design, params){
    current_design <- get_design(design)
    dds <- DESeqDataSetFromMatrix(countData = round(sd),
                                    colData   = as.data.frame(des),
                                    design    = current_design)
    if(params$filter_gene_counts){ # filter gene counts to decrease runtime. Not recommended for biomarker input!
        dds <- dds[rowSums(counts(dds)) > 1]
    }
    bpparam <- MulticoreParam(params$cpus)
    dds <- DESeq(dds, parallel = TRUE, BPPARAM = bpparam)
    return(dds)
}

# covariates are used to calculate within-group variability. Nuisance parameters are removed (e.g. for visualization) See:
# http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-after-vst-are-there-still-batches-in-the-pca-plot
regularize_data <- function(dds, design, covariates, nuisance, blind=FALSE){
  if (!is.na(nuisance)) {
    rld <- vst(dds, blind)
    mat <- assay(rld)
    if(!is.na(covariates)){
      condition <- formula(paste0("~",
                                  design,
                                  paste0(covariates[!covariates %in% nuisance],collapse = " + ")))
      
    } else {
      condition <- formula(paste0("~", design))
    }
    mm <- model.matrix(condition, colData(rld))
    mat <- limma::removeBatchEffect(mat, batch = rld[[nuisance]], design = mm)
    assay(rld) <- mat
  } else {
    rld <- vst(dds, blind)
  }
  return(rld)
}

get_DESeq_results <- function(dds, DESeqDesign, contrasts, design, params, current_group_filter, outdir){
    # Initial setup for DESeq2 contrasts
    resListAll <- list()
    resListFiltered <- list()
    resListDEGs <- list()
    Counts  <- counts(dds, normalized = TRUE)
    CPMdds  <- edgeR::cpm(counts(dds, normalized = TRUE))
    mergedDEGs <- c()

    for (x in 1:nrow(contrasts)) {  # For all comparisons to be done  
        condition1 <- contrasts[x, 2] # Control
        condition2 <- contrasts[x, 1] # Experimental

        contrast_string <- paste(condition2, "vs", condition1, sep = "_")

        message(contrast_string)

        DESeqDesign_subset <- as.matrix(DESeqDesign[DESeqDesign[, design] %in% c(condition1, condition2),])
        
        # sanity checks
        stopifnot(exprs = {
            nrow(DESeqDesign_subset) > 0
        })

        # run the DEseq p-value calculation
        message("Obtaining the DESeq2 results")

        currentContrast <- c(design, condition2, condition1)
        bpparam <- MulticoreParam(params$cpus)

        res <- DESeq2::results(dds,
                                parallel = TRUE,
                                BPPARAM = bpparam,
                                contrast = currentContrast,
                                pAdjustMethod = 'fdr',
                                cooksCutoff = params$cooks) # If Cooks cutoff disabled - manually inspect.

        res <- lfcShrink(dds,
                        contrast = currentContrast,
                        res = res,
                        type = "ashr")

        resListAll[[contrast_string]] <- res

        # Filter results using R-ODAF filters
        Filter <- matrix(data = NA, ncol = 3, nrow = nrow(Counts))
        rownames(Filter) <- rownames(Counts) # genes
        colnames(Filter) <- c("Low", "quantile", "spike")
        
        # Apply the "Relevance" condition
        message(paste0("Filtering genes: 75% of at least 1 group need to be above ", params$MinCount, " CPM"))
        SampPerGroup <- table(DESeqDesign_subset[, design])
        
        for (gene in 1:nrow(dds)) {
          CountsPass <- NULL
          for (group in 1:length(SampPerGroup)) { 
            sampleCols <- grep(dimnames(SampPerGroup)[[1]][group], DESeqDesign_subset[, design])
            Check <- sum(CPMdds[gene,sampleCols] >= params$MinCount) >= 0.75 * SampPerGroup[group]
            CountsPass <- c(CountsPass, Check)
          }
          if (sum(CountsPass) > 0) {Filter[gene, 1] <- 1 }  else {Filter[gene,1] <- 0 }
        }
        
        compte <- Counts[Filter[,1] == 1,]
        Filter <- Filter[rownames(Filter) %in% rownames(compte), , drop = F]
        
        message(paste0("Relevance filtering removed ", nrow(dds) - nrow(Filter),
                       " genes from the ", nrow(dds)," assessed. ",
                       nrow(Filter), " genes remaining"))
        
        message("Check median against third quantile" )
        message("AND")
        message("Check the presence of a spike" )

        for (gene in 1:nrow(compte)) {
            # Check the median against third quantile
            quantilePass <- NULL
            sampleColsg1 <- grep(dimnames(SampPerGroup)[[1]][1],DESeqDesign_subset[,design])
            sampleColsg2 <- grep(dimnames(SampPerGroup)[[1]][2],DESeqDesign_subset[,design])

            Check <- median(as.numeric(compte[gene, sampleColsg1])) > quantile(compte[gene, sampleColsg2], 0.75)[[1]]
            quantilePass <- c(quantilePass, Check)

            Check <- median(as.numeric(compte[gene, sampleColsg2])) > quantile(compte[gene, sampleColsg1], 0.75)[[1]]
            quantilePass <- c(quantilePass, Check)

            if (sum(quantilePass) > 0) {
                Filter[gene, 2] <- 1
            }  else {
                Filter[gene, 2] <- 0
            }

            # Check for spikes 
            spikePass <- NULL
            for (group in 1:length(SampPerGroup)) {
                sampleCols <- grep(dimnames(SampPerGroup)[[1]][group], DESeqDesign_subset[ ,design])
                if (max(compte[gene,sampleCols]) == 0) {Check <- FALSE} else {
                Check <- (max(compte[gene, sampleCols]) / sum(compte[gene, sampleCols])) >= 1.4 * (SampPerGroup[group])^(-0.66)
                spikePass <- c(spikePass, Check)
                }
            }
            if (sum(spikePass) > 1) {
                Filter[gene, 3] <- 0
            }  else {
                Filter[gene, 3] <- 1
            }
        }
        
        # Create output tables
        norm_data <- counts(dds[rownames(compte)], normalized = TRUE)
        DEsamples <- subset(res, res$padj < params$alpha & abs(res$log2FoldChange) > log2(params$linear_fc_filter) )
        if (nrow(DEsamples) == 0) {
          message("No significant results were found for this contrast. Moving on...")
          next
        }
        DECounts <- compte[rownames(compte) %in% rownames(DEsamples), , drop = F]
        Filter <- Filter[rownames(Filter) %in% rownames(DECounts), , drop = F]
        

        # Extract the final list of DEGs
        
        allCounts_all_filters <- res[rowSums(Filter) == 3 ,]
        DECounts_real <- DEsamples[rowSums(Filter) == 3 ,]
        DECounts_no_quant <- DEsamples[Filter[, 2] == 0 ,] # save these to output later 
        DECounts_spike <- DEsamples[Filter[, 3] == 0 ,] # save these to output later
        #TODO: output quantile rule failing and spike failing genes

        message(paste0("A total of ", nrow(DECounts_real),
                    " DEGs were selected (out of ", nrow(DECounts) ,"), after ", nrow(DECounts_no_quant),
                    " genes(s) removed by the quantile rule and ", nrow(DECounts_spike),
                    " gene(s) with a spike"))

        message("DESeq2 Done")

        colnames(norm_data) <- colData(dds)[, design]
        
        mergedDEGs <- c(mergedDEGs, rownames(DECounts_real))

        # TODO: the sapply at the end should be handling this, why doesn't it work?
        if (nrow(DECounts_real) > 0){
          resListDEGs[[contrast_string]] <- DECounts_real
        }
        if (nrow(DECounts_real) > 0){
          resListFiltered[[contrast_string]] <- allCounts_all_filters
        }
    }
    # If there are no significant results - remove the empty contrast from the list:
    resListAll <- resListAll[!sapply(resListAll, is.null)]
    resListFiltered <- resListFiltered[!sapply(resListFiltered, is.null)]
    resListDEGs <- resListDEGs[!sapply(resListDEGs, is.null)]
    
    mergedDEGs <- unique(mergedDEGs)
    
    return(list(resListAll=resListAll, resListFiltered=resListFiltered, resListDEGs=resListDEGs, mergedDEGs=mergedDEGs))
}


annotate_deseq_table <- function(deseq_results_list, params, filter_results = F) {
  x <- deseq_results_list
  annotated_results <- list()
  for (i in 1:length(x)) {
    print(i)
    deg_table <- x[[i]]
    # Add taxonomy
    if (is.null(deg_table)) {
      next
    } else if (nrow(deg_table) == 0) {
      next
    } else {
      deg_table <- cbind(Probe_Name = row.names(x[[i]]),
                         as(deg_table, "data.frame"),
                         contrast = gsub(pattern = paste0("log2.*", params$design, "\ "),
                                         replacement =  "",
                                         x = x[[i]]@elementMetadata[[2]][2]))
      deg_table <- dplyr::left_join(deg_table,
                                    params$biospyder,
                                    by = "Probe_Name")
      deg_table <- dplyr::mutate(deg_table, linearFoldChange = ifelse(log2FoldChange > 0,
                                                                      2 ^ log2FoldChange,
                                                                      -1 / (2 ^ log2FoldChange)))
      deg_table <- deg_table %>%
        dplyr::select(Probe_Name, Ensembl_Gene_ID, Gene_Symbol, baseMean, log2FoldChange, linearFoldChange, lfcSE, pvalue, padj, contrast)
      ## FILTERS ##
      if (filter_results == T) {
        deg_table <- deg_table[!is.na(deg_table$padj) & deg_table$padj < alpha & abs(deg_table$linearFoldChange) > linear_fc_filter, ]
      }
      annotated_results[[i]] <- deg_table %>% dplyr::distinct()
    }
  }
  y <- rbindlist(annotated_results)
  return(y)
}
