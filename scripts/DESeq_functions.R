get_design <- function(design){
  #return(formula(paste0("~", paste0(c(design), collapse = " + "))))
  if(!is.null(params$formula_override)) {
    des <- formula(paste0("~", params$formula_override))
  } else { des <- formula(paste0("~", design)) }
  return(des)
}


learn_deseq_model <- function(sd, des, design, params){
  current_design <- get_design(design)
  dds <- DESeqDataSetFromMatrix(countData = round(sd),
                                colData   = as.data.frame(des),
                                design    = current_design)
  
  if(params$filter_gene_counts){ # filter gene counts to decrease runtime. Not recommended for biomarker input!
    dds <- dds[rowSums(counts(dds)) > 1]
  }
  bpparam <- MulticoreParam(params$cpus)
  dds <- dds[rowSums(counts(dds)) > 1]
  dds <- DESeq(dds, parallel = TRUE, BPPARAM = bpparam)
  return(dds)
}

# covariates are used to calculate within-group variability. Batch is considered a nuisance parameter and is removed (e.g. for visualization) See:
# http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-after-vst-are-there-still-batches-in-the-pca-plot
regularize_data <- function(dds, design, covariates, batch_param, blind=FALSE){
  if (!is.na(batch_param)) {
    rld <- vst(dds, blind)

    mat <- assay(rld)
    if(!is.na(covariates)){
      condition <- formula(paste0("~",
                                  design,
                                  paste0(covariates[!covariates %in% batch_param],collapse = " + ")))
      
    } else {
      condition <- formula(paste0("~", design))
    }
    mm <- model.matrix(condition, colData(rld))
    if (length(unique(rld[[batch_param]])) > 1) {
      mat <- limma::removeBatchEffect(mat, batch = rld[[batch_param]], design = mm)
      assay(rld) <- mat
    }
  } else {
    if(nrow(dds) < 1000){
      rld <-varianceStabilizingTransformation(dds, blind)
    } else {
      rld <- vst(dds, blind)
    }
  }
  return(rld)
}

get_DESeq_results <- function(dds, exp_metadata, contrasts, design, params, current_group_filter, outdir){
    # Initial setup for DESeq2 contrasts
    resListAll <- list()
    resListFiltered <- list()
    resListDEGs <- list()
    filtered_table <- data.frame()
    Counts  <- counts(dds, normalized = TRUE)
    CPMdds  <- edgeR::cpm(counts(dds, normalized = TRUE))
    mergedDEGs <- c()

    for (x in 1:nrow(contrasts)) {  # For all comparisons to be done  
        condition1 <- contrasts[x, 2] # Control
        condition2 <- contrasts[x, 1] # Experimental

        contrast_string <- paste(condition2, "vs", condition1, sep = "_")

        message(contrast_string)

        exp_metadata_subset <- as.matrix(exp_metadata[exp_metadata[, design] %in% c(condition1, condition2),])
        
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
        SampPerGroup <- table(exp_metadata_subset[, design])
        if (!SampPerGroup[condition2] > 1) { next }
        
        for (gene in 1:nrow(dds)) {
          CountsPass <- NULL
          for (group in 1:length(SampPerGroup)) {
            sampleCols <- grep(names(SampPerGroup)[group], exp_metadata_subset[, design], fixed = T)
            sampleNamesGroup <- exp_metadata_subset[sampleCols, "original_names"]
            Check <- sum(CPMdds[gene,sampleNamesGroup] >= params$MinCount) >= 0.75 * SampPerGroup[group]
            CountsPass <- c(CountsPass, Check)
          }
          if (sum(CountsPass) > 0) {Filter[gene, 1] <- 1 }  else {Filter[gene,1] <- 0 }
          
        }

        compte <- Counts[Filter[,1] == 1,]
        Filter <- Filter[rownames(Filter) %in% rownames(compte), , drop = F]

        intitial_count <- nrow(dds)
        num_relevance_filtered <- nrow(dds) - nrow(Filter)
        message(paste0("Relevance filtering removed ", num_relevance_filtered,
                       " genes from the ", nrow(dds)," assessed. ",
                       nrow(Filter), " genes remaining"))
        
        # run the DEseq p-value calculation
        message("Obtaining the DESeq2 results")
        currentContrast <- c(design, condition2, condition1)
        bpparam <- MulticoreParam(params$cpus)
        res <- DESeq2::results(dds[rownames(compte),],
                               parallel = TRUE,
                               BPPARAM = bpparam,
                               contrast = currentContrast,
                               alpha = params$alpha,
                               pAdjustMethod = 'fdr',
                               cooksCutoff = params$cooks) # If Cooks cutoff disabled - manually inspect.

        res <- lfcShrink(dds,
                         contrast = currentContrast,
                         res = res,
                         BPPARAM = bpparam,
                         type = "ashr")
        resListAll[[contrast_string]] <- res
        
        DEsamples <<- subset(res, res$padj < params$alpha)
        if (nrow(DEsamples) == 0) {
          print("No significant results were found for this contrast. Moving on...")
          next
        }
        DECounts <- compte[rownames(compte) %in% rownames(DEsamples), , drop = F]
        Filter <- Filter[rownames(Filter) %in% rownames(DECounts), , drop = F]
        message("Check median against third quantile" )
        message("AND")
        message("Check the presence of a spike" )

        for (gene in 1:nrow(DECounts)) {
            # Check the median against third quantile
            quantilePass <- NULL
            sampleColsg1 <- grep(dimnames(SampPerGroup)[[1]][1],exp_metadata_subset[,design], fixed = T)
            sampleColsg2 <- grep(dimnames(SampPerGroup)[[1]][2],exp_metadata_subset[,design], fixed = T)
            
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
                sampleColsSpike <- grep(dimnames(SampPerGroup)[[1]][group], exp_metadata_subset[ ,design], fixed = T)
                sampleNamesSpike <- exp_metadata_subset[sampleColsSpike, "original_names"]
                if (max(DECounts[gene,sampleColsSpike]) == 0) {Check <- FALSE} else {
                  Check <- (max(DECounts[gene, sampleColsSpike]) / sum(DECounts[gene, sampleColsSpike])) >=
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

        # Extract the final list of DEGs
        
        message(paste0("Filtering by linear fold-change: linear FC needs to be above ", params$linear_fc_filter))
        
        allCounts_all_filters <- res[rowSums(Filter) == 3 ,]
        DECounts_real <- DEsamples[rowSums(Filter) == 3 & !is.na(DEsamples$padj) &  abs(DEsamples$log2FoldChange) > log2(params$linear_fc_filter) ,]
        DECounts_no_quant <- DEsamples[Filter[, 2] == 0 ,] # save these to output later 
        DECounts_spike <- DEsamples[Filter[, 3] == 0 ,] # save these to output later
        #TODO: output quantile rule failing and spike failing genes

        message(paste0("A total of ", nrow(DECounts_real),
                    " DEGs were selected (out of ", nrow(DECounts) ,"), after ", nrow(DECounts_no_quant),
                    " genes(s) removed by the quantile rule, ", nrow(DECounts_spike),
                    " gene(s) with a spike, and linear fold-change filtering was applied"))
        message("DESeq2 Done")

        mergedDEGs <- c(mergedDEGs, rownames(DECounts_real))
        
        filtered_table <- rbind(filtered_table, data.frame(facet = current_group_filter,
                                                           contrast = contrast_string,
                                                           initial = intitial_count,
                                                           relevance_filtered = num_relevance_filtered,
                                                           quantile_filtered = nrow(DECounts_no_quant),
                                                           spike_filtered = nrow(DECounts_spike),
                                                           passed_all_filters = nrow(DECounts_real)))

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
    return(list(resListAll=resListAll, resListFiltered=resListFiltered, resListDEGs=resListDEGs, mergedDEGs=mergedDEGs, filtered_table=filtered_table))
}


annotate_deseq_table <- function(deseq_results_list, params, filter_results = F, biosets_filter = F) {
  x <- deseq_results_list
  annotated_results <- list()
  for (i in 1:length(x)) {
    deg_table <- x[[i]]
    # Add taxonomy
    if (is.null(deg_table)) {
      next
    } else if (nrow(deg_table) == 0) {
      next
    } else {
      deg_table <- cbind(Feature_ID = row.names(x[[i]]),
                         as(deg_table, "data.frame"),
                         contrast = gsub(pattern = paste0("log2.*", params$design, "\ "),
                                         replacement =  "",
                                         x = x[[i]]@elementMetadata[[2]][2]))
      if(params$platform == "TempO-Seq"){
        deg_table <- dplyr::left_join(deg_table,
                                      params$biospyder, 
                                      by = c(Feature_ID = params$feature_id))
      } else{
        # need to catch a testForValidKeys error in the case where the only resulting genes have ensembl IDs that aren't in the AnnotationDBI database
        result = tryCatch({
          descriptions <- AnnotationDbi::select(get(params$species_data$orgdb),
                                                columns = c("ENSEMBL", "SYMBOL", "GENENAME"),
                                                keys = deg_table$Feature_ID,
                                                keytype="ENSEMBL") %>%
            distinct(ENSEMBL, .keep_all=TRUE)
          colnames(descriptions) <- c("Ensembl_Gene_ID","Gene_Symbol","description")
          descriptions$Feature_ID <- descriptions$Ensembl_Gene_ID
          deg_table <- dplyr::left_join(deg_table, descriptions, by="Feature_ID")
        }, error = function(e) {
          message("error")
        })
      }
      if(!("Gene_Symbol" %in% colnames(deg_table))){
        deg_table$Ensembl_Gene_ID <- deg_table$Feature_ID
        deg_table$Gene_Symbol <- NA
      }
      deg_table <- deg_table %>%
        mutate(linearFoldChange = ifelse(log2FoldChange > 0, 2 ^ log2FoldChange, -1 / (2 ^ log2FoldChange))) %>%
        mutate(Gene_Symbol_2 = coalesce(Gene_Symbol, Ensembl_Gene_ID)) %>%
        dplyr::select(Feature_ID, Ensembl_Gene_ID, Gene_Symbol = Gene_Symbol_2, baseMean, log2FoldChange, linearFoldChange, lfcSE, pvalue, padj, contrast)
      
      ## FILTERS ##
      if (biosets_filter == T) {
        # for biosets, filter on unadjusted p-value
        deg_table <- deg_table[!is.na(deg_table$pval) & deg_table$pval < params$alpha & abs(deg_table$linearFoldChange) > params$linear_fc_filter, ]
      }else if (filter_results == T) {
        deg_table <- deg_table[!is.na(deg_table$padj) & deg_table$padj < params$alpha & abs(deg_table$linearFoldChange) > params$linear_fc_filter, ]
      }
      annotated_results[[i]] <- deg_table %>% dplyr::distinct()
    }
  }
  y <- rbindlist(annotated_results)
  return(y)
}
