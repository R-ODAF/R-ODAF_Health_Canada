load(params$dataFile) # metadata, contrasts, counts, resultsList

# reformat data based on group_facet and display_group_facet

# case 1: no facet, no display facet
if(is.na(params$group_facet) && is.na(params$display_group_facet)){
  dds <- ddsList[['all']]
  resultsListAll <- overallResListAll[['all']]
  resultsListDEGs <- overallResListDEGs[['all']]
  rld <- rldList[['all']]
  mergedDEGs <- mergedDEGsList[['all']]
  
  # case 2: no facet, yes display facet
} else if(is.na(params$group_facet) && !is.na(params$display_group_facet)){
  # the data isn't already faceted but we need to extract the facet we need
  display_group_filter <- params$display_group_filter
  
  dds_all <- ddsList[['all']]
  resultsListAll_all <- overallResListAll[['all']]
  resultsListDEGs_all <- overallResListDEGs[['all']]
  rld_all <- rldList[['all']]
  mergedDEGs_all <- mergedDEGsList[['all']]
  
  metadata_subset <- subset_metadata(designList[['all']], design_to_use, contrasts, params$display_group_facet, display_group_filter)
  DESeqDesign_subset <- metadata_subset$DESeqDesign
  contrasts_subset <- metadata_subset$contrasts
  
  dds_subset <- subset_data(dds_all, DESeqDesign_subset)
  rld_subset <- subset_data(rld_all, DESeqDesign_subset)
  
  contrast_strings <- contrasts_subset %>% mutate(contrast_string = paste(V1,V2,sep="_vs_")) %>% pull(contrast_string)
  resultsListAll_subset <- resultsListAll_all[contrast_strings]
  resultsListDEGs_subset <- resultsListDEGs_all[contrast_strings]
  
  # note, in this case the merged DEGs will be for the whole experiment, not the display facet
  DESeqDesign <- DESeqDesign_subset
  contrasts <- contrasts_subset
  dds <- dds_subset
  resultsListAll <- resultsListAll_subset
  resultsListDEGs <- resultsListDEGs_subset
  rld <- rld_subset
  mergedDEGs <- mergedDEGs_all
  
  # TODO: add some tests here to make sure everything worked properly
  
  # case 3: yes facet, yes display facet
} else if(!is.na(params$group_facet) && !is.na(params$display_group_facet)){
  display_group_filter <- params$display_group_filter
  
  dds <- ddsList[[display_group_filter]]
  resultsListAll <- overallResListAll[[display_group_filter]]
  resultsListDEGs <- overallResListDEGs[[display_group_filter]]
  rld <- rldList[[display_group_filter]]
  mergedDEGs <- mergedDEGsList[[display_group_filter]]
  
  # case 4: yes facet, no display facet, this one doesn't make sense
} else {
  stop("Making a single report for faceted data not supported. Did you forget to set display_group_facet?")
}


# filter the regularized data a couple ways for different displays    
rld_DEGs <- rld[row.names(assay(rld)) %in% mergedDEGs]

rv = rowVars(assay(rld))
select = order(rv, decreasing = TRUE)[1:params$nBest]
rld_top <- rld[select,]

allResults <- annotate_deseq_table(resultsListAll, params, bs, filter_results = F)

