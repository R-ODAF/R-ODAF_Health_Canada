
source(here::here("scripts","data_functions.R"), local = TRUE)

# input data file
dataFile <- file.path(paths$RData, paste0(params$project_title, "_DEG_data.RData"))
#load(dataFile) # metadata, contrasts, counts, resultsList
attach(dataFile)
ddsList<-ddsList
overallResListAll<-overallResListAll
overallResListDEGs<-overallResListDEGs
rldList<-rldList
mergedDEGsList<-mergedDEGsList
designList<-designList
contrastsList<-contrastsList
design_to_use<-design_to_use
contrasts<-contrasts
filtered_table<-filtered_table
exp_metadata<-exp_metadata # original unfaceted design
exp_metadata_original <- exp_metadata
detach()

# reformat data based on group_facet and display_group_facet

# case 1: no facet, no display facet
if(is.na(params$group_facet) && is.na(params$display_group_facet)){
  dds <- ddsList[['all']]
  resultsListAll <- overallResListAll[['all']]
  resultsListDEGs <- overallResListDEGs[['all']]
  rld <- rldList[['all']]
  mergedDEGs <- mergedDEGsList[['all']]
  
  exp_metadata_subset<-exp_metadata
  contrasts_subset<-contrasts
  
  # case 2: no facet, yes display facet
} else if(is.na(params$group_facet) && !is.na(params$display_group_facet)){
  # the data isn't already faceted but we need to extract the facet we need
  display_group_filter <- params$display_group_filter
  
  dds_all <- ddsList[['all']]
  resultsListAll_all <- overallResListAll[['all']]
  resultsListDEGs_all <- overallResListDEGs[['all']]
  rld_all <- rldList[['all']]

  metadata_subset <- subset_metadata(designList[['all']], design_to_use, contrasts, params$display_group_facet, display_group_filter)
  exp_metadata_subset <- metadata_subset$exp_metadata
  contrasts_subset <- metadata_subset$contrasts
  
  dds_subset <- subset_data(dds_all, exp_metadata_subset)
  rld_subset <- subset_data(rld_all, exp_metadata_subset)
  contrast_strings <- contrasts_subset %>% mutate(contrast_string = paste(V1,V2,sep="_vs_")) %>% pull(contrast_string)
  resultsListAll_subset <- resultsListAll_all[contrast_strings]
  resultsListDEGs_subset <- resultsListDEGs_all[contrast_strings]

  exp_metadata <- exp_metadata_subset
  contrasts <- contrasts_subset
  dds <- dds_subset
  resultsListAll <- resultsListAll_subset
  resultsListDEGs <- resultsListDEGs_subset
  rld <- rld_subset
  # note, in this case the calculated merged DEGs will be for the whole experiment, not the display facet
  # so let's recalculate them
  mergedDEGs <- unique(unlist(sapply(resultsListDEGs,rownames),use.names=FALSE))
  # TODO: add some tests here to make sure everything worked properly
  
  # case 3: yes facet, yes display facet
} else if(!is.na(params$group_facet) && !is.na(params$display_group_facet)){
  if(params$group_facet != params$display_group_facet) {
    stop("Error: display_group_facet must match group_facet, otherwise DESeq2 results get mixed and matched.")
  }
  display_group_filter <- params$display_group_filter
  dds <- ddsList[[display_group_filter]]
  resultsListAll <- overallResListAll[[display_group_filter]]
  resultsListDEGs <- overallResListDEGs[[display_group_filter]]
  rld <- rldList[[display_group_filter]]
  mergedDEGs <- mergedDEGsList[[display_group_filter]]
  exp_metadata_subset <- designList[[display_group_filter]]
  contrasts_subset <- contrastsList[[display_group_filter]]
  # case 4: yes facet, no display facet, this one doesn't make sense
} else {
  stop("Making a single report for faceted data not supported. Did you forget to set display_group_facet?")
}


# filter the regularized data a couple ways for different displays    
rld_DEGs <- rld[row.names(assay(rld)) %in% mergedDEGs]

rv = rowVars(assay(rld))
select = order(rv, decreasing = TRUE)[1:params$nBest]
rld_top <- rld[select,]
select_heatmap = order(rv, decreasing = TRUE)[1:params$nHeatmapDEGs]
rld_top_heatmap <- rld[select_heatmap,]

allResults <- annotate_deseq_table(resultsListAll, params, filter_results = F)
significantResults <- annotate_deseq_table(resultsListDEGs, params, filter_results = F)

ordered_contrast_strings <- contrasts %>% mutate(contrast_string = paste(V1,'vs',V2,sep=" ")) %>% pull(contrast_string)


allResults$contrast <- factor(allResults$contrast, levels = ordered_contrast_strings)
significantResults$contrast <- factor(significantResults$contrast, levels = ordered_contrast_strings)
