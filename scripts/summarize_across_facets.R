library('tidyverse')
library('devtools')
library('DESeq2')
library('data.table')
library('yaml')
library('stringr')


# load params
source(here::here("scripts","file_functions.R"))
source(here::here("scripts","setup_functions.R"))
source(here::here("scripts","DESeq_functions.R"))


config <- yaml::read_yaml(file.path(here::here(),
                                    "config/config.yaml"),
                          eval.expr = T)

# Combine required params from config
params <- c(config$common, config$DESeq2)
# replace nulls in params with NA
params <- replace_nulls_in_config(params)
# If projectdir is not set, figure out current project root directory
projectdir <- params$projectdir
if (is.na(projectdir)) {
  projectdir <- here::here()
  params$projectdir <- projectdir
}

paths <- set_up_paths(params)
species_data <- load_species(params$species, params$wikipathways_filename, params$biospyder_manifest_file)
params$species_data <- species_data
params <- set_up_platform_params(params)

skip_extra <- c("DMSO","water") # Remove DMSO controls as a facet

# input data file
dataFile <- file.path(paths$RData, paste0(params$project_title, "_DEG_data.RData"))
load(dataFile) # metadata, contrasts, counts, resultsList

allResultsUnfaceted <- data.frame()
significantResultsUnfaceted <- data.frame()

ordered_contrasts <- paste0(contrasts$V1, " vs ", contrasts$V2)
ordered_levels <- expand.grid(facets,ordered_contrasts) %>% arrange(Var1) %>% mutate(asdf = paste0(Var1, ": ", Var2)) %>% pull(asdf)

for (f in facets){
  resultsListAll <- overallResListAll[[f]]
  resultsListDEGs <- overallResListDEGs[[f]]
  design <- designList[[f]]
  if (length(resultsListAll) > 0) {
    allResults <- annotate_deseq_table(resultsListAll, params, filter_results = F)
    allResultsUnfaceted <- allResults %>% mutate(facet=f) %>% rbind(allResultsUnfaceted)
  }
  if (length(resultsListDEGs) > 0) {
    significantResults <- annotate_deseq_table(resultsListDEGs, params, filter_results = F)
    significantResultsUnfaceted <- significantResults %>% mutate(facet=f) %>% rbind(significantResultsUnfaceted)
  }
}

prefix <- paste0(params$platform, "_",
                 params$project_title, "_",
                 paste(params$current_filter, collapse = "_"), "_",
                 format(Sys.time(),'%d-%m-%Y.%H.%M'))  


p1_data <- significantResultsUnfaceted %>%
    mutate(facet_contrast = factor(paste0(facet, ": ", contrast), levels = ordered_levels))


if (length(facets) < 10) {
  plot_size = 5 
} else {
  plot_size = round(sqrt(length(facets)))*1.5
}

if(length(facets) == 1){
  p1 = ggplot(p1_data, aes(x=contrast)) +
    geom_bar(aes(y=..count..)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, hjust=1)) +
    ylab("Number of DEGs") +
    xlab("Contrast")
} else if (length(facets) < 10) {
    p1 = ggplot(p1_data, aes(x=facet_contrast)) +
        geom_bar(aes(y=..count.., fill=facet)) +
        theme_bw() +
        theme(axis.text.x = element_text(angle=90, hjust=1)) +
        ylab("Number of DEGs") +
        xlab("Facet: contrast")
} else {
    p1 = ggplot(p1_data, aes(x=facet_contrast)) +
        geom_bar(aes(y=..count.., fill=facet)) +
        theme_bw() +
        theme(axis.text.x = element_blank(),
            legend.position = "none") +
        facet_wrap(~facet, scales = "free_x") +
        ylab("Number of DEGs") +
        xlab("Facet: contrast")
}
ggsave(file.path(report_dir,paste0(prefix,"_","DEG_summary_plot.png")),p1,
       width=plot_size, height=plot_size, units="in", dpi=300)


# plot filter stats
p2_data <- filtered_table %>%
  group_by(facet,contrast) %>%
  mutate(not_significant = initial-(relevance_filtered+quantile_filtered+spike_filtered)) %>%
  mutate(contrast = str_replace(contrast, "_vs_", " vs ")) %>%
  ungroup() %>%
  pivot_longer(cols=c(not_significant,relevance_filtered,quantile_filtered,spike_filtered,passed_all_filters)) %>%
  mutate(perc = value/initial) %>%
  mutate(name = factor(name, levels=c("relevance_filtered", "not_significant", "quantile_filtered", "spike_filtered", "passed_all_filters"))) %>%
  mutate(facet_contrast = factor(paste0(facet, ": ", contrast), levels = ordered_levels))

if(length(facets) == 1){
  p2 = ggplot(p2_data, aes(x=contrast,y=value)) +
    theme_bw() +
    geom_bar(stat="identity",position="dodge") +
    facet_wrap(~name, scales="free") +
    theme(axis.text.x = element_text(angle=90, hjust=1)) +
    ylab("Percent of all reads") +
    xlab("Contrast")
  
} else if (length(facets) < 10) {
  p2 = ggplot(p2_data, aes(x=facet_contrast,y=value,fill=facet)) +
    theme_bw() +
    geom_bar(stat="identity",position="dodge") +
    facet_wrap(~name, scales="free") +
    theme(axis.text.x = element_text(angle=90, hjust=1)) +
    ylab("Percent of all reads") +
    xlab("Facet: contrast")
} else {
  p2 = ggplot(p2_data, aes(x=facet_contrast,y=value,fill=facet)) +
    theme_bw() +
    geom_bar(stat="identity",position="dodge") +
    facet_wrap(~name, scales="free") +
    theme(axis.text.x = element_blank(),
          legend.position = "none") +
    ylab("Percent of all reads") +
    xlab("Facet: contrast")
}
  
ggsave(file.path(report_dir,paste0(prefix,"_","filter_summary_plot.png")),p2,
       width=plot_size, height=plot_size, units="in", dpi=300)


