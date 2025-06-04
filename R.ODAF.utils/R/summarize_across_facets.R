#' Process DESeq2 Results
#'
#' Annotates DESeq2 results and constructs a unified dataframe for all results.
#'
#' @param resultsListAll List of all results.
#' @param resultsListDEGs List of DEG results.
#' @param facets Vector of facets to process.
#' @param params List of parameters for annotation.
#' @return A list with combined results for all facets and DEGs.
#'
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import stringr
#' @export
summarize_across_facets <- function(overallResListAll, overallResListDEGs, filtered_table, facets, params) {
   allResultsUnfaceted <- data.frame()
   significantResultsUnfaceted <- data.frame()

   ordered_contrasts <- paste0(exp_contrasts$V1, " vs ", exp_contrasts$V2)
   ordered_levels <- expand.grid(facets, ordered_contrasts) %>%
      arrange(Var1) %>%
      mutate(ordered_levels = paste0(Var1, ": ", Var2)) %>%
      pull(ordered_levels)

   for (f in facets) {
      if (length(overallResListAll[[f]]) > 0) {
         allResults <- annotate_deseq_table(overallResListAll[[f]], params, filter_results = FALSE)
         allResultsUnfaceted <- rbind(allResultsUnfaceted, mutate(allResults, facet = f))
      }
      if (length(overallResListDEGs[[f]]) > 0) {
         significantResults <- annotate_deseq_table(overallResListDEGs[[f]], params, filter_results = FALSE)
         significantResultsUnfaceted <- rbind(significantResultsUnfaceted, mutate(significantResults, facet = f))
      }
   }

   prefix <- paste0(params$platform, "_",
      params$project_title, "_",
      paste(params$current_filter, collapse = "_"), "_",
      format(Sys.time(), '%d-%m-%Y.%H.%M'))

   mixedrank <- function(x) order(gtools::mixedorder(x))


   p1_data <- significantResultsUnfaceted %>%
    mutate(
    facet_contrast = paste0(facet, ": ", contrast) # Plain string
    )

   p1_data$facet_contrast <- factor(
    p1_data$facet_contrast,
    levels = mixedsort(unique(p1_data$facet_contrast))
   )

   if (length(facets) < 10) {
      plot_size <- 5
   } else {
      plot_size <- round(sqrt(length(facets))) * 1.5
   }

   if (length(facets) == 1) {
      p1 = ggplot(p1_data, aes(x = facet_contrast)) +
         geom_bar(aes(y = ..count..)) +
         theme_bw() +
         theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
         ylab("Number of DEGs") +
         xlab("Contrast")
   } else if (length(facets) < 10) {
      p1 = ggplot(p1_data, aes(x = facet_contrast)) +
         geom_bar(aes(y = ..count.., fill = facet)) +
         theme_bw() +
         theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
         ylab("Number of DEGs") +
         xlab("Facet: contrast")
   } else {
      p1 = ggplot(p1_data, aes(x = facet_contrast)) +
         geom_bar(aes(y = ..count.., fill = facet)) +
         theme_bw() +
         theme(axis.text.x = element_blank(),
            legend.position = "none") +
         facet_wrap(~facet, scales = "free_x") +
         ylab("Number of DEGs") +
         xlab("Facet: contrast")
   }
   # ggsave(file.path(paths$reports_dir, paste0(prefix, "_", "DEG_summary_plot.png")), p1,
   #    width = plot_size, height = plot_size, units = "in", dpi = 300)


   # plot filter stats


   p2_data <- filtered_table %>%
   group_by(facet, contrast) %>%
   mutate(not_significant = initial - (relevance_filtered + quantile_filtered + spike_filtered)) %>%
   mutate(contrast = str_replace(contrast, "_vs_", " vs ")) %>%
   ungroup() %>%
   tidyr::pivot_longer(
      cols = c(not_significant, relevance_filtered, quantile_filtered, spike_filtered, passed_all_filters)
   ) %>%
   mutate(perc = value / initial) %>%
   mutate(name = factor(
      name,
      levels = c("relevance_filtered", "not_significant", "quantile_filtered", "spike_filtered", "passed_all_filters")
   )) %>%
   mutate(facet_contrast = paste0(facet, ": ", contrast)) # facet_contrast remains *character* for now

   # Set the right levels using gtools::mixedsort()!
   ordered_levels <- mixedsort(unique(p2_data$facet_contrast))
   p2_data$facet_contrast <- factor(p2_data$facet_contrast, levels = ordered_levels)

   ordered_contrasts <- mixedsort(unique(p2_data$contrast))
   p2_data$contrast <- factor(p2_data$contrast, levels = ordered_contrasts)

   ## Now plotting code as before, referencing the factor columns
   if (length(facets) == 1) {
   p2 = ggplot(p2_data, aes(x = contrast, y = value)) +
      theme_bw() +
      geom_bar(stat = "identity", position = "dodge") +
      facet_wrap(~name, scales = "free") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      ylab("Number of genes") +
      xlab("Contrast")
   } else if (length(facets) < 10) {
   p2 = ggplot(p2_data, aes(x = facet_contrast, y = value, fill = facet)) +
      theme_bw() +
      geom_bar(stat = "identity", position = "dodge") +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      ylab("Number of genes") +
      xlab("Facet: contrast")
   } else {
   p2 = ggplot(p2_data, aes(x = facet_contrast, y = value, fill = facet)) +
      theme_bw() +
      geom_bar(stat = "identity", position = "dodge") +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      theme(axis.text.x = element_blank(),
            legend.position = "none") +
      ylab("Number of genes") +
      xlab("Facet: contrast")
   }

   # ggsave(file.path(paths$reports_dir, paste0(prefix, "_", "filter_summary_plot.png")), p2,
   #    width = plot_size, height = plot_size, units = "in", dpi = 300)
  return(list(p1 = p1, p2 = p2))
}
