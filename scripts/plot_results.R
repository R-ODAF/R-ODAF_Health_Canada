# Plot results from empirical FDR calculations

library(tidyverse)
theme_set(theme_bw())
raw_data <- list()
processed_data <- list()

processed_data$tPOds <- read.table(file.path(here::here(), "data", "processed",
                                             "PoD Meier FDR data - from Andrew, wrangled.txt"),
                                   sep = "\t", header = T, check.names = F) %>%
  pivot_longer(cols = c(5,7:8), names_to = "tPOD method", values_to = "FDR")

processed_data$tPOds <- processed_data$tPOds %>%
  mutate(Design = fct_reorder(Design, `Samples per group` , .fun = sum))

processed_data$tPOds <- processed_data$tPOds %>%
  mutate(display_category = paste0(`tPOD method`,", ",Doses," doses"))

processed_data$tPOds$display_category <- factor(processed_data$tPOds$display_category)

processed_data$tPOds$Filtering <- 
  factor(processed_data$tPOds$Filtering,
         levels = c("Default settings",
                    "Background Filtering",
                    "Williams Trend Test",
                    "Fold Change"
                    ))

processed_data$tPOds <- processed_data$tPOds  %>%
  mutate(nice_samples_per_group = paste0("n = ",`Samples per group`)) %>%
  mutate(nice_samples_per_group = fct_reorder(nice_samples_per_group, `Samples per group`))
  
  
  
  
  #########################

ggplot(processed_data$tPOds, aes(x = `tPOD method`, y = log10(FDR), fill = `tPOD method`)) +
  geom_jitter() +
  facet_grid(~`tPOD method`, scales = "free_x")

ggplot(processed_data$tPOds, aes(x = `tPOD method`, y = log10(FDR), fill = `tPOD method`)) +
  geom_jitter() +
  facet_grid(~`tPOD method`, scales = "free_x")

ggplot(processed_data$tPOds, aes(x = `Samples per group`, y = log10(FDR), fill = `tPOD method`)) +
  geom_jitter() +
  facet_grid(~`tPOD method`, scales = "free_x")

ggplot(processed_data$tPOds, aes(x = `Samples per group`, y = FDR, fill = `tPOD method`)) +
  geom_jitter() +
  facet_grid(Filtering~`tPOD method`, scales = "free")

ggplot(processed_data$tPOds, aes(x = `tPOD method`, y = FDR, shape = `tPOD method`)) +
  geom_jitter() +
  facet_grid(Filtering~`Samples per group`, scales = "free_x")

ggplot(processed_data$tPOds, aes(x = `tPOD method`, y = log2(FDR), shape = `tPOD method`)) +
  geom_jitter() +
  facet_grid(Filtering~`Samples per group`, scales = "free_x")

ggplot(processed_data$tPOds, aes(x = `Samples per group`, y = log2(FDR), shape = `tPOD method`)) +
  geom_jitter() +
  facet_grid(~Filtering, scales = "free_x")

ggplot(processed_data$tPOds, aes(x = `Samples per group`, y = log10(FDR), shape = `tPOD method`)) +
  geom_point() +
  facet_grid(~Filtering, scales = "free_x")

# Good plots
ggplot(processed_data$tPOds %>% dplyr::filter(!`tPOD method` %in% "Median BMCs"),
       aes(x = Design, y = FDR, fill = `tPOD method`)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept=0.05, col = "red") +
  facet_grid(~Filtering, scales = "free_x") +
  xlab("Group size (n)") +
  scale_x_discrete(labels = function(x) gsub(",", ",\n", x))

ggplot(processed_data$tPOds %>% dplyr::filter(!`tPOD method` %in% "Median BMCs"),
       aes(x = Design, y = FDR, fill = `tPOD method`)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept=0.05, col = "red", lty = "dashed") +
  facet_grid(~Filtering, scales = "free_x", switch = "x") +
  xlab(element_blank()) +
  scale_x_discrete(labels = function(x) gsub(",", ",\n", x)) +
  theme_minimal() +
  theme(legend.position = "bottom",
        #axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text.x = element_text(size = rel(1.5))) +
  ylab("FDR") +
  scale_fill_discrete(name = "tPOD Method")


ggplot(processed_data$tPOds %>% dplyr::filter(!`tPOD method` %in% "Median BMCs"),
       aes(x = Design, y = FDR, fill = Filtering)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept=0.05, col = "red") +
  facet_grid(~`tPOD method`, scales = "free_x") +
  xlab("Group size (n)") +
  scale_x_discrete(labels = function(x) gsub(",", ",\n", x))


custom.col <- c("#00AFBB", "#E7B800", "#FC4E07","#FFDB6D", "#C4961A", "#F4EDCA", 
                "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")

cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(processed_data$tPOds %>% dplyr::filter(!`tPOD method` %in% "Median BMCs"),
       aes(x = nice_samples_per_group, y = FDR, fill = Filtering)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept=0.05, col = "red", lty = "dashed") +
  facet_grid(`Cell line`~display_category, scales = "free_x", switch = "x") +
  xlab(element_blank()) +
  scale_x_discrete(labels = function(x) gsub(",", ",\n", x)) +
  theme_minimal() +
  theme(legend.position = "right",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = rel(1)),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text.x = element_text(size = rel(1))) +
  #ggtitle("Title") +
  ylab("FDR") +
  scale_fill_discrete(name = "Filtering",type = cbp1)
