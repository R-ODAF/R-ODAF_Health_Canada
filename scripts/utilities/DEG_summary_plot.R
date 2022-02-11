# Find DEG tables

files <- fs::dir_ls("./analysis/DEG_lists/", regexp = "DEG_summary.txt",
                    recurse = T)
files <- as.vector(files)

# Exclude the top-level summary...
files <- files[grep("DEG_lists/DEG_summary", files, invert = T)]

# Function to read in one  at a time and wrangle data
clean_data <- function(x) {
  df <- read.table(x, sep = "\t", header = F, fill = T)  
  df <- as.data.frame(t(df))
  df[1,1] <- "Condition"
  colnames(df) <- df[1,]
  df <- df[2:nrow(df),]
  # subset - shouldn't be needed but it is?
  df <- df[,1:2]
  final <-  reshape2::melt(df, id.vars=1)
  colnames(final) <- c("Condition","Sample","DEGs")
  final$DEGs <- as.numeric(final$DEGs)
  final$Condition <- factor(final$Condition, levels = final$Condition)
  return(final)
}

degs <- files %>%
  map_dfr(clean_data)

###################################################
# Optional... Data cleaning examples
# Some fun regex if you want to split things up...
# This takes the last space in the string as the delimiter
degs <- degs %>% tidyr::separate(Sample, c("Chemical","Timepoint"),
                                 sep = "\\s+(?=\\S*$)" )
#
###################################################

###################################################
# Plot
ggplot(degs, aes(x = Condition, y = DEGs, fill = Timepoint)) +
  geom_bar(stat="identity", position="dodge") +
  facet_wrap(~Chemical, scales = "free_x") +
  theme(axis.text.x = element_blank()) # element_text(angle = 90))

