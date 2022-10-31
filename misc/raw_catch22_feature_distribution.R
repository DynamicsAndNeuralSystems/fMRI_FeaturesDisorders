library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

# UCLA Schizophrenia
dataset_ID <- "UCLA_Schizophrenia"
data_path <- "/headnode1/abry4213/data/UCLA_Schizophrenia/"
rdata_path <- paste0(data_path, "processed_data/Rdata/")
plot_path <- "/headnode1/abry4213/data/UCLA_Schizophrenia/plots/"
noise_proc <- "AROMA+2P+GMR"

# Load raw catch22 data for UCLA schizophrenia
raw_catch22_data <- readRDS(paste0(rdata_path, "UCLA_Schizophrenia_catch22_filtered.Rds"))

# Plot feature distribution by each feature
raw_catch22_data %>%
  filter(Noise_Proc == noise_proc) %>%
  mutate(names = str_replace_all(names, "_", " ")) %>%
  ggplot(data=., mapping=aes(x = values, fill = names)) +
  geom_histogram() +
  facet_wrap(names ~ ., scales="free", ncol=4, 
             labeller = labeller(names = label_wrap_gen(22))) +
  ggtitle("catch22 Raw Feature Values for UCLA Schizophrenia") +
  ylab("# Regions/Subjects") +
  xlab("Feature Value") +
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5))
ggsave(paste0(plot_path, "catch22_raw_dist_UCLA_Schizophrenia_AROMA_2P_GMR.png"),
       width = 11, height = 11, units="in", dpi=300, background="white")