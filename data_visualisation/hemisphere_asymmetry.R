################################################################################
# Load libraries
################################################################################

# python_to_use <- "~/.conda/envs/pyspi/bin/python3"
python_to_use <- "/Users/abry4213/anaconda3/envs/pyspi/bin/python3"
pairwise_feature_set <- "pyspi14"
univariate_feature_set <- "catch24"
study_group_df <- data.frame(Study = c(rep("UCLA_CNP", 3)),
                             Noise_Proc = c(rep("AROMA+2P+GMR",3)),
                             Comparison_Group = c("Schizophrenia", "Bipolar", "ADHD"),
                             Group_Nickname = c("SCZ", "BPD", "ADHD"))
reticulate::use_python(python_to_use)

library(reticulate)
library(tidyverse)
library(icesTAF)
library(cowplot)
library(theft)
library(glue)
library(ggseg)
theme_set(theme_cowplot())
pyarrow_feather <- import("pyarrow.feather")

################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
source(glue("{github_dir}/data_visualisation/Manuscript_Draft_Visualisations_Helper.R"))
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/FigureS2/")
icesTAF::mkdir(plot_path)

# Load in univariate time-series feature info
TS_feature_info <- read.csv(glue("{github_dir}/data_visualisation/catch24_info.csv"))
UCLA_CNP_brain_region_info <- read.csv("~/data/UCLA_CNP/study_metadata/UCLA_CNP_Brain_Region_info.csv")

# Load catch24 data for UCLA CNP
UCLA_CNP_catch24 <- pyarrow_feather$read_feather("~/data/UCLA_CNP/processed_data/UCLA_CNP_AROMA_2P_GMR_catch24_filtered.feather")

# Load study metadata
UCLA_CNP_metadata <- pyarrow_feather$read_feather("~/data/UCLA_CNP/study_metadata/UCLA_CNP_sample_metadata.feather") 


# Calculate asymmetry index as in equation 4 in https://www-ncbi-nlm-nih-gov.ezproxy.library.sydney.edu.au/pmc/articles/PMC2726301/
asymmetry_index <- UCLA_CNP_catch24 %>%
  left_join(., UCLA_CNP_brain_region_info) %>%
  mutate(Hemisphere = case_when(str_detect(Brain_Region, "Left|lh-") ~ "Left",
                                str_detect(Brain_Region, "Right|rh-") ~ "Right")) %>%
  mutate(region = gsub("Left-|lh_|Right-|rh_", "", label)) %>%
  pivot_wider(id_cols = c(Sample_ID, region, names), 
              names_from = Hemisphere, values_from = values) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(AI = (Left - Right) / (abs(Left) + abs(Right)))


# Which features are generally most hemisphere-specific?
AI_by_feature <- asymmetry_index %>%
  left_join(UCLA_CNP_metadata) %>%
  group_by(names, Diagnosis) %>%
  summarise(mean_AI = mean(abs(AI), na.rm=T)) %>%
  arrange(desc(mean_AI))

# White brain regions generally exhibit the most hemispheric asymmetry?
AI_by_region <- asymmetry_index %>%
  left_join(UCLA_CNP_metadata) %>%
  group_by(region, Diagnosis) %>%
  summarise(mean_AI = mean(abs(AI), na.rm=T)) %>%
  arrange(desc(mean_AI))

# What is the global absolute AI for each condition?
AI_by_region %>%
  group_by(Diagnosis) %>%
  summarise(mean_global_AI = mean(mean_AI))

# Let's plot these per group
AI_by_region %>%
  ungroup() %>%
  mutate(label = paste0("lh_", region)) %>%
  dplyr::select(-region) %>%
  filter(Diagnosis == "Schizophrenia") %>%
  left_join(., dk %>% as_tibble) %>%
  ggseg(atlas = "dk", mapping = aes(fill = mean_AI),
        position = "stacked", colour = "gray70", hemisphere = "left") +
  theme_void() +
  theme(plot.title = element_blank()) 