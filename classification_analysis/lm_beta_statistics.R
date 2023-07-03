################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/t_statistics/")
TAF::mkdir(plot_path)

python_to_use <- "~/.conda/envs/pyspi/bin/python3"
# python_to_use <- "/Users/abry4213/anaconda3/envs/pyspi/bin/python3"
pairwise_feature_set <- "pyspi14"
univariate_feature_set <- "catch24"
data_path <- "~/data/TS_feature_manuscript"
study_group_df <- data.frame(Study = c(rep("UCLA_CNP", 3), "ABIDE_ASD"),
                             Noise_Proc = c(rep("AROMA+2P+GMR",3), "FC1000"),
                             Comparison_Group = c("Schizophrenia", "Bipolar", "ADHD", "ASD"),
                             Group_Nickname = c("SCZ", "BPD", "ADHD", "ASD"))

reticulate::use_python(python_to_use)

library(reticulate)

# Import pyarrow.feather as pyarrow_feather
pyarrow_feather <- import("pyarrow.feather")

################################################################################
# Load libraries
################################################################################
library(feather)
library(tidyverse)
library(glue)
library(icesTAF)
library(broom)

# Source visualisation script
source(glue("{github_dir}/data_visualisation/Manuscript_Draft_Visualisations_Helper.R"))

# Load in univariate time-series feature info
catch24_info <- read.csv(glue("{github_dir}/data_visualisation/catch24_info.csv"))
pyspi14_info <- read.csv(glue("{github_dir}/data_visualisation/SPI_info.csv"))

# Load study metadata
UCLA_CNP_metadata <- pyarrow_feather$read_feather("~/data/UCLA_CNP/study_metadata/UCLA_CNP_sample_metadata.feather") 
ABIDE_ASD_metadata <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_sample_metadata.feather") 

# Load brain region info
UCLA_CNP_brain_region_info <- read.csv("~/data/UCLA_CNP/study_metadata/UCLA_CNP_Brain_Region_info.csv")
ABIDE_ASD_brain_region_info <- read.table("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_Harvard_Oxford_cort_prob_2mm_ROI_lookup.txt", sep="\t", header=TRUE)

# Load SPI info
pyspi14_info <- read.csv(glue("{github_dir}/data_visualisation/SPI_info.csv"))

################################################################################
# Z-score helper function
z_score <- function(data) {
  data_z = (data - mean(data))/sd(data)
  return(data_z)
}

# Ridge plot for catch24 features' Mann-Whitney U-statistics across entire brain
lm_beta_stats_for_group <- function(comparison_group, input_data, study, group_nickname){
  res <- input_data %>%
    filter(Diagnosis %in% c(comparison_group, "Control"),
           Study == study) %>%
    mutate(Diagnosis = case_when(Diagnosis == "Schizophrenia" ~ "SCZ",
                                 Diagnosis == "Bipolar" ~ "BPD",
                                 T ~ Diagnosis)) %>%
    dplyr::select(Brain_Region, names, Diagnosis, values) %>%
    mutate(Diagnosis = factor(Diagnosis, levels = c("Control", group_nickname))) %>%
    group_by(Brain_Region, names) %>%
    mutate(values_z = z_score(values)) %>%
    nest() %>%
    mutate(
      fit = map(data, ~ lm(values_z ~ Diagnosis, data = .x)),
      tidied = map(fit, tidy)
    ) %>% 
    unnest(tidied) %>%
    dplyr::select(-data, -fit) %>%
    ungroup() %>%
    dplyr::rename("TS_Feature" = "names") %>% 
    filter(term != "(Intercept)") %>%
    mutate(Comparison_Group = group_nickname,
           Study = study)
  
  return(res)
}

if (!file.exists(glue("{data_path}/univariate_catch24_lm_beta_statistics_by_brain_region.feather"))) {
  # Load raw feature data
  UCLA_CNP_catch24 <- pyarrow_feather$read_feather("~/data/UCLA_CNP/processed_data/UCLA_CNP_AROMA_2P_GMR_catch24_filtered.feather")  %>%
    left_join(., UCLA_CNP_metadata) %>%
    mutate(Study = "UCLA_CNP")
  ABIDE_ASD_catch24 <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/processed_data/ABIDE_ASD_FC1000_catch24_filtered.feather")  %>%
    left_join(., ABIDE_ASD_metadata) %>%
    mutate(Study = "ABIDE_ASD")
  
  combined_univariate_data <- plyr::rbind.fill(UCLA_CNP_catch24, ABIDE_ASD_catch24)
  
  beta_stats_catch24_whole_brain <- 1:4 %>%
    purrr::map_df(~ lm_beta_stats_for_group(input_data = combined_univariate_data,
                                      comparison_group = study_group_df$Comparison_Group[.x],
                                      study = study_group_df$Study[.x],
                                      group_nickname = study_group_df$Group_Nickname[.x]))
  feather::write_feather(beta_stats_catch24_whole_brain, glue("{data_path}/univariate_catch24_lm_beta_statistics_by_brain_region.feather"))
} else {
  beta_stats_catch24_whole_brain <- feather::read_feather(glue("{data_path}/univariate_catch24_lm_beta_statistics_by_brain_region.feather"))
}

# Pairwise pyspi14 lm fitting
lm_beta_stats_for_group_pairwise <- function(comparison_group, input_data, study, group_nickname){
  res <- input_data %>%
    filter(Diagnosis %in% c(comparison_group, "Control"),
           Study == study) %>%
    mutate(Diagnosis = case_when(Diagnosis == "Schizophrenia" ~ "SCZ",
                                 Diagnosis == "Bipolar" ~ "BPD",
                                 T ~ Diagnosis)) %>%
    rowwise() %>%
    mutate(Region_Pair = paste0(brain_region_from, "_", brain_region_to)) %>%
    dplyr::select(Region_Pair, SPI, Diagnosis, value) %>%
    mutate(Diagnosis = factor(Diagnosis, levels = c("Control", group_nickname))) %>%
    group_by(Region_Pair, SPI) %>%
    mutate(values_z = z_score(value)) %>%
    nest() %>%
    mutate(
      fit = map(data, ~ possibly(lm, otherwise = NA_real_)(values_z ~ Diagnosis, data = .x)),
      tidied = map(fit, possibly(tidy, otherwise = NA_real_))
    ) %>% 
    unnest(tidied) %>%
    dplyr::select(-data, -fit) %>%
    ungroup() %>%
    filter(term != "(Intercept)") %>%
    mutate(Comparison_Group = group_nickname,
           Study = study)
  
  return(res)
}

if (!file.exists(glue("{data_path}/pairwise_pyspi14_lm_beta_statistics_by_region_pair.feather"))) {
  UCLA_CNP_pyspi14 <- pyarrow_feather$read_feather("~/data/UCLA_CNP/processed_data/UCLA_CNP_AROMA_2P_GMR_pyspi14_filtered.feather")  %>%
    left_join(., UCLA_CNP_metadata) %>%
    filter(!is.na(Diagnosis)) %>%
    mutate(Study = "UCLA_CNP")
  ABIDE_ASD_pyspi14 <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/processed_data/ABIDE_ASD_FC1000_pyspi14_filtered.feather")  %>%
    left_join(., ABIDE_ASD_metadata) %>%
    filter(!is.na(Diagnosis)) %>%
    mutate(Study = "ABIDE_ASD")
  
  beta_stats_pyspi14_whole_brain <- 1:4 %>%
    purrr::map_df(~ lm_beta_stats_for_group_pairwise(input_data = plyr::rbind.fill(UCLA_CNP_pyspi14,
                                                                    ABIDE_ASD_pyspi14),
                                      comparison_group = study_group_df$Comparison_Group[.x],
                                      study = study_group_df$Study[.x],
                                      group_nickname = study_group_df$Group_Nickname[.x]))

    feather::write_feather(beta_stats_pyspi14_whole_brain, glue("{data_path}/pairwise_pyspi14_lm_beta_statistics_by_region_pair.feather"))
} else {
  beta_stats_pyspi14_whole_brain <- feather::read_feather(glue("{data_path}/pairwise_pyspi14_lm_beta_statistics_by_region_pair.feather"))
}