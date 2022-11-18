################################################################################
# Load libraries
################################################################################

library(tidyverse)
library(icesTAF)
library(cowplot)
theme_set(theme_cowplot())
library(knitr)
library(kableExtra)

################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
source(paste0(github_dir, "helper_functions/classification/Linear_SVM.R"))
source(paste0(github_dir, "helper_functions/data_prep_and_QC/QC_functions_univariate.R"))
source(paste0(github_dir, "helper_functions/data_prep_and_QC/QC_functions_pairwise.R"))
plot_path <- paste0(github_dir, "plots/QC/")
icesTAF::mkdir(plot_path)

SCZ_data_path <- "~/data/UCLA_Schizophrenia/"
SCZ_rdata_path <- paste0(SCZ_data_path, "processed_data/Rdata/")
ASD_data_path <- "~/data/ABIDE_ASD/"
ASD_rdata_path <- paste0(ASD_data_path, "processed_data/Rdata/")

# Load subject metadata
SCZ_subject_metadata <- readRDS(paste0(SCZ_data_path, "UCLA_Schizophrenia_sample_metadata.Rds")) %>%
  dplyr::filter(Diagnosis %in% c("Control", "Schizophrenia"))
ASD_subject_metadata <- readRDS(paste0(ASD_data_path, "ABIDE_ASD_sample_metadata.Rds"))

univariate_feature_set <- "catch22"
pairwise_feature_set <- "pyspi14"

SCZ_noise_procs <- c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER")
ASD_noise_procs <- c("FC1000")

# Load unfiltered catch22 per dataset
SCZ_catch22 <- readRDS(paste0(SCZ_rdata_path, "UCLA_Schizophrenia_catch22.Rds")) %>%
  left_join(., SCZ_subject_metadata) %>%
  filter(Diagnosis %in% c("Control", "Schizophrenia"))
ASD_catch22 <- readRDS(paste0(ASD_rdata_path, "ABIDE_ASD_catch22.Rds")) %>%
  left_join(., ASD_subject_metadata) %>%
  filter(Diagnosis %in% c("Control", "ASD"))

################################################################################
# Univariate data
################################################################################

# Find SCZ subjects with NA values for one or more catch22 features
SCZ_NA_subjects <- find_univariate_sample_na(TS_feature_data = SCZ_catch22,
                                             dataset_ID = "UCLA_Schizophrenia",
                                             univariate_feature_set = "catch22")
SCZ_NA_subjects %>%
  left_join(., SCZ_subject_metadata %>% dplyr::select(Sample_ID, Diagnosis)) %>%
  dplyr::select(Sample_ID, Diagnosis, `AROMA+2P`:`AROMA+2P+GMR`) %>%
  kable() %>%
  kable_styling(full_width=F)

# Plot their time-series
plot_NA_sample_ts(dataset_ID = "UCLA_Schizophrenia",
                  grouping_var = "Brain_Region",
                  raw_TS_file = paste0(SCZ_data_path, "raw_data/UCLA_Schizophrenia_fMRI_TS.Rds"),
                  univariate_feature_set = univariate_feature_set,
                  NA_sample_IDs = SCZ_NA_subjects$Sample_ID,
                  noise_procs = SCZ_noise_procs)
ggsave(paste0(plot_path, 
              sprintf("SCZ_Raw_fMRI_Signal_NA_Subjects_Univariate_%s.png", 
                      univariate_feature_set)),
       width=7, height=6, units="in", dpi=300, bg="white")

# Find catch22 features that are missing for most/all SCZ participants
SCZ_NA_features <- find_univariate_feature_na(TS_feature_data = SCZ_catch22,
                                              dataset_ID = "UCLA_Schizophrenia",
                                              univariate_feature_set = "catch22")
SCZ_NA_features %>%
  kable() %>%
  kable_styling(full_width=F)

# Find SCZ subjects with NA values for one or more catch22 features
ASD_NA_subjects <- find_univariate_sample_na(TS_feature_data = ASD_catch22,
                                             dataset_ID = "ABIDE_ASD",
                                             univariate_feature_set = "catch22")
ASD_NA_subjects %>%
  kable() %>%
  kable_styling(full_width=F)

# Find catch22 features that are missing for most/all ASD participants
ASD_NA_features <- find_univariate_feature_na(TS_feature_data = ASD_catch22,
                                              dataset_ID = "ABIDE_ASD",
                                              univariate_feature_set = "catch22")
ASD_NA_features %>%
  kable() %>%
  kable_styling(full_width=F)
