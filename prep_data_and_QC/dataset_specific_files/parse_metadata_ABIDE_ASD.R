# Load needed libraries
library(tidyverse)
library(feather)
library(glue)

# Define paths specific to the UCLA CNP dataset
univariate_feature_set <- "catch24"
pairwise_feature_set <- "pyspi14"
github_dir <- "~/github/fMRI_FeaturesDisorders/"
data_path <- "~/data/ABIDE_ASD/"
output_data_path <- paste0(data_path, "processed_data/")
metadata <- read.csv(glue("{data_path}/study_metadata/Phenotypic_V1_0b_preprocessed1.csv"))
noise_proc <- "GSR"

############################## Save to feather ########################################

if (!file.exists(paste0(data_path,
                        "study_metadata/ABIDE_ASD_sample_metadata.feather"))) {
  # Load UCLA CNP sample metadata
  ABIDE_ASD_QC <- metadata %>%
    mutate(Diagnosis = ifelse(DX_GROUP == 1, "ASD", "Control"),
           Sex = ifelse(SEX == 1, "M", "F")) %>%
    dplyr::rename("Age" = "AGE_AT_SCAN",
                  "Sample_ID" = "FILE_ID",
                  "Site" = "SITE_ID") %>%
    filter(!is.na(Diagnosis)) %>%
    dplyr::select(Sample_ID, Diagnosis, Sex, Age, Site, func_mean_fd, qc_rater_1, qc_func_rater_2, qc_func_rater_3) %>%
    # Apply mean framewise displacement filter
    filter(func_mean_fd < 0.55) %>%
    pivot_longer(cols = c(qc_rater_1:qc_func_rater_3)) %>%
    group_by(Sample_ID) %>%
    # Filter to only subjects where two or more raters passed the functional QC data
    filter(sum(value=="fail") < 2) %>%
    pivot_wider(id_cols = c(Sample_ID:func_mean_fd), names_from=name, values_from=value)
  
  # Save to a joint feather file
  feather::write_feather(ABIDE_ASD_QC, 
                         paste0(data_path,
                                "study_metadata/ABIDE_ASD_sample_metadata.feather"))
} else {
  ABIDE_ASD_QC <- feather::read_feather(paste0(data_path,
                                               "study_metadata/ABIDE_ASD_sample_metadata.feather"))
}