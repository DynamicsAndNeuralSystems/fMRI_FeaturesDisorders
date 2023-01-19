# Define paths specific to this dataset
univariate_feature_set <- "catch22"
participant_csv <- "ABIDE_ASD_participants.csv"
github_dir <- "~/github/fMRI_FeaturesDisorders/"
data_path <- "~/data/ABIDE_ASD/"
dataset_ID <- "ABIDE_ASD"
noise_procs <- "FC1000"
raw_data_input_dir <- paste0(data_path, "raw_data/harvard_oxford_cort_prob_2mm/")

# Load needed libraries
library(tidyverse)
library(purrr)
library(feather)

# Define output directory for time-series .txt files
ts_output_dir <- paste0(data_path, "raw_data/time_series_files/FC1000/")
TAF::mkdir(ts_output_dir)

# Find list of subjects with rsfMRI data
subjects <- list.dirs(raw_data_input_dir, full.names = F, recursive = F)

# Function to read in data for each subject
output_csv_per_subject <- function(subject_ID) {
  subject_csv <- paste0(raw_data_input_dir, subject_ID, 
                        "/run_1/", subject_ID, "_task-Rest_confounds.csv")
  output_csv <- paste0(ts_output_dir, subject_ID, "_TS.csv")
  
  file.copy(subject_csv, output_csv)
}

subjects %>%
  purrr::map(~ output_csv_per_subject(.x))

# Save metadata
metadata <- read.csv(paste0(data_path, "study_metadata/", participant_csv),
                     colClasses = "character") %>%
  dplyr::rename("Sample_ID" = "subject_id",
                "Site" = "site",
                "Sex" = "sex",
                "Age" = "age",
                "ASD" = "asd")
feather::write_feather(metadata, paste0(data_path, sprintf("study_metadata/%s_sample_metadata.feather",
                                                         dataset_ID)))