# Define paths specific to this dataset
univariate_feature_set <- "catch24"
github_dir <- "~/github/fMRI_FeaturesDisorders/"
data_path <- "~/data/ABIDE_ASD/"
dataset_ID <- "ABIDE_ASD"
noise_procs <- "FC1000"
brain_region_lookup <- "ABIDE_ASD_Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv"
raw_data_input_dir <- paste0(data_path, "raw_data/preprocessed_connectomes_project/harvard_oxford_ROIs/")

# Load needed libraries
library(tidyverse)
library(purrr)
library(feather)

# Load QC-passing metadata
ABIDE_ASD_metadata <- feather::read_feather(paste0(data_path,
                                             "study_metadata/ABIDE_ASD_sample_metadata.feather"))

# Save brain region info to a feather file
ROI_info <- read.csv(paste0(data_path, "study_metadata/", brain_region_lookup))
feather::write_feather(ROI_info, paste0(data_path, sprintf("study_metadata/%s_Brain_Region_Lookup.feather",
                                                          dataset_ID)))

# Define output directory for time-series .txt files
ts_output_dir <- paste0(data_path, "raw_data/time_series_files/GSR/")
TAF::mkdir(ts_output_dir)

# Find list of subjects with rsfMRI data
subjects_with_ts_data <- list.files(raw_data_input_dir, full.names = F, recursive = F) %>%
  gsub("_rois_ho.1D", "", .)

subjects_passing_QC <- ABIDE_ASD_metadata$Sample_ID

subjects_to_parse <- subjects_with_ts_data[subjects_with_ts_data %in% subjects_passing_QC]

# Function to read in data for each subject
output_csv_per_subject <- function(subject_ID) {
  subject_ROI_ts_data <- read.table(glue("{raw_data_input_dir}/{subject_ID}_rois_ho.1D"),
                                    header = F) %>%
    dplyr::select(-V83) # github forums recommend removing this ROI
  output_csv <- paste0(ts_output_dir, subject_ID, "_TS.csv")
  
  write.table(subject_ROI_ts_data, output_csv, row.names = F, col.names = F, sep=",")
}

subjects_to_parse %>%
  purrr::map(~ output_csv_per_subject(.x))
