# to become argparse
github_dir <- "D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders/"
data_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/"
rdata_path <- paste0(data_path, "Rdata/")
set.seed(127)

# load libraries
library(theft)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

################################################################################
# UCLA data prep
################################################################################

### Source helper function
source(paste0(github_dir, "helper_functions/Data_preparation_univariate.R"))

### prep UCLA data from .mat file
mat_file <- paste0(data_path, "new/UCLA_time_series_four_groups.mat")
subject_csv <- paste0(data_path, "participants.csv")
load_mat_data(mat_file=mat_file, 
              subject_csv=subject_csv, 
              rdata_path=rdata_path, 
              overwrite=FALSE)

################################################################################
# Run catch22 with theft for UCLA dataset
################################################################################

### Source helper function
source(paste0(github_dir, "helper_functions/TS_feature_extraction.R"))

### Define noise-processing methods
noise_procs <- c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER")

catch22_all_regions(rdata_path = rdata_path, 
                    input_dataset = "UCLA",
                    noise_procs = noise_procs)

################################################################################
# Run catchaMouse16 with theft for UCLA dataset
################################################################################

# PLACEHOLDER

################################################################################
# Quality control for UCLA catch22 dataset
################################################################################

### Source helper function
source(paste0(github_dir, "helper_functions/QC_functions.R"))
input_dataset_name = "UCLA"
feature_set = "catch22"

# Find subjects with NA for all ROIs for all catch22 features
UCLA_catch22_NA_subjects <- find_univariate_subject_na(rdata_path = rdata_path,
                                                  input_dataset_name = input_dataset_name,
                                                  feature_set = feature_set,
                                                  noise_procs = noise_procs)
UCLA_catch22_NA_subjects

# We can load in the raw time-series datasets and 
# view the raw values for these subjects
plot_NA_subject_ts(rdata_path = rdata_path, 
                   input_dataset_name = input_dataset_name,
                   NA_subject_IDs = UCLA_catch22_NA_subjects$Subject_ID,
                   noise_procs = noise_procs)

# I will omit these six subjects from catch22 analysis going forward
remove_subjects_from_feature_matrix(rdata_path, 
                                     input_dataset_name = input_dataset_name,
                                     feature_set = feature_set,
                                     subject_IDs_to_drop = UCLA_catch22_NA_subjects$Subject_ID,
                                     noise_procs = c("AROMA+2P",
                                                     "AROMA+2P+GMR",
                                                     "AROMA+2P+DiCER"))

# We can save a dataframe containing the filtered subjects and their diagnoses:
noise_label = "AROMA_2P"

filtered_subject_info <- readRDS(paste0(rdata_path, 
                                        sprintf("%s_%s_%s_filtered.Rds",
                                                input_dataset_name,
                                                noise_label,
                                                feature_set))) %>%
  distinct(Subject_ID, group)
saveRDS(filtered_subject_info, file=paste0(rdata_path, 
                                           sprintf("%s_filtered_subject_info_%s.Rds",
                                                   input_dataset_name,
                                                   feature_set)))

# Lastly, we can z-score normalize the filtered catch22 feature matrix
z_score_feature_matrix(rdata_path = rdata_path,
                       input_dataset_name = input_dataset_name,
                       feature_set = feature_set,
                       noise_procs = noise_procs)

################################################################################
# Quality control for UCLA catchaMouse16 dataset
################################################################################
input_dataset_name = "UCLA"
feature_set = "catchaMouse16"

# First, look for feature(s) that are NA across all subjects:
UCLA_catchaMouse16_NA_features <- find_univariate_feature_na(rdata_path = rdata_path,
                                                             input_dataset_name = input_dataset_name,
                                                             feature_set = feature_set,
                                                             noise_procs = noise_procs)

UCLA_catchaMouse16_NA_features

# The feature `ST_LocalExtrema_n100_diffmaxabsmin` returned all NA values for 
# all subjects and brain regions, so this feature will be dropped.  
# I'll also check to see which subjects have NA for all catchaMouse16 features:

UCLA_catchaMouse16_NA_subjects <- find_univariate_subject_na(rdata_path = rdata_path,
                                                       input_dataset_name = input_dataset_name,
                                                       feature_set = feature_set,
                                                       noise_procs = noise_procs)
UCLA_catchaMouse16_NA_subjects

# We can load in the raw time-series datasets and 
# view the raw values for these subjects
plot_NA_subject_ts(rdata_path = rdata_path, 
                   input_dataset_name = input_dataset_name,
                   feature_set = "catchaMouse16",
                   NA_subject_IDs = UCLA_catchaMouse16_NA_subjects$Subject_ID,
                   noise_procs = noise_procs)

# I will omit these six subjects from catchaMouse16 analysis going forward
remove_subjects_from_feature_matrix(rdata_path, 
                                    input_dataset_name = input_dataset_name,
                                    feature_set = feature_set,
                                    subject_IDs_to_drop = UCLA_catchaMouse16_NA_subjects$Subject_ID,
                                    overwrite = T,
                                    noise_procs = c("AROMA+2P",
                                                    "AROMA+2P+GMR",
                                                    "AROMA+2P+DiCER"))

# Also remove the feature(s) that had no non-NA values for each noise-processing method
remove_feature_from_feature_matrix(rdata_path, 
                                   input_dataset_name = input_dataset_name,
                                   feature_set = feature_set,
                                   features_to_drop = unique(UCLA_catchaMouse16_NA_features$names),
                                   overwrite = T,
                                   noise_procs = c("AROMA+2P",
                                                   "AROMA+2P+GMR",
                                                   "AROMA+2P+DiCER"))

# Lastly, we can z-score normalize the filtered catchaMouse16 feature matrix
z_score_feature_matrix(rdata_path = rdata_path,
                       input_dataset_name = input_dataset_name,
                       feature_set = feature_set,
                       noise_procs = noise_procs)


# We can save a dataframe containing the filtered subjects and their diagnoses:
noise_label = "AROMA_2P"

filtered_subject_info <- readRDS(paste0(rdata_path, 
                                        sprintf("%s_%s_%s_filtered.Rds",
                                                input_dataset_name,
                                                noise_label,
                                                feature_set))) %>%
  distinct(Subject_ID, group)
saveRDS(filtered_subject_info, file=paste0(rdata_path, 
                                           sprintf("%s_filtered_subject_info_%s.Rds",
                                                   input_dataset_name,
                                                   feature_set)))


################################################################################
# Movement
################################################################################

# Compile or load in movement data
if (!file.exists(paste0(rdata_path, sprintf("%s_movement_data.Rds",
                                            input_dataset_name)))) {
  UCLA_movement_data <- compile_movement_data(fd_path=paste0(data_path, "movementData/"),
                                              input_dataset_name = "UCLA",
                                              subject_csv = paste0(data_path, "participants.csv"))
  saveRDS(UCLA_movement_data, file=paste0(rdata_path, sprintf("%s_movement_data.Rds",
                                                              input_dataset_name)))
} else {
  UCLA_movement_data <- readRDS(paste0(rdata_path, sprintf("%s_movement_data.Rds",
                                                           input_dataset_name)))
}

# Visualize the average FD by diagnosis in a boxplot:
plot_FD_vs_diagnosis(movement_data = UCLA_movement_data)

# Plot the number of subjects retained at each FD threshold:
plot_subjects_per_fd_threshold(movement_data = UCLA_movement_data)

# Plot the ratio of schizophrenia:control subjects retained at each FD threshold:
plot_schz_ctrl_ratio_per_fd_threshold(movement_data = UCLA_movement_data)
