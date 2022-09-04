# Define paths specific to this dataset
univariate_feature_set <- "catch22"
subject_csv <- "HCP100_subject_info.csv"
github_dir <- "/headnode1/abry4213/github/fMRI_FeaturesDisorders/"
data_path <- "/headnode1/abry4213/data/HCP100/"
dataset_ID <- "HCP100"
noise_procs <- c("AROMA+2P+GMR")

# Load needed libraries
require(plyr)
library(tidyverse)
library(R.matlab)
library(purrr)

# Iterate over noise procs
for (noise_proc in noise_procs) {
  # Get noise proc label
  noise_label <- gsub("\\+", "_", noise_proc)
  
  # Define output directory for time-series .txt files
  ts_output_dir <- paste0(data_path, "raw_data/time_series_files/", noise_label, "/")
  
  # Function to write DK atlas data to a CSV file per subject
  TS_data_to_CSV <- function(sample_ID, input_data_path, output_data_path) {
    sample_data <- R.matlab::readMat(paste0(input_data_path, sample_ID, "/cfg.mat"))
    TS_data <- sample_data[[2]]
    DK_TS_data <- as.data.frame(TS_data[[1]])
    # Write Desikan-Killiany parcellation data to a CSV file
    write.table(DK_TS_data, paste0(output_data_path, sample_ID, "_TS.csv"), 
                row.names = F, col.names = F, sep=",")
  }
  
  # Define samples to process
  samples_to_process <- list.dirs(paste0(data_path, "raw_data/cfg_data/"),
                                  recursive=F, full.names = F)
  
  # Iterate over samples to process
  samples_to_process %>%
    purrr::map_df( ~ TS_data_to_CSV(sample_ID = .x,
                                    input_data_path = paste0(data_path, "raw_data/cfg_data/"),
                                    output_data_path = ts_output_dir
    ))
}

# Save metadata
metadata <- read.csv(paste0(data_path, subject_csv)) %>%
  mutate(Sample_ID = gsub("_", "", Sample_ID)) 
saveRDS(metadata, file=paste0(data_path, sprintf("%s_sample_metadata.Rds",
                                                 dataset_ID)))
