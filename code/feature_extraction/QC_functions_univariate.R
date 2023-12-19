#------------------------------------
# This script compiles movement data across the UCLA subjects
#------------------------------------

#--------------------------------------
# Author: Annie Bryant, 28 August 2022
#--------------------------------------

library(tidyverse)
library(theft)
library(feather)
library(glue)

#-------------------------------------------------------------------------------
# Function to read in univariate TS feature data and return subjects with NA 
# values. For each noise-processing method, the number of TS features with all 
# NA for all brain regions are given.
#-------------------------------------------------------------------------------

find_univariate_sample_na <- function(TS_feature_data) {
  
  NA_sample_data <- TS_feature_data %>%
    group_by(Sample_ID, Noise_Proc, names) %>%
    filter(all(is.na(values))) %>%
    ungroup() %>%
    distinct(Sample_ID, Noise_Proc, names) %>%
    group_by(Sample_ID, Noise_Proc) %>%
    dplyr::summarise(num_na = n()) %>%
    # Only want to see subjects with NA for more than one feature
    filter(num_na > 1) %>%
    tidyr::pivot_wider(id_cols=Sample_ID, names_from=Noise_Proc, values_from=num_na) %>%
    mutate_all(~replace(., is.na(.), 0))
  
  return(NA_sample_data)
}

#-------------------------------------------------------------------------------
# Function to run quality control methods for univariate data
#-------------------------------------------------------------------------------

run_QC_for_univariate_dataset <- function(univariate_feature_set = "catch25",
                                          sample_metadata,
                                          catch24_results,
                                          fALFF_results,
                                          participants_to_drop_file) {
  
  
  # # Load TS feature data and subset by noise_proc
  # theft_processed_data <- pyarrow_feather$read_feather(glue("{data_path}/processed_data/{theft_processed_data_file}")) %>%
  #   filter(Noise_Proc == noise_proc) 
  # fALFF_processed_data <- pyarrow_feather$read_feather(glue("{data_path}/processed_data/{fALFF_processed_data_file}")) %>%
  #   mutate(Sample_ID = trimws(Sample_ID),
  #          Brain_Region = trimws(Brain_Region))

  # Merge theft and fALFF data
  TS_feature_data <- catch24_results %>%
    plyr::rbind.fill(fALFF_results) %>%
    dplyr::mutate(Feature_Set = univariate_feature_set)
    
  
  # Samples identified with missing data for all features:
  univar_NA_samples <- find_univariate_sample_na(TS_feature_data) %>%
    pull(Sample_ID)
  
  # Drop any samples shown above with NA features for 
  # one or more noise-processing methods:
  TS_feature_data_filtered <- subset(TS_feature_data, !(Sample_ID %in% univar_NA_samples))
  
  # Filter to samples in metadata
  TS_feature_data_filtered <- TS_feature_data_filtered %>%
    dplyr::filter(Sample_ID %in% sample_metadata$Sample_ID)

  # Find subjects to drop based on head movement
  lenient_FD_subjects_to_drop <- unname(unlist(read.table(participants_to_drop_file, header = TRUE, 
    sep = ",", stringsAsFactors = FALSE, colClasses = "character")[1]))
  TS_feature_data_filtered <- TS_feature_data_filtered %>%
    dplyr::filter(!(Sample_ID %in% lenient_FD_subjects_to_drop))

    # Return the filtered data
  return(TS_feature_data_filtered)
}
