#------------------------------------
# This script compiles movement data across the UCLA subjects
#------------------------------------

#--------------------------------------
# Author: Annie Bryant, 28 August 2022
#--------------------------------------

library(tidyverse)
library(theft)
library(feather)

#-------------------------------------------------------------------------------
# Function to read in univariate TS feature data and return subjects with NA 
# values. For each noise-processing method, the number of TS features with all 
# NA for all brain regions are given.
#-------------------------------------------------------------------------------

find_univariate_sample_na <- function(TS_feature_data, 
                                      dataset_ID = "UCLA_CNP") {
  
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
# Function to read in univariate TS feature data and return features with 
# NA for all subjects/brain regions
#-------------------------------------------------------------------------------

find_univariate_feature_na <- function(TS_feature_data, 
                                       dataset_ID = "UCLA_CNP",
                                       univariate_feature_set = "catch22") {
  
  NA_feature_data <- TS_feature_data %>%
    group_by(names, Noise_Proc) %>%
    filter(all(is.na(values))) %>%
    distinct(names, Noise_Proc) 
  
  return(NA_feature_data)
}

#-------------------------------------------------------------------------------
# Plot raw time-series data for subjects with NA values for all ROIs/features
#-------------------------------------------------------------------------------

plot_NA_sample_ts <- function(dataset_ID = "UCLA_CNP",
                              grouping_var = "Brain_Region",
                              raw_TS_file = "UCLA_CNP_AROMA_2P_GMR_fMRI_data.feather",
                              univariate_feature_set = "catch22",
                              NA_sample_IDs = c(),
                              noise_proc = "AROMA+2P+GMR") {
  
  if (length(NA_sample_IDs) > 0) {
    ts_data <- feather::read_feather(raw_TS_file)
    
    ts_data %>%
      filter(Sample_ID %in% NA_sample_IDs) %>%
      ggplot(data=., mapping=aes_string(x="timepoint", y="values", color=grouping_var)) +
      ggtitle(sprintf("Raw time-series for %s\nNA samples with %s",
                      gsub("_", " ", dataset_ID), univariate_feature_set)) +
      geom_line(alpha=0.6) +
      facet_grid(Sample_ID ~ Noise_Proc, switch="y") +
      theme(legend.position="none",
            strip.text.y.left = element_text(angle=0),
            plot.title = element_text(hjust=0.5))
  } else {
    cat("No NA data to show.\n")
  }
}

#-------------------------------------------------------------------------------
# Function to drop a list of samples from the given feature matrix
#-------------------------------------------------------------------------------

remove_samples_from_feature_matrix <- function(TS_feature_data, 
                                               sample_IDs_to_drop = c()) {
  
  cat("\nDropping samples:", paste(sample_IDs_to_drop, collapse=", "), "\n")
  
  TS_feature_data_filtered <- TS_feature_data %>%
    dplyr::filter(!(Sample_ID %in% sample_IDs_to_drop))
  
  return(TS_feature_data_filtered)
  
}

#-------------------------------------------------------------------------------
# Function to drop a feature(s) from the given feature matrix
#-------------------------------------------------------------------------------

remove_features_from_feature_matrix <- function(TS_feature_data, 
                                                features_to_drop = c()) {
  if (length(features_to_drop) > 0) {
    cat("\nDropping features:", paste(names, collapse=", "), "\n")
    
    TS_feature_data_filtered <- TS_feature_data %>%
      dplyr::filter(!(names %in% features_to_drop))
    
    return(TS_feature_data_filtered)
  } else {
    return(TS_feature_data)
  }
  
}

#-------------------------------------------------------------------------------
# Function to run quality control methods for univariate data
#-------------------------------------------------------------------------------

run_QC_for_univariate_dataset <- function(data_path, 
                                          proc_rdata_path,
                                          sample_metadata_file = "UCLA_CNP_sample_metadata.feather",
                                          dataset_ID = "UCLA_CNP",
                                          univariate_feature_set = "catch22",
                                          raw_TS_file = "UCLA_CNP_AROMA_2P_GMR_fMRI_data.feather",
                                          noise_proc = "AROMA+2P+GMR",
                                          plot_dir) {
  
  noise_label = gsub("\\+", "_", noise_proc)
  
  # Load sample metadata
  sample_metadata <- feather::read_feather(paste0(data_path, "study_metadata/", sample_metadata_file))
  
  # Load TS feature data and subset by noise_proc
  TS_feature_data <- feather::read_feather(paste0(rdata_path, dataset_ID, "_", noise_label, "_", 
                                                univariate_feature_set, ".feather")) %>%
    mutate(Feature_Set = univariate_feature_set)
  
  # Samples identified with missing data for all features:
  univar_NA_samples <- find_univariate_sample_na(TS_feature_data,
                                                 dataset_ID = dataset_ID) %>%
    pull(Sample_ID)
  
  # Drop any samples shown above with NA features for 
  # one or more noise-processing methods:
  TS_feature_data_filtered <- remove_samples_from_feature_matrix(TS_feature_data = TS_feature_data, 
                                                                 sample_IDs_to_drop = univar_NA_samples)
  
  # Features identified with missing data for all samples:
  univar_NA_features <- find_univariate_feature_na(TS_feature_data_filtered,
                                                   dataset_ID = dataset_ID,
                                                   univariate_feature_set = univariate_feature_set) %>%
    pull(names)
  
  # Drop any samples shown above with NA features for 
  # one or more noise-processing methods:
  TS_feature_data_filtered <- remove_features_from_feature_matrix(TS_feature_data = TS_feature_data_filtered, 
                                                                  features_to_drop = univar_NA_features)                                                               
  
  
  # Filter to samples in metadata
  TS_feature_data_filtered <- TS_feature_data_filtered %>%
    dplyr::filter(Sample_ID %in% sample_metadata$Sample_ID)
  
  # Save filtered data to feather
  feather::write_feather(TS_feature_data_filtered, paste0(proc_rdata_path,
                                                        sprintf("%s_%s_%s_filtered.feather",
                                                                dataset_ID,
                                                                noise_label,
                                                                univariate_feature_set)))
  
  # Save sample data post-filtering to a feather file:
  filtered_sample_info <- TS_feature_data_filtered %>%
    distinct(Sample_ID)                                            
  
  feather::write_feather(filtered_sample_info, paste0(proc_rdata_path, 
                                                    sprintf("%s_filtered_sample_info_%s_%s.feather",
                                                            dataset_ID,
                                                            noise_label,
                                                            univariate_feature_set)))
  
  cat("Sample info saved to:", paste0(proc_rdata_path, 
                                      sprintf("%s_filtered_sample_info_%s_%s.feather",
                                              dataset_ID,
                                              noise_label,
                                              univariate_feature_set)), "\n")
}
