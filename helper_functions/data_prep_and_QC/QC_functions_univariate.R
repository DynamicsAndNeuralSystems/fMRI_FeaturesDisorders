#------------------------------------
# This script compiles movement data across the UCLA subjects
#------------------------------------

#--------------------------------------
# Author: Annie Bryant, 28 August 2022
#--------------------------------------

library(tidyverse)
library(theft)

#-------------------------------------------------------------------------------
# Function to read in univariate TS feature data and return subjects with NA 
# values. For each noise-processing method, the number of TS features with all 
# NA for all brain regions are given.
#-------------------------------------------------------------------------------

find_univariate_sample_na <- function(TS_feature_data, 
                                      dataset_ID = "UCLA_CNP",
                                      univariate_feature_set = "catch22") {
  
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
# Function to read in a univariate TS feature matrix, z-score it, and save
# the z-scored matrix
#-------------------------------------------------------------------------------

z_score_feature_matrix <- function(noise_proc,
                                   TS_feature_data) {
  TS_feature_data_np <- subset(TS_feature_data, Noise_Proc == noise_proc)
  
  TS_feature_data_z <- normalise_feature_frame(TS_feature_data_np, 
                                               names_var = "names",
                                               values_var = "values", 
                                               method = "z-score")
  
  return(TS_feature_data_z)
}

z_score_all_noise_procs <- function(TS_feature_data,
                                    noise_procs = c("AROMA+2P",
                                                    "AROMA+2P+GMR",
                                                    "AROMA+2P+DiCER")) {
  
  TS_feature_data_z <- noise_procs %>%
    purrr::map_df( ~ z_score_feature_matrix(noise_proc = .x,
                                            TS_feature_data = TS_feature_data))
  
  return(TS_feature_data_z)
  
}



#-------------------------------------------------------------------------------
# Plot raw time-series data for subjects with NA values for all ROIs/features
#-------------------------------------------------------------------------------

plot_NA_sample_ts <- function(dataset_ID = "UCLA_CNP",
                              grouping_var = "Brain_Region",
                              raw_TS_file = "UCLA_CNP_fMRI_data.Rds",
                              univariate_feature_set = "catch22",
                              NA_sample_IDs = c(),
                              noise_procs = c("AROMA+2P",
                                              "AROMA+2P+GMR",
                                              "AROMA+2P+DiCER")) {
  
  if (length(NA_sample_IDs) > 0) {
    ts_data <- readRDS(raw_TS_file)
    
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
                                          sample_metadata_file = "UCLA_CNP_sample_metadata.Rds",
                                          dataset_ID = "UCLA_CNP",
                                          univariate_feature_set = "catch22",
                                          raw_TS_file = "UCLA_CNP_fMRI_data.Rds",
                                          add_catch2 = FALSE,
                                          noise_procs = c("AROMA+2P",
                                                          "AROMA+2P+GMR",
                                                          "AROMA+2P+DiCER"),
                                          plot_dir) {
  
  # Load sample metadata
  sample_metadata <- readRDS(paste0(data_path, "study_metadata/", sample_metadata_file))
  
  # Load TS feature data and subset by noise_proc
  TS_feature_data <- readRDS(paste0(rdata_path, dataset_ID, "_", 
                                    univariate_feature_set, ".Rds"))
  
  # Samples identified with missing data for all features:
  univar_NA_samples <- find_univariate_sample_na(TS_feature_data,
                                                 dataset_ID = dataset_ID,
                                                 univariate_feature_set = univariate_feature_set) %>%
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
  
  # Save filtered data to RDS
  saveRDS(TS_feature_data_filtered, file=paste0(proc_rdata_path,
                                                sprintf("%s_%s_filtered.Rds",
                                                        dataset_ID,
                                                        univariate_feature_set)))
  
  # Save sample data post-filtering to an `.Rds` file:
  filtered_sample_info <- TS_feature_data_filtered %>%
    distinct(Sample_ID)                                            
  
  saveRDS(filtered_sample_info, file=paste0(proc_rdata_path, 
                                            sprintf("%s_filtered_sample_info_%s.Rds",
                                                    dataset_ID,
                                                    univariate_feature_set)))
  
  cat("Sample info saved to:", paste0(proc_rdata_path, 
                                      sprintf("%s_filtered_sample_info_%s.Rds",
                                              dataset_ID,
                                              univariate_feature_set)), "\n")
  
  # Data normalisation: z-score the feature matrix as well. 
  TS_df_z <- z_score_all_noise_procs(TS_feature_data = TS_feature_data_filtered,
                                     noise_procs = noise_procs)
  
  saveRDS(TS_df_z, file = paste0(proc_rdata_path, sprintf("%s_%s_filtered_zscored.Rds",
                                                          dataset_ID,
                                                          univariate_feature_set)))
  
  cat("\nZ-scored data saved to:", paste0(proc_rdata_path, sprintf("%s_%s_filtered_zscored.Rds",
                                                                   dataset_ID,
                                                                   univariate_feature_set)),
      "\n")
  
  # OPTIONAL -- if user specifies to add mean and SD, save separately
  if (add_catch2) {
    TS_catch2_data <- readRDS(paste0(proc_rdata_path, dataset_ID, "_catch2.Rds"))
    
    # Filter to samples in metadata
    TS_catch2_filtered <- TS_catch2_data %>%
      dplyr::filter(Sample_ID %in% sample_metadata$Sample_ID)
    
    # Save filtered data to RDS
    saveRDS(TS_catch2_filtered, file=paste0(proc_rdata_path,
                                            sprintf("%s_catch2_filtered.Rds",
                                                    dataset_ID)))
    
    # Data normalisation: z-score the feature matrix as well. 
    TS_catch2_z <- z_score_all_noise_procs(TS_feature_data = TS_catch2_filtered,
                                           noise_procs = noise_procs)
    
    saveRDS(TS_catch2_z, file = paste0(proc_rdata_path, sprintf("%s_catch2_filtered_zscored.Rds",
                                                                dataset_ID)))
  }
}
