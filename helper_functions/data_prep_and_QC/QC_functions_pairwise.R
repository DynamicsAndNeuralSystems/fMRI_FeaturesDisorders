#--------------------------------------
# Author: Annie Bryant, 14 August 2022
#--------------------------------------

library(tidyverse)
library(theft)

#-------------------------------------------------------------------------------
# Function to read in pairwise TS feature data and return subjects with NA 
# values. For each noise-processing method, the number of TS features with all 
# NA for all brain region combinations are given.
#-------------------------------------------------------------------------------

find_pairwise_sample_na <- function(TS_feature_data, 
                                    dataset_ID = "UCLA_CNP",
                                    noise_proc = "AROMA+2P+GMR",
                                    pairwise_feature_set = "pyspi14") {

  # Load in pairwise data for this noise processing method
  TS_feature_data <- TS_feature_data %>%
    filter(Noise_Proc == noise_proc)
  
  NA_sample_data <- TS_feature_data %>%
    group_by(Sample_ID, Noise_Proc, SPI) %>%
    filter(all(is.na(value))) %>%
    ungroup() %>%
    distinct(Sample_ID, Noise_Proc, SPI) %>%
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

find_pairwise_feature_na <- function(TS_feature_data, 
                                     dataset_ID = "UCLA_CNP",
                                     noise_proc = "AROMA+2P+GMR",
                                     feature_set = "pyspi14") {
  
  # Load in pairwise data for this noise processing method
  TS_feature_data <- TS_feature_data %>%
    filter(Noise_Proc == noise_proc)
  
  NA_feature_data <- TS_feature_data %>%
    group_by(SPI, Noise_Proc) %>%
    filter(all(is.na(value))) %>%
    distinct(SPI, Noise_Proc) 
  
  return(NA_feature_data)
}

#-------------------------------------------------------------------------------
# Function to read in a pairwise TS feature matrix, z-score it, and save
# the z-scored matrix
#-------------------------------------------------------------------------------

z_score_pairwise_feature_matrix <- function(TS_feature_data_filtered, 
                                            proc_rdata_path,
                                            dataset_ID,
                                            pairwise_feature_set,
                                            noise_proc) {

  noise_label = gsub("\\+", "_", noise_proc)

  TS_feature_data_filtered <- subset(TS_feature_data_filtered, 
                                     Noise_Proc == noise_proc) %>%
    dplyr::rename("names"="SPI", "values"="value")
  
  TS_feature_data_z <- normalise_feature_frame(TS_feature_data_filtered, 
                                               names_var = "names",
                                               values_var = "values", 
                                               method = "z-score")
  
  saveRDS(TS_feature_data_z, file = paste0(proc_rdata_path, sprintf("%s_%s_%s_filtered_zscored.Rds",
                                                               dataset_ID,
                                                               noise_label,
                                                               pairwise_feature_set)))
  
  cat("\nZ-scored data saved to:", paste0(proc_rdata_path, sprintf("%s_%s_%s_filtered_zscored.Rds",
                                                                        dataset_ID,
                                                               noise_label,
                                                              pairwise_feature_set)),
      "\n")
}


#-------------------------------------------------------------------------------
# Function to drop a list of samples from the given feature matrix
#-------------------------------------------------------------------------------

remove_samples_from_feature_matrix <- function(TS_feature_data, 
                                               dataset_ID = "UCLA_CNP",
                                               pairwise_feature_set = "pyspi14",
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
                                                dataset_ID = "UCLA_CNP",
                                                pairwise_feature_set = "pyspi14",
                                                features_to_drop = c()) {
    
    cat("\nDropping samples:", paste(features_to_drop, collapse=", "), "\n")
  
  TS_feature_data_filtered <- TS_feature_data %>%
    dplyr::filter(!(SPI %in% features_to_drop))
  
  return(TS_feature_data_filtered)
  
}

#-------------------------------------------------------------------------------
# Function to run quality control methods for univariate data
#-------------------------------------------------------------------------------

run_QC_for_dataset <- function(data_path = "/headnode1/abry4213/data/UCLA_CNP_ABIDE_ASD/", 
                               proc_rdata_path,
                               sample_metadata,
                               dataset_ID = "UCLA_CNP",
                               pairwise_feature_set = "catch22",
                               noise_proc = "AROMA+2P+GMR") {
  
  noise_label = gsub("\\+", "_", noise_proc)

  # Processed rdata path
  if (is.null(proc_rdata_path)) {
    proc_rdata_path <- paste0(data_path, "processed_data/Rdata/")
  }
  
  # Load TS feature data and subset by noise_proc
  TS_feature_data <- readRDS(paste0(proc_rdata_path, dataset_ID, "_", noise_label, "_",
                                    pairwise_feature_set, ".Rds")) %>%
    filter(Noise_Proc == noise_proc)
  
  # Samples identified with missing data for all SPIs:
  pairwise_NA_samples <- find_pairwise_sample_na(TS_feature_data,
                                                 dataset_ID = dataset_ID,
                                                 pairwise_feature_set = pairwise_feature_set,
                                                 noise_proc = noise_proc) %>%
    pull(Sample_ID)
  
  # Drop any samples shown above with NA features for 
  # one or more noise-processing methods:
  TS_feature_data_filtered <- remove_samples_from_feature_matrix(TS_feature_data = TS_feature_data, 
                                                                 dataset_ID = dataset_ID,
                                                                 pairwise_feature_set = pairwise_feature_set,
                                                                 sample_IDs_to_drop = pairwise_NA_samples)
  
  # Samples identified with missing data for all SPIs:
  pairwise_NA_features <- find_pairwise_feature_na(TS_feature_data,
                                                   dataset_ID = dataset_ID,
                                                   noise_proc = noise_proc) %>%
    pull(SPI)
  # Drop any samples shown above with NA features for 
  # one or more noise-processing methods:
  TS_feature_data_filtered <- remove_features_from_feature_matrix(TS_feature_data = TS_feature_data_filtered, 
                                                                 dataset_ID = dataset_ID,
                                                                 pairwise_feature_set = pairwise_feature_set,
                                                                 features_to_drop = pairwise_NA_features)
  
  # Filter to samples in metadata
  TS_feature_data_filtered <- TS_feature_data_filtered %>%
    dplyr::filter(Sample_ID %in% sample_metadata$Sample_ID)
  
  # Save filtered data to RDS
  saveRDS(TS_feature_data_filtered, file=paste0(proc_rdata_path,
                                                sprintf("%s_%s_%s_filtered.Rds",
                                                        dataset_ID,
                                                        noise_label,
                                                        pairwise_feature_set)))

  # Save sample data post-filtering to an `.Rds` file:
  filtered_sample_info <- TS_feature_data_filtered %>%
    distinct(Sample_ID)
  
  saveRDS(filtered_sample_info, file=paste0(proc_rdata_path, 
                                            sprintf("%s_%s_filtered_sample_info_%s.Rds",
                                                    dataset_ID,
                                                    noise_label,
                                                    pairwise_feature_set)))
  
  cat("Sample info saved to:", paste0(proc_rdata_path, 
                                       sprintf("%s_%s_filtered_sample_info_%s.Rds",
                                               dataset_ID,
                                               noise_label,
                                               pairwise_feature_set)), "\n")
  
  # Data normalisation: z-score the feature matrix as well. 
  z_score_pairwise_feature_matrix(TS_feature_data_filtered = TS_feature_data_filtered,
                                  proc_rdata_path = proc_rdata_path, 
                                  dataset_ID = dataset_ID,
                                  pairwise_feature_set = pairwise_feature_set,
                                  noise_proc = noise_proc)
  
}
