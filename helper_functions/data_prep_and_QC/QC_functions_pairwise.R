#------------------------------------
# This script compiles movement data across the UCLA subjects
#------------------------------------

#--------------------------------------
# Author: Annie Bryant, 14 August 2022
#--------------------------------------

library(tidyverse)

#-------------------------------------------------------------------------------
# Function to read in pairwise TS feature data and return subjects with NA 
# values. For each noise-processing method, the number of TS features with all 
# NA for all brain region combinations are given.
#-------------------------------------------------------------------------------

find_pairwise_sample_na <- function(pydata_path, 
                                    input_dataset_name = "UCLA_Schizophrenia",
                                    noise_proc = "AROMA_2P_GMR",
                                    feature_set = "pyspi_19") {
  
  # Load in pairwise data for this noise processing method
  TS_feature_data <- readRDS(paste0(pydata_path, sprintf("%s_pairwise_%s.Rds",
                                                         noise_proc,
                                                         feature_set)))
  
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

find_pairwise_feature_na <- function(pydata_path, 
                                     input_dataset_name = "UCLA_Schizophrenia",
                                     noise_proc = "AROMA_2P_GMR",
                                     feature_set = "pyspi_19") {
  # Load in pairwise data for this noise processing method
  TS_feature_data <- readRDS(paste0(pydata_path, sprintf("%s_pairwise_%s.Rds",
                                                         noise_proc,
                                                         feature_set)))
  
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

z_score_pairwise_feature_matrix <- function(noise_proc,
                                   TS_df) {
  TS_df_noise_proc <- subset(TS_df, Noise_Proc == noise_proc)
  
  TS_feature_df_z <- normalise_feature_frame(TS_df_noise_proc, 
                                             names_var = "names",
                                             values_var = "values", 
                                             method = "z-score")
  
  return(TS_feature_df_z)
}

z_score_all_noise_procs <- function(rdata_path,
                                    input_dataset_name = "UCLA_Schizophrenia",
                                    feature_set = "catch22",
                                    noise_procs = c("AROMA+2P",
                                                    "AROMA+2P+GMR",
                                                    "AROMA+2P+DiCER")) {
  
  TS_df <- readRDS(paste0(rdata_path, sprintf("%s_%s_filtered.Rds",
                                              input_dataset_name,
                                              feature_set)))
  
  TS_df_z <- noise_procs %>%
    purrr::map_df( ~ z_score_feature_matrix(noise_proc = .x,
                                            TS_df = TS_df))
  
  saveRDS(TS_df_z, file = paste0(rdata_path, sprintf("%s_%s_filtered_zscored.Rds",
                                                             input_dataset_name,
                                                             feature_set)))
  
  cat("\nZ-scored data saved to:", paste0(rdata_path, sprintf("%s_%s_filtered_zscored.Rds",
                                                              input_dataset_name,
                                                              feature_set)),
      "\n")
}



#-------------------------------------------------------------------------------
# Plot raw time-series data for subjects with NA values for all ROIs/features
#-------------------------------------------------------------------------------

plot_NA_sample_ts <- function(rdata_path, 
                              input_dataset_name = "UCLA_Schizophrenia",
                              grouping_var = "Brain_Region",
                              raw_TS_file = "UCLA_Schizophrenia_fMRI_data",
                              feature_set = "catch22",
                              NA_sample_IDs = c(),
                              noise_procs = c("AROMA+2P",
                                              "AROMA+2P+GMR",
                                              "AROMA+2P+DiCER")) {
  
  if (length(NA_sample_IDs) > 0) {
    ts_data <- readRDS(raw_TS_file)
    
    ts_data %>%
      filter(Sample_ID %in% NA_sample_IDs) %>%
      ggplot(data=., mapping=aes_string(x="timepoint", y="value", color=grouping_var)) +
      ggtitle(sprintf("Raw time-series for %s\nNA samples with %s",
                      gsub("_", " ", input_dataset_name), feature_set)) +
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

remove_samples_from_feature_matrix <- function(rdata_path, 
                                                input_dataset_name = "UCLA_Schizophrenia",
                                                feature_set = "catch22",
                                                sample_IDs_to_drop = c()) {
  
  cat("\nDropping samples:", paste(sample_IDs_to_drop, collapse=", "), "\n")
  
  TS_df <- readRDS(paste0(rdata_path, sprintf("%s_%s.Rds", input_dataset_name, feature_set))) %>%
    dplyr::filter(!(Sample_ID %in% sample_IDs_to_drop))
  
  saveRDS(TS_df, paste0(rdata_path, sprintf("%s_%s_filtered.Rds", input_dataset_name, feature_set)))
}

#-------------------------------------------------------------------------------
# Function to drop a feature(s) from the given feature matrix
#-------------------------------------------------------------------------------

remove_feature_from_univariate_feature_matrix <- function(rdata_path, 
                                               input_dataset_name = "UCLA_Schizophrenia",
                                               feature_set = "catch22",
                                               features_to_drop = c(),
                                               overwrite = F,
                                               noise_procs = c("AROMA+2P",
                                                               "AROMA+2P+GMR",
                                                               "AROMA+2P+DiCER")) {
  for (noise_proc in noise_procs) {
    noise_label <- gsub("\\+", "_", noise_proc)
    if (!file.exists(paste0(rdata_path, 
                            sprintf("%s_%s_%s_filtered.Rds",
                                    input_dataset_name,
                                    noise_label,
                                    feature_set))) | overwrite) {
      TS_feature_df_filtered <- readRDS(paste0(rdata_path, sprintf("%s_%s_%s.Rds",
                                                                   input_dataset_name,
                                                                   noise_label,
                                                                   feature_set))) %>%
        dplyr::mutate(Noise_Proc = noise_proc) %>%
        dplyr::filter(!(names %in% features_to_drop))
      
      saveRDS(TS_feature_df_filtered, file = paste0(rdata_path, 
                                                    sprintf("%s_%s_%s_filtered.Rds",
                                                            input_dataset_name,
                                                            noise_label,
                                                            feature_set)))
    }
  }
}

#-------------------------------------------------------------------------------
# Function to run quality control methods for univariate data
#-------------------------------------------------------------------------------

run_QC_for_univariate_dataset <- function(rdata_path, 
                               dataset_ID = "UCLA_Schizophrenia",
                               univariate_feature_set = "catch22",
                               raw_TS_file = "UCLA_Schizophrenia_fMRI_data",
                               noise_procs = c("AROMA+2P",
                                               "AROMA+2P+GMR",
                                               "AROMA+2P+DiCER"),
                               plot_dir) {
  
  # Samples identified with missing data for one or more noise-processing methods:
  univar_NA_samples <- find_univariate_sample_na(rdata_path = rdata_path,
                                                   input_dataset_name = dataset_ID,
                                                   feature_set = univariate_feature_set) %>%
    pull(Sample_ID)
  
  # Plot the raw time-series data for these samples to confirm:
  tryCatch({
    plot_NA_sample_ts(rdata_path = rdata_path, 
                      input_dataset_name = dataset_ID,
                      raw_TS_file = raw_TS_file,
                      NA_sample_IDs = univar_NA_samples,
                      noise_procs = noise_procs)
    ggsave(paste0(plot_dir, dataset_ID, "_NA_TimeSeries.png"),
           width = 6, height = 6, units="in", dpi=300)
  }, error = function(e) cat("No NA time-series to plot.\n"))
  
  # Drop any samples shown above with NA features for 
  # one or more noise-processing methods:
  remove_samples_from_feature_matrix(rdata_path, 
                                     input_dataset_name = dataset_ID,
                                     feature_set = univariate_feature_set,
                                     sample_IDs_to_drop = univar_NA_samples)

  # Save sample data post-filtering to an `.Rds` file:
  filtered_sample_info <- readRDS(paste0(rdata_path, 
                                         sprintf("%s_%s_filtered.Rds",
                                                 dataset_ID,
                                                 univariate_feature_set))) %>%
    distinct(Sample_ID)
  saveRDS(filtered_sample_info, file=paste0(rdata_path, 
                                            sprintf("%s_filtered_sample_info_%s.Rds",
                                                    dataset_ID,
                                                    univariate_feature_set)))
  
  cat("Subject info saved to:", paste0(rdata_path, 
                                       sprintf("%s_filtered_sample_info_%s.Rds",
                                               dataset_ID,
                                               univariate_feature_set)), "\n")
  
  # Data normalisation: z-score the feature matrix as well. 
  z_score_all_noise_procs(rdata_path = rdata_path,
                          input_dataset_name = dataset_ID,
                          feature_set = univariate_feature_set,
                          noise_procs = noise_procs)
  
}

#-------------------------------------------------------------------------------
# Function to run quality control methods for univariate data
#-------------------------------------------------------------------------------

run_QC_for_pairwise_dataset <- function(rdata_path, 
                                          dataset_ID = "UCLA_Schizophrenia",
                                          univariate_feature_set = "catch22",
                                          raw_TS_file = "UCLA_Schizophrenia_fMRI_data",
                                          noise_procs = c("AROMA+2P",
                                                          "AROMA+2P+GMR",
                                                          "AROMA+2P+DiCER"),
                                          plot_dir) {
  
  # Samples identified with missing data for one or more noise-processing methods:
  univar_NA_samples <- find_univariate_sample_na(rdata_path = rdata_path,
                                                 input_dataset_name = dataset_ID,
                                                 feature_set = univariate_feature_set) %>%
    pull(Sample_ID)
  
  # Plot the raw time-series data for these samples to confirm:
  tryCatch({
    plot_NA_sample_ts(rdata_path = rdata_path, 
                      input_dataset_name = dataset_ID,
                      raw_TS_file = raw_TS_file,
                      NA_sample_IDs = univar_NA_samples,
                      noise_procs = noise_procs)
    ggsave(paste0(plot_dir, dataset_ID, "_NA_TimeSeries.png"),
           width = 6, height = 6, units="in", dpi=300)
  }, error = function(e) cat("No NA time-series to plot.\n"))
  
  # Drop any samples shown above with NA features for 
  # one or more noise-processing methods:
  remove_samples_from_feature_matrix(rdata_path, 
                                     input_dataset_name = dataset_ID,
                                     feature_set = univariate_feature_set,
                                     sample_IDs_to_drop = univar_NA_samples)
  
  # Save sample data post-filtering to an `.Rds` file:
  filtered_sample_info <- readRDS(paste0(rdata_path, 
                                         sprintf("%s_%s_filtered.Rds",
                                                 dataset_ID,
                                                 univariate_feature_set))) %>%
    distinct(Sample_ID)
  saveRDS(filtered_sample_info, file=paste0(rdata_path, 
                                            sprintf("%s_filtered_sample_info_%s.Rds",
                                                    dataset_ID,
                                                    univariate_feature_set)))
  
  cat("Subject info saved to:", paste0(rdata_path, 
                                       sprintf("%s_filtered_sample_info_%s.Rds",
                                               dataset_ID,
                                               univariate_feature_set)), "\n")
  
  # Data normalisation: z-score the feature matrix as well. 
  z_score_all_noise_procs(rdata_path = rdata_path,
                          input_dataset_name = dataset_ID,
                          feature_set = univariate_feature_set,
                          noise_procs = noise_procs)
  
}