#------------------------------------
# This script compiles movement data across the UCLA subjects
#------------------------------------

#--------------------------------------
# Author: Annie Bryant, 21 March 2022
#--------------------------------------

library(tidyverse)
library(theft)

#-------------------------------------------------------------------------------
# Function to read in univariate TS feature data and return subjects with NA 
# values. For each noise-processing method, the number of TS features with all 
# NA for all brain regions are given.
#-------------------------------------------------------------------------------

find_univariate_sample_na <- function(rdata_path, 
                                       input_dataset_name = "UCLA",
                                       feature_set = "catch22",
                                       noise_procs = c("AROMA+2P",
                                                       "AROMA+2P+GMR",
                                                       "AROMA+2P+DiCER")) {
  TS_feature_data <- readRDS(paste0(rdata_path, sprintf("%s_%s.Rds",
                                                        input_dataset_name,
                                                        feature_set)))
  
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

find_univariate_feature_na <- function(rdata_path, 
                                       input_dataset_name = "UCLA",
                                       feature_set = "catch22",
                                       noise_procs = c("AROMA+2P",
                                                       "AROMA+2P+GMR",
                                                       "AROMA+2P+DiCER")) {
  TS_feature_data <- readRDS(paste0(rdata_path, sprintf("%s_%s.Rds",
                                                        input_dataset_name,
                                                        feature_set)))
  
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
                                   TS_df) {
  TS_df_noise_proc <- subset(TS_df, Noise_Proc == noise_proc)
  
  TS_feature_df_z <- normalise_feature_frame(TS_df_noise_proc, 
                                             names_var = "names",
                                             values_var = "values", 
                                             method = "z-score")
  
  return(TS_feature_df_z)
}

z_score_all_noise_procs <- function(rdata_path,
                                    input_dataset_name = "UCLA",
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
                              input_dataset_name = "UCLA",
                              grouping_var = "Brain_Region",
                              raw_TS_file = "UCLA_fMRI_TimeSeries",
                              feature_set = "catch22",
                              NA_sample_IDs = c(),
                              noise_procs = c("AROMA+2P",
                                              "AROMA+2P+GMR",
                                              "AROMA+2P+DiCER")) {
  
  ts_data <- readRDS(raw_TS_file)
  
  p <- ts_data %>%
    filter(Subject_ID %in% NA_sample_IDs) %>%
    ggplot(data=., mapping=aes_string(x=timepoint, y=value, color=grouping_var)) +
    ggtitle(sprintf("Raw time-series for %s\nNA samples with %s",
                    input_dataset_name, feature_set)) +
    geom_line(alpha=0.6) +
    facet_grid(Subject_ID ~ Noise_Proc, switch="y") +
    theme(legend.position="none",
          strip.text.y.left = element_text(angle=0),
          plot.title = element_text(hjust=0.5))
  
  return(p)
}

#-------------------------------------------------------------------------------
# Function to drop a list of samples from the given feature matrix
#-------------------------------------------------------------------------------

remove_samples_from_feature_matrix <- function(rdata_path, 
                                                input_dataset_name = "UCLA",
                                                feature_set = "catch22",
                                                sample_IDs_to_drop = c()) {
  
  cat("\nDropping samples:", sample_IDs_to_drop, "\n")
  
  TS_df <- readRDS(paste0(rdata_path, sprintf("%s_%s.Rds", input_dataset_name, feature_set))) %>%
    dplyr::filter(!(Subject_ID %in% sample_IDs_to_drop))
  
  saveRDS(TS_df, paste0(rdata_path, sprintf("%s_%s_filtered.Rds", input_dataset_name, feature_set)))
}

#-------------------------------------------------------------------------------
# Function to drop a feature(s) from the given feature matrix
#-------------------------------------------------------------------------------

remove_feature_from_feature_matrix <- function(rdata_path, 
                                               input_dataset_name = "UCLA",
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
# Function to write a QC report
#-------------------------------------------------------------------------------

write_QC_report <- function(rdata_path, 
                            input_dataset_name = "UCLA",
                            univariate_feature_set = "catch22",
                            noise_procs = c("AROMA+2P",
                                            "AROMA+2P+GMR",
                                            "AROMA+2P+DiCER"),
                            plot_dir) {
  
  # Find subjects with NA for all ROIs for all catch22 features
  univar_NA_subjects <- find_univariate_subject_na(rdata_path = rdata_path,
                                                   input_dataset_name = dataset_ID,
                                                   feature_set = univariate_feature_set,
                                                   noise_procs = noise_procs)

  
}