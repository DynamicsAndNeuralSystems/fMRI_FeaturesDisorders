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

find_univariate_subject_na <- function(rdata_path, 
                                       input_dataset_name = "UCLA",
                                       feature_set = "catch22",
                                       noise_procs = c("AROMA+2P",
                                                       "AROMA+2P+GMR",
                                                       "AROMA+2P+DiCER")) {
  
  TS_feature_data_list <- list()
  for (noise_proc in noise_procs) {
    noise_label <- gsub("\\+", "_", noise_proc)
    TS_feature_df <- readRDS(paste0(rdata_path, sprintf("%s_%s_%s.Rds",
                                                        input_dataset_name,
                                                        noise_label,
                                                        feature_set))) %>%
      mutate(Noise_Proc = noise_proc)
    
    TS_feature_data_list[[noise_label]] <- TS_feature_df
  }
  TS_feature_data <- do.call(plyr::rbind.fill, TS_feature_data_list)
  
  NA_subjects_data <- TS_feature_data %>%
    group_by(Subject_ID, Noise_Proc, names) %>%
    filter(all(is.na(values))) %>%
    ungroup() %>%
    distinct(Subject_ID, Noise_Proc, names) %>%
    group_by(Subject_ID, Noise_Proc) %>%
    dplyr::summarise(num_na = n()) %>%
    # Only want to see subjects with NA for more than one feature
    filter(num_na > 1) %>%
    tidyr::pivot_wider(id_cols=Subject_ID, names_from=Noise_Proc, values_from=num_na) %>%
    mutate_all(~replace(., is.na(.), 0))
  
  return(NA_subjects_data)
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
  
  TS_feature_data_list <- list()
  for (noise_proc in noise_procs) {
    noise_label <- gsub("\\+", "_", noise_proc)
    TS_feature_df <- readRDS(paste0(rdata_path, sprintf("%s_%s_%s.Rds",
                                                        input_dataset_name,
                                                        noise_label,
                                                        feature_set))) %>%
      mutate(Noise_Proc = noise_proc)
    
    TS_feature_data_list[[noise_label]] <- TS_feature_df
  }
  TS_feature_data <- do.call(plyr::rbind.fill, TS_feature_data_list)
  
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

z_score_feature_matrix <- function(rdata_path, 
                                   input_dataset_name = "UCLA",
                                   feature_set = "catch22",
                                   noise_procs = c("AROMA+2P",
                                                   "AROMA+2P+GMR",
                                                   "AROMA+2P+DiCER")) {
  
  for (noise_proc in noise_procs) {
    noise_label <- gsub("\\+", "_", noise_proc)
    TS_feature_df <- readRDS(paste0(rdata_path, sprintf("%s_%s_%s_filtered.Rds",
                                                        input_dataset_name,
                                                        noise_label,
                                                        feature_set))) %>%
      mutate(Noise_Proc = noise_proc)
    
    TS_feature_df_z <- normalise_feature_frame(TS_feature_df, names_var = "names",
                                               values_var = "values", method = "z-score")
    
    saveRDS(TS_feature_df_z, file = paste0(rdata_path, sprintf("%s_%s_%s_filtered_zscored.Rds",
                                                               input_dataset_name,
                                                               noise_label,
                                                               feature_set)))
  }
}

#-------------------------------------------------------------------------------
# Plot raw time-series data for subjects with NA values for all ROIs/features
#-------------------------------------------------------------------------------

plot_NA_subject_ts <- function(rdata_path, 
                               input_dataset_name = "UCLA",
                               feature_set,
                               NA_subject_IDs,
                               noise_procs = c("AROMA+2P",
                                               "AROMA+2P+GMR",
                                               "AROMA+2P+DiCER")) {
  
  ts_data_list <- list()
  for (noise_proc in noise_procs) {
    noise_label <- gsub("\\+", "_", noise_proc)
    ts_df <- readRDS(paste0(rdata_path, sprintf("%s_%s.Rds",
                                                input_dataset_name,
                                                noise_label))) %>%
      filter(Subject_ID %in% NA_subject_IDs) %>%
      mutate(noise_proc = noise_proc)
    
    ts_data_list[[noise_label]] <- ts_df
  }
  ts_data <- do.call(plyr::rbind.fill, ts_data_list)
  
  p <- ts_data%>%
    filter(Subject_ID %in% NA_subject_IDs) %>%
    ggplot(data=., mapping=aes(x=timepoint, y=value, color=Brain_Region)) +
    ggtitle(sprintf("Raw BOLD Signal for %s\nNA subjects with %s",
                    input_dataset_name, feature_set)) +
    geom_line(alpha=0.6) +
    facet_grid(Subject_ID ~ noise_proc, switch="y") +
    theme(legend.position="none",
          strip.text.y.left = element_text(angle=0),
          plot.title = element_text(hjust=0.5))
  
  return(p)
}

#-------------------------------------------------------------------------------
# Function to drop a list of subjects from the given feature matrix
#-------------------------------------------------------------------------------

remove_subjects_from_feature_matrix <- function(rdata_path, 
                                                input_dataset_name = "UCLA",
                                                feature_set = "catch22",
                                                subject_IDs_to_drop,
                                                noise_procs = c("AROMA+2P",
                                                                "AROMA+2P+GMR",
                                                                "AROMA+2P+DiCER")) {
  for (noise_proc in noise_procs) {
    noise_label <- gsub("\\+", "_", noise_proc)
    if (!file.exists(paste0(rdata_path, 
                            sprintf("%s_%s_%s_filtered.Rds",
                                    input_dataset_name,
                                    noise_label,
                                    feature_set)))) {
      TS_feature_df_filtered <- readRDS(paste0(rdata_path, sprintf("%s_%s_%s.Rds",
                                                                   input_dataset_name,
                                                                   noise_label,
                                                                   feature_set))) %>%
        dplyr::mutate(Noise_Proc = noise_proc) %>%
        dplyr::filter(!(Subject_ID %in% subject_IDs_to_drop))
      
      saveRDS(TS_feature_df_filtered, file = paste0(rdata_path, 
                                                    sprintf("%s_%s_%s_filtered.Rds",
                                                            input_dataset_name,
                                                            noise_label,
                                                            feature_set)))
    }
  }
}

#-------------------------------------------------------------------------------
# Function to read in file with average fractional displacement (FD)
# As well as list of individual subject movement files
# And output a CSV containing the Subject ID, diagnosis, and FD
#-------------------------------------------------------------------------------

compile_movement_data <- function(fd_path, 
                                  input_dataset_name,
                                  subject_csv) {
  mov_data <- read.table(paste0(fd_path, sprintf("fdAvgs_%s.txt",
                                                 input_dataset_name)))
  colnames(mov_data)[1] <- "FD"
  
  mov_data_subjects <- list.files(fd_path, pattern="_movData.txt") %>%
    gsub("_movData.txt", "", .)
  
  mov_data$Subject_ID <- mov_data_subjects
  
  subject_info <- read.csv(subject_csv)
  colnames(subject_info)[1] <- "Subject_ID"
  
  mov_data <- mov_data %>%
    left_join(., subject_info) %>%
    dplyr::select(Subject_ID, diagnosis, FD)
  
  return(mov_data)
}

#-------------------------------------------------------------------------------
# Plot FD values by diagnosis boxplot
#-------------------------------------------------------------------------------
plot_FD_vs_diagnosis <- function(movement_data, input_dataset_name) {
  movement_data %>%
    ggplot(data=., mapping=aes(x=diagnosis, y=FD, fill=diagnosis)) +
    geom_boxplot() +
    ylab("Fractional Displacement (FD)") +
    xlab("Diagnosis") +
    ggtitle(sprintf("%s Subject Movement by Diagnosis",
                    input_dataset_name)) +
    theme(legend.position="none",
          plot.title=element_text(hjust=0.5)) 
}

#-------------------------------------------------------------------------------
# Plot the number of control vs. schizophrenia subjects retained per FD threshold
#-------------------------------------------------------------------------------
plot_subjects_per_fd_threshold <- function(movement_data,
                                           input_dataset_name) {
  fd_thresh_list <- list()
  for (fd_threshold in seq(0, 1, by=0.01)) {
    num_ctrl = nrow(subset(movement_data, FD <= fd_threshold & diagnosis=="CONTROL"))
    num_schz = nrow(subset(movement_data, FD <= fd_threshold & diagnosis=="SCHZ"))
    thresh_df <- data.frame(FD_Threshold = fd_threshold,
                            Control = num_ctrl,
                            Schizophrenia = num_schz)
    fd_thresh_list <- rlist::list.append(fd_thresh_list, thresh_df)
  }
  threshold_data <- do.call(plyr::rbind.fill, fd_thresh_list)
  
  threshold_data %>%
    pivot_longer(cols=c(-FD_Threshold),
                 names_to="Group",
                 values_to="# Subjects") %>%
    group_by(Group) %>%
    mutate(`% Subjects` = `# Subjects` / max(`# Subjects`, na.rm=T)) %>%
    pivot_longer(cols = c(`# Subjects`, `% Subjects`)) %>%
    ggplot(data=., mapping=aes(x=FD_Threshold, y=value, color=Group, group=Group)) +
    geom_line(size=2) +
    xlab("FD Threshold") +
    scale_x_reverse() +
    facet_grid(name ~ ., scales="free", switch="both") +
    ggtitle(sprintf("%s Subjects\nRetained by FD Threshold",
                    input_dataset_name)) +
    theme(legend.position="bottom",
          strip.placement = "outside",
          plot.title=element_text(hjust=0.5))
}

#-------------------------------------------------------------------------------
# Plot the ratio of schizophrenia:control subjects retained per FD threshold
#-------------------------------------------------------------------------------

plot_schz_ctrl_ratio_per_fd_threshold <- function(movement_data,
                                                  input_dataset_name) {
  fd_thresh_list <- list()
  for (fd_threshold in seq(0, 1, by=0.01)) {
    num_ctrl = nrow(subset(movement_data, 
                           FD <= fd_threshold & diagnosis=="CONTROL"))
    num_schz = nrow(subset(movement_data, 
                           FD <= fd_threshold & diagnosis=="SCHZ"))
    thresh_df <- data.frame(FD_Threshold = fd_threshold,
                            Control = num_ctrl,
                            Schizophrenia = num_schz)
    fd_thresh_list <- rlist::list.append(fd_thresh_list, thresh_df)
  }
  threshold_data <- do.call(plyr::rbind.fill, fd_thresh_list) %>%
    mutate(SCZ_to_CTRL = Schizophrenia / Control)
  
  threshold_data %>%
    ggplot(data=., mapping=aes(x=FD_Threshold, y=SCZ_to_CTRL)) +
    geom_line(size=2) +
    ylab("SCZ:CTRL Ratio") +
    xlab("FD Threshold") +
    scale_x_reverse() +
    ggtitle(sprintf("Ratio of SCZ:CTRL %s Subjects\nRetained by FD Threshold",
            input_dataset_name)) +
    theme(legend.position="bottom",
          plot.title=element_text(hjust=0.5))
}