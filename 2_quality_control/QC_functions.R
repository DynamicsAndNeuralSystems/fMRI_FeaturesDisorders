#------------------------------------
# This script compiles movement data across the UCLA subjects
#------------------------------------

#--------------------------------------
# Author: Annie Bryant, 21 March 2022
#--------------------------------------

library(tidyverse)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Function to read in catch22 data and return subjects with NA values

find_catch22_na <- function(rdata_path, noise_procs = c("AROMA+2P",
                                                        "AROMA+2P+GMR",
                                                        "AROMA+2P+DiCER")) {
  
  catch22_data_list <- list()
  for (noise_proc in noise_procs) {
    noise_label <- gsub("\\+", "_", noise_proc)
    catch22_df <- readRDS(paste0(rdata_path, "UCLA_", noise_label, "_catch22.Rds")) %>%
      mutate(noise_proc = noise_proc)
    
    catch22_data_list[[noise_label]] <- catch22_df
  }
  catch22_data <- do.call(plyr::rbind.fill, catch22_data_list)
  
  NA_subjects_data <- catch22_data %>%
    group_by(Subject_ID, noise_proc) %>%
    summarise(num_na = sum(is.na(values))) %>%
    filter(num_na > 0) %>%
    pivot_wider(id_cols=Subject_ID, names_from=noise_proc, values_from=num_na)
  
  return(NA_subjects_data)
}


#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Plot raw time-series data for subjects with catch22 NA values

plot_catch22_na_ts <- function(rdata_path, 
                               subject_IDs,
                               noise_procs = c("AROMA+2P",
                                                        "AROMA+2P+GMR",
                                                        "AROMA+2P+DiCER")) {
  
  ts_data_list <- list()
  for (noise_proc in noise_procs) {
    noise_label <- gsub("\\+", "_", noise_proc)
    ts_df <- readRDS(paste0(rdata_path, "UCLA_", noise_label, ".Rds")) %>%
      filter(Subject_ID %in% subject_IDs) %>%
      mutate(noise_proc = noise_proc)
    
    ts_data_list[[noise_label]] <- ts_df
  }
  ts_data <- do.call(plyr::rbind.fill, ts_data_list)
  
  p <- ts_data%>%
    filter(Subject_ID %in% NA_subjects_data$Subject_ID) %>%
    ggplot(data=., mapping=aes(x=timepoint, y=value, color=Brain_Region)) +
    geom_line(alpha=0.6) +
    facet_grid(Subject_ID ~ noise_proc, switch="y") +
    theme(legend.position="none",
          strip.text.y.left = element_text(angle=0))
  
  return(p)
}

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Function to read in file with average fractional displacement (FD)
# As well as list of individual subject movement files
# And output a CSV containing the Subject ID, diagnosis, and FD

compile_movement_data <- function(fd_path, subject_csv) {
  mov_data <- read.table(paste0(fd_path, "fdAvgs_UCLA.txt"))
  colnames(mov_data) <- "FD"
  
  mov_data_subjects <- list.files(fd_path, pattern="_movData.txt") %>%
    gsub("_movData.txt", "", .)
  
  mov_data$Subject_ID <- mov_data_subjects
  
  subject_info <- read.csv(subject_csv)
  colnames(subject_info)[1] <- "Subject_ID"
  
  mov_data %<>%
    left_join(., subject_info) %>%
    dplyr::select(Subject_ID, diagnosis, FD)
  
  return(mov_data)
}

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Plot FD values by diagnosis boxplot
plot_FD_vs_diagnosis <- function(movement_data) {
  movement_data %>%
    ggplot(data=., mapping=aes(x=diagnosis, y=FD, fill=diagnosis)) +
    geom_boxplot() +
    ylab("Fractional Displacement (FD)") +
    xlab("Diagnosis") +
    ggtitle("Subject Movement by Diagnosis") +
    theme(legend.position="none",
          plot.title=element_text(hjust=0.5)) 
}

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Plot the number of control vs. schizophrenia subjects retained per FD threshold
plot_subjects_per_fd_threshold <- function(movement_data) {
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
                 values_to="n") %>%
    ggplot(data=., mapping=aes(x=FD_Threshold, y=n, color=Group, group=Group)) +
    geom_line() +
    ylab("# Subjects") +
    xlab("FD Threshold") +
    scale_x_reverse() +
    ggtitle("Number of Subjects Retained\nby FD Threshold") +
    theme(legend.position="bottom",
          plot.title=element_text(hjust=0.5))
}
