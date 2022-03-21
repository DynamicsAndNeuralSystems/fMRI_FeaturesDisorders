#------------------------------------
# This script compiles movement data across the UCLA subjects
#------------------------------------

#--------------------------------------
# Author: Annie Bryant, 21 March 2022
#--------------------------------------

library(tidyverse)

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
  
  write.csv(mov_data, file=paste0(fd_path, "fdAvgs_UCLA.csv"),
            row.names = F)
}

