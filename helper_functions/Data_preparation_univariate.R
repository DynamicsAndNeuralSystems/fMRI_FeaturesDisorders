#------------------------------------
# This script sets out to produce a
# function for reading in matlab time
# series files into R
#------------------------------------

#--------------------------------------
# Author: Trent Henderson, 9 March 2021
# Updated: Annie Bryant, 1 August 2022
#--------------------------------------

require(plyr)
library(tidyverse)
library(R.matlab)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Function to load matlab .mat data for UCLA cohort
load_mat_data <- function(mat_file, subject_csv, rdata_path, overwrite=F) {
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Read in data
  cat("\nLoading in .mat file:", mat_file, "\n")
  mat_data <- readMat(mat_file)
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Load noise processing info
  Noise_Proc <- reshape2::melt(mat_data$noiseOptions) %>%
    dplyr::rename("noiseOptions" = "L1",
                  "Noise_Proc" = "value") %>%
    distinct(Noise_Proc, noiseOptions)
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Reshape data
  TS_data_long <- reshape2::melt(mat_data$time.series) %>%
    dplyr::rename(timepoint = Var1,
                  ROI_Index = Var2,
                  Subject_Index = Var3,
                  noiseOptions = Var4)
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Get unique IDs to join in identifiers
  ids <- reshape2::melt(mat_data$subject.list) %>% 
    group_by(value, L1) %>% 
    distinct() %>%
    dplyr::rename(Subject_ID = value,
                  Subject_Index = L1) %>%
    dplyr::select(c(Subject_ID, Subject_Index))
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Retrieve labels and clean up diagnosis names
  subject_info <- read.csv(subject_csv) %>%
    dplyr::rename(Subject_ID = 1) %>%
    distinct(Subject_ID, diagnosis, age, gender) %>%
    semi_join(., ids) %>%
    mutate(diagnosis = str_to_title(diagnosis)) %>%
    mutate(diagnosis = ifelse(diagnosis == "Adhd", "ADHD", diagnosis)) 
  
  # Save CSV file containing list of subjects with time-series data and diagnoses
  write.csv(subject_info, paste0(rdata_path, "UCLA_subjects_with_TS_data.csv"),
            row.names = F)
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Extract ROI names
  ROI_info <- reshape2::melt(mat_data$StructNames) %>%
    distinct(Var1, value) %>%
    dplyr::rename("ROI_Index"="Var1",
                  "Brain_Region"="value") %>%
    mutate(Brain_Region = gsub(" +", "", Brain_Region))
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Bring all data together
  TS_data_full <- inner_join(TS_data_long, ids, 
                             by=c("Subject_Index"="Subject_Index")) %>%
    inner_join(., subject_info, by=c("Subject_ID"="Subject_ID")) %>%
    inner_join(., Noise_Proc, by=c("noiseOptions"="noiseOptions")) %>%
    inner_join(., ROI_info, by=c("ROI_Index"="ROI_Index")) %>%
    dplyr::select(-noiseOptions, -Subject_Index, -ROI_Index) %>%
    filter(diagnosis %in% c("Schz", "Control"))
  
  if (!file.exists(paste0(rdata_path, "UCLA_fMRI_TimeSeries.Rds")) | overwrite) {
    cat("\nWriting UCLA fMRI time-series data to Rds object.", "\n")
    saveRDS(TS_data_full, file=paste0(rdata_path, "UCLA_fMRI_TimeSeries.Rds"))
    
  } else {
    cat("\nUCLA_fMRI_TimeSeries.Rds object already exists and --overwrite was not specified. Not writing new Rds object.\n")
  }
}

