# Define paths specific to this dataset
univariate_feature_set <- "catch22"
subject_csv <- "UCLA_CNP_participants.csv"
github_dir <- "~/github/fMRI_FeaturesDisorders/"
data_path <- "~/data/UCLA_CNP/"
dataset_ID <- "UCLA_CNP"
input_mat_file = "UCLA_time_series_four_groups.mat"
noise_procs <- c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER")
subject_csv = "UCLA_CNP_participants.csv"

# Load needed libraries
require(plyr)
library(tidyverse)
library(R.matlab)
library(arrow)

# Define output directory for time-series .txt files
ts_output_dir <- paste0(data_path, "raw_data/time_series_files/")

#-------------------------------------------------------------------------------
# Function to load matlab .mat data and output time-series data as .txt files
#-------------------------------------------------------------------------------
mat_data_into_TXT_files <- function(input_mat_file, 
                                    dataset_ID, 
                                    subject_csv, 
                                    data_path, 
                                    overwrite=F) {
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Read in data
  cat("\nLoading in .mat file:", input_mat_file, "\n")
  mat_data <- readMat(input_mat_file)
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Load noise processing info
  cat("Cleaning noise-processing info:\n")
  Noise_Proc <- data.frame(Noise_Proc = unlist(mat_data$noiseOptions)) %>%
    mutate(noiseOptions = 1:n())
  print(Noise_Proc)
  
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Reshape data
  cat("Reshaping data from wide to long.\n")
  TS_data_long <- reshape2::melt(mat_data$time.series) %>%
    dplyr::rename(timepoint = Var1,
                  Index = Var2,
                  Subject_Index = Var3,
                  noiseOptions = Var4)
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Get unique IDs to join in identifiers
  ids <- data.frame(Sample_ID = unlist(mat_data$subject.list)) %>% 
    mutate(Subject_Index = dplyr::row_number())
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Retrieve labels and clean up diagnosis names
  subject_info <- read.csv(subject_csv) %>%
    dplyr::rename(Sample_ID = 1) %>%
    distinct(Sample_ID, diagnosis, age, gender) %>%
    semi_join(., ids) %>%
    mutate(diagnosis = str_to_title(diagnosis)) %>%
    mutate(diagnosis = ifelse(diagnosis == "Adhd", "ADHD", diagnosis)) %>%
    dplyr::rename("Diagnosis" = "diagnosis",
                  "Age" = "age",
                  "Sex" = "gender") %>%
    dplyr::select(Sample_ID:Sex)
  
  # Save feather file containing list of subjects with time-series data and diagnoses
  arrow::write_feather(subject_info, paste0(data_path, sprintf("processed_data/%s_subjects_with_fMRI_TS_data.feather",
                                                               dataset_ID)))
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Extract ROI names
  ROI_info <- reshape2::melt(mat_data$StructNames) %>%
    distinct(Var1, value) %>%
    dplyr::rename("Index"="Var1",
                  "Brain_Region"="value") %>%
    mutate(Brain_Region = gsub(" +", "", Brain_Region))
  
  # Save feather mapping index to brain region name
  arrow::write_feather(ROI_info, paste0(data_path, sprintf("study_metadata/%s_Brain_Region_Lookup.feather",
                                                           dataset_ID)))
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Bring all data together
  TS_data_full <- inner_join(TS_data_long, ids, 
                             by=c("Subject_Index"="Subject_Index")) %>%
    inner_join(., subject_info, by=c("Sample_ID"="Sample_ID")) %>%
    inner_join(., Noise_Proc, by=c("noiseOptions"="noiseOptions")) %>%
    dplyr::select(-noiseOptions, -Subject_Index)
  
  # Separate data into TS versus metadata
  metadata <- TS_data_full %>%
    distinct(Sample_ID, Diagnosis, Age, Sex)
  arrow::write_feather(metadata, paste0(data_path, sprintf("study_metadata/%s_sample_metadata.feather",
                                                           dataset_ID)))
  
  for (noise_proc in Noise_Proc$Noise_Proc) {
    noise_label = gsub("\\+", "_", noise_proc)
    
    noise_proc_subset = TS_data_full %>%
      dplyr::filter(Noise_Proc == noise_proc)
    
    # Make output directory
    np_output_dir = paste0(ts_output_dir, noise_label, "/")
    TAF::mkdir(np_output_dir)
    
    # Save a .csv file per sample
    for (sample in unique(noise_proc_subset$Sample_ID)) {
      if (!file.exists(paste0(np_output_dir, sample, "_TS.csv"))) {
        noise_proc_subset %>%
          filter(Sample_ID == sample) %>%
          dplyr::select(timepoint, Index, value) %>%
          pivot_wider(names_from=Index, values_from=value) %>%
          dplyr::select(-timepoint) %>%
          write.csv(., 
                    file = paste0(np_output_dir, sample, "_TS.csv"),
                    col.names = F, row.names=F)
      }
    }
    
  }
}

#-------------------------------------------------------------------------------
# Prep data from .mat file
#-------------------------------------------------------------------------------
mat_data_into_TXT_files(input_mat_file=paste0(data_path, 
                                              "raw_data/",
                                              input_mat_file), 
                        dataset_ID = dataset_ID,
                        subject_csv = paste0(data_path, "study_metadata/", subject_csv), 
                        data_path = data_path)