# Command-line arguments to parse
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")
parser$add_argument("--github_dir", default="/project/hctsa/annie/github/fMRI_FeaturesDisorders/")
parser$add_argument("--data_path", default="/project/hctsa/annie/data/UCLA_Schizophrenia/")
parser$add_argument("--input_mat_file", default="", nargs="?")
parser$add_argument("--subject_csv", default="participants.csv")
parser$add_argument("--noise_procs", default=c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER"), nargs='*', action='append')
parser$add_argument("--dataset_ID", default="UCLA_Schizophrenia")

# Parse input arguments
args <- parser$parse_args()
data_path <- args$data_path
input_mat_file <- args$input_mat_file
subject_csv <- args$subject_csv
noise_procs <- args$noise_procs
dataset_ID <- args$dataset_ID
github_dir <- args$github_dir

rdata_path <- paste0(data_path, "Rdata/")
plot_dir <- paste0(github_dir, "plots/")

# Load needed libraries
require(plyr)
library(tidyverse)
library(R.matlab)

#-------------------------------------------------------------------------------
# Function to load matlab .mat data for UCLA cohort
#-------------------------------------------------------------------------------
load_mat_data <- function(mat_file, dataset_ID, subject_csv, rdata_path, overwrite=F) {
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Read in data
  cat("\nLoading in .mat file:", mat_file, "\n")
  mat_data <- readMat(mat_file)
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Load noise processing info
  cat("Cleaning noise-processing info:\n")
  Noise_Proc <- reshape2::melt(mat_data$noiseOptions) %>%
    dplyr::rename("noiseOptions" = "L1",
                  "Noise_Proc" = "value") %>%
    distinct(Noise_Proc, noiseOptions)
  print(Noise_Proc)
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Reshape data
  cat("Reshaping data from wide to long.\n")
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
    dplyr::rename(Sample_ID = value,
                  Subject_Index = L1) %>%
    dplyr::select(c(Sample_ID, Subject_Index))
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Retrieve labels and clean up diagnosis names
  subject_info <- read.csv(subject_csv) %>%
    dplyr::rename(Sample_ID = 1) %>%
    distinct(Sample_ID, diagnosis, age, gender) %>%
    semi_join(., ids) %>%
    mutate(diagnosis = str_to_title(diagnosis)) %>%
    mutate(diagnosis = ifelse(diagnosis == "Adhd", "ADHD", diagnosis)) 
  
  # Save CSV file containing list of subjects with time-series data and diagnoses
  write.csv(subject_info, paste0(rdata_path, sprintf("%s_subjects_with_TS_data.csv",
                                                     dataset_ID)),
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
    inner_join(., subject_info, by=c("Sample_ID"="Sample_ID")) %>%
    inner_join(., Noise_Proc, by=c("noiseOptions"="noiseOptions")) %>%
    inner_join(., ROI_info, by=c("ROI_Index"="ROI_Index")) %>%
    dplyr::select(-noiseOptions, -Subject_Index, -ROI_Index) %>%
    filter(diagnosis %in% c("Schz", "Control"))
  
  # Separate data into TS versus metadata
  metadata <- TS_data_full %>%
    distinct(Sample_ID, diagnosis, age, gender)
  saveRDS(metadata, file=paste0(rdata_path, sprintf("%s_subject_metadata.Rds",
                                                    dataset_ID)))
  
  TS_data_for_analysis <- TS_data_full %>%
    dplyr::select(Sample_ID, Brain_Region, Noise_Proc, timepoint, value) 
  
  if (!file.exists(paste0(rdata_path, sprintf("%s_fMRI_data.Rds",
                                              dataset_ID))) | overwrite) {
    cat("\nWriting", dataset_ID, "fMRI time-series data to Rds object.", "\n")
    saveRDS(TS_data_for_analysis, file=paste0(rdata_path, sprintf("%s_fMRI_data.Rds",
                                                                  dataset_ID)))
    
  } else {
    cat("\nfMRI .Rds object already exists and --overwrite was not specified. Not writing new Rds object.\n")
  }
}

#-------------------------------------------------------------------------------
# Prep data from .mat file
#-------------------------------------------------------------------------------
load_mat_data(mat_file=paste0(data_path, input_mat_file), 
              dataset_ID = dataset_ID,
              subject_csv=paste0(data_path, subject_csv), 
              rdata_path=rdata_path, 
              overwrite=TRUE)