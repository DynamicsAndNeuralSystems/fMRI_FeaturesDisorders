#------------------------------------
# This script sets out to produce a
# function for reading in matlab time
# series files into R
#------------------------------------

#--------------------------------------
# Author: Trent Henderson, 9 March 2021
# Updated: Annie Bryant, 15 March 2022
#--------------------------------------

require(plyr)
library(tidyverse)
library(R.matlab)
library(argparse)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Parse command-line arguments
parser <- ArgumentParser(description='Read in matlab time-series data and convert to an R data object.')
parser$add_argument("--mat_file", help=".mat file containing the time-series data and other metadata.")
parser$add_argument("--label_metadata", help="CSV file containing sample metadata info.",
                    default="D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/participants.csv")
parser$add_argument("--data_path", help="File path to store resulting Rdata objects.",
                    default="D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/")
parser$add_argument("--overwrite", help="Should the Rds object be overwritten if it already exists? Default is F.",
                    action="store_true", default=FALSE)

# Parse arguments
args <- parser$parse_args()
mat_file <- args$mat_file
label_metadata <- args$label_metadata
data_path <- args$data_path
overwrite <- args$overwrite
rdata_path <- paste0(data_path, "Rdata/")

# Create Rdata directory within data_path if it doesn't already exist
if (!dir.exists(rdata_path)) {
  dir.create(rdata_path, showWarnings = FALSE)
}

# DEBUG ONLY
# data_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/"
# rdata_path <- paste0(data_path, "Rdata/")
# mat_file <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/new/UCLA_time_series_four_groups.mat"
# label_metadata <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/participants.csv"
# overwrite <- T

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Function to load matlab .mat data for UCLA cohort

load_mat_data <- function(mat_file, label_metadata, rdata_path, overwrite=F) {
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Read in data
  cat("\nLoading in .mat file:", mat_file, "\n")
  mat_data <- readMat(mat_file)
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Load noise processing info
  noise_proc <- reshape2::melt(mat_data$noiseOptions) %>%
    dplyr::rename("noiseOptions" = "L1",
                  "noise_proc" = "value") %>%
    distinct(noise_proc, noiseOptions)
  
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
  labels <- read.csv(label_metadata) %>%
    dplyr::rename(Subject_ID = 1) %>%
    distinct(Subject_ID, diagnosis, age, gender) %>%
    semi_join(., ids) %>%
    mutate(diagnosis = str_to_title(diagnosis)) %>%
    mutate(diagnosis = ifelse(diagnosis == "Adhd", "ADHD", diagnosis)) 
  
  # Save CSV file containing list of subjects with time-series data and diagnoses
  write.csv(labels, paste0(rdata_path, "UCLA_subjects_with_TS_data.csv"),
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
    inner_join(., labels, by=c("Subject_ID"="Subject_ID")) %>%
    inner_join(., noise_proc, by=c("noiseOptions"="noiseOptions")) %>%
    inner_join(., ROI_info, by=c("ROI_Index"="ROI_Index")) %>%
    dplyr::select(-noiseOptions, -Subject_Index, -ROI_Index) %>%
    filter(diagnosis %in% c("Schz", "Control"))
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Split data by noise processing
  TS_data_split <- split(TS_data_full, TS_data_full$noise_proc)
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Save each dataframe to an Rds object
  for (noise_proc in names(TS_data_split)) {
    tmp <- TS_data_split[[noise_proc]]
    noise_label <- gsub("\\+", "_", noise_proc)
    if (!file.exists(paste0(rdata_path, 
                            sprintf("UCLA_%s.Rds", noise_label))) | 
        overwrite) {
      cat("\nWriting UCLA", noise_proc, "data to Rds object.", "\n")
      saveRDS(tmp, 
              file=paste0(rdata_path, sprintf("UCLA_%s.Rds", noise_label)))
      
    } else {
      cat("\nUCLA", noise_proc, "Rds object already exists and --overwrite was not specified. Not writing new Rds object.\n")
    }
  }
}

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Call function using command-line inputs
load_mat_data(mat_file=mat_file, 
              label_metadata=label_metadata, 
              rdata_path=rdata_path, 
              overwrite=overwrite)