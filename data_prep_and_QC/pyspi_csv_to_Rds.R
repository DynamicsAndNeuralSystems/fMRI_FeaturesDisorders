# Command-line arguments to parse
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")
parser$add_argument("--data_path", default="/project/hctsa/annie/data/UCLA_Schizophrenia/")
parser$add_argument("--brain_region_lookup", default="", nargs="?")
parser$add_argument("--pairwise_feature_set", default="pyspi_19")
parser$add_argument("--noise_procs", default=c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER"), nargs='*', action='append')
parser$add_argument("--dataset_ID", default="UCLA_Schizophrenia")

# Parse input arguments
args <- parser$parse_args()
data_path <- args$data_path
brain_region_lookup <- args$brain_region_lookup
pairwise_feature_set <- args$pairwise_feature_set
dataset_ID <- args$dataset_ID
noise_procs <- args$noise_procs

pydata_path <- paste0(data_path, "pydata/")

# univariate_feature_set <- "catch22"
# subject_csv <- "participants.csv"
# project_path <- "D:/Virtual_Machines/Shared_Folder/github/"

# UCLA schizophrenia
# pydata_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/UCLA_Schizophrenia/pydata/"
# data_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/UCLA_Schizophrenia/"
# dataset_ID <- "UCLA_Schizophrenia"
# noise_procs <- c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER")
# pairwise_feature_set <- "pyspi_19"
# brain_region_lookup <- "Brain_Region_info.csv"

# Load packages
library(tidyverse)

# Read in brain region lookup table
brain_region_LUT <- read.csv(paste0(data_path, brain_region_lookup))

# Function to read in a subject's pyspi data from a CSV and output a dataframe
read_subject_csv <- function(subject_csv, subject_ID) {
  tryCatch({
    subject_data <- read.csv(subject_csv) %>%
      dplyr::mutate(Sample_ID = subject_ID)
    
    return(subject_data)
  }, error = function(e) {
    cat("Could not process data for", subject_ID, "\n")
  })
}

# Iterate over each noise-processing method
for (noise_proc in noise_procs) {
  noise_label = gsub("\\+", "_", noise_proc)
  if (!(file.exists(paste0(pydata_path, noise_label, "_pairwise_", 
                           pairwise_feature_set, ".Rds")))) {
    # Get list of subjects with processed pyspi data
    subjects <- list.dirs(paste0(data_path, "pydata/", noise_label), 
                          recursive = F,
                          full.names = F)
    
    # Iterate over each subject and store pyspi data
    res <- subjects[1:3] %>%
      purrr::map_df(~ read_subject_csv(subject_ID = .x,
                                       subject_csv = paste0(data_path,
                                                            "pydata/", noise_label, 
                                                            "/", .x, "/calc.csv"))) %>%
      dplyr::mutate(Noise_Proc = noise_proc)
    
    saveRDS(res, paste0(pydata_path, noise_label, "_pairwise_", 
                        pairwise_feature_set, ".Rds"))
  }
}






