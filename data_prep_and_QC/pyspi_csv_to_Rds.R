# Command-line arguments to parse
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")
parser$add_argument("--data_path", default="/project/hctsa/annie/data/UCLA_Schizophrenia/")
parser$add_argument("--pairwise_feature_set", default="pyspi_19")
parser$add_argument("--sample_metadata", default="participants.csv")
parser$add_argument("--brain_region_lookup", default="Brain_Region_info.csv")
parser$add_argument("--label_vars", default=c("CONTROL", "SCHIZOPHRENIA"), nargs="*", action="append")
parser$add_argument("--noise_procs", default=c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER"), nargs='*', action='append')
parser$add_argument("--dataset_ID", default="UCLA_Schizophrenia")
parser$add_argument("--overwrite", default=FALSE, action="store_true")

# Parse input arguments
args <- parser$parse_args()
data_path <- args$data_path
label_vars <- args$label_vars
sample_metadata <- args$sample_metadata
brain_region_lookup <- args$brain_region_lookup
pairwise_feature_set <- args$pairwise_feature_set
dataset_ID <- args$dataset_ID
noise_procs <- args$noise_procs
overwrite <- args$overwrite

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

# ABIDE ASD
# pydata_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/ABIDE_ASD/pydata/"
# data_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/ABIDE_ASD/"
# dataset_ID <- "ABIDE_ASD"
# noise_procs <- c("FC1000")
# pairwise_feature_set <- "pyspi_19"
# brain_region_lookup <- "Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv"

# Load packages
library(tidyverse)
library(purrr)

# Read in brain region lookup table
brain_region_LUT <- read.csv(paste0(data_path, brain_region_lookup)) %>%
  mutate(Index = as.numeric(Index))

# Load sample info
sample_metadata_df <- read.csv(paste0(data_path, sample_metadata)) %>%
  mutate(Sample_ID == gsub("_", "", Sample_ID))

# Function to read in a subject's pyspi data from a CSV and output a dataframe
read_sample_pyspi_data <- function(sample_data_file, sample_ID) {
  tryCatch({
    sample_data <- read.csv(sample_data) %>%
      dplyr::mutate(Sample_ID = subject_ID,
                    brain_region_1 = 1+as.numeric(gsub("proc-", "", brain_region_1)),
                    brain_region_2 = 1+as.numeric(gsub("proc-", "", brain_region_2))) %>%
      dplyr::rename("Index" = "brain_region_1") %>%
      left_join(., brain_region_LUT) %>%
      dplyr::select(-Index) %>%
      dplyr::rename("brain_region_1" = "Brain_Region",
                    "Index" = "brain_region_2") %>%
      left_join(., brain_region_LUT) %>%
      dplyr::select(-Index) %>%
      dplyr::rename("brain_region_2" = "Brain_Region")
    
    return(sample_data)
  }, error = function(e) {
    cat("Could not process data for", sample_ID, "\n")
    message(e)
  })
}

# Iterate over each noise-processing method
for (noise_proc in noise_procs) {
  tryCatch({noise_label = gsub("\\+", "_", noise_proc)
  if (!(file.exists(paste0(pydata_path, noise_label, "_pairwise_", 
                           pairwise_feature_set, ".Rds"))) | overwrite) {
    # Get list of subjects with processed pyspi data
    subjects <- list.dirs(paste0(data_path, "pydata/", noise_label), 
                          recursive = F,
                          full.names = F)
    
    # Iterate over each subject and store pyspi data
    res <- subjects %>%
      purrr::map_df(~ read_sample_pyspi_data(sample_ID = .x,
                                       sample_data_file = paste0(data_path,
                                                            "pydata/", noise_label, 
                                                            "/", .x, "/calc.csv"))) %>%
      dplyr::mutate(Noise_Proc = noise_proc)
    
    saveRDS(res, paste0(pydata_path, noise_label, "_pairwise_", 
                        pairwise_feature_set, ".Rds"))
  }
  }, error = function(e) message(e)
  )
}

# Merge all pairwise data for subjects in user-specified groups
if (!(file.exists(paste0(pydata_path, dataset_ID, "_pairwise_", 
                         pairwise_feature_set, ".Rds")))) {
  full_res <- noise_procs %>%
    purrr::map_df(~ readRDS(paste0(pydata_path, gsub("\\+", "_", .x), 
                                   "_pairwise_", 
                                   pairwise_feature_set, ".Rds"))) %>%
    left_join(., sample_metadata_df) %>%
    dplyr::filter(Diagnosis %in% label_vars)
  
  saveRDS(full_res, 
          paste0(pydata_path, dataset_ID, "_pairwise_", 
                 pairwise_feature_set, ".Rds"))
}






