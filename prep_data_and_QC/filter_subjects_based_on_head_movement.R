#-------------------------------------------------------------------------------
# Define paths
#-------------------------------------------------------------------------------
# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--data_path", default="~/data/UCLA_CNP/")
parser$add_argument("--sample_metadata_file", default="UCLA_CNP_sample_metadata.feather")
parser$add_argument("--dataset_ID", default="UCLA_CNP")

# Parse input arguments
args <- parser$parse_args()
data_path <- args$data_path
sample_metadata_file <- args$sample_metadata_file
dataset_ID <- args$dataset_ID

# Load libraries
library(tidyverse)
library(glue)

# Load mean framewise displacement data
mean_FD <- read.table(glue("{data_path}/movement_data/fmriprep/{dataset_ID}_mFD.txt"),
                      sep=",", colClasses = "character")

colnames(mean_FD) <- c("Sample_ID", "Jenkinson", "Power", "VanDijk")

mean_FD <- mean_FD %>%
  dplyr::select(Sample_ID, Power) %>%
  mutate(Power = as.numeric(Power))

# Find subjects to drop based on the "lenient" routine from Parkes et al 2018 -- mFD Power < 0.55 for all retained subjects
subjects_to_drop <- mean_FD %>%
  filter(Power > 0.55) %>%
  pull(Sample_ID)

# Save this list of subjects to a .txt file
write.table(subjects_to_drop, 
            file=glue("{data_path}/movement_data/fmriprep/{dataset_ID}_subjects_to_drop_lenient.txt"),
            row.names = F,
            col.names = F,
            quote=F)