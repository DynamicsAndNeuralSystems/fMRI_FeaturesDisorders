# Command-line arguments to parse
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")
parser$add_argument("--github_dir", default="/project/hctsa/annie/github/fMRI_FeaturesDisorders/")
parser$add_argument("--data_path", default="/project/hctsa/annie/data/ABIDE_ASD/")
parser$add_argument("--input_mat_file", default="")
parser$add_argument("--subject_csv", default="participants.csv")
parser$add_argument("--noise_procs", default="FC1000", nargs='?', action='append')
parser$add_argument("--brain_region_lookup", default="Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv")
parser$add_argument("--parcellation_name", default="harvard_oxford_cort_prob_2mm")
parser$add_argument("--dataset_ID", default="ABIDE_ASD")

# Parse input arguments
args <- parser$parse_args()
data_path <- args$data_path
input_mat_file <- args$input_mat_file
subject_csv <- args$subject_csv
noise_procs <- args$noise_procs
parcellation_name <- args$parcellation_name
brain_region_lookup <- args$brain_region_lookup
dataset_ID <- args$dataset_ID
github_dir <- args$github_dir

rdata_path <- paste0(data_path, "Rdata/")
plot_dir <- paste0(data_path, "plots/")

# Set parameters
options(scipen = 999)
set.seed(127)

# load libraries
library(theft)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

# Load Harvard-Oxford atlas lookup table
brain_region_lookup <- read.csv(paste0(data_path, brain_region_lookup))

# Load subject metadata
if (!file.exists(paste0(rdata_path, sprintf("%s_subject_metadata.Rds",
                                            dataset_ID)))) {
  subject_metadata <- read.csv(paste0(data_path, subject_csv)) %>%
    mutate(subject_id = as.character(subject_id)) %>%
    mutate(subject_id = str_replace_all(subject_id, "_", ""))
  
  saveRDS(subject_metadata, file=paste0(rdata_path, sprintf("%s_subject_metadata.Rds",
                                                            dataset_ID)))
} else {
  subject_metadata <- readRDS(paste0(rdata_path, sprintf("%s_subject_metadata.Rds",
                                                         dataset_ID)))
}

# Find subjects with parcellated fMRI time-series data
subjects_with_fMRI <- list.dirs(paste0(data_path, parcellation_name),
                                full.names = F) %>%
  subset(. %in% subject_metadata$subject_id)


# Define function that reads in fMRI time-series data for a given subject
# And reshapes to output data in a tidy (long) format,
# With one row per timepoint/brain region
read_subject_fMRI_data <- function(subject) {
  subject_fMRI_TS <- read.csv(paste0(data_path, sprintf("%s/%s/run_1/%s_task-Rest_confounds.csv",
                                                        parcellation_name, subject, subject)),
                              header = F)
  colnames(subject_fMRI_TS) <- brain_region_lookup$Brain_Region
  subject_fMRI_TS$timepoint <- 1:nrow(subject_fMRI_TS)
  subject_fMRI_TS$Sample_ID <- subject
  subject_fMRI_TS$Noise_Proc <- unique(noise_procs)

  subject_fMRI_TS_long <- subject_fMRI_TS %>%
    pivot_longer(cols=c(-Sample_ID, -timepoint, -Noise_Proc),
                 names_to = "Brain_Region",
                 values_to = "value") %>%
    mutate(Brain_Region = str_replace_all(Brain_Region, " |, ", "_"))

  return(subject_fMRI_TS_long)
}

if (!file.exists(paste0(rdata_path, sprintf("%s_fMRI_data.Rds",
                                            dataset_ID)))) {
  fMRI_data <- subjects_with_fMRI %>%
    purrr::map_df(~ read_subject_fMRI_data(.x))

  saveRDS(fMRI_data, file=paste0(rdata_path, sprintf("%s_fMRI_data.Rds",
                                                     dataset_ID)))
}
