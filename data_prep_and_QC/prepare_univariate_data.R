#-------------------------------------------------------------------------------
# Define paths
#-------------------------------------------------------------------------------
# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--github_dir", default="/headnode1/abry4213/github/")
parser$add_argument("--data_path", default="/headnode1/abry4213/data/UCLA_Schizophrenia/")
parser$add_argument("--univariate_feature_set", default="catch22")
parser$add_argument("--sample_metadata_file", default="UCLA_Schizophrenia_sample_metadata.Rds")
parser$add_argument("--brain_region_lookup", default="", nargs='?')
parser$add_argument("--noise_procs", default=c(""))
parser$add_argument("--dataset_ID", default="UCLA_Schizophrenia")
parser$add_argument("--run_number", nargs='?')

# Parse input arguments
args <- parser$parse_args()
github_dir <- args$github_dir
data_path <- args$data_path
univariate_feature_set <- args$univariate_feature_set
sample_metadata_file <- args$sample_metadata_file
brain_region_lookup <- args$brain_region_lookup
noise_procs <- args$noise_procs
dataset_ID <- args$dataset_ID
run_number <- args$run_number

# univariate_feature_set <- "catch22"
# github_dir <- "/headnode1/abry4213/github/"

# UCLA schizophrenia
# data_path <- "/headnode1/abry4213/data/UCLA_Schizophrenia/"
# dataset_ID <- "UCLA_Schizophrenia"
# sample_metadata_file <- "UCLA_Schizophrenia_sample_metadata.Rds"
# noise_procs <- c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER")
# brain_region_lookup <- "Brain_Region_info.csv"

# ABIDE ASD
# data_path <- "/headnode1/abry4213/data/ABIDE_ASD/"
# dataset_ID <- "ABIDE_ASD"
# sample_metadata_file <- "ABIDE_ASD_sample_metadata.Rds"
# noise_procs <- c("FC1000")
# brain_region_lookup <- "Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv"

# HCP100
# data_path <- "/headnode1/abry4213/data/HCP100/"
# dataset_ID <- "HCP100"
# sample_metadata_file <- "HCP100_sample_metadata.Rds"
# noise_procs <- c("AROMA+2P+GMR")
# brain_region_lookup <- "Brain_Region_info.csv"

if (!is.null(run_number)) {
  rdata_path <- paste0(data_path, "processed_data_run", run_number, "/Rdata/")
  plot_dir <- paste0(data_path, "plots_run", run_number, "/")
} else {
  rdata_path <- paste0(data_path, "processed_data/Rdata/")
  plot_dir <- paste0(data_path, "plots/")
}

icesTAF::mkdir(plot_dir)
icesTAF::mkdir(rdata_path)

# Set the seed
set.seed(127)

# Load tidyverse
library(tidyverse)
library(theft)

# Unlist noise-processing methods
tryCatch({
  noise_procs <- stringr::str_split(noise_procs, ";")[[1]]
  noise_procs <- unlist(noise_procs)
}, error = function(e) {
  message(e)
})


#-------------------------------------------------------------------------------
# Source helper scripts
#-------------------------------------------------------------------------------
# Set working directory to file location
helper_script_dir = paste0(github_dir, "fMRI_FeaturesDisorders/helper_functions/")
source(paste0(helper_script_dir, "data_prep_and_QC/TS_feature_extraction.R"))
source(paste0(helper_script_dir, "data_prep_and_QC/QC_functions_univariate.R"))


# DIY rlist::list.append
list.append <- function (.data, ...) 
{
  if (is.list(.data)) {
    c(.data, list(...))
  }
  else {
    c(.data, ..., recursive = FALSE)
  }
}


#-------------------------------------------------------------------------------
# Prepare data using dataset-specific script
#-------------------------------------------------------------------------------
system(sprintf("Rscript %s/data_prep_and_QC/dataset_specific_files/prepare_%s_time_series.R",
               github_dir, dataset_ID))

#-------------------------------------------------------------------------------
# Read in TS data for dataset, partitioned by noise-processing method
#-------------------------------------------------------------------------------

# Load brain region lookup table
brain_region_lookup_table <- read.csv(paste0(data_path, brain_region_lookup))
read_in_sample_TS_data <- function(sample_ID, noise_proc,
                                   brain_region_lookup_table) {
  noise_label <- gsub("\\+", "_", noise_proc)
  TS_data <- read.csv(paste0(data_path,
                             "raw_data/time_series_files/",
                             noise_label, "/",
                             sample_ID, "_TS.csv"),
                      header=F) %>%
    mutate(timepoint = 1:nrow(.)) %>%
    pivot_longer(cols = c(-timepoint),
                 names_to = "Index",
                 values_to = "values") %>%
    mutate(Sample_ID = sample_ID,
           Noise_Proc = noise_proc,
           Index = as.numeric(gsub("X|V", "", Index))) %>%
    left_join(., brain_region_lookup_table) %>%
    dplyr::select(Sample_ID, Noise_Proc, Brain_Region, timepoint, values)
}

if (!file.exists(paste0(data_path, "raw_data/", dataset_ID, "_fMRI_TS.Rds"))) {
  noise_proc_TS_data_list <- list()
  for (noise_proc in noise_procs) {
    noise_label <- gsub("\\+", "_", noise_proc)
    sample_IDs <- list.files(paste0(data_path, "raw_data/time_series_files/", noise_label)) %>%
      gsub("_TS.csv", "", .)
    np_TS_data <- sample_IDs %>%
      purrr::map_df(~ read_in_sample_TS_data(sample_ID = .x,
                                             noise_proc = noise_proc,
                                             brain_region_lookup_table = brain_region_lookup_table))
    noise_proc_TS_data_list <- list.append(noise_proc_TS_data_list, np_TS_data)
  }
  
  full_TS_data <- do.call(plyr::rbind.fill, noise_proc_TS_data_list)
  saveRDS(full_TS_data, paste0(data_path, "raw_data/",
                               dataset_ID, "_fMRI_TS.Rds"))
} else {
  full_TS_data <- readRDS(paste0(data_path, "raw_data/",
                                 dataset_ID, "_fMRI_TS.Rds"))
}


#-------------------------------------------------------------------------------
# Run catch22
#-------------------------------------------------------------------------------
catch22_all_samples(full_TS_data = full_TS_data,
                    rdata_path = rdata_path,
                    dataset_ID = dataset_ID,
                    unique_columns = c("Sample_ID", "Brain_Region", "Noise_Proc"),
                    output_column_names = c("Sample_ID", "Brain_Region", "Noise_Proc"))

#-------------------------------------------------------------------------------
# Perform QC for catch22 data
#-------------------------------------------------------------------------------
run_QC_for_univariate_dataset(data_path = data_path, 
                              proc_rdata_path = rdata_path,
                              sample_metadata_file = sample_metadata_file,
                              dataset_ID = dataset_ID,
                              univariate_feature_set = univariate_feature_set,
                              raw_TS_file = paste0(data_path, "raw_data/",
                                                   dataset_ID, "_fMRI_TS.Rds"),
                              noise_procs = noise_procs,
                              plot_dir = plot_dir)
