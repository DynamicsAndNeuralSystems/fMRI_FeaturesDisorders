#-------------------------------------------------------------------------------
# Define paths
#-------------------------------------------------------------------------------
# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--github_dir", default="~/github/")
parser$add_argument("--data_path", default="~/data/UCLA_CNP/")
parser$add_argument("--univariate_feature_set", default="catch24")
parser$add_argument("--sample_metadata_file", default="UCLA_CNP_sample_metadata.feather")
parser$add_argument("--brain_region_lookup", default="", nargs='?')
parser$add_argument("--noise_proc", default="")
parser$add_argument("--dataset_ID", default="UCLA_CNP")

# Parse input arguments
args <- parser$parse_args()
github_dir <- args$github_dir
data_path <- args$data_path
univariate_feature_set <- args$univariate_feature_set
sample_metadata_file <- args$sample_metadata_file
brain_region_lookup <- args$brain_region_lookup
noise_proc <- args$noise_proc
dataset_ID <- args$dataset_ID

# univariate_feature_set <- "catch25"
# github_dir <- "~/github/"
# github_dir <- "~/Library/CloudStorage/OneDrive-TheUniversityofSydney(Students)/github/"

# # UCLA CNP
# data_path <- "/headnode1/abry4213/data/UCLA_CNP/"
# dataset_ID <- "UCLA_CNP"
# sample_metadata_file <- "UCLA_CNP_sample_metadata.feather"
# noise_proc <- "AROMA+2P+GMR"
# brain_region_lookup <- "UCLA_CNP_Brain_Region_Lookup.feather"

# # ABIDE ASD
# data_path <- "~/data/ABIDE_ASD/"
# dataset_ID <- "ABIDE_ASD"
# sample_metadata_file <- "ABIDE_ASD_sample_metadata.feather"
# noise_proc <- "FC1000"
# brain_region_lookup <- "ABIDE_ASD_Harvard_Oxford_cort_prob_2mm_ROI_lookup.txt"

noise_label <- gsub("\\+", "_", noise_proc)
plot_dir <- paste0(data_path, "plots/")

TAF::mkdir(plot_dir)

# Set the seed
set.seed(127)

python_to_use <- "~/.conda/envs/pyspi/bin/python3"
# python_to_use <- "/Users/abry4213/anaconda3/envs/pyspi/bin/python3"
reticulate::use_python(python_to_use)

# Load tidyverse
library(tidyverse)
library(theft)
library(feather)
library(reticulate)

# Import pyarrow.feather as pyarrow_feather
pyarrow_feather <- import("pyarrow.feather")

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
# system(sprintf("Rscript %s/fMRI_FeaturesDisorders/prep_data_and_QC/dataset_specific_files/prepare_%s_time_series.R",
#                github_dir, dataset_ID))

#-------------------------------------------------------------------------------
# Read in TS data for dataset, partitioned by noise-processing method
#-------------------------------------------------------------------------------

# Load brain region lookup table
if (str_detect(brain_region_lookup, ".feather")) {
  brain_region_lookup_table <- pyarrow_feather$read_feather(paste0(data_path, "study_metadata/", brain_region_lookup))
} else {
  brain_region_lookup_table <- read.table(paste0(data_path, "study_metadata/", brain_region_lookup),
                                          sep = "\t", header = T)
}

read_in_sample_TS_data <- function(sample_ID, dataset_ID, noise_proc,
                                   brain_region_lookup_table) {
  noise_label <- gsub("\\+", "_", noise_proc)
  tryCatch({TS_data <- read.csv(paste0(data_path,
                                       "raw_data/time_series_files/",
                                       noise_label, "/",
                                       sample_ID, "_TS.csv"),
                                header=T) %>%
    mutate(timepoint = 1:nrow(.)) %>%
    pivot_longer(cols = c(-timepoint),
                 names_to = "Index",
                 values_to = "values") %>%
    mutate(Sample_ID = sample_ID,
           Noise_Proc = noise_proc,
           Index = as.numeric(gsub("X|V", "", Index))) %>%
    left_join(., brain_region_lookup_table) %>%
    dplyr::select(Sample_ID, Noise_Proc, Brain_Region, timepoint, values)
  return(TS_data)},
  error = function(e) {
    cat("\nError for subject", sample_ID, "\n")
    message(e)
  })
}

if (!file.exists(paste0(data_path, "raw_data/",
                        dataset_ID, "_", noise_label, "_fMRI_TS.feather"))) {
  noise_proc_TS_data_list <- list()
  noise_label <- gsub("\\+", "_", noise_proc)
  sample_IDs <- list.files(paste0(data_path, "raw_data/time_series_files/", noise_label)) %>%
    gsub("_TS.csv", "", .)
  np_TS_data <- sample_IDs %>%
    purrr::map_df(~ read_in_sample_TS_data(sample_ID = .x,
                                           dataset_ID = dataset_ID,
                                           noise_proc = noise_proc,
                                           brain_region_lookup_table = brain_region_lookup_table))
  
  if (nrow(np_TS_data) > 0) {
    feather::write_feather(np_TS_data, paste0(data_path, "raw_data/",
                                              dataset_ID, "_", noise_label, "_fMRI_TS.feather"))
  } else {
    stop("No time-series data available to save.")
  }
  
} else {
  np_TS_data <- pyarrow_feather$read_feather(paste0(data_path, "raw_data/",
                                           dataset_ID, "_", noise_label,  "_fMRI_TS.feather"))
}

#-------------------------------------------------------------------------------
# Run catch24
#-------------------------------------------------------------------------------
# catch24_all_samples(full_TS_data = np_TS_data,
#                     data_path = data_path,
#                     noise_proc = noise_proc,
#                     dataset_ID = dataset_ID,
#                     unique_columns = c("Sample_ID", "Brain_Region", "Noise_Proc"),
#                     output_column_names = c("Sample_ID", "Brain_Region", "Noise_Proc"))

#-------------------------------------------------------------------------------
# Compute fALFF to be 25th feature for all samples
#-------------------------------------------------------------------------------
# file_base <- glue("{data_path}/raw_data/{dataset_ID}_{noise_label}_fMRI_TS")

# # Convert mat file to feather
# system(sprintf("python3 %s/fMRI_FeaturesDisorders/prep_data_and_QC/feather_to_mat.py %s mat",
#                github_dir, file_base))
               
# # Run matlab fALFF script
# SP_fALFF_script_path <- glue("{github_dir}/fMRI_FeaturesDisorders/prep_data_and_QC/")
# TS_mat_file <- glue("{data_path}/raw_data/{dataset_ID}_{noise_label}_fMRI_TS.mat")
# output_mat_file <- glue("{data_path}/processed_data/{dataset_ID}_{noise_label}_fALFF.mat")
# system(sprintf("/usr/physics/matlab/matlab2019b/bin/matlab -nodisplay -singleCompThread -r \"compute_regional_fALFF %s %s %s %s; exit\"",
#                  data_path, SP_fALFF_script_path, TS_mat_file, output_mat_file))

# # Convert resulting mat file back to feather
# fALFF_file_base =  glue("{data_path}/processed_data/{dataset_ID}_{noise_label}_fALFF")
# samples_to_use = glue("{data_path}/processed_data/{dataset_ID}_filtered_sample_info_{noise_label}_final.feather")
# system(sprintf("python3 %s/fMRI_FeaturesDisorders/prep_data_and_QC/feather_to_mat.py %s feather %s %s %s",
#                github_dir, fALFF_file_base, dataset_ID, noise_proc, samples_to_use))


#-------------------------------------------------------------------------------
# Perform QC for catch25 data
#-------------------------------------------------------------------------------
univariate_feature_set <- "catch25"
run_QC_for_univariate_dataset(data_path = data_path, 
                              sample_metadata_file = sample_metadata_file,
                              dataset_ID = dataset_ID,
                              univariate_feature_set = univariate_feature_set,
                              raw_TS_file = paste0(data_path, "raw_data/",
                                                    dataset_ID, "_", noise_label, "_fMRI_TS.feather"),
                              theft_processed_data_file = glue("{dataset_ID}_{noise_label}_catch24.feather"),
                              fALFF_processed_data_file = glue("{dataset_ID}_{noise_label}_fALFF.feather"),
                              noise_proc = noise_proc)

