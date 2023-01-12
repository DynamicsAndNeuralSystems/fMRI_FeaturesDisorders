#-------------------------------------------------------------------------------
# Define paths
#-------------------------------------------------------------------------------
# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--github_dir", default="~/github/")
parser$add_argument("--data_path", default="~/data/UCLA_CNP/")
parser$add_argument("--univariate_feature_set", default="catch22")
parser$add_argument("--sample_metadata_file", default="UCLA_CNP_sample_metadata.Rds")
parser$add_argument("--brain_region_lookup", default="", nargs='?')
parser$add_argument("--noise_proc", default="")
parser$add_argument("--dataset_ID", default="UCLA_CNP")
parser$add_argument("--add_catch2", action="store_true", default=FALSE)

# Parse input arguments
args <- parser$parse_args()
github_dir <- args$github_dir
data_path <- args$data_path
univariate_feature_set <- args$univariate_feature_set
sample_metadata_file <- args$sample_metadata_file
brain_region_lookup <- args$brain_region_lookup
noise_proc <- args$noise_proc
dataset_ID <- args$dataset_ID
add_catch2 <- args$add_catch2

# univariate_feature_set <- "catch22"
# github_dir <- "~/github/"
# add_catch2 <- TRUE

# # UCLA CNP
# data_path <- "~/data/UCLA_CNP_ABIDE_ASD/"
# dataset_ID <- "UCLA_CNP"
# sample_metadata_file <- "UCLA_CNP_sample_metadata.Rds"
# noise_proc <- "AROMA+2P+GMR"
# brain_region_lookup <- "UCLA_CNP_Brain_Region_info.csv"

# # ABIDE ASD
# data_path <- "~/data/UCLA_CNP_ABIDE_ASD/"
# dataset_ID <- "ABIDE_ASD"
# sample_metadata_file <- "ABIDE_ASD_sample_metadata.Rds"
# noise_proc <- "FC1000"
# brain_region_lookup <- "ABIDE_ASD_Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv"

noise_label <- gsub("\\+", "_", noise_proc)
rdata_path <- paste0(data_path, "processed_data/Rdata/")
plot_dir <- paste0(data_path, "plots/")

TAF::mkdir(plot_dir)
TAF::mkdir(rdata_path)

# Set the seed
set.seed(127)

# Load tidyverse
library(tidyverse)
library(theft)

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
system(sprintf("Rscript %s/fMRI_FeaturesDisorders/prep_data_and_QC/dataset_specific_files/prepare_%s_time_series.R",
               github_dir, dataset_ID))

#-------------------------------------------------------------------------------
# Read in TS data for dataset, partitioned by noise-processing method
#-------------------------------------------------------------------------------

# Load brain region lookup table
brain_region_lookup_table <- read.csv(paste0(data_path, "study_metadata/", brain_region_lookup))
read_in_sample_TS_data <- function(sample_ID, dataset_ID, noise_proc,
                                   brain_region_lookup_table) {
  noise_label <- gsub("\\+", "_", noise_proc)
  TS_data <- read.csv(paste0(data_path,
                             "raw_data/", dataset_ID, "/time_series_files/",
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

if (!file.exists(paste0(data_path, "raw_data/",
                        dataset_ID, "_", noise_label, "_fMRI_TS.Rds"))) {
  noise_proc_TS_data_list <- list()
  noise_label <- gsub("\\+", "_", noise_proc)
  sample_IDs <- list.files(paste0(data_path, "raw_data/", dataset_ID, "/time_series_files/", noise_label)) %>%
    gsub("_TS.csv", "", .)
  np_TS_data <- sample_IDs %>%
    purrr::map_df(~ read_in_sample_TS_data(sample_ID = .x,
                                           dataset_ID = dataset_ID,
                                            noise_proc = noise_proc,
                                            brain_region_lookup_table = brain_region_lookup_table))

  saveRDS(np_TS_data, paste0(data_path, "raw_data/",
                               dataset_ID, "_", noise_label, "_fMRI_TS.Rds"))
} else {
  np_TS_data <- readRDS(paste0(data_path, "raw_data/",
                                 dataset_ID, "_", noise_label,  "_fMRI_TS.Rds"))
}


#-------------------------------------------------------------------------------
# Run catch22
#-------------------------------------------------------------------------------
catch22_all_samples(full_TS_data = np_TS_data,
                    rdata_path = rdata_path,
                    noise_proc = noise_proc,
                    dataset_ID = dataset_ID,
                    unique_columns = c("Sample_ID", "Brain_Region", "Noise_Proc"),
                    output_column_names = c("Sample_ID", "Brain_Region", "Noise_Proc"),
                    add_mean_SD = add_catch2,
                    overwrite = FALSE)

#-------------------------------------------------------------------------------
# Perform QC for catch22 data
#-------------------------------------------------------------------------------
run_QC_for_univariate_dataset(data_path = data_path, 
                              proc_rdata_path = rdata_path,
                              sample_metadata_file = sample_metadata_file,
                              dataset_ID = dataset_ID,
                              univariate_feature_set = univariate_feature_set,
                              add_catch2 = add_catch2,
                              raw_TS_file = paste0(data_path, "raw_data/",
                                                   dataset_ID, "_", noise_label, "_fMRI_TS.Rds"),
                              noise_proc = noise_proc,
                              plot_dir = plot_dir)

# save filtered catch22 data to a CSV
catch22_filtered_zscored <- readRDS(paste0(rdata_path, sprintf("%s_%s_%s_filtered_zscored.Rds",
                                                          dataset_ID,
                                                          noise_label,
                                                          univariate_feature_set)))
                                                          
write_feather(catch22_filtered_zscored, paste0(rdata_path, sprintf("%s_%s_%s_filtered_zscored.feather",
                                                                   dataset_ID,
                                                                   noise_label,
                                                                   univariate_feature_set)))                                                     