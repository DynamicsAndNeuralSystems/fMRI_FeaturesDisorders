#-------------------------------------------------------------------------------
# Define paths
#-------------------------------------------------------------------------------
# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--project_path", default="/project/hctsa/annie/")
parser$add_argument("--github_dir", default="/project/hctsa/annie/github/fMRI_FeaturesDisorders/")
parser$add_argument("--data_path", default="/project/hctsa/annie/data/UCLA_Schizophrenia/")
parser$add_argument("--rdata_path", default="/project/hctsa/annie/data/UCLA_Schizophrenia/Rdata/")
parser$add_argument("--univariate_feature_set", default="catch22")
parser$add_argument("--parcellation_name", default="aparc+aseg", nargs='?')
parser$add_argument("--input_mat_file", default="")
parser$add_argument("--subject_csv", default="participants.csv")
parser$add_argument("--brain_region_lookup", default="", nargs='?')
parser$add_argument("--noise_procs", default=c(""), nargs="*", action="append")
parser$add_argument("--dataset_ID", default="UCLA_Schizophrenia")

# univariate_feature_set <- "catch22"
# subject_csv <- "participants.csv"
# project_path <- "D:/Virtual_Machines/Shared_Folder/github/"
# github_dir <- "D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders/"

# UCLA schizophrenia
# rdata_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/UCLA_Schizophrenia/Rdata/"
# data_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/UCLA_Schizophrenia/"
# dataset_ID <- "UCLA_Schizophrenia"
# input_mat_file = "new/UCLA_time_series_four_groups.mat"
# noise_procs <- c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER")
# parcellation_name <- "aparc+aseg"
# brain_region_lookup <- ""

# ABIDE ASD
# rdata_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/ABIDE_ASD/Rdata/"
# data_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/ABIDE_ASD/"
# dataset_ID <- "ABIDE_ASD"
# noise_procs <- c("FC1000")
# input_mat_file = "NA"
# parcellation_name <- "harvard_oxford_cort_prob_2mm"
# brain_region_lookup <- "Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv"

# Parse input arguments
args <- parser$parse_args()
project_path <- args$project_path
github_dir <- args$github_dir
data_path <- args$data_path
rdata_path <- args$rdata_path
univariate_feature_set <- args$univariate_feature_set
input_mat_file <- args$input_mat_file
subject_csv <- args$subject_csv
parcellation_name <- args$parcellation_name
noise_procs <- args$noise_procs
dataset_ID <- args$dataset_ID
plot_dir <- args$plot_dir
brain_region_lookup <- args$brain_region_lookup

plot_dir <- paste0(data_path, "plots/")
icesTAF::mkdir(plot_dir)

# Set the seed
set.seed(127)

# Load tidyverse
library(tidyverse)

# Unlist noise-processing methods
tryCatch({
  noise_procs <- unlist(noise_procs)
}, error = function(e) {})


#-------------------------------------------------------------------------------
# Source helper scripts
#-------------------------------------------------------------------------------
# Set working directory to file location
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
helper_script_dir = "helper_functions/"
files.sources = list.files(helper_script_dir, pattern=".R$", full.names = T) %>% .[!(str_detect(., "cluster"))]
invisible(sapply(files.sources, source))

#-------------------------------------------------------------------------------
# Prepare data using dataset-specific script
#-------------------------------------------------------------------------------
system(sprintf("Rscript %s/data_prep_and_QC/prepare_%s_univariate_data.R --data_path %s --input_mat_file %s --subject_csv %s --noise_procs %s --dataset_ID %s --github_dir %s --parcellation_name %s --brain_region_lookup %s",
               github_dir, dataset_ID, data_path, input_mat_file,
               subject_csv, paste(noise_procs, collapse=" "), dataset_ID, github_dir,
               parcellation_name, brain_region_lookup))

#-------------------------------------------------------------------------------
# Run catch22
#-------------------------------------------------------------------------------
catch22_all_samples(TS_data_file = paste0(rdata_path, sprintf("%s_fMRI_data.Rds",
                                                              dataset_ID)),
                    rdata_path,
                    input_dataset = dataset_ID,
                    unique_columns = c("Sample_ID", "Brain_Region", "Noise_Proc"),
                    output_column_names = c("Sample_ID", "Brain_Region", "Noise_Proc"))

#-------------------------------------------------------------------------------
# Perform QC for catch22 data
#-------------------------------------------------------------------------------
run_QC_for_univariate_dataset(rdata_path = rdata_path, 
                              dataset_ID = dataset_ID,
                              univariate_feature_set = univariate_feature_set,
                              raw_TS_file = paste0(rdata_path, dataset_ID, "_fMRI_data.Rds"),
                              noise_procs = noise_procs,
                              plot_dir = plot_dir)