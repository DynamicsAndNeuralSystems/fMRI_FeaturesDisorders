#-------------------------------------------------------------------------------
# Define paths
#-------------------------------------------------------------------------------
# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--project_path", default="/project/hctsa/annie/")
parser$add_argument("--github_dir", default="/project/hctsa/annie/github/fMRI_FeaturesDisorders/")
parser$add_argument("--data_path", default="/project/hctsa/annie/data/scz/UCLA/")
parser$add_argument("--rdata_path", default="/project/hctsa/annie/data/scz/UCLA/Rdata/")
parser$add_argument("--univariate_feature_set", default="catch22")
parser$add_argument("--input_mat_file", default="new/UCLA_time_series_four_groups.mat")
parser$add_argument("--subject_csv", default="participants.csv")
parser$add_argument("--noise_procs", default=c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER"))
parser$add_argument("--dataset_ID", default="UCLA")
parser$add_argument("--plot_dir", default="/project/hctsa/annie/github/fMRI_FeaturesDisorders/plots/")
# project_path <- "D:/Virtual_Machines/Shared_Folder/github/"
# github_dir <- "D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders/"
# rdata_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/Rdata/"
# data_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/"
# univariate_feature_set <- "catch22"
# dataset_ID <- "UCLA"
# input_mat_file = "new/UCLA_time_series_four_groups.mat"
# subject_csv <- "participants.csv"
# noise_procs <- c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER")
# plot_dir <- paste0(github_dir, "plots/")

# Parse input arguments
args <- parser$parse_args()
project_path <- args$project_path
github_dir <- args$github_dir
data_path <- args$data_path
rdata_path <- args$rdata_path
univariate_feature_set <- args$univariate_feature_set
input_mat_file <- args$input_mat_file
subject_csv <- args$subject_csv
noise_procs <- args$noise_procs
dataset_ID <- args$dataset_ID
plot_dir <- args$plot_dir

# Set the seed
set.seed(127)

# Load tidyverse
library(tidyverse)

#-------------------------------------------------------------------------------
# Source helper scripts
#-------------------------------------------------------------------------------
# Set working directory to file location
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
helper_script_dir = "helper_functions/"
files.sources = list.files(helper_script_dir, pattern=".R", full.names = T) %>% .[!(str_detect(., "cluster"))]
sapply(files.sources, source)

#-------------------------------------------------------------------------------
# Prep univariate data
#-------------------------------------------------------------------------------
load_mat_data(mat_file=paste0(data_path, input_mat_file), 
              subject_csv=paste0(data_path, subject_csv), 
              rdata_path=rdata_path, 
              overwrite=TRUE)

#-------------------------------------------------------------------------------
# Run catch22
#-------------------------------------------------------------------------------
catch22_all_samples(TS_data_file = paste0(rdata_path, sprintf("%s_fMRI_TimeSeries.Rds",
                                                              dataset_ID)), 
                    rdata_path,
                    input_dataset = dataset_ID,
                    output_column_names = c("Sample_ID", "Brain_Region", "Noise_Proc"),
                    noise_procs = noise_procs)

#-------------------------------------------------------------------------------
# Perform QC for catch22 data
#-------------------------------------------------------------------------------
rmarkdown::render(input = paste0(helper_script_dir, "QC_report_template.Rmd"),
                  output_file = paste0(plot_dir, "QC_report_catch22.html"),
                  params = list(rdata_path = rdata_path,
                                dataset_ID = dataset_ID,
                                univariate_feature_set = univariate_feature_set,
                                noise_procs = noise_procs,
                                raw_TS_file = paste0(rdata_path, "UCLA_fMRI_TimeSeries.Rds")))

#-------------------------------------------------------------------------------
# Prep univariate data
#-------------------------------------------------------------------------------
load_mat_data(mat_file=paste0(data_path, input_mat_file), 
              subject_csv=paste0(data_path, subject_csv), 
              rdata_path=rdata_path, 
              overwrite=TRUE)