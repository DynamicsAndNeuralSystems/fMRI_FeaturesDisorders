#-------------------------------------------------------------------------------
# Define paths
#-------------------------------------------------------------------------------
# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--project_path", default="/project/hctsa/annie/")
parser$add_argument("--github_dir", default="/project/hctsa/annie/github/")
parser$add_argument("--data_path", default="/project/hctsa/annie/data/scz/UCLA/")
parser$add_argument("--rdata_path", default="/project/hctsa/annie/data/scz/UCLA/Rdata/")
parser$add_argument("--feature_set", default="catch22")
parser$add_argument("--input_mat_file", default="new/UCLA_time_series_four_groups.mat")
parser$add_argument("--subject_csv", default="participants.csv")
# project_path <- "D:/Virtual_Machines/Shared_Folder/github/"
# github_dir <- "D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders/"
# rdata_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/Rdata/"
# univariate_feature_set <- "catch22"
# input_mat_file = paste0(data_path, "new/UCLA_time_series_four_groups.mat")
# subject_csv <- paste0(data_path, "participants.csv")

# Parse input arguments
args <- parser$parse_args()
project_path <- args$project_path
github_dir <- args$github_dir
rdata_path <- args$rdata_path
univariate_feature_set <- args$univariate_feature_set
input_mat_file <- args$input_mat_file
subject_csv <- args$subject_csv

# Set the seed
set.seed(127)

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