library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--project_path", default="/project/hctsa/annie/")
parser$add_argument("--github_dir", default="/project/hctsa/annie/github/")
parser$add_argument("--dataset_ID", default="UCLA_Schizophrenia")
parser$add_argument("--data_path", default="/project/hctsa/annie/data/UCLA_Schizophrenia/")
parser$add_argument("--rdata_path", default="/project/hctsa/annie/data/UCLA_Schizophrenia/Rdata/")
parser$add_argument("--plot_dir", default="/project/hctsa/annie/data/UCLA_Schizophrenia/plots/Misclassification_Analysis/")
parser$add_argument("--univariate_feature_set", default="catch22")

# Parse input arguments
args <- parser$parse_args()
project_path <- args$project_path
dataset_ID <- args$dataset_ID
github_dir <- args$github_dir
data_path <- args$data_path
rdata_path <- args$rdata_path
univariate_feature_set <- args$univariate_feature_set

# project_path <- "D:/Virtual_Machines/Shared_Folder/github/"
# dataset_ID <- "UCLA_Schizophrenia"
# github_dir <- "D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders/"
# data_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/UCLA_Schizophrenia/"
# plot_dir <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/UCLA_Schizophrenia/plots/Misclassification_Analysis/"
# rdata_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/UCLA_Schizophrenia/Rdata/"
# univariate_feature_set <- "catch22"

# Load functions
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
set.seed(127)