#-------------------------------------------------------------------------------
# Define paths
#-------------------------------------------------------------------------------
# Parse arguments

library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--github_dir", default="~/github/")
parser$add_argument("--data_path", default="~/data/UCLA_CNP/")
parser$add_argument("--python_to_use", default="~/.conda/envs/pyspi/bin/python3")
parser$add_argument("--movement_file_path", default="~/data/UCLA_CNP/movement_data/movement_files/")
parser$add_argument("--sample_metadata_file", default="UCLA_CNP_sample_metadata.feather")
parser$add_argument("--dataset_ID", default="UCLA_CNP")

# Parse input arguments
args <- parser$parse_args()
github_dir <- args$github_dir
data_path <- args$data_path
python_to_use <- args$python_to_use
sample_metadata_file <- args$sample_metadata_file
movement_file_path <- args$movement_file_path
dataset_ID <- args$dataset_ID

# python_to_use <- "~/.conda/envs/pyspi/bin/python3"
# github_dir <- "~/github/"
# data_path <- "~/data/UCLA_CNP/"
# sample_metadata_file <- "UCLA_CNP_sample_metadata.feather"
# movement_file_path <- "~/data/UCLA_CNP/movement_data/movement_files/"
# dataset_ID <- "UCLA_CNP"

reticulate::use_python(python_to_use)
library(feather)
library(tidyverse)
library(reticulate)
library(glue)

movement_data <- read.table(glue("{movement_file_path}/{dataset_ID}_mFD.txt"), sep=",")
colnames(movement_data) <- c("Sample_ID", "FD_Jenkinson", "FD_Power", "FD_VanDijk")

