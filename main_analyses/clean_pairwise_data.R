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


# Load libraries
library()