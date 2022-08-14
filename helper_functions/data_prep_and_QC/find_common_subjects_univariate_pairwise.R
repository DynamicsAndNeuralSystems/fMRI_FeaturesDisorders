# Command-line arguments to parse
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")
parser$add_argument("--data_path", default="/project/hctsa/annie/data/UCLA_Schizophrenia/")
parser$add_argument("--univariate_feature_set", default="catch22")
parser$add_argument("--pairwise_feature_set", default="pyspi_19")
parser$add_argument("--noise_procs", default=c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER"), nargs='*', action='append')
parser$add_argument("--dataset_ID", default="UCLA_Schizophrenia")

# Parse input arguments
args <- parser$parse_args()
data_path <- args$data_path
brain_region_lookup <- args$brain_region_lookup
univariate_feature_set <- args$univariate_feature_set
pairwise_feature_set <- args$pairwise_feature_set
dataset_ID <- args$dataset_ID
noise_procs <- args$noise_procs
overwrite <- args$overwrite

pydata_path <- paste0(data_path, "pydata/")
rdata_path <- paste0(data_path, "Rdata/")
set.seed(127)

library(tidyverse)

# Load info on univariate (theft) subjects
univariate_subject_info <- readRDS(paste0(rdata_path, "UCLA_filtered_subject_info_catch22.Rds"))

# 
pairwise_subject_info <- readRDS(paste0(rdata_path, "Filtered_subject_info_pyspi_19.Rds"))

intersection <- inner_join(pairwise_subject_info, univariate_subject_info)
saveRDS(intersection, file=paste0(rdata_path, "UCLA_Subjects_with_Univariate_and_Pairwise.Rds"))
