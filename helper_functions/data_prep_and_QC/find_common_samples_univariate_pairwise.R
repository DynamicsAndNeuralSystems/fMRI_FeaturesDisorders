# Command-line arguments to parse
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")
parser$add_argument("--data_path", default="/headnode1/abry4213/data/UCLA_Schizophrenia/")
parser$add_argument("--univariate_feature_set", default="catch22")
parser$add_argument("--pairwise_feature_set", default="pyspi14")
parser$add_argument("--dataset_ID", default="UCLA_Schizophrenia")

# Parse input arguments
args <- parser$parse_args()
data_path <- args$data_path
univariate_feature_set <- args$univariate_feature_set
pairwise_feature_set <- args$pairwise_feature_set
dataset_ID <- args$dataset_ID

library(tidyverse)

rdata_path <- paste0(data_path, "processed_data/Rdata/")

# Read in univariate data
univariate_samples <- readRDS(paste0(rdata_path, 
                                     dataset_ID,
                                     "_filtered_sample_info_",
                                     univariate_feature_set,
                                     ".Rds")) %>%
  dplyr::select(Sample_ID)

# Read in pairwise data CSV
pairwise_samples <- readRDS(paste0(rdata_path,
                                   dataset_ID,
                                   "_filtered_sample_info_",
                                   pairwise_feature_set,
                                   ".Rds"))%>%
  dplyr::select(Sample_ID)

# Find the intersection
intersection <- inner_join(univariate_samples, pairwise_samples)

# Write the intersection results to an Rds file
saveRDS(intersection,
          paste0(rdata_path, dataset_ID, "_samples_with_univariate_",
                 univariate_feature_set,
                 "_and_pairwise_", pairwise_feature_set,
                 "_filtered.Rds"))