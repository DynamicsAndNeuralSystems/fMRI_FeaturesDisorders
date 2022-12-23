# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--github_dir", default="/headnode1/abry4213/github/")
parser$add_argument("--data_path", default="/headnode1/abry4213/data/UCLA_Schizophrenia/")
parser$add_argument("--sample_metadata_file", default="UCLA_Schizophrenia_sample_metadata.Rds")
parser$add_argument("--feature_set", default="catch22")
parser$add_argument("--dataset_ID", default="UCLA_Schizophrenia")

# Parse input arguments
args <- parser$parse_args()
github_dir <- args$github_dir
sample_metadata_file <- args$sample_metadata_file
data_path <- args$data_path
feature_set <- args$feature_set
dataset_ID <- args$dataset_ID

# feature_set <- "catch22"
# github_dir <- "~/github/"
# data_path <- "~/data/UCLA_CNP_ABIDE_ASD/"
# dataset_ID <- "UCLA_CNP_ABIDE_ASD"
# sample_metadata_file <- "UCLA_CNP_ABIDE_ASD_sample_metadata.Rds"

rdata_path <- paste0(data_path, "processed_data/Rdata/")
plot_dir <- paste0(data_path, "plots/")

#-------------------------------------------------------------------------------
# Source helper scripts
#-------------------------------------------------------------------------------
helper_script_dir = paste0(github_dir, "fMRI_FeaturesDisorders/helper_functions/classification/")
source(paste0(helper_script_dir, "Linear_SVM.R"))
source(paste0(helper_script_dir, "Null_distributions.R"))

# Set the seed
set.seed(127)

# Load tidyverse
library(tidyverse)

# Use e1071 SVM with a linear kernel
kernel = "linear"


################################################################################
# Read in 10-fold 10-repeat CV SVM results
################################################################################

SVM_balanced_accuracy_across_repeats <- readRDS(paste0(rdata_path, sprintf("%s_univariate_CV_linear_SVM_%s_and_catch2_balanced_accuracy.Rds",
                                                                           dataset_ID,
                                                                           feature_set))) %>%
  filter(univariate_feature_set == feature_set)

################################################################################
# Read in null model results and calculat p-values
################################################################################

weighting_name <- "inv_prob"
grouping_types <- c("ROI", "Feature", "All")
grouping_vars <- c("Brain_Region", "Feature", "Combo")

# Function to read in model permutatino null data for files in directory
read_in_files_from_dir <- function(grouping_type, dataset_ID, univariate_feature_set, weighting_name) {
  # Define data directory
  output_data_dir <- paste0(rdata_path, sprintf("%s_%s_wise_%s_%s_null_model_fits/",
                                                dataset_ID,
                                                grouping_type, 
                                                feature_set, 
                                                weighting_name))
  
  model_permutation_null_weighting <- list.files(output_data_dir, pattern="Rds") %>%
    purrr::map_df(~ readRDS(paste0(output_data_dir, .x))) %>%
    mutate(Analysis_Type = grouping_type)
  
  return(model_permutation_null_weighting)
}

# Generate null-model fits distribution
if (!file.exists(paste0(rdata_path, sprintf("%s_CV_SVM_model_permutation_null_%s_%s.Rds",
                                            dataset_ID,
                                            feature_set,
                                            weighting_name)))) {
  
  model_permutation_null_weighting <- grouping_type %>%
    purrr::map_df(~ read_in_files_from_dir(grouping_type = .x,
                                           dataset_ID = dataset_ID,
                                           univariate_feature_set = feature_set,
                                           weighting_name = weighting_name))
  
  saveRDS(model_permutation_null_weighting, paste0(rdata_path, sprintf("%s_CV_SVM_model_permutation_null_%s_%s.Rds",
                                                                       dataset_ID,
                                                                       feature_set,
                                                                       weighting_name)))
} else {
  model_permutation_null_weighting <- readRDS(paste0(rdata_path, sprintf("%s_CV_SVM_model_permutation_null_%s_%s.Rds",
                                                                         dataset_ID,
                                                                         feature_set,
                                                                         weighting_name)))
}

# Empirically derive p-values based on null model fits distribution
if (!file.exists(paste0(rdata_path, sprintf("%s_CV_SVM_permutation_null_%s_%s_pvals.Rds",
                                            dataset_ID,
                                            feature_set,
                                            weighting_name)))) {
  cat("Now calculating p-values\n")
  pvalues <- calc_empirical_nulls(class_res = subset(SVM_balanced_accuracy_across_repeats, Analysis_Type=="ROI"),
                                  null_data = model_permutation_null_weighting,
                                  feature_set = feature_set,
                                  use_pooled_null = TRUE,
                                  is_main_data_averaged = TRUE,
                                  grouping_var = grouping_var)
  
  saveRDS(pvalues, file=paste0(rdata_path, sprintf("%s_CV_SVM_permutation_null_%s_%s_pvals.Rds",
                                                   dataset_ID,
                                                   feature_set,
                                                   weighting_name)))
  
}

