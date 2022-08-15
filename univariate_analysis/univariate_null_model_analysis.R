# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--github_dir", default="/project/hctsa/annie/github/")
parser$add_argument("--rdata_path", default="/project/hctsa/annie/data/scz/UCLA/Rdata/")
parser$add_argument("--univariate_feature_set", default="catch22")
parser$add_argument("--noise_proc", default="AROMA+2P+GMR")
parser$add_argument("--dataset_ID", default="UCLA_Schizophrenia")

# Parse input arguments
args <- parser$parse_args()
github_dir <- args$github_dir
rdata_path <- args$rdata_path
noise_proc <- args$noise_proc
univariate_feature_set <- args$univariate_feature_set
dataset_ID <- args$dataset_ID

# univariate_feature_set <- "catch22"
# pairwise_feature_set <- "pyspi_19"
# subject_csv <- "participants.csv"
# github_dir <- "/headnode1/abry4213/github/fMRI_FeaturesDisorders/"

# UCLA schizophrenia
# rdata_path <- "/headnode1/abry4213/data/UCLA_Schizophrenia/Rdata/"
# data_path <- "/headnode1/abry4213/data/UCLA_Schizophrenia/"
# dataset_ID <- "UCLA_Schizophrenia"
# noise_proc <- "AROMA+2P+GMR"

# ABIDE ASD
# rdata_path <- "/headnode1/abry4213/data/ABIDE_ASD/Rdata/"
# data_path <- "/headnode1/abry4213/data/ABIDE_ASD/"
# dataset_ID <- "ABIDE_ASD"
# noise_proc <- "FC1000"


#-------------------------------------------------------------------------------
# Source helper scripts
#-------------------------------------------------------------------------------
# Set working directory to file location
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
helper_script_dir = "../helper_functions/classification/"
source(paste0(helper_script_dir, "Linear_SVM.R"))
source(paste0(helper_script_dir, "Null_distributions.R"))

# Set the seed
set.seed(127)

# Load tidyverse
library(tidyverse)

# Unlist noise-processing methods
tryCatch({
  noise_procs <- unlist(noise_procs)
}, error = function(e) {})

# Use e1071 SVM with a linear kernel
kernel = "linear"

################################################################################
# Define weighting parameters
################################################################################

grouping_param_df <- data.frame(grouping_type = c("ROI", "Feature", "Combo"),
                                grouping_var = c("Brain_Region", "Feature", "Combo"),
                                SVM_feature_var = c("Feature", "Brain_Region", "Combo")) 

weighting_param_df <- data.frame(name = c("inv_prob"),
                                 use_inv_prob_weighting = c(TRUE))

for (i in 1:nrow(grouping_param_df)) {
  grouping_type = grouping_param_df$grouping_type[i]
  grouping_var = grouping_param_df$grouping_var[i]
  SVM_feature_var = grouping_param_df$SVM_feature_var[i]
  
  #### 10-fold linear SVM with different weights
  # Iterate over weighting_param_df 
  for (j in 1:nrow(weighting_param_df)) {
    weighting_name <- weighting_param_df$name[j]
    use_inv_prob_weighting <- weighting_param_df$use_inv_prob_weighting[j]
    
    # Generate null-model fits distribution
    if (!file.exists(paste0(rdata_path, sprintf("%s_%s_wise_model_permutation_null_%s_%s.Rds",
                                                dataset_ID,
                                                grouping_type,
                                                univariate_feature_set,
                                                weighting_name)))) {
      # Define data directory
      output_data_dir <- paste0(rdata_path, sprintf("%s_%s_wise_%s_%s_null_model_fits/",
                                                    dataset_ID,
                                                    grouping_type, 
                                                    univariate_feature_set, 
                                                    weighting_name))
      
      
      ## Concatenate null results and save to RDS file
      model_permutation_null_weighting <- list.files(output_data_dir, pattern="Rds") %>%
        purrr::map_df(~ readRDS(paste0(output_data_dir, .x)))
      saveRDS(model_permutation_null_weighting, paste0(rdata_path, sprintf("%s_%s_wise_model_permutation_null_%s_%s.Rds",
                                                                           dataset_ID,
                                                                           grouping_type,
                                                                           univariate_feature_set,
                                                                           weighting_name)))
    } else {
      model_permutation_null_weighting <- readRDS(paste0(rdata_path, sprintf("%s_%s_wise_model_permutation_null_%s_%s.Rds",
                                                                             dataset_ID,
                                                                             grouping_type,
                                                                             univariate_feature_set,
                                                                             weighting_name)))
    }
    
    # Empirically derive p-values based on null model fits distribution
    if (!file.exists(paste0(rdata_path, sprintf("%s_wise_CV_linear_SVM_model_permutation_null_%s_%s_pvals.Rds",
                                                grouping_type,
                                                univariate_feature_set,
                                                weighting_name)))) {
      group_wise_SVM_balanced_accuracy <- readRDS(paste0(rdata_path, sprintf("%s_wise_CV_linear_SVM_%s_%s_balacc.Rds",
                                                                             grouping_type, 
                                                                             univariate_feature_set, 
                                                                             weighting_name)))
      
      # Calculate p-values
      pvalues <- calc_empirical_nulls(class_res = group_wise_SVM_balanced_accuracy,
                                      null_data = model_permutation_null_weighting,
                                      feature_set = univariate_feature_set,
                                      noise_proc = noise_proc,
                                      use_pooled_null = TRUE,
                                      is_main_data_averaged = TRUE,
                                      grouping_var = grouping_var)
      
      saveRDS(pvalues, file=paste0(rdata_path, sprintf("%s_wise_CV_linear_SVM_model_permutation_null_%s_%s_pvals.Rds",
                                                       grouping_type,
                                                       univariate_feature_set,
                                                       weighting_name)))
    }
  }
}
