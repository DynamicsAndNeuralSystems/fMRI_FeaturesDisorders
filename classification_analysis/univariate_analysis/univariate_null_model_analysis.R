# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--github_dir", default="/headnode1/abry4213/github/")
parser$add_argument("--data_path", default="/headnode1/abry4213/data/UCLA_Schizophrenia/")
parser$add_argument("--sample_metadata_file", default="UCLA_Schizophrenia_sample_metadata.Rds")
parser$add_argument("--univariate_feature_set", default="catch22")
parser$add_argument("--noise_proc", default="AROMA+2P+GMR")
parser$add_argument("--dataset_ID", default="UCLA_Schizophrenia")
parser$add_argument("--run_number", nargs="?")

# Parse input arguments
args <- parser$parse_args()
github_dir <- args$github_dir
sample_metadata_file <- args$sample_metadata_file
data_path <- args$data_path
noise_proc <- args$noise_proc
univariate_feature_set <- args$univariate_feature_set
dataset_ID <- args$dataset_ID
run_number <- args$run_number

# univariate_feature_set <- "catch22"
# pairwise_feature_set <- "pyspi14"
# github_dir <- "/headnode1/abry4213/github/"

# UCLA schizophrenia
# data_path <- "/headnode1/abry4213/data/UCLA_Schizophrenia/"
# dataset_ID <- "UCLA_Schizophrenia"
# sample_metadata_file <- "UCLA_Schizophrenia_sample_metadata.Rds"
# noise_proc <- "AROMA+2P+GMR"
# run_number <- 2

# ABIDE ASD
# data_path <- "/headnode1/abry4213/data/ABIDE_ASD/"
# sample_metadata_file <- "ABIDE_ASD_sample_metadata.Rds"
# dataset_ID <- "ABIDE_ASD"
# noise_proc <- "FC1000"

if (!is.null(run_number)) {
  rdata_path <- paste0(data_path, "processed_data_run", run_number, "/Rdata/")
  plot_dir <- paste0(data_path, "plots_run", run_number, "/")
} else {
  rdata_path <- paste0(data_path, "processed_data/Rdata/")
  plot_dir <- paste0(data_path, "plots/")
}


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
      cat("Now merging null model fits data\n")
      # Define data directory
      output_data_dir <- paste0(rdata_path, sprintf("%s_%s_wise_%s_%s_null_model_fits/",
                                                    dataset_ID,
                                                    grouping_type, 
                                                    univariate_feature_set, 
                                                    weighting_name))
      
      
      ## Concatenate null results and save to RDS file
      model_permutation_null_weighting <- list.files(output_data_dir, pattern="Rds") %>%
        purrr::map_df(~ readRDS(paste0(output_data_dir, .x)))
      cat("Now saving null model fits data\n")
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
      cat("Now calculating p-values\n")
      group_wise_SVM_balanced_accuracy_across_repeats <- readRDS(paste0(rdata_path, sprintf("%s_wise_CV_linear_SVM_%s_%s_balacc_across_repeats.Rds",
                                                                                            grouping_type, 
                                                                                            univariate_feature_set, 
                                                                                            weighting_name)))
      
      # Calculate p-values
      pvalues <- calc_empirical_nulls(class_res = group_wise_SVM_balanced_accuracy_across_repeats,
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
