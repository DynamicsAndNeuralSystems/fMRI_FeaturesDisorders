#-------------------------------------------------------------------------------
# Define paths
#-------------------------------------------------------------------------------
# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--github_dir", default="/headnode1/abry4213/github/fMRI_FeaturesDisorders/")
parser$add_argument("--data_path", default="/headnode1/abry4213/data/UCLA_Schizophrenia/")
parser$add_argument("--rdata_path", default="/headnode1/abry4213/data/UCLA_Schizophrenia/Rdata/")
parser$add_argument("--pairwise_feature_set", default="pyspi_19")
parser$add_argument("--univariate_feature_set", default="catch22")
parser$add_argument("--noise_procs", default=c(""), nargs="*", action="append")
parser$add_argument("--dataset_ID", default="UCLA_Schizophrenia")

# univariate_feature_set <- "catch22"
# pairwise_feature_set <- "pyspi_19"
# subject_csv <- "participants.csv"
# github_dir <- "/headnode1/abry4213/github/fMRI_FeaturesDisorders/"

# UCLA schizophrenia
# rdata_path <- "/headnode1/abry4213/data/UCLA_Schizophrenia/Rdata/"
# data_path <- "/headnode1/abry4213/data/UCLA_Schizophrenia/"
# dataset_ID <- "UCLA_Schizophrenia"
# noise_procs <- c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER")

# ABIDE ASD
# rdata_path <- "/headnode1/abry4213/data/ABIDE_ASD/Rdata/"
# data_path <- "/headnode1/abry4213/data/ABIDE_ASD/"
# dataset_ID <- "ABIDE_ASD"
# noise_procs <- c("FC1000")

# Parse input arguments
args <- parser$parse_args()
github_dir <- args$github_dir
data_path <- args$data_path
rdata_path <- args$rdata_path
univariate_feature_set <- args$univariate_feature_set
noise_procs <- args$noise_procs
dataset_ID <- args$dataset_ID

plot_dir <- paste0(data_path, "plots/")
icesTAF::mkdir(plot_dir)

# Set the seed
set.seed(127)

# Load tidyverse
library(tidyverse)

# Unlist noise-processing methods
tryCatch({
  noise_procs <- unlist(noise_procs)
}, error = function(e) {})


#-------------------------------------------------------------------------------
# Source helper scripts
#-------------------------------------------------------------------------------
# Set working directory to file location
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
helper_script_dir = "../helper_functions/"
source(paste0(helper_script_dir, "Linear_SVM.R"))
source(paste0(helper_script_dir, "Null_distributions.R"))

################################################################################
# Generate model-free shuffle null distribution
################################################################################
# if (!file.exists(paste0(rdata_path, sprintf("%s_Null_Model_Free_Shuffles_%s.Rds",
#                                             dataset_ID, univariate_feature_set)))) {
#   model_free_shuffle_null_res <- run_model_free_n_shuffles(num_shuffles = 100000,
#                                                            feature_set = univariate_feature_set,
#                                                            dataset_ID = dataset_ID,
#                                                            rdata_path = rdata_path)
#   saveRDS(model_free_shuffle_null_res, file = paste0(rdata_path, sprintf("%s_Null_Model_Free_Shuffles_%s.Rds",
#                                                                          dataset_ID,
#                                                                          univariate_feature_set)))
# } else {
#   model_free_shuffle_null_res <- readRDS(paste0(rdata_path, sprintf("%s_Null_Model_Free_Shuffles_%s.Rds",
#                                                                     dataset_ID,
#                                                                     univariate_feature_set)))
# }

################################################################################
# Create ten folds to use for all analyses
################################################################################
subjects_to_use <- read.csv(paste0(data_path, sprintf("%s_samples_with_univariate_%s_and_pairwise_%s.csv",
                                                      dataset_ID,
                                                      univariate_feature_set,
                                                      pairwise_feature_set)))
if (!file.exists(paste0(rdata_path, dataset_ID, "_samples_per_10_folds.Rds"))) {
  # Make folds
  set.seed(127)
  k = 10
  sample_folds <- caret::createFolds(subjects_to_use$group, k = k, list = TRUE, returnTrain = FALSE)
  
  # Save to Rds file
  saveRDS(sample_folds, file=paste0(rdata_path, dataset_ID, "_samples_per_10_folds.Rds"))
} else {
  sample_folds <- readRDS(paste0(rdata_path, dataset_ID, "_samples_per_10_folds.Rds"))
}

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
    
    # Run given weighting for 10-fold CV linear SVM
    if (!file.exists(paste0(rdata_path, sprintf("%s_wise_CV_linear_SVM_%s_%s.Rds",
                                                grouping_type, feature_set, weighting_name)))) {
      group_wise_SVM_CV_weighting <- run_univariate_cv_svm_by_input_var(rdata_path = rdata_path,
                                                                        feature_set = feature_set,
                                                                        test_package = test_package,
                                                                        svm_kernel = kernel,
                                                                        grouping_var = grouping_var,
                                                                        flds = sample_folds,
                                                                        svm_feature_var = SVM_feature_var,
                                                                        out_of_sample_only = TRUE,
                                                                        use_inv_prob_weighting = use_inv_prob_weighting,
                                                                        noise_procs = noise_procs)
      saveRDS(group_wise_SVM_CV_weighting, file=paste0(rdata_path, 
                                                       sprintf("%s_wise_CV_linear_SVM_%s_%s.Rds",
                                                               grouping_type,
                                                               feature_set, 
                                                               weighting_name)))
    }
    
    #### Calculate balanced accuracy across all folds
    if (!file.exists(paste0(rdata_path, sprintf("%s_wise_CV_linear_SVM_%s_%s_balacc.Rds",
                                                grouping_type, feature_set, weighting_name)))) {
      group_wise_SVM_balanced_accuracy <- group_wise_SVM_CV_weighting %>%
        group_by(grouping_var, Noise_Proc, Sample_Type) %>%
        summarise(accuracy = sum(Prediction_Correct) / n(),
                  balanced_accuracy = caret::confusionMatrix(data = Predicted_Diagnosis,
                                                             reference = Actual_Diagnosis)$byClass[["Balanced Accuracy"]])
      
      saveRDS(group_wise_SVM_balanced_accuracy, file=paste0(rdata_path, sprintf("%s_wise_CV_linear_SVM_%s_%s_balacc.Rds",
                                                                                grouping_type, feature_set, weighting_name)))
    }
  }
}