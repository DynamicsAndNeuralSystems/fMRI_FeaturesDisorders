#-------------------------------------------------------------------------------
# Define paths
#-------------------------------------------------------------------------------
# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--github_dir", default="/headnode1/abry4213/github/")
parser$add_argument("--data_path", default="/headnode1/abry4213/data/UCLA_Schizophrenia/")
parser$add_argument("--sample_metadata_file", default="UCLA_Schizophrenia_sample_metadata.Rds")
parser$add_argument("--pairwise_feature_set", default="pyspi14")
parser$add_argument("--univariate_feature_set", default="catch22")
parser$add_argument("--noise_procs", default=c(""))
parser$add_argument("--noise_proc_for_null", default=c(""))
parser$add_argument("--dataset_ID", default="UCLA_Schizophrenia")
parser$add_argument("--email")
parser$add_argument("--add_catch2", action="store_true", default=FALSE)

# Parse input arguments
args <- parser$parse_args()
github_dir <- args$github_dir
data_path <- args$data_path
rdata_path <- args$rdata_path
pairwise_feature_set <- args$pairwise_feature_set
univariate_feature_set <- args$univariate_feature_set
noise_procs <- args$noise_procs
noise_proc_for_null <- args$noise_proc_for_null
dataset_ID <- args$dataset_ID
sample_metadata_file <- args$sample_metadata_file
email <- args$email
add_catch2 <- args$add_catch2

cat("noise_procs:", noise_procs, "\n")
cat("noise_proc_for_null:", noise_proc_for_null, "\n")

# univariate_feature_set <- "catch22"
# pairwise_feature_set <- "pyspi14"
# github_dir <- "~/github/"
# email <- "abry4213@uni.sydney.edu.au"

# UCLA schizophrenia
# data_path <- "~/data/UCLA_Schizophrenia/"
# dataset_ID <- "UCLA_Schizophrenia"
# sample_metadata_file <- "UCLA_Schizophrenia_sample_metadata.Rds"
# noise_procs <- "AROMA+2P;AROMA+2P+GMR;AROMA+2P+DiCER"
# noise_proc_for_null <- "AROMA+2P+GMR"

# ABIDE ASD
# data_path <- "~/data/ABIDE_ASD/"
# sample_metadata_file <- "ABIDE_ASD_sample_metadata.Rds"
# dataset_ID <- "ABIDE_ASD"
# noise_procs <- c("FC1000")
# noise_proc_for_null <- "FC1000"

rdata_path <- paste0(data_path, "processed_data/Rdata/")
plot_dir <- paste0(data_path, "plots/")

icesTAF::mkdir(plot_dir)

# Set the seed
set.seed(127)

# Load tidyverse
library(tidyverse)

# Unlist noise-processing methods
tryCatch({
  noise_procs <- stringr::str_split(noise_procs, ";")[[1]]
  noise_procs <- unlist(noise_procs)
}, error = function(e) {})


#-------------------------------------------------------------------------------
# Source helper scripts
#-------------------------------------------------------------------------------
helper_script_dir = paste0(github_dir, "fMRI_FeaturesDisorders/helper_functions/classification/")
source(paste0(helper_script_dir, "Linear_SVM.R"))
source(paste0(helper_script_dir, "Null_distributions.R"))

# Load sample metadata
sample_metadata <- readRDS(paste0(data_path, sample_metadata_file))

################################################################################
# Create 10 repeats of 10 folds to use for all analyses
################################################################################
subjects_to_use <- readRDS(paste0(rdata_path, sprintf("%s_samples_with_univariate_%s_and_pairwise_%s_filtered.Rds",
                                                      dataset_ID,
                                                      univariate_feature_set,
                                                      pairwise_feature_set)))

if (!file.exists(paste0(rdata_path, dataset_ID, "_samples_per_10_folds_10_repeats.Rds"))) {
  sample_folds <- list()
  
  for (i in 1:10) {
    # Make folds
    k = 10
    samples_with_diagnosis <- subjects_to_use %>%
      left_join(., sample_metadata)
    sample_folds_i <- caret::createFolds(samples_with_diagnosis$Diagnosis, k = k, list = TRUE, returnTrain = FALSE)
    sample_folds[[i]] <- sample_folds_i
  }
  
  # Save RDS file
  saveRDS(sample_folds, paste0(rdata_path, dataset_ID, "_samples_per_10_folds_10_repeats.Rds"))
} else {
  sample_folds <- readRDS(paste0(rdata_path, dataset_ID, "_samples_per_10_folds_10_repeats.Rds"))
}

################################################################################
# Define weighting parameters
################################################################################

grouping_param_df <- data.frame(grouping_type = c("ROI", "Feature", "Combo"),
                                grouping_var = c("Brain_Region", "Feature", "Combo"),
                                SVM_feature_var = c("Feature", "Brain_Region", "Combo")) 

# Use a linear kernel
kernel = "linear"
weighting_name <- "inv_prob"
use_inv_prob_weighting <- TRUE

for (i in 1:nrow(grouping_param_df)) {
  grouping_type = grouping_param_df$grouping_type[i]
  grouping_var = grouping_param_df$grouping_var[i]
  SVM_feature_var = grouping_param_df$SVM_feature_var[i]
  
  #### 10-fold linear SVM with different weights
  
  # Run given weighting for 10-fold CV linear SVM
  if (!file.exists(paste0(rdata_path, sprintf("%s_wise_CV_linear_SVM_%s_%s.Rds",
                                              grouping_type, 
                                              univariate_feature_set, 
                                              weighting_name)))) {
    
    group_wise_SVM_CV_weighting <- 1:length(sample_folds) %>%
      purrr::map_df( ~ run_univariate_cv_svm_by_input_var(data_path = data_path,
                                                           rdata_path = rdata_path,
                                                           dataset_ID = dataset_ID,
                                                           sample_metadata = sample_metadata,
                                                           univariate_feature_set = univariate_feature_set,
                                                           pairwise_feature_set = pairwise_feature_set,
                                                           svm_kernel = kernel,
                                                           grouping_var = grouping_var,
                                                           flds = sample_folds[[.x]],
                                                           repeat_number = .x,
                                                           svm_feature_var = SVM_feature_var,
                                                           out_of_sample_only = TRUE,
                                                           use_inv_prob_weighting = use_inv_prob_weighting,
                                                           noise_procs = noise_procs))
    
    saveRDS(group_wise_SVM_CV_weighting, file=paste0(rdata_path, 
                                                     sprintf("%s_wise_CV_linear_SVM_%s_%s.Rds",
                                                             grouping_type,
                                                             univariate_feature_set, 
                                                             weighting_name)))
  } else {
    group_wise_SVM_CV_weighting <- readRDS(paste0(rdata_path, 
                                                  sprintf("%s_wise_CV_linear_SVM_%s_%s.Rds",
                                                          grouping_type,
                                                          univariate_feature_set, 
                                                          weighting_name)))
  }
    
    #### Calculate averaged balanced accuracy across all folds and repeats
    if (!file.exists(paste0(rdata_path, sprintf("%s_wise_CV_linear_SVM_%s_%s_balacc_across_repeats.Rds",
                                                grouping_type, 
                                                univariate_feature_set, 
                                                weighting_name)))) {
      # First find balanced accuracy per repeat across folds
      group_wise_SVM_balanced_accuracy <- group_wise_SVM_CV_weighting %>%
        group_by(grouping_var, Noise_Proc, Sample_Type, fold_number, repeat_number) %>%
        summarise(accuracy = sum(Prediction_Correct) / n(),
                  balanced_accuracy = caret::confusionMatrix(data = Predicted_Diagnosis,
                                                             reference = Actual_Diagnosis)$byClass[["Balanced Accuracy"]])
      saveRDS(group_wise_SVM_balanced_accuracy, file=paste0(rdata_path, sprintf("%s_wise_CV_linear_SVM_%s_%s_balacc.Rds",
                                                                                grouping_type, 
                                                                                univariate_feature_set, 
                                                                                weighting_name))) 
      
      # Then find averaged balanced accuracy across all repeats
      group_wise_SVM_balanced_accuracy_across_repeats <- group_wise_SVM_balanced_accuracy %>%
        group_by(grouping_var, Noise_Proc, Sample_Type) %>%
        summarise(mean_accuracy = mean(accuracy, na.rm=T),
                  mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T)) %>%
        dplyr::rename("accuracy" = "mean_accuracy",
                      "balanced_accuracy" = "mean_balanced_accuracy")
      saveRDS(group_wise_SVM_balanced_accuracy_across_repeats, file=paste0(rdata_path, sprintf("%s_wise_CV_linear_SVM_%s_%s_balacc_across_repeats.Rds",
                                                                                grouping_type, 
                                                                                univariate_feature_set, 
                                                                                weighting_name)))
    } else {
      group_wise_SVM_balanced_accuracy_across_repeats <- readRDS(paste0(rdata_path, 
                                                                        sprintf("%s_wise_CV_linear_SVM_%s_%s_balacc_across_repeats.Rds",
                                                                             grouping_type, 
                                                                             univariate_feature_set, 
                                                                             weighting_name)))
    }
}
