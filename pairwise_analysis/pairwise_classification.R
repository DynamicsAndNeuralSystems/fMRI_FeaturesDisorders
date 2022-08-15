#-------------------------------------------------------------------------------
# Define paths
#-------------------------------------------------------------------------------
# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--github_dir", default="/headnode1/abry4213/github/fMRI_FeaturesDisorders/")
parser$add_argument("--data_path", default="/headnode1/abry4213/data/UCLA_Schizophrenia/")
parser$add_argument("--num_null_permutations", default=1000)
parser$add_argument("--pairwise_feature_set", default="pyspi_19")
parser$add_argument("--univariate_feature_set", default="catch22")
parser$add_argument("--noise_procs", default=c(""), nargs="*", action="append")
parser$add_argument("--noise_proc_for_null", default=c(""))
parser$add_argument("--dataset_ID", default="UCLA_Schizophrenia")

# univariate_feature_set <- "catch22"
# pairwise_feature_set <- "pyspi_19"
# subject_csv <- "participants.csv"
# github_dir <- "/headnode1/abry4213/github/fMRI_FeaturesDisorders/"

# UCLA schizophrenia
# rdata_path <- "/headnode1/abry4213/data/UCLA_Schizophrenia/Rdata/"
# pydata_path <- "/headnode1/abry4213/data/UCLA_Schizophrenia/pydata/"
# data_path <- "/headnode1/abry4213/data/UCLA_Schizophrenia/"
# dataset_ID <- "UCLA_Schizophrenia"
# noise_procs <- c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER")
# noise_proc_for_null <- "AROMA+2P+GMR"

# ABIDE ASD
# rdata_path <- "/headnode1/abry4213/data/ABIDE_ASD/Rdata/"
# pydata_path <- "/headnode1/abry4213/data/ABIDE_ASD/pydata/"
# data_path <- "/headnode1/abry4213/data/ABIDE_ASD/"
# dataset_ID <- "ABIDE_ASD"
# noise_procs <- c("FC1000")
# noise_proc_for_null <- "FC1000"

# Parse input arguments
args <- parser$parse_args()
github_dir <- args$github_dir
data_path <- args$data_path
num_null_permutations <- args$num_null_permutations
pairwise_feature_set <- args$pairwise_feature_set
univariate_feature_set <- args$univariate_feature_set
noise_procs <- args$noise_procs
noise_proc_for_null <- args$noise_proc_for_null
dataset_ID <- args$dataset_ID

plot_dir <- paste0(data_path, "plots/")
icesTAF::mkdir(plot_dir)

rdata_path <- paste0(data_path, "Rdata/")
pydata_path <- paste0(data_path, "pydata/")

# Set the seed
set.seed(127)

# Load tidyverse
library(tidyverse)

# Unlist noise-processing methods
tryCatch({
  noise_procs <- unlist(noise_procs)
}, error = function(e) {})


################################################################################
# Source helper scripts
################################################################################
# Set working directory to file location
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
helper_script_dir = "../helper_functions/classification/"
source(paste0(helper_script_dir, "Linear_SVM.R"))
source(paste0(helper_script_dir, "Null_distributions.R"))

################################################################################
# Load pyspi data
################################################################################
pyspi_data_file <- sprintf("%s/%s_pairwise_%s.Rds",
                           pydata_path, dataset_ID, pairwise_feature_set) 
pyspi_data <- readRDS(pyspi_data_file) %>%
  dplyr::filter(Noise_Proc %in% noise_proc_for_null) %>%
  dplyr::mutate(Diagnosis = stringr::str_to_sentence(Diagnosis))

SPI_directionality <- read.csv(paste0(github_dir, "pairwise_analysis/SPI_Direction_Info.csv"))

################################################################################
# Create ten folds to use for all analyses
################################################################################
subjects_to_use <- readRDS(paste0(data_path, sprintf("%s_samples_with_univariate_%s_and_pairwise_%s.Rds",
                                                     dataset_ID,
                                                     univariate_feature_set,
                                                     pairwise_feature_set)))


if (!file.exists(paste0(rdata_path, dataset_ID, "_samples_per_10_folds.Rds"))) {
  # Make folds
  set.seed(127)
  k = 10
  sample_folds <- caret::createFolds(subjects_to_use$Diagnosis, k = k, list = TRUE, returnTrain = FALSE)
  
  # Save to Rds file
  saveRDS(sample_folds, file=paste0(rdata_path, dataset_ID, "_samples_per_10_folds.Rds"))
} else {
  sample_folds <- readRDS(paste0(rdata_path, dataset_ID, "_samples_per_10_folds.Rds"))
}

################################################################################
# Define weighting parameters
################################################################################

weighting_param_df <- data.frame(name = c("inv_prob"),
                                 use_inv_prob_weighting = c(TRUE))

SVM_grouping_params <- data.frame(grouping_var = c("SPI"),
                                  SVM_feature_var = c("region_pair"))

# Use a linear kernel
kernel = "linear"

for (i in 1:nrow(SVM_grouping_params)) {
  grouping_var = SVM_grouping_params$grouping_var[i]
  SVM_feature_var = SVM_grouping_params$SVM_feature_var[i]
  
  #### 10-fold linear SVM with different weights
  # Iterate over weighting_param_df 
  for (j in 1:nrow(weighting_param_df)) {
    weighting_name <- weighting_param_df$name[j]
    use_inv_prob_weighting <- weighting_param_df$use_inv_prob_weighting[j]
    
    # Run given weighting for 10-fold CV linear SVM
    if (!file.exists(paste0(rdata_path, sprintf("%s_pairwise_CV_linear_SVM_%s_%s.Rds",
                                                grouping_var, 
                                                pairwise_feature_set, 
                                                weighting_name)))) {
      tryCatch({group_wise_SVM_CV_weighting <- run_pairwise_cv_svm_by_input_var(pairwise_data = pyspi_data,
                                                                                data_path = data_path,
                                                                                SPI_directionality = SPI_directionality,
                                                                                svm_kernel = "linear",
                                                                                num_k_folds = 10,
                                                                                flds = sample_folds,
                                                                                grouping_var = grouping_var,
                                                                                svm_feature_var = SVM_feature_var,
                                                                                noise_proc = noise_proc_for_null,
                                                                                out_of_sample_only = TRUE,
                                                                                use_inv_prob_weighting = use_inv_prob_weighting,
                                                                                shuffle_labels = FALSE)
      saveRDS(group_wise_SVM_CV_weighting, file=paste0(rdata_path, 
                                                       sprintf("%s_pairwise_CV_linear_SVM_%s_%s.Rds",
                                                               grouping_var,
                                                               pairwise_feature_set, 
                                                               weighting_name)))
      }, error = function(e) {
        cat("Could not run linear SVM for", grouping_var, pairwise_feature_set, ".\n")
        message(e)
      })
    } else {
      group_wise_SVM_CV_weighting <- readRDS(paste0(rdata_path, 
                                                    sprintf("%s_pairwise_CV_linear_SVM_%s_%s.Rds",
                                                            grouping_var,
                                                            pairwise_feature_set, 
                                                            weighting_name)))
    }
    
    #### Calculate balanced accuracy across all folds
    if (!file.exists(paste0(rdata_path, sprintf("%s_pairwise_CV_linear_SVM_%s_%s_balacc.Rds",
                                                grouping_var, 
                                                pairwise_feature_set, 
                                                weighting_name)))) {
      group_wise_SVM_balanced_accuracy <- group_wise_SVM_CV_weighting %>%
        group_by(grouping_var, Noise_Proc, Sample_Type) %>%
        summarise(accuracy = sum(Prediction_Correct) / n(),
                  balanced_accuracy = caret::confusionMatrix(data = Predicted_Diagnosis,
                                                             reference = Actual_Diagnosis)$byClass[["Balanced Accuracy"]])
      
      saveRDS(group_wise_SVM_balanced_accuracy, file=paste0(rdata_path, sprintf("%s_pairwise_CV_linear_SVM_%s_%s_balacc.Rds",
                                                                                grouping_var, 
                                                                                pairwise_feature_set, 
                                                                                weighting_name)))
    } else {
      group_wise_SVM_balanced_accuracy <- readRDS(paste0(rdata_path, sprintf("%s_pairwise_CV_linear_SVM_%s_%s_balacc.Rds",
                                                                             grouping_var, 
                                                                             pairwise_feature_set, 
                                                                             weighting_name)))
    }
  }
}
