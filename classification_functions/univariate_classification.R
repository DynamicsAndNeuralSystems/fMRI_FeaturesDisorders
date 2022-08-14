#-------------------------------------------------------------------------------
# Define paths
#-------------------------------------------------------------------------------
# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--github_dir", default="/headnode1/abry4213/github/fMRI_FeaturesDisorders/")
parser$add_argument("--data_path", default="/headnode1/abry4213/data/UCLA_Schizophrenia/")
parser$add_argument("--rdata_path", default="/headnode1/abry4213/data/UCLA_Schizophrenia/Rdata/")
parser$add_argument("--univariate_feature_set", default="catch22")
parser$add_argument("--noise_procs", default=c(""), nargs="*", action="append")
parser$add_argument("--dataset_ID", default="UCLA_Schizophrenia")

# univariate_feature_set <- "catch22"
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
files.sources = list.files(helper_script_dir, pattern=".R$", full.names = T) %>% .[!(str_detect(., "cluster"))]
invisible(sapply(files.sources, source))

################################################################################
# Generate model-free shuffle null distribution
################################################################################
# if (!file.exists(paste0(rdata_path, sprintf("Null_Model_Free_Shuffles_%s.Rds",
#                                             univariate_feature_set)))) {
#   model_free_shuffle_null_res <- run_model_free_n_shuffles(num_shuffles = 100000,
#                                                            feature_set = univariate_feature_set,
#                                                            dataset_ID = dataset_ID,
#                                                            rdata_path = rdata_path)
#   saveRDS(model_free_shuffle_null_res, file = paste0(rdata_path, sprintf("Null_Model_Free_Shuffles_%s.Rds",
#                                                                          univariate_feature_set)))
# } else {
#   model_free_shuffle_null_res <- readRDS(paste0(rdata_path, sprintf("Null_Model_Free_Shuffles_%s.Rds",
#                                                                     univariate_feature_set)))
# }

################################################################################
# Create ten folds to use for all analyses
################################################################################

subjects_to_use <- readRDS(paste0(rdata_path, dataset_ID, "_Subjects_with_Univariate_and_Pairwise.Rds"))

if (!file.exists(paste0(rdata_path, "Subjects_per_10_folds.Rds"))) {
  # Make folds
  set.seed(127)
  k = 10
  subject_folds <- caret::createFolds(subjects_to_use$group, k = k, list = TRUE, returnTrain = FALSE)
  
  # Save to Rds file
  saveRDS(subject_folds, file=paste0(rdata_path, "Subjects_per_10_folds.Rds"))
} else {
  subject_folds <- readRDS(paste0(rdata_path, "Subjects_per_10_folds.Rds"))
}

################################################################################
# Define weighting parameters
################################################################################

grouping_param_df <- data.frame(grouping_type = c("ROI", "Feature", "Combo"),
                                grouping_var = c("Brain_Region", "Feature", "Combo"),
                                SVM_feature_var = c("Feature", "Brain_Region", "Combo")) 

weighting_param_df <- data.frame(name = c("inv_prob"),
                                 use_inv_prob_weighting = c(TRUE),
                                 use_SMOTE = c(FALSE))

