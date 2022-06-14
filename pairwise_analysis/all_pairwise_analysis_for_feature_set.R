# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--github_dir", default="/home/abry4213/github/")
parser$add_argument("--rdata_path", default="/home/abry4213/data/scz/UCLA/Rdata/")
parser$add_argument("--feature_set", default="catch22")
# github_dir <- "D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders/"
# rdata_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/Rdata/"
# study <- "D:/Virtual_Machines/Shared_Folder/PhD_work/"
# data_path <- paste0(study, "data/scz/UCLA/")
# pydata_path <- paste0(study, "data/scz/UCLA/pydata/")
# feature_set <- "catch22"

# Parse input arguments
args <- parser$parse_args()
github_dir <- args$github_dir
rdata_path <- args$rdata_path
feature_set <- args$feature_set

### Source functions
# Main
source(paste0(github_dir, "helper_functions/Linear_SVM.R"))
source(paste0(github_dir, "helper_functions/Visualization.R"))
source(paste0(github_dir, "helper_functions/Null_distributions.R"))

set.seed(127)

# Compare just AROMA+2P+GMR
noise_proc = "AROMA+2P+GMR"

# Use e1071 SVM with a linear kernel
test_package = "e1071"
kernel = "linear"

################################################################################
# Load pyspi data
################################################################################
pyspi_data <- readRDS(paste0(pydata_path, "UCLA_all_subject_pyspi_AROMA_2P_GMR_filtered_zscored.Rds")) %>%
  mutate(group = stringr::str_to_sentence(group))
# pyspi_data_pearson <- subset(pyspi_data, SPI=="cov_EmpiricalCovariance")
SPI_directionality <- read.csv("SPI_Direction_Info.csv")

################################################################################
# Generate model-free shuffle null distribution
################################################################################
if (!file.exists(paste0(rdata_path, sprintf("Null_Model_Free_Shuffles_%s.Rds",
                                            feature_set)))) {
  model_free_shuffle_null_res <- run_model_free_n_shuffles(num_shuffles = 100000,
                                                           feature_set = feature_set,
                                                           rdata_path = rdata_path)
  saveRDS(model_free_shuffle_null_res, file = paste0(rdata_path, sprintf("Null_Model_Free_Shuffles_%s.Rds",
                                                                         feature_set)))
} else {
  model_free_shuffle_null_res <- readRDS(paste0(rdata_path, sprintf("Null_Model_Free_Shuffles_%s.Rds",
                                                                    feature_set)))
}


################################################################################
# Define weighting parameters
################################################################################

weighting_param_df <- data.frame(name = c("unweighted", "inv_prob", "SMOTE"),
                                 use_inv_prob_weighting = c(FALSE, TRUE, FALSE),
                                 use_SMOTE = c(FALSE, FALSE, TRUE))

################################################################################
# Per-feature, all ROI pair combinations
################################################################################
feature_set <- "pyspi_19"

#### 10-fold linear SVM with different weights
# Iterate over weighting_param_df 
for (i in 1:nrow(weighting_param_df)) {
  weighting_name <- weighting_param_df$name[i]
  use_inv_prob_weighting <- weighting_param_df$use_inv_prob_weighting[i]
  use_SMOTE <- weighting_param_df$use_SMOTE[i]
  
  # Run given weighting for 10-fold CV linear SVM
  if (!file.exists(paste0(rdata_path, sprintf("pyspi_ROI_pairwise_CV_linear_SVM_%s_%s.Rds",
                                              feature_set, weighting_name)))) {
    pyspi_SPI_pairwise_SVM_CV_weighting <- run_pairwise_cv_svm_by_input_var(pairwise_data = pyspi_data,
                                                                            SPI_directionality = SPI_directionality,
                                                                            svm_kernel = "linear",
                                                                            grouping_var = "SPI",
                                                                            svm_feature_var = "region_pair",
                                                                            test_package = "e1071",
                                                                            noise_proc = "AROMA+2P+GMR",
                                                                            return_all_fold_metrics = TRUE,
                                                                            use_inv_prob_weighting = use_inv_prob_weighting,
                                                                            use_SMOTE = use_SMOTE,
                                                                            shuffle_labels = FALSE)
    saveRDS(pyspi_SPI_pairwise_SVM_CV_weighting, file=paste0(rdata_path, 
                                                           sprintf("pyspi_SPI_pairwise_CV_linear_SVM_%s_%s.Rds",
                                                                   feature_set, weighting_name)))
  }
}

# #### Calculate p values from model-free shuffle null distribution
# for (weighting_name in unique(weighting_param_df$name)) {
#   if (!file.exists(paste0(rdata_path, sprintf("pyspi_SPI_pairwise_CV_linear_SVM_%s_%s_model_free_shuffle_pvals.Rds",
#                                               feature_set, weighting_name)))) {
#     pyspi_SPI_pairwise_SVM_CV_weighting <- readRDS(paste0(rdata_path, 
#                                                    sprintf("pyspi_SPI_pairwise_CV_linear_SVM_%s_%s.Rds",
#                                                            feature_set, weighting_name)))
#     
#     # Calculate p-values
#     pvalues <- calc_empirical_nulls(class_res = pyspi_SPI_pairwise_SVM_CV_weighting,
#                                     null_data = model_free_shuffle_null_res,
#                                     is_data_averaged = FALSE,
#                                     grouping_var = "SPI")
#     
#     saveRDS(pvalues, file=paste0(rdata_path, sprintf("pyspi_SPI_pairwise_CV_linear_SVM_%s_%s_model_free_shuffle_pvals.Rds",
#                                                      feature_set, weighting_name)))
#   }
# }
# 
# #### Generate empirical null model distributions per SPI
# for (i in 1:nrow(weighting_param_df)) {
#   weighting_name <- weighting_param_df$name[i]
#   use_inv_prob_weighting <- weighting_param_df$use_inv_prob_weighting[i]
#   use_SMOTE <- weighting_param_df$use_SMOTE[i]
#   
#   # Generate null-model fits distribution
#   if (!file.exists(paste0(rdata_path, sprintf("pyspi_SPI_pairwise_model_permutation_null_%s_%s.Rds",
#                                               feature_set, weighting_name)))) {
#     model_permutation_null_weighting <- run_null_model_n_permutations_pairwise(pairwise_data = pyspi_data,
#                                                                                noise_proc = "AROMA+2P+GMR",
#                                                                                feature_set = feature_set,
#                                                                                test_package = "e1071",
#                                                                                svm_kernel = "linear",
#                                                                                grouping_var = "SPI",
#                                                                                svm_feature_var = "region_pair",
#                                                                                SPI_directionality = SPI_directionality,
#                                                                                num_permutations = 5,
#                                                                                use_inv_prob_weighting = use_inv_prob_weighting,
#                                                                                use_SMOTE = use_SMOTE)
#     
#     saveRDS(model_permutation_null_weighting, file=paste0(rdata_path, sprintf("pyspi_SPI_pairwise_model_permutation_null_%s_%s.Rds",
#                                                                               feature_set, weighting_name)))
#   } else {
#     model_permutation_null_weighting <- readRDS(paste0(rdata_path, sprintf("pyspi_SPI_pairwise_model_permutation_null_%s_%s.Rds",
#                                                                            feature_set, weighting_name)))
#   }
#   
#   # Empirically derive p-values based on null model fits distribution
#   if (!file.exists(paste0(rdata_path, sprintf("pyspi_SPI_pairwise_CV_linear_SVM_%s_%s_null_model_fit_pvals.Rds",
#                                               feature_set, weighting_name)))) {
#     pyspi_SPI_pairwise_SVM_CV_weighting <- readRDS(paste0(rdata_path, 
#                                                    sprintf("pyspi_SPI_pairwise_CV_linear_SVM_%s_%s.Rds",
#                                                            feature_set, weighting_name)))
#     
#     # Calculate p-values
#     pvalues <- calc_empirical_nulls(class_res = pyspi_SPI_pairwise_SVM_CV_weighting,
#                                     null_data = model_permutation_null_weighting,
#                                     is_data_averaged = FALSE,
#                                     grouping_var = "SPI")
#     
#     saveRDS(pvalues, file=paste0(rdata_path, sprintf("pyspi_SPI_pairwise_CV_linear_SVM_%s_%s_null_model_fit_pvals.Rds",
#                                                      feature_set, weighting_name)))
#   }
# }
# 
# 
# ################################################################################
# # TS Feature-wise analysis
# ################################################################################
# 
# #### 10-fold linear SVM with different weights
# # Iterate over weighting_param_df 
# for (i in 1:nrow(weighting_param_df)) {
#   weighting_name <- weighting_param_df$name[i]
#   use_inv_prob_weighting <- weighting_param_df$use_inv_prob_weighting[i]
#   use_SMOTE <- weighting_param_df$use_SMOTE[i]
#   
#   # Run given weighting for 10-fold CV linear SVM
#   if (!file.exists(paste0(rdata_path, sprintf("Feature_wise_CV_linear_SVM_%s_%s.Rds",
#                                               feature_set, weighting_name)))) {
#     feature_wise_SVM_CV_weighting <- run_cv_svm_by_input_var(rdata_path = rdata_path,
#                                                             feature_set = feature_set,
#                                                             test_package = test_package,
#                                                             svm_kernel = kernel,
#                                                             grouping_var = "Feature",
#                                                             svm_feature_var = "Brain_Region",
#                                                             use_inv_prob_weighting = use_inv_prob_weighting,
#                                                             use_SMOTE = use_SMOTE,
#                                                             noise_procs = noise_procs)
#     saveRDS(feature_wise_SVM_CV_weighting, file=paste0(rdata_path, 
#                                                       sprintf("Feature_wise_CV_linear_SVM_%s_%s.Rds",
#                                                               feature_set, weighting_name)))
#   }
# }
# 
# #### Calculate p values from model-free shuffle null distribution
# for (weighting_name in unique(weighting_param_df$name)) {
#   if (!file.exists(paste0(rdata_path, sprintf("Feature_wise_CV_linear_SVM_%s_%s_model_free_shuffle_pvals.Rds",
#                                               feature_set, weighting_name)))) {
#     feature_wise_SVM_CV_weighting <- readRDS(paste0(rdata_path, 
#                                                    sprintf("Feature_wise_CV_linear_SVM_%s_%s.Rds",
#                                                            feature_set, weighting_name)))
#     
#     # Calculate p-values
#     pvalues <- calc_empirical_nulls(class_res = feature_wise_SVM_CV_weighting,
#                                     null_data = model_free_shuffle_null_res,
#                                     grouping_var = "Brain_Region")
#     
#     saveRDS(pvalues, file=paste0(rdata_path, sprintf("Feature_wise_CV_linear_SVM_%s_%s_model_free_shuffle_pvals.Rds",
#                                                      feature_set, weighting_name)))
#   }
# }
# 
# #### Generate empirical null model distributions per brain region
# for (i in 1:nrow(weighting_param_df)) {
#   weighting_name <- weighting_param_df$name[i]
#   use_inv_prob_weighting <- weighting_param_df$use_inv_prob_weighting[i]
#   use_SMOTE <- weighting_param_df$use_SMOTE[i]
#   
#   # Generate null-model fits distribution
#   if (!file.exists(paste0(rdata_path, sprintf("Feature_wise_model_permutation_null_%s_%s.Rds",
#                                               feature_set, weighting_name)))) {
#     model_permutation_null_weighting <- run_null_model_n_permutations(rdata_path,
#                                                                       feature_set = feature_set,
#                                                                       noise_procs = noise_procs,
#                                                                       grouping_var = "Feature",
#                                                                       svm_feature_var = "Brain_Region",
#                                                                       num_permutations = 40,
#                                                                       use_inv_prob_weighting = use_inv_prob_weighting,
#                                                                       use_SMOTE = use_SMOTE)
#     
#     saveRDS(model_permutation_null_weighting, file=paste0(rdata_path, sprintf("Feature_wise_model_permutation_null_%s_%s.Rds",
#                                                                               feature_set, weighting_name)))
#   }
#   
#   # Empirically derive p-values based on null model fits distribution
#   if (!file.exists(paste0(rdata_path, sprintf("Feature_wise_CV_linear_SVM_%s_%s_null_model_fit_pvals.Rds",
#                                               feature_set, weighting_name)))) {
#     feature_wise_SVM_CV_weighting <- readRDS(paste0(rdata_path, 
#                                                    sprintf("Feature_wise_CV_linear_SVM_%s_%s.Rds",
#                                                            feature_set, weighting_name)))
#     
#     # Calculate p-values
#     pvalues <- calc_empirical_nulls(class_res = feature_wise_SVM_CV_weighting,
#                                     null_data = model_permutation_null_weighting,
#                                     grouping_var = "Brain_Region")
#     
#     saveRDS(pvalues, file=paste0(rdata_path, sprintf("Feature_wise_CV_linear_SVM_%s_%s_null_model_fit_pvals.Rds",
#                                                      feature_set, weighting_name)))
#   }
# }
# 
# ################################################################################
# # ROI/Feature Combo-wise analysis
# ################################################################################
# 
# #### 10-fold linear SVM with different weights
# # Iterate over weighting_param_df 
# for (i in 1:nrow(weighting_param_df)) {
#   weighting_name <- weighting_param_df$name[i]
#   use_inv_prob_weighting <- weighting_param_df$use_inv_prob_weighting[i]
#   use_SMOTE <- weighting_param_df$use_SMOTE[i]
#   
#   # Run given weighting for 10-fold CV linear SVM
#   if (!file.exists(paste0(rdata_path, sprintf("Combo_wise_CV_linear_SVM_%s_%s.Rds",
#                                               feature_set, weighting_name)))) {
#     combo_wise_SVM_CV_weighting <- run_cv_svm_by_input_var(rdata_path = rdata_path,
#                                                              feature_set = feature_set,
#                                                              test_package = test_package,
#                                                              svm_kernel = kernel,
#                                                              grouping_var = "Combo",
#                                                              svm_feature_var = "Combo",
#                                                              use_inv_prob_weighting = use_inv_prob_weighting,
#                                                              use_SMOTE = use_SMOTE,
#                                                              noise_procs = noise_procs)
#     saveRDS(combo_wise_SVM_CV_weighting, file=paste0(rdata_path, 
#                                                        sprintf("Combo_wise_CV_linear_SVM_%s_%s.Rds",
#                                                                feature_set, weighting_name)))
#   }
# }
# 
# #### Calculate p values from model-free shuffle null distribution
# for (weighting_name in unique(weighting_param_df$name)) {
#   if (!file.exists(paste0(rdata_path, sprintf("Combo_wise_CV_linear_SVM_%s_%s_model_free_shuffle_pvals.Rds",
#                                               feature_set, weighting_name)))) {
#     combo_wise_SVM_CV_weighting <- readRDS(paste0(rdata_path, 
#                                                     sprintf("Combo_wise_CV_linear_SVM_%s_%s.Rds",
#                                                             feature_set, weighting_name)))
#     
#     # Calculate p-values
#     pvalues <- calc_empirical_nulls(class_res = combo_wise_SVM_CV_weighting,
#                                     null_data = model_free_shuffle_null_res,
#                                     grouping_var = "Brain_Region")
#     
#     saveRDS(pvalues, file=paste0(rdata_path, sprintf("Combo_wise_CV_linear_SVM_%s_%s_model_free_shuffle_pvals.Rds",
#                                                      feature_set, weighting_name)))
#   }
# }
# 
# #### Generate empirical null model distributions per brain region
# for (i in 1:nrow(weighting_param_df)) {
#   weighting_name <- weighting_param_df$name[i]
#   use_inv_prob_weighting <- weighting_param_df$use_inv_prob_weighting[i]
#   use_SMOTE <- weighting_param_df$use_SMOTE[i]
#   
#   # Generate null-model fits distribution
#   if (!file.exists(paste0(rdata_path, sprintf("Combo_wise_model_permutation_null_%s_%s.Rds",
#                                               feature_set, weighting_name)))) {
#     model_permutation_null_weighting <- run_null_model_n_permutations(rdata_path,
#                                                                       feature_set = feature_set,
#                                                                       noise_procs = noise_procs,
#                                                                       grouping_var = "Combo",
#                                                                       svm_feature_var = "Combo",
#                                                                       num_permutations = 100,
#                                                                       use_inv_prob_weighting = use_inv_prob_weighting,
#                                                                       use_SMOTE = use_SMOTE)
#     
#     saveRDS(model_permutation_null_weighting, file=paste0(rdata_path, sprintf("Combo_wise_model_permutation_null_%s_%s.Rds",
#                                                                               feature_set, weighting_name)))
#   }
#   
#   # Empirically derive p-values based on null model fits distribution
#   if (!file.exists(paste0(rdata_path, sprintf("Combo_wise_CV_linear_SVM_%s_%s_null_model_fit_pvals.Rds",
#                                               feature_set, weighting_name)))) {
#     combo_wise_SVM_CV_weighting <- readRDS(paste0(rdata_path, 
#                                                     sprintf("Combo_wise_CV_linear_SVM_%s_%s.Rds",
#                                                             feature_set, weighting_name)))
#     
#     # Calculate p-values
#     pvalues <- calc_empirical_nulls(class_res = combo_wise_SVM_CV_weighting,
#                                     null_data = model_permutation_null_weighting,
#                                     grouping_var = "Brain_Region")
#     
#     saveRDS(pvalues, file=paste0(rdata_path, sprintf("Combo_wise_CV_linear_SVM_%s_%s_null_model_fit_pvals.Rds",
#                                                      feature_set, weighting_name)))
#   }
# }