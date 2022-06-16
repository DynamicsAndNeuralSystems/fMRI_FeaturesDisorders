# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--github_dir", default="/project/hctsa/annie/github/")
parser$add_argument("--rdata_path", default="/project/hctsa/annie/data/scz/UCLA/Rdata/")
parser$add_argument("--feature_set", default="catch22")
# github_dir <- "D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders/"
# rdata_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/Rdata/"
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

# Compare all three noise processing methods
noise_procs = c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER")

# Use e1071 SVM with a linear kernel
test_package = "e1071"
kernel = "linear"


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
# Brain Region-wise analysis
################################################################################

#### 10-fold linear SVM with different weights
# Iterate over weighting_param_df 
for (i in 1:nrow(weighting_param_df)) {
  weighting_name <- weighting_param_df$name[i]
  use_inv_prob_weighting <- weighting_param_df$use_inv_prob_weighting[i]
  use_SMOTE <- weighting_param_df$use_SMOTE[i]
  
  # Run given weighting for 10-fold CV linear SVM
  if (!file.exists(paste0(rdata_path, sprintf("ROI_wise_CV_linear_SVM_%s_%s.Rds",
                                              feature_set, weighting_name)))) {
    region_wise_SVM_CV_weighting <- run_univariate_cv_svm_by_input_var(rdata_path = rdata_path,
                                                             feature_set = feature_set,
                                                             test_package = test_package,
                                                             svm_kernel = kernel,
                                                             grouping_var = "Brain_Region",
                                                             svm_feature_var = "Feature",
                                                             use_inv_prob_weighting = use_inv_prob_weighting,
                                                             use_SMOTE = use_SMOTE,
                                                             noise_procs = noise_procs)
    saveRDS(region_wise_SVM_CV_weighting, file=paste0(rdata_path, 
                                                       sprintf("ROI_wise_CV_linear_SVM_%s_%s.Rds",
                                                               feature_set, weighting_name)))
  }
}

#### Calculate p values from model-free shuffle null distribution
for (weighting_name in unique(weighting_param_df$name)) {
  if (!file.exists(paste0(rdata_path, sprintf("ROI_wise_CV_linear_SVM_%s_%s_model_free_shuffle_pvals.Rds",
                                              feature_set, weighting_name)))) {
    region_wise_SVM_CV_weighting <- readRDS(paste0(rdata_path, 
                                                   sprintf("ROI_wise_CV_linear_SVM_%s_%s.Rds",
                                                           feature_set, weighting_name)))
    
    # Calculate p-values
    pvalues <- calc_empirical_nulls(class_res = region_wise_SVM_CV_weighting,
                                    null_data = model_free_shuffle_null_res,
                                    grouping_var = "Brain_Region")
    
    saveRDS(pvalues, file=paste0(rdata_path, sprintf("ROI_wise_CV_linear_SVM_%s_%s_model_free_shuffle_pvals.Rds",
                                                     feature_set, weighting_name)))
  }
}

#### Generate empirical null model distributions per brain region
for (i in 1:nrow(weighting_param_df)) {
  weighting_name <- weighting_param_df$name[i]
  use_inv_prob_weighting <- weighting_param_df$use_inv_prob_weighting[i]
  use_SMOTE <- weighting_param_df$use_SMOTE[i]
  
  # Generate null-model fits distribution
  if (!file.exists(paste0(rdata_path, sprintf("ROI_wise_model_permutation_null_%s_%s.Rds",
                                              feature_set, weighting_name)))) {
    model_permutation_null_weighting <- run_null_model_n_permutations(rdata_path,
                                                                       feature_set = feature_set,
                                                                       noise_procs = noise_procs,
                                                                       grouping_var = "Brain_Region",
                                                                       svm_feature_var = "Feature",
                                                                       num_permutations = 10,
                                                                       use_inv_prob_weighting = use_inv_prob_weighting,
                                                                       use_SMOTE = use_SMOTE)
    
    saveRDS(model_permutation_null_weighting, file=paste0(rdata_path, sprintf("ROI_wise_model_permutation_null_%s_%s.Rds",
                                                                               feature_set, weighting_name)))
  }
  
  # Empirically derive p-values based on null model fits distribution
  if (!file.exists(paste0(rdata_path, sprintf("ROI_wise_CV_linear_SVM_%s_%s_null_model_fit_pvals.Rds",
                                              feature_set, weighting_name)))) {
    region_wise_SVM_CV_weighting <- readRDS(paste0(rdata_path, 
                                                   sprintf("ROI_wise_CV_linear_SVM_%s_%s.Rds",
                                                           feature_set, weighting_name)))
    
    # Calculate p-values
    pvalues <- calc_empirical_nulls(class_res = region_wise_SVM_CV_weighting,
                                    null_data = model_permutation_null_weighting,
                                    grouping_var = "Brain_Region")
    
    saveRDS(pvalues, file=paste0(rdata_path, sprintf("ROI_wise_CV_linear_SVM_%s_%s_null_model_fit_pvals.Rds",
                                                     feature_set, weighting_name)))
  }
}


################################################################################
# TS Feature-wise analysis
################################################################################

#### 10-fold linear SVM with different weights
# Iterate over weighting_param_df 
for (i in 1:nrow(weighting_param_df)) {
  weighting_name <- weighting_param_df$name[i]
  use_inv_prob_weighting <- weighting_param_df$use_inv_prob_weighting[i]
  use_SMOTE <- weighting_param_df$use_SMOTE[i]
  
  # Run given weighting for 10-fold CV linear SVM
  if (!file.exists(paste0(rdata_path, sprintf("Feature_wise_CV_linear_SVM_%s_%s.Rds",
                                              feature_set, weighting_name)))) {
    feature_wise_SVM_CV_weighting <- run_univariate_cv_svm_by_input_var(rdata_path = rdata_path,
                                                            feature_set = feature_set,
                                                            test_package = test_package,
                                                            svm_kernel = kernel,
                                                            grouping_var = "Feature",
                                                            svm_feature_var = "Brain_Region",
                                                            use_inv_prob_weighting = use_inv_prob_weighting,
                                                            use_SMOTE = use_SMOTE,
                                                            noise_procs = noise_procs)
    saveRDS(feature_wise_SVM_CV_weighting, file=paste0(rdata_path, 
                                                      sprintf("Feature_wise_CV_linear_SVM_%s_%s.Rds",
                                                              feature_set, weighting_name)))
  }
}

#### Calculate p values from model-free shuffle null distribution
for (weighting_name in unique(weighting_param_df$name)) {
  if (!file.exists(paste0(rdata_path, sprintf("Feature_wise_CV_linear_SVM_%s_%s_model_free_shuffle_pvals.Rds",
                                              feature_set, weighting_name)))) {
    feature_wise_SVM_CV_weighting <- readRDS(paste0(rdata_path, 
                                                   sprintf("Feature_wise_CV_linear_SVM_%s_%s.Rds",
                                                           feature_set, weighting_name)))
    
    # Calculate p-values
    pvalues <- calc_empirical_nulls(class_res = feature_wise_SVM_CV_weighting,
                                    null_data = model_free_shuffle_null_res,
                                    grouping_var = "Brain_Region")
    
    saveRDS(pvalues, file=paste0(rdata_path, sprintf("Feature_wise_CV_linear_SVM_%s_%s_model_free_shuffle_pvals.Rds",
                                                     feature_set, weighting_name)))
  }
}

#### Generate empirical null model distributions per brain region
for (i in 1:nrow(weighting_param_df)) {
  weighting_name <- weighting_param_df$name[i]
  use_inv_prob_weighting <- weighting_param_df$use_inv_prob_weighting[i]
  use_SMOTE <- weighting_param_df$use_SMOTE[i]
  
  # Generate null-model fits distribution
  if (!file.exists(paste0(rdata_path, sprintf("Feature_wise_model_permutation_null_%s_%s.Rds",
                                              feature_set, weighting_name)))) {
    model_permutation_null_weighting <- run_null_model_n_permutations(rdata_path,
                                                                      feature_set = feature_set,
                                                                      noise_procs = noise_procs,
                                                                      grouping_var = "Feature",
                                                                      svm_feature_var = "Brain_Region",
                                                                      num_permutations = 40,
                                                                      use_inv_prob_weighting = use_inv_prob_weighting,
                                                                      use_SMOTE = use_SMOTE)
    
    saveRDS(model_permutation_null_weighting, file=paste0(rdata_path, sprintf("Feature_wise_model_permutation_null_%s_%s.Rds",
                                                                              feature_set, weighting_name)))
  }
  
  # Empirically derive p-values based on null model fits distribution
  if (!file.exists(paste0(rdata_path, sprintf("Feature_wise_CV_linear_SVM_%s_%s_null_model_fit_pvals.Rds",
                                              feature_set, weighting_name)))) {
    feature_wise_SVM_CV_weighting <- readRDS(paste0(rdata_path, 
                                                   sprintf("Feature_wise_CV_linear_SVM_%s_%s.Rds",
                                                           feature_set, weighting_name)))
    
    # Calculate p-values
    pvalues <- calc_empirical_nulls(class_res = feature_wise_SVM_CV_weighting,
                                    null_data = model_permutation_null_weighting,
                                    grouping_var = "Brain_Region")
    
    saveRDS(pvalues, file=paste0(rdata_path, sprintf("Feature_wise_CV_linear_SVM_%s_%s_null_model_fit_pvals.Rds",
                                                     feature_set, weighting_name)))
  }
}

################################################################################
# ROI/Feature Combo-wise analysis
################################################################################

#### 10-fold linear SVM with different weights
# Iterate over weighting_param_df 
for (i in 1:nrow(weighting_param_df)) {
  weighting_name <- weighting_param_df$name[i]
  use_inv_prob_weighting <- weighting_param_df$use_inv_prob_weighting[i]
  use_SMOTE <- weighting_param_df$use_SMOTE[i]
  
  # Run given weighting for 10-fold CV linear SVM
  if (!file.exists(paste0(rdata_path, sprintf("Combo_wise_CV_linear_SVM_%s_%s.Rds",
                                              feature_set, weighting_name)))) {
    combo_wise_SVM_CV_weighting <- run_univariate_cv_svm_by_input_var(rdata_path = rdata_path,
                                                             feature_set = feature_set,
                                                             test_package = test_package,
                                                             svm_kernel = kernel,
                                                             grouping_var = "Combo",
                                                             svm_feature_var = "Combo",
                                                             use_inv_prob_weighting = use_inv_prob_weighting,
                                                             use_SMOTE = use_SMOTE,
                                                             noise_procs = noise_procs)
    saveRDS(combo_wise_SVM_CV_weighting, file=paste0(rdata_path, 
                                                       sprintf("Combo_wise_CV_linear_SVM_%s_%s.Rds",
                                                               feature_set, weighting_name)))
  }
}

#### Calculate p values from model-free shuffle null distribution
for (weighting_name in unique(weighting_param_df$name)) {
  if (!file.exists(paste0(rdata_path, sprintf("Combo_wise_CV_linear_SVM_%s_%s_model_free_shuffle_pvals.Rds",
                                              feature_set, weighting_name)))) {
    combo_wise_SVM_CV_weighting <- readRDS(paste0(rdata_path, 
                                                    sprintf("Combo_wise_CV_linear_SVM_%s_%s.Rds",
                                                            feature_set, weighting_name)))
    
    # Calculate p-values
    pvalues <- calc_empirical_nulls(class_res = combo_wise_SVM_CV_weighting,
                                    null_data = model_free_shuffle_null_res,
                                    grouping_var = "Brain_Region")
    
    saveRDS(pvalues, file=paste0(rdata_path, sprintf("Combo_wise_CV_linear_SVM_%s_%s_model_free_shuffle_pvals.Rds",
                                                     feature_set, weighting_name)))
  }
}

#### Generate empirical null model distributions per brain region
for (i in 1:nrow(weighting_param_df)) {
  weighting_name <- weighting_param_df$name[i]
  use_inv_prob_weighting <- weighting_param_df$use_inv_prob_weighting[i]
  use_SMOTE <- weighting_param_df$use_SMOTE[i]
  
  # Generate null-model fits distribution
  if (!file.exists(paste0(rdata_path, sprintf("Combo_wise_model_permutation_null_%s_%s.Rds",
                                              feature_set, weighting_name)))) {
    model_permutation_null_weighting <- run_null_model_n_permutations(rdata_path,
                                                                      feature_set = feature_set,
                                                                      noise_procs = noise_procs,
                                                                      grouping_var = "Combo",
                                                                      svm_feature_var = "Combo",
                                                                      num_permutations = 100,
                                                                      use_inv_prob_weighting = use_inv_prob_weighting,
                                                                      use_SMOTE = use_SMOTE)
    
    saveRDS(model_permutation_null_weighting, file=paste0(rdata_path, sprintf("Combo_wise_model_permutation_null_%s_%s.Rds",
                                                                              feature_set, weighting_name)))
  }
  
  # Empirically derive p-values based on null model fits distribution
  if (!file.exists(paste0(rdata_path, sprintf("Combo_wise_CV_linear_SVM_%s_%s_null_model_fit_pvals.Rds",
                                              feature_set, weighting_name)))) {
    combo_wise_SVM_CV_weighting <- readRDS(paste0(rdata_path, 
                                                    sprintf("Combo_wise_CV_linear_SVM_%s_%s.Rds",
                                                            feature_set, weighting_name)))
    
    # Calculate p-values
    pvalues <- calc_empirical_nulls(class_res = combo_wise_SVM_CV_weighting,
                                    null_data = model_permutation_null_weighting,
                                    grouping_var = "Brain_Region")
    
    saveRDS(pvalues, file=paste0(rdata_path, sprintf("Combo_wise_CV_linear_SVM_%s_%s_null_model_fit_pvals.Rds",
                                                     feature_set, weighting_name)))
  }
}