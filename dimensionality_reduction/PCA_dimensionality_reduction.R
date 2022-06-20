github_dir <- "D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders/"
source(paste0(github_dir, "helper_functions/Linear_SVM.R"))
source(paste0(github_dir, "helper_functions/Visualization.R"))
source(paste0(github_dir, "helper_functions/Null_distributions.R"))
rdata_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/Rdata/"
library(penalizedSVM)
library(factoextra)
library(FactoMineR)

noise_proc = "AROMA+2P"
noise_label = "AROMA_2P"
univariate_feature_set = "catch22"

################################################################################
# Univariate PCA
################################################################################

# Load catch22 data for current noise processing method
catch22_feature_matrix <- readRDS(paste0(rdata_path, 
                                         sprintf("UCLA_%s_catch22_filtered_zscored.Rds", 
                                                 noise_label))) 

# ROI-wise and Feature-wise analyses are not prone to overfitting since
# the number of variables (22 or 82, respectively) is less than 
# the number of subjects (166)

# However, the combination of ROIs + Features --> 1,804 variables

# Generate vector of groups for PCA+SVM
group_vector <- catch22_feature_matrix %>%
  distinct(Subject_ID, group) %>%
  pull(group)

###### ROI-wise ###### 
# Run PCA
ROI_wise_PCA_list <- run_PCA_by_group_var(feature_matrix = catch22_feature_matrix,
                                          grouping_variable = "Brain_Region",
                                          feature_var = "names")

# Run linear SVM with increasing # PCs
weighting_name = "inv_prob"
if (!file.exists(paste0(rdata_path, sprintf("ROI_wise_PCA_CV_linear_SVM_%s_%s.Rds",
                                            univariate_feature_set, weighting_name)))) {
  ROI_wise_PCA_linear_SVM_res <- run_SVM_from_PCA(list_of_PCA_res = ROI_wise_PCA_list,
                                                  subject_dx_list = group_vector,
                                                  c_values = 1,
                                                  interval = 2,
                                                  use_inv_prob_weighting = TRUE,
                                                  use_SMOTE = FALSE,
                                                  return_all_fold_metrics = FALSE) 
  
  save(ROI_wise_PCA_linear_SVM_res, file = paste0(rdata_path, sprintf("ROI_wise_PCA_CV_linear_SVM_%s_%s.Rds",
                                                                      univariate_feature_set, weighting_name)))
} else {
  ROI_wise_PCA_linear_SVM_res <- readRDS(paste0(rdata_path, sprintf("ROI_wise_PCA_CV_linear_SVM_%s_%s.Rds",
                                                                    univariate_feature_set, weighting_name)))
}

###### Feature-wise ###### 

# Run PCA
Feature_wise_PCA_list <- run_PCA_by_group_var(feature_matrix = catch22_feature_matrix,
                                              grouping_variable = "names",
                                              feature_var = "Brain_Region")

# Run linear SVM with increasing # PCs
weighting_name = "inv_prob"
if (!file.exists(paste0(rdata_path, sprintf("Feature_wise_PCA_CV_linear_SVM_%s_%s.Rds",
                                            univariate_feature_set, weighting_name)))) {
  Feature_wise_PCA_linear_SVM_res <- run_SVM_from_PCA(list_of_PCA_res = Feature_wise_PCA_list,
                                                  subject_dx_list = group_vector,
                                                  c_values = 1,
                                                  interval = 2,
                                                  use_inv_prob_weighting = TRUE,
                                                  use_SMOTE = FALSE,
                                                  return_all_fold_metrics = FALSE) 
  
  save(Feature_wise_PCA_linear_SVM_res, file = paste0(rdata_path, sprintf("Feature_wise_PCA_CV_linear_SVM_%s_%s.Rds",
                                                                      univariate_feature_set, weighting_name)))
} else {
  Feature_wise_PCA_linear_SVM_res <- readRDS(paste0(rdata_path, sprintf("Feature_wise_PCA_CV_linear_SVM_%s_%s.Rds",
                                                                    univariate_feature_set, weighting_name)))
}

###### Combo-wise ###### 
Combo_wise_PCA_list <- run_PCA_by_group_var(feature_matrix = catch22_feature_matrix,
                                              grouping_variable = "Combo",
                                              feature_var = "Combo")

# Run linear SVM with increasing # PCs
weighting_name = "inv_prob"
if (!file.exists(paste0(rdata_path, sprintf("Combo_wise_PCA_CV_linear_SVM_%s_%s.Rds",
                                            univariate_feature_set, weighting_name)))) {
  Combo_wise_PCA_linear_SVM_res <- run_SVM_from_PCA(list_of_PCA_res = Combo_wise_PCA_list,
                                                      subject_dx_list = group_vector,
                                                      c_values = 1,
                                                      interval = 2,
                                                      use_inv_prob_weighting = TRUE,
                                                      use_SMOTE = FALSE,
                                                      return_all_fold_metrics = FALSE) 
  
  save(Combo_wise_PCA_linear_SVM_res, file = paste0(rdata_path, sprintf("Combo_wise_PCA_CV_linear_SVM_%s_%s.Rds",
                                                                          univariate_feature_set, weighting_name)))
} else {
  Combo_wise_PCA_linear_SVM_res <- readRDS(paste0(rdata_path, sprintf("Combo_wise_PCA_CV_linear_SVM_%s_%s.Rds",
                                                                        univariate_feature_set, weighting_name)))
}
