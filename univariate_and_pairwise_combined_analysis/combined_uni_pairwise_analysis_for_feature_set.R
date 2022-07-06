# Parse arguments
library(argparse)
library(factoextra)
library(FactoMineR)

parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--github_dir", default="/project/hctsa/annie/github/")
parser$add_argument("--rdata_path", default="/project/hctsa/annie/data/scz/UCLA/Rdata/")
parser$add_argument("--pydata_path", default="/project/hctsa/annie/data/scz/UCLA/pydata/")
parser$add_argument("--noise_proc", default="AROMA+2P+GMR")
parser$add_argument("--univariate_feature_set", default="catch22")
parser$add_argument("--pairwise_feature_set", default="pyspi_19")
# github_dir <- "D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders/"
# rdata_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/Rdata/"
# pydata_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/pydata/"
# noise_proc <- "AROMA+2P+GMR"
# univariate_feature_set <- "catch22"
# pairwise_feature_set <- "pyspi_19"

# Parse input arguments
args <- parser$parse_args()
github_dir <- args$github_dir
rdata_path <- args$rdata_path
pydata_path <- args$pydata_path
noise_proc <- args$noise_proc
univariate_feature_set <- args$univariate_feature_set
pairwise_feature_set <- args$pairwise_feature_set

### Source functions
# Main
source(paste0(github_dir, "helper_functions/Linear_SVM.R"))
source(paste0(github_dir, "helper_functions/Visualization.R"))
source(paste0(github_dir, "helper_functions/Null_distributions.R"))
source(paste0(github_dir, "helper_functions/PCA_functions.R"))

set.seed(127)
noise_label = gsub("\\+", "_", noise_proc)

# Use e1071 SVM with a linear kernel
test_package = "e1071"
kernel = "linear"

# ###############################################################################
# Load data
# ###############################################################################
univariate_data <- readRDS(paste0(rdata_path, sprintf("UCLA_%s_%s_filtered_zscored.Rds",
                                                      noise_label,
                                                      univariate_feature_set)))


pairwise_data <- readRDS(paste0(pydata_path, sprintf("UCLA_all_subject_%s_%s_filtered_zscored.Rds",
                                                     pairwise_feature_set, noise_label)))

SPI_directionality <- read.csv(paste0(github_dir, "pairwise_analysis/SPI_Direction_Info.csv"))

# Merge subjects
univariate_subjects <- univariate_data %>% ungroup() %>% distinct(Subject_ID, group)
pairwise_subjects <- pairwise_data %>% ungroup() %>% distinct(Subject_ID, group)
combined_subjects <- inner_join(univariate_subjects, pairwise_subjects)

univariate_data <- univariate_data %>% semi_join(combined_subjects)
pairwise_data <- pairwise_data %>% semi_join(combined_subjects)

if (!file.exists(paste0(rdata_path, "Filtered_subject_info_", 
                        univariate_feature_set, "_",
                        pairwise_feature_set, ".Rds"))) {
  saveRDS(combined_subjects, file=paste0(rdata_path, "Filtered_subject_info_", 
                                         univariate_feature_set, "_",
                                         pairwise_feature_set, ".Rds"))
}


################################################################################
# Generate model-free shuffle null distribution
################################################################################
if (!file.exists(paste0(rdata_path, sprintf("Null_Model_Free_Shuffles_combined_%s_%s.Rds",
                                            univariate_feature_set, pairwise_feature_set)))) {
  model_free_shuffle_null_res <- run_model_free_n_shuffles(num_shuffles = 100000,
                                                           feature_set = paste0(univariate_feature_set, "_",
                                                                                pairwise_feature_set),
                                                           rdata_path = rdata_path)
  saveRDS(model_free_shuffle_null_res, file = paste0(rdata_path, sprintf("Null_Model_Free_Shuffles_combined_%s_%s.Rds",
                                                                         univariate_feature_set, pairwise_feature_set)))
} else {
  model_free_shuffle_null_res <- readRDS(paste0(rdata_path, sprintf("Null_Model_Free_Shuffles_combined_%s_%s.Rds",
                                                                    univariate_feature_set, pairwise_feature_set)))
}


################################################################################
# Define weighting parameters
################################################################################
weighting_param_df <- data.frame(name = c("unweighted", "inv_prob", "SMOTE"),
                                 use_inv_prob_weighting = c(FALSE, TRUE, FALSE),
                                 use_SMOTE = c(FALSE, FALSE, TRUE))

################################################################################
# All catch22 features + pearson correlation
################################################################################

#### 10-fold linear SVM with different weights
# Iterate over weighting_param_df
for (i in 1:nrow(weighting_param_df)) {
  weighting_name <- weighting_param_df$name[i]
  use_inv_prob_weighting <- weighting_param_df$use_inv_prob_weighting[i]
  use_SMOTE <- weighting_param_df$use_SMOTE[i]
  
  # Run given weighting for 10-fold CV linear SVM
  if (!file.exists(paste0(rdata_path, sprintf("Univariate_%s_Pairwise_%s_CV_linear_SVM_%s.Rds",
                                              univariate_feature_set, pairwise_feature_set, weighting_name)))) {
    univariate_pairwise_SVM_CV_weighting <- run_combined_uni_pairwise_cv_svm_by_input_var(univariate_data = univariate_data,
                                                                                          univariate_feature_set = univariate_feature_set,
                                                                                          pairwise_data = pairwise_data,
                                                                                          pairwise_feature_set = pairwise_feature_set,
                                                                                          SPI_directionality = SPI_directionality,
                                                                                          svm_kernel = "linear",
                                                                                          test_package = "e1071",
                                                                                          noise_proc = "AROMA+2P+GMR",
                                                                                          return_all_fold_metrics = TRUE,
                                                                                          use_inv_prob_weighting = use_inv_prob_weighting,
                                                                                          use_SMOTE = use_SMOTE,
                                                                                          shuffle_labels = FALSE)
    saveRDS(univariate_pairwise_SVM_CV_weighting, file=paste0(rdata_path,
                                                              sprintf("Univariate_%s_Pairwise_%s_CV_linear_SVM_%s.Rds",
                                                                      univariate_feature_set, pairwise_feature_set, weighting_name)))
  }
}

#### Calculate p values from model-free shuffle null distribution
for (weighting_name in unique(weighting_param_df$name)) {
  if (!file.exists(paste0(rdata_path, sprintf("Univariate_%s_Pairwise_%s_CV_linear_SVM_%s_model_free_shuffle_pvals.Rds",
                                              univariate_feature_set, pairwise_feature_set, weighting_name)))) {
    univariate_pairwise_SVM_CV_weighting <- readRDS(paste0(rdata_path,
                                                           sprintf("Univariate_%s_Pairwise_%s_CV_linear_SVM_%s.Rds",
                                                                   univariate_feature_set, pairwise_feature_set, weighting_name)))

    combined_feature_set = paste0(univariate_feature_set, "_", pairwise_feature_set)
    # Calculate p-values
    pvalues <- calc_empirical_nulls(class_res = univariate_pairwise_SVM_CV_weighting,
                                    null_data = model_free_shuffle_null_res,
                                    feature_set = combined_feature_set,
                                    is_main_data_averaged = FALSE,
                                    grouping_var = "SPI")

    saveRDS(pvalues, file=paste0(rdata_path, sprintf("Univariate_%s_Pairwise_%s_CV_linear_SVM_%s_model_free_shuffle_pvals.Rds",
                                                     univariate_feature_set, pairwise_feature_set, weighting_name)))
  }
}

#### Generate empirical null model distributions per SPI
for (i in 1:nrow(weighting_param_df)) {
  weighting_name <- weighting_param_df$name[i]
  use_inv_prob_weighting <- weighting_param_df$use_inv_prob_weighting[i]
  use_SMOTE <- weighting_param_df$use_SMOTE[i]

  # Generate null-model fits distribution
  if (!file.exists(paste0(rdata_path, sprintf("Univariate_%s_Pairwise_%s_CV_linear_SVM_%s_model_permutation_null.Rds",
                                              univariate_feature_set, pairwise_feature_set, weighting_name)))) {
    model_permutation_null_weighting <- run_null_model_n_permutations_univariate_pairwise_combo(univariate_data = univariate_data,
                                                                                                univariate_feature_set = univariate_feature_set,
                                                                                                pairwise_data = pairwise_data,
                                                                                                pairwise_feature_set = pairwise_feature_set,
                                                                                                SPI_directionality = SPI_directionality,
                                                                                                svm_kernel = "linear",
                                                                                                test_package = "e1071",
                                                                                                noise_proc = "AROMA+2P+GMR",
                                                                                                num_permutations = 5,
                                                                                                return_all_fold_metrics = TRUE,
                                                                                                use_inv_prob_weighting = use_inv_prob_weighting,
                                                                                                use_SMOTE = use_SMOTE)

    saveRDS(model_permutation_null_weighting, file=paste0(rdata_path, sprintf("Univariate_%s_Pairwise_%s_CV_linear_SVM_%s_model_permutation_null.Rds",
                                                                              univariate_feature_set, pairwise_feature_set, weighting_name)))
  } else {
    model_permutation_null_weighting <- readRDS(paste0(rdata_path, sprintf("Univariate_%s_Pairwise_%s_CV_linear_SVM_%s_model_permutation_null.Rds",
                                                                           univariate_feature_set, pairwise_feature_set, weighting_name)))
  }

  # # Empirically derive p-values based on null model fits distribution
  # if (!file.exists(paste0(rdata_path, sprintf("pyspi_SPI_pairwise_CV_linear_SVM_null_model_fit_pvals_%s_%s.Rds",
  #                                             feature_set, weighting_name)))) {
  #   pyspi_SPI_pairwise_SVM_CV_weighting <- readRDS(paste0(rdata_path,
  #                                                  sprintf("pyspi_SPI_pairwise_CV_linear_SVM_%s_%s.Rds",
  #                                                          feature_set, weighting_name)))
  # 
  #   # Calculate p-values
  #   pvalues <- calc_empirical_nulls(class_res = pyspi_SPI_pairwise_SVM_CV_weighting,
  #                                   null_data = model_permutation_null_weighting,
  #                                   feature_set = feature_set,
  #                                   is_main_data_averaged = FALSE,
  #                                   grouping_var = "SPI")
  # 
  #   saveRDS(pvalues, file=paste0(rdata_path, sprintf("pyspi_SPI_pairwise_CV_linear_SVM_null_model_fit_pvals_%s_%s.Rds",
  #                                                    feature_set, weighting_name)))
  # }
}

#### PCA Analysis
if (!file.exists(paste0(rdata_path, sprintf("Univariate_%s_Pairwise_%s_PCA.Rds",
                                            univariate_feature_set, pairwise_feature_set)))) {
  combined_uni_pairwise_PCA_list <- run_PCA_for_uni_pairwise_combo(univariate_data = univariate_data,
                                                                   pairwise_data = pairwise_data)
  saveRDS(combined_uni_pairwise_PCA_list, file=paste0(rdata_path, sprintf("Univariate_%s_Pairwise_%s_PCA.Rds",
                                                                          univariate_feature_set, pairwise_feature_set)))
  
} else {
  combined_uni_pairwise_PCA_list <- readRDS(paste0(rdata_path, sprintf("Univariate_%s_Pairwise_%s_PCA.Rds",
                                                                       univariate_feature_set, pairwise_feature_set)))
}

# Run linear SVM with increasing # PCs
weighting_name = "inv_prob"
if (!file.exists(paste0(rdata_path, sprintf("Univariate_%s_Pairwise_%s_PCA_CV_linear_SVM_%s.Rds",
                                            univariate_feature_set, pairwise_feature_set,
                                            weighting_name)))) {
  combined_uni_pairwise_PCA_linear_SVM_res <- run_SVM_from_PCA(list_of_PCA_res = combined_uni_pairwise_PCA_list,
                                                               subject_dx_list = combined_subjects$group,
                                                               c_values = 1,
                                                               interval = 2,
                                                               use_inv_prob_weighting = TRUE,
                                                               use_SMOTE = FALSE,
                                                               return_all_fold_metrics = FALSE) 
  
  save(combined_uni_pairwise_PCA_linear_SVM_res, file = paste0(rdata_path, sprintf("Univariate_%s_Pairwise_%s_PCA_CV_linear_SVM_%s.Rds",
                                                                      univariate_feature_set, pairwise_feature_set,
                                                                      weighting_name)))
}