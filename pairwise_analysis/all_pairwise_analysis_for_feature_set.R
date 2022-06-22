# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--project_path", default="/project/hctsa/annie/")
parser$add_argument("--github_dir", default="/project/hctsa/annie/github/")
parser$add_argument("--rdata_path", default="/project/hctsa/annie/data/scz/UCLA/Rdata/")
parser$add_argument("--pydata_path", default="/project/hctsa/annie/data/scz/UCLA/pydata/")
parser$add_argument("--feature_set", default="pyspi_19")
# github_dir <- "D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders/"
# rdata_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/Rdata/"
# pydata_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/pydata/"
# feature_set <- "pyspi_19"

# Parse input arguments
args <- parser$parse_args()
project_path <- args$project_path
github_dir <- args$github_dir
rdata_path <- args$rdata_path
pydata_path <- args$pydata_path
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

# # ###############################################################################
# # Load pyspi data
# # ###############################################################################
# pyspi_data_file <- paste0(pydata_path, "UCLA_all_subject_pyspi_19_AROMA_2P_GMR_filtered_zscored.Rds")
# pyspi_data <- readRDS(pyspi_data_file) %>%
#   mutate(group = stringr::str_to_sentence(group))
# 
# SPI_directionality <- read.csv(paste0(github_dir, "pairwise_analysis/SPI_Direction_Info.csv"))
# 
# ################################################################################
# # Generate model-free shuffle null distribution
# ################################################################################
# if (!file.exists(paste0(rdata_path, sprintf("Null_Model_Free_Shuffles_%s.Rds",
#                                             feature_set)))) {
#   model_free_shuffle_null_res <- run_model_free_n_shuffles(num_shuffles = 100000,
#                                                            feature_set = feature_set,
#                                                            rdata_path = rdata_path)
#   saveRDS(model_free_shuffle_null_res, file = paste0(rdata_path, sprintf("Null_Model_Free_Shuffles_%s.Rds",
#                                                                          feature_set)))
# } else {
#   model_free_shuffle_null_res <- readRDS(paste0(rdata_path, sprintf("Null_Model_Free_Shuffles_%s.Rds",
#                                                                     feature_set)))
# }
# 
# 
# ################################################################################
# # Define weighting parameters
# ################################################################################
# 
# weighting_param_df <- data.frame(name = c("inv_prob"),
#                                  use_inv_prob_weighting = c(TRUE),
#                                  use_SMOTE = c(FALSE))
# 
# ################################################################################
# # Per ROI pair
# ################################################################################
# 
# 
# 
# ################################################################################
# # Per-feature, all ROI pair combinations
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
#   if (!file.exists(paste0(rdata_path, sprintf("pyspi_SPI_pairwise_CV_linear_SVM_%s_%s.Rds",
#                                               feature_set, weighting_name)))) {
#     pyspi_SPI_pairwise_SVM_CV_weighting <- run_pairwise_cv_svm_by_input_var(pairwise_data = pyspi_data,
#                                                                             SPI_directionality = SPI_directionality,
#                                                                             svm_kernel = "linear",
#                                                                             grouping_var = "SPI",
#                                                                             svm_feature_var = "region_pair",
#                                                                             test_package = "e1071",
#                                                                             noise_proc = "AROMA+2P+GMR",
#                                                                             return_all_fold_metrics = TRUE,
#                                                                             use_inv_prob_weighting = use_inv_prob_weighting,
#                                                                             use_SMOTE = use_SMOTE,
#                                                                             shuffle_labels = FALSE)
#     saveRDS(pyspi_SPI_pairwise_SVM_CV_weighting, file=paste0(rdata_path,
#                                                            sprintf("pyspi_SPI_pairwise_CV_linear_SVM_%s_%s.Rds",
#                                                                    feature_set, weighting_name)))
#   }
# }
# 
# #### Calculate p values from model-free shuffle null distribution
# for (weighting_name in unique(weighting_param_df$name)) {
#   if (!file.exists(paste0(rdata_path, sprintf("pyspi_SPI_pairwise_CV_linear_SVM_model_free_shuffle_pvals_%s_%s.Rds",
#                                               feature_set, weighting_name)))) {
#     pyspi_SPI_pairwise_SVM_CV_weighting <- readRDS(paste0(rdata_path,
#                                                    sprintf("pyspi_SPI_pairwise_CV_linear_SVM_%s_%s.Rds",
#                                                            feature_set, weighting_name)))
# 
#     # Calculate p-values
#     pvalues <- calc_empirical_nulls(class_res = pyspi_SPI_pairwise_SVM_CV_weighting,
#                                     null_data = model_free_shuffle_null_res,
#                                     feature_set = feature_set,
#                                     is_data_averaged = FALSE,
#                                     grouping_var = "SPI")
# 
#     saveRDS(pvalues, file=paste0(rdata_path, sprintf("pyspi_SPI_pairwise_CV_linear_SVM_model_free_shuffle_pvals_%s_%s.Rds",
#                                                      feature_set, weighting_name)))
#   }
# }

# template file
num_permutations <- 5
template_pbs_file <- paste0(github_dir, "pairwise_analysis/template_null_model_fit.pbs")

lookup_list <- list("PROJECT_NAME" = "hctsa", 
                    "NAME" = "pyspi_SPIwise_null_model_fit",
                    "MEMNUM" = "20",
                    "NCPUS" = "1",
                    "GITHUB_DIR" = github_dir,
                    "PROJECT_DIR" = project_path,
                    "EMAIL" = "abry4213@uni.sydney.edu.au",
                    "PBS_NOTIFY" = "abe",
                    "WALL_HRS" = "1",
                    "PAIRWISE_DATA_FILE" = paste0(pydata_path, "UCLA_all_subject_pyspi_19_AROMA_2P_GMR_filtered_zscored.Rds"),
                    "SPI_DIRECTIONALITY_FILE" = paste0(github_dir, "pairwise_analysis/SPI_Direction_Info.csv"),
                    "FEATURE_SET" = "pyspi_19",
                    "GROUPING_VAR" = "SPI",
                    "SVM_FEATURE_VAR" = "region_pair",
                    "NOISE_PROC" = noise_proc)

to_be_replaced <- names(lookup_list)
replacement_values <- unlist(unname(lookup_list))

#### Generate empirical null model distributions per SPI
for (i in 1:nrow(weighting_param_df)) {
  weighting_name <- weighting_param_df$name[i]
  use_inv_prob_weighting <- weighting_param_df$use_inv_prob_weighting[i]
  use_SMOTE <- weighting_param_df$use_SMOTE[i]
  
  output_data_dir <- paste0(rdata_path, sprintf("Pairwise_%s_%s_null_model_fits/",
                                           weighting_name, feature_set))
  
  output_scripts_dir <- paste0(github_dir, sprintf("pairwise_analysis/Pairwise_%s_%s_null_model_fits/",
                                                   weighting_name, feature_set))
  icesTAF::mkdir(output_data_dir)
  icesTAF::mkdir(output_scripts_dir)

  for (j in 1:num_permutations) {
    new_pbs_file <- readLines(template_pbs_file)
    
    # Replace file paths
    pbs_text_replaced <- mgsub::mgsub(new_pbs_file,
                                      to_be_replaced,
                                      replacement_values)
    
    # Replace null iteration number
    pbs_text_replaced <- gsub("iteri", j, pbs_text_replaced)
    
    # Write updated PBS script to file
    output_pbs_file <- writeLines(pbs_text_replaced, 
                                  paste0(output_scripts_dir, 
                                         "null_iter_", j, ".pbs"))
  }
}

# ################################################################################
# # ROI Pair/SPI Combo-wise analysis
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
#   if (!file.exists(paste0(rdata_path, sprintf("pyspi_Combo_pairwise_CV_linear_SVM_%s_%s.Rds",
#                                               feature_set, weighting_name)))) {
#     pyspi_combo_pairwise_SVM_CV_weighting <- run_pairwise_cv_svm_by_input_var(pairwise_data = pyspi_data,
#                                                                               SPI_directionality = SPI_directionality,
#                                                                               svm_kernel = "linear",
#                                                                               grouping_var = "Combo",
#                                                                               svm_feature_var = "Combo",
#                                                                               test_package = "e1071",
#                                                                               noise_proc = "AROMA+2P+GMR",
#                                                                               return_all_fold_metrics = TRUE,
#                                                                               use_inv_prob_weighting = use_inv_prob_weighting,
#                                                                               use_SMOTE = use_SMOTE,
#                                                                               shuffle_labels = FALSE)
#     saveRDS(pyspi_combo_pairwise_SVM_CV_weighting, file=paste0(rdata_path,
#                                                        sprintf("pyspi_Combo_pairwise_CV_linear_SVM_%s_%s.Rds",
#                                                                feature_set, weighting_name)))
#   }
# }
# 
# #### Calculate p values from model-free shuffle null distribution
# for (weighting_name in unique(weighting_param_df$name)) {
#   if (!file.exists(paste0(rdata_path, sprintf("pyspi_Combo_pairwise_CV_linear_SVM_model_free_shuffle_pvals_%s_%s.Rds",
#                                               feature_set, weighting_name)))) {
#     pyspi_combo_pairwise_SVM_CV_weighting <- readRDS(paste0(rdata_path,
#                                                             sprintf("pyspi_Combo_pairwise_CV_linear_SVM_%s_%s.Rds",
#                                                                     feature_set, weighting_name)))
# 
#     # Calculate p-values
#     pvalues <- calc_empirical_nulls(class_res = pyspi_combo_pairwise_SVM_CV_weighting,
#                                     null_data = model_free_shuffle_null_res,
#                                     feature_set = feature_set,
#                                     is_data_averaged = FALSE,
#                                     grouping_var = "Combo")
# 
#     saveRDS(pvalues, file=paste0(rdata_path, sprintf("pyspi_Combo_pairwise_CV_linear_SVM_model_free_shuffle_pvals_%s_%s.Rds",
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
#   if (!file.exists(paste0(rdata_path, sprintf("pyspi_Combo_pairwise_model_permutation_null_%s_%s_500perm.Rds",
#                                               feature_set, weighting_name)))) {
#     model_permutation_null_weighting <- run_null_model_n_permutations_pairwise(pairwise_data = pyspi_data,
#                                                                                noise_proc = "AROMA+2P+GMR",
#                                                                                feature_set = feature_set,
#                                                                                test_package = "e1071",
#                                                                                svm_kernel = "linear",
#                                                                                grouping_var = "Combo",
#                                                                                svm_feature_var = "Combo",
#                                                                                return_all_fold_metrics = FALSE,
#                                                                                SPI_directionality = SPI_directionality,
#                                                                                num_permutations = 500,
#                                                                                use_inv_prob_weighting = use_inv_prob_weighting,
#                                                                                use_SMOTE = use_SMOTE)
# 
#     saveRDS(model_permutation_null_weighting, file=paste0(rdata_path, sprintf("pyspi_Combo_pairwise_model_permutation_null_%s_%s_500perm.Rds",
#                                                                               feature_set, weighting_name)))
#   } else {
#     model_permutation_null_weighting <- readRDS(paste0(rdata_path, sprintf("pyspi_Combo_pairwise_model_permutation_null_%s_%s_500perm.Rds",
#                                                                            feature_set, weighting_name)))
#   }
# 
#   # Empirically derive p-values based on null model fits distribution
#   if (!file.exists(paste0(rdata_path, sprintf("pyspi_Combo_pairwise_CV_linear_SVM_model_permutation_null_pvals_%s_%s_500perm.Rds",
#                                               feature_set, weighting_name)))) {
#     combo_wise_SVM_CV_weighting <- readRDS(paste0(rdata_path,
#                                                     sprintf("pyspi_Combo_pairwise_CV_linear_SVM_%s_%s.Rds",
#                                                             feature_set, weighting_name)))
# 
#     # Calculate p-values
#     pvalues <- calc_empirical_nulls(class_res = combo_wise_SVM_CV_weighting,
#                                     null_data = model_permutation_null_weighting,
#                                     feature_set = feature_set,
#                                     is_data_averaged = FALSE,
#                                     grouping_var = "Combo")
# 
#     saveRDS(pvalues, file=paste0(rdata_path, sprintf("pyspi_Combo_pairwise_CV_linear_SVM_model_permutation_null_pvals_%s_%s_500perm.Rds",
#                                                      feature_set, weighting_name)))
#   }
# }