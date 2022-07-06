# Parse arguments
library(argparse)
library(factoextra)
library(FactoMineR)

parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--project_path", default="/project/hctsa/annie/")
parser$add_argument("--github_dir", default="/project/hctsa/annie/github/")
parser$add_argument("--rdata_path", default="/project/hctsa/annie/data/scz/UCLA/Rdata/")
parser$add_argument("--pydata_path", default="/project/hctsa/annie/data/scz/UCLA/pydata/")
parser$add_argument("--noise_proc", default="AROMA+2P+GMR")
parser$add_argument("--univariate_feature_set", default="catch22")
parser$add_argument("--pairwise_feature_set", default="pyspi_19")
# project_path <- ""
# github_dir <- "D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders/"
# rdata_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/Rdata/"
# pydata_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/pydata/"
# noise_proc <- "AROMA+2P+GMR"
# univariate_feature_set <- "catch22"
# pairwise_feature_set <- "pyspi_19"

# Parse input arguments
args <- parser$parse_args()
project_path <- args$project_path
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
# weighting_param_df <- data.frame(name = c("unweighted", "inv_prob", "SMOTE"),
#                                  use_inv_prob_weighting = c(FALSE, TRUE, FALSE),
#                                  use_SMOTE = c(FALSE, FALSE, TRUE))
weighting_param_df <- data.frame(name = c("inv_prob"),
                                 use_inv_prob_weighting = c(TRUE),
                                 use_SMOTE = c(FALSE))

grouping_df <- data.frame(grouping_var = "SPI",
                          SVM_feature_var = "region_pair")

################################################################################
# All catch22 features with each SPI individually
################################################################################

#### 10-fold linear SVM with different weights
# Iterate over weighting_param_df
for (i in 1:nrow(weighting_param_df)) {
  weighting_name <- weighting_param_df$name[i]
  use_inv_prob_weighting <- weighting_param_df$use_inv_prob_weighting[i]
  use_SMOTE <- weighting_param_df$use_SMOTE[i]
  
  grouping_var <- grouping_df$grouping_var
  SVM_feature_var <- grouping_df$SVM_feature_var
  
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
  
  #### Calculate p values from model-free shuffle null distribution
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
  
  #### Generate empirical null model distributions per SPI
  weighting_null_dist_file <- paste0(rdata_path, sprintf("Univariate_%s_Pairwise_%s_CV_linear_SVM_%s_model_permutation_null.Rds",
                                                         univariate_feature_set, pairwise_feature_set, weighting_name))
  

  # Generate null-model fits distribution
  # if (!file.exists(weighting_null_dist_file)) {
    # Output script dir
    output_data_dir <- paste0(rdata_path, sprintf("univariate_and_SPI_pairwise_%s_null_model_fits/",
                                                  weighting_name))
    output_scripts_dir <- paste0(github_dir, sprintf("univariate_and_pairwise_combined_analysis/univariate_and_SPI_wise_%s_null_model_fits/",
                                                     weighting_name))
    
    # save preprocessed univariate and pairwise data to files
    univariate_data_file = paste0(rdata_path, "univariate_data_for_combined_uni_pairwise.Rds")
    if (!file.exists(univariate_data_file)) {
      saveRDS(univariate_data, file = univariate_data_file)
    }
    
    pairwise_data_file = paste0(rdata_path, "pairwise_data_for_combined_uni_pairwise.Rds")
    if (!file.exists(pairwise_data_file)) {
      saveRDS(pairwise_data, file = pairwise_data_file)
    }
    
    icesTAF::mkdir(output_data_dir)
    icesTAF::mkdir(output_scripts_dir)
    
    num_permutations <- 180
    nperm_per_iter <- 3
    num_k_folds <- 10
    template_pbs_file <- paste0(github_dir, "univariate_and_pairwise_combined_analysis/template_null_model_fit.pbs")
    
    lookup_list <- list("PROJECT_NAME" = "hctsa",
                        "NAME" = "univariate_pairwise_combined_SPI_wise_null_model_fit",
                        "MEMNUM" = "20",
                        "NCPUS" = "1",
                        "GITHUB_DIR" = github_dir,
                        "PROJECT_DIR" = project_path,
                        "EMAIL" = "abry4213@uni.sydney.edu.au",
                        "PBS_NOTIFY" = "abe",
                        "WALL_HRS" = "8",
                        "UNIVARIATE_DATA_FILE" = univariate_data_file,
                        "PAIRWISE_DATA_FILE" = pairwise_data_file,
                        "SPI_DIRECTIONALITY_FILE" = paste0(github_dir, "pairwise_analysis/SPI_Direction_Info.csv"),
                        "NUM_K_FOLDS" = num_k_folds,
                        "NUM_PERMS_PER_ITER" = nperm_per_iter,
                        "OUTPUT_DATA_DIR" = output_data_dir,
                        "FEATURE_SET" = "pyspi_19",
                        "GROUPING_VAR" = grouping_var,
                        "SVM_FEATURE_VAR" = SVM_feature_var,
                        "NOISE_PROC" = noise_proc)
    
    to_be_replaced <- names(lookup_list)
    replacement_values <- unlist(unname(lookup_list))
    
    for (p in 1:num_permutations) {
      
      
      # Run command if null file doesn't exist
      if (!file.exists(sprintf("%s/Univariate_%s_Pairwise_%s_CV_linear_SVM_%s_null_model_fit_iter_%s.Rds",
                               output_data_dir, univariate_feature_set, 
                               pairwise_feature_set, weighting_name, p))) {
        cat("\nNow running null perms for iteration", p, "\n")
        new_pbs_file <- readLines(template_pbs_file)
        
        # Replace file paths
        pbs_text_replaced <- mgsub::mgsub(new_pbs_file,
                                          to_be_replaced,
                                          replacement_values)
        
        # Replace null iteration number
        pbs_text_replaced <- gsub("iterj", p, pbs_text_replaced)
        
        # Write updated PBS script to file
        output_pbs_file <- writeLines(pbs_text_replaced,
                                      paste0(output_scripts_dir,
                                             "null_iter_", p, ".pbs"))
        
        system(paste0("qsub ", output_scripts_dir, "null_iter_", p, ".pbs"))
        
      }
    }
    
    # saveRDS(model_permutation_null_weighting, file=paste0(rdata_path, sprintf("Univariate_%s_Pairwise_%s_CV_linear_SVM_%s_model_permutation_null.Rds",
    #                                                                           univariate_feature_set, pairwise_feature_set, weighting_name)))
  # } else {
  #   model_permutation_null_weighting <- readRDS(paste0(rdata_path, sprintf("Univariate_%s_Pairwise_%s_CV_linear_SVM_%s_model_permutation_null.Rds",
  #                                                                          univariate_feature_set, pairwise_feature_set, weighting_name)))
  # }

  # # Empirically derive p-values based on null model fits distribution
  # if (!file.exists(paste0(rdata_path, sprintf("Univariate_%s_Pairwise_%s_CV_linear_SVM_%s_null_model_fit_pvals.Rds",
  #                                             univariate_feature_set, pairwise_feature_set, weighting_name)))) {
  #   univariate_pairwise_SVM_CV_weighting <- readRDS(paste0(rdata_path,
  #                                                          sprintf("Univariate_%s_Pairwise_%s_CV_linear_SVM_%s.Rds",
  #                                                                  univariate_feature_set, pairwise_feature_set, weighting_name)))
  # 
  #   # Calculate p-values
  #   pvalues <- calc_empirical_nulls(class_res = univariate_pairwise_SVM_CV_weighting,
  #                                   null_data = model_permutation_null_weighting,
  #                                   feature_set = feature_set,
  #                                   is_main_data_averaged = FALSE,
  #                                   grouping_var = "SPI")
  # 
  #   saveRDS(pvalues, file=paste0(rdata_path, sprintf("Univariate_%s_Pairwise_%s_CV_linear_SVM_%s_null_model_fit_pvals.Rds",
  #                                                    univariate_feature_set, pairwise_feature_set, weighting_name)))
  # }
  
}

################################################################################
# PCA dimensionality reduction
################################################################################

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