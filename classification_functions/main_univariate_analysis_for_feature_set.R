# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--project_path", default="/project/hctsa/annie/")
parser$add_argument("--github_dir", default="/project/hctsa/annie/github/")
parser$add_argument("--rdata_path", default="/project/hctsa/annie/data/scz/UCLA/Rdata/")
parser$add_argument("--feature_set", default="catch22")
# project_path <- "D:/Virtual_Machines/Shared_Folder/github/"
# github_dir <- "D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders/"
# rdata_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/Rdata/"
# feature_set <- "catch22"

# Parse input arguments
args <- parser$parse_args()
project_path <- args$project_path
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
# Create ten folds to use for all analyses
################################################################################

subjects_to_use <- readRDS(paste0(rdata_path, "UCLA_Subjects_with_Univariate_and_Pairwise.Rds"))

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

# weighting_param_df <- data.frame(name = c("unweighted", "inv_prob"),
#                                  use_inv_prob_weighting = c(FALSE, TRUE),
#                                  use_SMOTE = c(FALSE, FALSE))


for (i in 1:nrow(grouping_param_df)) {
  grouping_type = grouping_param_df$grouping_type[i]
  grouping_var = grouping_param_df$grouping_var[i]
  SVM_feature_var = grouping_param_df$SVM_feature_var[i]
  
  #### 10-fold linear SVM with different weights
  # Iterate over weighting_param_df 
  for (j in 1:nrow(weighting_param_df)) {
    weighting_name <- weighting_param_df$name[j]
    use_inv_prob_weighting <- weighting_param_df$use_inv_prob_weighting[j]
    use_SMOTE <- weighting_param_df$use_SMOTE[j]
    
    # Run given weighting for 10-fold CV linear SVM
    if (!file.exists(paste0(rdata_path, sprintf("%s_wise_CV_linear_SVM_%s_%s.Rds",
                                                grouping_type, feature_set, weighting_name)))) {
      group_wise_SVM_CV_weighting <- run_univariate_cv_svm_by_input_var(rdata_path = rdata_path,
                                                                        feature_set = feature_set,
                                                                        test_package = test_package,
                                                                        svm_kernel = kernel,
                                                                        grouping_var = grouping_var,
                                                                        flds = subject_folds,
                                                                        svm_feature_var = SVM_feature_var,
                                                                        out_of_sample_only = TRUE,
                                                                        use_inv_prob_weighting = use_inv_prob_weighting,
                                                                        use_SMOTE = use_SMOTE,
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
    
    #### Calculate p values from model-free shuffle null distribution
    if (!file.exists(paste0(rdata_path, sprintf("%s_wise_CV_linear_SVM_%s_%s_model_free_shuffle_pvals.Rds",
                                                grouping_type,
                                                feature_set, 
                                                weighting_name)))) {
      
      group_wise_SVM_balanced_accuracy <- readRDS(paste0(rdata_path, sprintf("%s_wise_CV_linear_SVM_%s_%s_balacc.Rds",
                                                                             grouping_type, feature_set, weighting_name)))
      
      # Calculate p-values
      pvalues <- calc_empirical_nulls(class_res = group_wise_SVM_balanced_accuracy,
                                      null_data = model_free_shuffle_null_res,
                                      feature_set = feature_set,
                                      use_pooled_null = TRUE,
                                      is_main_data_averaged = TRUE,
                                      grouping_var = grouping_var)
      
      saveRDS(pvalues, file=paste0(rdata_path, sprintf("%s_wise_CV_linear_SVM_%s_%s_model_free_shuffle_pvals.Rds",
                                                       grouping_type,
                                                       feature_set, 
                                                       weighting_name)))
    }
    
    # Generate null-model fits distribution
    if (!file.exists(paste0(rdata_path, sprintf("%s_wise_model_permutation_null_%s_%s.Rds",
                                                grouping_type,
                                                feature_set,
                                                weighting_name)))) {
      
      # template file
      if (grouping_type %in% c("ROI", "Feature")) {
        num_permutations <- 50
        nperm_per_iter <- 20
        wall_hrs <- "6"
      } else if (grouping_type == "Combo") {
        num_permutations <- 75
        nperm_per_iter <- 150
        wall_hrs <- "12"
      }
      num_k_folds <- 10
      template_pbs_file <- paste0(github_dir, "univariate_analysis/template_null_model_fit.pbs")
      
      output_data_dir <- paste0(rdata_path, sprintf("%s_wise_%s_%s_null_model_fits/",
                                                    grouping_type, feature_set, weighting_name))
      
      output_scripts_dir <- paste0(github_dir, sprintf("univariate_analysis/%s_wise_%s_%s_null_model_fits/",
                                                       grouping_type, feature_set, weighting_name))
      icesTAF::mkdir(output_data_dir)
      icesTAF::mkdir(output_scripts_dir)
      
      lookup_list <- list("PROJECT_NAME" = "hctsa",
                          "NAME" = sprintf("univariate_%s_wise_null_model_fit",
                                           grouping_type),
                          "MEMNUM" = "20",
                          "NCPUS" = "1",
                          "GITHUB_DIR" = github_dir,
                          "PROJECT_DIR" = project_path,
                          "EMAIL" = "abry4213@uni.sydney.edu.au",
                          "PBS_NOTIFY" = "a",
                          "WALL_HRS" = wall_hrs,
                          "NUM_K_FOLDS" = num_k_folds,
                          "NUM_PERMS_PER_ITER" = nperm_per_iter,
                          "OUTPUT_DATA_DIR" = output_data_dir,
                          "FEATURE_SET" = "catch22",
                          "GROUPING_VAR" = grouping_var,
                          "SVM_FEATURE_VAR" = SVM_feature_var,
                          "WEIGHTING_NAME" = weighting_name)
      
      to_be_replaced <- names(lookup_list)
      replacement_values <- unlist(unname(lookup_list))
      
      for (p in 1:num_permutations) {
        
        # Run command if null file doesn't exist
        if (!file.exists(sprintf("%s/%s_wise_%s_%s_null_model_fit_iter_%s.Rds",
                                 output_data_dir, grouping_var, feature_set, weighting_name, p))) {
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
      
    }
    
  }
}
