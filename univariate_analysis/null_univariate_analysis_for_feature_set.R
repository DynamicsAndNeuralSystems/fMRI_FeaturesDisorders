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
# feature_set <- "catchaMouse16"

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
# Define weighting parameters
################################################################################

grouping_param_df <- data.frame(grouping_type = c("ROI", "Feature", "Combo"),
                                grouping_var = c("Brain_Region", "Feature", "Combo"),
                                SVM_feature_var = c("Feature", "Brain_Region", "Combo")) 

weighting_param_df <- data.frame(name = c("unweighted", "inv_prob"),
                                 use_inv_prob_weighting = c(FALSE, TRUE),
                                 use_SMOTE = c(FALSE, FALSE))


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
    
    
    # Generate null-model fits distribution
    if (!file.exists(paste0(rdata_path, sprintf("%s_wise_model_permutation_null_%s_%s.Rds",
                                                grouping_type,
                                                feature_set,
                                                weighting_name)))) {
      # Define data directory
      output_data_dir <- paste0(rdata_path, sprintf("%s_wise_%s_%s_null_model_fits/",
                                                    grouping_type, feature_set, weighting_name))
      
      ## Concatenate null results and save to RDS file
      model_permutation_null_weighting <- list.files(output_data_dir, pattern="Rds") %>%
        purrr::map_df(~ readRDS(paste0(output_data_dir, .x)))
      saveRDS(model_permutation_null_weighting, paste0(rdata_path, sprintf("%s_wise_model_permutation_null_%s_%s.Rds",
                                                                           grouping_type,
                                                                           feature_set,
                                                                           weighting_name)))
    } else {
      model_permutation_null_weighting <- readRDS(paste0(rdata_path, sprintf("%s_wise_model_permutation_null_%s_%s.Rds",
                                                                             grouping_type,
                                                                             feature_set,
                                                                             weighting_name)))
    }
    
    # Empirically derive p-values based on null model fits distribution
    if (!file.exists(paste0(rdata_path, sprintf("%s_wise_CV_linear_SVM_model_permutation_null_%s_%s_pvals.Rds",
                                                grouping_type,
                                                feature_set,
                                                weighting_name)))) {
      group_wise_SVM_balanced_accuracy <- readRDS(paste0(rdata_path, 
                                                         sprintf("%s_wise_CV_linear_SVM_%s_%s_balacc.Rds",
                                                                             grouping_type, feature_set, weighting_name)))
      
      # Calculate p-values
      pvalues <- calc_empirical_nulls(class_res = group_wise_SVM_balanced_accuracy,
                                      null_data = model_permutation_null_weighting,
                                      feature_set = feature_set,
                                      use_pooled_null = TRUE,
                                      is_main_data_averaged = TRUE,
                                      grouping_var = grouping_var)
      
      saveRDS(pvalues, file=paste0(rdata_path, sprintf("%s_wise_CV_linear_SVM_model_permutation_null_%s_%s_pvals.Rds",
                                                       grouping_type,
                                                       feature_set,
                                                       weighting_name)))
    }
  }
}
