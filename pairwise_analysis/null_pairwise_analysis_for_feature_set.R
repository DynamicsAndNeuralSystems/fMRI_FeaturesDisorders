# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--project_path", default="/project/hctsa/annie/")
parser$add_argument("--github_dir", default="/project/hctsa/annie/github/fMRI_FeaturesDisorders/")
parser$add_argument("--rdata_path", default="/project/hctsa/annie/data/scz/UCLA/Rdata/")
parser$add_argument("--pydata_path", default="/project/hctsa/annie/data/scz/UCLA/pydata/")
parser$add_argument("--feature_set", default="pyspi_19")
# project_path <- "D:/Virtual_Machines/Shared_Folder/"
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

# ###############################################################################
# Load pyspi data
# ###############################################################################

SPI_directionality <- read.csv(paste0(github_dir, "pairwise_analysis/SPI_Direction_Info.csv"))


################################################################################
# Define weighting parameters
################################################################################

weighting_param_df <- data.frame(name = c("inv_prob"),
                                 use_inv_prob_weighting = c(TRUE),
                                 use_SMOTE = c(FALSE))

SVM_grouping_params <- data.frame(grouping_var = c("SPI"),
                                  SVM_feature_var = c("region_pair"))


################################################################################
# Run linear SVM for each grouping var
################################################################################

for (i in 1:nrow(weighting_param_df)) {
  weighting_name <- weighting_param_df$name[i]
  use_inv_prob_weighting <- weighting_param_df$use_inv_prob_weighting[i]
  use_SMOTE <- weighting_param_df$use_SMOTE[i]
  
  for (j in 1:nrow(SVM_grouping_params)) {
    grouping_var = SVM_grouping_params$grouping_var[j]
    SVM_feature_var = SVM_grouping_params$SVM_feature_var[j]
    
    weighting_null_dist_file <- paste0(rdata_path, sprintf("pyspi_%s_pairwise_%s_%s_null_model_fits.Rds",
                                                           grouping_var, feature_set, weighting_name))
    
    # Run null perm iterations if overall null distribution data file doesn't exist
    if (!file.exists(weighting_null_dist_file)) {

      # Output script dir
      output_data_dir <- paste0(rdata_path, sprintf("pyspi_%s_pairwise_%s_%s_null_model_fits/",
                                                    grouping_var, feature_set, weighting_name))
      output_scripts_dir <- paste0(github_dir, sprintf("pairwise_analysis/pyspi_%s_pairwise_%s_%s_null_model_fits/",
                                                       grouping_var, weighting_name, feature_set))
      
      ## Concatenate null results and save to RDS file
      null_model_fit_res <- list.files(output_data_dir, pattern="Rds") %>%
        purrr::map_df(~ readRDS(paste0(output_data_dir, .x)))
      saveRDS(null_model_fit_res, paste0(rdata_path, sprintf("pyspi_%s_pairwise_%s_%s_null_model_fits.Rds",
                                                             grouping_var, feature_set, weighting_name)))
    } else {
      null_model_fit_res <- readRDS(paste0(rdata_path, sprintf("pyspi_%s_pairwise_%s_%s_null_model_fits.Rds",
                                                               grouping_var, feature_set, weighting_name)))
    }
    
    #### Calculate p-values from empirical model null distributions
    if (!file.exists(paste0(rdata_path, sprintf("pyspi_%s_pairwise_CV_linear_SVM_null_model_fits_pvals_%s_%s.Rds",
                                                grouping_var, feature_set, weighting_name)))) {
      group_wise_SVM_balanced_accuracy <- readRDS(paste0(rdata_path,
                                                         sprintf("pyspi_%s_pairwise_CV_linear_SVM_%s_%s_balacc.Rds",
                                                                 grouping_var,
                                                                 feature_set, 
                                                                 weighting_name)))

      # Calculate p-values
      pvalues <- calc_empirical_nulls(class_res = group_wise_SVM_balanced_accuracy,
                                      null_data = null_model_fit_res,
                                      feature_set = feature_set,
                                      use_pooled_null = TRUE,
                                      is_main_data_averaged = TRUE,
                                      is_null_data_averaged = TRUE,
                                      grouping_var = grouping_var)

      saveRDS(pvalues, file=paste0(rdata_path, sprintf("pyspi_%s_pairwise_CV_linear_SVM_null_model_fits_pvals_%s_%s.Rds",
                                                       grouping_var, feature_set, weighting_name)))
    }
  }
}