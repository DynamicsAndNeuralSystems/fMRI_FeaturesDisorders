#-------------------------------------------------------------------------------
# Define paths
#-------------------------------------------------------------------------------
# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--github_dir", default="/headnode1/abry4213/github/")
parser$add_argument("--data_path", default="/headnode1/abry4213/data/UCLA_Schizophrenia/")
parser$add_argument("--sample_metadata_file", default="UCLA_Schizophrenia_sample_metadata.Rds")
parser$add_argument("--pairwise_feature_set", default="pyspi14")
parser$add_argument("--univariate_feature_set", default="catch22")
parser$add_argument("--noise_proc_for_null", default=c(""))
parser$add_argument("--dataset_ID", default="UCLA_Schizophrenia")
parser$add_argument("--email")
parser$add_argument("--run_number", nargs='?')

# Parse input arguments
args <- parser$parse_args()
github_dir <- args$github_dir
data_path <- args$data_path
rdata_path <- args$rdata_path
pairwise_feature_set <- args$pairwise_feature_set
univariate_feature_set <- args$univariate_feature_set
noise_proc_for_null <- args$noise_proc_for_null
dataset_ID <- args$dataset_ID
sample_metadata_file <- args$sample_metadata_file
email <- args$email
run_number <- args$run_number

# 
# univariate_feature_set <- "catch22"
# pairwise_feature_set <- "pyspi14"
# github_dir <- "/headnode1/abry4213/github/"
# email <- "abry4213@uni.sydney.edu.au"

# UCLA schizophrenia
# data_path <- "/headnode1/abry4213/data/UCLA_Schizophrenia/"
# dataset_ID <- "UCLA_Schizophrenia"
# sample_metadata_file <- "UCLA_Schizophrenia_sample_metadata.Rds"
# noise_procs <- "AROMA+2P;AROMA+2P+GMR;AROMA+2P+DiCER"
# noise_proc_for_null <- "AROMA+2P+GMR"

# ABIDE ASD
# data_path <- "/headnode1/abry4213/data/ABIDE_ASD/"
# sample_metadata_file <- "ABIDE_ASD_sample_metadata.Rds"
# dataset_ID <- "ABIDE_ASD"
# noise_procs <- c("FC1000")
# noise_proc_for_null <- "FC1000"

if (!is.null(run_number)) {
  rdata_path <- paste0(data_path, "processed_data_run", run_number, "/Rdata/")
  plot_dir <- paste0(data_path, "plots_run", run_number, "/")
} else {
  rdata_path <- paste0(data_path, "processed_data/Rdata/")
  plot_dir <- paste0(data_path, "plots/")
}

icesTAF::mkdir(plot_dir)

# Set the seed
set.seed(127)

# Load tidyverse
library(tidyverse)


#-------------------------------------------------------------------------------
# Source helper scripts
#-------------------------------------------------------------------------------
helper_script_dir = paste0(github_dir, "fMRI_FeaturesDisorders/helper_functions/classification/")
source(paste0(helper_script_dir, "Linear_SVM.R"))
source(paste0(helper_script_dir, "Null_distributions.R"))

# Load sample metadata
sample_metadata <- readRDS(paste0(data_path, sample_metadata_file))

kernel = "linear"
noise_label = gsub("\\+", "_", noise_proc_for_null)


################################################################################
# Define grouping parameters
################################################################################
weighting_name="inv_prob"

grouping_param_df <- data.frame(grouping_var = "SPI",
                                SVM_feature_var = "region_pair")

#### 10-fold linear SVM with different weights
# Iterate over grouping_param_df
for (i in 1:nrow(grouping_param_df)) {
  grouping_var = grouping_param_df$grouping_var[i]
  SVM_feature_var = grouping_param_df$SVM_feature_var[i]
  
  # Generate null-model fits distribution
  if (!file.exists(paste0(rdata_path, sprintf("%s_univariate_%s_pairwise_%s_model_permutation_null_%s.Rds",
                                              dataset_ID,
                                              univariate_feature_set,
                                              pairwise_feature_set,
                                              weighting_name)))) {
    cat("Now merging null model fits data\n")
    # Define data directory
    run_number <- ifelse(is.null(run_number), "", run_number)
    output_data_dir <- paste0(rdata_path, sprintf("%s_univariate_%s_pairwise_%s_%s_null_model_fits%s/",
                                                  dataset_ID,
                                                  univariate_feature_set,
                                                  pairwise_feature_set,
                                                  weighting_name,
                                                  run_number))
    
    
    ## Concatenate null results and save to RDS file
    model_permutation_null_weighting <- list.files(output_data_dir, pattern="Rds") %>%
      purrr::map_df(~ readRDS(paste0(output_data_dir, .x)))
    cat("Now saving null model fits data\n")
    saveRDS(model_permutation_null_weighting, paste0(rdata_path, sprintf("%s_univariate_%s_pairwise_%s_model_permutation_null_%s.Rds",
                                                                         dataset_ID,
                                                                         univariate_feature_set,
                                                                         pairwise_feature_set,
                                                                         weighting_name)))
  } else {
    model_permutation_null_weighting <- readRDS(paste0(rdata_path, sprintf("%s_univariate_%s_pairwise_%s_model_permutation_null_%s.Rds",
                                                                           dataset_ID,
                                                                           univariate_feature_set,
                                                                           pairwise_feature_set,
                                                                           weighting_name)))
  }
  

  
  # Empirically derive p-values based on null model fits distribution
  if (!file.exists(paste0(rdata_path, sprintf("univariate_%s_pairwise_%s_CV_linear_SVM_model_permutation_null_%s_pvals.Rds",
                                              univariate_feature_set,
                                              pairwise_feature_set,
                                              weighting_name)))) {
    cat("Now calculating p-values\n")
    univariate_pairwise_SVM_balanced_accuracy_across_repeats <- readRDS(paste0(rdata_path, sprintf("univariate_%s_pairwise_%s_CV_linear_SVM_%s_balacc_across_repeats.Rds",
                                                                                    univariate_feature_set, 
                                                                                    pairwise_feature_set,
                                                                                    weighting_name)))
    
    # Calculate p-values
    pvalues <- calc_empirical_nulls(class_res = univariate_pairwise_SVM_balanced_accuracy_across_repeats,
                                    null_data = model_permutation_null_weighting,
                                    feature_set = pairwise_feature_set,
                                    noise_proc = noise_proc_for_null,
                                    use_pooled_null = TRUE,
                                    is_main_data_averaged = TRUE,
                                    grouping_var = grouping_var)
    
    saveRDS(pvalues, file=paste0(rdata_path, sprintf("univariate_%s_pairwise_%s_CV_linear_SVM_model_permutation_null_%s_pvals.Rds",
                                                     univariate_feature_set,
                                                     pairwise_feature_set,
                                                     weighting_name)))
  }
}