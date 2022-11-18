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

# Parse input arguments
args <- parser$parse_args()
github_dir <- args$github_dir
data_path <- args$data_path
sample_metadata_file <- args$sample_metadata_file
pairwise_feature_set <- args$pairwise_feature_set
univariate_feature_set <- args$univariate_feature_set
noise_proc_for_null <- args$noise_proc_for_null
dataset_ID <- args$dataset_ID
email <- args$email

cat("noise_proc_for_null:", noise_proc_for_null, "\n")

# univariate_feature_set <- "catch22"
# pairwise_feature_set <- "pyspi14"
# github_dir <- "/headnode/abry4213/github/"
# email <- "abry4213@uni.sydney.edu.au"

# UCLA schizophrenia
# data_path <- "/headnode1/abry4213/data/UCLA_Schizophrenia/"
# dataset_ID <- "UCLA_Schizophrenia"
# sample_metadata_file <- "UCLA_Schizophrenia_sample_metadata.Rds"
# noise_proc_for_null <- "AROMA+2P+GMR"

# ABIDE ASD
# data_path <- "/headnode1/abry4213/data/ABIDE_ASD/"
# sample_metadata_file <- "ABIDE_ASD_sample_metadata.Rds"
# dataset_ID <- "ABIDE_ASD"
# noise_proc_for_null <- "FC1000"

rdata_path <- paste0(data_path, "processed_data/Rdata/")
plot_dir <- paste0(data_path, "plots/")

icesTAF::mkdir(plot_dir)

# Set the seed
set.seed(127)

# Load tidyverse
library(tidyverse)

# Unlist noise-processing methods
tryCatch({
  noise_procs <- stringr::str_split(noise_procs, ";")[[1]]
  noise_procs <- unlist(noise_procs)
}, error = function(e) {})

grouping_param_df <- data.frame(grouping_var = c("SPI"),
                                SVM_feature_var = c("region_pair"))

# Use a linear kernel
kernel = "linear"
weighting_name <- "inv_prob"
use_inv_prob_weighting <- TRUE

for (i in 1:nrow(grouping_param_df)) {
  grouping_var = grouping_param_df$grouping_var[i]
  SVM_feature_var = grouping_param_df$SVM_feature_var[i]
  
  ############################################################################
  # Null model fits
  ############################################################################
  
  # We want to run 1,000 null model fits, and we can run 10 permutations per iteration
  num_permutations <- 200
  nperm_per_iter <- 5
  wall_hrs <- "24"
  # Use 10-fold cross-validation
  num_k_folds <- 10
  # Define the univariate template PBS script
  template_pbs_file <- paste0(github_dir, "fMRI_FeaturesDisorders/helper_functions/classification/template_pairwise_null_model_fit.pbs")
  
  # Where to store null model fit results
  output_data_dir <- paste0(rdata_path, sprintf("%s_%s_wise_%s_%s_null_model_fits/",
                                                dataset_ID,
                                                grouping_var,
                                                pairwise_feature_set,
                                                weighting_name))
  
  # Where to save PBS script to
  output_scripts_dir <- paste0(github_dir, sprintf("fMRI_FeaturesDisorders/classification_analysis/pairwise_analysis/null_pbs_scripts/%s_%s_wise_%s_%s_null_model_fits/",
                                                   dataset_ID,
                                                   grouping_var,
                                                   pairwise_feature_set,
                                                   weighting_name))
  
  cat("\nNow generating null PBS scripts for", grouping_var, "\n")
  cat("Script location:", output_scripts_dir, "\n")
  
  # Make these directories
  icesTAF::mkdir(output_data_dir)
  icesTAF::mkdir(output_scripts_dir)
  
  
  # Lookup table for PBS script
  lookup_list <- list("NAME" = sprintf("pairwise_%s_wise_null_model_fit",
                                       grouping_var),
                      "MEMNUM" = "20",
                      "NCPUS" = "1",
                      "DATASET_ID" = dataset_ID,
                      "GITHUB_DIR" = github_dir,
                      "DATA_PATH" = data_path,
                      "RDATA_PATH" = rdata_path,
                      "PAIRWISE_DATA_FILE" = pyspi_data_file,
                      "SPI_DIRECTIONALITY_FILE" = SPI_directionality_file,
                      "EMAIL" = email,
                      "PBS_NOTIFY" = "a",
                      "WALL_HRS" = wall_hrs,
                      "NOISE_PROC_FOR_NULL" = noise_proc_for_null,
                      "NUM_K_FOLDS" = num_k_folds,
                      "NUM_PERMS_PER_ITER" = nperm_per_iter,
                      "OUTPUT_DATA_DIR" = output_data_dir,
                      "SAMPLE_METADATA_FILE" = sample_metadata_file,
                      "UNIVARIATE_FEATURE_SET" = univariate_feature_set,
                      "PAIRWISE_FEATURE_SET" = pairwise_feature_set,
                      "GROUPING_VAR" = grouping_var,
                      "SVM_FEATURE_VAR" = SVM_feature_var,
                      "WEIGHTING_NAME" = weighting_name)
  
  to_be_replaced <- names(lookup_list)
  replacement_values <- unlist(unname(lookup_list))
  
  # Create a PBS script per iteration
  for (p in 1:num_permutations) {
    
    # Run command if null file doesn't exist
    if (!file.exists(sprintf("%s/%s_wise_%s_%s_null_model_fit_iter_%s.Rds",
                             output_data_dir, grouping_var, pairwise_feature_set,
                             weighting_name, p))) {
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
      
    }
  }
}