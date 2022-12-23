#-------------------------------------------------------------------------------
# Define paths
#-------------------------------------------------------------------------------
# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--github_dir", default="/headnode1/abry4213/github/")
parser$add_argument("--data_path", default="/headnode1/abry4213/data/UCLA_CNP_ABIDE_ASD/")
parser$add_argument("--sample_metadata_file", default="UCLA_Schizophrenia_sample_metadata.Rds")
parser$add_argument("--pairwise_feature_set", default="pyspi14")
parser$add_argument("--univariate_feature_set", default="catch22")
parser$add_argument("--dataset_ID", default="UCLA_CNP_ABIDE_ASD")
parser$add_argument("--email")
parser$add_argument("--add_catch2", action="store_true", default=FALSE)

# Parse input arguments
args <- parser$parse_args()
github_dir <- args$github_dir
data_path <- args$data_path
rdata_path <- args$rdata_path
pairwise_feature_set <- args$pairwise_feature_set
univariate_feature_set <- args$univariate_feature_set
dataset_ID <- args$dataset_ID
sample_metadata_file <- args$sample_metadata_file
email <- args$email
add_catch2 <- args$add_catch2

# univariate_feature_set <- "catch22"
# pairwise_feature_set <- "pyspi14"
# github_dir <- "~/github/"
# email <- "abry4213@uni.sydney.edu.au"
# data_path <- "~/data/UCLA_CNP_ABIDE_ASD/"
# dataset_ID <- "UCLA_CNP_ABIDE_ASD"
# sample_metadata_file <- "UCLA_CNP_ABIDE_ASD_sample_metadata.Rds"
# add_catch2 <- TRUE

rdata_path <- paste0(data_path, "processed_data/Rdata/")
plot_dir <- paste0(data_path, "plots/")

cat("\nAdd catch2?", add_catch2, "\n")

if (add_catch2) {
  feature_sets <- c(univariate_feature_set, "catch2")
} else {
  feature_sets <- c(univariate_feature_set)
}

cat("\nPrepping null model scripts for the following feature sets:\n")
print(feature_sets)

TAF::mkdir(plot_dir)

# Set the seed
set.seed(127)

# Load tidyverse
library(tidyverse)

grouping_param_df <- data.frame(grouping_type = c("ROI", "Feature", "Combo"),
                                grouping_var = c("Brain_Region", "Feature", "Combo"),
                                SVM_feature_var = c("Feature", "Brain_Region", "Combo")) 

# Use a linear kernel
svm_kernel = "linear"
weighting_name <- "inv_prob"
use_inv_prob_weighting <- TRUE

for (i in 1:nrow(grouping_param_df)) {
  grouping_type = grouping_param_df$grouping_type[i]
  grouping_var = grouping_param_df$grouping_var[i]
  SVM_feature_var = grouping_param_df$SVM_feature_var[i]
  
  ############################################################################
  # Null model fits
  ############################################################################
  
  # We want to run 1,000 null model fits, and we can run 10 permutations per iteration
  if (grouping_type == "Combo") {
    num_permutations <- 100
    nperm_per_iter <- 50
    wall_hrs <- "24"
  } else {
    num_permutations <- 100
    nperm_per_iter <- 10
    wall_hrs <- "12"
  }
  
  # Use 10-fold cross-validation
  num_k_folds <- 10
  # Define the univariate template PBS script
  template_pbs_file <- paste0(github_dir, "fMRI_FeaturesDisorders/helper_functions/classification/template_univariate_null_model_fit.pbs")
  
  for (feature_set in feature_sets) {
        # Where to store null model fit results
    output_data_dir <- paste0(rdata_path, sprintf("%s_%s_wise_%s_%s_null_model_fits/",
                                                  dataset_ID,
                                                  grouping_type,
                                                  feature_set,
                                                  weighting_name))
    
    # Where to save PBS script to
    output_scripts_dir <- paste0(github_dir, sprintf("fMRI_FeaturesDisorders/classification_analysis/univariate_analysis/null_pbs_scripts/%s_%s_wise_%s_%s_null_model_fits/",
                                                    dataset_ID,
                                                    grouping_type,
                                                    feature_set,
                                                    weighting_name))
    
    cat("\nNow generating null PBS scripts for", grouping_type, "\n")
    cat("Script location:", output_scripts_dir, "\n")
    
    # Make these directories
    TAF::mkdir(output_data_dir)
    TAF::mkdir(output_scripts_dir)
    
    # Lookup table for PBS script
    lookup_list <- list("NAME" = sprintf("univariate_%s_wise_%s",
                                        grouping_type, weighting_name),
                        "MEMNUM" = "20",
                        "NCPUS" = "1",
                        "DATASET_ID" = dataset_ID,
                        "GITHUB_DIR" = github_dir,
                        "DATA_PATH" = data_path,
                        "RDATA_PATH" = rdata_path,
                        "EMAIL" = email,
                        "PBS_NOTIFY" = "a",
                        "WALL_HRS" = wall_hrs,
                        "NUM_K_FOLDS" = num_k_folds,
                        "NUM_PERMS_PER_ITER" = nperm_per_iter,
                        "OUTPUT_DATA_DIR" = output_data_dir,
                        "SAMPLE_METADATA_FILE" = sample_metadata_file,
                        "UNIVARIATE_FEATURE_SET" = feature_set,
                        "PAIRWISE_FEATURE_SET" = pairwise_feature_set,
                        "GROUPING_VAR" = grouping_var,
                        "SVM_FEATURE_VAR" = SVM_feature_var,
                        "WEIGHTING_NAME" = weighting_name)
    
    to_be_replaced <- names(lookup_list)
    replacement_values <- unlist(unname(lookup_list))
    
    # Create a PBS script per iteration
    for (p in 1:num_permutations) {
      
      # Run command if null file doesn't exist
      if (!file.exists(sprintf("%s/%s_%s_wise_%s_%s_null_model_fit_iter_%s.Rds",
                              output_data_dir, dataset_ID,
                              grouping_var, feature_set,
                              weighting_name, p))) {
        new_pbs_file <- readLines(template_pbs_file)
        
        # Replace file paths
        pbs_text_replaced <- mgsub::mgsub(new_pbs_file,
                                          to_be_replaced,
                                          replacement_values)
        
        # Replace null iteration number
        pbs_text_replaced <- gsub("iterj", p, pbs_text_replaced)
        
        # Write updated PBS script to fil
        output_pbs_file <- writeLines(pbs_text_replaced,
                                      paste0(output_scripts_dir,
                                            "null_iter_", p, ".pbs"))
        
      }
    }
  }
}