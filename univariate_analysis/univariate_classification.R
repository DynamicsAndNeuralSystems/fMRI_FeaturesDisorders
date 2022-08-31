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
parser$add_argument("--noise_procs", default=c(""), nargs="*", action="append")
parser$add_argument("--noise_proc_for_null", default=c(""))
parser$add_argument("--dataset_ID", default="UCLA_Schizophrenia")
parser$add_argument("--run_number")

# Parse input arguments
args <- parser$parse_args()
github_dir <- args$github_dir
data_path <- args$data_path
rdata_path <- args$rdata_path
pairwise_feature_set <- args$pairwise_feature_set
univariate_feature_set <- args$univariate_feature_set
noise_procs <- args$noise_procs
noise_proc_for_null <- args$noise_proc_for_null
dataset_ID <- args$dataset_ID
sample_metadata_file <- args$sample_metadata_file
run_number <- args$run_number
# 
# univariate_feature_set <- "catch22"
# pairwise_feature_set <- "pyspi14"
# github_dir <- "/headnode1/abry4213/github/"

# UCLA schizophrenia
# data_path <- "/headnode1/abry4213/data/UCLA_Schizophrenia/"
# dataset_ID <- "UCLA_Schizophrenia"
# sample_metadata_file <- "UCLA_Schizophrenia_sample_metadata.Rds"
# noise_procs <- "AROMA+2P;AROMA+2P+GMR;AROMA+2P+DiCER"
# noise_proc_for_null <- "AROMA+2P+GMR"
# run_number <- 1

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

# Unlist noise-processing methods
tryCatch({
  noise_procs <- stringr::str_split(noise_procs, ";")[[1]]
  noise_procs <- unlist(noise_procs)
}, error = function(e) {})


#-------------------------------------------------------------------------------
# Source helper scripts
#-------------------------------------------------------------------------------
# Set working directory to file location
helper_script_dir = paste0(github_dir, "fMRI_FeaturesDisorders/helper_functions/classification/")
source(paste0(helper_script_dir, "Linear_SVM.R"))
source(paste0(helper_script_dir, "Null_distributions.R"))

# Load sample metadata
sample_metadata <- readRDS(paste0(data_path, sample_metadata_file))

################################################################################
# Create ten folds to use for all analyses
################################################################################
subjects_to_use <- readRDS(paste0(rdata_path, sprintf("%s_samples_with_univariate_%s_and_pairwise_%s_filtered.Rds",
                                                     dataset_ID,
                                                     univariate_feature_set,
                                                     pairwise_feature_set)))


if (!file.exists(paste0(rdata_path, dataset_ID, "_samples_per_10_folds.Rds"))) {
  # Make folds
  set.seed(127)
  k = 10
  samples_with_diagnosis <- subjects_to_use %>%
    left_join(., sample_metadata)
  sample_folds <- caret::createFolds(samples_with_diagnosis$Diagnosis, k = k, list = TRUE, returnTrain = FALSE)
  
  # Save to Rds file
  saveRDS(sample_folds, file=paste0(rdata_path, dataset_ID, "_samples_per_10_folds.Rds"))
} else {
  sample_folds <- readRDS(paste0(rdata_path, dataset_ID, "_samples_per_10_folds.Rds"))
}

################################################################################
# Define weighting parameters
################################################################################

grouping_param_df <- data.frame(grouping_type = c("ROI", "Feature", "Combo"),
                                grouping_var = c("Brain_Region", "Feature", "Combo"),
                                SVM_feature_var = c("Feature", "Brain_Region", "Combo")) 

weighting_param_df <- data.frame(name = c("inv_prob"),
                                 use_inv_prob_weighting = c(TRUE))

# Use a linear kernel
kernel = "linear"

for (i in 1:nrow(grouping_param_df)) {
  grouping_type = grouping_param_df$grouping_type[i]
  grouping_var = grouping_param_df$grouping_var[i]
  SVM_feature_var = grouping_param_df$SVM_feature_var[i]
  
  #### 10-fold linear SVM with different weights
  # Iterate over weighting_param_df 
  for (j in 1:nrow(weighting_param_df)) {
    weighting_name <- weighting_param_df$name[j]
    use_inv_prob_weighting <- weighting_param_df$use_inv_prob_weighting[j]
    
    # Run given weighting for 10-fold CV linear SVM
    if (!file.exists(paste0(rdata_path, sprintf("%s_wise_CV_linear_SVM_%s_%s.Rds",
                                                grouping_type, 
                                                univariate_feature_set, 
                                                weighting_name)))) {
      group_wise_SVM_CV_weighting <- run_univariate_cv_svm_by_input_var(data_path = data_path,
                                                                        dataset_ID = dataset_ID,
                                                                        sample_metadata = sample_metadata,
                                                                        univariate_feature_set = univariate_feature_set,
                                                                        pairwise_feature_set = pairwise_feature_set,
                                                                        svm_kernel = kernel,
                                                                        grouping_var = grouping_var,
                                                                        flds = sample_folds,
                                                                        svm_feature_var = SVM_feature_var,
                                                                        out_of_sample_only = TRUE,
                                                                        use_inv_prob_weighting = use_inv_prob_weighting,
                                                                        noise_procs = noise_procs)
      saveRDS(group_wise_SVM_CV_weighting, file=paste0(rdata_path, 
                                                       sprintf("%s_wise_CV_linear_SVM_%s_%s.Rds",
                                                               grouping_type,
                                                               univariate_feature_set, 
                                                               weighting_name)))
    } else {
      group_wise_SVM_CV_weighting <- readRDS(paste0(rdata_path, 
                                                    sprintf("%s_wise_CV_linear_SVM_%s_%s.Rds",
                                                            grouping_type,
                                                            univariate_feature_set, 
                                                            weighting_name)))
    }
    
    #### Calculate balanced accuracy across all folds
    if (!file.exists(paste0(rdata_path, sprintf("%s_wise_CV_linear_SVM_%s_%s_balacc.Rds",
                                                grouping_type, 
                                                univariate_feature_set, 
                                                weighting_name)))) {
      group_wise_SVM_balanced_accuracy <- group_wise_SVM_CV_weighting %>%
        group_by(grouping_var, Noise_Proc, Sample_Type) %>%
        summarise(accuracy = sum(Prediction_Correct) / n(),
                  balanced_accuracy = caret::confusionMatrix(data = Predicted_Diagnosis,
                                                             reference = Actual_Diagnosis)$byClass[["Balanced Accuracy"]])
      
      saveRDS(group_wise_SVM_balanced_accuracy, file=paste0(rdata_path, sprintf("%s_wise_CV_linear_SVM_%s_%s_balacc.Rds",
                                                                                grouping_type, 
                                                                                univariate_feature_set, 
                                                                                weighting_name)))
    } else {
      group_wise_SVM_balanced_accuracy <- readRDS(paste0(rdata_path, sprintf("%s_wise_CV_linear_SVM_%s_%s_balacc.Rds",
                                                                             grouping_type, 
                                                                             univariate_feature_set, 
                                                                             weighting_name)))
    }
    
    ############################################################################
    # Null model fits
    ############################################################################
    
    # We want to run 1,000 null model fits, and we can run 10 permutations per iteration
    num_permutations <- 100
    nperm_per_iter <- 10
    wall_hrs <- "12"
    # Use 10-fold cross-validation
    num_k_folds <- 10
    # Define the univariate template PBS script
    template_pbs_file <- paste0(github_dir, "fMRI_FeaturesDisorders/helper_functions/classification/template_univariate_null_model_fit.pbs")

    # Where to store null model fit results
    output_data_dir <- paste0(rdata_path, sprintf("%s_%s_wise_%s_%s_null_model_fits/",
                                                  dataset_ID,
                                                  grouping_type,
                                                  univariate_feature_set,
                                                  weighting_name))

    # Where to save PBS script to
    output_scripts_dir <- paste0(github_dir, sprintf("fMRI_FeaturesDisorders/univariate_analysis/%s_%s_wise_%s_%s_null_model_fits/",
                                                     dataset_ID,
                                                     grouping_type,
                                                     univariate_feature_set,
                                                     weighting_name))
    
    # Make these directories
    icesTAF::mkdir(output_data_dir)
    icesTAF::mkdir(output_scripts_dir)

    # Lookup table for PBS script
    lookup_list <- list("NAME" = sprintf("univariate_%s_wise_null_model_fit",
                                         grouping_type),
                        "MEMNUM" = "20",
                        "NCPUS" = "1",
                        "DATASET_ID" = dataset_ID,
                        "GITHUB_DIR" = github_dir,
                        "DATA_PATH" = data_path,
                        "EMAIL" = "abry4213@uni.sydney.edu.au",
                        "PBS_NOTIFY" = "a",
                        "WALL_HRS" = wall_hrs,
                        "NOISE_PROCS" = noise_proc_for_null,
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
                               output_data_dir, grouping_var, univariate_feature_set,
                               weighting_name, p))) {
        cat("\nNow creating pbs script for for iteration", p, "\n")
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
}
