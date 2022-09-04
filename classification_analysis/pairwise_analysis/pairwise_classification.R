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
sample_metadata_file <- args$sample_metadata_file
pairwise_feature_set <- args$pairwise_feature_set
univariate_feature_set <- args$univariate_feature_set
noise_proc_for_null <- args$noise_proc_for_null
dataset_ID <- args$dataset_ID
email <- args$email
run_number <- args$run_number

cat("noise_proc_for_null:", noise_proc_for_null, "\n")
cat("run_number:", run_number, "\n")
# 
# univariate_feature_set <- "catch22"
# pairwise_feature_set <- "pyspi14"
# github_dir <- "/headnode1/abry4213/github/"

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

################################################################################
# Load pyspi data
################################################################################
pyspi_data_file <- sprintf("%s/%s_%s_filtered_zscored.Rds",
                           rdata_path, dataset_ID, pairwise_feature_set) 

pyspi_data <- readRDS(pyspi_data_file) %>%
  dplyr::filter(Noise_Proc %in% noise_proc_for_null) 

SPI_directionality_file <- paste0(github_dir, "fMRI_FeaturesDisorders/classification_analysis/pairwise_analysis/SPI_Direction_Info.csv")
SPI_directionality <- read.csv(SPI_directionality_file)

# Load sample metadata
sample_metadata <- readRDS(paste0(data_path, sample_metadata_file))

################################################################################
# Create ten folds to use for all analyses
################################################################################
subjects_to_use <- readRDS(paste0(rdata_path, sprintf("%s_samples_with_univariate_%s_and_pairwise_%s_filtered.Rds",
                                                      dataset_ID,
                                                      univariate_feature_set,
                                                      pairwise_feature_set)))

pyspi_data <- pyspi_data %>%
  semi_join(., subjects_to_use)

if (!file.exists(paste0(rdata_path, dataset_ID, "_samples_per_10_folds_10_repeats.Rds"))) {
  sample_folds <- list()
  for (i in 1:10) {
    # Make folds
    k = 10
    samples_with_diagnosis <- subjects_to_use %>%
      left_join(., sample_metadata)
    sample_folds_i <- caret::createFolds(samples_with_diagnosis$Diagnosis, k = k, list = TRUE, returnTrain = FALSE)
    sample_folds[[i]] <- sample_folds_i
  }
  
  # Save RDS file
  saveRDS(sample_folds, paste0(rdata_path, dataset_ID, "_samples_per_10_folds_10_repeats.Rds"))
} else {
  sample_folds <- readRDS(paste0(rdata_path, dataset_ID, "_samples_per_10_folds_10_repeats.Rds"))
}


################################################################################
# Define weighting parameters
################################################################################
# Use a linear kernel
kernel = "linear"
weighting_name <- "inv_prob"
use_inv_prob_weighting <- TRUE

grouping_param_df <- data.frame(grouping_var = c("SPI"),
                                SVM_feature_var = c("region_pair"))


for (i in 1:nrow(grouping_param_df)) {
  grouping_var = grouping_param_df$grouping_var[i]
  SVM_feature_var = grouping_param_df$SVM_feature_var[i]
  
  #### 10-fold linear SVM with different weights
  
  # Run given weighting for 10-fold CV linear SVM
  if (!file.exists(paste0(rdata_path, sprintf("%s_wise_CV_linear_SVM_%s_%s.Rds",
                                              grouping_var, 
                                              pairwise_feature_set, 
                                              weighting_name)))) {
    tryCatch({ 
      group_wise_SVM_CV_weighting_list <- list()
      for (idx in 1:length(sample_folds)) {
        tryCatch({
          if (!file.exists(paste0(rdata_path, sprintf("%s_wise_CV_linear_SVM_%s_%s_repeat%s.Rds",
                                                      grouping_var, pairwise_feature_set,
                                                      weighting_name, idx)))) {
            repeat_res <- run_pairwise_cv_svm_by_input_var(pairwise_data = pyspi_data,
                                                           dataset_ID = dataset_ID,
                                                           data_path = data_path,
                                                           rdata_path = rdata_path,
                                                           sample_metadata = sample_metadata,
                                                           SPI_directionality = SPI_directionality,
                                                           svm_kernel = kernel,
                                                           num_k_folds = 10,
                                                           flds = sample_folds[[idx]],
                                                           repeat_number = idx,
                                                           grouping_var = grouping_var,
                                                           svm_feature_var = SVM_feature_var,
                                                           noise_proc = noise_proc_for_null,
                                                           out_of_sample_only = TRUE,
                                                           use_inv_prob_weighting = use_inv_prob_weighting,
                                                           shuffle_labels = FALSE)
            saveRDS(repeat_res, file=paste0(rdata_path, sprintf("%s_wise_CV_linear_SVM_%s_%s_repeat%s.Rds",
                                                                grouping_var, pairwise_feature_set,
                                                                weighting_name, idx)))
          } else {
            repeat_res <- readRDS(paste0(rdata_path, sprintf("%s_wise_CV_linear_SVM_%s_%s_repeat%s.Rds",
                                                             grouping_var, pairwise_feature_set,
                                                             weighting_name, idx)))
          }
          group_wise_SVM_CV_weighting_list <- list.append(group_wise_SVM_CV_weighting_list, repeat_res)
        }, error = function(e) {
          cat("Error for repeat number:", idx, "\n")
          message(e)
        })
        
      }
      
      group_wise_SVM_CV_weighting <- do.call(plyr::rbind.fill, group_wise_SVM_CV_weighting_list)
      saveRDS(group_wise_SVM_CV_weighting, file=paste0(rdata_path,
                                                       sprintf("%s_wise_CV_linear_SVM_%s_%s.Rds",
                                                               grouping_var,
                                                               pairwise_feature_set,
                                                               weighting_name)))
    }, error = function(e) {
      cat("Could not run linear SVM for", grouping_var, pairwise_feature_set, ".\n")
      message(e)
    })
  } else {
    group_wise_SVM_CV_weighting <- readRDS(paste0(rdata_path, 
                                                  sprintf("%s_wise_CV_linear_SVM_%s_%s.Rds",
                                                          grouping_var,
                                                          pairwise_feature_set, 
                                                          weighting_name)))
  }
  
  #### Calculate balanced accuracy across all folds
  if (!file.exists(paste0(rdata_path, sprintf("%s_wise_CV_linear_SVM_%s_%s_balacc.Rds",
                                              grouping_var, 
                                              pairwise_feature_set, 
                                              weighting_name)))) {
    group_wise_SVM_balanced_accuracy <- group_wise_SVM_CV_weighting %>%
      group_by(grouping_var, Noise_Proc, Sample_Type, repeat_number) %>%
      summarise(accuracy = sum(Prediction_Correct) / n(),
                balanced_accuracy = caret::confusionMatrix(data = Predicted_Diagnosis,
                                                           reference = Actual_Diagnosis)$byClass[["Balanced Accuracy"]]) %>%
      group_by(grouping_var, Noise_Proc, Sample_Type) %>%
      summarise(mean_accuracy = mean(accuracy, na.rm=T),
                mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T)) %>%
      dplyr::rename("accuracy" = "mean_accuracy",
                    "balanced_accuracy" = "mean_balanced_accuracy")
    
    saveRDS(group_wise_SVM_balanced_accuracy, file=paste0(rdata_path, sprintf("%s_wise_CV_linear_SVM_%s_%s_balacc.Rds",
                                                                              grouping_var, 
                                                                              pairwise_feature_set, 
                                                                              weighting_name)))
  } else {
    group_wise_SVM_balanced_accuracy <- readRDS(paste0(rdata_path, sprintf("%s_wise_CV_linear_SVM_%s_%s_balacc.Rds",
                                                                           grouping_var, 
                                                                           pairwise_feature_set, 
                                                                           weighting_name)))
  }
  
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
  
  run_number = ifelse(is.null(run_number), "", run_number)
  
  # Where to save PBS script to
  output_scripts_dir <- paste0(github_dir, sprintf("fMRI_FeaturesDisorders/classification_analysis/pairwise_analysis/null_pbs_scripts/%s_%s_wise_%s_%s_null_model_fits%s/",
                                                   dataset_ID,
                                                   grouping_var,
                                                   pairwise_feature_set,
                                                   weighting_name,
                                                   run_number))
  
  cat("\nNow generating null PBS scripts for", grouping_var, "\n")
  cat("Script location:", output_scripts_dir, "\n")
  
  # Make these directories
  icesTAF::mkdir(output_data_dir)
  icesTAF::mkdir(output_scripts_dir)
  
  
  # Lookup table for PBS script
  lookup_list <- list("NAME" = sprintf("pairwise_%s_wise_null_model_fit%s",
                                       grouping_var, run_number),
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
