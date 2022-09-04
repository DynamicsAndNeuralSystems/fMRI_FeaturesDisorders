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

# ###############################################################################
# Load data
# ###############################################################################
noise_label = gsub("\\+", "_", noise_proc_for_null)

univariate_data_file <- paste0(rdata_path, sprintf("%s_%s_filtered_zscored.Rds",
                                                   dataset_ID, univariate_feature_set))
univariate_feature_data <- readRDS(univariate_data_file) %>%
  dplyr::filter(Noise_Proc == noise_proc_for_null)


pairwise_data_file <- paste0(rdata_path, sprintf("%s_%s_filtered_zscored.Rds",
                                                 dataset_ID, pairwise_feature_set))
pairwise_feature_data <- readRDS(pairwise_data_file) %>%
  dplyr::filter(Noise_Proc == noise_proc_for_null)

SPI_directionality_file <- paste0(github_dir, "fMRI_FeaturesDisorders/classification_analysis/pairwise_analysis/SPI_Direction_Info.csv")
SPI_directionality <- read.csv(SPI_directionality_file)

# Filter subjects to only those with data available for both
univariate_feature_data <- univariate_feature_data %>% semi_join(subjects_to_use)
pairwise_feature_data <- pairwise_feature_data %>% semi_join(subjects_to_use)

################################################################################
# Define weighting parameters
################################################################################

weighting_param_df <- data.frame(name = c("inv_prob"),
                                 use_inv_prob_weighting = c(TRUE))

grouping_param_df <- data.frame(grouping_var = "SPI",
                                SVM_feature_var = "region_pair")

################################################################################
# All catch22 features with each SPI individually
################################################################################

#### 10-fold linear SVM with different weights
# Iterate over weighting_param_df
for (i in 1:nrow(grouping_param_df)) {
  grouping_var = grouping_param_df$grouping_var[i]
  SVM_feature_var = grouping_param_df$SVM_feature_var[i]
  
  #### 10-fold linear SVM with different weights
  # Iterate over weighting_param_df 
  for (j in 1:nrow(weighting_param_df)) {
    weighting_name <- weighting_param_df$name[j]
    use_inv_prob_weighting <- weighting_param_df$use_inv_prob_weighting[j]
    
    # Run given weighting for 10-fold CV linear SVM
    if (!file.exists(paste0(rdata_path, sprintf("univariate_%s_pairwise_%s_CV_linear_SVM_%s.Rds",
                                                univariate_feature_set, pairwise_feature_set, weighting_name)))) {
      tryCatch({univariate_pairwise_SVM_CV_weighting <- run_combined_uni_pairwise_cv_svm_by_input_var(dataset_ID = dataset_ID,
                                                                                                      data_path = data_path,
                                                                                                      rdata_path = rdata_path,
                                                                                                      univariate_data = univariate_feature_data,
                                                                                                      univariate_feature_set = univariate_feature_set,
                                                                                                      pairwise_data = pairwise_feature_data,
                                                                                                      pairwise_feature_set = pairwise_feature_set,
                                                                                                      SPI_directionality = SPI_directionality,
                                                                                                      flds = sample_folds,
                                                                                                      num_k_folds = 10,
                                                                                                      svm_kernel = kernel,
                                                                                                      noise_proc = noise_proc_for_null,
                                                                                                      out_of_sample_only = TRUE,
                                                                                                      use_inv_prob_weighting = use_inv_prob_weighting,
                                                                                                      shuffle_labels = FALSE)
      saveRDS(univariate_pairwise_SVM_CV_weighting, file=paste0(rdata_path,
                                                                sprintf("univariate_%s_pairwise_%s_CV_linear_SVM_%s.Rds",
                                                                        univariate_feature_set, pairwise_feature_set, weighting_name)))
      }, error = function(e) {
        message(e)
      })
    } else {
      univariate_pairwise_SVM_CV_weighting <- readRDS(paste0(rdata_path,
                                                             sprintf("univariate_%s_pairwise_%s_CV_linear_SVM_%s.Rds",
                                                                     univariate_feature_set, pairwise_feature_set, weighting_name)))
    }
    
    #### Calculate balanced accuracy across all folds
    if (!file.exists(paste0(rdata_path, sprintf("univariate_%s_pairwise_%s_CV_linear_SVM_%s_balacc.Rds",
                                                univariate_feature_set, 
                                                pairwise_feature_set,
                                                weighting_name)))) {
      univariate_pairwise_SVM_balanced_accuracy <- univariate_pairwise_SVM_CV_weighting %>%
        group_by(SPI, Noise_Proc, Sample_Type) %>%
        summarise(accuracy = sum(Prediction_Correct) / n(),
                  balanced_accuracy = caret::confusionMatrix(data = Predicted_Diagnosis,
                                                             reference = Actual_Diagnosis)$byClass[["Balanced Accuracy"]])
      
      saveRDS(univariate_pairwise_SVM_balanced_accuracy, file=paste0(rdata_path, sprintf("univariate_%s_pairwise_%s_CV_linear_SVM_%s_balacc.Rds",
                                                                                         univariate_feature_set, 
                                                                                         pairwise_feature_set,
                                                                                         weighting_name)))
    } else {
      univariate_pairwise_SVM_balanced_accuracy <- readRDS(paste0(rdata_path, sprintf("univariate_%s_pairwise_%s_CV_linear_SVM_%s_balacc.Rds",
                                                                                      univariate_feature_set, 
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
    # Define the univariate+pairwise combined template PBS script
    template_pbs_file <- paste0(github_dir, "fMRI_FeaturesDisorders/helper_functions/classification/template_combined_univariate_pairwise_null_model_fit.pbs")
    
    run_number = ifelse(is.null(run_number), "", run_number)
    
    # Where to store null model fit results
    output_data_dir <- paste0(rdata_path, sprintf("%s_univariate_%s_pairwise_%s_%s_null_model_fits%s/",
                                                  dataset_ID,
                                                  univariate_feature_set,
                                                  pairwise_feature_set,
                                                  weighting_name,
                                                  run_number))
    
    
    
    # Where to save PBS script to
    output_scripts_dir <- paste0(github_dir, sprintf("fMRI_FeaturesDisorders/classification_analysis/combined_univariate_pairwise/null_pbs_scripts/%s_univariate_%s_pairwise_%s_%s_null_model_fits%s/",
                                                     dataset_ID,
                                                     univariate_feature_set,
                                                     pairwise_feature_set,
                                                     weighting_name,
                                                     run_number))
    
    cat("\nNow generating null PBS scripts for.\n")
    cat("Script location:", output_scripts_dir, "\n")
    
    # Make these directories
    icesTAF::mkdir(output_data_dir)
    icesTAF::mkdir(output_scripts_dir)
    
    
    # Lookup table for PBS script
    lookup_list <- list("NAME" = sprintf("univariate_%s_pairwise_%s_null_model_fit%s",
                                         univariate_feature_set, 
                                         pairwise_feature_set,
                                         run_number),
                        "MEMNUM" = "20",
                        "NCPUS" = "1",
                        "DATASET_ID" = dataset_ID,
                        "GITHUB_DIR" = github_dir,
                        "DATA_PATH" = data_path,
                        "RDATA_PATH" = rdata_path,
                        "UNIVARIATE_DATA_FILE" = univariate_data_file,
                        "UNIVARIATE_FEATURE_SET" = univariate_feature_set,
                        "PAIRWISE_DATA_FILE" = pairwise_data_file,
                        "PAIRWISE_FEATURE_SET" = pairwise_feature_set,
                        "SPI_DIRECTIONALITY_FILE" = SPI_directionality_file,
                        "EMAIL" = email,
                        "PBS_NOTIFY" = "a",
                        "WALL_HRS" = wall_hrs,
                        "NOISE_PROCS" = noise_proc_for_null,
                        "NUM_K_FOLDS" = num_k_folds,
                        "NUM_PERMS_PER_ITER" = nperm_per_iter,
                        "OUTPUT_DATA_DIR" = output_data_dir,
                        "SAMPLE_METADATA_FILE" = sample_metadata_file,
                        "GROUPING_VAR" = grouping_var,
                        "SVM_FEATURE_VAR" = SVM_feature_var,
                        "WEIGHTING_NAME" = weighting_name)
    
    to_be_replaced <- names(lookup_list)
    replacement_values <- unlist(unname(lookup_list))
    
    # Create a PBS script per iteration
    for (p in 1:num_permutations) {
      
      # Run command if null file doesn't exist
      if (!file.exists(sprintf("%s/univariate_%s_pairwise_%s_CV_linear_SVM_%s_null_model_fit_iter_%s.Rds",
                               output_data_dir, univariate_feature_set, 
                               pairwise_feature_set, weighting_name, p))) {
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
  
  saveRDS(combined_uni_pairwise_PCA_linear_SVM_res, file = paste0(rdata_path, sprintf("Univariate_%s_Pairwise_%s_PCA_CV_linear_SVM_%s.Rds",
                                                                                      univariate_feature_set, pairwise_feature_set,
                                                                                      weighting_name)))
} else {
  combined_uni_pairwise_PCA_linear_SVM_res <- readRDS(file=paste0(rdata_path, sprintf("Univariate_%s_Pairwise_%s_PCA_CV_linear_SVM_%s.Rds",
                                                                                      univariate_feature_set, pairwise_feature_set,
                                                                                      weighting_name)))
}