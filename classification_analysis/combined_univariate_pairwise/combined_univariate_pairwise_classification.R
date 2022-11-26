#-------------------------------------------------------------------------------
# Define paths
#-------------------------------------------------------------------------------
# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--github_dir", default="/headnode1/abry4213/github/")
parser$add_argument("--data_path", default="/headnode1/abry4213/data/UCLA_Schizophrenia/")
parser$add_argument("--sample_metadata_file", default="UCLA_Schizophrenia_sample_metadata.Rds")
parser$add_argument("--pairwise_feature_set", default="pyspi14_corrected")
parser$add_argument("--univariate_feature_set", default="catch22")
parser$add_argument("--noise_procs", default=c(""))
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
noise_procs <- args$noise_procs
noise_proc_for_null <- args$noise_proc_for_null
dataset_ID <- args$dataset_ID
sample_metadata_file <- args$sample_metadata_file
email <- args$email
run_number <- args$run_number

# univariate_feature_set <- "catch22"
# pairwise_feature_set <- "pyspi14_corrected"
# github_dir <- "~/github/"
# email <- "abry4213@uni.sydney.edu.au"

# UCLA schizophrenia
# data_path <- "~/data/UCLA_Schizophrenia/"
# dataset_ID <- "UCLA_Schizophrenia"
# sample_metadata_file <- "UCLA_Schizophrenia_sample_metadata.Rds"
# noise_procs <- "AROMA+2P;AROMA+2P+GMR;AROMA+2P+DiCER"
# noise_procs <- "AROMA+2P+GMR"
# noise_proc_for_null <- "AROMA+2P+GMR"

# ABIDE ASD
# data_path <- "~/data/ABIDE_ASD/"
# sample_metadata_file <- "ABIDE_ASD_sample_metadata.Rds"
# dataset_ID <- "ABIDE_ASD"
# noise_procs <- c("FC1000")
# noise_proc_for_null <- "FC1000"

rdata_path <- paste0(data_path, "processed_data/Rdata/")
plot_dir <- paste0(data_path, "plots/")

TAF::mkdir(plot_dir)

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
# Use a linear kernel
kernel = "linear"
weighting_name <- "inv_prob"
use_inv_prob_weighting <- TRUE

################################################################################
# All catch22 features with each SPI individually
################################################################################

#### 10-fold linear SVM with different weights
# Run given weighting for 10-fold CV linear SVM
if (!file.exists(paste0(rdata_path, sprintf("univariate_%s_pairwise_%s_CV_linear_SVM_%s.Rds",
                                            univariate_feature_set, pairwise_feature_set, weighting_name)))) {
  tryCatch({
    univariate_pairwise_SVM_CV_weighting_list <- list()
    for (idx in 1:length(sample_folds)) {
      repeat_res <- run_combined_uni_pairwise_cv_svm_by_input_var(dataset_ID = dataset_ID,
                                                                  data_path = data_path,
                                                                  rdata_path = rdata_path,
                                                                  univariate_data = univariate_feature_data,
                                                                  univariate_feature_set = univariate_feature_set,
                                                                  pairwise_data = pairwise_feature_data,
                                                                  pairwise_feature_set = pairwise_feature_set,
                                                                  SPI_directionality = SPI_directionality,
                                                                  flds = sample_folds[[idx]],
                                                                  repeat_number = idx,
                                                                  num_k_folds = 10,
                                                                  svm_kernel = kernel,
                                                                  noise_proc = noise_proc_for_null,
                                                                  out_of_sample_only = TRUE,
                                                                  use_inv_prob_weighting = use_inv_prob_weighting,
                                                                  shuffle_labels = FALSE,
                                                                  drop_NaN = FALSE,
                                                                  impute_NaN = TRUE)
      
      univariate_pairwise_SVM_CV_weighting_list <- list.append(univariate_pairwise_SVM_CV_weighting_list,
                                                               repeat_res)
      
    }
    univariate_pairwise_SVM_CV_weighting <- do.call(plyr::rbind.fill, 
                                                    univariate_pairwise_SVM_CV_weighting_list)
    saveRDS(univariate_pairwise_SVM_CV_weighting, file=paste0(rdata_path,
                                                              sprintf("univariate_%s_pairwise_%s_CV_linear_SVM_%s.Rds",
                                                                      univariate_feature_set, pairwise_feature_set, weighting_name)))
  }, error = function(e) {
    cat("Could not run linear SVM for", pairwise_feature_set, ".\n")
    message(e)
  })
} else {
  univariate_pairwise_SVM_CV_weighting <- readRDS(paste0(rdata_path,
                                                         sprintf("univariate_%s_pairwise_%s_CV_linear_SVM_%s.Rds",
                                                                 univariate_feature_set, pairwise_feature_set, weighting_name)))
}

#### Calculate balanced accuracy across all folds
if (!file.exists(paste0(rdata_path, sprintf("univariate_%s_pairwise_%s_CV_linear_SVM_%s_balacc_across_repeats.Rds",
                                            univariate_feature_set, 
                                            pairwise_feature_set,
                                            weighting_name)))) {
  
  # First find balanced accuracy per repeat across folds
  univariate_pairwise_SVM_balanced_accuracy <- univariate_pairwise_SVM_CV_weighting %>%
    group_by(SPI, Noise_Proc, Sample_Type, fold_number, repeat_number) %>%
    summarise(accuracy = sum(Prediction_Correct) / n(),
              balanced_accuracy = caret::confusionMatrix(data = Predicted_Diagnosis,
                                                         reference = Actual_Diagnosis)$byClass[["Balanced Accuracy"]])
  
  saveRDS(univariate_pairwise_SVM_balanced_accuracy, file=paste0(rdata_path, sprintf("univariate_%s_pairwise_%s_CV_linear_SVM_%s_balacc.Rds",
                                                                                     univariate_feature_set, 
                                                                                     pairwise_feature_set,
                                                                                     weighting_name)))
  
  # Then find averaged balanced accuracy across all repeats
  univariate_pairwise_SVM_balanced_accuracy_across_repeats <- univariate_pairwise_SVM_balanced_accuracy %>%
    group_by(SPI, Noise_Proc, Sample_Type) %>%
    summarise(mean_accuracy = mean(accuracy, na.rm=T),
              mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T)) %>%
    dplyr::rename("accuracy" = "mean_accuracy",
                  "balanced_accuracy" = "mean_balanced_accuracy")
  
  saveRDS(univariate_pairwise_SVM_balanced_accuracy_across_repeats, file=paste0(rdata_path, sprintf("univariate_%s_pairwise_%s_CV_linear_SVM_%s_balacc_across_repeats.Rds",
                                                                                                    univariate_feature_set, 
                                                                                                    pairwise_feature_set,
                                                                                                    weighting_name)))
  
  
} else {
  univariate_pairwise_SVM_balanced_accuracy_across_repeats <- readRDS(paste0(rdata_path, sprintf("univariate_%s_pairwise_%s_CV_linear_SVM_%s_balacc_across_repeats.Rds",
                                                                                                 univariate_feature_set, 
                                                                                                 pairwise_feature_set,
                                                                                                 weighting_name)))
}

