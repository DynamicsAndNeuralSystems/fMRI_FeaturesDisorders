# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

# Data directories
parser$add_argument("--pairwise_data_file", default = "")
parser$add_argument("--univariate_data_file", default = "")
parser$add_argument("--dataset_ID", default = "UCLA_Schizophrenia")
parser$add_argument("--sample_metadata_file", default="UCLA_Schizophrenia_sample_metadata.Rds")
parser$add_argument("--SPI_directionality_file", default="/headnode1/abry4213/github/fMRI_FeaturesDisorders/pairwise_analysis/SPI_Direction_Info.csv")
parser$add_argument("--data_path", default="/headnode1/abry4213/data/UCLA_Schizophrenia/")
parser$add_argument("--rdata_path", default="/headnode1/abry4213/data/UCLA_Schizophrenia/processed_data/Rdata/")
parser$add_argument("--output_data_dir", default="/headnode1/abry4213/data/UCLA_Schizophrenia/Rdata/Pairwise_pyspi_19_inv_prob_null_model_fits/")
parser$add_argument("--github_dir", default="/headnode1/abry4213/github/")

# Permutation arguments
parser$add_argument("--num_k_folds", default=10)
parser$add_argument("--null_iter_number", default=1)
parser$add_argument("--num_perms_for_iter", default=1)
parser$add_argument("--svm_kernel", default="linear")

# Feature sets
parser$add_argument("--univariate_feature_set", default="catch22")
parser$add_argument("--pairwise_feature_set", default="pyspi_19")
parser$add_argument("--grouping_var", default="SPI")
parser$add_argument("--svm_feature_var", default="region_pair")
parser$add_argument("--noise_proc", default="AROMA+2P+GMR")

# Output metrics
parser$add_argument("--return_all_fold_metrics", action='store_true', default=FALSE)
parser$add_argument("--weighting_name", default="unweighted")
parser$add_argument("--use_inv_prob_weighting", action='store_true', default=FALSE)

# Indicate whether analysis is univariate, pariwise, or both
parser$add_argument("--univariate", action='store_true', default=FALSE)
parser$add_argument("--pairwise", action='store_true', default=FALSE)
parser$add_argument("--combined_univariate_pairwise", action='store_true', default=FALSE)

# Parse input arguments
args <- parser$parse_args()

# Data directories
dataset_ID <- args$dataset_ID
sample_metadata_file <- args$sample_metadata_file
pairwise_data_file <- args$pairwise_data_file
univariate_data_file <- args$univariate_data_file
SPI_directionality_file <- args$SPI_directionality_file
data_path <- args$data_path
rdata_path <- args$rdata_path
output_data_dir <- args$output_data_dir
github_dir <- args$github_dir

# Permutation arguments
num_k_folds <- as.numeric(args$num_k_folds)
null_iter_number <- args$null_iter_number
num_perms_for_iter <- args$num_perms_for_iter
svm_kernel <- args$svm_kernel

# Feature sets
grouping_var <- args$grouping_var
svm_feature_var <- args$svm_feature_var
noise_proc <- args$noise_proc
pairwise_feature_set <- args$pairwise_feature_set
univariate_feature_set <- args$univariate_feature_set

# Output metrics
return_all_fold_metrics <- args$return_all_fold_metrics
weighting_name <- args$weighting_name
use_inv_prob_weighting <- args$use_inv_prob_weighting

# Indicate whether analysis is univariate, pariwise, or both
univariate <- args$univariate
pairwise <- args$pairwise
combined_univariate_pairwise <- args$combined_univariate_pairwise


################### Debugging ######################
# univariate_feature_set <- "catch22"
# pairwise_feature_set <- "pyspi14"
# github_dir <- "/headnode1/abry4213/github/"
# email <- "abry4213@uni.sydney.edu.au"
# num_k_folds <- 10
# null_iter_number <- 1
# num_perms_for_iter <- 5
# svm_kernel <- "linear"
# SPI_directionality_file <- paste0(github_dir, "fMRI_FeaturesDisorders/classification_analysis/pairwise_analysis/SPI_Direction_Info.csv")

# UCLA schizophrenia
# dataset_ID <- "UCLA_Schizophrenia"
# sample_metadata_file <- "UCLA_Schizophrenia_sample_metadata.Rds"
# data_path <- "/headnode1/abry4213/data/UCLA_Schizophrenia/"
# noise_proc <- "AROMA+2P+GMR"

# ABIDE ASD
# dataset_ID <- "ABIDE_ASD"
# sample_metadata_file <- "ABIDE_ASD_sample_metadata.Rds"
# data_path <- "/headnode1/abry4213/data/ABIDE_ASD/"
# noise_proc <- "FC1000"

# All
# use_inv_prob_weighting = TRUE
# rdata_path <- paste0(data_path, "processed_data/Rdata/")
# output_data_dir <- paste0(rdata_path, dataset_ID, "_univariate_",
#                           univariate_feature_set, "_pairwise_",
#                           pairwise_feature_set, "_inv_prob_null_model_fits/")
# univariate_data_file <- paste0(rdata_path, sprintf("%s_%s_filtered_zscored.Rds",
#                                                    dataset_ID, univariate_feature_set))
# pairwise_data_file <- paste0(rdata_path, sprintf("%s_%s_filtered_zscored.Rds",
#                                                  dataset_ID, pairwise_feature_set))

####################################################

if (is.null(rdata_path)) {
  rdata_path <- paste0(data_path, "processed_data/Rdata/")
}

# Load sample metadata
sample_metadata <- readRDS(paste0(data_path, sample_metadata_file))

# Source linear SVM functions
source(paste0(github_dir, "fMRI_FeaturesDisorders/helper_functions/classification/Linear_SVM.R"))

icesTAF::mkdir(output_data_dir)

cat("\nNumber of k-folds:", num_k_folds, "\n")
cat("\nNum permutations per iteration:", num_perms_for_iter, "\n")

if (univariate & !file.exists(sprintf("%s/%s_wise_%s_%s_null_model_fit_iter_%s.Rds",
                                      output_data_dir, grouping_var, univariate_feature_set, 
                                      weighting_name, null_iter_number))) {
  # Run null iteration
  null_out <- 1:num_perms_for_iter  %>%
    purrr::map_df( ~ run_univariate_cv_svm_by_input_var(data_path = data_path,
                                                        rdata_path = rdata_path,
                                                        sample_metadata = sample_metadata,
                                                        dataset_ID = dataset_ID,
                                                        svm_kernel = svm_kernel,
                                                        pairwise_feature_set = pairwise_feature_set,
                                                        univariate_feature_set = univariate_feature_set,
                                                        grouping_var = grouping_var,
                                                        svm_feature_var = svm_feature_var,
                                                        noise_procs = noise_proc,
                                                        num_k_folds = num_k_folds,
                                                        out_of_sample_only = TRUE,
                                                        use_inv_prob_weighting = use_inv_prob_weighting,
                                                        shuffle_labels = T) %>%
                     # Keep track of which null iteration this is
                     mutate(Null_Iter_Number = .x + (.x * (as.numeric(null_iter_number) - 1))))
  
  null_out <- null_out %>% 
    group_by(grouping_var, Noise_Proc, Sample_Type, Null_Iter_Number) %>%
    summarise(accuracy = sum(Prediction_Correct) / n(),
              balanced_accuracy = caret::confusionMatrix(data = Predicted_Diagnosis,
                                                         reference = Actual_Diagnosis)$byClass[["Balanced Accuracy"]])
  
  # Save null results to RDS
  saveRDS(null_out, file=sprintf("%s/%s_wise_%s_%s_null_model_fit_iter_%s.Rds",
                                 output_data_dir, grouping_var, univariate_feature_set, 
                                 weighting_name, null_iter_number))
}

if (pairwise & !(file.exists(sprintf("%s/%s_wise_%s_%s_null_model_fit_iter_%s.Rds",
                                 output_data_dir, grouping_var, pairwise_feature_set, 
                                 weighting_name, null_iter_number)))) {
  # Load data
  pairwise_data <- readRDS(pairwise_data_file)
  SPI_directionality <- read.csv(SPI_directionality_file)
  
  
  cat("\nHead of pairwise data:\n")
  print(head(pairwise_data))
  cat("\nHead of SPI directionality data:\n")
  print(head(SPI_directionality))
  
  # Run null iteration
  null_out <- 1:num_perms_for_iter %>%
    purrr::map_df( ~ run_pairwise_cv_svm_by_input_var(data_path = data_path,
                                                      rdata_path = rdata_path,
                                                      dataset_ID = dataset_ID,
                                                      sample_metadata = sample_metadata,
                                                      pairwise_data = pairwise_data,
                                                      SPI_directionality = SPI_directionality,
                                                      svm_kernel = svm_kernel,
                                                      grouping_var = grouping_var,
                                                      svm_feature_var = svm_feature_var,
                                                      noise_proc = noise_proc,
                                                      num_k_folds = num_k_folds,
                                                      out_of_sample_only = TRUE,
                                                      use_inv_prob_weighting = use_inv_prob_weighting,
                                                      shuffle_labels = T) %>%
                     # Keep track of which null iteration this is
                     mutate(Null_Iter_Number = .x + (.x * (as.numeric(null_iter_number) - 1))))
  
  # Save null results to RDS
  saveRDS(null_out, file=sprintf("%s/%s_wise_%s_%s_null_model_fit_iter_%s.Rds",
                                 output_data_dir, grouping_var, pairwise_feature_set, 
                                 weighting_name, null_iter_number))
}

if (combined_univariate_pairwise & !file.exists(sprintf("%s/univariate_%s_pairwise_%s_CV_linear_SVM_%s_null_model_fit_iter_%s.Rds",
                                                        output_data_dir, 
                                                        univariate_feature_set, 
                                                        pairwise_feature_set, 
                                                        weighting_name, 
                                                        null_iter_number))) {
  univariate_feature_data <- readRDS(univariate_data_file) %>%
    dplyr::filter(Noise_Proc == noise_proc)
    
  pairwise_feature_data <- readRDS(pairwise_data_file) %>%
    dplyr::filter(Noise_Proc == noise_proc)
  
  SPI_directionality <- read.csv(SPI_directionality_file)
  
  # Run null iteration
  null_out <- 1:num_perms_for_iter %>%
    purrr::map_df( ~ run_combined_uni_pairwise_cv_svm_by_input_var(dataset_ID = dataset_ID,
                                                                   data_path = data_path,
                                                                   rdata_path = rdata_path,
                                                                   univariate_data = univariate_feature_data,
                                                                   univariate_feature_set = univariate_feature_set,
                                                                   pairwise_data = pairwise_feature_data,
                                                                   pairwise_feature_set = pairwise_feature_set,
                                                                   SPI_directionality = SPI_directionality,
                                                                   num_k_folds = 10,
                                                                   svm_kernel = svm_kernel,
                                                                   noise_proc = noise_proc,
                                                                   out_of_sample_only = TRUE,
                                                                   use_inv_prob_weighting = use_inv_prob_weighting,
                                                                   shuffle_labels = TRUE) %>%
                     # Keep track of which null iteration this is
                     mutate(Null_Iter_Number = .x + (.x * (as.numeric(null_iter_number) - 1))))
  
  # Save null results to RDS
  saveRDS(null_out, file=sprintf("%s/univariate_%s_pairwise_%s_CV_linear_SVM_%s_null_model_fit_iter_%s.Rds",
                                 output_data_dir, univariate_feature_set, 
                                 pairwise_feature_set, weighting_name, null_iter_number))
}

