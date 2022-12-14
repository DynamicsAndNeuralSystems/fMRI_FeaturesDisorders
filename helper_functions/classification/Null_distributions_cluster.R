# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

# Data directories
parser$add_argument("--pairwise_data_file", default = "")
parser$add_argument("--univariate_data_file", default = "")
parser$add_argument("--dataset_ID", default = "UCLA_CNP_ABIDE_ASD")
parser$add_argument("--sample_metadata_file", default="UCLA_CNP_ABIDE_ASD_sample_metadata.Rds")
parser$add_argument("--SPI_directionality_file", default="/headnode1/abry4213/github/fMRI_FeaturesDisorders/pairwise_analysis/SPI_Direction_Info.csv")
parser$add_argument("--data_path", default="/headnode1/abry4213/data/UCLA_CNP_ABIDE_ASD/")
parser$add_argument("--rdata_path", default="/headnode1/abry4213/data/UCLA_CNP_ABIDE_ASD/processed_data/Rdata/")
parser$add_argument("--output_data_dir", default="/headnode1/abry4213/data/UCLA_CNP_ABIDE_ASD/Rdata/Pairwise_pyspi_19_inv_prob_null_model_fits/")
parser$add_argument("--github_dir", default="/headnode1/abry4213/github/")

# Permutation arguments
parser$add_argument("--num_k_folds", default=10)
parser$add_argument("--null_iter_number", default=1)
parser$add_argument("--num_perms_for_iter", default=1)
parser$add_argument("--svm_kernel", default="linear")

# Feature sets
parser$add_argument("--univariate_feature_set", default="catch22")
parser$add_argument("--pairwise_feature_set", default="pyspi14")
parser$add_argument("--grouping_var", default="SPI")
parser$add_argument("--svm_feature_var", default="region_pair")

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

# Load sample metadata
sample_metadata <- readRDS(paste0(data_path, sample_metadata_file))

# Load contrasts to control participants
study_comparisons_to_control <- readRDS(paste0(rdata_path, dataset_ID, "_comparisons_to_control.Rds"))


# Source linear SVM functions
source(paste0(github_dir, "fMRI_FeaturesDisorders/helper_functions/classification/Linear_SVM.R"))

TAF::mkdir(output_data_dir)

cat("\nNumber of k-folds:", num_k_folds, "\n")
cat("\nNum permutations per iteration:", num_perms_for_iter, "\n")

if (univariate & !file.exists(sprintf("%s/%s_%s_wise_%s_%s_null_model_fit_iter_%s.Rds",
                                      output_data_dir, dataset_ID,
                                      grouping_var, univariate_feature_set, 
                                      weighting_name, null_iter_number))) {
  # Load univariate features
  univariate_features <- readRDS(paste0(rdata_path, dataset_ID, "_", univariate_feature_set,
                                               "_and_catch2_filtered_zscored.Rds")) %>%
    left_join(., sample_metadata) 
  
  univariate_null_res_list <- list()
  
  for (i in 1:nrow(study_comparisons_to_control)) {
    study <- study_comparisons_to_control$Study[i]
    group_to_compare <- study_comparisons_to_control$Group_to_Compare[i] 
    
    # Get subjects to compare
    sample_groups <- sample_metadata %>%
      filter(Sample_ID %in% subjects_to_use,
             Study == study,
             Diagnosis %in% c("Control", group_to_compare)) %>%
      distinct(Sample_ID, Diagnosis)
    
    # Combine control data with diagnosis group data
    group_data_for_SVM <- univariate_features %>%
      filter(Study == study,
             Diagnosis %in% c(group_to_compare, "Control")) 
    
    # Filter to univariate feature set
    group_data_for_SVM_feature_set <- group_data_for_SVM %>%
      filter(feature_set == univariate_feature_set)
    
    # Run null iteration
    null_out <- 1:num_perms_for_iter  %>%
      purrr::map_df( ~ run_univariate_cv_svm_by_input_var(feature_matrix = group_data_for_SVM_feature_set,
                                                          dataset_ID = dataset_ID,
                                                          svm_kernel = svm_kernel,
                                                          univariate_feature_set = univariate_feature_set,
                                                          sample_groups = sample_groups,
                                                          grouping_var = grouping_var,
                                                          svm_feature_var = svm_feature_var,
                                                          num_k_folds = num_k_folds,
                                                          out_of_sample_only = TRUE,
                                                          use_inv_prob_weighting = use_inv_prob_weighting,
                                                          shuffle_labels = T) %>%
                       # Keep track of which null iteration this is
                       mutate(Null_Iter_Number = .x + (.x * (as.numeric(null_iter_number) - 1)))) 
    
    # Add study and group info
    null_out$Study <- study
    null_out$Group_to_Compare <- group_to_compare
    
    univariate_null_res_list <- list.append(univariate_null_res_list, null_out)
  }
  
  # Concatenate results
  univariate_null_res <- do.call(plyr::rbind.fill, univariate_null_res_list) %>% 
    group_by(grouping_var, Study, Group_to_Compare, Sample_Type, fold_number, Null_Iter_Number) %>%
    # First find accuracy and balanced accuracy by fold
    summarise(accuracy = sum(Prediction_Correct) / n(),
              balanced_accuracy = caret::confusionMatrix(data = Predicted_Diagnosis,
                                                         reference = Actual_Diagnosis)$byClass[["Balanced Accuracy"]]) %>%
    # Then take average acc/balacc across all ten folds per iteration
    group_by(grouping_var, Study, Group_to_Compare, Sample_Type, Null_Iter_Number) %>%
    summarise(mean_accuracy = mean(accuracy, na.rm=T),
              mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T)) %>%
    dplyr::rename("accuracy" = "mean_accuracy",
                  "balanced_accuracy" = "mean_balanced_accuracy")
  
  # Save null results to RDS
  saveRDS(null_out, file=sprintf("%s/%s_%s_wise_%s_%s_null_model_fit_iter_%s.Rds",
                                 output_data_dir, dataset_ID, grouping_var, univariate_feature_set, 
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

  
  null_out <- null_out %>% 
    group_by(grouping_var, Noise_Proc, Sample_Type, fold_number, Null_Iter_Number) %>%
    # First find accuracy and balanced accuracy by fold
    summarise(accuracy = sum(Prediction_Correct) / n(),
              balanced_accuracy = caret::confusionMatrix(data = Predicted_Diagnosis,
                                                         reference = Actual_Diagnosis)$byClass[["Balanced Accuracy"]]) %>%
    # Then take average acc/balacc across all ten folds per iteration
    group_by(grouping_var, Noise_Proc, Sample_Type, Null_Iter_Number) %>%
    summarise(mean_accuracy = mean(accuracy, na.rm=T),
              mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T)) %>%
    dplyr::rename("accuracy" = "mean_accuracy",
                  "balanced_accuracy" = "mean_balanced_accuracy")
  
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
  
  null_out <- null_out %>% 
    group_by(SPI, Noise_Proc, Sample_Type, fold_number, Null_Iter_Number) %>%
    # First find accuracy and balanced accuracy by fold
    summarise(accuracy = sum(Prediction_Correct) / n(),
              balanced_accuracy = caret::confusionMatrix(data = Predicted_Diagnosis,
                                                         reference = Actual_Diagnosis)$byClass[["Balanced Accuracy"]]) %>%
    # Then take average acc/balacc across all ten folds per iteration
    group_by(SPI, Noise_Proc, Sample_Type, Null_Iter_Number) %>%
    summarise(mean_accuracy = mean(accuracy, na.rm=T),
              mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T)) %>%
    dplyr::rename("accuracy" = "mean_accuracy",
                  "balanced_accuracy" = "mean_balanced_accuracy")
  
  # Save null results to RDS
  saveRDS(null_out, file=sprintf("%s/univariate_%s_pairwise_%s_CV_linear_SVM_%s_null_model_fit_iter_%s.Rds",
                                 output_data_dir, univariate_feature_set, 
                                 pairwise_feature_set, weighting_name, null_iter_number))
}

