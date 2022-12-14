#-------------------------------------------------------------------------------
# Define paths
#-------------------------------------------------------------------------------
# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--github_dir", default="/headnode1/abry4213/github/")
parser$add_argument("--data_path", default="/headnode1/abry4213/data/UCLA_CNP_ABIDE_ASD/")
parser$add_argument("--sample_metadata_file", default="UCLA_CNP_ABIDE_ASD_sample_metadata.Rds")
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
# 
# univariate_feature_set <- "catch22"
# pairwise_feature_set <- "pyspi14"
# github_dir <- "~/github/"
# email <- "abry4213@uni.sydney.edu.au"
# add_catch2 <- TRUE
# data_path <- "~/data/UCLA_CNP_ABIDE_ASD/"
# dataset_ID <- "UCLA_CNP_ABIDE_ASD"
# sample_metadata_file <- "UCLA_CNP_ABIDE_ASD_sample_metadata.Rds"

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

subjects_to_use <- readRDS(paste0(rdata_path, sprintf("%s_samples_with_univariate_%s_and_pairwise_%s_filtered.Rds",
                                                      dataset_ID,
                                                      univariate_feature_set,
                                                      pairwise_feature_set)))

# Load sample metadata
sample_metadata <- readRDS(paste0(data_path, "study_metadata/", sample_metadata_file)) %>%
  filter(Sample_ID %in% subjects_to_use)

################################################################################
# Create 10 repeats of 10 folds to use for all analyses
################################################################################

# Load contrasts to control participants
study_comparisons_to_control <- readRDS(paste0(rdata_path, dataset_ID, "_comparisons_to_control.Rds"))

# Use 10-fold CV
k <- 10

# Iterate over each contrast to controls per dataset
if (!file.exists(paste0(rdata_path, dataset_ID, "_samples_per_10_folds_10_repeats.Rds"))) {
  sample_folds_list <- list()
  for (i in 1:nrow(study_comparisons_to_control)) {
    study <- study_comparisons_to_control$Study[i]
    group_to_compare <- study_comparisons_to_control$Group_to_Compare[i]
    
    group_folds <- list()
    
    control_data <- sample_metadata %>%
      filter(Study == study, 
             Diagnosis == "Control") %>%
      distinct(Sample_ID, Diagnosis) 
    
    control_vector <- control_data %>%
      pull(Diagnosis)
    
    group_data <- sample_metadata %>%
      filter(Study == study,
             Diagnosis == group_to_compare) %>%
      distinct(Sample_ID, Diagnosis) 
    
    group_vector <- group_data %>%
      pull(Diagnosis)
    
    merged_data <- plyr::rbind.fill(control_data,
                                    group_data)
    
    # 10 repeats
    for (r in 1:10) {
      # 10 folds per repeat for 10-fold CV
      group_folds_r <- caret::createFolds(c(control_vector, group_vector),
                                          k = k,
                                          list = TRUE,
                                          returnTrain = FALSE)
      
      # Map indices to sample IDs to avoid ambiguity
      group_folds_samples_r <- purrr::map(group_folds_r, ~ (merged_data[.x,]) %>% pull(Sample_ID))
      
      # Append to list
      group_folds[[r]] <- group_folds_samples_r
    }
    
    # Combine group splits
    group_df <- as.data.frame(cbind(Study = rep(study, 10),
                                    Group_to_Compare = rep(group_to_compare, 10),
                                    Repeat_Number = 1:10,
                                    Folds = group_folds))
    
    # Append to dataframe
    sample_folds_list <- list.append(sample_folds_list, group_df)
  }
  
  # Combine fold results into a dataframe
  sample_folds <- do.call(plyr::rbind.fill, sample_folds_list)
  
  # Save to an Rds file
  saveRDS(sample_folds, paste0(rdata_path, dataset_ID, 
                               sprintf("_samples_per_%s_folds_10_repeats.Rds",
                                       k)))
} else {
  sample_folds <- readRDS(paste0(rdata_path, dataset_ID, 
                                 sprintf("_samples_per_%s_folds_10_repeats.Rds",
                                         k)))
}

################################################################################
# Define grouping_parameters
################################################################################

grouping_param_df <- data.frame(grouping_type = c("ROI", "Feature", "Combo"),
                                grouping_var = c("Brain_Region", "Feature", "Combo"),
                                SVM_feature_var = c("Feature", "Brain_Region", "Combo")) 

# Use a linear kernel
svm_kernel = "linear"
weighting_name <- "inv_prob"
use_inv_prob_weighting <- TRUE

if (add_catch2) {
  output_file_RDS <- paste0(rdata_path, 
                            sprintf("%s_univariate_CV_linear_SVM_%s_and_catch2_%s.Rds",
                                    dataset_ID,
                                    univariate_feature_set, 
                                    weighting_name))
} else {
  output_file_RDS <- paste0(rdata_path, 
                            sprintf("%s_univariate_CV_linear_SVM_%s_%s.Rds",
                                    dataset_ID,
                                    univariate_feature_set, 
                                    weighting_name))
}

if (!file.exists(output_file_RDS)) {
  # Load non-control data 
  univariate_features_groups <- readRDS(paste0(rdata_path, "UCLA_CNP_ABIDE_ASD_", univariate_feature_set,
                                               "_and_catch2_filtered_zscored.Rds")) %>%
    left_join(., sample_metadata) %>%
    filter(Diagnosis != "Control")
  
  # Load control data
  univariate_features_controls <- readRDS(paste0(rdata_path, "UCLA_CNP_ABIDE_ASD_", univariate_feature_set,
                                                 "_and_catch2_filtered_zscored.Rds")) %>%
    left_join(., sample_metadata) %>%
    filter(Diagnosis == "Control")
  
  univariate_class_res_list <- list()
  
  for (i in 1:nrow(study_comparisons_to_control)) {
    study <- study_comparisons_to_control$Study[i]
    group_to_compare <- study_comparisons_to_control$Group_to_Compare[i] 
    
    # Get subjects to compare
    sample_groups <- sample_metadata %>%
      filter(Sample_ID %in% subjects_to_use,
             Study == study,
             Diagnosis %in% c("Control", group_to_compare)) %>%
      distinct(Sample_ID, Diagnosis)
    
    # Isolate group folds
    group_folds <- sample_folds %>%
      filter(Study == study,
             Group_to_Compare == group_to_compare) %>%
      pull(Folds)
    
    # Combine control data with diagnosis group data
    group_data_for_SVM <- univariate_features_groups %>%
      filter(Study == study,
             Diagnosis == group_to_compare) %>%
      plyr::rbind.fill(., univariate_features_controls %>%
                         filter(Study == study)) 
    
    group_data_for_SVM_feature_set <- group_data_for_SVM %>%
      filter(feature_set == univariate_feature_set)
    
    # If add catch2 flag is on, also prepare SVM data for catch2
    if (add_catch2) {
      group_data_for_SVM_catch2 <- subset(group_data_for_SVM, feature_set=="catch2")
    }
    
    # Run SVM
    for (j in 1:nrow(grouping_param_df)) {
      grouping_type = grouping_param_df$grouping_type[j]
      grouping_var = grouping_param_df$grouping_var[j]
      SVM_feature_var = grouping_param_df$SVM_feature_var[j]
      
      group_wise_SVM_CV_feature_set <- 1:length(group_folds) %>%
        purrr::map_df( ~ run_univariate_cv_svm_by_input_var(feature_matrix = group_data_for_SVM_feature_set,
                                                            dataset_ID = dataset_ID,
                                                            univariate_feature_set = univariate_feature_set,
                                                            sample_groups = sample_groups,
                                                            svm_kernel = svm_kernel,
                                                            grouping_var = grouping_var,
                                                            flds = group_folds[[.x]],
                                                            repeat_number = .x,
                                                            svm_feature_var = SVM_feature_var,
                                                            out_of_sample_only = TRUE,
                                                            use_inv_prob_weighting = use_inv_prob_weighting)) %>%
        mutate(Study = study,
               Analysis_Type = grouping_type,
               Group_to_Compare = group_to_compare)
      
      univariate_class_res_list <- list.append(univariate_class_res_list, group_wise_SVM_CV_feature_set)
      
      # Check if add_catch2 flag is set
      if (add_catch2) {
        group_wise_SVM_CV_catch2 <- 1:length(group_folds) %>%
          purrr::map_df( ~ run_univariate_cv_svm_by_input_var(feature_matrix = group_data_for_SVM_catch2,
                                                              dataset_ID = dataset_ID,
                                                              univariate_feature_set = "catch2",
                                                              sample_groups = sample_groups,
                                                              svm_kernel = svm_kernel,
                                                              grouping_var = grouping_var,
                                                              flds = group_folds[[.x]],
                                                              repeat_number = .x,
                                                              svm_feature_var = SVM_feature_var,
                                                              out_of_sample_only = TRUE,
                                                              use_inv_prob_weighting = use_inv_prob_weighting)) %>%
          mutate(Study = study,
                 Analysis_Type = grouping_type,
                 Group_to_Compare = group_to_compare,
                 univariate_feature_set = "catch2")
        
        
        univariate_class_res_list <- list.append(univariate_class_res_list, group_wise_SVM_CV_catch2)
      }
    }
  }
  univariate_class_res <- do.call(plyr::rbind.fill, univariate_class_res_list)
  saveRDS(univariate_class_res, file = output_file_RDS)
} else {
  univariate_class_res <- readRDS(file = output_file_RDS)
}


################################################################################
# Calculate balanced accuracy across folds/repeats

# Function to calculate the repeat-wise balanced accuracy for the  given group and study
find_balanced_accuracy_by_repeat <- function(group, study, univariate_class_res) {
  repeat_wise_balanced_accuracy <- univariate_class_res %>%
    # Filter by study and group to compare against controls
    filter(Study == study, Group_to_Compare == group) %>%
    
    # Drop unused factor levels
    droplevels() %>%
    group_by(grouping_var, univariate_feature_set, num_SVM_features, 
             Study, Analysis_Type, Group_to_Compare, repeat_number) %>%
    
    # calculate accuracy and balanced accuracy
    summarise(accuracy = sum(Prediction_Correct) / n(),
              balanced_accuracy = caret::confusionMatrix(data = Predicted_Diagnosis,
                                                         reference = Actual_Diagnosis)$byClass[["Balanced Accuracy"]])
  
  return(repeat_wise_balanced_accuracy)
}

output_file_name <- ifelse(add_catch2,
                           paste0(rdata_path, sprintf("%s_univariate_CV_linear_SVM_%s_and_catch2_balanced_accuracy_by_repeats.Rds",
                                                      dataset_ID,
                                                      univariate_feature_set)),
                           paste0(rdata_path, sprintf("%s_univariate_CV_linear_SVM_%s_balanced_accuracy_by_repeats.Rds",
                                                      dataset_ID,
                                                      univariate_feature_set)))

if (!file.exists(output_file_name)) { 
  # Find the balanced accuracy for each study and comparison group combo
  balanced_accuracy_by_repeat <- 1:nrow(study_comparisons_to_control) %>%
    purrr::map_df(~ find_balanced_accuracy_by_repeat(group = study_comparisons_to_control$Group_to_Compare[.x],
                                                     study = study_comparisons_to_control$Study[.x],
                                                     univariate_class_res = univariate_class_res))
  
  # Save the repeat-wise balanced accuracy
  saveRDS(balanced_accuracy_by_repeat, 
          file = output_file_name)
} else {
  balanced_accuracy_by_repeat <- readRDS(file = output_file_name)
}

################################################################################
# Calculate balanced accuracy across folds/repeats

find_balanced_accuracy_by_fold <- function(group, study, univariate_class_res) {
  full_balanced_accuracy <- univariate_class_res %>%
    # Filter by study and group to compare against controls
    filter(Study == study, Group_to_Compare == group) %>%
    
    # Drop unused factor levels
    droplevels() %>%
    group_by(grouping_var, univariate_feature_set, num_SVM_features, 
             Study, Analysis_Type, Group_to_Compare, fold_number, repeat_number) %>%
    
    # calculate accuracy and balanced accuracy per fold/repeat
    summarise(mean_accuracy = sum(Prediction_Correct) / n(),
              mean_balanced_accuracy = caret::confusionMatrix(data = Predicted_Diagnosis,
                                                              reference = Actual_Diagnosis)$byClass[["Balanced Accuracy"]]) %>%
    
    # take the average and SD across all folds/repeats
    ungroup() %>%
    group_by(grouping_var, univariate_feature_set, num_SVM_features, 
             Study, Analysis_Type, Group_to_Compare) %>%
    summarise(balanced_accuracy = mean(mean_balanced_accuracy, na.rm=T),
              balanced_accuracy_SD = sd(mean_balanced_accuracy, na.rm=T),
              accuracy = mean(mean_accuracy, na.rm=T),
              accuracy_SD = sd(mean_accuracy, na.rm=T))
  
  
  return(full_balanced_accuracy)
}


# Save balanced accuracy results aggregated across repeats
output_file_name <- ifelse(add_catch2,
                           paste0(rdata_path, sprintf("%s_univariate_CV_linear_SVM_%s_and_catch2_balanced_accuracy.Rds",
                                                      dataset_ID,
                                                      univariate_feature_set)),
                           paste0(rdata_path, sprintf("%s_univariate_CV_linear_SVM_%s_balanced_accuracy.Rds",
                                                      dataset_ID,
                                                      univariate_feature_set)))

if (!file.exists(output_file_name)) {
  balanced_accuracy_results <- 1:nrow(study_comparisons_to_control) %>%
    purrr::map_df(~ find_balanced_accuracy_by_fold(group = study_comparisons_to_control$Group_to_Compare[.x],
                                                   study = study_comparisons_to_control$Study[.x],
                                                   univariate_class_res = univariate_class_res))
  
  saveRDS(balanced_accuracy_results, 
          file = output_file_name)
} else {
  balanced_accuracy_results <- readRDS(file = output_file_name)
}

