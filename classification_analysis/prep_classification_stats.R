#-------------------------------------------------------------------------------
# Define paths
#-------------------------------------------------------------------------------

python_to_use <- "~/.conda/envs/pyspi/bin/python3"
# python_to_use <- "/Users/abry4213/opt/anaconda3/envs/pyspi/bin/python3"
pairwise_feature_set <- "pyspi14"
github_dir <- "~/github/"
data_path <- "~/data/"
scaler <- "robustsigmoid"
UCLA_CNP_sample_metadata_file <- "UCLA_CNP_sample_metadata.feather"
UCLA_CNP_noise_proc <- "AROMA+2P+GMR"
ABIDE_ASD_sample_metadata_file <- "ABIDE_ASD_sample_metadata.feather"
ABIDE_ASD_noise_proc <- "FC1000"
output_data_path <- glue::glue("{data_path}/TS_feature_manuscript/")
study_group_df <- data.frame(Study = c(rep("UCLA_CNP", 3), "ABIDE_ASD"),
                             Noise_Proc = c(rep("AROMA+2P+GMR",3), "FC1000"),
                             Comparison_Group = c("Schizophrenia", "ADHD", "Bipolar", "ASD"))
univariate_feature_sets <- c("catch22")

# DIY rlist::list.append
list.append <- function (.data, ...) 
{
  if (is.list(.data)) {
    c(.data, list(...))
  }
  else {
    c(.data, ..., recursive = FALSE)
  }
}

reticulate::use_python(python_to_use)
library(feather)
library(tidyverse)
library(reticulate)
library(glue)

# Import pyarrow.feather as pyarrow_feather
pyarrow_feather <- import("pyarrow.feather")

################################################################################
# Function definitions
################################################################################

# Functions to calculate empirical p-value
compare_main_and_null <- function(main_df_iter, null_distribution_df) {
  # Filter null to this data -- keep all grouping vars in this analysis type
  null_distribution_df <- null_distribution_df %>%
    dplyr::select(-index, -group_var, -Null_Iter_Number) %>%
    semi_join(., main_df_iter)
  
  # Compare main balanced accuracy with that of the empirical null distribution
  main_balanced_accuracy_value <- main_df_iter$Balanced_Accuracy_Across_Repeats
  null_balanced_accuracy_values <- null_distribution_df$Null_Balanced_Accuracy
  
  # Find proportion of iterations for which the main balanced accuracy is greater
  prop_main_greater <- sum(main_balanced_accuracy_value > null_balanced_accuracy_values)/length(null_balanced_accuracy_values)
  
  # Find p-value
  p_value <- 1 - prop_main_greater
  
  # Organise into dataframe
  main_df_iter$p_value <- p_value
  return(main_df_iter)
}

################################################################################
# Compile balanced accuracy results
################################################################################

#### Univariate ####
# Load balanced accuracy data, or construct if needed
if (!file.exists(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_univariate_{scaler}_scaler_balanced_accuracy_all_folds.feather"))) {
  univariate_balanced_accuracy_all_folds_list <- list()
  # First iterate over each study/comparison group
  for (i in 1:nrow(study_group_df)) {
    dataset_ID <- study_group_df$Study[i]
    noise_proc <- study_group_df$Noise_Proc[i]
    noise_label = gsub("\\+", "_", noise_proc)
    comparison_group <- study_group_df$Comparison_Group[i]
    
    # Now iterate over each univariate feature set
    for (featset in univariate_feature_sets) {
      balacc_across_folds <- pyarrow_feather$read_feather(glue("{data_path}/{dataset_ID}/processed_data/{dataset_ID}_{comparison_group}_Univariate_{featset}_{scaler}_scaler_SVM_balanced_accuracy.feather"))
      balacc_across_folds$Study <- dataset_ID
      balacc_across_folds$Univariate_Feature_Set <- featset
      # Append to list
      univariate_balanced_accuracy_all_folds_list <- list.append(univariate_balanced_accuracy_all_folds_list, balacc_across_folds)
    }
  }
  
  # Combine the list results into a dataframe
  univariate_balanced_accuracy_all_folds <- do.call(plyr::rbind.fill, 
                                                    univariate_balanced_accuracy_all_folds_list) 
  
  # Save to feather file
  feather::write_feather(univariate_balanced_accuracy_all_folds,
                         glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_univariate_{scaler}_scaler_balanced_accuracy_all_folds.feather"))
} else {
  univariate_balanced_accuracy_all_folds <- feather::read_feather(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_univariate_{scaler}_scaler_balanced_accuracy_all_folds.feather"))
}

# Aggregate the main results across folds and then across repeats
univariate_balanced_accuracy <- univariate_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Analysis_Type, group_var, Repeat_Number) %>%
  reframe(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
            Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T)) %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Analysis_Type, group_var) %>%
  reframe(Balanced_Accuracy_Across_Repeats = mean(Balanced_Accuracy_Across_Folds, na.rm=T),
          Balanced_Accuracy_Across_Repeats_SD = sd(Balanced_Accuracy_Across_Folds, na.rm=T))

#### Pairwise ####            
# Load balanced accuracy data, or construct if needed
if (!file.exists(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_pairwise_{scaler}_scaler_balanced_accuracy_all_folds.feather"))) {
  pairwise_balanced_accuracy_all_folds_list <- list()
  # First iterate over each study/comparison group
  for (i in 1:nrow(study_group_df)) {
    dataset_ID <- study_group_df$Study[i]
    noise_proc <- study_group_df$Noise_Proc[i]
    noise_label = gsub("\\+", "_", noise_proc)
    comparison_group <- study_group_df$Comparison_Group[i]
    
    balacc_across_folds <- pyarrow_feather$read_feather(glue("{data_path}/{dataset_ID}/processed_data/{dataset_ID}_{comparison_group}_Pairwise_{pairwise_feature_set}_{scaler}_scaler_SVM_balanced_accuracy.feather"))
    balacc_across_folds$Study <- dataset_ID
    balacc_across_folds$Pairwise_Feature_Set <- pairwise_feature_set
    # Append to list
    pairwise_balanced_accuracy_all_folds_list <- list.append(pairwise_balanced_accuracy_all_folds_list, balacc_across_folds)
  }
  
  # Combine the list results into a dataframe
  pairwise_balanced_accuracy_all_folds <- do.call(plyr::rbind.fill, 
                                                    pairwise_balanced_accuracy_all_folds_list) %>%
    mutate(Pairwise_Feature_Set = Pairwise_Feature_Set)
  
  # Save to feather file
  feather::write_feather(pairwise_balanced_accuracy_all_folds,
                         glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_pairwise_{scaler}_scaler_balanced_accuracy_all_folds.feather"))
} else {
  pairwise_balanced_accuracy_all_folds <- feather::read_feather(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_pairwise_{scaler}_scaler_balanced_accuracy_all_folds.feather"))
}

# Aggregate the main results across folds and then across repeats
pairwise_balanced_accuracy <- pairwise_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Pairwise_Feature_Set, Analysis_Type, group_var, Repeat_Number) %>%
  reframe(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
          Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T)) %>%
  group_by(Study, Comparison_Group, Pairwise_Feature_Set, Analysis_Type, group_var) %>%
  reframe(Balanced_Accuracy_Across_Repeats = mean(Balanced_Accuracy_Across_Folds, na.rm=T),
          Balanced_Accuracy_Across_Repeats_SD = sd(Balanced_Accuracy_Across_Folds, na.rm=T))

################################################################################
# Compile individual subject predictions
################################################################################

# Univariate
if (!file.exists(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_univariate_{scaler}_scaler_subject_class_predictions.feather"))) {
  univariate_subject_predictions_list <- list()
  # First iterate over each study/comparison group
  for (i in 1:nrow(study_group_df)) {
    dataset_ID <- study_group_df$Study[i]
    noise_proc <- study_group_df$Noise_Proc[i]
    noise_label = gsub("\\+", "_", noise_proc)
    comparison_group <- study_group_df$Comparison_Group[i]
    
    # Now iterate over each univariate feature set
    for (featset in univariate_feature_sets) {
      subject_predictions <- pyarrow_feather$read_feather(glue("{data_path}/{dataset_ID}/processed_data/{dataset_ID}_{comparison_group}_Univariate_{featset}_{scaler}_scaler_SVM_sample_predictions.feather"))
      subject_predictions$Study <- dataset_ID
      subject_predictions$Univariate_Feature_Set <- featset
      # Append to list
      univariate_subject_predictions_list <- list.append(univariate_subject_predictions_list, subject_predictions)
    }
  }
  
  # Combine the list results into a dataframe
  univariate_subject_predictions <- do.call(plyr::rbind.fill, 
                                                    univariate_subject_predictions_list) %>%
    mutate(Analysis_Type = case_when(str_detect(group_var, "_") ~ "TS_Feature",
                                     group_var == "Combo" ~ "Combo",
                                     T ~ "Brain_Region"))
  
  # Save to feather file
  feather::write_feather(univariate_subject_predictions,
                         glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_univariate_{scaler}_scaler_subject_class_predictions.feather"))
} else {
  univariate_subject_predictions <- feather::read_feather(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_univariate_{scaler}_scaler_subject_class_predictions.feather"))
}

# Pairwise
if (!file.exists(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_pairwise_{scaler}_scaler_subject_class_predictions.feather"))) {
  pairwise_subject_predictions_list <- list()
  # First iterate over each study/comparison group
  for (i in 1:nrow(study_group_df)) {
    dataset_ID <- study_group_df$Study[i]
    noise_proc <- study_group_df$Noise_Proc[i]
    noise_label = gsub("\\+", "_", noise_proc)
    comparison_group <- study_group_df$Comparison_Group[i]
    
    subject_predictions <- pyarrow_feather$read_feather(glue("{data_path}/{dataset_ID}/processed_data/{dataset_ID}_{comparison_group}_Pairwise_{pairwise_feature_set}_{scaler}_scaler_SVM_sample_predictions.feather"))
    subject_predictions$Study <- dataset_ID
    subject_predictions$Pairwise_Feature_Set <- pairwise_feature_set
    # Append to list
    pairwise_subject_predictions_list <- list.append(pairwise_subject_predictions_list, subject_predictions)
  }
  
  # Combine the list results into a dataframe
  pairwise_subject_predictions <- do.call(plyr::rbind.fill, 
                                                    pairwise_subject_predictions_list) %>%
    mutate(Analysis_Type = case_when(str_detect(group_var, "_") ~ "TS_Feature",
                                     group_var == "Combo" ~ "Combo",
                                     T ~ "Brain_Region"))
  
  # Save to feather file
  feather::write_feather(pairwise_subject_predictions,
                         glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_pairwise_{scaler}_scaler_subject_class_predictions.feather"))
} else {
  pairwise_subject_predictions <- feather::read_feather(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_pairwise_{scaler}_scaler_subject_class_predictions.feather"))
}

################################################################################
# Compile fold assignments
################################################################################

# Univariate
if (!file.exists(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_univariate_{scaler}_scaler_fold_assignments.feather"))) {
  univariate_fold_assignments_list <- list()
  # First iterate over each study/comparison group
  for (i in 1:nrow(study_group_df)) {
    dataset_ID <- study_group_df$Study[i]
    noise_proc <- study_group_df$Noise_Proc[i]
    noise_label = gsub("\\+", "_", noise_proc)
    comparison_group <- study_group_df$Comparison_Group[i]
    
    # Now iterate over each univariate feature set
    for (featset in univariate_feature_sets) {
      fold_assignments <- pyarrow_feather$read_feather(glue("{data_path}/{dataset_ID}/processed_data/{dataset_ID}_{comparison_group}_Univariate_{featset}_{scaler}_scaler_SVM_fold_assignments.feather"))
      fold_assignments$Study <- dataset_ID
      fold_assignments$Univariate_Feature_Set <- featset
      # Append to list
      univariate_fold_assignments_list <- list.append(univariate_fold_assignments_list, fold_assignments)
    }
  }
  
  # Combine the list results into a dataframe
  univariate_fold_assignments <- do.call(plyr::rbind.fill, 
                                                    univariate_fold_assignments_list) %>%
    mutate(Analysis_Type = case_when(str_detect(group_var, "_") ~ "TS_Feature",
                                     group_var == "Combo" ~ "Combo",
                                     T ~ "Brain_Region"))
  
  # Save to feather file
  feather::write_feather(univariate_fold_assignments,
                         glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_univariate_{scaler}_scaler_fold_assignments.feather"))
} else {
  univariate_fold_assignments <- feather::read_feather(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_univariate_{scaler}_scaler_fold_assignments.feather"))
}

# Pairwise
if (!file.exists(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_pairwise_{scaler}_scaler_fold_assignments.feather"))) {
  pairwise_fold_assignments_list <- list()
  # First iterate over each study/comparison group
  for (i in 1:nrow(study_group_df)) {
    dataset_ID <- study_group_df$Study[i]
    noise_proc <- study_group_df$Noise_Proc[i]
    noise_label = gsub("\\+", "_", noise_proc)
    comparison_group <- study_group_df$Comparison_Group[i]
    
    fold_assignments <- pyarrow_feather$read_feather(glue("{data_path}/{dataset_ID}/processed_data/{dataset_ID}_{comparison_group}_Pairwise_{pairwise_feature_set}_{scaler}_scaler_SVM_fold_assignments.feather"))
    fold_assignments$Study <- dataset_ID
    fold_assignments$Pairwise_Feature_Set <- pairwise_feature_set
    # Append to list
    pairwise_fold_assignments_list <- list.append(pairwise_fold_assignments_list, fold_assignments)
  }
  
  # Combine the list results into a dataframe
  pairwise_fold_assignments <- do.call(plyr::rbind.fill, 
                                                    pairwise_fold_assignments_list) %>%
    mutate(Analysis_Type = case_when(str_detect(group_var, "_") ~ "TS_Feature",
                                     group_var == "Combo" ~ "Combo",
                                     T ~ "Brain_Region"))
  
  # Save to feather file
  feather::write_feather(pairwise_fold_assignments,
                         glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_pairwise_{scaler}_scaler_fold_assignments.feather"))
} else {
  pairwise_fold_assignments <- feather::read_feather(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_pairwise_{scaler}_scaler_fold_assignments.feather"))
}


################################################################################
# Compile SVM coefficient results
################################################################################

# Univariate
if (!file.exists(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_univariate_{scaler}_scaler_SVM_coefficients.feather"))) {
  univariate_SVM_coefs_list <- list()
  # First iterate over each study/comparison group
  for (i in 1:nrow(study_group_df)) {
    dataset_ID <- study_group_df$Study[i]
    noise_proc <- study_group_df$Noise_Proc[i]
    noise_label = gsub("\\+", "_", noise_proc)
    comparison_group <- study_group_df$Comparison_Group[i]
    
    # Now iterate over each univariate feature set
    for (featset in univariate_feature_sets) {
      SVM_coefs <- pyarrow_feather$read_feather(glue("{data_path}/{dataset_ID}/processed_data/{dataset_ID}_{comparison_group}_Univariate_{featset}_{scaler}_scaler_SVM_fold_SVM_coefficients.feather")) %>%
        group_by(`Feature Name`, Analysis_Type, group_var, Comparison_Group, Scaling_Type) %>%
        summarise(CoefficientM = mean(Coefficient, na.rm=T),
                  Coefficient_SD = sd(Coefficient, na.rm=T)) %>%
        dplyr::rename("Coefficient" = "CoefficientM",
                      "Feature_Name" = "Feature Name")
      SVM_coefs$Study <- dataset_ID
      SVM_coefs$Univariate_Feature_Set <- featset
      # Append to list
      univariate_SVM_coefs_list <- list.append(univariate_SVM_coefs_list, SVM_coefs)
    }
  }
  
  # Combine the list results into a dataframe
  univariate_SVM_coefs <- do.call(plyr::rbind.fill, 
                                              univariate_SVM_coefs_list) %>%
    mutate(Analysis_Type = case_when(str_detect(group_var, "_") ~ "TS_Feature",
                                     group_var == "Combo" ~ "Combo",
                                     T ~ "Brain_Region"))
  
  # Save to feather file
  feather::write_feather(univariate_SVM_coefs,
                         glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_univariate_{scaler}_scaler_SVM_coefficients.feather"))
} else {
  univariate_SVM_coefs <- feather::read_feather(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_univariate_{scaler}_scaler_SVM_coefficients.feather"))
}

# Pairwise
if (!file.exists(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_pairwise_{scaler}_scaler_SVM_coefficients.feather"))) {
  pairwise_SVM_coefs_list <- list()
  # First iterate over each study/comparison group
  for (i in 1:nrow(study_group_df)) {
    dataset_ID <- study_group_df$Study[i]
    noise_proc <- study_group_df$Noise_Proc[i]
    noise_label = gsub("\\+", "_", noise_proc)
    comparison_group <- study_group_df$Comparison_Group[i]
    
    SVM_coefs <- pyarrow_feather$read_feather(glue("{data_path}/{dataset_ID}/processed_data/{dataset_ID}_{comparison_group}_Pairwise_{pairwise_feature_set}_{scaler}_scaler_SVM_fold_SVM_coefficients.feather")) %>%
      group_by(`Feature Name`, Analysis_Type, group_var, Comparison_Group, Scaling_Type) %>%
      summarise(CoefficientM = mean(Coefficient, na.rm=T),
                Coefficient_SD = sd(Coefficient, na.rm=T)) %>%
      dplyr::rename("Coefficient" = "CoefficientM",
                    "Feature_Name" = "Feature Name")
    SVM_coefs$Study <- dataset_ID
    SVM_coefs$Pairwise_Feature_Set <- pairwise_feature_set
    # Append to list
    pairwise_SVM_coefs_list <- list.append(pairwise_SVM_coefs_list, SVM_coefs)
  }
  
  # Combine the list results into a dataframe
  pairwise_SVM_coefs <- do.call(plyr::rbind.fill, 
                                              pairwise_SVM_coefs_list) %>%
    mutate(Analysis_Type = case_when(str_detect(group_var, "_") ~ "TS_Feature",
                                     group_var == "Combo" ~ "Combo",
                                     T ~ "Brain_Region"))
  
  # Save to feather file
  feather::write_feather(pairwise_SVM_coefs,
                         glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_pairwise_{scaler}_scaler_SVM_coefficients.feather"))
} else {
  pairwise_SVM_coefs <- feather::read_feather(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_pairwise_{scaler}_scaler_SVM_coefficients.feather"))
}

################################################################################
# Compile null results
################################################################################

# Univariate
if (!file.exists(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_univariate_{scaler}_scaler_null_balanced_accuracy_distributions.feather"))) {
  univariate_null_balanced_accuracy_list <- list()
  # First iterate over each study/comparison group
  for (i in 1:nrow(study_group_df)) {
    dataset_ID <- study_group_df$Study[i]
    noise_proc <- study_group_df$Noise_Proc[i]
    noise_label = gsub("\\+", "_", noise_proc)
    comparison_group <- study_group_df$Comparison_Group[i]
    
    # Now iterate over each univariate feature set
    for (featset in univariate_feature_sets) {
      null_dist <- pyarrow_feather$read_feather(glue("{data_path}/{dataset_ID}/processed_data/{dataset_ID}_{comparison_group}_Univariate_{featset}_{scaler}_scaler_SVM_null_balanced_accuracy_distributions.feather"))
      null_dist$Study <- dataset_ID
      null_dist$Univariate_Feature_Set <- featset
      # Append to list
      univariate_null_balanced_accuracy_list <- list.append(univariate_null_balanced_accuracy_list, null_dist)
    }
  }
  
  # Combine the list results into a dataframe
  univariate_null_balanced_accuracy <- do.call(plyr::rbind.fill, 
                                               univariate_null_balanced_accuracy_list) %>%
    mutate(Analysis_Type = case_when(str_detect(group_var, "_") ~ "Univariate_TS_Feature",
                                     group_var == "Combo" ~ "Univariate_Combo",
                                     T ~ "Univariate_Brain_Region"))
  
  # Save to feather file
  feather::write_feather(univariate_null_balanced_accuracy,
                         glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_univariate_{scaler}_scaler_null_balanced_accuracy_distributions.feather"))
} else {
  univariate_null_balanced_accuracy <- feather::read_feather(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_univariate_{scaler}_scaler_null_balanced_accuracy_distributions.feather"))
}

# Pairwise
if (!file.exists(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_pairwise_{scaler}_scaler_null_balanced_accuracy_distributions.feather"))) {
  pairwise_null_balanced_accuracy_list <- list()
  # First iterate over each study/comparison group
  for (i in 1:nrow(study_group_df)) {
    dataset_ID <- study_group_df$Study[i]
    noise_proc <- study_group_df$Noise_Proc[i]
    noise_label = gsub("\\+", "_", noise_proc)
    comparison_group <- study_group_df$Comparison_Group[i]
    
    null_dist <- pyarrow_feather$read_feather(glue("{data_path}/{dataset_ID}/processed_data/{dataset_ID}_{comparison_group}_Pairwise_{pairwise_feature_set}_{scaler}_scaler_SVM_null_balanced_accuracy_distributions.feather"))
    null_dist$Study <- dataset_ID
    null_dist$Pairwise_Feature_Set <- pairwise_feature_set
    # Append to list
    pairwise_null_balanced_accuracy_list <- list.append(pairwise_null_balanced_accuracy_list, null_dist)
  }
  
  # Combine the list results into a dataframe
  pairwise_null_balanced_accuracy <- do.call(plyr::rbind.fill, 
                                               pairwise_null_balanced_accuracy_list) %>%
    mutate(Analysis_Type = "Pairwise_SPI")
  
  # Save to feather file
  feather::write_feather(pairwise_null_balanced_accuracy,
                         glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_pairwise_{scaler}_scaler_null_balanced_accuracy_distributions.feather"))
} else {
  pairwise_null_balanced_accuracy <- feather::read_feather(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_pairwise_{scaler}_scaler_null_balanced_accuracy_distributions.feather"))
}

################################################################################
# Compile p-value results
################################################################################

# Univariate
if (!file.exists(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_univariate_{scaler}_scaler_empirical_p_values.feather"))) {
  univariate_split <- univariate_balanced_accuracy %>%
    group_by(Study, Comparison_Group, Analysis_Type, Univariate_Feature_Set, group_var) %>%
    group_split()
  
  univariate_p_values <- univariate_split %>%
    purrr::map_df(~ compare_main_and_null(main_df_iter = .x,
                                          null_distribution_df = univariate_null_balanced_accuracy))
  
  # Adjust p-values by group
  univariate_p_values <- univariate_p_values %>%
    group_by(Study, Comparison_Group, Univariate_Feature_Set, Analysis_Type) %>%
    mutate(p_value_BH = p.adjust(p_value, method="BH"),
           p_value_Bonferroni = p.adjust(p_value, method="bonferroni"))
  
  feather::write_feather(univariate_p_values, glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_univariate_{scaler}_scaler_empirical_p_values.feather"))
} else {
  univariate_p_values <- feather::read_feather(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_univariate_{scaler}_scaler_empirical_p_values.feather"))
}

# Pairwise
if (!file.exists(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_pairwise_{scaler}_scaler_empirical_p_values.feather"))) {
  pairwise_split <- pairwise_balanced_accuracy %>%
    group_by(Study, Comparison_Group, Analysis_Type, Pairwise_Feature_Set, group_var) %>%
    group_split()
  
  pairwise_p_values <- pairwise_split %>%
    purrr::map_df(~ compare_main_and_null(main_df_iter = .x,
                                          null_distribution_df = pairwise_null_balanced_accuracy))
  
  # Adjust p-values by group
  pairwise_p_values <- pairwise_p_values %>%
    group_by(Study, Comparison_Group, Pairwise_Feature_Set, Analysis_Type) %>%
    mutate(p_value_BH = p.adjust(p_value, method="BH"), 
           p_value_Bonferroni = p.adjust(p_value, method="bonferroni"))
  
  feather::write_feather(pairwise_p_values, glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_pairwise_{scaler}_scaler_empirical_p_values.feather"))
} else {
  pairwise_p_values <- feather::read_feather(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_pairwise_{scaler}_scaler_empirical_p_values.feather"))
}