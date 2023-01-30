#-------------------------------------------------------------------------------
# Define paths
#-------------------------------------------------------------------------------

# python_to_use <- "~/.conda/envs/pyspi/bin/python3"
python_to_use <- "/Users/abry4213/opt/anaconda3/envs/pyspi/bin/python3"
pairwise_feature_set <- "pyspi14"
github_dir <- "~/github/"
data_path <- "~/data/"
UCLA_CNP_sample_metadata_file <- "UCLA_CNP_sample_metadata.feather"
UCLA_CNP_noise_proc <- "AROMA+2P+GMR"
ABIDE_ASD_sample_metadata_file <- "ABIDE_ASD_sample_metadata.feather"
ABIDE_ASD_noise_proc <- "FC1000"
output_data_path <- glue::glue("{data_path}/TS_feature_manuscript/")
study_group_df <- data.frame(Study = c(rep("UCLA_CNP", 3), "ABIDE_ASD"),
                             Noise_Proc = c(rep("AROMA+2P+GMR",3), "FC1000"),
                             Comparison_Group = c("Schizophrenia", "ADHD", "Bipolar", "ASD"))
univariate_feature_sets <- c("catch22", "catch2", "catch24")

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
    dplyr::select(-index, -group_var) %>%
    semi_join(., main_df_iter)
  
  # Compare main balanced accuracy with that of the empirical null distribution
  main_balanced_accuracy_value <- main_df_iter$Balanced_Accuracy_Across_Folds
  null_balanced_accuracy_values <- null_distribution_df$Null_Balanced_Accuracy
  
  # Find proportion of iterations for which the main balanced accuracy is greater
  prop_main_greater <- sum(main_balanced_accuracy_value > null_balanced_accuracy_values)/length(null_balanced_accuracy_values)
  
  # Find p-value
  p_value <- 1 - prop_main_greater
  
  # Organise into dataframe
  main_df_iter$p_value <- p_value
  return(main_df_iter)
}

calculate_empirical_p_values <- function(main_balanced_accuracy_split, 
                                         null_balanced_accuracy) {
  
  # Iterate over splits
  p_values <- main_balanced_accuracy_split %>%
    purrr::map_df(~ compare_main_and_null(main_df_iter = .x,
                                          null_distribution_df = null_balanced_accuracy))
  
  return(p_values)
}

################################################################################
# Sanity checks
################################################################################

univariate_fold_assignments_list <- list()

# First iterate over each study/comparison group
for (i in 1:nrow(study_group_df)) {
  dataset_ID <- study_group_df$Study[i]
  noise_proc <- study_group_df$Noise_Proc[i]
  noise_label = gsub("\\+", "_", noise_proc)
  comparison_group <- study_group_df$Comparison_Group[i]
  
  # Now iterate over each univariate feature set
  for (featset in univariate_feature_sets) {
    fold_assignments_iter <- pyarrow_feather$read_feather(glue("{data_path}/{dataset_ID}/processed_data/{dataset_ID}_{comparison_group}_Univariate_{featset}_SVM_fold_assignments.feather"))
    fold_assignments_iter$Study <- dataset_ID
    fold_assignments_iter$Univariate_Feature_Set <- featset
    # Append to list
    univariate_fold_assignments_list <- list.append(univariate_fold_assignments_list, fold_assignments_iter)
  }
}

# Combine the list results into a dataframe
all_univariate_fold_assignments <- do.call(plyr::rbind.fill, 
                                           univariate_fold_assignments_list) %>%
  mutate(Analysis_Type = case_when(str_detect(group_var, "_") ~ "TS_Feature",
                                   group_var == "Combo" ~ "Combo",
                                   T ~ "Brain_Region"))

# Confirm that for the first repeat, each sample is assigned to the same fold every time
all_univariate_fold_assignments %>%
  filter(Repeat == 1) %>%
  group_by(Study, Comparison_Group, Sample_Index, Univariate_Feature_Set) %>%
  summarise(num_folds = length(unique(Fold))) %>%
  filter(num_folds != 1)


################################################################################
# Compile balanced accuracy results
################################################################################

# Load balanced accuracy data, or construct if needed
if (!file.exists(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_univariate_balanced_accuracy_all_folds.feather"))) {
  univariate_balanced_accuracy_all_folds_list <- list()
  # First iterate over each study/comparison group
  for (i in 1:nrow(study_group_df)) {
    dataset_ID <- study_group_df$Study[i]
    noise_proc <- study_group_df$Noise_Proc[i]
    noise_label = gsub("\\+", "_", noise_proc)
    comparison_group <- study_group_df$Comparison_Group[i]
    
    # Now iterate over each univariate feature set
    for (featset in univariate_feature_sets) {
      balacc_across_folds <- pyarrow_feather$read_feather(glue("{data_path}/{dataset_ID}/processed_data/{dataset_ID}_{comparison_group}_Univariate_{featset}_SVM_balanced_accuracy.feather"))
      balacc_across_folds$Study <- dataset_ID
      balacc_across_folds$Univariate_Feature_Set <- featset
      # Append to list
      univariate_balanced_accuracy_all_folds_list <- list.append(univariate_balanced_accuracy_all_folds_list, balacc_across_folds)
    }
  }
  
  # Combine the list results into a dataframe
  univariate_balanced_accuracy_all_folds <- do.call(plyr::rbind.fill, 
                                                    univariate_balanced_accuracy_all_folds_list) %>%
    mutate(Analysis_Type = case_when(str_detect(group_var, "_") ~ "TS_Feature",
                                     group_var == "Combo" ~ "Combo",
                                     T ~ "Brain_Region"))
  
  # Save to feather file
  feather::write_feather(univariate_balanced_accuracy_all_folds,
                         glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_univariate_balanced_accuracy_all_folds.feather"))
} else {
  univariate_balanced_accuracy_all_folds <- feather::read_feather(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_univariate_balanced_accuracy_all_folds.feather"))
}

# Aggregate the main results by repeat
univariate_balanced_accuracy_by_repeats <- univariate_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Analysis_Type, group_var, Repeat_Number) %>%
  summarise(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
            Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T))

# Aggregate the main results across all folds, independent of repeat
univariate_balanced_accuracy <- univariate_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Analysis_Type, group_var) %>%
  summarise(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
            Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T))

# Null results
if (!file.exists(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_univariate_null_balanced_accuracy_distributions.feather"))) {
  univariate_null_balanced_accuracy_list <- list()
  # First iterate over each study/comparison group
  for (i in 1:nrow(study_group_df)) {
    dataset_ID <- study_group_df$Study[i]
    noise_proc <- study_group_df$Noise_Proc[i]
    noise_label = gsub("\\+", "_", noise_proc)
    comparison_group <- study_group_df$Comparison_Group[i]
    
    # Now iterate over each univariate feature set
    for (featset in univariate_feature_sets) {
      null_dist <- pyarrow_feather$read_feather(glue("{data_path}/{dataset_ID}/processed_data/{dataset_ID}_{comparison_group}_Univariate_{featset}_SVM_null_balanced_accuracy_distributions.feather"))
      null_dist$Study <- dataset_ID
      null_dist$Univariate_Feature_Set <- featset
      # Append to list
      univariate_null_balanced_accuracy_list <- list.append(univariate_null_balanced_accuracy_list, null_dist)
    }
  }
  
  # Combine the list results into a dataframe
  univariate_null_balanced_accuracy <- do.call(plyr::rbind.fill, 
                                               univariate_null_balanced_accuracy_list) %>%
    mutate(Analysis_Type = case_when(str_detect(group_var, "_") ~ "TS_Feature",
                                     group_var == "Combo" ~ "Combo",
                                     T ~ "Brain_Region"))
  
  # Save to feather file
  feather::write_feather(univariate_null_balanced_accuracy,
                         glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_univariate_null_balanced_accuracy_distributions.feather"))
} else {
  univariate_null_balanced_accuracy <- feather::read_feather(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_univariate_null_balanced_accuracy_distributions.feather"))
}

# Calculate p-values based on empirical nulls
# Sanity check the p-value correction
if (!file.exists(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_univariate_empirical_p_values.feather"))) {
  univariate_split <- univariate_balanced_accuracy %>%
    group_by(Study, Comparison_Group, Analysis_Type, Univariate_Feature_Set, group_var) %>%
    group_split()
  
  univariate_p_values <- calculate_empirical_p_values(main_balanced_accuracy_split = univariate_split,
                                                      null_balanced_accuracy = univariate_null_balanced_accuracy)
  
  # Adjust p-values by group
  univariate_p_values <- univariate_p_values %>%
    group_by(Study, Comparison_Group, Univariate_Feature_Set, Analysis_Type) %>%
    mutate(p_value_BH = p.adjust(p_value, method="BH"))
  
  feather::write_feather(univariate_p_values, glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_univariate_empirical_p_values.feather"))
} else {
  univariate_p_values <- feather::read_feather(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_univariate_empirical_p_values.feather"))
}

# PAIRWISE
# Load balanced accuracy data, or construct if needed
if (!file.exists(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_pairwise_balanced_accuracy_all_folds.feather"))) {
  pairwise_balanced_accuracy_all_folds_list <- list()
  # First iterate over each study/comparison group
  for (i in 1:nrow(study_group_df)) {
    dataset_ID <- study_group_df$Study[i]
    noise_proc <- study_group_df$Noise_Proc[i]
    noise_label = gsub("\\+", "_", noise_proc)
    comparison_group <- study_group_df$Comparison_Group[i]
    balacc_across_folds <- pyarrow_feather$read_feather(glue("{data_path}/{dataset_ID}/processed_data/{dataset_ID}_{comparison_group}_Pairwise_{pairwise_feature_set}_SVM_balanced_accuracy.feather"))
    balacc_across_folds$Study <- dataset_ID
    balacc_across_folds$Pairwise_Feature_Set <- pairwise_feature_set
    # Append to list
    pairwise_balanced_accuracy_all_folds_list <- list.append(pairwise_balanced_accuracy_all_folds_list, balacc_across_folds)
  }
  
  # Combine the list results into a dataframe
  pairwise_balanced_accuracy_all_folds <- do.call(plyr::rbind.fill, 
                                                  pairwise_balanced_accuracy_all_folds_list) %>%
    mutate(Analysis_Type = "SPI")
  
  # Save to feather file
  feather::write_feather(pairwise_balanced_accuracy_all_folds,
                         glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_pairwise_balanced_accuracy_all_folds.feather"))
} else {
  pairwise_balanced_accuracy_all_folds <- feather::read_feather(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_pairwise_balanced_accuracy_all_folds.feather"))
}

# Aggregate the main results by repeat
pairwise_balanced_accuracy_by_repeats <- pairwise_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Pairwise_Feature_Set, Analysis_Type, group_var, Repeat_Number) %>%
  summarise(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
            Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T))

# Aggregate the main results across all folds, independent of repeat
pairwise_balanced_accuracy <- pairwise_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Pairwise_Feature_Set, Analysis_Type, group_var) %>%
  summarise(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
            Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T))

# Null results
if (!file.exists(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_pairwise_null_balanced_accuracy_distributions.feather"))) {
  pairwise_null_balanced_accuracy_list <- list()
  # First iterate over each study/comparison group
  for (i in 1:nrow(study_group_df)) {
    dataset_ID <- study_group_df$Study[i]
    noise_proc <- study_group_df$Noise_Proc[i]
    noise_label = gsub("\\+", "_", noise_proc)
    comparison_group <- study_group_df$Comparison_Group[i]
    
    null_dist <- pyarrow_feather$read_feather(glue("{data_path}/{dataset_ID}/processed_data/{dataset_ID}_{comparison_group}_Pairwise_{pairwise_feature_set}_SVM_null_balanced_accuracy_distributions.feather"))
    null_dist$Study <- dataset_ID
    null_dist$Pairwise_Feature_Set <- pairwise_feature_set
    # Append to list
    pairwise_null_balanced_accuracy_list <- list.append(pairwise_null_balanced_accuracy_list, null_dist)
  }
  
  # Combine the list results into a dataframe
  pairwise_null_balanced_accuracy <- do.call(plyr::rbind.fill, 
                                             pairwise_null_balanced_accuracy_list) %>%
    mutate(Analysis_Type = "SPI")
  
  # Save to feather file
  feather::write_feather(pairwise_null_balanced_accuracy,
                         glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_pairwise_null_balanced_accuracy_distributions.feather"))
} else {
  pairwise_null_balanced_accuracy <- feather::read_feather(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_pairwise_null_balanced_accuracy_distributions.feather"))
}

# Calculate p-values based on empirical nullss
if (!file.exists(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_pairwise_empirical_p_values.feather"))) {
  pairwise_split <- pairwise_balanced_accuracy %>%
    group_by(Study, Comparison_Group, Analysis_Type, Pairwise_Feature_Set, group_var) %>%
    group_split()
  
  pairwise_p_values <- calculate_empirical_p_values(main_balanced_accuracy_split = pairwise_split,
                                                    null_balanced_accuracy = pairwise_null_balanced_accuracy)
  
  # Adjust p-values by group
  pairwise_p_values <- pairwise_p_values %>%
    group_by(Study, Comparison_Group, Pairwise_Feature_Set, Analysis_Type) %>%
    mutate(p_value_BH = p.adjust(p_value, method="BH"))
  
  feather::write_feather(pairwise_p_values, glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_pairwise_empirical_p_values.feather"))
}

# UNIVARIATE + PAIRWISE COMBO
# Load balanced accuracy data, or construct if needed
if (!file.exists(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_combo_univariate_pairwise_balanced_accuracy_all_folds.feather"))) {
  combo_univariate_pairwise_balanced_accuracy_all_folds_list <- list()
  # First iterate over each study/comparison group
  for (i in 1:nrow(study_group_df)) {
    dataset_ID <- study_group_df$Study[i]
    noise_proc <- study_group_df$Noise_Proc[i]
    noise_label = gsub("\\+", "_", noise_proc)
    comparison_group <- study_group_df$Comparison_Group[i]
    
    # Now iterate over each univariate feature set
    for (featset in univariate_feature_sets) {
      balacc_across_folds <- pyarrow_feather$read_feather(glue("{data_path}/{dataset_ID}/processed_data/{dataset_ID}_{comparison_group}_Univariate_{featset}_Pairwise_{pairwise_feature_set}_SVM_balanced_accuracy.feather"))
      balacc_across_folds$Study <- dataset_ID
      balacc_across_folds$Univariate_Feature_Set <- featset
      balacc_across_folds$Pairwise_Feature_Set <- pairwise_feature_set
      # Append to list
      combo_univariate_pairwise_balanced_accuracy_all_folds_list <- list.append(combo_univariate_pairwise_balanced_accuracy_all_folds_list, balacc_across_folds)
    }
  }
  
  # Combine the list results into a dataframe
  combo_univariate_pairwise_balanced_accuracy_all_folds <- do.call(plyr::rbind.fill, 
                                                                   combo_univariate_pairwise_balanced_accuracy_all_folds_list) %>%
    mutate(Analysis_Type = "SPI_Combo")
  
  # Save to feather file
  feather::write_feather(combo_univariate_pairwise_balanced_accuracy_all_folds,
                         glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_combo_univariate_pairwise_balanced_accuracy_all_folds.feather"))
} else {
  combo_univariate_pairwise_balanced_accuracy_all_folds <- feather::read_feather(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_combo_univariate_pairwise_balanced_accuracy_all_folds.feather"))
}

# Aggregate the main results by repeat
combo_univariate_pairwise_balanced_accuracy_by_repeats <- combo_univariate_pairwise_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Pairwise_Feature_Set, Analysis_Type, group_var, Repeat_Number) %>%
  summarise(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
            Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T))

# Aggregate the main results across all folds, independent of repeat
combo_univariate_pairwise_balanced_accuracy <- combo_univariate_pairwise_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Pairwise_Feature_Set, Analysis_Type, group_var) %>%
  summarise(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
            Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T))

# Null results
if (!file.exists(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_combo_univariate_pairwise_null_balanced_accuracy_distributions.feather"))) {
  combo_univariate_pairwise_null_balanced_accuracy_list <- list()
  # First iterate over each study/comparison group
  for (i in 1:nrow(study_group_df)) {
    dataset_ID <- study_group_df$Study[i]
    noise_proc <- study_group_df$Noise_Proc[i]
    noise_label = gsub("\\+", "_", noise_proc)
    comparison_group <- study_group_df$Comparison_Group[i]
    
    # Now iterate over each univariate feature set
    for (featset in univariate_feature_sets) {
      null_dist <- pyarrow_feather$read_feather(glue("{data_path}/{dataset_ID}/processed_data/{dataset_ID}_{comparison_group}_Univariate_{featset}_Pairwise_{pairwise_feature_set}_SVM_null_balanced_accuracy_distributions.feather"))
      null_dist$Study <- dataset_ID
      null_dist$Univariate_Feature_Set <- featset
      null_dist$Pairwise_Feature_Set <- pairwise_feature_set
      # Append to list
      combo_univariate_pairwise_null_balanced_accuracy_list <- list.append(combo_univariate_pairwise_null_balanced_accuracy_list, null_dist)
    }
  }
  
  # Combine the list results into a dataframe
  combo_univariate_pairwise_null_balanced_accuracy <- do.call(plyr::rbind.fill, 
                                                              combo_univariate_pairwise_null_balanced_accuracy_list) %>%
    mutate(Analysis_Type = "SPI_Combo")
  
  # Save to feather file
  feather::write_feather(combo_univariate_pairwise_null_balanced_accuracy,
                         glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_combo_univariate_pairwise_null_balanced_accuracy_distributions.feather"))
} else {
  combo_univariate_pairwise_null_balanced_accuracy <- feather::read_feather(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_combo_univariate_pairwise_null_balanced_accuracy_distributions.feather"))
}

# Calculate p-values based on empirical nulls
if (!file.exists(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_combo_univariate_pairwise_empirical_p_values.feather"))) {
  combo_univariate_pairwise_split <- combo_univariate_pairwise_balanced_accuracy %>%
    group_by(Study, Comparison_Group, Analysis_Type, Univariate_Feature_Set, Pairwise_Feature_Set, group_var) %>%
    group_split()
  
  combo_univariate_pairwise_p_values <- calculate_empirical_p_values(main_balanced_accuracy_split = combo_univariate_pairwise_split,
                                                                     null_balanced_accuracy = combo_univariate_pairwise_null_balanced_accuracy)
  
  # Adjust p-values by group
  combo_univariate_pairwise_p_values <- combo_univariate_pairwise_p_values %>%
    group_by(Study, Comparison_Group, Univariate_Feature_Set, Pairwise_Feature_Set, Analysis_Type) %>%
    mutate(p_value_BH = p.adjust(p_value, method="BH"))
  
  feather::write_feather(combo_univariate_pairwise_p_values, glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_combo_univariate_pairwise_empirical_p_values.feather"))
}

bivar <- pyarrow_feather$read_feather(glue("{data_path}/ABIDE_ASD/processed_data/ABIDE_ASD_FC1000_pyspi14_filtered_zscored.feather"))

ABIDE_metadata <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_sample_metadata.feather")

data_to_analyse <- bivar %>%
  distinct(Sample_ID) %>%
  left_join(., ABIDE_metadata)

data_to_analyse %>%
  group_by(Diagnosis) %>%
  summarise(N = n(),
            Num_Female = sum(Sex=="F"),
            Perc_Female = round(100*Num_Female/N, 1),
            Age_Mean = round(mean(as.numeric(Age), na.rm=T), 1),
            Age_SD = round(sd(as.numeric(Age)), 1))
