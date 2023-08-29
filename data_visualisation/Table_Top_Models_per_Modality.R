################################################################################
# Load libraries
################################################################################

python_to_use <- "/Users/abry4213/anaconda3/envs/pyspi/bin/python3"
reticulate::use_python(python_to_use)
library(reticulate)

# Import pyarrow.feather as pyarrow_feather
pyarrow_feather <- import("pyarrow.feather")

library(tidyverse)
library(knitr)
library(kableExtra)
library(glue)
library(rstatix)

################################################################################
# Define study/data paths
################################################################################

# Filter to just those samples used in each study
UCLA_CNP_subjects_to_keep <- pyarrow_feather$read_feather("~/data/UCLA_CNP/processed_data/UCLA_CNP_filtered_sample_info_AROMA_2P_GMR_catch24_pyspi14.feather")
ABIDE_ASD_subjects_to_keep <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/processed_data/ABIDE_ASD_filtered_sample_info_FC1000_catch24_pyspi14.feather")

UCLA_CNP_study_metadata <- pyarrow_feather$read_feather("~/data/UCLA_CNP/study_metadata/UCLA_CNP_sample_metadata.feather") %>%
  semi_join(., UCLA_CNP_subjects_to_keep) %>%
  mutate(Age = as.numeric(Age),
         Study = "UCLA_CNP")
ABIDE_ASD_study_metadata <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_sample_metadata.feather") %>%
  semi_join(., ABIDE_ASD_subjects_to_keep) %>%
  mutate(Age = as.numeric(Age),
         Study = "ABIDE")

# Load univariate classification results across all folds
univariate_balanced_accuracy_AUC_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_mixedsigmoid_scaler_balanced_accuracy_AUC_all_folds.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set, kernel==SVM_kernel)
# Compute mean + SD performance across all folds
univariate_balanced_accuracy <- univariate_balanced_accuracy_AUC_all_folds %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Analysis_Type, group_var, kernel) %>%
  reframe(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
          Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T),
          ROC_AUC_Across_Folds = mean(ROC_AUC, na.rm=T),
          ROC_AUC_Across_Folds_SD = sd(ROC_AUC, na.rm=T))

# Load pairwise stats data
pairwise_balanced_accuracy_AUC_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_mixedsigmoid_scaler_balanced_accuracy_AUC_all_folds.feather"))
# Compute mean + SD performance across all folds
pairwise_balanced_accuracy <- pairwise_balanced_accuracy_AUC_all_folds %>%
  group_by(Study, Comparison_Group, Pairwise_Feature_Set, Analysis_Type, group_var) %>%
  reframe(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
          Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T),
          ROC_AUC_Across_Folds = mean(ROC_AUC, na.rm=T),
          ROC_AUC_Across_Folds_SD = sd(ROC_AUC, na.rm=T))

# Load combined pairwise + univariate data
combo_univariate_pairwise_balanced_accuracy_AUC_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combined_univariate_pairwise_mixedsigmoid_scaler_balanced_accuracy_AUC_all_folds.feather")) %>%
  mutate(Analysis_Type = "SPI_Univariate_Combo")
# Compute mean + SD performance across all folds
combo_univariate_pairwise_balanced_accuracy <- combo_univariate_pairwise_balanced_accuracy_AUC_all_folds %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Pairwise_Feature_Set, Analysis_Type, group_var) %>%
  reframe(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
          Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T),
          ROC_AUC_Across_Folds = mean(ROC_AUC, na.rm=T),
          ROC_AUC_Across_Folds_SD = sd(ROC_AUC, na.rm=T))


# Combine all the datasets into one
full_balacc_res <- do.call(plyr::rbind.fill, list(univariate_balanced_accuracy,
                                                  pairwise_balanced_accuracy,
                                                  combo_univariate_pairwise_balanced_accuracy))


full_balacc_res %>%
  dplyr::group_by(Study, Comparison_Group, Analysis_Type) %>%
  dplyr::filter(Balanced_Accuracy_Across_Folds == max(Balanced_Accuracy_Across_Folds)) %>%
  ungroup() %>%
  mutate(val_to_show = glue("{group_var}: {100*round(Balanced_Accuracy_Across_Folds,3)} +/- {100*round(Balanced_Accuracy_Across_Folds_SD,3)}")) %>%
  pivot_wider(id_cols = c(Analysis_Type),
              names_from = Comparison_Group,
              values_from = val_to_show)