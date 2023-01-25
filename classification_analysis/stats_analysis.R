python_to_use <- "/headnode1/abry4213/.conda/envs/pyspi/bin/python3"
reticulate::use_python(python_to_use)
library(feather)
library(tidyverse)
library(reticulate)

# Import pyarrow.feather as pyarrow_feather
pyarrow_feather <- import("pyarrow.feather")

################################################################################
# Univariate catch22 results
################################################################################

# Main results

# UCLA CNP ADHD
UCLA_CNP_ADHD_univariate_balanced_accuracy_folds_repeats <- pyarrow_feather$read_feather("~/data/UCLA_CNP/processed_data/UCLA_CNP_ADHD_Univariate_catch22_SVM_balanced_accuracy.feather")
UCLA_CNP_ADHD_univariate_balanced_accuracy_folds_repeats$Study <- "UCLA_CNP"

# UCLA CNP Bipolar
UCLA_CNP_Bipolar_univariate_balanced_accuracy_folds_repeats <- pyarrow_feather$read_feather("~/data/UCLA_CNP/processed_data/UCLA_CNP_Bipolar_Univariate_catch22_SVM_balanced_accuracy.feather")
UCLA_CNP_Bipolar_univariate_balanced_accuracy_folds_repeats$Study <- "UCLA_CNP"

# UCLA CNP Schizophrenia
UCLA_CNP_Schizophrenia_univariate_balanced_accuracy_folds_repeats <- pyarrow_feather$read_feather("~/data/UCLA_CNP/processed_data/UCLA_CNP_Schizophrenia_Univariate_catch22_SVM_balanced_accuracy.feather")
UCLA_CNP_Schizophrenia_univariate_balanced_accuracy_folds_repeats$Study <- "UCLA_CNP"

# ABIDE ASD
ABIDE_ASD_univariate_balanced_accuracy_folds_repeats <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/processed_data/ABIDE_ASD_ASD_Univariate_catch22_SVM_balanced_accuracy.feather")
ABIDE_ASD_univariate_balanced_accuracy_folds_repeats$Study <- "ABIDE_ASD"

# Merge the main results
univariate_balanced_accuracy_folds_repeats <- do.call(plyr::rbind.fill, list(UCLA_CNP_ADHD_univariate_balanced_accuracy_folds_repeats,
                                                                             UCLA_CNP_Bipolar_univariate_balanced_accuracy_folds_repeats,
                                                                             UCLA_CNP_Schizophrenia_univariate_balanced_accuracy_folds_repeats,
                                                                             ABIDE_ASD_univariate_balanced_accuracy_folds_repeats)) %>%
  dplyr::select(-index)

# Aggregate the main results by repeat
univariate_balanced_accuracy_repeats <- univariate_balanced_accuracy_folds_repeats %>%
  group_by(Study, Comparison_Group, group_var, Repeat_Number) %>%
  summarise(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
            Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T))

# Aggregate the main results across all folds, independent of repeat
univariate_balanced_accuracy <- univariate_balanced_accuracy_folds_repeats %>%
  group_by(Study, Comparison_Group, group_var) %>%
  summarise(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
            Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T))