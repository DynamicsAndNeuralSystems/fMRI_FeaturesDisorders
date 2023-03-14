#-------------------------------------------------------------------------------
# Define paths
#-------------------------------------------------------------------------------

python_to_use <- "~/.conda/envs/pyspi/bin/python3"
# python_to_use <- "/Users/abry4213/opt/anaconda3/envs/pyspi/bin/python3"
pairwise_feature_set <- "pyspi14"
github_dir <- "~/github/"
data_path <- "~/data/"
scaler <- "standard"
sample_metadata_file <- "UCLA_CNP_sample_metadata.feather"
noise_proc <- "AROMA+2P+GMR"
output_data_path <- glue::glue("{data_path}/TS_feature_manuscript/FTM_vs_catch24/")
TAF::mkdir(output_data_path)
univariate_feature_sets <- c("catch22", "FTM", "catch24")

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
library(cowplot)
library(patchwork)
theme_set(theme_cowplot())

# Import pyarrow.feather as pyarrow_feather
pyarrow_feather <- import("pyarrow.feather")

# Prep metadata for just the 166 participants included in this analysis
metadata <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP/study_metadata/{sample_metadata_file}")) %>%
  dplyr::select(-Study) %>%
  mutate(Age = as.numeric(Age))

# Load catch24 feature values
catch24_feature_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP/processed_data/UCLA_CNP_AROMA_2P_GMR_catch24_filtered.feather")) %>%
  left_join(., metadata) %>%
  filter(Diagnosis %in% c("Schizophrenia", "Control"))

# Identify subjects to use
subjects_to_use <- unique(catch24_feature_values$Sample_ID)

# Filter metadata and save
metadata <- metadata %>% filter(Sample_ID %in% subjects_to_use)
save(metadata, file=glue("{output_data_path}/metadata.Rda"))

# Load FTM feature values
FTM_feature_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP/processed_data/UCLA_CNP_AROMA_2P_GMR_FTM_filtered.feather"))%>%
  left_join(., metadata) %>%
  filter(Diagnosis %in% c("Schizophrenia", "Control"),
         Sample_ID %in% subjects_to_use) 

# Merge FTM + catch24 feature values and save
all_feature_values <- plyr::rbind.fill(catch24_feature_values, FTM_feature_values) %>%
  filter(!is.na(Diagnosis))
save(all_feature_values, file=glue("{output_data_path}/all_feature_values.Rda"))

################################################################################
# Balanced accuracy by repeat, averaging first by test fold balanced accuracy
################################################################################

# Load catch24 SVM balanced accuracy by repeat
catch24_SVM_balanced_accuracy_by_repeat <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP/processed_data/UCLA_CNP_Schizophrenia_Univariate_catch24_{scaler}_scaler_SVM_balanced_accuracy.feather")) %>%
  filter(Comparison_Group == "Schizophrenia") %>%
  mutate(Univariate_Feature_Set = "catch24",
         Repeat_Number = Repeat_Number + 1) %>%
  group_by(Analysis_Type, Univariate_Feature_Set, group_var, Repeat_Number) %>%
  summarise(Repeat_Balanced_Accuracy = mean(Balanced_Accuracy, na.rm=T))

# Load catch22 SVM balanced accuracy by repeat
catch22_SVM_balanced_accuracy_by_repeat <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP/processed_data/UCLA_CNP_Schizophrenia_Univariate_catch22_{scaler}_scaler_SVM_balanced_accuracy.feather")) %>%
  filter(Comparison_Group == "Schizophrenia") %>%
  mutate(Univariate_Feature_Set = "catch22",
         Repeat_Number = Repeat_Number + 1) %>%
  group_by(Analysis_Type, Univariate_Feature_Set, group_var, Repeat_Number) %>%
  summarise(Repeat_Balanced_Accuracy = mean(Balanced_Accuracy, na.rm=T))

# Load FTM SVM balanced accuracy
FTM_SVM_balanced_accuracy_by_repeat <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP/processed_data/UCLA_CNP_Schizophrenia_Univariate_FTM_{scaler}_scaler_SVM_balanced_accuracy.feather")) %>%
  filter(Comparison_Group == "Schizophrenia") %>%
  mutate(Univariate_Feature_Set = "FTM",
         Repeat_Number = Repeat_Number + 1) %>%
  group_by(Analysis_Type, Univariate_Feature_Set, group_var, Repeat_Number) %>%
  summarise(Repeat_Balanced_Accuracy = mean(Balanced_Accuracy, na.rm=T))

# Merge FTM + catch22 + catch24 balanced accuracy
balanced_accuracy_by_repeats <- plyr::rbind.fill(catch24_SVM_balanced_accuracy_by_repeat,
                                                 catch22_SVM_balanced_accuracy_by_repeat,
                                                 FTM_SVM_balanced_accuracy_by_repeat)
save(balanced_accuracy_by_repeats, file=glue("{output_data_path}/balanced_accuracy_by_repeats.Rda"))



################################################################################
# Null distributions
################################################################################

# Load catch24 null balanced accuracy res
catch24_null_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP/processed_data/UCLA_CNP_Schizophrenia_Univariate_catch24_{scaler}_scaler_SVM_null_balanced_accuracy_distributions.feather"))
catch24_null_balanced_accuracy_all_folds$Univariate_Feature_Set <- "catch24"

# Load catch22 null balanced accuracy res
catch22_null_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP/processed_data/UCLA_CNP_Schizophrenia_Univariate_catch22_{scaler}_scaler_SVM_null_balanced_accuracy_distributions.feather"))
catch22_null_balanced_accuracy_all_folds$Univariate_Feature_Set <- "catch22"

# Load FTM null balanced accuracy res
FTM_null_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP/processed_data/UCLA_CNP_Schizophrenia_Univariate_FTM_{scaler}_scaler_SVM_null_balanced_accuracy_distributions.feather"))
FTM_null_balanced_accuracy_all_folds$Univariate_Feature_Set <- "FTM"

# Merge FTM + catch22 + catch24 null balanced accuracy res
null_balanced_accuracy_all_folds <- plyr::rbind.fill(catch24_null_balanced_accuracy_all_folds,
                                                     catch22_null_balanced_accuracy_all_folds,
                                                    FTM_null_balanced_accuracy_all_folds) %>%
  dplyr::select(-Comparison_Group, -Scaling_Type)
save(null_balanced_accuracy_all_folds, file=glue("{output_data_path}/Null_SVM_balanced_accuracy.Rda"))
