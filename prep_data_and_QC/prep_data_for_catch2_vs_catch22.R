#-------------------------------------------------------------------------------
# Define paths
#-------------------------------------------------------------------------------

python_to_use <- "~/.conda/envs/pyspi/bin/python3"
# python_to_use <- "/Users/abry4213/opt/anaconda3/envs/pyspi/bin/python3"
pairwise_feature_set <- "pyspi14"
github_dir <- "~/github/"
data_path <- "~/data/"
scaler <- "robustsigmoid"
sample_metadata_file <- "UCLA_CNP_sample_metadata.feather"
noise_proc <- "AROMA+2P+GMR"
output_data_path <- glue::glue("{data_path}/TS_feature_manuscript/catch2_catch22/")
TAF::mkdir(output_data_path)
univariate_feature_sets <- c("catch22", "catch2")

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
library(R.matlab)

# Import pyarrow.feather as pyarrow_feather
pyarrow_feather <- import("pyarrow.feather")

# Prep metadata for just the 166 participants included in this analysis
metadata <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP/study_metadata/{sample_metadata_file}")) %>%
  dplyr::select(-Study) %>%
  mutate(Age = as.numeric(Age))

# Load catch22 feature values
catch22_feature_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP/processed_data/UCLA_CNP_AROMA_2P_GMR_catch22_filtered.feather")) %>%
  left_join(., metadata) %>%
  filter(Diagnosis %in% c("Schizophrenia", "Control")) 

# Identify subjects to use
subjects_to_use <- unique(catch22_feature_values$Sample_ID)

# Filter metadata and save
metadata <- metadata %>% filter(Sample_ID %in% subjects_to_use)
save(metadata, file=glue("{output_data_path}/metadata.Rda"))

# Prep mean framewise displacement data
mFD_data <- read.table(glue("{data_path}/UCLA_CNP/movement_data/UCLA_CNP_mFD.txt"),
                       sep=",")
colnames(mFD_data) <- c("Sample_ID", "mFD_Jenkinson", "mFD_Power", "mFD_VanDijk")
# Filter to samples, link with metadata, and save
mFD_data <- mFD_data %>%
  filter(Sample_ID %in% subjects_to_use) %>%
  left_join(., metadata) %>%
  dplyr::select(-mFD_Jenkinson, -mFD_VanDijk)
save(mFD_data, file=glue("{output_data_path}/movement_data.Rda"))

# Prep number of subjects retained with each movement threshold
FD_data <- as.data.frame(R.matlab::readMat(glue("{data_path}/UCLA_CNP/movement_data/UCLA_CNP_all_FD.mat"))[[1]])
colnames(FD_data) <- c("Sample_ID", "FD_Jenkinson", "FD_Power", "FD_VanDijk")

FD_data <- unnest(FD_data) %>%
  unnest(cols = c(Sample_ID, FD_Jenkinson, FD_Power, FD_VanDijk)) %>%
  # Filter to subjects used here
  filter(Sample_ID %in% subjects_to_use) %>%
  group_by(Sample_ID) %>%
  mutate(Frame_Number = row_number()) %>%
  ungroup()
colnames(FD_data) <- c("Sample_ID", "FD_Jenkinson", "FD_Power", "FD_VanDijk", "Frame_Number")

FD_data <- left_join(FD_data, metadata)

# Filter to subjects used here

# No threshold
no_threshold_counts <- FD_data %>%
  distinct(Sample_ID, Diagnosis) %>%
  group_by(Diagnosis) %>%
  count() %>%
  mutate(Threshold = "None")

# Lenient threshold
lenient_threshold_counts <- mFD_data %>%
  filter(mFD_Power < 0.55) %>%
  group_by(Diagnosis) %>%
  count() %>%
  mutate(Threshold = "Lenient")

# Stringent threshold
stringent_threshold_counts <- FD_data %>%
  # No frames with FD > 5mm
  group_by(Sample_ID) %>%
  filter(max(FD_Power) < 5) %>%
  # <20% of frames above 0.2mm
  filter(sum(FD_Power > 0.2) < 0.2*n()) %>%
  # Mean FD < 0.25mm
  filter(mean(FD_Power) < 0.25) %>%
  ungroup() %>%
  distinct(Sample_ID, Diagnosis) %>%
  group_by(Diagnosis) %>%
  count() %>%
  mutate(Threshold = "Stringent")

movement_threshold_counts <- plyr::rbind.fill(list(no_threshold_counts,
                                                   lenient_threshold_counts,
                                                   stringent_threshold_counts))

save(movement_threshold_counts, file=glue("{output_data_path}/movement_threshold_counts.Rda"))


# Load catch2 feature values
catch2_feature_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP/processed_data/UCLA_CNP_AROMA_2P_GMR_catch2_filtered.feather"))%>%
  left_join(., metadata) %>%
  filter(Diagnosis %in% c("Schizophrenia", "Control")) 

# Merge catch2 + catch22 feature values and save
all_feature_values <- plyr::rbind.fill(catch22_feature_values, catch2_feature_values) %>%
  filter(!is.na(Diagnosis))
save(all_feature_values, file=glue("{output_data_path}/all_feature_values.Rda"))

# Load catch22 SVM balanced accuracy by repeat
catch22_SVM_balanced_accuracy_by_repeat <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP/processed_data/UCLA_CNP_Schizophrenia_Univariate_catch22_robustsigmoid_scaler_SVM_sample_predictions.feather")) %>%
  filter(Comparison_Group == "Schizophrenia") %>%
  mutate(Univariate_Feature_Set = "catch22",
         CV_Predicted_Diagnosis = factor(CV_Predicted_Diagnosis, levels = c("Control", "Schizophrenia")),
         Diagnosis = factor(Diagnosis, levels = c("Control", "Schizophrenia"))) %>%
  group_by(Analysis_Type, Univariate_Feature_Set, group_var, Repeat_Number) %>%
  summarise(Repeat_Balanced_Accuracy = caret::confusionMatrix(data = CV_Predicted_Diagnosis,
                                                              reference = Diagnosis)$byClass[["Balanced Accuracy"]])

# Load catch2 SVM individual subject predictions
catch2_SVM_balanced_accuracy_by_repeat <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP/processed_data/UCLA_CNP_Schizophrenia_Univariate_catch2_robustsigmoid_scaler_SVM_sample_predictions.feather")) %>%
  filter(Comparison_Group == "Schizophrenia") %>%
  mutate(Univariate_Feature_Set = "FTM",
         CV_Predicted_Diagnosis = factor(CV_Predicted_Diagnosis, levels = c("Control", "Schizophrenia")),
         Diagnosis = factor(Diagnosis, levels = c("Control", "Schizophrenia"))) %>%
  group_by(Analysis_Type, Univariate_Feature_Set, group_var, Repeat_Number) %>%
  summarise(Repeat_Balanced_Accuracy = caret::confusionMatrix(data = CV_Predicted_Diagnosis,
                                                              reference = Diagnosis)$byClass[["Balanced Accuracy"]])

# Merge catch2 + catch22 subject predictions
balanced_accuracy_by_repeats <- plyr::rbind.fill(catch22_SVM_balanced_accuracy_by_repeat,
                                                 catch2_SVM_balanced_accuracy_by_repeat)
save(balanced_accuracy_by_repeats, file=glue("{output_data_path}/balanced_accuracy_by_repeats.Rda"))

# Load catch22 SVM coefficients
catch22_SVM_coefficients <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP/processed_data/UCLA_CNP_Schizophrenia_Univariate_catch22_robustsigmoid_scaler_SVM_fold_SVM_coefficients.feather")) %>%
  filter(Comparison_Group == "Schizophrenia") %>%
  mutate(Univariate_Feature_Set = "catch22") %>%
  group_by(Analysis_Type, Univariate_Feature_Set, group_var, `Feature Name`, `Repeat Number`) %>%
  summarise(mean_coef = mean(Coefficient, na.rm=T)) %>%
  group_by(Analysis_Type, Univariate_Feature_Set, group_var, `Feature Name`) %>%
  summarise(coef = mean(mean_coef, na.rm=T))

# Load catch2 SVM coefficients
catch2_SVM_coefficients <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP/processed_data/UCLA_CNP_Schizophrenia_Univariate_catch2_robustsigmoid_scaler_SVM_fold_SVM_coefficients.feather")) %>%
  filter(Comparison_Group == "Schizophrenia") %>%
  mutate(Univariate_Feature_Set = "FTM") %>%
  group_by(Analysis_Type, Univariate_Feature_Set, group_var, `Feature Name`, `Repeat Number`) %>%
  summarise(mean_coef = mean(Coefficient, na.rm=T)) %>%
  group_by(Analysis_Type, Univariate_Feature_Set, group_var, `Feature Name`) %>%
  summarise(coef = mean(mean_coef, na.rm=T))

# Merge catch2 + catch22 SVM coefficients
SVM_coefficients <- plyr::rbind.fill(catch22_SVM_coefficients,
                                                    catch2_SVM_coefficients)
save(SVM_coefficients, file=glue("{output_data_path}/SVM_coefficients.Rda"))

# Load catch22 null balanced accuracy res
catch22_null_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP/processed_data/UCLA_CNP_Schizophrenia_Univariate_catch22_robustsigmoid_scaler_SVM_null_balanced_accuracy_distributions.feather"))
catch22_null_balanced_accuracy_all_folds$Univariate_Feature_Set <- "catch22"

# Load catch2 null balanced accuracy res
catch2_null_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP/processed_data/UCLA_CNP_Schizophrenia_Univariate_catch2_robustsigmoid_scaler_SVM_null_balanced_accuracy_distributions.feather"))
catch2_null_balanced_accuracy_all_folds$Univariate_Feature_Set <- "FTM"

# Merge catch2 + catch22 null balanced accuracy res
null_balanced_accuracy_all_folds <- plyr::rbind.fill(catch22_null_balanced_accuracy_all_folds,
                                                    catch2_null_balanced_accuracy_all_folds) %>%
  dplyr::select(-Comparison_Group, -Scaling_Type)
save(null_balanced_accuracy_all_folds, file=glue("{output_data_path}/Null_SVM_balanced_accuracy.Rda"))

# Load catch22 movement-based SVM results
catch22_movement_SVM_res <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP/processed_data/UCLA_CNP_Schizophrenia_Univariate_catch22_robustsigmoid_scaler_movement_SVM_balanced_accuracy.feather"))
catch22_movement_SVM_res$Univariate_Feature_Set <- "catch22"

# Load catch2 movement-based SVM results
catch2_movement_SVM_res <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP/processed_data/UCLA_CNP_Schizophrenia_Univariate_catch2_robustsigmoid_scaler_movement_SVM_balanced_accuracy.feather"))
catch2_movement_SVM_res$Univariate_Feature_Set <- "FTM"

# Merge movement-based SVM results and save
all_movement_SVM_res <- plyr::rbind.fill(catch22_movement_SVM_res,
                                         catch2_movement_SVM_res) %>%
  dplyr::select(Threshold_Type, Univariate_Feature_Set, Analysis_Type, group_var, Repeat_Number, Fold, Balanced_Accuracy)
save(all_movement_SVM_res, file=glue("{output_data_path}/movement_based_SVM_results.Rda"))