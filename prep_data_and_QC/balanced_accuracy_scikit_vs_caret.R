#-------------------------------------------------------------------------------
# Define paths
#-------------------------------------------------------------------------------

# python_to_use <- "~/.conda/envs/pyspi/bin/python3"
python_to_use <- "/Users/abry4213/opt/anaconda3/envs/pyspi/bin/python3"
univariate_feature_sets <- "catch22"
pairwise_feature_set <- "pyspi14"
github_dir <- "~/github/"
data_path <- "~/data/"
scaler <- "robustsigmoid"
sample_metadata_file <- "UCLA_CNP_sample_metadata.feather"
noise_proc <- "AROMA+2P+GMR"

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

# Identify subjects to use
catch24_feature_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP/processed_data/UCLA_CNP_AROMA_2P_GMR_catch24_filtered.feather"))
subjects_to_use <- unique(catch24_feature_values$Sample_ID)

# Filter metadata and save
metadata <- metadata %>% filter(Sample_ID %in% subjects_to_use)

################################################################################
# Fold-wise
################################################################################

# Balanced accuracy by fold, scikit-learn
catch22_SVM_balanced_accuracy_by_fold_scikit <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP/processed_data/UCLA_CNP_Schizophrenia_Univariate_catch22_mixedsigmoid_scaler_SVM_balanced_accuracy.feather")) %>%
  filter(Comparison_Group == "Schizophrenia") %>%
  mutate(Univariate_Feature_Set = "catch22",
         Repeat_Number = Repeat_Number + 1,
         Fold_Accuracy_Method = "scikit_learn")

# Balanced accuracy by fold, caret
fold_wise_sample_predictions <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP/processed_data/UCLA_CNP_Schizophrenia_Univariate_catch22_mixedsigmoid_scaler_SVM_sample_predictions.feather")) %>%
  filter(Comparison_Group == "Schizophrenia") %>%
  mutate(Univariate_Feature_Set = "catch22",
         CV_Predicted_Diagnosis = factor(CV_Predicted_Diagnosis, levels = c("Control", "Schizophrenia")),
         Diagnosis = factor(Diagnosis, levels = c("Control", "Schizophrenia"))) %>%
  dplyr::select(Sample_ID, Repeat_Number, Analysis_Type, Univariate_Feature_Set, group_var, Diagnosis, CV_Predicted_Diagnosis)

fold_assignments <-  pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP/processed_data/UCLA_CNP_Schizophrenia_Univariate_catch22_mixedsigmoid_scaler_SVM_fold_assignments.feather")) %>%
  dplyr::select(Sample_ID, Repeat, Fold, Analysis_Type, group_var) %>%
  dplyr::rename("Repeat_Number" = "Repeat") 

catch22_SVM_balanced_accuracy_by_fold_caret <- left_join(fold_wise_sample_predictions,
                                                         fold_assignments) %>%
  group_by(Analysis_Type, Univariate_Feature_Set, group_var, Repeat_Number, Fold) %>%
  summarise(Balanced_Accuracy = caret::confusionMatrix(data = CV_Predicted_Diagnosis,
                                                       reference = Diagnosis)$byClass[["Balanced Accuracy"]],
            Accuracy = sum(CV_Predicted_Diagnosis==Diagnosis)/n()) %>%
  mutate(Univariate_Feature_Set = "catch22",
         Fold_Accuracy_Method = "caret")

################################################################################
# Repeat-wise
################################################################################

# Balanced accuracy by repeat, averaging first by test fold balanced accuracy automatically within scikit-learn
catch22_SVM_balanced_accuracy_by_fold_repeat_scikit <- catch22_SVM_balanced_accuracy_by_fold_scikit %>%
  group_by(Analysis_Type, Univariate_Feature_Set, group_var, Repeat_Number) %>%
  summarise(Repeat_Balanced_Accuracy = mean(Balanced_Accuracy, na.rm=T)) %>%
  mutate(Aggregate_Type = "fold_repeat_scikit")

# Balanced accuracy by repeat, averaging first by test fold balanced accuracy manually with caret
catch22_SVM_balanced_accuracy_by_fold_repeat_caret <- catch22_SVM_balanced_accuracy_by_fold_caret %>%
  group_by(Analysis_Type, Univariate_Feature_Set, group_var, Repeat_Number) %>%
  summarise(Repeat_Balanced_Accuracy = mean(Balanced_Accuracy),
            Repeat_Accuracy = mean(Accuracy)) %>%
  mutate(Aggregate_Type = "fold_repeat_caret")

# Balanced accuracy by repeat, without averaging by repeat
catch22_SVM_balanced_accuracy_by_repeat_caret <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP/processed_data/UCLA_CNP_Schizophrenia_Univariate_catch22_mixedsigmoid_scaler_SVM_sample_predictions.feather")) %>%
  filter(Comparison_Group == "Schizophrenia") %>%
  mutate(Univariate_Feature_Set = "catch22",
         CV_Predicted_Diagnosis = factor(CV_Predicted_Diagnosis, levels = c("Control", "Schizophrenia")),
         Diagnosis = factor(Diagnosis, levels = c("Control", "Schizophrenia"))) %>%
  group_by(Analysis_Type, Univariate_Feature_Set, group_var, Repeat_Number) %>%
  summarise(Repeat_Balanced_Accuracy = caret::confusionMatrix(data = CV_Predicted_Diagnosis,
                                                              reference = Diagnosis)$byClass[["Balanced Accuracy"]],
            Repeat_Accuracy = sum(CV_Predicted_Diagnosis==Diagnosis)/n()) %>%
  mutate(Aggregate_Type = "repeat_caret")

################################################################################
# Comparisons
################################################################################

# Q1: Does scikit-learn calculate balanced accuracy differently from caret?
plyr::rbind.fill(catch22_SVM_balanced_accuracy_by_fold_scikit,
                 catch22_SVM_balanced_accuracy_by_fold_caret) %>%
  dplyr::select(-index, -Comparison_Group, -Scaling_Type) %>%
  pivot_wider(id_cols = c(Fold, Repeat_Number, Analysis_Type, 
                          Univariate_Feature_Set, group_var),
              names_from = Fold_Accuracy_Method,
              values_from = Balanced_Accuracy) %>%
  filter(scikit_learn != caret)
# A1: nope, all values are the same. Woohoo!

# Q2: Does regular accuracy differ using caret when you aggregate by fold first versus go straight to repeat-wise?
plyr::rbind.fill(catch22_SVM_balanced_accuracy_by_fold_repeat_caret,
                 catch22_SVM_balanced_accuracy_by_repeat_caret) %>%
  pivot_wider(id_cols = c(Repeat_Number, Analysis_Type, 
                          Univariate_Feature_Set, group_var),
              names_from = Aggregate_Type,
              values_from = Repeat_Accuracy) %>%
  filter(fold_repeat_caret != repeat_caret) 
# Yep

# Q3: Does balanced accuracy differ using caret when you aggregate by fold first versus go straight to repeat-wise?
plyr::rbind.fill(catch22_SVM_balanced_accuracy_by_fold_repeat_caret,
                 catch22_SVM_balanced_accuracy_by_repeat_caret)%>%
  pivot_wider(id_cols = c(Repeat_Number, Analysis_Type, 
                          Univariate_Feature_Set, group_var),
              names_from = Aggregate_Type,
              values_from = Repeat_Balanced_Accuracy) %>%
  filter(fold_repeat_caret != repeat_caret) 
# Yep

# Compare accuracy with vs without folds
# Let's zoom in to e.g. right amygdala for schizophrenia
right_amyg_res <- fold_wise_sample_predictions %>%
  filter(group_var == "Right-Amygdala") %>%
  left_join(., fold_assignments) %>%
  group_by(Repeat_Number, Fold) %>%
  summarise(Num_Correct = sum(Diagnosis == CV_Predicted_Diagnosis),
            Num_Total = n())

# Pivot based on number of correct predictions
right_amyg_res %>%
  pivot_wider(id_cols = Fold, names_from = Repeat_Number, values_from = Num_Correct)

spaghetti_plot_acc <- plyr::rbind.fill(catch22_SVM_balanced_accuracy_by_fold_repeat_caret,
                                       catch22_SVM_balanced_accuracy_by_repeat_caret) %>%
  mutate(Aggregate_Type = ifelse(Aggregate_Type == "fold_repeat_caret", "Fold then repeat", "Overall repeat")) %>%
  dplyr::group_by(Analysis_Type, Univariate_Feature_Set, group_var, Aggregate_Type) %>%
  summarise(mean_acc = 100*mean(Repeat_Accuracy)) %>%
  mutate(line_var = paste0(Univariate_Feature_Set, "_", group_var)) %>%
  ggplot(data=., mapping=aes(x=Aggregate_Type, y=mean_acc)) +
  scale_x_discrete(expand=c(0.05,0.05)) +
  geom_line(aes(group = line_var), alpha=0.3) +
  xlab("Accuracy Aggregation Type") +
  ylab("Mean Accuracy\nAcross Resamples")

# Compare balanced accuracy with vs without folds
spaghetti_plot_balacc <- plyr::rbind.fill(catch22_SVM_balanced_accuracy_by_fold_repeat_caret,
                                         catch22_SVM_balanced_accuracy_by_repeat_caret) %>%
  mutate(Aggregate_Type = ifelse(Aggregate_Type == "fold_repeat_caret", "Fold then repeat", "Overall repeat")) %>%
  dplyr::group_by(Analysis_Type, Univariate_Feature_Set, group_var, Aggregate_Type) %>%
  summarise(mean_balacc = 100*mean(Repeat_Balanced_Accuracy)) %>%
  mutate(line_var = paste0(Univariate_Feature_Set, "_", group_var)) %>%
  ggplot(data=., mapping=aes(x=Aggregate_Type, y=mean_balacc)) +
  scale_x_discrete(expand=c(0.05,0.05)) +
  geom_line(aes(group = line_var), alpha=0.3) +
  xlab("Balanced Accuracy Aggregation Type") +
  ylab("Mean Balanced Accuracy\nAcross Resamples")

# Plot distribution of fold-wise minus direct repeat-wise accuracy
histogram_plot_acc <- plyr::rbind.fill(catch22_SVM_balanced_accuracy_by_fold_repeat_caret,
                                         catch22_SVM_balanced_accuracy_by_repeat_caret) %>%
  dplyr::group_by(Analysis_Type, Univariate_Feature_Set, group_var, Aggregate_Type) %>%
  summarise(mean_acc = 100*mean(Repeat_Accuracy)) %>%
  pivot_wider(id_cols = c(Analysis_Type, Univariate_Feature_Set, group_var),
              names_from=Aggregate_Type,
              values_from=mean_acc) %>%
  mutate(acc_diff = fold_repeat_caret - repeat_caret) %>%
  ggplot(data=., mapping=aes(x=acc_diff)) +
  scale_x_continuous(limits=c(-8,8), 
                     breaks=c(-8,-4,0,4,8)) +
  ylab("# Models") +
  xlab("Fold-wise vs. overall aggregate\n accuracy (\u0394%)") +
  geom_histogram(fill="gray70") +
  geom_vline(xintercept = 0)

# Plot distribution of fold-wise minus direct repeat-wise balanced accuracy
histogram_plot_balacc <- plyr::rbind.fill(catch22_SVM_balanced_accuracy_by_fold_repeat_caret,
                                          catch22_SVM_balanced_accuracy_by_repeat_caret) %>%
  dplyr::group_by(Analysis_Type, Univariate_Feature_Set, group_var, Aggregate_Type) %>%
  summarise(mean_balacc = 100*mean(Repeat_Balanced_Accuracy)) %>%
  pivot_wider(id_cols = c(Analysis_Type, Univariate_Feature_Set, group_var),
              names_from=Aggregate_Type,
              values_from=mean_balacc) %>%
  mutate(balacc_diff = fold_repeat_caret - repeat_caret) %>%
  ggplot(data=., mapping=aes(x=balacc_diff)) +
  scale_x_continuous(limits=c(-8,8), 
                     breaks=c(-8,-4,0,4,8)) +
  ylab("# Models") +
  xlab("Fold-wise vs. overall aggregate\nbalanced accuracy (\u0394%)") +
  geom_histogram(fill="gray70") +
  geom_vline(xintercept = 0)

(spaghetti_plot_acc + histogram_plot_acc) / (spaghetti_plot_balacc + histogram_plot_balacc) + 
  plot_layout(widths = c(1, 1.5))
