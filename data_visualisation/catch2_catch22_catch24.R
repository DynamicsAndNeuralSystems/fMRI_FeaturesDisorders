################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/catch2_misc/")
TAF::mkdir(plot_path)

# python_to_use <- "~/.conda/envs/pyspi/bin/python3"
python_to_use <- "/Users/abry4213/opt/anaconda3/envs/pyspi/bin/python3"
pairwise_feature_set <- "pyspi14"
univariate_feature_set <- "catch22"
UCLA_CNP_data_path <- "~/data/UCLA_CNP/"
ABIDE_ASD_data_path <- "~/data/ABIDE_ASD/"
data_path <- "~/data/TS_feature_manuscript/"
study_group_df <- data.frame(Study = c(rep("UCLA_CNP", 3), "ABIDE_ASD"),
                             Noise_Proc = c(rep("AROMA+2P+GMR",3), "FC1000"),
                             Comparison_Group = c("Schizophrenia", "ADHD", "Bipolar", "ASD"))
univariate_feature_sets <- c("catch22", "catch2", "catch24")

ABIDE_ASD_brain_region_info <- read.csv("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv")

reticulate::use_python(python_to_use)

library(reticulate)

# Import pyarrow.feather as pyarrow_feather
pyarrow_feather <- import("pyarrow.feather")

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

################################################################################
# Load libraries
################################################################################
library(feather)
library(tidyverse)
library(glue)
library(icesTAF)
library(cowplot)
library(ggseg)
library(ggsegHO)
library(knitr)
library(kableExtra)
library(RColorBrewer)
library(correctR)
theme_set(theme_cowplot())

# Load in univariate time-series feature info
TS_feature_info <- read.csv(glue("{github_dir}/data_visualisation/manuscript_figures/catch24_info.csv"))
# Load data
UCLA_CNP_metadata <- pyarrow_feather$read_feather(glue("{UCLA_CNP_data_path}/study_metadata/UCLA_CNP_sample_metadata.feather"))
ABIDE_ASD_metadata <- pyarrow_feather$read_feather(glue("{ABIDE_ASD_data_path}/study_metadata/ABIDE_ASD_sample_metadata.feather"))
metadata <- plyr::rbind.fill(UCLA_CNP_metadata, ABIDE_ASD_metadata)
univariate_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_robustsigmoid_scaler_balanced_accuracy_all_folds.feather"))

univariate_SVM_coefficients <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_robustsigmoid_scaler_fold_SVM_coefficients.feather"))

univariate_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_robustsigmoid_scaler_empirical_p_values.feather"))

univariate_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_null_balanced_accuracy_distributions.feather"))

univariate_fold_dist <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_robustsigmoid_scaler_fold_assignments.feather"))

# TO COMMENT OUT
SCZ_univariate_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{UCLA_CNP_data_path}/processed_data/UCLA_CNP_Schizophrenia_Univariate_catch22_robustsigmoid_scaler_SVM_balanced_accuracy.feather"))

# Aggregate the main results across all folds, independent of repeat
SCZ_univariate_balanced_accuracy <- SCZ_univariate_balanced_accuracy_all_folds %>%
  mutate(Study = "UCLA_CNP", Univariate_Feature_Set = "catch22") %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Analysis_Type, group_var) %>%
  summarise(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
            Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T))


SCZ_univariate_SVM_coefficients <- pyarrow_feather$read_feather(glue("{UCLA_CNP_data_path}/processed_data/UCLA_CNP_Schizophrenia_Univariate_catch22_robustsigmoid_scaler_SVM_fold_SVM_coefficients.feather"))

SCZ_univariate_null_distribution <- pyarrow_feather$read_feather(glue("{UCLA_CNP_data_path}/processed_data/UCLA_CNP_Schizophrenia_Univariate_catch22_robustsigmoid_scaler_SVM_null_balanced_accuracy_distributions.feather"))

univariate_split <- SCZ_univariate_balanced_accuracy %>%
  group_by(Study, Comparison_Group, Analysis_Type, Univariate_Feature_Set, group_var) %>%
  group_split()

SCZ_univariate_p_values <- calculate_empirical_p_values(main_balanced_accuracy_split = univariate_split,
                                                    null_balanced_accuracy = SCZ_univariate_null_distribution)

# Adjust p-values by group
SCZ_univariate_p_values <- SCZ_univariate_p_values %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Analysis_Type) %>%
  mutate(p_value_BH = p.adjust(p_value, method="BH"))
# TO COMMENT OUT

sample_sizes <- univariate_fold_dist %>%
  distinct(Study, Sample_ID) %>%
  left_join(., metadata) %>%
  group_by(Study, Diagnosis) %>%
  count()

# Brain Region-wise
univariate_balanced_accuracy_all_folds %>%
  filter(Analysis_Type == "Brain_Region") %>%
  group_by(Study, Comparison_Group, Analysis_Type, Scaling_Type, Univariate_Feature_Set, group_var) %>%
  summarise(BalAcc_Mean = mean(Balanced_Accuracy, na.rm=T),
            BalAcc_SD = sd(Balanced_Accuracy, na.rm=T)) %>%
  ggplot(data=., mapping=aes(x = Univariate_Feature_Set, y = BalAcc_Mean, 
                             group = group_var, color = Comparison_Group)) +
  scale_colour_brewer(palette = "Dark2") +
  geom_line(alpha=0.4) +
  facet_wrap(Comparison_Group ~ .) +
  ggtitle("Mean Balanced Accuracy for Brain Regions\nAcross CV Folds + Repeats") +
  ylab("Mean Balanced Accuracy") +
  xlab("Feature Set") +
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5))
ggsave(glue("{plot_path}/BalAcc_Brain_Region_Feature_Set.png"),
       width = 6, height = 5, units="in", dpi=300)

# Find which regions are better/worse/no difference in catch2 vs catch22
schizophrenia_region_data <- univariate_balanced_accuracy_all_folds %>%
  filter(Study == "UCLA_CNP", Comparison_Group == "Schizophrenia", 
         Univariate_Feature_Set %in% c("catch2", "catch22"),
         Analysis_Type == "Brain_Region") %>%
  group_by(Study, group_var, Univariate_Feature_Set) %>%
  summarise(Mean_BalAcc = mean(Balanced_Accuracy, na.rm=T)) %>%
  pivot_wider(id_cols = c(Study, group_var), 
              names_from = Univariate_Feature_Set, 
              values_from = Mean_BalAcc)

schizophrenia_region_data %>%
  mutate(Status = case_when(catch2 - catch22 > 0.05 ~ "catch2 better",
                            catch22 - catch2 > 0.05 ~ "catch22 better",
                            T ~ "neither")) %>%
  filter(Status != "neither") %>%
  pivot_longer(cols = c(catch2, catch22), names_to="Feature_Set", values_to="values") %>%
  ggplot(data = ., mapping=aes(x = Feature_Set, y=values, group=group_var, color = Status)) +
  scale_colour_brewer(palette = "Dark2") +
  ggtitle("Feature Set Performance by Brain Region\nin UCLA CNP Schizophrenia Dataset") +
  ylab("Mean Balanced Accuracy") +
  xlab("Feature Set") +
  facet_grid(. ~ Status) +
  geom_line() +
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5, size=12))
ggsave(glue("{plot_path}/BalAcc_Brain_Region_Comparison_5percent.png"),
       width = 6, height = 4, units="in", dpi=300)

# TS Feature
univariate_balanced_accuracy_all_folds %>%
  filter(Analysis_Type == "TS_Feature") %>%
  group_by(Study, Comparison_Group, Analysis_Type, Scaling_Type, Univariate_Feature_Set, group_var) %>%
  summarise(BalAcc_Mean = mean(Balanced_Accuracy, na.rm=T),
            BalAcc_SD = sd(Balanced_Accuracy, na.rm=T)) %>%
  ggplot(data=., mapping=aes(x = Univariate_Feature_Set, y = BalAcc_Mean,
                             fill = Comparison_Group)) +
  geom_violin() +
  facet_wrap(Comparison_Group ~ .) +
  ggtitle("Mean Balanced Accuracy for TS Features\nAcross CV Folds + Repeats") +
  scale_colour_brewer(palette = "Dark2") +
  ylab("Mean Balanced Accuracy") +
  xlab("Feature Set") +
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5))
ggsave(glue("{plot_path}/BalAcc_TSFeature_Feature_Set.png"),
       width = 6, height = 5, units="in", dpi=300)

# Schizophrenia only -- TS Feature
schizophrenia_feature_data <- univariate_balanced_accuracy_all_folds %>%
  filter(Study == "UCLA_CNP", Comparison_Group == "Schizophrenia", 
         Univariate_Feature_Set %in% c("catch2", "catch22"),
         Analysis_Type == "TS_Feature") %>%
  group_by(Study, group_var, Univariate_Feature_Set) %>%
  summarise(Mean_BalAcc = mean(Balanced_Accuracy, na.rm=T))

schizophrenia_feature_data %>%
  ggplot(data=., mapping=aes(x = Mean_BalAcc, color = Univariate_Feature_Set)) +
  geom_vline(aes(xintercept = Mean_BalAcc, color = Univariate_Feature_Set), size=1.5) +
  ggtitle("Mean Balanced Accuracy for TS Features\nAcross CV Folds + Repeats") +
  labs(color = "Feature Set") +
  scale_colour_brewer(palette = "Dark2") +
  xlab("Mean Balanced Accuracy") +
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        plot.title = element_text(hjust=0.5,margin=margin(0,0,20,0)))
ggsave(glue("{plot_path}/BalAcc_TSFeature_Feature_Set_Schizophrenia.png"),
       width = 6, height = 3, units="in", dpi=300)


# Combo-wise
univariate_balanced_accuracy_all_folds %>%
  filter(Analysis_Type == "Combo") %>%
  group_by(Study, Comparison_Group, Analysis_Type, Scaling_Type, Univariate_Feature_Set, group_var) %>%
  summarise(BalAcc_Mean = mean(Balanced_Accuracy, na.rm=T),
            BalAcc_SD = sd(Balanced_Accuracy, na.rm=T)) %>%
  ggplot(data=., mapping=aes(x = Univariate_Feature_Set, y = BalAcc_Mean, 
                             group = group_var, color = Comparison_Group)) +
  geom_line() +
  geom_ribbon(aes(ymin = BalAcc_Mean - BalAcc_SD, 
                  ymax = BalAcc_Mean + BalAcc_SD,
                  fill = Comparison_Group), alpha=0.4, color=NA) +
  facet_wrap(Comparison_Group ~ .) +
  ggtitle("Mean Balanced Accuracy for Univariate Combo\nAcross CV Folds + Repeats") +
  scale_fill_brewer(palette = "Dark2") +
  scale_colour_brewer(palette = "Dark2") +
  ylab("Mean Balanced Accuracy") +
  xlab("Feature Set") +
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5))
ggsave(glue("{plot_path}/BalAcc_Combo_Feature_Set.png"),
       width = 6, height = 5, units="in", dpi=300)

# Schizophrenia only -- Combo
schizophrenia_combo_data <- univariate_balanced_accuracy_all_folds %>%
  filter(Study == "UCLA_CNP", Comparison_Group == "Schizophrenia", 
         Univariate_Feature_Set %in% c("catch2", "catch22"),
         Analysis_Type == "Combo") %>%
  group_by(Study, group_var, Univariate_Feature_Set, Repeat_Number) %>%
  summarise(Mean_BalAcc = mean(Balanced_Accuracy, na.rm=T)) 

schizophrenia_combo_data %>%
  ggplot(data=., mapping=aes(x = Univariate_Feature_Set, 
                             y = Mean_BalAcc,
                             fill = Univariate_Feature_Set)) +
  ggtitle("Mean Balanced Accuracy for Univariate Combo\nAcross CV Folds + Repeats") +
  scale_fill_brewer(palette = "Dark2") +
  labs(fill = "Feature Set") +
  geom_violin() +
  geom_boxplot(fill=NA, width=0.1) +
  ylab("Mean Balanced Accuracy") +
  xlab("Feature Set") +
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5))
ggsave(glue("{plot_path}/BalAcc_Combo_Feature_Set_Schizophrenia.png"),
       width = 6, height = 5, units="in", dpi=300)

# Use correctR to test for difference across repeated k-fold CV
univariate_balanced_accuracy_all_folds %>%
  filter(Study == "UCLA_CNP", Comparison_Group == "Schizophrenia", 
         Univariate_Feature_Set %in% c("catch2", "catch22"),
         Analysis_Type == "Combo") %>%
  dplyr::rename("model" = "Univariate_Feature_Set",
                k = "Fold",
                r = "Repeat_Number",
                values = "Balanced_Accuracy") %>%
  repkfold_ttest(data = .,
                 n1 = ceiling(0.9*num_samples),
                 n2 = floor(0.1*num_samples),
                 k = 10,
                 r = 10)
