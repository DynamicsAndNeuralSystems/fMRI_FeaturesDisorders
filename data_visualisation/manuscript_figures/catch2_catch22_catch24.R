################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/FigureS2/")
TAF::mkdir(plot_path)

# python_to_use <- "~/.conda/envs/pyspi/bin/python3"
python_to_use <- "/Users/abry4213/opt/anaconda3/envs/pyspi/bin/python3"
pairwise_feature_set <- "pyspi14"
univariate_feature_set <- "catch22"
data_path <- "~/data/TS_feature_manuscript"
study_group_df <- data.frame(Study = c(rep("UCLA_CNP", 3), "ABIDE_ASD"),
                             Noise_Proc = c(rep("AROMA+2P+GMR",3), "FC1000"),
                             Comparison_Group = c("Schizophrenia", "ADHD", "Bipolar", "ASD"))
univariate_feature_sets <- c("catch22", "catch2", "catch24")

ABIDE_ASD_brain_region_info <- read.csv("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv")

reticulate::use_python(python_to_use)

library(reticulate)

# Import pyarrow.feather as pyarrow_feather
pyarrow_feather <- import("pyarrow.feather")

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
theme_set(theme_cowplot())

# Load in univariate time-series feature info
TS_feature_info <- read.csv(glue("{github_dir}/data_visualisation/manuscript_figures/catch24_info.csv"))
# Load data
univariate_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_robustsigmoid_scaler_balanced_accuracy_all_folds.feather"))

univariate_SVM_coefficients <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_robustsigmoid_scaler_balanced_accuracy_all_folds.feather"))

univariate_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_robustsigmoid_scaler_empirical_p_values.feather"))

univariate_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_null_balanced_accuracy_distributions.feather"))

# Brain Region-wise
univariate_balanced_accuracy_all_folds %>%
  filter(Analysis_Type == "Brain_Region") %>%
  group_by(Study, Comparison_Group, Analysis_Type, Scaling_Type, Univariate_Feature_Set, group_var) %>%
  summarise(BalAcc_Mean = mean(Balanced_Accuracy, na.rm=T),
            BalAcc_SD = sd(Balanced_Accuracy, na.rm=T)) %>%
  ggplot(data=., mapping=aes(x = Univariate_Feature_Set, y = BalAcc_Mean, 
                             group = group_var, color = Comparison_Group)) +
  geom_line(alpha=0.4) +
  facet_wrap(Comparison_Group ~ .) +
  ggtitle("Mean Balanced Accuracy for Brain Regions\nAcross CV Folds + Repeats") +
  ylab("Mean Balanced Accuracy") +
  xlab("Feature Set") +
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5))
ggsave(glue("{plot_path}/BalAcc_Brain_Region_Feature_Set.png"),
       width = 6, height = 5, units="in", dpi=300)

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
  ylab("Mean Balanced Accuracy") +
  xlab("Feature Set") +
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5))
ggsave(glue("{plot_path}/BalAcc_TSFeature_Feature_Set.png"),
       width = 6, height = 5, units="in", dpi=300)

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
  ylab("Mean Balanced Accuracy") +
  xlab("Feature Set") +
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5))
ggsave(glue("{plot_path}/BalAcc_Combo_Feature_Set.png"),
       width = 6, height = 5, units="in", dpi=300)