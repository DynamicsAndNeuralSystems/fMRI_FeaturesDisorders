################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/univariate_results/")
TAF::mkdir(plot_path)

# python_to_use <- "~/.conda/envs/pyspi/bin/python3"
python_to_use <- "/Users/abry4213/anaconda3/envs/pyspi/bin/python3"
pairwise_feature_set <- "pyspi14"
univariate_feature_set <- "catch24"
data_path <- "~/data/TS_feature_manuscript"
study_group_df <- data.frame(Study = c(rep("UCLA_CNP", 3), "ABIDE_ASD"),
                             Noise_Proc = c(rep("AROMA+2P+GMR",3), "FC1000"),
                             Comparison_Group = c("Schizophrenia", "Bipolar", "ADHD", "ASD"),
                             Group_Nickname = c("SCZ", "BPD", "ADHD", "ASD"))
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
library(ggsegDefaultExtra)
library(knitr)
library(kableExtra)
library(patchwork)
library(broom)
library(colorspace)
library(see)
library(ggridges)
library(splitstackshape)
library(DescTools)
library(metR)
library(scales)
library(correctR)
theme_set(theme_cowplot())

# Source visualisation script
source(glue("{github_dir}/data_visualisation/Manuscript_Draft_Visualisations_Helper.R"))

UCLA_CNP_brain_region_info <- read.csv("~/data/UCLA_CNP/study_metadata/UCLA_CNP_Brain_Region_info.csv")
ABIDE_ASD_brain_region_info <- read.table("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_Harvard_Oxford_cort_prob_2mm_ROI_lookup.txt", sep=";", header = T) %>%
  mutate(Brain_Region = ifelse(Index==45, "Heschl's Gyrus (includes H1 and H2)", Brain_Region))

# Load in univariate time-series feature info
TS_feature_info <- read.csv(glue("{github_dir}/data_visualisation/catch24_info.csv"))
# Load data
univariate_balanced_accuracy_AUC_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_mixedsigmoid_scaler_balanced_accuracy_AUC_all_folds.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)
univariate_balanced_accuracy <- univariate_balanced_accuracy_AUC_all_folds %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Analysis_Type, group_var, Repeat_Number, kernel) %>%
  reframe(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
          ROC_AUC_Across_Folds = mean(ROC_AUC, na.rm=T)) %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Analysis_Type, group_var, kernel) %>%
  reframe(Balanced_Accuracy_Across_Repeats = mean(Balanced_Accuracy_Across_Folds, na.rm=T),
          Balanced_Accuracy_Across_Repeats_SD = sd(Balanced_Accuracy_Across_Folds, na.rm=T),
          ROC_AUC_Across_Repeats = mean(ROC_AUC_Across_Folds, na.rm=T),
          ROC_AUC_Across_Repeats_SD = sd(ROC_AUC_Across_Folds, na.rm=T))


# Plot linear vs. RBF kernel for each clinical group
univariate_balanced_accuracy %>%
  mutate(Analysis_Type = case_when(Analysis_Type == "Univariate_Brain_Region" ~ "Brain Regions",
                                   Analysis_Type == "Univariate_TS_Feature" ~ "catch24 Features",
                                   T ~ "Combination")) %>%
  ggplot(data=., mapping=aes(x=kernel, y=100*Balanced_Accuracy_Across_Repeats, 
                             color=Analysis_Type, group=group_var)) +
  geom_line(alpha=0.5) +
  ylab("Mean Balanced Accuracy (%)") +
  xlab("SVM Kernel Type") +
  facet_grid(Comparison_Group ~ Analysis_Type, scales="free", space="free", switch = "y") +
  theme(legend.position = "none",
        strip.placement = "outside")
ggsave(glue("{plot_path}/univariate_linear_vs_RBF_kernel.pdf"),
       width=6, height=4, units="in", dpi=300)