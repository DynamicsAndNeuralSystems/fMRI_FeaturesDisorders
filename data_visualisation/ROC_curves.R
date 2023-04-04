################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/univariate_results/")
TAF::mkdir(plot_path)

# python_to_use <- "~/.conda/envs/pyspi/bin/python3"
python_to_use <- "/Users/abry4213/opt/anaconda3/envs/pyspi/bin/python3"
pairwise_feature_set <- "pyspi14"
univariate_feature_set <- "catch24"
data_path <- "~/data/TS_feature_manuscript"
study_group_df <- data.frame(Study = c(rep("UCLA_CNP", 3), "ABIDE_ASD"),
                             Noise_Proc = c(rep("AROMA+2P+GMR",3), "FC1000"),
                             Comparison_Group = c("Schizophrenia", "Bipolar", "ADHD", "ASD"),
                             Group_Nickname = c("SCZ", "BPD", "ADHD", "ASD"))


UCLA_CNP_brain_region_info <- read.csv("~/data/UCLA_CNP/study_metadata/UCLA_CNP_Brain_Region_info.csv")
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
library(MESS)
library(broom)
library(colorspace)
theme_set(theme_cowplot())

# Source visualisation script
source(glue("{github_dir}/data_visualisation/Manuscript_Draft_Visualisations_Helper.R"))

# Load in univariate time-series feature info
TS_feature_info <- read.csv(glue("{github_dir}/data_visualisation/catch24_info.csv"))

# Load data
univariate_ROC <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_mixedsigmoid_scaler_ROC.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)
pairwise_ROC <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_mixedsigmoid_scaler_ROC.feather")) 
combined_univariate_pairwise_ROC <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combined_univariate_pairwise_mixedsigmoid_scaler_ROC.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)

temp <- pyarrow_feather$read_feather("~/Desktop/ABIDE_ASD_GSR_catch2_filtered.feather")


################################################################################
# ROC Curves + AUC
################################################################################

# Univariate combo
univariate_ROC %>%
  filter(group_var == "Combo") %>%
  filter(Comparison_Group=="BPD") %>%
  # group_by(Comparison_Group) %>%
  summarise(AUC = auc(fpr, tpr))



<- AUC(univariate_ROC %>% filter(group_var=="Combo", Comparison_Group=="Schizophrenia") %>% pull(fpr),
                            univariate_ROC %>% filter(group_var=="Combo", Comparison_Group=="Schizophrenia") %>% pull(tpr))

univariate_ROC %>%
  filter(group_var=="Combo", Comparison_Group=="Schizophrenia") %>%
  ggplot(data=., mapping=aes(x=fpr,y=tpr)) +
  geom_abline(slope=1, linetype=2) +
  geom_smooth(se=T) +
  xlab("False positive rate (FPR)") +
  ylab("True positive rate (TPR)") +
  coord_equal() +
  annotate(geom="text", label=glue("AUC: {round(univariate_combo_AUC, 3)}"),
           x = 0.7, y=0.25, size=8)