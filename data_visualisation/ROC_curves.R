################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/ROC_AUC/")
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
univariate_ROC <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_mixedsigmoid_scaler_ROC_TPR_FPR.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)
pairwise_ROC <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_mixedsigmoid_scaler_ROC_TPR_FPR.feather")) 
combined_univariate_pairwise_ROC <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combined_univariate_pairwise_mixedsigmoid_scaler_ROC_TPR_FPR.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)

UCLA_CNP_brain_region_info <- read.csv("~/data/UCLA_CNP/study_metadata/UCLA_CNP_Brain_Region_info.csv")
ABIDE_ASD_brain_region_info <- read.table("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_Harvard_Oxford_cort_prob_2mm_ROI_lookup.txt", sep=";", header = T) %>%
  mutate(Brain_Region = ifelse(Index==45, "Heschl's Gyrus (includes H1 and H2)", Brain_Region))

################################################################################
# ROC Curves + AUC
################################################################################

# Univariate top-performing region per group
top_region_AUC <- univariate_ROC %>%
  filter(Analysis_Type == "Brain_Region") %>%
  group_by(Comparison_Group, Study, group_var, Repeat_Number) %>%
  summarise(AUC = auc(fpr, tpr)) %>%
  group_by(Comparison_Group, Study, group_var) %>%
  summarise(Mean_AUC = mean(AUC)) %>%
  group_by(Comparison_Group) %>%
  filter(Mean_AUC == max(Mean_AUC)) %>%
  mutate(Group_Nickname =  case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                     Comparison_Group == "Bipolar" ~ "BPD",
                                     T ~ Comparison_Group))
top_region_AUC

univariate_ROC %>%
  semi_join(top_region_AUC) %>%
  ggplot(data=.) +
  geom_abline(slope=1, linetype=2) +
  geom_smooth(se=T, aes(color=Comparison_Group, x=fpr,y=tpr)) +
  xlab("False positive rate (FPR)") +
  ylab("True positive rate (TPR)") +
  coord_equal() +
  geom_text(data = top_region_AUC,
            aes(label=glue("{Group_Nickname}: {round(Mean_AUC, 2)}"),
                color=Comparison_Group),
            x = 1, y=c(0.1, 0.2, 0.3, 0.4), 
            size=4.5, hjust=1) +
  theme(legend.position = "none") +
  scale_color_manual(values=c("Control" = "#5BB67B", 
                              "Schizophrenia" = "#573DC7", 
                              "Bipolar" = "#D5492A", 
                              "ADHD" = "#0F9EA9", 
                              "ASD" = "#C47B2F"))
ggsave(glue("{plot_path}/univariate_top_region_ROC.png"),
       width=3, height=3, units="in", dpi=300)

# Univariate top-performing feature per group
top_feature_AUC <- univariate_ROC %>%
  filter(Analysis_Type == "TS_Feature") %>%
  group_by(Comparison_Group, Study, group_var, Repeat_Number) %>%
  summarise(AUC = auc(fpr, tpr)) %>%
  group_by(Comparison_Group, Study, group_var) %>%
  summarise(Mean_AUC = mean(AUC)) %>%
  group_by(Comparison_Group) %>%
  filter(Mean_AUC == max(Mean_AUC)) %>%
  mutate(Group_Nickname =  case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                     Comparison_Group == "Bipolar" ~ "BPD",
                                     T ~ Comparison_Group))
top_feature_AUC

univariate_ROC %>%
  semi_join(top_feature_AUC) %>%
  ggplot(data=.) +
  geom_abline(slope=1, linetype=2) +
  geom_smooth(se=T, aes(color=Comparison_Group, x=fpr,y=tpr)) +
  xlab("False positive rate (FPR)") +
  ylab("True positive rate (TPR)") +
  coord_equal() +
  geom_text(data = top_feature_AUC,
            aes(label=glue("{Group_Nickname}: {round(Mean_AUC, 2)}"),
                color=Comparison_Group),
            x = 1, y=c(0.1, 0.2, 0.3, 0.4), 
            size=4.5, hjust=1) +
  theme(legend.position = "none") +
  scale_color_manual(values=c("Control" = "#5BB67B", 
                              "Schizophrenia" = "#573DC7", 
                              "Bipolar" = "#D5492A", 
                              "ADHD" = "#0F9EA9", 
                              "ASD" = "#C47B2F"))
ggsave(glue("{plot_path}/univariate_top_feature_ROC.png"),
       width=3, height=3, units="in", dpi=300)



# Univariate combo
mean_AUC <- univariate_ROC %>%
  filter(group_var == "Combo") %>%
  group_by(Comparison_Group, Study, Repeat_Number) %>%
  summarise(AUC = auc(fpr, tpr)) %>%
  group_by(Comparison_Group, Study) %>%
  summarise(Mean_AUC = mean(AUC)) %>%
  mutate(Group_Nickname =  case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                       Comparison_Group == "Bipolar" ~ "BPD",
                                       T ~ Comparison_Group))

univariate_ROC %>%
  mutate(Group_Nickname =  case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                     Comparison_Group == "Bipolar" ~ "BPD",
                                     T ~ Comparison_Group)) %>%
  filter(group_var=="Combo") %>%
  ggplot(data=.) +
  geom_abline(slope=1, linetype=2) +
  geom_smooth(se=T, aes(color=Comparison_Group, x=fpr,y=tpr)) +
  xlab("False positive rate (FPR)") +
  ylab("True positive rate (TPR)") +
  coord_equal() +
  geom_text(data = mean_AUC,
            aes(label=glue("{Group_Nickname}: {round(Mean_AUC, 2)}"),
                color=Comparison_Group),
            x = 1, y=c(0.1, 0.2, 0.3, 0.4), 
            size=4.5, hjust=1) +
  theme(legend.position = "none") +
  scale_color_manual(values=c("Control" = "#5BB67B", 
                              "Schizophrenia" = "#573DC7", 
                              "Bipolar" = "#D5492A", 
                              "ADHD" = "#0F9EA9", 
                              "ASD" = "#C47B2F"))
ggsave(glue("{plot_path}/univariate_catch24_combo_ROC.png"),
       width=3, height=3, units="in", dpi=300)


# Pairwise top-performing SPI with univariate info per group
top_feature_AUC <- combined_univariate_pairwise_ROC %>%
  filter(Analysis_Type == "SPI_Univariate_Combo") %>%
  group_by(Comparison_Group, Study, group_var, Repeat_Number) %>%
  summarise(AUC = auc(fpr, tpr)) %>%
  group_by(Comparison_Group, Study, group_var) %>%
  summarise(Mean_AUC = mean(AUC)) %>%
  group_by(Comparison_Group) %>%
  filter(Mean_AUC == max(Mean_AUC)) %>%
  mutate(Group_Nickname =  case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                     Comparison_Group == "Bipolar" ~ "BPD",
                                     T ~ Comparison_Group))
top_feature_AUC

combined_univariate_pairwise_ROC %>%
  semi_join(top_feature_AUC) %>%
  ggplot(data=.) +
  geom_abline(slope=1, linetype=2) +
  geom_smooth(se=T, aes(color=Comparison_Group, x=fpr,y=tpr)) +
  xlab("False positive rate (FPR)") +
  ylab("True positive rate (TPR)") +
  coord_equal() +
  geom_text(data = top_feature_AUC,
            aes(label=glue("{Group_Nickname}: {round(Mean_AUC, 2)}"),
                color=Comparison_Group),
            x = 1, y=c(0.1, 0.2, 0.3, 0.4), 
            size=4.5, hjust=1) +
  theme(legend.position = "none") +
  scale_color_manual(values=c("Control" = "#5BB67B", 
                              "Schizophrenia" = "#573DC7", 
                              "Bipolar" = "#D5492A", 
                              "ADHD" = "#0F9EA9", 
                              "ASD" = "#C47B2F"))
ggsave(glue("{plot_path}/univariate_top_SPI_univariate_pairwise_combo_ROC.png"),
       width=3, height=3, units="in", dpi=300)


# Pairwise all SPIs
mean_AUC <- pairwise_all_SPIs_ROC %>%
  filter(group_var == "Combo") %>%
  group_by(Comparison_Group, Study, Repeat_Number) %>%
  summarise(AUC = auc(fpr, tpr)) %>%
  group_by(Comparison_Group, Study) %>%
  summarise(Mean_AUC = mean(AUC)) %>%
  mutate(Group_Nickname =  case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                     Comparison_Group == "Bipolar" ~ "BPD",
                                     T ~ Comparison_Group))

pairwise_all_SPIs_ROC %>%
  mutate(Group_Nickname =  case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                     Comparison_Group == "Bipolar" ~ "BPD",
                                     T ~ Comparison_Group)) %>%
  filter(group_var == "Combo") %>%
  ggplot(data=.) +
  geom_abline(slope=1, linetype=2) +
  geom_smooth(se=T, aes(color=Comparison_Group, x=fpr,y=tpr)) +
  xlab("False positive rate (FPR)") +
  ylab("True positive rate (TPR)") +
  coord_equal() +
  geom_text(data = mean_AUC,
            aes(label=glue("{Group_Nickname}: {round(Mean_AUC, 2)}"),
                color=Comparison_Group),
            x = 1, y=c(0.1, 0.2, 0.3, 0.4), 
            size=4.5, hjust=1) +
  theme(legend.position = "none") +
  scale_color_manual(values=c("Control" = "#5BB67B", 
                              "Schizophrenia" = "#573DC7", 
                              "Bipolar" = "#D5492A", 
                              "ADHD" = "#0F9EA9", 
                              "ASD" = "#C47B2F"))
ggsave(glue("{plot_path}/pairwise_all_SPIs_combo_ROC.png"),
       width=3, height=3, units="in", dpi=300)