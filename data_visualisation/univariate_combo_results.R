################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/Library/CloudStorage/OneDrive-TheUniversityofSydney(Students)/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/univariate_results/")
TAF::mkdir(plot_path)

# python_to_use <- "~/.conda/envs/pyspi/bin/python3"
python_to_use <- "/Users/abry4213/anaconda3/envs/pyspi/bin/python3"
pairwise_feature_set <- "pyspi14"
univariate_feature_set <- "catch25"
SVM_kernel <- "Linear"
data_path <- "~/data/TS_feature_manuscript"
study_group_df <- data.frame(Study = c(rep("UCLA_CNP", 3), "ABIDE_ASD"),
                             Noise_Proc = c(rep("AROMA+2P+GMR",3), "FC1000"),
                             Comparison_Group = c("Schizophrenia", "Bipolar", "ADHD", "ASD"),
                             Group_Nickname = c("SCZ", "BP", "ADHD", "ASD"))
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
library(patchwork)
library(broom)
library(colorspace)
library(scales)
library(correctR)
library(splitstackshape)
library(FactoMineR)
library("factoextra")
library(ggpubr)

theme_set(theme_cowplot())

# Source visualisation script
source(glue("{github_dir}/data_visualisation/Manuscript_Draft_Visualisations_Helper.R"))

UCLA_CNP_brain_region_info <- read.csv("~/data/UCLA_CNP/study_metadata/UCLA_CNP_Brain_Region_info.csv")
ABIDE_ASD_brain_region_info <- read.table("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_Harvard_Oxford_cort_prob_2mm_ROI_lookup.txt", sep=";", header = T) %>%
  mutate(Brain_Region = ifelse(Index==45, "Heschl's Gyrus (includes H1 and H2)", Brain_Region))

# Load in univariate time-series feature info
TS_feature_info <- read.csv(glue("{github_dir}/data_visualisation/catch25_info.csv"))

# Load univariate classification results across all folds
univariate_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_balanced_accuracy_all_folds.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set, kernel==SVM_kernel)
# Compute mean + SD performance across all folds
univariate_balanced_accuracy <- univariate_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Analysis_Type, group_var, kernel) %>%
  reframe(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
          Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T))

# Load p-values
univariate_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_empirical_p_values.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set) %>%
  dplyr::select(-Balanced_Accuracy_Across_Folds) %>%
  left_join(., univariate_balanced_accuracy)

# Load null distribution
univariate_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_null_balanced_accuracy_distributions.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)

# Load study metadata
UCLA_CNP_metadata <- pyarrow_feather$read_feather("~/data/UCLA_CNP/study_metadata/UCLA_CNP_sample_metadata.feather") %>%
  mutate(Study="UCLA_CNP")
ABIDE_ASD_metadata <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_sample_metadata.feather")  %>%
  mutate(Study="ABIDE_ASD")

################################################################################
# Bowtie plot comparing each brain region and feature to combo-wise performance
num_comparisons <- univariate_p_values %>%
  filter(Analysis_Type != "Univariate_Combo") %>%
  group_by(Study, Comparison_Group) %>%
  count()

univariate_p_values %>%
  filter(Analysis_Type == "Univariate_Combo") %>% 
  left_join(., num_comparisons) %>%
  expandRows("n") %>%
  group_by(Study, Comparison_Group) %>%
  mutate(group_ID = paste0(Comparison_Group, "_", row_number())) %>%
  plyr::rbind.fill(., univariate_p_values %>%
                     filter(Analysis_Type != "Univariate_Combo") %>%
                     group_by(Study, Comparison_Group) %>%
                     mutate(group_ID = paste0(Comparison_Group, "_", row_number()))) %>%
  rowwise() %>%
  mutate(Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                      Comparison_Group == "Bipolar" ~ "BP",
                                      T ~ Comparison_Group),
         Analysis_Label = case_when(Analysis_Type == "Univariate_Combo" ~ "All\nRegions \u00D7\nFeatures",
                                    Analysis_Type == "Univariate_TS_Feature" ~ "Individual\ncatch25\nFeatures",
                                    T ~ "Individual\nBrain\nRegions")) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BP", "ADHD", "ASD")),
         Analysis_Label = factor(Analysis_Label, levels = c("Individual\ncatch25\nFeatures",
                                                            "All\nRegions \u00D7\nFeatures",
                                                            "Individual\nBrain\nRegions"))) %>%
  group_by(Study, Comparison_Group, group_ID) %>%
  mutate(Comparison_Sig = paste0(Comparison_Group, "_", p_value_Bonferroni[Analysis_Type != "Univariate_Combo"]<0.05)) %>% 
  ggplot(data=., mapping=aes(x=Analysis_Label, y=100*Balanced_Accuracy_Across_Folds, 
                             group=group_ID, color=Comparison_Sig)) +
  geom_line(alpha=0.5) +
  geom_point(aes(size=Analysis_Type)) +
  ylab("Mean Balanced Accuracy Across Repeats (%)") +
  xlab("Classifier Inputs") +
  geom_hline(yintercept = 50, linetype=2) +
  scale_x_discrete(expand=c(0.05,0.05,0.05,0.05)) +
  scale_size_manual(values=c("Univariate_Combo"=4,
                             "Univariate_TS_Feature"=1.5,
                             "Univariate_Brain_Region" = 1.5)) +
  scale_color_manual(values = c("SCZ_TRUE"="#573DC7", 
                                "BP_TRUE"="#D5492A", 
                                "ADHD_TRUE"="#0F9EA9",
                                "ASD_TRUE"="#C47B2F",
                                "SCZ_FALSE"="gray70", 
                                "BP_FALSE"="gray70", 
                                "ADHD_FALSE"="gray70",
                                "ASD_FALSE"="gray70")) +
  facet_grid(. ~ Comparison_Group, scales="free_x", space="free") +
  theme(legend.position = "none",
        plot.margin = margin(1,10,1,1),
        strip.background = element_blank(),
        strip.text = element_text(face="bold")) + 
  coord_flip()
ggsave(glue("{plot_path}/univariate_bowtie_balanced_accuracy.svg"),
       width=9, height=2.5, units="in", dpi=300)


################################################################################
# How do samples compare in PC space?

# UCLA CNP first
data_for_PCA_UCLA_CNP <- UCLA_CNP_catch25 %>%
  left_join(., TS_feature_info, by=c("names"="feature_name")) %>%
  mutate(unique_ID = paste0(Figure_name, "__", Brain_Region), .keep="unused") %>%
  dplyr::select(unique_ID, Sample_ID, Diagnosis, values) %>%
  pivot_wider(id_cols=c(Sample_ID, Diagnosis), names_from=unique_ID, values_from=values) %>%
  mutate(Diagnosis=factor(Diagnosis, levels=c("Control",
                                              "Schizophrenia",
                                              "Bipolar",
                                              "ADHD")))

UCLA_CNP_pca_res <- PCA(select(data_for_PCA_UCLA_CNP, c(-Sample_ID, -Diagnosis)), graph = FALSE, scale.unit = TRUE)
UCLA_CNP_pca_scores <- as.data.frame(UCLA_CNP_pca_res$ind$coord) %>%
  mutate(Sample_ID = data_for_PCA_UCLA_CNP$Sample_ID,
         Diagnosis = data_for_PCA_UCLA_CNP$Diagnosis)

# How do ABIDE samples compare in PC space?
data_for_PCA_ABIDE <- ABIDE_ASD_catch25 %>%
  left_join(., TS_feature_info, by=c("names"="feature_name")) %>%
  mutate(unique_ID = paste0(Figure_name, "__", Brain_Region), .keep="unused") %>%
  dplyr::select(unique_ID, Sample_ID, Diagnosis, values) %>%
  pivot_wider(id_cols=c(Sample_ID, Diagnosis), names_from=unique_ID, values_from=values) %>%
  mutate(Diagnosis=factor(Diagnosis, levels=c("Control",
                                              "ASD")))

ABIDE_pca_res <- PCA(select(data_for_PCA_ABIDE, c(-Sample_ID, -Diagnosis)), graph = FALSE, scale.unit = TRUE)
ABIDE_pca_scores <- as.data.frame(ABIDE_pca_res$ind$coord) %>%
  mutate(Sample_ID = data_for_PCA_ABIDE$Sample_ID,
         Diagnosis = data_for_PCA_ABIDE$Diagnosis)

UCLA_CNP_biplot <- UCLA_CNP_pca_scores %>%
  mutate(Diagnosis = case_when(Diagnosis == "Schizophrenia" ~ "SCZ",
                               Diagnosis == "Bipolar" ~ "BP",
                               T ~ Diagnosis)) %>%
  ggplot(data=., mapping=aes(x=`Dim.1`, y=`Dim.2`, color=Diagnosis)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype=2) +
  geom_hline(yintercept = 0, linetype=2) +
  stat_chull(aes(color = Diagnosis, fill = Diagnosis), 
             alpha = 0.1, geom = "polygon")+
  ylab("PC2") +
  xlab("PC1") +
  scale_fill_manual(values = c("SCZ"="#573DC7", 
                                "BP"="#D5492A", 
                                "ADHD"="#0F9EA9",
                                "Control" = "#5BB67B")) +
  scale_color_manual(values = c("SCZ"="#573DC7", 
                                "BP"="#D5492A", 
                                "ADHD"="#0F9EA9",
                                "Control" = "#5BB67B")) +
  theme(legend.position = "bottom")

ABIDE_biplot <- ABIDE_pca_scores %>%
  ggplot(data=., mapping=aes(x=`Dim.1`, y=`Dim.2`, color=Diagnosis)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype=2) +
  geom_hline(yintercept = 0, linetype=2) +
  stat_chull(aes(color = Diagnosis, fill = Diagnosis), 
             alpha = 0.1, geom = "polygon")+
  ylab("PC2") +
  xlab("PC1") +
  scale_fill_manual(values = c("ASD" = "#C47B2F",
                                "Control" = "#5BB67B")) +
  scale_color_manual(values = c("ASD" = "#C47B2F",
                                "Control" = "#5BB67B")) +
  theme(legend.position = "bottom")

# Combine biplots
UCLA_CNP_biplot + ABIDE_biplot
ggsave(glue("{plot_path}/Univariate_combo_PCA_biplots.svg"),
       width=8, height=4.5, units="in", dpi=300)

################################################################################
# Compare performance with L1-regularized SVM balanced accuracy implemented with LinearSVC()
SVM_L1_balanced_accuracy_by_folds <- pyarrow_feather$read_feather(glue("{data_path}/SVM_L1_Regularized_Balanced_Accuracy.feather"))
# Aggregate balanced accuracy by repeats
SVM_L1_balanced_accuracy_by_repeats <- SVM_L1_balanced_accuracy_by_folds %>%
  group_by(Study, Comparison_Group, Analysis_Type, Repeat_Number) %>%
  summarise(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
            Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T)) 
SVM_L1_balanced_accuracy <- SVM_L1_balanced_accuracy_by_repeats %>%
  group_by(Study, Comparison_Group, Analysis_Type) %>%
  summarise(Balanced_Accuracy_Across_Repeats = mean(Balanced_Accuracy_Across_Folds, na.rm=T),
            Balanced_Accuracy_Across_Repeats_SD = sd(Balanced_Accuracy_Across_Folds, na.rm=T)) 

SVM_L1_regularized_coefficients_by_folds <- pyarrow_feather$read_feather(glue("{data_path}/SVM_L1_Regularized_Coefficients.feather")) %>%
  dplyr::rename("Feature_Name" = "Feature Name", "Repeat_Number" = "Repeat Number")
SVM_L1_regularized_coefficients_by_repeats <- SVM_L1_regularized_coefficients_by_folds %>%
  group_by(Study, Comparison_Group, Analysis_Type, Feature_Name, Repeat_Number) %>%
  summarise(Coefficient_Across_Folds = mean(Coefficient, na.rm=T),
            Coefficient_Across_Folds_SD = sd(Coefficient, na.rm=T)) 
SVM_L1_regularized_coefficients <- SVM_L1_regularized_coefficients_by_repeats %>%
  group_by(Study, Comparison_Group, Analysis_Type, Feature_Name) %>%
  summarise(Coefficient_Across_Repeats = mean(Coefficient_Across_Folds, na.rm=T),
            Coefficient_Across_Repeats_SD = sd(Coefficient_Across_Folds, na.rm=T)) 

# Find number of zero vs non-zero components per group
SVM_L1_regularized_coefficients %>%
  group_by(Study, Comparison_Group) %>%
  summarise(num_zero = sum(Coefficient_Across_Repeats == 0),
            num_nonzero = sum(Coefficient_Across_Repeats != 0),
            total = n())

# Plot performance in normal SVM vs. L1-regularized SVM
SVM_L1_balanced_accuracy_by_repeats %>%
  dplyr::select(-Balanced_Accuracy_Across_Folds_SD) %>%
  plyr::rbind.fill(univariate_balanced_accuracy_by_repeats %>% filter(Analysis_Type=="Univariate_Combo")) %>%
  left_join(study_group_df) %>%
  mutate(Group_Nickname = factor(Group_Nickname, levels=c("SCZ", "BP", "ADHD", "ASD"))) %>%
  ggplot(data=., mapping=aes(x=Group_Nickname, y=Balanced_Accuracy_Across_Folds)) +
  geom_violin(aes(fill = Analysis_Type), position = position_dodge(width = 0.75),
              scale="width", width=0.6) +
  xlab("Comparison Group") +
  ylab("Balanced Accuracy\nby Repeat") +
  scale_fill_manual(values = c("cadetblue2", "darkcyan")) +
  geom_boxplot(aes(fill = Analysis_Type), color="black", width=0.1, position = position_dodge(width = 0.75)) +
  theme(legend.position = "none")
ggsave(glue("{plot_path}/Univariate_Combo_with_vs_without_Regularization.svg"),
       width=4, height=2.5, units="in", dpi=300)
