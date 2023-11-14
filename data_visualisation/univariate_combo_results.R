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
library(see)
library(ggpp)

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
# BOWTIE PLOT 

# Summarise performance metrics
univariate_p_values %>%
  filter(Analysis_Type == "Univariate_Combo") %>% 
  mutate(summarymean = round(Balanced_Accuracy_Across_Folds*100,1),
         summarysd = round(Balanced_Accuracy_Across_Folds_SD*100,1)) %>%
  dplyr::select(Comparison_Group, summarymean, summarysd) %>%
  arrange(desc(summarymean))

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
                                    T ~ "Individual\nBrain\nRegions"),
         Analysis_Sig = paste0(Comparison_Group, "_", p_value_HolmBonferroni < 0.05)) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BP", "ADHD", "ASD")),
         Analysis_Label = factor(Analysis_Label, levels = c("Individual\ncatch25\nFeatures",
                                                            "All\nRegions \u00D7\nFeatures",
                                                            "Individual\nBrain\nRegions"))) %>%
  group_by(Study, Comparison_Group, group_ID) %>%
  mutate(Comparison_Sig = paste0(Comparison_Group, "_", p_value_HolmBonferroni[Analysis_Type != "Univariate_Combo"]<0.05)) %>% 
  ggplot(data=., mapping=aes(x=Analysis_Label, y=100*Balanced_Accuracy_Across_Folds, 
                             group=group_ID, color=Comparison_Sig)) +
  geom_line(alpha=0.5) +
  geom_point(aes(size=Analysis_Type, color=Analysis_Sig)) +
  ylab("Mean Balanced Accuracy Across Repeats (%)") +
  xlab("Classifier Inputs") +
  geom_hline(yintercept = 50, linetype=2) +
  scale_x_discrete(expand=c(0.05,0.05,0.05,0.05)) +
  scale_size_manual(values=c("Univariate_Combo"=4,
                             "Univariate_TS_Feature"=1.5,
                             "Univariate_Brain_Region" = 1.5)) +
  scale_color_manual(values = c("SCZ_TRUE"="#9d60a8", 
                                "BP_TRUE"="#2F77C0", 
                                "ADHD_TRUE"="#e45075",
                                "ASD_TRUE"="#E28328",
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
# Take only the top-performing region and time-series feature from training data
univariate_balanced_accuracy_best_in_sample <- univariate_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Analysis_Type, group_var, kernel) %>%
  reframe(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
          Training_Balanced_Accuracy_Across_Folds = mean(Training_Balanced_Accuracy, na.rm=T),
          Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T)) %>%
  group_by(Study, Comparison_Group, Analysis_Type) %>%
  filter(Training_Balanced_Accuracy_Across_Folds == max(Training_Balanced_Accuracy_Across_Folds))

univariate_p_values_best_in_sample <- univariate_balanced_accuracy_best_in_sample %>%
  filter(group_var!="Combo") %>%
  left_join(., univariate_p_values) %>%
  group_by(Study, Comparison_Group) %>%
  mutate(group_ID = paste0(Comparison_Group, "_", row_number()))

# Bowtie plot comparing each brain region and feature to combo-wise performance
num_comparisons <- univariate_balanced_accuracy_best_in_sample %>%
  filter(Analysis_Type != "Univariate_Combo") %>%
  group_by(Study, Comparison_Group) %>%
  count()

univariate_p_values %>%
  filter(Analysis_Type == "Univariate_Combo") %>% 
  left_join(., num_comparisons) %>%
  expandRows("n") %>%
  group_by(Study, Comparison_Group) %>%
  mutate(group_ID = paste0(Comparison_Group, "_", row_number())) %>%
  plyr::rbind.fill(., univariate_p_values_best_in_sample) %>%
  rowwise() %>%
  mutate(Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                      Comparison_Group == "Bipolar" ~ "BP",
                                      T ~ Comparison_Group),
         Analysis_Label = case_when(Analysis_Type == "Univariate_Combo" ~ "All\nRegions \u00D7\nFeatures",
                                    Analysis_Type == "Univariate_TS_Feature" ~ "Individual\ncatch25\nFeatures",
                                    T ~ "Individual\nBrain\nRegions"),
         Analysis_Sig = paste0(Comparison_Group, "_", p_value_HolmBonferroni < 0.05)) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BP", "ADHD", "ASD")),
         Analysis_Label = factor(Analysis_Label, levels = c("Individual\ncatch25\nFeatures",
                                                            "All\nRegions \u00D7\nFeatures",
                                                            "Individual\nBrain\nRegions"))) %>%
  group_by(Study, Comparison_Group, group_ID) %>%
  mutate(Comparison_Sig = paste0(Comparison_Group, "_", p_value_HolmBonferroni[Analysis_Type != "Univariate_Combo"]<0.05)) %>% 
  ggplot(data=., mapping=aes(x=Analysis_Label, y=100*Balanced_Accuracy_Across_Folds, 
                             group=group_ID, color=Comparison_Sig)) +
  geom_line(alpha=0.5) +
  geom_point(aes(size=Analysis_Type, color=Analysis_Sig)) +
  ylab("Mean Balanced Accuracy Across Repeats (%)") +
  xlab("Classifier Inputs") +
  geom_hline(yintercept = 50, linetype=2) +
  scale_x_discrete(expand=c(0.05,0.05,0.05,0.05)) +
  scale_size_manual(values=c("Univariate_Combo"=4,
                             "Univariate_TS_Feature"=1.5,
                             "Univariate_Brain_Region" = 1.5)) +
  scale_color_manual(values = c("SCZ_TRUE"="#9d60a8", 
                                "BP_TRUE"="#2F77C0", 
                                "ADHD_TRUE"="#e45075",
                                "ASD_TRUE"="#E28328",
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
ggsave(glue("{plot_path}/univariate_bowtie_balanced_accuracy_supplement.svg"),
       width=9, height=2.5, units="in", dpi=300)

################################################################################
# Compare performance with L1-regularized SVM balanced accuracy implemented with LinearSVC()
SVM_L1_balanced_accuracy_by_folds <- pyarrow_feather$read_feather(glue("{data_path}/SVM_L1_Regularized_Balanced_Accuracy.feather"))
SVM_L1_coefficients <- pyarrow_feather$read_feather(glue("{data_path}/SVM_L1_Regularized_Coefficients.feather"))

# Find number of unique features by disorder
SVM_L1_coefficients %>%
  group_by(Comparison_Group) %>%
  distinct(`Feature Name`) %>%
  count()

# Aggregate balanced accuracy by repeats
SVM_L1_balanced_accuracy <- SVM_L1_balanced_accuracy_by_folds %>%
  group_by(Study, Comparison_Group, Analysis_Type) %>%
  summarise(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
            Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T)) 

SVM_L1_regularized_coefficients_by_folds <- pyarrow_feather$read_feather(glue("{data_path}/SVM_L1_Regularized_Coefficients.feather")) %>%
  dplyr::rename("Feature_Name" = "Feature Name", "Repeat_Number" = "Repeat Number")
SVM_L1_regularized_coefficients <- SVM_L1_regularized_coefficients_by_folds %>%
  group_by(Study, Comparison_Group, Analysis_Type, Feature_Name) %>%
  summarise(Coefficient_Across_Folds = mean(Coefficient, na.rm=T),
            Coefficient_Across_Folds_SD = sd(Coefficient, na.rm=T)) 

# Find number of zero vs non-zero components per group
SVM_L1_regularized_coefficients %>%
  group_by(Study, Comparison_Group, Analysis_Type) %>%
  summarise(num_zero = sum(Coefficient_Across_Folds == 0),
            num_nonzero = sum(Coefficient_Across_Folds != 0),
            total = n())

# Plot performance in normal SVM vs. L1-regularized SVM
SVM_L1_balanced_accuracy_by_folds %>%
  plyr::rbind.fill(univariate_balanced_accuracy_all_folds %>% filter(Analysis_Type=="Univariate_Combo")) %>%
  left_join(study_group_df) %>%
  mutate(Group_Nickname = factor(Group_Nickname, levels=c("SCZ", "BP", "ADHD", "ASD"))) %>%
  ggplot(data=., mapping=aes(x=Analysis_Type, y=100*Balanced_Accuracy)) +
  xlab("SVM Type") +
  ylab("Balanced Accuracy (%)") +
  geom_violinhalf(aes(fill = Analysis_Type), 
                  position = position_nudge(x=0.2),
                  scale="width", width=0.6)  +
  geom_boxplot(width=0.1, notch=FALSE, notchwidth = 0.4, outlier.shape = NA,
               fill=NA, color="white",
               position = position_nudge(x=0.27), coef = 0) +
  geom_hline(yintercept = 50, linetype=2, alpha=0.5) +
  facet_grid(. ~ Group_Nickname) +
  geom_point(aes(color = Analysis_Type), position = position_jitter(width=0.1),
             size = 1.75, alpha=0.6, stroke=0) +
  labs(fill="", color="") +
  scale_fill_manual(values = c("#319719", "#8dc57f")) +
  scale_color_manual(values = c("#319719", "#8dc57f")) +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(face="bold"),
        strip.background = element_blank())
ggsave(glue("{plot_path}/Univariate_Combo_with_vs_without_Regularization.svg"),
       width=6, height=3.5, units="in", dpi=300)

################################################################################
# How do samples compare in PC space?

# Start list for PCA score results
pca_scores_list <- list()

for (i in 1:nrow(study_group_df)) {
  dataset_ID <- study_group_df$Study[i]
  comparison_group <- study_group_df$Comparison_Group[i]
  group_nickname <- study_group_df$Group_Nickname
  
  if (dataset_ID == "UCLA_CNP") {
    input_catch25_data = UCLA_CNP_catch25
  } else {
    input_catch25_data = ABIDE_ASD_catch25
  }
  
  # Prep data for PCA
  data_for_PCA <- input_catch25_data %>%
    filter(Diagnosis %in% c("Control", comparison_group)) %>%
    left_join(., TS_feature_info, by=c("names"="feature_name")) %>%
    mutate(unique_ID = paste0(Figure_name, "__", Brain_Region), .keep="unused") %>%
    dplyr::select(unique_ID, Sample_ID, Diagnosis, values) %>%
    pivot_wider(id_cols=c(Sample_ID, Diagnosis), names_from=unique_ID, values_from=values) %>%
    mutate(Diagnosis=factor(Diagnosis, levels=c("Control", comparison_group)))
  
  # Compute PCA
  pca_res <- PCA(select(data_for_PCA, c(-Sample_ID, -Diagnosis)), ncp=25, graph = FALSE, scale.unit = TRUE)
  pca_scores <- as.data.frame(pca_res$ind$coord) %>%
    mutate(Sample_ID = data_for_PCA$Sample_ID,
           Diagnosis = data_for_PCA$Diagnosis,
           Comparison_Group = comparison_group,
           Study = dataset_ID)
  
  pca_scores_list <- list.append(pca_scores_list, pca_scores)
}

pca_scores <- do.call(plyr::rbind.fill, pca_scores_list)

pyarrow_feather$write_feather(pca_scores, "~/data/TS_feature_manuscript/Univariate_combo_first25_PCs.feather")

pca_scores %>%
  mutate(Diagnosis = case_when(Diagnosis == "Schizophrenia" ~ "SCZ",
                               Diagnosis == "Bipolar" ~ "BP",
                               T ~ Diagnosis),
         Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                      Comparison_Group == "Bipolar" ~ "BP",
                                      T ~ Comparison_Group)) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BP", "ADHD", "ASD"))) %>%
  ggplot(data=., mapping=aes(x=`Dim.1`, y=`Dim.2`, color=Diagnosis)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype=2) +
  geom_hline(yintercept = 0, linetype=2) +
  stat_chull(aes(color = Diagnosis, fill = Diagnosis), 
             alpha = 0.1, geom = "polygon")+
  ylab("PC2") +
  xlab("PC1") +
  facet_wrap(Comparison_Group ~ ., scales="free", nrow=1) +
  scale_fill_manual(values = c("SCZ"="#9d60a8", 
                                "BP"="#2F77C0", 
                                "ADHD"="#e45075",
                                "ASD" = "#E28328",
                                "Control" = "#439E47")) +
  scale_color_manual(values = c("SCZ"="#9d60a8", 
                                "BP"="#2F77C0", 
                                "ADHD"="#e45075",
                                "ASD" = "#E28328",
                                "Control" = "#439E47")) +
  theme(legend.position = "bottom",
        strip.text = element_text(face="bold"),
        strip.background = element_blank())
ggsave(glue("{plot_path}/Univariate_combo_PCA_biplots.svg"),
       width=11, height=3.5, units="in", dpi=300)

################################################################################
# Comparing classification performance using the first 25 PCs per disorder

univariate_combo_first25_PCs_balanced_accuracy <- pyarrow_feather$read_feather("~/data/TS_feature_manuscript/Univariate_combo_first25_PCs_balanced_accuracy.feather")

univariate_combo_first25_PCs_balanced_accuracy %>%
  plyr::rbind.fill(univariate_balanced_accuracy_all_folds %>% filter(Analysis_Type=="Univariate_Combo")) %>%
  left_join(study_group_df) %>%
  mutate(Group_Nickname = factor(Group_Nickname, levels=c("SCZ", "BP", "ADHD", "ASD"))) %>%
  ggplot(data=., mapping=aes(x=Analysis_Type, y=100*Balanced_Accuracy)) +
  xlab("SVM Type") +
  ylab("Balanced Accuracy (%)") +
  geom_violinhalf(aes(fill = Analysis_Type), 
                  position = position_nudge(x=0.2),
                  scale="width", width=0.6)  +
  geom_boxplot(width=0.1, notch=FALSE, notchwidth = 0.4, outlier.shape = NA,
               fill=NA, color="white",
               position = position_nudge(x=0.27), coef = 0) +
  geom_hline(yintercept = 50, linetype=2, alpha=0.5) +
  facet_grid(. ~ Group_Nickname) +
  geom_point(aes(color = Analysis_Type), position = position_jitter(width=0.1),
             size = 1.75, alpha=0.6, stroke=0) +
  labs(fill="", color="") +
  scale_fill_manual(values = c("#a86ba3", "#cfafcd")) +
  scale_color_manual(values = c("#a86ba3", "#cfafcd")) +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(face="bold"),
        strip.background = element_blank())
ggsave(glue("{plot_path}/Univariate_Combo_vs_First_25_PCs.svg"),
       width=6, height=3.5, units="in", dpi=300)