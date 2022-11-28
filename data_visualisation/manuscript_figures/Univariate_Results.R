################################################################################
# Load libraries
################################################################################

library(tidyverse)
library(icesTAF)
library(cowplot)
theme_set(theme_cowplot())

################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/Figure2/")
TAF::mkdir(plot_path)

source(paste0(github_dir, "helper_functions/data_prep_and_QC/QC_functions_univariate.R"))
source(paste0(github_dir, "helper_functions/Visualization.R"))

SCZ_data_path <- "~/data/UCLA_Schizophrenia/"
SCZ_rdata_path <- paste0(SCZ_data_path, "processed_data/Rdata/")
ASD_data_path <- "~/data/ABIDE_ASD/"
ASD_rdata_path <- paste0(ASD_data_path, "processed_data/Rdata/")

noise_proc_SCZ <- "AROMA+2P+GMR"
noise_proc_ASD <- "FC1000"
univariate_feature_set <- "catch22"
pairwise_feature_set <- "pyspi14_corrected"

ASD_brain_region_info <- read.csv(paste0(ASD_data_path, 
                                           "Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv"))

################################################################################
# Figure 2A univariate region-wise results
################################################################################

### Load univariate region-wise p-value data
# UCLA Schizophrenia
SCZ_ROI_pvals <- readRDS(paste0(SCZ_rdata_path, "ROI_wise_CV_linear_SVM_model_permutation_null_catch22_inv_prob_pvals.Rds")) %>%
  filter(Noise_Proc == noise_proc_SCZ) %>%
  ungroup() %>%
  mutate(bal_acc_p_adj_bonf = p.adjust(bal_acc_p, method="bonferroni")) %>%
  dplyr::select(grouping_var, bal_acc_p, bal_acc_p_adj, bal_acc_p_adj_bonf)
# ABIDE ASD
ASD_ROI_pvals <- readRDS(paste0(ASD_rdata_path, "ROI_wise_CV_linear_SVM_model_permutation_null_catch22_inv_prob_pvals.Rds")) %>%
  filter(Noise_Proc == noise_proc_ASD) %>%
  ungroup() %>%
  mutate(bal_acc_p_adj_bonf = p.adjust(bal_acc_p, method="bonferroni")) %>%
  dplyr::select(grouping_var, bal_acc_p, bal_acc_p_adj, bal_acc_p_adj_bonf)

### Identify significant brain regions
# UCLA Schizophrenia
sig_regions_SCZ_univar_ROI <- SCZ_ROI_pvals %>%
  filter(bal_acc_p_adj < 0.05) %>%
  pull(grouping_var)
# ABIDE ASD
sig_regions_ASD_univar_ROI <- ASD_ROI_pvals %>%
  filter(bal_acc_p_adj < 0.05) %>%
  pull(grouping_var)

### Load full univariate region-wise data
# UCLA Schizophrenia
SCZ_ROI_main_full <- readRDS(paste0(SCZ_rdata_path, "ROI_wise_CV_linear_SVM_catch22_inv_prob_balacc.Rds")) %>%
  filter(Noise_Proc == noise_proc_SCZ)
# ABIDE ASD
ASD_ROI_main_full <- readRDS(paste0(ASD_rdata_path, "ROI_wise_CV_linear_SVM_catch22_inv_prob_balacc.Rds")) %>%
  filter(Noise_Proc == noise_proc_ASD)

### Aggregate univariate region-wise data per repeat,  
# and filter to significant regions only
# UCLA Schizophrenia
SCZ_ROI_main_repeats <- SCZ_ROI_main_full %>%
  filter(grouping_var %in% sig_regions_SCZ_univar_ROI) %>%
  group_by(grouping_var, Sample_Type, Noise_Proc, repeat_number) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., SCZ_ROI_pvals) %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))
# ABIDE ASD
ASD_ROI_main_repeats <- ASD_ROI_main_full %>%
  filter(grouping_var %in% sig_regions_ASD_univar_ROI) %>%
  group_by(grouping_var, Sample_Type, Noise_Proc, repeat_number) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., ASD_ROI_pvals) %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))

### Aggregate all univariate region-wise data across repeats
# UCLA Schizophrenia
SCZ_ROI_main <- SCZ_ROI_main_repeats %>%
  group_by(grouping_var, Sample_Type, Noise_Proc) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T),
            mean_SD = sd(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., SCZ_ROI_pvals) %>%
  ungroup() %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))
# ABIDE ASD
ASD_ROI_main <- ASD_ROI_main_repeats %>%
  group_by(grouping_var, Sample_Type, Noise_Proc) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T),
            mean_SD = sd(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., ASD_ROI_pvals) %>%
  ungroup() %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))

### Load univariate region-wise empirical null distributions
# UCLA Schizophrenia
SCZ_ROI_null <- readRDS(paste0(SCZ_rdata_path, "UCLA_Schizophrenia_ROI_wise_model_permutation_null_catch22_inv_prob.Rds"))
# ABIDE ASD
ASD_ROI_null <- readRDS(paste0(ASD_rdata_path, "ABIDE_ASD_ROI_wise_model_permutation_null_catch22_inv_prob.Rds"))

# Plot boxplot with shaded null region per dataset

### UCLA boxplot with shaded null region
plot_boxplot_shaded_null(dataset_ID = "UCLA_Schizophrenia",
                         grouping_var_name = "Brain Region",
                         main_data_by_repeat = SCZ_ROI_main_repeats,
                         fill_color = "#F0224B",
                         null_mean_value = mean(SCZ_ROI_null$balanced_accuracy, na.rm=T),
                         null_SD_value = sd(SCZ_ROI_null$balanced_accuracy, na.rm=T))
ggsave(paste0(plot_path, "UCLA_Schizophrenia_Brain_Region_sig_boxplot.png"),
       width = 3.7, height = 2, units="in", dpi=300)

### ABIDE boxplot with shaded null region
plot_boxplot_shaded_null(dataset_ID = "ABIDE_ASD",
                         grouping_var_name = "Brain Region",
                         main_data_by_repeat = ASD_ROI_main_repeats,
                         fill_color = "#F0224B",
                         null_mean_value = mean(ASD_ROI_null$balanced_accuracy, na.rm=T),
                         null_SD_value = sd(ASD_ROI_null$balanced_accuracy, na.rm=T))
ggsave(paste0(plot_path, "ABIDE_ASD_Brain_Region_sig_boxplot.png"),
       width = 3.7, height = 2, units="in", dpi=300)

### UCLA Schizophrenia cortex and subcortex data in brain
SCZ_ROI_main$Type = "Main"
SCZ_ROI_null$Type = "Null"

SCZ_ROI_main_for_ggseg <- SCZ_ROI_main %>%
  dplyr::select(grouping_var) %>%
  dplyr::rename("label" = "grouping_var") %>%
  distinct() %>%
  mutate(fillyes = T) %>%
  mutate(label = ifelse(str_detect(label, "ctx-"),
                        gsub("-", "_", label),
                        as.character(label))) %>%
  mutate(label = gsub("ctx_", "", label))

# cortical data
plot_significant_regions_ggseg(dataset_ID = "UCLA_Schizophrenia",
                               atlas_name = "dk",
                               atlas_data = dk,
                               main_data_for_ggseg = SCZ_ROI_main_for_ggseg,
                               fill_color = "#F0224B")
ggsave(paste0(plot_path, "UCLA_Schizophrenia_ROI_wise_cortex_balanced_accuracy.png"),
       width = 4, height = 1.5, units = "in", dpi = 300)
# subcortical data
plot_significant_regions_ggseg(dataset_ID = "UCLA_Schizophrenia",
                               atlas_name = "aseg",
                               atlas_data = aseg,
                               main_data_for_ggseg = SCZ_ROI_main_for_ggseg,
                               fill_color = "#F0224B")
ggsave(paste0(plot_path, "UCLA_ROI_wise_subcortex_balanced_accuracy.png"),
       width = 4, height = 4, units = "in", dpi = 300)


### ABIDE ASD cortex and subcortex data in brain
ASD_ROI_main$Type = "Main"
ASD_ROI_null$Type = "Null"

ASD_ROI_main_for_ggseg <- ASD_ROI_main %>%
  dplyr::select(grouping_var) %>%
  dplyr::rename("Brain_Region" = "grouping_var") %>%
  distinct() %>%
  left_join(., ASD_brain_region_info) %>% 
  dplyr::rename("region" = "ggseg") %>%
  mutate(fillyes = T)

# cortical data
plot_significant_regions_ggseg(dataset_ID = "ABIDE_ASD",
                               atlas_name = "hoCort",
                               atlas_data = hoCort,
                               main_data_for_ggseg = ASD_ROI_main_for_ggseg,
                               fill_color = "#F0224B")
ggsave(paste0(plot_path, "ABIDE_ASD_ROI_wise_cortex_balanced_accuracy.png"),
       width = 4, height = 1.5, units = "in", dpi = 300)

################################################################################
# Misc: confirm no significant univariate features per dataset
################################################################################
# UCLA Schizophrenia
SCZ_feature_pvals <- readRDS(paste0(SCZ_rdata_path, "Feature_wise_CV_linear_SVM_model_permutation_null_catch22_inv_prob_pvals.Rds")) %>%
  filter(Noise_Proc == noise_proc_SCZ) %>%
  ungroup() %>%
  mutate(bal_acc_p_adj_bonf = p.adjust(bal_acc_p, method="bonferroni")) %>%
  dplyr::select(grouping_var, bal_acc_p, bal_acc_p_adj, bal_acc_p_adj_bonf)

SCZ_feature_pvals %>% filter(bal_acc_p_adj < 0.05)

# ABIDE ASD
ASD_feature_pvals <- readRDS(paste0(ASD_rdata_path, "Feature_wise_CV_linear_SVM_model_permutation_null_catch22_inv_prob_pvals.Rds")) %>%
  filter(Noise_Proc == noise_proc_ASD) %>%
  ungroup() %>%
  mutate(bal_acc_p_adj_bonf = p.adjust(bal_acc_p, method="bonferroni")) %>%
  dplyr::select(grouping_var, bal_acc_p, bal_acc_p_adj, bal_acc_p_adj_bonf)

ASD_feature_pvals %>% filter(bal_acc_p_adj < 0.05)

################################################################################
# Figure 2B univariate combo-wise results
################################################################################

### Load univariate combo-wise p-value data
# UCLA Schizophrenia
SCZ_combo_pvals <- readRDS(paste0(SCZ_rdata_path, "Combo_wise_CV_linear_SVM_model_permutation_null_catch22_inv_prob_pvals.Rds")) %>%
  filter(Noise_Proc == noise_proc_SCZ) %>%
  ungroup() %>%
  mutate(bal_acc_p_adj_bonf = p.adjust(bal_acc_p, method="bonferroni")) %>%
  dplyr::select(grouping_var, bal_acc_p, bal_acc_p_adj, bal_acc_p_adj_bonf)
# ABIDE ASD
ASD_combo_pvals <- readRDS(paste0(ASD_rdata_path, "Combo_wise_CV_linear_SVM_model_permutation_null_catch22_inv_prob_pvals.Rds")) %>%
  filter(Noise_Proc == noise_proc_ASD) %>%
  ungroup() %>%
  mutate(bal_acc_p_adj_bonf = p.adjust(bal_acc_p, method="bonferroni")) %>%
  dplyr::select(grouping_var, bal_acc_p, bal_acc_p_adj, bal_acc_p_adj_bonf)

### Load full univariate combo-wise data
# UCLA Schizophrenia
SCZ_combo_main_full <- readRDS(paste0(SCZ_rdata_path, 
                                       "Combo_wise_CV_linear_SVM_catch22_inv_prob_balacc.Rds")) %>%
  filter(Noise_Proc == noise_proc_SCZ) 
# ABIDE ASD
ASD_combo_main_full <- readRDS(paste0(ASD_rdata_path, 
                                        "Combo_wise_CV_linear_SVM_catch22_inv_prob_balacc.Rds")) %>%
  filter(Noise_Proc == noise_proc_ASD)

### Aggregate univariate combo-wise data per repeat
# UCLA Schizophrenia
SCZ_combo_main_repeats <- SCZ_combo_main_full %>%
  group_by(grouping_var, Sample_Type, Noise_Proc, repeat_number) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., SCZ_combo_pvals) %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))
# ABIDE ASD
ASD_combo_main_repeats <- ASD_combo_main_full %>%
  group_by(grouping_var, Sample_Type, Noise_Proc, repeat_number) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., ASD_combo_pvals) %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))

### Aggregate all univariate combo-wise data across repeats
# UCLA Schizophrenia
SCZ_combo_main <- SCZ_combo_main_repeats %>%
  group_by(grouping_var, Sample_Type, Noise_Proc) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T),
            mean_SD = sd(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., SCZ_combo_pvals) %>%
  ungroup() %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))
# ABIDE ASD
ASD_combo_main <- ASD_combo_main_repeats %>%
  group_by(grouping_var, Sample_Type, Noise_Proc) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T),
            mean_SD = sd(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., ASD_combo_pvals) %>%
  ungroup() %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))

### Load univariate combo-wise empirical null distributions
# UCLA Schizophrenia
SCZ_combo_null <- readRDS(paste0(SCZ_rdata_path, "UCLA_Schizophrenia_Combo_wise_model_permutation_null_catch22_inv_prob.Rds"))
# ABIDE ASD
ASD_combo_null <- readRDS(paste0(ASD_rdata_path, "ABIDE_ASD_Combo_wise_model_permutation_null_catch22_inv_prob.Rds"))

### Violin plot comparison for univariate results
# UCLA Schizophrenia
plot_univar_region_vs_combo_violin(dataset_ID = "UCLA_Schizophrenia",
                                   ROI_data_full = SCZ_ROI_main_full,
                                   ROI_pvals = SCZ_ROI_pvals,
                                   combo_data = SCZ_combo_main,
                                   region_fill_color = "#F0224B",
                                   combo_fill_color = "#9B51B4")
ggsave(paste0(plot_path, "UCLA_Schizophrenia_Univariate_Comparison_Boxplot.png"),
       width = 3, height = 2.75, units = "in", dpi = 300)
# ABIDE ASD
plot_univar_region_vs_combo_violin(dataset_ID = "ABIDE_ASD",
                                   ROI_data_full = ASD_ROI_main_full,
                                   ROI_pvals = ASD_ROI_pvals,
                                   combo_data = ASD_combo_main,
                                   region_fill_color = "#F0224B",
                                   combo_fill_color = "#9B51B4")
ggsave(paste0(plot_path, "ABIDE_ASD_Univariate_Comparison_Boxplot.png"),
       width = 3, height = 2.75, units = "in", dpi = 300)
