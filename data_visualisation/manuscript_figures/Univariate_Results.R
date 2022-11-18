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
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/")
icesTAF::mkdir(plot_path)

source(paste0(github_dir, "helper_functions/data_prep_and_QC/QC_functions_univariate.R"))
source(paste0(github_dir, "helper_functions/Visualization.R"))

UCLA_data_path <- "~/data/UCLA_Schizophrenia/"
UCLA_rdata_path <- paste0(UCLA_data_path, "processed_data/Rdata/")
ABIDE_data_path <- "~/data/ABIDE_ASD/"
ABIDE_rdata_path <- paste0(ABIDE_data_path, "processed_data/Rdata/")


ABIDE_brain_region_info <- read.csv(paste0(ABIDE_data_path, 
                                           "Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv"))

icesTAF::mkdir(paste0(plot_path, "Figure2/"))

################################################################################
# Figure 2A univariate region-wise results
################################################################################
icesTAF::mkdir(paste0(plot_path, "Figure2/"))

### Load univariate region-wise p-value data
# UCLA Schizophrenia
UCLA_ROI_pvals <- readRDS(paste0(UCLA_rdata_path, "ROI_wise_CV_linear_SVM_model_permutation_null_catch22_inv_prob_pvals.Rds")) %>%
  filter(Noise_Proc == "AROMA+2P+GMR") %>%
  ungroup() %>%
  mutate(bal_acc_p_adj_bonf = p.adjust(bal_acc_p, method="bonferroni")) %>%
  dplyr::select(grouping_var, bal_acc_p, bal_acc_p_adj, bal_acc_p_adj_bonf)
# ABIDE ASD
ABIDE_ROI_pvals <- readRDS(paste0(ABIDE_rdata_path, "ROI_wise_CV_linear_SVM_model_permutation_null_catch22_inv_prob_pvals.Rds")) %>%
  filter(Noise_Proc == "FC1000") %>%
  ungroup() %>%
  mutate(bal_acc_p_adj_bonf = p.adjust(bal_acc_p, method="bonferroni")) %>%
  dplyr::select(grouping_var, bal_acc_p, bal_acc_p_adj, bal_acc_p_adj_bonf)

### Identify significant brain regions
# UCLA Schizophrenia
sig_regions_UCLA_univar_ROI <- UCLA_ROI_pvals %>%
  filter(bal_acc_p_adj < 0.05) %>%
  pull(grouping_var)
# ABIDE ASD
sig_regions_ABIDE_univar_ROI <- ABIDE_ROI_pvals %>%
  filter(bal_acc_p_adj < 0.05) %>%
  pull(grouping_var)

### Load full univariate region-wise data
# UCLA Schizophrenia
UCLA_ROI_main_full <- readRDS(paste0(UCLA_rdata_path, "ROI_wise_CV_linear_SVM_catch22_inv_prob_balacc.Rds")) %>%
  filter(Noise_Proc == "AROMA+2P+GMR")
# ABIDE ASD
ABIDE_ROI_main_full <- readRDS(paste0(ABIDE_rdata_path, "ROI_wise_CV_linear_SVM_catch22_inv_prob_balacc.Rds")) %>%
  filter(Noise_Proc == "FC1000")

### Aggregate univariate region-wise data per repeat,  
# and filter to significant regions only
# UCLA Schizophrenia
UCLA_ROI_main_repeats <- UCLA_ROI_main_full %>%
  filter(grouping_var %in% sig_regions_UCLA_univar_ROI) %>%
  group_by(grouping_var, Sample_Type, Noise_Proc, repeat_number) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., UCLA_ROI_pvals) %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))
# ABIDE ASD
ABIDE_ROI_main_repeats <- ABIDE_ROI_main_full %>%
  filter(grouping_var %in% sig_regions_ABIDE_univar_ROI) %>%
  group_by(grouping_var, Sample_Type, Noise_Proc, repeat_number) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., ABIDE_ROI_pvals) %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))

### Aggregate all univariate region-wise data across repeats
# UCLA Schizophrenia
UCLA_ROI_main <- UCLA_ROI_main_repeats %>%
  group_by(grouping_var, Sample_Type, Noise_Proc) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T),
            mean_SD = sd(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., UCLA_ROI_pvals) %>%
  ungroup() %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))
# ABIDE ASD
ABIDE_ROI_main <- ABIDE_ROI_main_repeats %>%
  group_by(grouping_var, Sample_Type, Noise_Proc) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T),
            mean_SD = sd(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., ABIDE_ROI_pvals) %>%
  ungroup() %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))

### Load univariate region-wise empirical null distributions
# UCLA Schizophrenia
UCLA_ROI_null <- readRDS(paste0(UCLA_rdata_path, "UCLA_Schizophrenia_ROI_wise_model_permutation_null_catch22_inv_prob.Rds"))
# ABIDE ASD
ABIDE_ROI_null <- readRDS(paste0(ABIDE_rdata_path, "ABIDE_ASD_ROI_wise_model_permutation_null_catch22_inv_prob.Rds"))

# Plot boxplot with shaded null region per dataset

### UCLA boxplot with shaded null region
plot_boxplot_shaded_null(dataset_ID = "UCLA_Schizophrenia",
                         grouping_var_name = "Brain Region",
                         main_data_by_repeat = UCLA_ROI_main_repeats,
                         fill_color = "#F0224B",
                         null_mean_value = mean(UCLA_ROI_null$balanced_accuracy, na.rm=T),
                         null_SD_value = sd(UCLA_ROI_null$balanced_accuracy, na.rm=T))
ggsave(paste0(plot_path, "Figure2/UCLA_Schizophrenia_Brain_Region_sig_boxplot.png"),
       width = 3.7, height = 2, units="in", dpi=300)

### ABIDE boxplot with shaded null region
plot_boxplot_shaded_null(dataset_ID = "ABIDE_ASD",
                         grouping_var_name = "Brain Region",
                         main_data_by_repeat = ABIDE_ROI_main_repeats,
                         fill_color = "#F0224B",
                         null_mean_value = mean(ABIDE_ROI_null$balanced_accuracy, na.rm=T),
                         null_SD_value = sd(ABIDE_ROI_null$balanced_accuracy, na.rm=T))
ggsave(paste0(plot_path, "Figure2/ABIDE_ASD_Brain_Region_sig_boxplot.png"),
       width = 3.7, height = 1.6, units="in", dpi=300)

### UCLA Schizophrenia cortex and subcortex data in brain
UCLA_ROI_main$Type = "Main"
UCLA_ROI_null$Type = "Null"

UCLA_ROI_main_for_ggseg <- UCLA_ROI_main %>%
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
                               main_data_for_ggseg = UCLA_ROI_main_for_ggseg,
                               fill_color = "#F0224B")
ggsave(paste0(plot_path, "Figure2/UCLA_Schizophrenia_ROI_wise_cortex_balanced_accuracy.png"),
       width = 4, height = 1.5, units = "in", dpi = 300)
# subcortical data
plot_significant_regions_ggseg(dataset_ID = "UCLA_Schizophrenia",
                               atlas_name = "aseg",
                               atlas_data = aseg,
                               main_data_for_ggseg = UCLA_ROI_main_for_ggseg,
                               fill_color = "#F0224B")
ggsave(paste0(plot_path, "Figure2/UCLA_ROI_wise_subcortex_balanced_accuracy.png"),
       width = 4, height = 4, units = "in", dpi = 300)


### ABIDE ASD cortex and subcortex data in brain
ABIDE_ROI_main$Type = "Main"
ABIDE_ROI_null$Type = "Null"

ABIDE_ROI_main_for_ggseg <- ABIDE_ROI_main %>%
  dplyr::select(grouping_var) %>%
  dplyr::rename("Brain_Region" = "grouping_var") %>%
  distinct() %>%
  left_join(., ABIDE_brain_region_info) %>% 
  dplyr::rename("region" = "ggseg") %>%
  mutate(fillyes = T)

# cortical data
plot_significant_regions_ggseg(dataset_ID = "ABIDE_ASD",
                               atlas_name = "hoCort",
                               atlas_data = hoCort,
                               main_data_for_ggseg = ABIDE_ROI_main_for_ggseg,
                               fill_color = "#F0224B")
ggsave(paste0(plot_path, "Figure2/ABIDE_ASD_ROI_wise_cortex_balanced_accuracy.png"),
       width = 4, height = 1.5, units = "in", dpi = 300)


################################################################################
# Figure 2B univariate combo-wise results
################################################################################

### Load univariate combo-wise p-value data
# UCLA Schizophrenia
UCLA_combo_pvals <- readRDS(paste0(UCLA_rdata_path, "Combo_wise_CV_linear_SVM_model_permutation_null_catch22_inv_prob_pvals.Rds")) %>%
  filter(Noise_Proc == "AROMA+2P+GMR") %>%
  ungroup() %>%
  mutate(bal_acc_p_adj_bonf = p.adjust(bal_acc_p, method="bonferroni")) %>%
  dplyr::select(grouping_var, bal_acc_p, bal_acc_p_adj, bal_acc_p_adj_bonf)
# ABIDE ASD
ABIDE_combo_pvals <- readRDS(paste0(ABIDE_rdata_path, "Combo_wise_CV_linear_SVM_model_permutation_null_catch22_inv_prob_pvals.Rds")) %>%
  filter(Noise_Proc == "FC1000") %>%
  ungroup() %>%
  mutate(bal_acc_p_adj_bonf = p.adjust(bal_acc_p, method="bonferroni")) %>%
  dplyr::select(grouping_var, bal_acc_p, bal_acc_p_adj, bal_acc_p_adj_bonf)

### Load full univariate combo-wise data
# UCLA Schizophrenia
UCLA_combo_main_full <- readRDS(paste0(UCLA_rdata_path, 
                                       "Combo_wise_CV_linear_SVM_catch22_inv_prob_balacc.Rds")) %>%
  filter(Noise_Proc == "AROMA+2P+GMR") 
# ABIDE ASD
ABIDE_combo_main_full <- readRDS(paste0(ABIDE_rdata_path, 
                                        "Combo_wise_CV_linear_SVM_catch22_inv_prob_balacc.Rds")) %>%
  filter(Noise_Proc == "FC1000") 

### Aggregate univariate combo-wise data per repeat
# UCLA Schizophrenia
UCLA_combo_main_repeats <- UCLA_combo_main_full %>%
  group_by(grouping_var, Sample_Type, Noise_Proc, repeat_number) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., UCLA_combo_pvals) %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))
# ABIDE ASD
ABIDE_combo_main_repeats <- ABIDE_combo_main_full %>%
  group_by(grouping_var, Sample_Type, Noise_Proc, repeat_number) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., ABIDE_combo_pvals) %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))

### Aggregate all univariate combo-wise data across repeats
# UCLA Schizophrenia
UCLA_combo_main <- UCLA_combo_main_repeats %>%
  group_by(grouping_var, Sample_Type, Noise_Proc) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T),
            mean_SD = sd(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., UCLA_combo_pvals) %>%
  ungroup() %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))
# ABIDE ASD
ABIDE_combo_main <- ABIDE_combo_main_repeats %>%
  group_by(grouping_var, Sample_Type, Noise_Proc) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T),
            mean_SD = sd(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., ABIDE_combo_pvals) %>%
  ungroup() %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))

### Load univariate combo-wise empirical null distributions
# UCLA Schizophrenia
UCLA_combo_null <- readRDS(paste0(UCLA_rdata_path, "UCLA_Schizophrenia_Combo_wise_model_permutation_null_catch22_inv_prob.Rds"))
# ABIDE ASD
ABIDE_combo_null <- readRDS(paste0(ABIDE_rdata_path, "ABIDE_ASD_Combo_wise_model_permutation_null_catch22_inv_prob.Rds"))

### Violin plot comparison for univariate results
# UCLA Schizophrenia
plot_univar_region_vs_combo_violin(dataset_ID = "UCLA_Schizophrenia",
                                   ROI_data_full = UCLA_ROI_main_full,
                                   ROI_pvals = UCLA_ROI_pvals,
                                   combo_data = UCLA_combo_main,
                                   region_fill_color = "#F0224B",
                                   combo_fill_color = "#9B51B4")
ggsave(paste0(plot_path, "Figure2/UCLA_Univariate_Comparison_Boxplot.png"),
       width = 3, height = 2.75, units = "in", dpi = 300)
# ABIDE ASD
plot_univar_region_vs_combo_violin(dataset_ID = "ABIDE_ASD",
                                   ROI_data_full = ABIDE_ROI_main_full,
                                   ROI_pvals = ABIDE_ROI_pvals,
                                   combo_data = ABIDE_combo_main,
                                   region_fill_color = "#F0224B",
                                   combo_fill_color = "#9B51B4")
ggsave(paste0(plot_path, "Figure2/ABIDE_Univariate_Comparison_Boxplot.png"),
       width = 3, height = 2.75, units = "in", dpi = 300)
