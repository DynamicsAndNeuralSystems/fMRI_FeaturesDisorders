################################################################################
# Load libraries
################################################################################

library(tidyverse)
library(icesTAF)
# library(ggseg)
library(tidyverse)
library(plotly)
# library(ggseg3d)
# library(ggsegHO)
library(reshape2)
library(patchwork)
library(cowplot)
library(broom)
theme_set(theme_cowplot())



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

################################################################################
# Figure 3A pairwise SPI-wise analysis results
################################################################################
icesTAF::mkdir(paste0(plot_path, "Figure3/"))

### Load pairwise SPI-wise p-value data
# UCLA Schizophrenia
UCLA_SPI_pvals <- readRDS(paste0(UCLA_rdata_path, "SPI_wise_CV_linear_SVM_model_permutation_null_pyspi14_inv_prob_pvals.Rds")) %>%
  filter(Noise_Proc == "AROMA+2P+GMR") %>%
  ungroup() %>%
  mutate(bal_acc_p_adj_bonf = p.adjust(bal_acc_p, method="bonferroni")) %>%
  dplyr::select(grouping_var, bal_acc_p, bal_acc_p_adj, bal_acc_p_adj_bonf)
# ABIDE ASD
ABIDE_SPI_pvals <- readRDS(paste0(ABIDE_rdata_path, "SPI_wise_CV_linear_SVM_model_permutation_null_pyspi14_inv_prob_pvals.Rds")) %>%
  filter(Noise_Proc == "FC1000") %>%
  ungroup() %>%
  mutate(bal_acc_p_adj_bonf = p.adjust(bal_acc_p, method="bonferroni")) %>%
  dplyr::select(grouping_var, bal_acc_p, bal_acc_p_adj, bal_acc_p_adj_bonf)

### Load full pairwise SPI-wise data
# UCLA Schizophrenia
UCLA_SPI_main_full <- readRDS(paste0(UCLA_rdata_path, "SPI_wise_CV_linear_SVM_pyspi14_inv_prob_balacc.Rds"))%>%
  filter(Noise_Proc == "AROMA+2P+GMR") 
# ABIDE ASD
ABIDE_SPI_main_full <- readRDS(paste0(ABIDE_rdata_path, "SPI_wise_CV_linear_SVM_pyspi14_inv_prob_balacc.Rds")) %>%
  filter(Noise_Proc == "FC1000") 

### Aggregate pairwise SPI-wise data per repeat 
# And filter by BH adjusted p-value < 0.05
# UCLA Schizophrenia
UCLA_SPI_main_repeats <- UCLA_SPI_main_full %>%
  group_by(grouping_var, Sample_Type, Noise_Proc, repeat_number) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., UCLA_SPI_pvals) %>%
  ungroup() %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))

UCLA_SPI_main_repeats_sig <- UCLA_SPI_main_repeats %>%
  filter(bal_acc_p_adj < 0.05)

# ABIDE ASD
ABIDE_SPI_main_repeats <- ABIDE_SPI_main_full %>%
  group_by(grouping_var, Sample_Type, Noise_Proc, repeat_number) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., ABIDE_SPI_pvals) %>%
  ungroup() %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))
ABIDE_SPI_main_repeats_sig <- ABIDE_SPI_main_repeats %>%
  filter(bal_acc_p_adj < 0.05)

### Aggregate all pairwise SPI-wise data across repeats
# UCLA Schizophrenia
UCLA_SPI_main <- UCLA_SPI_main_repeats %>%
  group_by(grouping_var, Sample_Type, Noise_Proc) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T),
            mean_SD = sd(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., UCLA_SPI_pvals) %>%
  ungroup() %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))

UCLA_SPI_main_sig <- UCLA_SPI_main %>%
  filter(bal_acc_p_adj < 0.05)

# ABIDE ASD
ABIDE_SPI_main <- ABIDE_SPI_main_repeats %>%
  group_by(grouping_var, Sample_Type, Noise_Proc) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T),
            mean_SD = sd(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., ABIDE_SPI_pvals) %>%
  ungroup() %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))

ABIDE_SPI_main_sig <- ABIDE_SPI_main %>%
  filter(bal_acc_p_adj < 0.05)

### Load pairwise SPI-wise empirical null distributions
# UCLA Schizophrenia
UCLA_SPI_null <- readRDS(paste0(UCLA_rdata_path, "UCLA_Schizophrenia_SPI_wise_model_permutation_null_pyspi14_inv_prob.Rds"))
# ABIDE ASD
ABIDE_SPI_null <- readRDS(paste0(ABIDE_rdata_path, "ABIDE_ASD_SPI_wise_model_permutation_null_pyspi14_inv_prob.Rds"))

### SPI boxplot with shaded null region
# UCLA Schizophrenia
plot_boxplot_shaded_null(dataset_ID = "UCLA_Schizophrenia",
                         grouping_var_name = "Pairwise SPI",
                         main_data_by_repeat = UCLA_SPI_main_repeats_sig,
                         fill_color = "chartreuse3",
                         wrap_length = 100,
                         null_mean_value = mean(UCLA_SPI_null$balanced_accuracy, na.rm=T),
                         null_SD_value = sd(UCLA_SPI_null$balanced_accuracy, na.rm=T))
ggsave(paste0(plot_path, "Figure3/UCLA_Schizophrenia_SPI_sig_boxplot.png"),
       width = 8, height = 2.5, units="in", dpi=300)
# ABIDE ASD
plot_boxplot_shaded_null(dataset_ID = "ABIDE_ASD",
                         grouping_var_name = "Pairwise SPI",
                         main_data_by_repeat = ABIDE_SPI_main_repeats_sig,
                         fill_color = "chartreuse3",
                         wrap_length = 100,
                         null_mean_value = mean(ABIDE_SPI_null$balanced_accuracy, na.rm=T),
                         null_SD_value = sd(ABIDE_SPI_null$balanced_accuracy, na.rm=T))
ggsave(paste0(plot_path, "Figure3/ABIDE_ASD_SPI_sig_boxplot.png"),
       width = 8, height = 2.5, units="in", dpi=300)


################################################################################
# Figure 3B violin plot comparing all balanced accuracy values
################################################################################

### Load combo univariate and pairwise SPI-wise p-value data
# UCLA Schizophrenia
UCLA_SPI_combo_pvals <- readRDS(paste0(UCLA_rdata_path, "univariate_catch22_pairwise_pyspi14_CV_linear_SVM_model_permutation_null_inv_prob_pvals.Rds")) %>%
  filter(Noise_Proc == "AROMA+2P+GMR") %>%
  ungroup() %>%
  mutate(bal_acc_p_adj_bonf = p.adjust(bal_acc_p, method="bonferroni")) %>%
  dplyr::select(grouping_var, bal_acc_p, bal_acc_p_adj, bal_acc_p_adj_bonf)
# ABIDE ASD
ABIDE_SPI_combo_pvals <- readRDS(paste0(ABIDE_rdata_path, "univariate_catch22_pairwise_pyspi14_CV_linear_SVM_model_permutation_null_inv_prob_pvals.Rds")) %>%
  filter(Noise_Proc == "FC1000") %>%
  ungroup() %>%
  mutate(bal_acc_p_adj_bonf = p.adjust(bal_acc_p, method="bonferroni")) %>%
  dplyr::select(grouping_var, bal_acc_p, bal_acc_p_adj, bal_acc_p_adj_bonf)

### Load full combo univariate and pairwise SPI-wise data
# UCLA Schizophrenia
UCLA_SPI_combo_main_full <- readRDS(paste0(UCLA_rdata_path, "univariate_catch22_pairwise_pyspi14_CV_linear_SVM_inv_prob_balacc.Rds"))%>%
  filter(Noise_Proc == "AROMA+2P+GMR") %>%
  dplyr::rename("grouping_var" = "SPI")
# ABIDE ASD
ABIDE_SPI_combo_main_full <- readRDS(paste0(ABIDE_rdata_path, "univariate_catch22_pairwise_pyspi14_CV_linear_SVM_inv_prob_balacc.Rds")) %>%
  filter(Noise_Proc == "FC1000") %>%
  dplyr::rename("grouping_var" = "SPI")

### Aggregate combo univariate and pairwise SPI-wise data per repeat 
# And filter by BH adjusted p-value < 0.05
# UCLA Schizophrenia
UCLA_SPI_combo_main_repeats <- UCLA_SPI_combo_main_full %>%
  group_by(grouping_var, Sample_Type, Noise_Proc, repeat_number) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., UCLA_SPI_combo_pvals) %>%
  ungroup() %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))

UCLA_SPI_combo_main_repeats_sig <- UCLA_SPI_combo_main_repeats %>%
  filter(bal_acc_p_adj < 0.05)

# ABIDE ASD
ABIDE_SPI_combo_main_repeats <- ABIDE_SPI_combo_main_full %>%
  group_by(grouping_var, Sample_Type, Noise_Proc, repeat_number) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., ABIDE_SPI_combo_pvals) %>%
  ungroup() %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))

ABIDE_SPI_combo_main_repeats_sig <- ABIDE_SPI_combo_main_repeats %>%
  filter(bal_acc_p_adj < 0.05)

### Aggregate all combo univariate pairwise SPI-wise data across repeats
# UCLA Schizophrenia
UCLA_SPI_combo_main <- UCLA_SPI_combo_main_repeats %>%
  group_by(grouping_var, Sample_Type, Noise_Proc) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T),
            mean_SD = sd(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., UCLA_SPI_combo_pvals) %>%
  ungroup() %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))

UCLA_SPI_combo_main_sig <- UCLA_SPI_combo_main %>%
  filter(bal_acc_p_adj < 0.05)

# ABIDE ASD
ABIDE_SPI_combo_main <- ABIDE_SPI_combo_main_repeats %>%
  group_by(grouping_var, Sample_Type, Noise_Proc) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T),
            mean_SD = sd(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., ABIDE_SPI_combo_pvals) %>%
  ungroup() %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))

ABIDE_SPI_combo_main_sig <- ABIDE_SPI_combo_main %>%
  filter(bal_acc_p_adj < 0.05)

### Violin plots that show spaghetti lines colored by whether difference is significant
# Find which SPIs are statistically different with vs without univariate
# combo-wise information
# UCLA Schizophrenia
plot_SPI_with_without_univar(dataset_ID,
                             SPI_data_main_repeats = UCLA_SPI_main_repeats,
                             SPI_combo_data_main_repeats = UCLA_SPI_combo_main_repeats,
                             SPI_only_color = "chartreuse3",
                             SPI_univar_color = "darkgoldenrod2")

ggsave(paste0(plot_path, "Figure3/UCLA_SPI_With_vs_Without_Univariate.png"),
       width = 7, height = 3, units = "in", dpi = 300)

# ABIDE ASD
plot_SPI_with_without_univar(dataset_ID,
                             SPI_data_main_repeats = ABIDE_SPI_main_repeats,
                             SPI_combo_data_main_repeats = ABIDE_SPI_combo_main_repeats,
                             SPI_only_color = "chartreuse3",
                             SPI_univar_color = "darkgoldenrod2")

ggsave(paste0(plot_path, "Figure3/ABIDE_SPI_With_vs_Without_Univariate.png"),
       width = 7, height = 3, units = "in", dpi = 300)