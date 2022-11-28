################################################################################
# Load libraries
################################################################################

library(tidyverse)
library(icesTAF)
library(plotly)
library(cowplot)
theme_set(theme_cowplot())

################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
source(paste0(github_dir, "data_visualisation/manuscript_figures/Manuscript_Draft_Visualisations_Helper.R"))
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/Figure3/")
TAF::mkdir(plot_path)

SCZ_data_path <- "~/data/UCLA_Schizophrenia/"
SCZ_rdata_path <- paste0(SCZ_data_path, "processed_data/Rdata/")
ASD_data_path <- "~/data/ABIDE_ASD/"
ASD_rdata_path <- paste0(ASD_data_path, "processed_data/Rdata/")

noise_proc_SCZ <- "AROMA+2P+GMR"
noise_proc_ASD <- "FC1000"

pairwise_feature_set <- "pyspi14_corrected"

################################################################################
# Figure 3A pairwise SPI-wise analysis results
################################################################################

### Load subject-wise predictions
SCZ_SPI_subject_preds <- readRDS(paste0(SCZ_rdata_path,
                                        sprintf("SPI_wise_CV_linear_SVM_%s_inv_prob.Rds",
                                                pairwise_feature_set)))
ASD_SPI_subject_preds <- readRDS(paste0(ASD_rdata_path,
                                        sprintf("SPI_wise_CV_linear_SVM_%s_inv_prob.Rds",
                                                pairwise_feature_set)))

### Load pairwise SPI-wise p-value data
# UCLA Schizophrenia
SCZ_SPI_pvals <- readRDS(paste0(SCZ_rdata_path, 
                                sprintf("SPI_wise_CV_linear_SVM_model_permutation_null_%s_inv_prob_pvals.Rds",
                                        pairwise_feature_set))) %>%
  filter(Noise_Proc == noise_proc_SCZ) %>%
  ungroup() %>%
  mutate(bal_acc_p_adj_bonf = p.adjust(bal_acc_p, method="bonferroni")) %>%
  dplyr::select(grouping_var, bal_acc_p, bal_acc_p_adj, bal_acc_p_adj_bonf)

# ABIDE ASD
ASD_SPI_pvals <- readRDS(paste0(ASD_rdata_path, 
                                sprintf("SPI_wise_CV_linear_SVM_model_permutation_null_%s_inv_prob_pvals.Rds",
                                        pairwise_feature_set))) %>%
  filter(Noise_Proc == noise_proc_ASD) %>%
  ungroup() %>%
  mutate(bal_acc_p_adj_bonf = p.adjust(bal_acc_p, method="bonferroni")) %>%
  dplyr::select(grouping_var, bal_acc_p, bal_acc_p_adj, bal_acc_p_adj_bonf)

### Load full pairwise SPI-wise data
# UCLA Schizophrenia
SCZ_SPI_main_full <- readRDS(paste0(SCZ_rdata_path, 
                                    sprintf("SPI_wise_CV_linear_SVM_%s_inv_prob_balacc.Rds",
                                            pairwise_feature_set))) %>%
  filter(Noise_Proc == noise_proc_SCZ) 

# ABIDE ASD
ASD_SPI_main_full <- readRDS(paste0(ASD_rdata_path, 
                                      sprintf("SPI_wise_CV_linear_SVM_%s_inv_prob_balacc.Rds",
                                              pairwise_feature_set))) %>%
  filter(Noise_Proc == noise_proc_ASD) 

### Aggregate pairwise SPI-wise data per repeat 
# And filter by BH adjusted p-value < 0.05
# UCLA Schizophrenia
SCZ_SPI_main_repeats <- SCZ_SPI_main_full %>%
  group_by(grouping_var, Sample_Type, Noise_Proc, repeat_number) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., SCZ_SPI_pvals) %>%
  ungroup() %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))

SCZ_SPI_main_repeats_sig <- SCZ_SPI_main_repeats %>%
  filter(bal_acc_p_adj < 0.05)

# ABIDE ASD
ASD_SPI_main_repeats <- ASD_SPI_main_full %>%
  group_by(grouping_var, Sample_Type, Noise_Proc, repeat_number) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., ASD_SPI_pvals) %>%
  ungroup() %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))
ASD_SPI_main_repeats_sig <- ASD_SPI_main_repeats %>%
  filter(bal_acc_p_adj < 0.05)

### Aggregate all pairwise SPI-wise data across repeats
# UCLA Schizophrenia
SCZ_SPI_main <- SCZ_SPI_main_repeats %>%
  group_by(grouping_var, Sample_Type, Noise_Proc) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T),
            mean_SD = sd(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., SCZ_SPI_pvals) %>%
  ungroup() %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))

SCZ_SPI_main_sig <- SCZ_SPI_main %>%
  filter(bal_acc_p_adj < 0.05)

# ABIDE ASD
ASD_SPI_main <- ASD_SPI_main_repeats %>%
  group_by(grouping_var, Sample_Type, Noise_Proc) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T),
            mean_SD = sd(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., ASD_SPI_pvals) %>%
  ungroup() %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))

ASD_SPI_main_sig <- ASD_SPI_main %>%
  filter(bal_acc_p_adj < 0.05)

### Load pairwise SPI-wise empirical null distributions
# UCLA Schizophrenia
SCZ_SPI_null <- readRDS(paste0(SCZ_rdata_path, 
                               sprintf("UCLA_Schizophrenia_SPI_wise_model_permutation_null_%s_inv_prob.Rds",
                                       pairwise_feature_set)))
# ABIDE ASD
ASD_SPI_null <- readRDS(paste0(ASD_rdata_path, 
                               sprintf("ABIDE_ASD_SPI_wise_model_permutation_null_%s_inv_prob.Rds",
                                       pairwise_feature_set)))

### SPI boxplot with shaded null region
# UCLA Schizophrenia
plot_boxplot_shaded_null(dataset_ID = "UCLA_Schizophrenia",
                         grouping_var_name = "Pairwise SPI",
                         main_data_by_repeat = SCZ_SPI_main_repeats_sig,
                         fill_color = "chartreuse3",
                         wrap_length = 100,
                         null_mean_value = mean(SCZ_SPI_null$balanced_accuracy, na.rm=T),
                         null_SD_value = sd(SCZ_SPI_null$balanced_accuracy, na.rm=T))
ggsave(paste0(plot_path, "UCLA_Schizophrenia_SPI_sig_boxplot.png"),
       width = 8, height = 2.5, units="in", dpi=300)
# ABIDE ASD
plot_boxplot_shaded_null(dataset_ID = "ABIDE_ASD",
                         grouping_var_name = "Pairwise SPI",
                         main_data_by_repeat = ASD_SPI_main_repeats_sig,
                         fill_color = "chartreuse3",
                         wrap_length = 100,
                         null_mean_value = mean(ASD_SPI_null$balanced_accuracy, na.rm=T),
                         null_SD_value = sd(ASD_SPI_null$balanced_accuracy, na.rm=T))
ggsave(paste0(plot_path, "ABIDE_ASD_SPI_sig_boxplot.png"),
       width = 8, height = 2.5, units="in", dpi=300)


