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
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/Figure4/")
icesTAF::mkdir(plot_path)

SCZ_data_path <- "~/data/UCLA_Schizophrenia/"
SCZ_rdata_path <- paste0(SCZ_data_path, "processed_data/Rdata/")
ASD_data_path <- "~/data/ABIDE_ASD/"
ASD_rdata_path <- paste0(ASD_data_path, "processed_data/Rdata/")


################################################################################
# Figure 4 violin plot comparing all balanced accuracy values
################################################################################

### Load combo univariate and pairwise SPI-wise p-value data
# UCLA Schizophrenia
UCLA_SPI_combo_pvals <- readRDS(paste0(UCLA_rdata_path, 
                                       "univariate_catch22_pairwise_pyspi14_CV_linear_SVM_model_permutation_null_inv_prob_pvals.Rds")) %>%
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