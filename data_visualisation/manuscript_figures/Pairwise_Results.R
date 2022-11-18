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
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/Figure3/")
icesTAF::mkdir(plot_path)

SCZ_data_path <- "~/data/UCLA_Schizophrenia/"
SCZ_rdata_path <- paste0(SCZ_data_path, "processed_data/Rdata/")
ASD_data_path <- "~/data/ABIDE_ASD/"
ASD_rdata_path <- paste0(ASD_data_path, "processed_data/Rdata/")

noise_proc_SCZ <- "AROMA+2P+GMR"
noise_proc_ASD <- "FC1000"

raw_SCZ_pyspi_res <- readRDS(paste0(SCZ_rdata_path, "UCLA_Schizophrenia_pyspi14_filtered.Rds"))
raw_ASD_pyspi_res <- readRDS(paste0(ASD_rdata_path, "ABIDE_ASD_pyspi14_filtered.Rds"))

################################################################################
# Figure 3A pairwise SPI-wise analysis results
################################################################################

raw_SCZ_pyspi_res %>%
  filter(brain_region_1 != brain_region_2) %>%
  filter(is.na(value) | is.nan(value)) %>%
  mutate(SPI = str_replace_all(SPI, "_", " ")) %>%
  group_by(SPI, brain_region_1) %>%
  summarise(num_NA = n()) %>%
  ggplot(data=., mapping=aes(x = brain_region_1,
                             y = SPI,
                             fill = num_NA)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="red") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.4),
        legend.position = "none") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 15))

raw_SCZ_pyspi_res %>%
  filter(brain_region_1 != brain_region_2) %>%
  filter(is.na(value) | is.nan(value)) %>%
  mutate(SPI = str_replace_all(SPI, "_", " ")) %>%
  group_by(SPI, Sample_ID) %>%
  summarise(num_NA = n()) %>%
  ggplot(data=., mapping=aes(x = Sample_ID,
                             y = SPI,
                             fill = num_NA)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="red") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.4),
        legend.position = "none") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 15))


### Load pairwise SPI-wise p-value data
# UCLA Schizophrenia
UCLA_SPI_pvals <- readRDS(paste0(UCLA_rdata_path, "SPI_wise_CV_linear_SVM_model_permutation_null_pyspi14_inv_prob_pvals.Rds")) %>%
  filter(Noise_Proc == noise_proc_SCZ) %>%
  ungroup() %>%
  mutate(bal_acc_p_adj_bonf = p.adjust(bal_acc_p, method="bonferroni")) %>%
  dplyr::select(grouping_var, bal_acc_p, bal_acc_p_adj, bal_acc_p_adj_bonf)
# ABIDE ASD
ABIDE_SPI_pvals <- readRDS(paste0(ABIDE_rdata_path, "SPI_wise_CV_linear_SVM_model_permutation_null_pyspi14_inv_prob_pvals.Rds")) %>%
  filter(Noise_Proc == noise_proc_ASD) %>%
  ungroup() %>%
  mutate(bal_acc_p_adj_bonf = p.adjust(bal_acc_p, method="bonferroni")) %>%
  dplyr::select(grouping_var, bal_acc_p, bal_acc_p_adj, bal_acc_p_adj_bonf)

### Load full pairwise SPI-wise data
# UCLA Schizophrenia
UCLA_SPI_main_full <- readRDS(paste0(UCLA_rdata_path, "SPI_wise_CV_linear_SVM_pyspi14_inv_prob_balacc.Rds"))%>%
  filter(Noise_Proc == noise_proc_SCZ) 
# ABIDE ASD
ABIDE_SPI_main_full <- readRDS(paste0(ABIDE_rdata_path, "SPI_wise_CV_linear_SVM_pyspi14_inv_prob_balacc.Rds")) %>%
  filter(Noise_Proc == noise_proc_ASD) 

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

