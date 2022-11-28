################################################################################
# Load libraries
################################################################################

library(tidyverse)
library(icesTAF)
library(plotly)
library(cowplot)
library(broom)
theme_set(theme_cowplot())

################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/Figure4/")
TAF::mkdir(plot_path)

SCZ_data_path <- "~/data/UCLA_Schizophrenia/"
SCZ_rdata_path <- paste0(SCZ_data_path, "processed_data/Rdata/")
ASD_data_path <- "~/data/ABIDE_ASD/"
ASD_rdata_path <- paste0(ASD_data_path, "processed_data/Rdata/")

noise_proc_SCZ <- "AROMA+2P+GMR"
noise_proc_ASD <- "FC1000"

univariate_feature_set <- "catch22"
pairwise_feature_set <- "pyspi14_corrected"

################################################################################
# Figure 4 violin plot comparing all balanced accuracy values
################################################################################

### Load combo univariate and pairwise SPI-wise p-value data
# UCLA Schizophrenia
SCZ_SPI_pvals <- readRDS(paste0(SCZ_rdata_path, 
                                sprintf("SPI_wise_CV_linear_SVM_model_permutation_null_%s_inv_prob_pvals.Rds",
                                        pairwise_feature_set))) %>%
  filter(Noise_Proc == noise_proc_SCZ) %>%
  ungroup() %>%
  mutate(bal_acc_p_adj_bonf = p.adjust(bal_acc_p, method="bonferroni")) %>%
  dplyr::select(grouping_var, bal_acc_p, bal_acc_p_adj, bal_acc_p_adj_bonf)

SCZ_SPI_combo_pvals <- readRDS(paste0(SCZ_rdata_path, 
                                       sprintf("univariate_%s_pairwise_%s_CV_linear_SVM_model_permutation_null_inv_prob_pvals.Rds",
                                               univariate_feature_set,
                                               pairwise_feature_set))) %>%
  filter(Noise_Proc == noise_proc_SCZ) %>%
  ungroup() %>%
  mutate(bal_acc_p_adj_bonf = p.adjust(bal_acc_p, method="bonferroni")) %>%
  dplyr::select(grouping_var, bal_acc_p, bal_acc_p_adj, bal_acc_p_adj_bonf)

# ABIDE ASD
ASD_SPI_combo_pvals <- readRDS(paste0(ASD_rdata_path,
                                        sprintf("univariate_%s_pairwise_%s_CV_linear_SVM_model_permutation_null_inv_prob_pvals.Rds",
                                                univariate_feature_set,
                                                pairwise_feature_set))) %>%
  filter(Noise_Proc == "FC1000") %>%
  ungroup() %>%
  mutate(bal_acc_p_adj_bonf = p.adjust(bal_acc_p, method="bonferroni")) %>%
  dplyr::select(grouping_var, bal_acc_p, bal_acc_p_adj, bal_acc_p_adj_bonf)

ASD_SPI_pvals <- readRDS(paste0(ASD_rdata_path, 
                                sprintf("SPI_wise_CV_linear_SVM_model_permutation_null_%s_inv_prob_pvals.Rds",
                                        pairwise_feature_set))) %>%
  filter(Noise_Proc == noise_proc_ASD) %>%
  ungroup() %>%
  mutate(bal_acc_p_adj_bonf = p.adjust(bal_acc_p, method="bonferroni")) %>%
  dplyr::select(grouping_var, bal_acc_p, bal_acc_p_adj, bal_acc_p_adj_bonf)

### Load full combo univariate and pairwise SPI-wise data
# UCLA Schizophrenia
SCZ_SPI_combo_main_full <- readRDS(paste0(SCZ_rdata_path, sprintf("univariate_%s_pairwise_%s_CV_linear_SVM_inv_prob_balacc.Rds",
                                                                  univariate_feature_set,
                                                                  pairwise_feature_set))) %>%
  filter(Noise_Proc == noise_proc_SCZ) %>%
  dplyr::rename("grouping_var" = "SPI")

SCZ_SPI_main_full <- readRDS(paste0(SCZ_rdata_path, 
                                    sprintf("SPI_wise_CV_linear_SVM_%s_inv_prob_balacc.Rds",
                                            pairwise_feature_set))) %>%
  filter(Noise_Proc == noise_proc_SCZ) 

# ABIDE ASD
ASD_SPI_combo_main_full <- readRDS(paste0(ASD_rdata_path, sprintf("univariate_%s_pairwise_%s_CV_linear_SVM_inv_prob_balacc.Rds",
                                                                  univariate_feature_set,
                                                                  pairwise_feature_set))) %>%
  filter(Noise_Proc == "FC1000") %>%
  dplyr::rename("grouping_var" = "SPI")

ASD_SPI_main_full <- readRDS(paste0(ASD_rdata_path, 
                                    sprintf("SPI_wise_CV_linear_SVM_%s_inv_prob_balacc.Rds",
                                            pairwise_feature_set))) %>%
  filter(Noise_Proc == noise_proc_ASD) 

### Aggregate combo univariate and pairwise SPI-wise data per repeat 
# And filter by BH adjusted p-value < 0.05
# UCLA Schizophrenia
SCZ_SPI_combo_main_repeats <- SCZ_SPI_combo_main_full %>%
  group_by(grouping_var, Sample_Type, Noise_Proc, repeat_number) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., SCZ_SPI_combo_pvals) %>%
  ungroup() %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))

SCZ_SPI_combo_main_repeats_sig <- SCZ_SPI_combo_main_repeats %>%
  filter(bal_acc_p_adj < 0.05)

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
ASD_SPI_combo_main_repeats <- ASD_SPI_combo_main_full %>%
  group_by(grouping_var, Sample_Type, Noise_Proc, repeat_number) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., ASD_SPI_combo_pvals) %>%
  ungroup() %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))

ASD_SPI_combo_main_repeats_sig <- ASD_SPI_combo_main_repeats %>%
  filter(bal_acc_p_adj < 0.05)

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


### Aggregate all combo univariate pairwise SPI-wise data across repeats
# UCLA Schizophrenia
SCZ_SPI_combo_main <- SCZ_SPI_combo_main_repeats %>%
  group_by(grouping_var, Sample_Type, Noise_Proc) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T),
            mean_SD = sd(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., SCZ_SPI_combo_pvals) %>%
  ungroup() %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))

SCZ_SPI_combo_main_sig <- SCZ_SPI_combo_main %>%
  filter(bal_acc_p_adj < 0.05)

# ABIDE ASD
ASD_SPI_combo_main <- ASD_SPI_combo_main_repeats %>%
  group_by(grouping_var, Sample_Type, Noise_Proc) %>%
  summarise(mean_balanced_accuracy = mean(balanced_accuracy, na.rm=T),
            mean_SD = sd(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balanced_accuracy") %>%
  left_join(., ASD_SPI_combo_pvals) %>%
  ungroup() %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))

ASD_SPI_combo_main_sig <- ASD_SPI_combo_main %>%
  filter(bal_acc_p_adj < 0.05)

### Load pairwise SPI-wise empirical null distributions
# UCLA Schizophrenia
SCZ_SPI_combo_null <- readRDS(paste0(SCZ_rdata_path, 
                               sprintf("UCLA_Schizophrenia_univariate_%s_pairwise_%s_model_permutation_null_inv_prob.Rds",
                                       univariate_feature_set,
                                       pairwise_feature_set)))
# ABIDE ASD
ASD_SPI_combo_null <- readRDS(paste0(ASD_rdata_path, 
                                     sprintf("ABIDE_ASD_univariate_%s_pairwise_%s_model_permutation_null_inv_prob.Rds",
                                             univariate_feature_set,
                                             pairwise_feature_set)))

### Boxplot for SPI+combo significant features
# UCLA Schizophrenia
plot_boxplot_shaded_null(dataset_ID = "UCLA_Schizophrenia",
                         grouping_var_name = "SPI +\nUnivariate Combo",
                         main_data_by_repeat = SCZ_SPI_combo_main_repeats_sig,
                         fill_color = "darkgoldenrod2",
                         wrap_length = 100,
                         null_mean_value = mean(SCZ_SPI_combo_null$balanced_accuracy, na.rm=T),
                         null_SD_value = sd(SCZ_SPI_combo_null$balanced_accuracy, na.rm=T))
ggsave(paste0(plot_path, "UCLA_Schizophrenia_SPI_combo_sig_boxplot.png"),
       width = 7, height = 3, units="in", dpi=300)

# ABIDE ASD
plot_boxplot_shaded_null(dataset_ID = "ABIDE_ASD",
                         grouping_var_name = "SPI +\nUnivariate Combo",
                         main_data_by_repeat = SCZ_SPI_combo_main_repeats_sig,
                         fill_color = "darkgoldenrod2",
                         wrap_length = 100,
                         null_mean_value = mean(SCZ_SPI_combo_null$balanced_accuracy, na.rm=T),
                         null_SD_value = sd(SCZ_SPI_combo_null$balanced_accuracy, na.rm=T))
ggsave(paste0(plot_path, "ABIDE_ASD_SPI_combo_sig_boxplot.png"),
       width = 7, height = 3, units="in", dpi=300)


### corrected re-sampled t-test statistic
# currently just using standard student's t-test, need to update
merged_SCZ_SPI_combo_data <- SCZ_SPI_main_repeats %>% 
  mutate(Analysis = "SPI") %>%
  plyr::rbind.fill(., SCZ_SPI_combo_main_repeats %>% 
                     mutate(Analysis = "SPI+combo")) %>%
  mutate(Analysis = factor(Analysis, levels = c("SPI",
                                                "SPI+combo")))

# Find SPIs with a significant difference with vs without univariate combo info
# Using paired two-way Student's t-test
t_test_res_SCZ_SPIs <- merged_SCZ_SPI_combo_data %>%
  nest(data = -grouping_var) %>%
  mutate(test = map(data, ~ t.test(balanced_accuracy ~ Analysis, data = .x, paired=T)),
         tidied = map(test, broom::tidy)) %>%
  unnest(tidied) %>%
  ungroup() %>%
  mutate(t_test_p_adj = p.adjust(p.value, method="BH")) %>%
  dplyr::select(grouping_var, statistic, t_test_p_adj)

merged_SCZ_SPI_combo_data %>%
  left_join(., t_test_res_SCZ_SPIs) %>%
  mutate(grouping_var_char = gsub("_", " ", as.character(grouping_var))) %>%
  mutate(grouping_var = factor(grouping_var, levels = rev(levels(grouping_var)))) %>%
  mutate(Analysis = as.character(Analysis)) %>%
  rowwise() %>%
  mutate(tofill_violin = ifelse(bal_acc_p_adj < 0.05, Analysis, "no"),
         tofill_bg = ifelse(t_test_p_adj < 0.05, "yes", "no")) %>%
  ggplot(data=., mapping=aes(x = grouping_var,
                             y = balanced_accuracy)) +
  geom_rect(aes(fill = tofill_bg), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.3) +
  geom_violin(aes(fill = tofill_violin))  +
  facet_wrap(grouping_var ~ ., scales="free_x", nrow=1) +
  scale_fill_manual(values = c("gray90", "chartreuse3", "darkgoldenrod2", "gray90")) +
  ylab("Balanced Accuracy") +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow=2)) +
  theme(axis.text.x = element_text(angle=90),
        strip.text = element_blank())


### Violin plots that show spaghetti lines colored by whether difference is significant
# Find which SPIs are statistically different with vs without univariate
# combo-wise information
# UCLA Schizophrenia
plot_SPI_with_without_univar(dataset_ID = "UCLA_Schizophrenia",
                             SPI_data_main_repeats = SCZ_SPI_main_repeats,
                             SPI_combo_data_main_repeats = SCZ_SPI_combo_main_repeats,
                             SPI_only_color = "chartreuse3",
                             SPI_univar_color = "darkgoldenrod2")

ggsave(paste0(plot_path, "UCLA_Schizophrenia_SPI_With_vs_Without_Univariate.png"),
       width = 7, height = 3, units = "in", dpi = 300)

# ABIDE ASD
plot_SPI_with_without_univar(dataset_ID = "ABIDE_ASD",
                             SPI_data_main_repeats = ASD_SPI_main_repeats,
                             SPI_combo_data_main_repeats = ASD_SPI_combo_main_repeats,
                             SPI_only_color = "chartreuse3",
                             SPI_univar_color = "darkgoldenrod2")

ggsave(paste0(plot_path, "ABIDE_ASD_SPI_With_vs_Without_Univariate.png"),
       width = 7, height = 3, units = "in", dpi = 300)