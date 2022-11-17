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
# Figure 2: Univariate results
################################################################################
icesTAF::mkdir(paste0(plot_path, "Figure2/"))

# Load UCLA p-values
UCLA_ROI_pvals <- readRDS(paste0(UCLA_rdata_path, "ROI_wise_CV_linear_SVM_model_permutation_null_catch22_inv_prob_pvals.Rds")) %>%
  filter(Noise_Proc == "AROMA+2P+GMR") %>%
  ungroup() %>%
  mutate(bal_acc_p_adj_bonf = p.adjust(bal_acc_p, method="bonferroni")) %>%
  dplyr::select(grouping_var, bal_acc_p, bal_acc_p_adj, bal_acc_p_adj_bonf)

# Load UCLA main data
UCLA_ROI_main <- readRDS(paste0(UCLA_rdata_path, "ROI_wise_CV_linear_SVM_catch22_inv_prob_balacc.Rds")) %>%
  filter(Noise_Proc == "AROMA+2P+GMR") %>%
  group_by(grouping_var, repeat_number) %>%
  summarise(mean_balacc = mean(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balacc") %>%
  left_join(., UCLA_ROI_pvals) %>%
  ungroup() %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))

# Load UCLA null data
UCLA_ROI_null <- readRDS(paste0(UCLA_rdata_path, "UCLA_Schizophrenia_ROI_wise_model_permutation_null_catch22_inv_prob.Rds"))

# Load ABIDE p-values
ABIDE_ROI_pvals <- readRDS(paste0(ABIDE_rdata_path, "ROI_wise_CV_linear_SVM_model_permutation_null_catch22_inv_prob_pvals.Rds")) %>%
  filter(Noise_Proc == "FC1000") %>%
  ungroup() %>%
  mutate(bal_acc_p_adj_bonf = p.adjust(bal_acc_p, method="bonferroni")) %>%
  dplyr::select(grouping_var, bal_acc_p, bal_acc_p_adj, bal_acc_p_adj_bonf)

# Load ABIDE main data
ABIDE_ROI_main <- readRDS(paste0(ABIDE_rdata_path, "ROI_wise_CV_linear_SVM_catch22_inv_prob_balacc.Rds")) %>%
  filter(Noise_Proc == "FC1000") %>%
  group_by(grouping_var, repeat_number) %>%
  summarise(mean_balacc = mean(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balacc") %>%
  left_join(., ABIDE_ROI_pvals)

# Load ABIDE null data
ABIDE_ROI_null <- readRDS(paste0(ABIDE_rdata_path, "ABIDE_ASD_ROI_wise_model_permutation_null_catch22_inv_prob.Rds"))

# UCLA boxplot with shaded null region
UCLA_null_mean <- mean(UCLA_ROI_null$balanced_accuracy)
UCLA_null_SD <- sd(UCLA_ROI_null$balanced_accuracy)
ggplot() +
  geom_boxplot(data = UCLA_ROI_main %>% filter(bal_acc_p_adj < 0.05) %>%
                 mutate(grouping_var = str_replace_all(grouping_var, "ctx-lh-|Left-", "Left ")) %>%
                 mutate(grouping_var = str_replace_all(grouping_var, "ctx-rh-|Right-", "Right ")) %>%
                 mutate(grouping_var = fct_reorder(grouping_var, 
                                                   balanced_accuracy,
                                                   .fun = mean)), 
               aes(y=grouping_var, x=balanced_accuracy),
               fill = "#F0224B", color="white") +
  geom_rect(data = UCLA_ROI_main, 
            xmin = UCLA_null_mean - UCLA_null_SD,
            xmax = UCLA_null_mean + UCLA_null_SD,
            ymin=0, ymax=Inf, fill="gray85") +
  geom_vline(xintercept = UCLA_null_mean) +
  scale_x_continuous(limits = c(0.95*(UCLA_null_mean - UCLA_null_SD),
                                1.05*max(UCLA_ROI_main %>% 
                                          filter(bal_acc_p_adj < 0.05) %>% 
                                          pull(balanced_accuracy)))) +
  xlab("Balanced Accuracy") +
  ylab("Brain Region")

ggsave(paste0(plot_path, "Figure2/UCLA_Schizophrenia_Brain_Region_sig_boxplot.png"),
       width = 3.7, height = 2, units="in", dpi=300)

# ABIDE boxplot with shaded null region
ABIDE_null_mean <- mean(ABIDE_ROI_null$balanced_accuracy)
ABIDE_null_SD <- sd(ABIDE_ROI_null$balanced_accuracy)
ggplot() +
  geom_boxplot(data = ABIDE_ROI_main %>% filter(bal_acc_p_adj < 0.05) %>%
                 mutate(grouping_var = str_replace_all(grouping_var, "ctx-lh-|Left-", "Left ")) %>%
                 mutate(grouping_var = str_replace_all(grouping_var, "ctx-rh-|Right-", "Right ")) %>%
                 mutate(grouping_var = fct_reorder(grouping_var, 
                                                   balanced_accuracy,
                                                   .fun = mean)), 
               aes(y=grouping_var, x=balanced_accuracy),
               fill = "#F0224B", color="white") +
  geom_rect(data = ABIDE_ROI_main, 
            xmin = ABIDE_null_mean - ABIDE_null_SD,
            xmax = ABIDE_null_mean + ABIDE_null_SD,
            ymin=0, ymax=Inf, fill="gray85") +
  geom_vline(xintercept = ABIDE_null_mean) +
  scale_x_continuous(limits = c(0.95*(ABIDE_null_mean - ABIDE_null_SD),
                                1.05*max(ABIDE_ROI_main %>% 
                                           filter(bal_acc_p_adj < 0.05) %>% 
                                           pull(balanced_accuracy)))) +
  xlab("Balanced Accuracy") +
  ylab("Brain Region") +
  theme(axis.line = element_line(color="white"),
        axis.text = element_text(color="white"),
        axis.title = element_text(color="white"),
        axis.ticks = element_line(color="white")) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 20))

ggsave(paste0(plot_path, "ABIDE_ASD_Brain_Region_sig_boxplot.png"),
       width = 3.7, height = 1.6, units="in", dpi=300)

# Univariate region-wise
UCLA_ROI_main$Type = "Main"
UCLA_ROI_null$Type = "Null"
ABIDE_ROI_main$Type = "Main"
ABIDE_ROI_null$Type = "Null"

# UCLA ggseg to highlight significant brain regions
UCLA_significant_brain_regions <- UCLA_ROI_main %>%
  filter(bal_acc_p_adj < 0.05) %>%
  dplyr::select(grouping_var) %>%
  dplyr::rename("label" = "grouping_var") %>%
  distinct() %>%
  mutate(fillyes = T) %>%
  mutate(label = ifelse(str_detect(label, "ctx-"),
                        gsub("-", "_", label),
                        label)) %>%
  mutate(label = gsub("ctx_", "", label)) %>%
  left_join(., dk %>% as_tibble()) %>%
  filter(!is.na(region)) %>%
  ungroup() %>%
  dplyr::select(-label)

UCLA_significant_brain_regions %>%
  ggseg(atlas = "dk", mapping = aes(fill = fillyes),
        position = "dispersed", colour = "darkgrey") +
  scale_fill_manual(values = c("#F0224B"), na.value = NA) +
  labs(fill = "Balanced Accuracy") +
  theme_void() +
  theme(plot.title = element_blank(),
        legend.position = "none") 
ggsave(paste0(plot_path, "UCLA_Schizophrenia_ROI_wise_cortex_balanced_accuracy.png"),
       width = 4, height = 1.5, units = "in", dpi = 300)

# UCLA Subcortical data
aseg_lookup <- aseg$data

UCLA_significant_subcortex <- UCLA_ROI_main %>%
  filter(bal_acc_p_adj < 0.05) %>%
  dplyr::select(grouping_var) %>%
  dplyr::rename("label" = "grouping_var") %>%
  distinct() %>%
  mutate(fillyes = T) %>%
  mutate(label = ifelse(str_detect(label, "ctx-"),
                        gsub("-", "_", label),
                        label)) %>%
  mutate(label = gsub("ctx_", "", label)) %>%
  left_join(., aseg %>% as_tibble()) %>%
  filter(!is.na(region)) %>%
  ungroup() %>%
  dplyr::select(-label)

UCLA_significant_subcortex %>%
  ggseg(atlas="aseg", mapping=aes(fill = fillyes),
        position="stacked", colour="darkgrey", side = "coronal") +
  scale_fill_manual(values = c("#F0224B"), na.value = NA) +
  labs(fill = "Balanced Accuracy") +
  theme_cowplot() +
  theme_void() +
  theme(plot.title = element_blank(),
        legend.position = "none") 
ggsave(paste0(plot_path, "UCLA_ROI_wise_subcortex_balanced_accuracy.png"),
       width = 4, height = 4, units = "in", dpi = 300)

# ABIDE
ABIDE_significant_brain_regions <- ABIDE_ROI_main %>%
  filter(bal_acc_p_adj < 0.05) %>%
  dplyr::select(grouping_var) %>%
  dplyr::rename("Brain_Region" = "grouping_var") %>%
  distinct() %>%
  left_join(., ABIDE_brain_region_info) %>% 
  dplyr::rename("region" = "ggseg") %>%
  mutate(fillyes = T) %>%
  left_join(., hoCort %>% as_tibble()) %>%
  filter(!is.na(region)) 

ABIDE_significant_brain_regions %>%
  ggseg(atlas = "hoCort", mapping = aes(fill = fillyes),
        position = "dispersed", colour = "darkgrey") +
  scale_fill_manual(values = c("#F0224B"), na.value = NA) +
  labs(fill = "Balanced Accuracy") +
  theme_void() +
  theme(plot.title = element_blank(),
        legend.position = "none") 
ggsave(paste0(plot_path, "ABIDE_ASD_ROI_wise_cortex_balanced_accuracy.png"),
       width = 4, height = 1.5, units = "in", dpi = 300)

################################################################################
# Slide 7 Combo-wise analysis
################################################################################

# Load UCLA p-values
UCLA_combo_pvals <- readRDS(paste0(UCLA_rdata_path, "Combo_wise_CV_linear_SVM_model_permutation_null_catch22_inv_prob_pvals.Rds")) %>%
  filter(Noise_Proc == "AROMA+2P+GMR") %>%
  ungroup() %>%
  mutate(bal_acc_p_adj_bonf = p.adjust(bal_acc_p, method="bonferroni")) %>%
  dplyr::select(grouping_var, bal_acc_p, bal_acc_p_adj, bal_acc_p_adj_bonf)

# Load UCLA main data
UCLA_combo_main <- readRDS(paste0(UCLA_rdata_path, "Combo_wise_CV_linear_SVM_catch22_inv_prob_balacc.Rds")) %>%
  filter(Noise_Proc == "AROMA+2P+GMR") %>%
  group_by(grouping_var, repeat_number) %>%
  summarise(mean_balacc = mean(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balacc") %>%
  left_join(., UCLA_combo_pvals)

# Load UCLA null data
UCLA_combo_null <- readRDS(paste0(UCLA_rdata_path, "UCLA_Schizophrenia_Combo_wise_model_permutation_null_catch22_inv_prob.Rds"))

# UCLA Violin plot comparison for univariate results
UCLA_uni_region <- UCLA_ROI_main %>%
  mutate(Analysis = "Region-wise",
         Variable_Type = "Univariate") %>%
  group_by(grouping_var, Analysis, Variable_Type) %>%
  summarise(mean_balacc = mean(balanced_accuracy)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balacc")

UCLA_uni_combo <- UCLA_combo_main %>%
  mutate(Analysis = "Combo-wise",
         Variable_Type = "Univariate") %>%
  group_by(grouping_var, Analysis, Variable_Type) %>%
  summarise(mean_balacc = mean(balanced_accuracy))  %>%
  dplyr::rename("balanced_accuracy" = "mean_balacc")

UCLA_merged_data <- do.call(plyr::rbind.fill, list(UCLA_uni_region, 
                                                   UCLA_uni_combo)) %>%
  mutate( Analysis = factor(Analysis, levels = c("Region-wise", 
                                                 "Combo-wise"))) 

UCLA_merged_data %>%
  ggplot(data=., mapping = aes(x = Analysis, y = balanced_accuracy)) +
  geom_violin(aes(fill = Analysis)) +
  stat_summary(data = subset(UCLA_merged_data, Analysis == "Combo-wise" & Variable_Type == "Univariate"),
               geom = "crossbar", fun = "mean", aes(color=Analysis), size=1,
               color = "#9B51B4") +
  geom_boxplot(fill=NA, width=0.1, color="black") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
  scale_fill_manual(values = c("#F0224B")) +
  ylab("Balanced Accuracy") +
  theme(legend.position = "none",
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.y.left = element_blank(),
        axis.text.y = element_text(color="white"),
        axis.text.x = element_text(angle=45, hjust=1, color="white"),
        axis.title.x = element_blank(),
        axis.line = element_line(color="white"),
        axis.title.y = element_text(color="white"),
        axis.ticks = element_line(color="white"))

ggsave(paste0(plot_path, "UCLA_Univariate_Comparison_Boxplot.png"),
       width = 3, height = 2.75, units = "in", dpi = 300)

# Load ABIDE p-values
ABIDE_combo_pvals <- readRDS(paste0(ABIDE_rdata_path, "Combo_wise_CV_linear_SVM_model_permutation_null_catch22_inv_prob_pvals.Rds")) %>%
  filter(Noise_Proc == "FC1000") %>%
  ungroup() %>%
  mutate(bal_acc_p_adj_bonf = p.adjust(bal_acc_p, method="bonferroni")) %>%
  dplyr::select(grouping_var, bal_acc_p, bal_acc_p_adj, bal_acc_p_adj_bonf)

# Load ABIDE main data
ABIDE_combo_main <- readRDS(paste0(ABIDE_rdata_path, "Combo_wise_CV_linear_SVM_catch22_inv_prob_balacc.Rds")) %>%
  filter(Noise_Proc == "FC1000") %>%
  group_by(grouping_var, repeat_number) %>%
  summarise(mean_balacc = mean(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balacc") %>%
  left_join(., ABIDE_combo_pvals)

# Load ABIDE null data
ABIDE_combo_null <- readRDS(paste0(ABIDE_rdata_path, "ABIDE_ASD_Combo_wise_model_permutation_null_catch22_inv_prob.Rds"))

# ABIDE Violin plot comparison for univariate results
ABIDE_uni_region <- ABIDE_ROI_main %>%
  mutate(Analysis = "Region-wise",
         Variable_Type = "Univariate") %>%
  group_by(grouping_var, Analysis, Variable_Type) %>%
  summarise(mean_balacc = mean(balanced_accuracy))  %>%
  dplyr::rename("balanced_accuracy" = "mean_balacc")

ABIDE_uni_combo <- ABIDE_combo_main %>%
  mutate(Analysis = "Combo-wise",
         Variable_Type = "Univariate") %>%
  group_by(grouping_var, Analysis, Variable_Type) %>%
  summarise(mean_balacc = mean(balanced_accuracy))  %>%
  dplyr::rename("balanced_accuracy" = "mean_balacc")

ABIDE_merged_data <- do.call(plyr::rbind.fill, list(ABIDE_uni_region, 
                                                   ABIDE_uni_combo)) %>%
  mutate( Analysis = factor(Analysis, levels = c("Region-wise", 
                                                 "Combo-wise"))) 

ABIDE_merged_data %>%
  ggplot(data=., mapping = aes(x = Analysis, y = balanced_accuracy)) +
  geom_violin(aes(fill = Analysis)) +
  stat_summary(data = subset(ABIDE_merged_data, Analysis == "Combo-wise" & Variable_Type == "Univariate"),
               geom = "crossbar", fun = "mean", aes(color=Analysis), size=1,
               color = "#9B51B4") +
  geom_boxplot(fill=NA, width=0.1, color="black") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
  scale_fill_manual(values = c("#F0224B")) +
  ylab("Balanced Accuracy") +
  theme(legend.position = "none",
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.y.left = element_blank(),
        axis.text.y = element_text(color="white"),
        axis.text.x = element_text(angle=45, hjust=1, color="white"),
        axis.title.x = element_blank(),
        axis.line = element_line(color="white"),
        axis.title.y = element_text(color="white"),
        axis.ticks = element_line(color="white"))

ggsave(paste0(plot_path, "ABIDE_Univariate_Comparison_Boxplot.png"),
       width = 3, height = 2.75, units = "in", dpi = 300)


################################################################################
# Slide 8 Pairwise SPI-wise analysis
################################################################################

# Load UCLA p-values
UCLA_SPI_pvals <- readRDS(paste0(UCLA_rdata_path, "SPI_wise_CV_linear_SVM_model_permutation_null_pyspi14_inv_prob_pvals.Rds")) %>%
  filter(Noise_Proc == "AROMA+2P+GMR") %>%
  ungroup() %>%
  mutate(bal_acc_p_adj_bonf = p.adjust(bal_acc_p, method="bonferroni")) %>%
  dplyr::select(grouping_var, bal_acc_p, bal_acc_p_adj, bal_acc_p_adj_bonf)

# Load UCLA main data
UCLA_SPI_main <- readRDS(paste0(UCLA_rdata_path, "SPI_wise_CV_linear_SVM_pyspi14_inv_prob_balacc.Rds")) %>%
  filter(Noise_Proc == "AROMA+2P+GMR") %>%
  group_by(grouping_var, repeat_number) %>%
  summarise(mean_balacc = mean(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balacc") %>%
  left_join(., UCLA_SPI_pvals) %>%
  ungroup() %>%
  mutate(grouping_var = str_replace_all(grouping_var, "_", " ")) %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))

# Load UCLA null data
UCLA_SPI_null <- readRDS(paste0(UCLA_rdata_path, "UCLA_Schizophrenia_SPI_wise_model_permutation_null_pyspi14_inv_prob.Rds"))

# UCLA boxplot with shaded null region
UCLA_null_mean <- mean(UCLA_SPI_null$balanced_accuracy)
UCLA_null_SD <- sd(UCLA_SPI_null$balanced_accuracy)
ggplot() +
  geom_boxplot(data = UCLA_SPI_main %>% filter(bal_acc_p_adj < 0.05), 
               aes(y=grouping_var, x=balanced_accuracy),
               fill = "chartreuse3", color="white") +
  geom_rect(data = UCLA_SPI_main, 
            xmin = UCLA_null_mean - UCLA_null_SD,
            xmax = UCLA_null_mean + UCLA_null_SD,
            ymin=0, ymax=Inf, fill="gray85") +
  geom_vline(xintercept = UCLA_null_mean) +
  scale_x_continuous(limits = c(0.98*(UCLA_null_mean - UCLA_null_SD),
                                1.05*max(UCLA_SPI_main %>% 
                                           filter(bal_acc_p_adj < 0.05) %>% 
                                           pull(balanced_accuracy)))) +
  xlab("Balanced Accuracy") +
  ylab("Pairwise SPI") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
  theme(axis.line = element_line(color="white"),
        axis.text = element_text(color="white"),
        axis.title = element_text(color="white"),
        axis.ticks = element_line(color="white"))

ggsave(paste0(plot_path, "UCLA_Schizophrenia_SPI_sig_boxplot.png"),
       width = 7, height = 2.5, units="in", dpi=300)


# Load ABIDE p-values
ABIDE_SPI_pvals <- readRDS(paste0(ABIDE_rdata_path, "SPI_wise_CV_linear_SVM_model_permutation_null_pyspi14_inv_prob_pvals.Rds")) %>%
  filter(Noise_Proc == "FC1000") %>%
  ungroup() %>%
  mutate(bal_acc_p_adj_bonf = p.adjust(bal_acc_p, method="bonferroni")) %>%
  dplyr::select(grouping_var, bal_acc_p, bal_acc_p_adj, bal_acc_p_adj_bonf)

# Load ABIDE main data
ABIDE_SPI_main <- readRDS(paste0(ABIDE_rdata_path, "SPI_wise_CV_linear_SVM_pyspi14_inv_prob_balacc.Rds")) %>%
  filter(Noise_Proc == "FC1000") %>%
  group_by(grouping_var, repeat_number) %>%
  summarise(mean_balacc = mean(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balacc") %>%
  left_join(., ABIDE_SPI_pvals) %>%
  ungroup() %>%
  mutate(grouping_var = str_replace_all(grouping_var, "_", " ")) %>%
  mutate(grouping_var = fct_reorder(grouping_var, 
                                    balanced_accuracy,
                                    .fun = mean))

# Load ABIDE null data
ABIDE_SPI_null <- readRDS(paste0(ABIDE_rdata_path, "ABIDE_ASD_SPI_wise_model_permutation_null_pyspi14_inv_prob.Rds"))

# ABIDE boxplot with shaded null region
ABIDE_null_mean <- mean(ABIDE_SPI_null$balanced_accuracy)
ABIDE_null_SD <- sd(ABIDE_SPI_null$balanced_accuracy)
ggplot() +
  geom_boxplot(data = ABIDE_SPI_main %>% filter(bal_acc_p_adj < 0.05), 
               aes(y=grouping_var, x=balanced_accuracy),
               fill = "chartreuse3", color="white") +
  geom_rect(data = ABIDE_SPI_main, 
            xmin = ABIDE_null_mean - ABIDE_null_SD,
            xmax = ABIDE_null_mean + ABIDE_null_SD,
            ymin=0, ymax=Inf, fill="gray85") +
  geom_vline(xintercept = ABIDE_null_mean) +
  scale_x_continuous(limits = c(0.98*(ABIDE_null_mean - ABIDE_null_SD),
                                1.05*max(ABIDE_SPI_main %>% 
                                           filter(bal_acc_p_adj < 0.05) %>% 
                                           pull(balanced_accuracy)))) +
  xlab("Balanced Accuracy") +
  ylab("Pairwise SPI") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
  theme(axis.line = element_line(color="white"),
        axis.text = element_text(color="white"),
        axis.title = element_text(color="white"),
        axis.ticks = element_line(color="white"))

ggsave(paste0(plot_path, "ABIDE_ASD_SPI_sig_boxplot.png"),
       width = 7, height = 2.5, units="in", dpi=300)


# Stacked violin plot univariate + pairwise + combo

# Load UCLA SPI combo p-values
UCLA_SPI_combo_pvals <- readRDS(paste0(UCLA_rdata_path, "univariate_catch22_pairwise_pyspi14_CV_linear_SVM_model_permutation_null_inv_prob_pvals.Rds")) %>%
  filter(Noise_Proc == "AROMA+2P+GMR") %>%
  ungroup() %>%
  mutate(bal_acc_p_adj_bonf = p.adjust(bal_acc_p, method="bonferroni")) %>%
  dplyr::select(grouping_var, bal_acc_p, bal_acc_p_adj, bal_acc_p_adj_bonf) %>%
  dplyr::rename("SPI" = "grouping_var")

# Load SPI combo data
UCLA_SPI_combo_main <- readRDS(paste0(UCLA_rdata_path, "univariate_catch22_pairwise_pyspi14_CV_linear_SVM_inv_prob_balacc.Rds")) %>%
  filter(Noise_Proc == "AROMA+2P+GMR") %>%
  group_by(SPI, repeat_number) %>%
  summarise(mean_balacc = mean(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balacc") %>%
  left_join(., UCLA_SPI_combo_pvals) %>%
  ungroup() %>%
  mutate(SPI = str_replace_all(SPI, "_", " ")) %>%
  mutate(SPI = fct_reorder(SPI, 
                           balanced_accuracy,
                           .fun = mean))

UCLA_uni_region$Analysis <- "Univariate Region-wise"
UCLA_uni_region$Variable_Type <- "Univariate"
UCLA_uni_combo$Analysis <- "Univariate Combo-wise"
UCLA_uni_combo$Variable_Type <- "Univariate"
  
UCLA_pair_feature <- UCLA_SPI_main %>%
  group_by(grouping_var) %>%
  summarise(mean_balacc = mean(balanced_accuracy)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balacc") %>%
  mutate(Analysis = "Pairwise SPI-wise",
         Variable_Type = "Pairwise")

UCLA_pair_combo <- UCLA_SPI_combo_main %>%
  dplyr::rename("grouping_var" = "SPI") %>%
  group_by(grouping_var) %>%
  summarise(mean_balacc = mean(balanced_accuracy)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balacc") %>%
  mutate(Analysis = "Pairwise SPI + Univariate Combo-wise",
         Variable_Type = "Pairwise")

UCLA_merged_data <- do.call(plyr::rbind.fill, list(UCLA_uni_region, UCLA_uni_combo, 
                                              UCLA_pair_feature, UCLA_pair_combo)) %>%
  mutate(Variable_Type = factor(Variable_Type, levels = c("Univariate",
                                                          "Pairwise")),
         Analysis = factor(Analysis, levels = c(
           "Univariate Region-wise", 
           "Univariate Combo-wise", 
           "Pairwise SPI-wise",
           "Pairwise SPI + Univariate Combo-wise")))

UCLA_merged_data %>%
  ggplot(data=., mapping = aes(x = Analysis, y = balanced_accuracy)) +
  geom_violin(aes(fill = Analysis), color="white") +
  stat_summary(data = subset(UCLA_merged_data, Analysis == "Univariate Combo-wise"),
               geom = "crossbar", fun = "mean", aes(color=Analysis), size=1,
               color = "#9B51B4") +
  geom_boxplot(fill=NA, width=0.1, color="white") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
  geom_line(data = subset(UCLA_merged_data, Variable_Type=="Pairwise"),
            aes(x=Analysis, y=balanced_accuracy, group=grouping_var),
            color="white", alpha=0.5) +
  scale_fill_manual(values = c("#F0224B", "chartreuse3", "darkgoldenrod2")) +
  ylab("Balanced Accuracy") +
  theme(legend.position = "none",
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.y.left = element_blank(),
        axis.title.x = element_blank()) +
  theme(axis.line = element_line(color="white"),
        axis.text = element_text(color="white"),
        strip.text = element_text(color="white"),
        axis.title = element_text(color="white"),
        axis.ticks = element_line(color="white"))
ggsave(paste0(plot_dir, "UCLA_Final_Comparison_Violin.png"),
       width = 7, height = 3, units = "in", dpi = 300)

# Load ABIDE SPI combo p-values
ABIDE_SPI_combo_pvals <- readRDS(paste0(ABIDE_rdata_path, "univariate_catch22_pairwise_pyspi14_CV_linear_SVM_model_permutation_null_inv_prob_pvals.Rds")) %>%
  filter(Noise_Proc == "FC1000") %>%
  ungroup() %>%
  mutate(bal_acc_p_adj_bonf = p.adjust(bal_acc_p, method="bonferroni")) %>%
  dplyr::select(grouping_var, bal_acc_p, bal_acc_p_adj, bal_acc_p_adj_bonf) %>%
  dplyr::rename("SPI" = "grouping_var")

# Load SPI combo data
ABIDE_SPI_combo_main <- readRDS(paste0(ABIDE_rdata_path, "univariate_catch22_pairwise_pyspi14_CV_linear_SVM_inv_prob_balacc.Rds")) %>%
  filter(Noise_Proc == "FC1000") %>%
  group_by(SPI, repeat_number) %>%
  summarise(mean_balacc = mean(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balacc") %>%
  left_join(., ABIDE_SPI_combo_pvals) %>%
  ungroup() %>%
  mutate(SPI = str_replace_all(SPI, "_", " ")) %>%
  mutate(SPI = fct_reorder(SPI, 
                           balanced_accuracy,
                           .fun = mean))

ABIDE_uni_region$Analysis <- "Univariate Region-wise"
ABIDE_uni_region$Variable_Type <- "Univariate"
ABIDE_uni_combo$Analysis <- "Univariate Combo-wise"
ABIDE_uni_combo$Variable_Type <- "Univariate"

ABIDE_pair_feature <- ABIDE_SPI_main %>%
  group_by(grouping_var) %>%
  summarise(mean_balacc = mean(balanced_accuracy)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balacc") %>%
  mutate(Analysis = "Pairwise SPI-wise",
         Variable_Type = "Pairwise")

ABIDE_pair_combo <- ABIDE_SPI_combo_main %>%
  dplyr::rename("grouping_var" = "SPI") %>%
  group_by(grouping_var) %>%
  summarise(mean_balacc = mean(balanced_accuracy)) %>%
  dplyr::rename("balanced_accuracy" = "mean_balacc") %>%
  mutate(Analysis = "Pairwise SPI + Univariate Combo-wise",
         Variable_Type = "Pairwise")

ABIDE_merged_data <- do.call(plyr::rbind.fill, list(ABIDE_uni_region, ABIDE_uni_combo, 
                                                   ABIDE_pair_feature, ABIDE_pair_combo)) %>%
  mutate(Variable_Type = factor(Variable_Type, levels = c("Univariate",
                                                          "Pairwise")),
         Analysis = factor(Analysis, levels = c(
           "Univariate Region-wise", 
           "Univariate Combo-wise", 
           "Pairwise SPI-wise",
           "Pairwise SPI + Univariate Combo-wise")))

ABIDE_merged_data %>%
  ggplot(data=., mapping = aes(x = Analysis, y = balanced_accuracy)) +
  geom_violin(aes(fill = Analysis), color="white") +
  stat_summary(data = subset(ABIDE_merged_data, Analysis == "Univariate Combo-wise"),
               geom = "crossbar", fun = "mean", aes(color=Analysis), size=1,
               color = "#9B51B4") +
  geom_boxplot(fill=NA, width=0.1, color="white") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
  geom_line(data = subset(ABIDE_merged_data, Variable_Type=="Pairwise"),
            aes(x=Analysis, y=balanced_accuracy, group=grouping_var),
            color="white", alpha=0.5) +
  scale_fill_manual(values = c("#F0224B", "chartreuse3", "darkgoldenrod2")) +
  ylab("Balanced Accuracy") +
  theme(legend.position = "none",
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.y.left = element_blank(),
        axis.title.x = element_blank()) +
  theme(axis.line = element_line(color="white"),
        axis.text = element_text(color="white"),
        strip.text = element_text(color="white"),
        axis.title = element_text(color="white"),
        axis.ticks = element_line(color="white"))
ggsave(paste0(plot_dir, "ABIDE_Final_Comparison_Violin.png"),
       width = 7, height = 3, units = "in", dpi = 300)
