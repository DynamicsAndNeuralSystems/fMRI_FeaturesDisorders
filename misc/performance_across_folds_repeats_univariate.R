library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

# # UCLA Schizophrenia
# dataset_ID <- "UCLA_Schizophrenia"
# data_path <- "/headnode1/abry4213/data/UCLA_Schizophrenia/"
# rdata_path <- paste0(data_path, "processed_data/Rdata/")
# plot_path <- "/headnode1/abry4213/data/UCLA_Schizophrenia/plots/"
# brain_region_interest <- "ctx-rh-postcentral"
# feature_interest <- "SP_Summaries_welch_rect_centroid"
# noise_proc <- "AROMA+2P+GMR"

# ABIDE ASD
dataset_ID <- "ABIDE_ASD"
data_path <- "/headnode1/abry4213/data/ABIDE_ASD/"
rdata_path <- paste0(data_path, "processed_data/Rdata/")
plot_path <- "/headnode1/abry4213/data/ABIDE_ASD/plots/"
brain_region_interest <- "Superior Frontal Gyrus"
noise_proc <- "FC1000"
feature_interest <- "FC_LocalSimple_mean3_stderr"

############### ROI-wise ############### 
# Load in data containing balanced accuracy per fold per repeat
ROI_wise_balanced_accuracy_across_folds <- readRDS(paste0(rdata_path, "ROI_wise_CV_linear_SVM_catch22_inv_prob_balacc.Rds"))

# Find mean and SD for balanced accuracy across folds per each repeat
repeat_avg <- ROI_wise_balanced_accuracy_across_folds %>%
  filter(grouping_var == brain_region_interest,
         Noise_Proc == noise_proc) %>%
  group_by(repeat_number) %>%
  summarise(mean = mean(balanced_accuracy, na.rm=T),
            SD = sd(balanced_accuracy, na.rm=T)) %>%
  pivot_longer(cols=c(mean, SD),
               names_to="Metric",
               values_to="balanced_accuracy") %>%
  mutate(fold_number = ifelse(Metric=="mean", "Mean", "SD"),
         label = round(100*balanced_accuracy, 2)) %>%
  mutate(fold_number = factor(fold_number, levels=rev(c(1:10, "Mean", "SD"))))

# Plot balanced accuracy by fold X repeat
ROI_wise_balanced_accuracy_across_folds %>%
  filter(grouping_var == brain_region_interest,
         Noise_Proc == noise_proc) %>%
  mutate(repeat_number = as.factor(repeat_number),
         fold_number = factor(fold_number, levels=rev(c(1:10, "Mean", "SD")))) %>%
  ggplot(data=., mapping=aes(x=repeat_number, 
                             y=fold_number)) +
  # Fill in tiles per fold/repeat
  geom_tile(aes(fill=100*balanced_accuracy)) +
  # Mean tiles
  geom_tile(data=repeat_avg %>% filter(Metric=="mean"), 
            mapping=aes(fill=100*balanced_accuracy)) +
  # Standard deviation tiles (gray)
  geom_tile(data=repeat_avg %>% filter(Metric=="SD"), fill="gray70") +
  scale_fill_gradient(low="white", high="#F0224B") +
  ylab("Fold Number") +
  xlab("Repeat Number") +
  ggtitle(paste0("Balanced Accuracy across 10 Folds/Repeats\nfor ", dataset_ID, " ", brain_region_interest)) +
  geom_text(data=repeat_avg, aes(label=label)) +
  labs(fill="Balanced Accuracy %") +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust=0.5),
        legend.key.width = unit(1.5,"cm"))
ggsave(paste0(plot_path, dataset_ID, "_Balanced_Accuracy_Fold_Repeat_",
              brain_region_interest, ".png"),
       width=7, height=5, units="in", dpi=300, bg="white")


############### Feature-wise ############### 
# Load in data containing balanced accuracy per fold per repeat
Feature_wise_balanced_accuracy_across_folds <- readRDS(paste0(rdata_path, 
                                                          "Feature_wise_CV_linear_SVM_catch22_inv_prob_balacc.Rds"))

# Find mean and SD for balanced accuracy across folds per each repeat
repeat_avg <- Feature_wise_balanced_accuracy_across_folds %>%
  filter(grouping_var == feature_interest,
         Noise_Proc == noise_proc) %>%
  group_by(repeat_number) %>%
  summarise(mean = mean(balanced_accuracy, na.rm=T),
            SD = sd(balanced_accuracy, na.rm=T)) %>%
  pivot_longer(cols=c(mean, SD),
               names_to="Metric",
               values_to="balanced_accuracy") %>%
  mutate(fold_number = ifelse(Metric=="mean", "Mean", "SD"),
         label = round(100*balanced_accuracy, 2)) %>%
  mutate(fold_number = factor(fold_number, levels=rev(c(1:10, "Mean", "SD"))))

# Plot balanced accuracy by fold X repeat
Feature_wise_balanced_accuracy_across_folds %>%
  filter(grouping_var == feature_interest,
         Noise_Proc == noise_proc) %>%
  mutate(repeat_number = as.factor(repeat_number),
         fold_number = factor(fold_number, levels=rev(c(1:10, "Mean", "SD")))) %>%
  ggplot(data=., mapping=aes(x=repeat_number, 
                             y=fold_number)) +
  # Fill in tiles per fold/repeat
  geom_tile(aes(fill=100*balanced_accuracy)) +
  # Mean tiles
  geom_tile(data=repeat_avg %>% filter(Metric=="mean"), 
            mapping=aes(fill=100*balanced_accuracy)) +
  # Standard deviation tiles (gray)
  geom_tile(data=repeat_avg %>% filter(Metric=="SD"), fill="gray70") +
  scale_fill_gradient(low="white", high="#5C86C6") +
  ylab("Fold Number") +
  xlab("Repeat Number") +
  ggtitle(paste0("Balanced Accuracy across 10 Folds/Repeats\nfor ", 
                 dataset_ID, "\n", feature_interest)) +
  geom_text(data=repeat_avg, aes(label=label)) +
  labs(fill="Balanced Accuracy %") +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust=0.5),
        legend.key.width = unit(1.5,"cm"))
ggsave(paste0(plot_path, dataset_ID, 
              "_Balanced_Accuracy_Fold_Repeat_",
              feature_interest, ".png"),
       width=7, height=6, units="in", dpi=300, bg="white")

############### Combo-wise ############### 
# Load in data containing balanced accuracy per fold per repeat
Combo_wise_balanced_accuracy_across_folds <- readRDS(paste0(rdata_path, 
                                                              "Combo_wise_CV_linear_SVM_catch22_inv_prob_balacc.Rds"))

# Find mean and SD for balanced accuracy across folds per each repeat
repeat_avg <- Combo_wise_balanced_accuracy_across_folds %>%
  filter(Noise_Proc == noise_proc) %>%
  group_by(repeat_number) %>%
  summarise(mean = mean(balanced_accuracy, na.rm=T),
            SD = sd(balanced_accuracy, na.rm=T)) %>%
  pivot_longer(cols=c(mean, SD),
               names_to="Metric",
               values_to="balanced_accuracy") %>%
  mutate(fold_number = ifelse(Metric=="mean", "Mean", "SD"),
         label = round(100*balanced_accuracy, 2)) %>%
  mutate(fold_number = factor(fold_number, levels=rev(c(1:10, "Mean", "SD"))))

# Plot balanced accuracy by fold X repeat
Combo_wise_balanced_accuracy_across_folds %>%
  filter(Noise_Proc == noise_proc) %>%
  mutate(repeat_number = as.factor(repeat_number),
         fold_number = factor(fold_number, levels=rev(c(1:10, "Mean", "SD")))) %>%
  ggplot(data=., mapping=aes(x=repeat_number, 
                             y=fold_number)) +
  # Fill in tiles per fold/repeat
  geom_tile(aes(fill=100*balanced_accuracy)) +
  # Mean tiles
  geom_tile(data=repeat_avg %>% filter(Metric=="mean"), 
            mapping=aes(fill=100*balanced_accuracy)) +
  # Standard deviation tiles (gray)
  geom_tile(data=repeat_avg %>% filter(Metric=="SD"), fill="gray70") +
  scale_fill_gradient(low="white", high="#9B51B4") +
  ylab("Fold Number") +
  xlab("Repeat Number") +
  ggtitle(paste0("Balanced Accuracy across 10 Folds/Repeats\nfor ",
          dataset_ID, "\nUnivariate Combo-wise Analysis")) +
  geom_text(data=repeat_avg, aes(label=label)) +
  labs(fill="Balanced Accuracy %") +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust=0.5),
        legend.key.width = unit(1.5,"cm"))
ggsave(paste0(plot_path, dataset_ID, "_Balanced_Accuracy_Fold_Repeat_Univariate_Combo_Wise.png"),
       width=7, height=6, units="in", dpi=300, bg="white")