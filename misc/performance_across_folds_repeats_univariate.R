library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

# Define data path
data_path <- "/headnode1/abry4213/data/UCLA_Schizophrenia/"

# Output figure path
plot_path <- "/headnode1/abry4213/data/UCLA_Schizophrenia/plots/"

# Load in data containing balanced accuracy per fold per repeat
balanced_accuracy_across_folds <- readRDS(paste0(data_path, 
                                                 "processed_data/Rdata/ROI_wise_CV_linear_SVM_catch22_inv_prob_balacc.Rds"))


# Use right postcentral gyrus as an example region

# Find mean and SD for balanced accuracy across folds per each repeat
repeat_avg <- balanced_accuracy_across_folds %>%
  filter(grouping_var == "ctx-rh-postcentral",
         Noise_Proc == "AROMA+2P+GMR") %>%
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
balanced_accuracy_across_folds %>%
  filter(grouping_var == "ctx-rh-postcentral",
         Noise_Proc == "AROMA+2P+GMR") %>%
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
  ylab("Fold Number") +
  xlab("Repeat Number") +
  ggtitle("Balanced Accuracy across 10 Folds/Repeats\nfor Right Postcentral Gyrus") +
  geom_text(data=repeat_avg, aes(label=label)) +
  labs(fill="Balanced Accuracy %") +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust=0.5),
        legend.key.width = unit(1.5,"cm"))
ggsave(paste0(plot_path, "UCLA_Schizophrenia_Balanced_Accuracy_Fold_Repeat_Right_Postcentral.png"),
       width=7, height=5, units="in", dpi=300, bg="white")