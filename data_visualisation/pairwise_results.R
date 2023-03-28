################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/pairwise_results/")
TAF::mkdir(plot_path)

# python_to_use <- "~/.conda/envs/pyspi/bin/python3"
python_to_use <- "/Users/abry4213/opt/anaconda3/envs/pyspi/bin/python3"
pairwise_feature_set <- "pyspi14"
univariate_feature_set <- "catch24"
data_path <- "~/data/TS_feature_manuscript"
study_group_df <- data.frame(Study = c(rep("UCLA_CNP", 3), "ABIDE_ASD"),
                             Noise_Proc = c(rep("AROMA+2P+GMR",3), "FC1000"),
                             Comparison_Group = c("Schizophrenia", "ADHD", "Bipolar", "ASD"))

ABIDE_ASD_brain_region_info <- read.csv("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv")

reticulate::use_python(python_to_use)

library(reticulate)

# Import pyarrow.feather as pyarrow_feather
pyarrow_feather <- import("pyarrow.feather")

################################################################################
# Load libraries
################################################################################
library(feather)
library(tidyverse)
library(glue)
library(icesTAF)
library(cowplot)
library(ggnewscale)
library(scales)
theme_set(theme_cowplot())

# Source visualisation script
source(glue("{github_dir}/data_visualisation/Manuscript_Draft_Visualisations_Helper.R"))

# Load in SPI info
SPI_info <- read.csv(glue("{github_dir}/data_visualisation/SPI_info.csv"))

# Load stats data
pairwise_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_mixedsigmoid_scaler_balanced_accuracy_all_folds.feather"))
pairwise_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_mixedsigmoid_scaler_empirical_p_values.feather"))
pairwise_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_mixedsigmoid_scaler_null_balanced_accuracy_distributions.feather"))

combo_univariate_pairwise_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combined_univariate_pairwise_mixedsigmoid_scaler_balanced_accuracy_all_folds.feather"))
# combo_univariate_pairwise_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_mixedsigmoid_scaler_empirical_p_values.feather"))
# combo_univariate_pairwise_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_mixedsigmoid_scaler_null_balanced_accuracy_distributions.feather"))

pairwise_all_SPIs_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_all_SPIs_mixedsigmoid_scaler_balanced_accuracy_all_folds.feather"))
# pairwise_all_SPIs_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_all_SPIs_mixedsigmoid_scaler_empirical_p_values.feather"))
# pairwise_all_SPIs_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_all_SPIs_pairwise_mixedsigmoid_scaler_null_balanced_accuracy_distributions.feather"))

# Aggregate the main results across folds and then across repeats
pairwise_balanced_accuracy <- pairwise_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Pairwise_Feature_Set, Analysis_Type, group_var, Repeat_Number) %>%
  reframe(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
          Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T)) %>%
  group_by(Study, Comparison_Group, Pairwise_Feature_Set, Analysis_Type, group_var) %>%
  reframe(Balanced_Accuracy_Across_Repeats = mean(Balanced_Accuracy_Across_Folds, na.rm=T),
          Balanced_Accuracy_Across_Repeats_SD = sd(Balanced_Accuracy_Across_Folds, na.rm=T))

combo_univariate_pairwise_balanced_accuracy <- combo_univariate_pairwise_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Pairwise_Feature_Set, Analysis_Type, group_var, Repeat_Number) %>%
  reframe(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
          Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T)) %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Pairwise_Feature_Set, Analysis_Type, group_var) %>%
  reframe(Balanced_Accuracy_Across_Repeats = mean(Balanced_Accuracy_Across_Folds, na.rm=T),
          Balanced_Accuracy_Across_Repeats_SD = sd(Balanced_Accuracy_Across_Folds, na.rm=T))

# Aggregate balanced accuracy by repeats
pairwise_balanced_accuracy_by_repeats <- pairwise_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Pairwise_Feature_Set, Analysis_Type, group_var, Repeat_Number) %>%
  summarise(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
            Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T)) %>%
  left_join(., pairwise_p_values %>% dplyr::select(Study:group_var, p_value:p_value_BH))

combo_univariate_pairwise_balanced_accuracy_by_repeats <- combo_univariate_pairwise_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Pairwise_Feature_Set, Analysis_Type, group_var, Repeat_Number) %>%
  summarise(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
            Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T)) %>%
  left_join(., pairwise_p_values %>% dplyr::select(Study:group_var, p_value:p_value_BH))

################################################################################
# SPI-wise SVM results
################################################################################

# Annotation bar with SPI type
pairwise_p_values %>%
  filter(Pairwise_Feature_Set == pairwise_feature_set,
         Analysis_Type == "Pairwise_SPI",
         p_value_Bonferroni < 0.05) %>%
  dplyr::rename("SPI" = "group_var") %>%
  left_join(., SPI_info) %>%
  mutate(Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                      Comparison_Group == "Bipolar" ~ "BPD",
                                      T ~ Comparison_Group),
         Balanced_Accuracy_Across_Repeats = 100*Balanced_Accuracy_Across_Repeats) %>%
  mutate(Nickname = fct_reorder(Nickname, Balanced_Accuracy_Across_Repeats, .fun=mean),
         Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD"))) %>%
  group_by(SPI, Category) %>%
  summarise(Balacc_Sum = sum(Balanced_Accuracy_Across_Repeats)) %>%
  ungroup() %>%
  mutate(SPI = fct_reorder(SPI, Balacc_Sum),
         Category = fct_reorder(Category, Balacc_Sum, .fun=sum, .desc=T)) %>%
  ggplot(data=., mapping=aes(x=0, y=SPI, fill=Category)) +
  geom_tile() +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text=element_text(size=14)) +
  guides(fill = guide_legend(title.position = "top", 
                             ncol = 2,
                             byrow=T,
                             title.hjust = 0.5)) 
ggsave(glue("{plot_path}/SPI_wise_colorbar.png"),
       width=6, height=6, units="in", dpi=300)

# Actual heatmap
pairwise_p_values %>%
  filter(Pairwise_Feature_Set == pairwise_feature_set,
         Analysis_Type == "Pairwise_SPI",
         p_value_Bonferroni < 0.05) %>%
  dplyr::rename("SPI" = "group_var") %>%
  left_join(., SPI_info) %>%
  mutate(Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                      Comparison_Group == "Bipolar" ~ "BPD",
                                      T ~ Comparison_Group),
         Balanced_Accuracy_Across_Repeats = 100*Balanced_Accuracy_Across_Repeats) %>%
  mutate(Nickname = fct_reorder(Nickname, Balanced_Accuracy_Across_Repeats, .fun=mean),
         Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD"))) %>%
  ggplot(data=., mapping=aes(x=Comparison_Group, y=Nickname, 
                             fill=Balanced_Accuracy_Across_Repeats)) +
  geom_tile()+
  geom_text(aes(label = round(Balanced_Accuracy_Across_Repeats, 1))) +
  scale_fill_gradientn(colors=c(alpha("#AC77BD", 0.3), "#AC77BD"), 
                       na.value=NA)  + 
  scale_y_discrete(labels = wrap_format(28)) +
  labs(fill = "Mean Balanced Accuracy (%)") +
  xlab("Clinical Group") +
  ylab("Pairwise SPI") +
  theme(legend.position="bottom")  +
  guides(fill = guide_colorbar(title.position = "top", 
                               nrow = 1,
                               barwidth = 12, 
                               barheight = 1,
                               title.hjust = 0.5)) 
ggsave(glue("{plot_path}/SPI_wise_results.png"),
       width=5, height=4.5, units="in", dpi=300)

################################################################################
# Region pair-wise T-tests for Pearson
demo_data <- 

################################################################################
# Plot significant results relative to their respective null distributions
null_data_for_plot <- pairwise_null_distribution %>%
  filter(Analysis_Type == "Pairwise_SPI") %>%
  dplyr::select(Study, Comparison_Group, Null_Balanced_Accuracy) %>%
  mutate(Type = "Null") %>%
  dplyr::rename("Balanced_Accuracy_Across_Repeats" = "Null_Balanced_Accuracy") %>%
  plyr::rbind.fill(., pairwise_p_values %>%
                     filter(Pairwise_Feature_Set == pairwise_feature_set,
                            Analysis_Type == "Pairwise_SPI",
                            p_value_Bonferroni < 0.05) %>%
                     dplyr::select(Study, Comparison_Group, Balanced_Accuracy_Across_Repeats) %>%
                     mutate(Type = "Main")) %>%
  mutate(Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                      Comparison_Group == "Bipolar" ~ "BPD",
                                      T ~ Comparison_Group)) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD")))
null_data_for_plot %>%
  ggplot(data=., mapping=aes(x=100*Balanced_Accuracy_Across_Repeats)) +
  geom_histogram(data = subset(null_data_for_plot, Type=="Null"),
                 fill="gray80", bins=50) +
  geom_vline(data = subset(null_data_for_plot, Type=="Main"),
             aes(xintercept = 100*Balanced_Accuracy_Across_Repeats,
                 color = Comparison_Group),
             linewidth=0.5) +
  scale_color_manual(values=c("Control" = "#5BB67B", 
                              "SCZ" = "#573DC7", 
                              "BPD" = "#D5492A", 
                              "ADHD" = "#0F9EA9", 
                              "ASD" = "#C47B2F")) +
  facet_wrap(Comparison_Group ~ ., ncol=2, scales="free") +
  xlab("Balanced Accuracy Across Repeats (%)") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.placement = "outside") +
  theme(legend.position="none")
ggsave(glue("{plot_path}/pairwise_SPI_main_vs_null_balanced_acc.png"),
       width=6, height=3.5, units="in", dpi=300)


################################################################################
# Compare SPIs with vs without univariate info
################################################################################

plyr::rbind.fill(pairwise_p_values,
                 combo_univariate_pairwise_balanced_accuracy) %>%
  mutate(Analysis_Type = ifelse(Analysis_Type == "Pairwise_SPI", "SPI Only", "SPI + Univariate"),
         sig = ifelse(p_value_Bonferroni < 0.05, "Significant", "Not significant")) %>%
  mutate(Analysis_Type = factor(Analysis_Type, levels=c("SPI Only", "SPI + Univariate"))) %>%
  mutate(Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                      Comparison_Group == "Bipolar" ~ "BPD",
                                      T ~ Comparison_Group)) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD"))) %>%
  ggplot(data=., mapping=aes(x=Analysis_Type, y=Balanced_Accuracy_Across_Repeats,
                             group = group_var)) +
  geom_line(aes(color = Comparison_Group), show.legend = FALSE) +
  scale_color_manual(values=c("Control" = "#5BB67B", 
                              "SCZ" = "#573DC7", 
                              "BPD" = "#D5492A", 
                              "ADHD" = "#0F9EA9", 
                              "ASD" = "#C47B2F")) +
  new_scale_colour() +  # start a new scale
  geom_point(aes(color = sig)) +
  scale_color_manual(values = c("gray60", "#5BB67B")) +
  scale_x_discrete(labels = wrap_format(7)) +
  xlab("Analysis Type") +
  ylab("Mean Balanced Accuracy (%)") +
  facet_wrap(Comparison_Group ~ ., ncol=2, scales="free_y") +
  theme(legend.position = "bottom",
        legend.title = element_blank())
ggsave(glue("{plot_path}/SPI_with_vs_without_univariate_spaghetti.png"),
       width=5, height=5.5, units="in", dpi=300)