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
library(scales)
theme_set(theme_cowplot())

# Source visualisation script
source(glue("{github_dir}/data_visualisation/Manuscript_Draft_Visualisations_Helper.R"))

# Load in SPI info
SPI_info <- read.csv(glue("{github_dir}/data_visualisation/SPI_info.csv"))

# Load data
pairwise_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_mixedsigmoid_scaler_balanced_accuracy_all_folds.feather"))
pairwise_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_mixedsigmoid_scaler_empirical_p_values.feather"))
pairwise_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_mixedsigmoid_scaler_null_balanced_accuracy_distributions.feather"))

# Aggregate the main results across folds and then across repeats
pairwise_balanced_accuracy <- pairwise_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Pairwise_Feature_Set, Analysis_Type, group_var, Repeat_Number) %>%
  reframe(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
          Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T)) %>%
  group_by(Study, Comparison_Group, Pairwise_Feature_Set, Analysis_Type, group_var) %>%
  reframe(Balanced_Accuracy_Across_Repeats = mean(Balanced_Accuracy_Across_Folds, na.rm=T),
          Balanced_Accuracy_Across_Repeats_SD = sd(Balanced_Accuracy_Across_Folds, na.rm=T))

# Aggregate balanced accuracy by repeats
pairwise_balanced_accuracy_by_repeats <- pairwise_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Pairwise_Feature_Set, Analysis_Type, group_var, Repeat_Number) %>%
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
