################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/pairwise_results/")
TAF::mkdir(plot_path)

# python_to_use <- "~/.conda/envs/pyspi/bin/python3"
python_to_use <- "/Users/abry4213/opt/anaconda3/envs/pyspi/bin/python3"
pairwise_feature_set <- "pyspi14"
data_path <- "~/data/TS_feature_manuscript"
study_group_df <- data.frame(Study = c(rep("UCLA_CNP", 3), "ABIDE_ASD"),
                             Noise_Proc = c(rep("AROMA+2P+GMR",3), "FC1000"),
                             Comparison_Group = c("Schizophrenia", "ADHD", "Bipolar", "ASD"))
univariate_feature_sets <- c("catch22", "catch2", "catch24")

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
theme_set(theme_cowplot())

# Source visualisation script
source(glue("{github_dir}/data_visualisation/Manuscript_Draft_Visualisations_Helper.R"))

# Load in SPI info
SPI_info <- read.csv(glue("{github_dir}/data_visualisation/SPI_info.csv"))

pairwise_null_distribution %>%
  filter(group_var == "cov_EmpiricalCovariance") %>%
  ggplot(data=., mapping=aes(x=Null_Balanced_Accuracy)) +
  geom_histogram()

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
  mutate(Nickname = fct_reorder(Nickname, Balanced_Accuracy_Across_Repeats, .fun=sum),
         Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD"))) %>%
  ggplot(data=., mapping=aes(x=Comparison_Group, y=Nickname, 
                             fill=Balanced_Accuracy_Across_Repeats)) +
  geom_tile()+
  scale_fill_gradientn(colors=c(alpha("#4C7FC0", 0.3), "#4C7FC0"), 
                       na.value=NA, 
                       limits=c(54, 68),
                       breaks=seq(54, 68, by=4)) +
  labs(fill = "Mean Balanced Accuracy (%)") +
  xlab("Clinical Group") +
  ylab("Pairwise SPI") +
  theme(legend.position="bottom")  +
  guides(fill = guide_colorbar(title.position = "top", 
                               nrow = 1,
                               barwidth = 12, 
                               barheight = 1,
                               title.hjust = 0.5)) 
ggsave(glue("{plot_path}/Feature_wise_results.png"),
       width=5.5, height=5.5, units="in", dpi=300)

for (i in 1:nrow(study_group_df)) {
  dataset_ID <- study_group_df$Study[i]
  comparison_group <- study_group_df$Comparison_Group[i]
  
  significant_SPIs <- pairwise_p_values %>%
    filter(Study == dataset_ID,
           Comparison_Group == comparison_group,
           Analysis_Type == "SPI") %>%
    filter(p_value_BH < 0.05) %>%
    pull(group_var)
  
  # Only move forward if 1+ significant brain regions was detected 
  if (length(significant_SPIs) > 0) {
    # Pull out relevant null data
    null_data_to_plot <- pairwise_null_distribution %>%
      dplyr::rename("SPI" = "group_var") %>%
      filter(Study == dataset_ID,
             Comparison_Group == comparison_group,
             Analysis_Type == "SPI") %>%
      left_join(., SPI_info) %>%
      dplyr::rename("group_var" = "Nickname")
    
    # Pull out data for repeats
    repeat_data_to_plot <- pairwise_balanced_accuracy_by_repeats %>%
      dplyr::rename("SPI" = "group_var") %>%
      filter(Study == dataset_ID,
             Comparison_Group == comparison_group,
             SPI %in% significant_SPIs,
             Analysis_Type == "SPI")  %>%
      left_join(., SPI_info) %>%
      dplyr::rename("group_var" = "Nickname")
    
    ### UCLA boxplot with shaded null region
    plot_boxplot_shaded_null(dataset_ID = dataset_ID,
                             grouping_var_name = "",
                             main_data_by_repeat = repeat_data_to_plot,
                             fill_color = "chartreuse3",
                             wrap_length=50,
                             null_mean_value = mean(null_data_to_plot$Null_Balanced_Accuracy, na.rm=T),
                             null_SD_value = sd(null_data_to_plot$Null_Balanced_Accuracy, na.rm=T))
    ggsave(glue("{plot_path}/{dataset_ID}_{comparison_group}_{pairwise_feature_set}_SPI_sig_boxplot.png"),
           width=max(5, 3.5+sqrt(length(significant_SPIs))), 
           height=max(1.25, sqrt(length(significant_SPIs))), units="in", dpi=300)
  }
}
