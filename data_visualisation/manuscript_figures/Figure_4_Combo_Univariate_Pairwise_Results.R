################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/Figure4/")
TAF::mkdir(plot_path)

python_to_use <- "~/.conda/envs/pyspi/bin/python3"
python_to_use <- "/Users/abry4213/opt/anaconda3/envs/pyspi/bin/python3"
pairwise_feature_set <- "pyspi14"
data_path <- "~/data/TS_feature_manuscript"
UCLA_CNP_sample_metadata <- "~/data/UCLA_CNP/study_metadata/UCLA_CNP_sample_metadata.feather"
ABIDE_ASD_sample_metadata <- "~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_sample_metadata.feather"
study_group_df <- data.frame(Study = c(rep("UCLA_CNP", 3), "ABIDE_ASD"),
                             Num_Samples = c(166, 157, 167, 1150), 
                             Noise_Proc = c(rep("AROMA+2P+GMR",3), "FC1000"),
                             Comparison_Group = c("Schizophrenia", "ADHD", "Bipolar", "ASD"))
univariate_feature_sets <- c("catch22", "catch2", "catch24")

ABIDE_ASD_brain_region_info <- read.csv("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv")

reticulate::use_python(python_to_use)

library(reticulate)

# Import pyarrow.feather as pyarrow_feather
pyarrow_feather <- import("pyarrow.feather")
# DIY rlist::list.append
list.append <- function (.data, ...) 
{
  if (is.list(.data)) {
    c(.data, list(...))
  }
  else {
    c(.data, ..., recursive = FALSE)
  }
}

################################################################################
# Load libraries
################################################################################
library(feather)
library(tidyverse)
library(glue)
library(icesTAF)
library(cowplot)
library(correctR)
library(scales)
theme_set(theme_cowplot())

# Source visualisation script
source(glue("{github_dir}/data_visualisation/manuscript_figures/Manuscript_Draft_Visualisations_Helper.R"))

# Load in metadata
UCLA_CNP_metadata <- pyarrow_feather$read_feather(UCLA_CNP_sample_metadata)
ABIDE_ASD_metadata <- pyarrow_feather$read_feather(ABIDE_ASD_sample_metadata)

# Load in SPI info
SPI_info <- read.csv(glue("{github_dir}/data_visualisation/manuscript_figures/SPI_info.csv"))

# Load data
pairwise_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_robustsigmoid_scaler_balanced_accuracy_all_folds.feather"))
pairwise_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_robustsigmoid_scaler_empirical_p_values.feather"))
pairwise_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_robustsigmoid_scaler_null_balanced_accuracy_distributions.feather"))
combo_univariate_pairwise_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combo_univariate_pairwise_robustsigmoid_scaler_balanced_accuracy_all_folds.feather"))
combo_univariate_pairwise_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combo_univariate_pairwise_robustsigmoid_scaler_empirical_p_values.feather"))
combo_univariate_pairwise_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combo_univariate_pairwise_robustsigmoid_scaler_null_balanced_accuracy_distributions.feather"))

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
  left_join(., combo_univariate_pairwise_p_values %>% dplyr::select(Study:group_var, p_value:p_value_BH))

# Aggregate the main results across all folds, independent of repeat
pairwise_balanced_accuracy <- pairwise_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Pairwise_Feature_Set, Analysis_Type, group_var) %>%
  summarise(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
            Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T))

combo_univariate_pairwise_balanced_accuracy <- combo_univariate_pairwise_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Pairwise_Feature_Set, Analysis_Type, group_var) %>%
  summarise(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
            Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T))

################################################################################
# Figure 4A SPI- and univariate combo-wise SVM results
################################################################################

for (i in 1:nrow(study_group_df)) {
  dataset_ID <- study_group_df$Study[i]
  comparison_group <- study_group_df$Comparison_Group[i]
  
  significant_SPIs <- combo_univariate_pairwise_p_values %>%
    filter(Study == dataset_ID,
           Comparison_Group == comparison_group,
           Univariate_Feature_Set == univariate_feature_set,
           Analysis_Type == "SPI_Combo") %>%
    filter(p_value_BH < 0.05) %>%
    pull(group_var)
  
  # Only move forward if 1+ significant brain regions was detected 
  if (length(significant_SPIs) > 0) {
    # Pull out relevant null data
    null_data_to_plot <- combo_univariate_pairwise_null_distribution %>%
      dplyr::rename("SPI" = "group_var") %>%
      filter(Study == dataset_ID,
             Univariate_Feature_Set == univariate_feature_set,
             Comparison_Group == comparison_group,
             Analysis_Type == "SPI_Combo") %>%
      left_join(., SPI_info) %>%
      dplyr::rename("group_var" = "Nickname")
    
    # Pull out data for repeats
    repeat_data_to_plot <- combo_univariate_pairwise_balanced_accuracy_by_repeats %>%
      dplyr::rename("SPI" = "group_var") %>%
      filter(Study == dataset_ID,
             Comparison_Group == comparison_group,
             Univariate_Feature_Set == univariate_feature_set,
             SPI %in% significant_SPIs,
             Analysis_Type == "SPI_Combo")  %>%
      left_join(., SPI_info) %>%
      dplyr::rename("group_var" = "Nickname")
    
    ### UCLA boxplot with shaded null region
    plot_boxplot_shaded_null(dataset_ID = dataset_ID,
                             grouping_var_name = "",
                             main_data_by_repeat = repeat_data_to_plot,
                             fill_color = "darkgoldenrod2",
                             wrap_length=50,
                             null_mean_value = mean(null_data_to_plot$Null_Balanced_Accuracy, na.rm=T),
                             null_SD_value = sd(null_data_to_plot$Null_Balanced_Accuracy, na.rm=T))
    ggsave(glue("{plot_path}/{dataset_ID}_{comparison_group}_combo_{pairwise_feature_set}_{univariate_feature_set}_sig_boxplot.png"),
           width=max(5, 3.5+sqrt(length(significant_SPIs))), 
           height=max(1.25, sqrt(length(significant_SPIs))), units="in", dpi=300)
    
  }
}

################################################################################
# Figure 4B Comparison of each SPI with vs. without univariate info
################################################################################

# Merge the pairwise and combo balanced accuracy datasets
merged_balanced_accuracy_all_folds <- plyr::rbind.fill(pairwise_balanced_accuracy_all_folds,
                                                       combo_univariate_pairwise_balanced_accuracy_all_folds) %>%
  mutate(Analysis_Type = factor(Analysis_Type, levels = c("SPI",
                                                          "SPI_Combo")))


# Run repeated k-fold cross-validation correction
rkcv_list <- list()
for (i in 1:nrow(study_group_df)) {
  dataset_ID <- study_group_df$Study[i]
  comparison_group <- study_group_df$Comparison_Group[i]
  num_samples <- study_group_df$Num_Samples[i]
  
  
  data_for_correction <- merged_balanced_accuracy_all_folds %>%
    filter(Comparison_Group == comparison_group) %>%
    dplyr::rename("model" = "Analysis_Type",
                  "values" = "Balanced_Accuracy",
                  "k" = "Fold",
                  "r" = "Repeat_Number") %>%
    mutate(r = r + 1)
  
  SPI_only_data <- data_for_correction %>% 
    filter(model == "SPI")
  SPI_combo_data <- data_for_correction %>% 
    filter(model == "SPI_Combo", 
           Univariate_Feature_Set == univariate_feature_set)
  merged_data_for_correction <- plyr::rbind.fill(SPI_only_data,
                                                 SPI_combo_data)
  
  # Iterate over each SPI
  for (this_SPI in unique(merged_data_for_correction$group_var)) {
    tryCatch({SPI_rkcv <- repkfold_ttest(data = subset(merged_data_for_correction, group_var == this_SPI),
                                         n1 = ceiling(0.9*num_samples), 
                                         n2 = floor(0.1*num_samples), 
                                         k = 10, 
                                         r = 10)
    
    SPI_rkcv$SPI <- this_SPI
    SPI_rkcv$Study <- dataset_ID
    SPI_rkcv$Comparison_Group <- comparison_group
    SPI_rkcv$Univariate_Feature_Set <- univariate_feature_set
    
    # Append to list
    rkcv_list <- list.append(rkcv_list, SPI_rkcv)
    }, error = function(e) {
      print(e)
      cat("Error for SPI:", this_SPI, "\n")
    })
    
  }
}

rkcv_results <- do.call(plyr::rbind.fill, rkcv_list)

# Plot each SPI with vs without univariate info as a violin
for (i in 1:nrow(study_group_df)) {
  dataset_ID <- study_group_df$Study[i]
  comparison_group <- study_group_df$Comparison_Group[i]
  
  SPI_data <- pairwise_balanced_accuracy_by_repeats %>%
    filter(Comparison_Group == comparison_group)
  SPI_combo_data <- combo_univariate_pairwise_balanced_accuracy_by_repeats %>%
    filter(Comparison_Group == comparison_group,
           Univariate_Feature_Set == univariate_feature_set)
  
  merged_data_to_plot <- plyr::rbind.fill(SPI_data, SPI_combo_data) %>%
    mutate(fill_color = case_when(Analysis_Type == "SPI" & p_value_BH < 0.05 ~ "chartreuse3",
                                  Analysis_Type == "SPI_Combo" & p_value_BH < 0.05 ~ "darkgoldenrod2",
                                  Analysis_Type == "SPI" & p_value_BH > 0.05 ~ "gray90_1",
                                  T ~ "gray90_2")) %>%
    dplyr::rename("SPI" = "group_var")  %>%
    left_join(., SPI_info) %>%
    mutate(Nickname = fct_reorder(Nickname, Balanced_Accuracy_Across_Folds, 
                                  .fun=mean, .desc=T)) 
  
  merged_data_to_plot %>%
    ggplot(data=., mapping=aes(x = Nickname, y = Balanced_Accuracy_Across_Folds,
                               fill = fill_color)) +
    geom_violin() +
    facet_wrap(Nickname ~ ., scales="free_x", nrow=1) +
    xlab("Pairwise Statistic (SPI)") + 
    ylab("Balanced Accuracy\nby Repeat CV") +
    scale_fill_manual(values = c("chartreuse3", "darkgoldenrod2", "gray90", "gray90")) +
    scale_x_discrete(labels = wrap_format(20)) +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.4),
          strip.text = element_blank()) +
    theme(legend.position = "none")
  ggsave(glue("{plot_path}/{dataset_ID}_{comparison_group}_combo_{pairwise_feature_set}_{univariate_feature_set}_violins.png"),
         width = 8.5, height = 4, units = "in", dpi = 300)
}
