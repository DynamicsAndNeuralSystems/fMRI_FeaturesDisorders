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
library(ggseg)
library(ggsegHO)
theme_set(theme_cowplot())

# Source visualisation script
source(glue("{github_dir}/data_visualisation/manuscript_figures/Manuscript_Draft_Visualisations_Helper.R"))

# Load in SPI info
SPI_info <- read.csv(glue("{github_dir}/data_visualisation/manuscript_figures/SPI_info.csv"))

# Load data
combo_univariate_pairwise_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combo_univariate_pairwise_balanced_accuracy_all_folds.feather"))
combo_univariate_pairwise_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combo_univariate_pairwise_empirical_p_values.feather"))
combo_univariate_pairwise_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combo_univariate_pairwise_null_balanced_accuracy_distributions.feather"))

# Aggregate balanced accuracy by repeats
combo_univariate_pairwise_balanced_accuracy_by_repeats <- combo_univariate_pairwise_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Pairwise_Feature_Set, Analysis_Type, group_var, Repeat_Number) %>%
  summarise(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
            Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T)) %>%
  left_join(., combo_univariate_pairwise_p_values %>% dplyr::select(Study:group_var, p_value:p_value_BH))

################################################################################
# Figure 4 SPI- and univariate combo-wise SVM results
################################################################################

for (i in 1:nrow(study_group_df)) {
  dataset_ID <- study_group_df$Study[i]
  comparison_group <- study_group_df$Comparison_Group[i]
  
  for (featset in c("catch2", "catch22", "catch24")) {
    significant_SPIs <- combo_univariate_pairwise_p_values %>%
      filter(Study == dataset_ID,
             Comparison_Group == comparison_group,
             Univariate_Feature_Set == featset,
             Analysis_Type == "SPI_Combo") %>%
      filter(p_value_BH < 0.05) %>%
      pull(group_var)
    
    # Only move forward if 1+ significant brain regions was detected 
    if (length(significant_SPIs) > 0) {
      # Pull out relevant null data
      null_data_to_plot <- combo_univariate_pairwise_null_distribution %>%
        dplyr::rename("SPI" = "group_var") %>%
        filter(Study == dataset_ID,
               Univariate_Feature_Set == featset,
               Comparison_Group == comparison_group,
               Analysis_Type == "SPI_Combo") %>%
        left_join(., SPI_info) %>%
        dplyr::rename("group_var" = "Nickname")
      
      # Pull out data for repeats
      repeat_data_to_plot <- combo_univariate_pairwise_balanced_accuracy_by_repeats %>%
        dplyr::rename("SPI" = "group_var") %>%
        filter(Study == dataset_ID,
               Comparison_Group == comparison_group,
               Univariate_Feature_Set == featset,
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
      ggsave(glue("{plot_path}/{dataset_ID}_{comparison_group}_combo_{pairwise_feature_set}_{featset}_sig_boxplot.png"),
             width=max(5, 3.5+sqrt(length(significant_SPIs))), 
             height=max(1.25, sqrt(length(significant_SPIs))), units="in", dpi=300)
    }
  }
}
