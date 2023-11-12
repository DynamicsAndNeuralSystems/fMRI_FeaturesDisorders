################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/Library/CloudStorage/OneDrive-TheUniversityofSydney(Students)/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/univariate_results/")
TAF::mkdir(plot_path)

# python_to_use <- "~/.conda/envs/pyspi/bin/python3"
python_to_use <- "/Users/abry4213/anaconda3/envs/pyspi/bin/python3"
pairwise_feature_set <- "pyspi14"
univariate_feature_set <- "catch25"
SVM_kernel <- "Linear"
data_path <- "~/data/TS_feature_manuscript"
study_group_df <- data.frame(Study = c(rep("UCLA_CNP", 3), "ABIDE_ASD"),
                             Noise_Proc = c(rep("AROMA+2P+GMR",3), "FC1000"),
                             Comparison_Group = c("Schizophrenia", "Bipolar", "ADHD", "ASD"),
                             Group_Nickname = c("SCZ", "BP", "ADHD", "ASD"))
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
library(ggsegDefaultExtra)
library(patchwork)
library(broom)
library(colorspace)
library(scales)
library(correctR)
library(splitstackshape)
theme_set(theme_cowplot())

# Source visualisation script
source(glue("{github_dir}/data_visualisation/Manuscript_Draft_Visualisations_Helper.R"))

UCLA_CNP_brain_region_info <- read.csv("~/data/UCLA_CNP/study_metadata/UCLA_CNP_Brain_Region_info.csv")
ABIDE_ASD_brain_region_info <- read.table("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_Harvard_Oxford_cort_prob_2mm_ROI_lookup.txt", sep=";", header = T) %>%
  mutate(Brain_Region = ifelse(Index==45, "Heschl's Gyrus (includes H1 and H2)", Brain_Region))

# Load in univariate time-series feature info
TS_feature_info <- read.csv(glue("{github_dir}/data_visualisation/catch25_info.csv"))

# Load univariate classification results across all folds
univariate_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_balanced_accuracy_all_folds.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set, kernel==SVM_kernel)
# Compute mean + SD performance across all folds
univariate_balanced_accuracy <- univariate_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Analysis_Type, group_var, kernel) %>%
  reframe(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
          Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T))

# Load p-values
univariate_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_empirical_p_values.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set) %>%
  dplyr::select(-Balanced_Accuracy_Across_Folds) %>%
  left_join(., univariate_balanced_accuracy)

# Load null distribution
univariate_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_null_balanced_accuracy_distributions.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)

# Load study metadata
UCLA_CNP_metadata <- pyarrow_feather$read_feather("~/data/UCLA_CNP/study_metadata/UCLA_CNP_sample_metadata.feather") %>%
  mutate(Study="UCLA_CNP")
ABIDE_ASD_metadata <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_sample_metadata.feather")  %>%
  mutate(Study="ABIDE_ASD")

################################################################################
# Univariate region-wise analysis
################################################################################
         
# Plot balanced accuracy in the brain
plot_balacc_in_brain <- function(significant_univariate_region_wise_results, 
                                 color_palette=c("#FFEE75", "#FCA769", "#fb6555", "#D32345", "#401057"),
                                 bin_seq_range=seq(50,75,by=5)) {

  # Find max fill and min fill values
  min_fill <- floor(min(significant_univariate_region_wise_results$Balanced_Accuracy_Across_Folds))
  max_fill <- ceiling(max(significant_univariate_region_wise_results$Balanced_Accuracy_Across_Folds))
  
  # Initialize list of ggseg plots
  ggseg_plot_list <- list()
  
  # First plot within brain using ggseg
  for (i in 1:nrow(study_group_df)) {
    dataset_ID <- study_group_df$Study[i]
    comparison_group <- study_group_df$Comparison_Group[i]
    
    # Define atlas by study
    atlas <- ifelse(dataset_ID == "UCLA_CNP", "dk", "hoCort")
    
    # If dataset is ABIDE ASD, convert regions to ggseg regions
    if (dataset_ID == "ABIDE_ASD") {
      significant_data_for_ggseg <- significant_univariate_region_wise_results %>%
        filter(Study == dataset_ID,
               Comparison_Group == comparison_group) %>%
        dplyr::rename("Brain_Region" = "group_var") %>%
        left_join(., ABIDE_ASD_brain_region_info)
      
    } else {
      # Extract sig results to plot
      significant_data_for_ggseg <- significant_univariate_region_wise_results %>%
        filter(Study == dataset_ID,
               Comparison_Group == comparison_group) %>%
        distinct() %>%
        mutate(label = ifelse(str_detect(group_var, "ctx-"),
                              gsub("-", "_", group_var),
                              as.character(group_var))) %>%
        mutate(label = gsub("ctx_", "", label))
    }
    
    
    # Plot balanced accuracy data in cortex
    dataset_ggseg <- plot_data_with_ggseg_discrete(dataset_ID = dataset_ID,
                                                   atlas_name = atlas,
                                                   atlas_data = get(atlas) %>% as_tibble(),
                                                   data_to_plot = significant_data_for_ggseg,
                                                   fill_variable = "Balanced_Accuracy_Across_Folds",
                                                   fill_colors = color_palette,
                                                   bin_seq = bin_seq_range,
                                                   line_color = "gray30",
                                                   na_color = "white") +
      labs(fill = "Mean Balanced Accuracy (%)") +
      theme(plot.title = element_blank())
    
    # Append to list
    ggseg_plot_list <- list.append(ggseg_plot_list, dataset_ggseg)
    
    # Add subcortical data for UCLA CNP
    if (dataset_ID == "UCLA_CNP") {
      dataset_ggseg_subctx <- plot_data_with_ggseg_discrete(dataset_ID = dataset_ID,
                                                            atlas_name = "aseg",
                                                            atlas_data = aseg %>% as_tibble(),
                                                            data_to_plot = significant_data_for_ggseg,
                                                            fill_variable = "Balanced_Accuracy_Across_Folds",
                                                            fill_colors = color_palette,
                                                            bin_seq = bin_seq_range,
                                                            line_color = "gray30",
                                                            na_color = "white") +
        labs(fill = "Mean Balanced Accuracy (%)") +
        theme(plot.title = element_blank()) 
      # Append to list
      ggseg_plot_list <- list.append(ggseg_plot_list, dataset_ggseg_subctx)
    }
  }
  return(ggseg_plot_list)
}

# catch25 features
# Define dataset with univariate region-wise results
significant_univariate_region_wise_results <- univariate_p_values %>%
  filter(Analysis_Type == "Univariate_Brain_Region") %>%
  filter(p_value_HolmBonferroni < 0.05) %>%
  mutate(Balanced_Accuracy_Across_Folds = 100*Balanced_Accuracy_Across_Folds)

catch25_region_wise_balacc_plot_list <- plot_balacc_in_brain(significant_univariate_region_wise_results,
                                                             bin_seq_range=seq(50,75,by=5))
wrap_plots(catch25_region_wise_balacc_plot_list, 
           ncol=2, 
           byrow=T) + 
  plot_layout(guides = "collect")
ggsave(glue("{plot_path}/Region_wise_balacc_ggseg.svg"),
       width=5, height=7, units="in", dpi=300)


# Get region count by group
significant_univariate_region_wise_results %>%
  count(Comparison_Group)

# Save regional results to a CSV file as a table
univariate_p_values %>%
  filter(Analysis_Type=="Univariate_Brain_Region") %>%
  dplyr::select(Comparison_Group, group_var, p_value_HolmBonferroni, Balanced_Accuracy_Across_Folds, Balanced_Accuracy_Across_Folds_SD) %>%
  dplyr::rename("Brain_Region" = "group_var",
                "Disorder" = "Comparison_Group") %>%
  mutate(Balanced_Accuracy_Across_Folds = round(100*Balanced_Accuracy_Across_Folds,1),
         Balanced_Accuracy_Across_Folds_SD = round(100*Balanced_Accuracy_Across_Folds_SD,1),
         Disorder = case_when(Disorder == "Bipolar" ~ "BP",
                              Disorder == "Schizophrenia" ~ "SCZ",
                              T ~ Disorder)) %>%
  write.csv(., glue("{plot_path}/Region_wise_performance_results.csv"),
            row.names = F)