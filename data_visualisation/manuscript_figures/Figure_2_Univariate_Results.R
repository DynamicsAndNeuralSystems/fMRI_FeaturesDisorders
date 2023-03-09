################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/Figure2/")
TAF::mkdir(plot_path)

# python_to_use <- "~/.conda/envs/pyspi/bin/python3"
python_to_use <- "/Users/abry4213/opt/anaconda3/envs/pyspi/bin/python3"
pairwise_feature_set <- "pyspi14"
univariate_feature_set <- "catch22"
data_path <- "~/data/TS_feature_manuscript"
study_group_df <- data.frame(Study = c(rep("UCLA_CNP", 3), "ABIDE_ASD"),
                             Noise_Proc = c(rep("AROMA+2P+GMR",3), "FC1000"),
                             Comparison_Group = c("Schizophrenia", "ADHD", "Bipolar", "ASD"))
univariate_feature_set <- "catch22"

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
library(knitr)
library(kableExtra)
library(patchwork)
theme_set(theme_cowplot())

# Source visualisation script
source(glue("{github_dir}/data_visualisation/manuscript_figures/Manuscript_Draft_Visualisations_Helper.R"))

# Load in univariate time-series feature info
TS_feature_info <- read.csv(glue("{github_dir}/data_visualisation/manuscript_figures/catch24_info.csv"))
# Load data
univariate_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_robustsigmoid_scaler_balanced_accuracy_all_folds.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)
univariate_SVM_coefficients <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_robustsigmoid_scaler_balanced_accuracy_all_folds.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)
univariate_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_robustsigmoid_scaler_empirical_p_values.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)
univariate_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_robustsigmoid_scaler_null_balanced_accuracy_distributions.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)

# Aggregate balanced accuracy by repeats
univariate_balanced_accuracy_by_repeats <- univariate_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Analysis_Type, group_var, Repeat_Number) %>%
  summarise(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
            Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T)) %>%
  left_join(., univariate_p_values %>% dplyr::select(Study:group_var, p_value:p_value_Bonferroni))

################################################################################
# Figure 2A univariate region-wise results
################################################################################

# Define dataset with univariate region-wise results
significant_univariate_region_wise_results <- univariate_p_values %>%
  filter(Univariate_Feature_Set == univariate_feature_set,
         Analysis_Type == "Univariate_Brain_Region") %>%
  filter(p_value_Bonferroni < 0.05) %>%
  mutate(Balanced_Accuracy_Across_Repeats = 100*Balanced_Accuracy_Across_Repeats)

# Find max fill and min fill values
min_fill <- floor(min(significant_univariate_region_wise_results$Balanced_Accuracy_Across_Repeats))
max_fill <- ceiling(max(significant_univariate_region_wise_results$Balanced_Accuracy_Across_Repeats))

# Initialize list of ggseg plots
ggseg_plot_list <- list()

# First plot within brain using ggseg
for (i in 1:nrow(study_group_df)) {
  dataset_ID <- study_group_df$Study[i]
  comparison_group <- study_group_df$Comparison_Group[i]
  
  # Define atlas by study
  atlas <- ifelse(dataset_ID == "UCLA_CNP", "dk", "hoCort")
  
  # Extract sig results to plot
  significant_data_for_ggseg <- significant_univariate_region_wise_results %>%
    filter(Study == dataset_ID,
           Comparison_Group == comparison_group) %>%
    distinct() %>%
    mutate(label = ifelse(str_detect(group_var, "ctx-"),
                          gsub("-", "_", group_var),
                          as.character(group_var))) %>%
    mutate(label = gsub("ctx_", "", label))
  
  # Plot balanced accuracy data in cortex
  dataset_ggseg <- plot_significant_regions_ggseg(dataset_ID = dataset_ID,
                                                  atlas_name = atlas,
                                                  atlas_data = get(atlas),
                                                  main_data_for_ggseg = significant_data_for_ggseg,
                                                  min_fill = min_fill,
                                                  max_fill = max_fill,
                                                  fill_color = "#F0224B")
  # Append to list
  ggseg_plot_list <- list.append(ggseg_plot_list, dataset_ggseg)
  
  # Add subcortical data for UCLA CNP
  if (dataset_ID == "UCLA_CNP") {
    dataset_ggseg_subctx <- plot_significant_regions_ggseg(dataset_ID = dataset_ID,
                                                    atlas_name = "aseg",
                                                    atlas_data = aseg,
                                                    main_data_for_ggseg = significant_data_for_ggseg,
                                                    min_fill = min_fill,
                                                    max_fill = max_fill,
                                                    fill_color = "#F0224B")
    # Append to list
    ggseg_plot_list <- list.append(ggseg_plot_list, dataset_ggseg_subctx)
  }
}

# Combine plots
wrap_plots(ggseg_plot_list, nrow=2)

for (i in 1:nrow(study_group_df)) {
  dataset_ID <- study_group_df$Study[i]
  comparison_group <- study_group_df$Comparison_Group[i]
  
  # Define atlas by study
  atlas <- ifelse(dataset_ID == "UCLA_CNP", "dk", "hoCort")
  
  significant_results <- univariate_p_values %>%
    filter(Study == dataset_ID,
           Comparison_Group == comparison_group,
           Univariate_Feature_Set == univariate_feature_set,
           Analysis_Type == "Univariate_Brain_Region") %>%
    filter(p_value_Bonferroni < 0.05)
  significant_brain_regions <- significant_results$group_var
  
  # Only move forward if 1+ significant brain regions was detected 
  if (length(significant_brain_regions) > 0) {
    # Pull out relevant null data
    null_data_to_plot <- univariate_null_distribution %>%
      filter(Study == dataset_ID,
             Comparison_Group == comparison_group,
             Univariate_Feature_Set == univariate_feature_set,
             Analysis_Type == "Univariate_Brain_Region") 
    
    # Pull out data for repeats
    repeat_data_to_plot <- univariate_balanced_accuracy_by_repeats %>%
      filter(Study == dataset_ID,
             Comparison_Group == comparison_group,
             Univariate_Feature_Set == univariate_feature_set,
             group_var %in% significant_brain_regions,
             Analysis_Type == "Univariate_Brain_Region") 
    
    ### UCLA boxplot with shaded null region
    plot_boxplot_shaded_null(dataset_ID = "UCLA_CNP",
                             grouping_var_name = "Brain Region",
                             main_data_by_repeat = repeat_data_to_plot,
                             fill_color = "#F0224B",
                             wrap_length=15,
                             null_mean_value = mean(null_data_to_plot$Null_Balanced_Accuracy, na.rm=T),
                             null_SD_value = sd(null_data_to_plot$Null_Balanced_Accuracy, na.rm=T))
    ggsave(glue("{plot_path}/{dataset_ID}_{comparison_group}_{univariate_feature_set}_Brain_Region_sig_boxplot.png"),
           width=4, height=2, units="in", dpi=300)
    
    # Plot in the brain
    significant_data_for_ggseg <- significant_results %>%
      distinct() %>%
      mutate(label = ifelse(str_detect(group_var, "ctx-"),
                            gsub("-", "_", group_var),
                            as.character(group_var))) %>%
      mutate(label = gsub("ctx_", "", label))
    plot_significant_regions_ggseg(dataset_ID = dataset_ID,
                                   atlas_name = atlas,
                                   atlas_data = get(atlas),
                                   main_data_for_ggseg = significant_data_for_ggseg,
                                   fill_color = "#F0224B")
    ggsave(glue("{plot_path}/{dataset_ID}_{comparison_group}_cortex_significant_univariate_{univariate_feature_set}.png"),
           width=4, height=2.5, units="in", dpi=300)
    
    # For UCLA_CNP, also plot subcortex
    if (dataset_ID == "UCLA_CNP") {
      plot_significant_regions_ggseg(dataset_ID = dataset_ID,
                                     atlas_name = "aseg",
                                     atlas_data = aseg,
                                     main_data_for_ggseg = significant_data_for_ggseg,
                                     fill_color = "#F0224B")
      ggsave(glue("{plot_path}/{dataset_ID}_{comparison_group}_subcortex_significant_univariate_{univariate_feature_set}.png"),
             width=2.5, height=4, units="in", dpi=300)
    }
  }
  
}

################################################################################
# Figure 2B feature-based
################################################################################

for (i in 1:nrow(study_group_df)) {
    dataset_ID <- study_group_df$Study[i]
    comparison_group <- study_group_df$Comparison_Group[i]
    
    significant_TS_features <- univariate_p_values %>%
      filter(Study == dataset_ID,
             Comparison_Group == comparison_group,
             Univariate_Feature_Set == univariate_feature_set,
             Analysis_Type == "TS_Feature") %>%
      filter(p_value_BH < 0.05) %>%
      pull(group_var)
    
    # Only move forward if 1+ significant brain regions was detected 
    if (length(significant_TS_features) > 0) {
      # Pull out relevant null data
      null_data_to_plot <- univariate_null_distribution %>%
        dplyr::rename("TS_Feature" = "group_var") %>%
        filter(Study == dataset_ID,
               Comparison_Group == comparison_group,
               Univariate_Feature_Set == univariate_feature_set,
               Analysis_Type == "TS_Feature")  %>%
        left_join(., TS_feature_info) %>%
        dplyr::rename("group_var" = "Nickname")
      
      # Pull out data for repeats
      repeat_data_to_plot <- univariate_balanced_accuracy_by_repeats %>%
        dplyr::rename("TS_Feature" = "group_var") %>%
        filter(Study == dataset_ID,
               Comparison_Group == comparison_group,
               Univariate_Feature_Set == univariate_feature_set,
               TS_Feature %in% significant_TS_features,
               Analysis_Type == "TS_Feature") %>%
        left_join(., TS_feature_info) %>%
        dplyr::rename("group_var" = "Nickname")
      
      ### UCLA boxplot with shaded null region
      plot_boxplot_shaded_null(dataset_ID = "UCLA_CNP",
                               grouping_var_name = "",
                               main_data_by_repeat = repeat_data_to_plot,
                               fill_color = "#4A7CBB",
                               wrap_length=50,
                               null_mean_value = mean(null_data_to_plot$Null_Balanced_Accuracy, na.rm=T),
                               null_SD_value = sd(null_data_to_plot$Null_Balanced_Accuracy, na.rm=T))
      ggsave(glue("{plot_path}/{dataset_ID}_{comparison_group}_{univariate_feature_set}_TS_Feature_sig_boxplot.png"),
             width=3.5+sqrt(length(significant_TS_features)), 
             height=2.75, units="in", dpi=300)
     }
}

################################################################################
# Figure 2C combo-based
################################################################################

group_order <- univariate_balanced_accuracy_by_repeats %>%
  filter(Univariate_Feature_Set==univariate_feature_set,
         Analysis_Type=="Combo") %>%
  ungroup() %>%
  mutate(is_sig = p_value_BH < 0.05,
         Comparison_Group = fct_reorder(Comparison_Group, Balanced_Accuracy_Across_Folds, .fun=mean)) %>%
  pull(Comparison_Group) %>%
  levels()

null_data <- univariate_null_distribution %>%
  filter(Univariate_Feature_Set==univariate_feature_set,
         Analysis_Type == "Combo") %>%
  group_by(Comparison_Group, Univariate_Feature_Set) %>%
  summarise(null_balacc_mean = mean(Null_Balanced_Accuracy, na.rm=T),
            null_balacc_SD = sd(Null_Balanced_Accuracy, na.rm=T)) %>%
  ungroup() %>%
  mutate(Comparison_Group = as.numeric(factor(Comparison_Group, levels = group_order)))

main_data <- univariate_balanced_accuracy_by_repeats %>%
  filter(Univariate_Feature_Set==univariate_feature_set,
         Analysis_Type=="Combo") %>%
  ungroup() %>%
  mutate(is_sig = p_value_BH < 0.05,
         Comparison_Group = as.numeric(fct_reorder(Comparison_Group, Balanced_Accuracy_Across_Folds, .fun=mean)))

ggplot() +
  geom_rect(data = null_data, 
            aes(ymin = null_balacc_mean - null_balacc_SD,
                ymax = null_balacc_mean + null_balacc_SD,
                xmin = Comparison_Group - 0.5, 
                xmax = Comparison_Group + 0.5),
            fill="gray80") +
  geom_crossbar(data = null_data, aes(x = Comparison_Group,
                                      y = null_balacc_mean,  
                                      xmin = Comparison_Group - 0.5, 
                                      xmax = Comparison_Group + 0.5)) +
  geom_boxplot(data = main_data, 
               aes(x = Comparison_Group, 
                   y = Balanced_Accuracy_Across_Folds, 
                   fill = is_sig,
                   group = Comparison_Group))  +
  scale_fill_manual(values=c("gray85", "#9B51B4")) +
  # Add labels for group names
  scale_x_continuous(breaks = 1:4,
                     labels = group_order) +
  scale_y_continuous(expand=c(0,0,0.1,0)) +
  theme(legend.position = "none") +
  coord_flip()  +
  ylab("Balanced Accuracy Across Repeats") +
  xlab("Comparison\nGroup")

ggsave(glue("{plot_path}/All_{univariate_feature_set}_Combo_boxplot.png"),
       width=5,
       height=2.25)