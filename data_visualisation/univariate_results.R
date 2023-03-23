################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/univariate_results/")
TAF::mkdir(plot_path)

# python_to_use <- "~/.conda/envs/pyspi/bin/python3"
python_to_use <- "/Users/abry4213/opt/anaconda3/envs/pyspi/bin/python3"
pairwise_feature_set <- "pyspi14"
univariate_feature_set <- "catch24"
data_path <- "~/data/TS_feature_manuscript"
study_group_df <- data.frame(Study = c(rep("UCLA_CNP", 3), "ABIDE_ASD"),
                             Noise_Proc = c(rep("AROMA+2P+GMR",3), "FC1000"),
                             Comparison_Group = c("Schizophrenia", "Bipolar", "ADHD", "ASD"))

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
source(glue("{github_dir}/data_visualisation/Manuscript_Draft_Visualisations_Helper.R"))

# Load in univariate time-series feature info
TS_feature_info <- read.csv(glue("{github_dir}/data_visualisation/catch24_info.csv"))
# Load data
univariate_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_mixedsigmoid_scaler_balanced_accuracy_all_folds.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)
univariate_SVM_coefficients <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_mixedsigmoid_scaler_balanced_accuracy_all_folds.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)
univariate_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_mixedsigmoid_scaler_empirical_p_values.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)
univariate_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_mixedsigmoid_scaler_null_balanced_accuracy_distributions.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)

# Aggregate balanced accuracy by repeats
univariate_balanced_accuracy_by_repeats <- univariate_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Analysis_Type, group_var, Repeat_Number) %>%
  summarise(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
            Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T)) %>%
  left_join(., univariate_p_values %>% dplyr::select(Study:group_var, p_value:p_value_Bonferroni))

################################################################################
# Figure 2B univariate region-wise results
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
  
  # If dataset is ABIDE ASD, convert regions to ggseg regions
  if (dataset_ID == "ABIDE_ASD") {
    significant_data_for_ggseg <- significant_univariate_region_wise_results %>%
      filter(Study == dataset_ID,
             Comparison_Group == comparison_group) %>%
      dplyr::rename("Brain_Region" = "group_var") %>%
      left_join(., ABIDE_ASD_brain_region_info) %>%
      dplyr::rename("region" = "ggseg")
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
wrap_plots(ggseg_plot_list, nrow=2, widths = c(0.25, 0.1, 0.25, 0.1)) + 
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom')
ggsave(glue("{plot_path}/Region_wise_results.png"),
       width=8, height=6, units="in", dpi=300)

################################################################################
# Figure 2C feature-based
################################################################################

# Annotation bar with feature type
univariate_p_values %>%
  filter(Univariate_Feature_Set == univariate_feature_set,
         Analysis_Type == "Univariate_TS_Feature",
         p_value_Bonferroni < 0.05) %>%
  dplyr::rename("TS_Feature" = "group_var") %>%
  left_join(., TS_feature_info) %>%
  group_by(TS_Feature, Category) %>%
  summarise(Balacc_Sum = sum(Balanced_Accuracy_Across_Repeats)) %>%
  ungroup() %>%
  mutate(TS_Feature = fct_reorder(TS_Feature, Balacc_Sum),
         Category = fct_reorder(Category, Balacc_Sum, .fun=sum, .desc=T)) %>%
  ggplot(data=., mapping=aes(x=0, y=TS_Feature, fill=Category)) +
  geom_tile() +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text=element_text(size=14)) +
  guides(fill = guide_legend(title.position = "top", 
                             ncol = 2,
                             byrow=T,
                             title.hjust = 0.5)) 
ggsave(glue("{plot_path}/Feature_wise_colorbar.png"),
       width=6, height=6, units="in", dpi=300)

# Actual heatmap
univariate_p_values %>%
  filter(Univariate_Feature_Set == univariate_feature_set,
         Analysis_Type == "Univariate_TS_Feature",
         p_value_Bonferroni < 0.05) %>%
  dplyr::rename("TS_Feature" = "group_var") %>%
  left_join(., TS_feature_info) %>%
  mutate(Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                      Comparison_Group == "Bipolar" ~ "BPD",
                                      T ~ Comparison_Group),
         Balanced_Accuracy_Across_Repeats = 100*Balanced_Accuracy_Across_Repeats) %>%
  mutate(Figure_Name = fct_reorder(Figure_Name, Balanced_Accuracy_Across_Repeats, .fun=sum),
         Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD"))) %>%
  ggplot(data=., mapping=aes(x=Comparison_Group, y=Figure_Name, 
                             fill=Balanced_Accuracy_Across_Repeats)) +
  geom_tile()+
  geom_text(aes(label = round(Balanced_Accuracy_Across_Repeats, 1))) +
  scale_fill_gradientn(colors=c(alpha("#4C7FC0", 0.3), "#4C7FC0"), 
                       na.value=NA) +
  labs(fill = "Mean Balanced Accuracy (%)") +
  xlab("Clinical Group") +
  ylab("Univariate time-series feature") +
  theme(legend.position="bottom")  +
  guides(fill = guide_colorbar(title.position = "top", 
                               nrow = 1,
                               barwidth = 12, 
                               barheight = 1,
                               title.hjust = 0.5)) 
ggsave(glue("{plot_path}/Feature_wise_results.png"),
       width=5.5, height=5.5, units="in", dpi=300)

################################################################################
# Figure 2D combo-based
################################################################################

univariate_balanced_accuracy_by_repeats %>%
  filter(Univariate_Feature_Set == univariate_feature_set,
         Analysis_Type == "Univariate_Combo") %>%
  mutate(sig_fill = p_value_Bonferroni < 0.05,
         Balanced_Accuracy_Across_Folds = 100*Balanced_Accuracy_Across_Folds,
         Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                      Comparison_Group == "Bipolar" ~ "BPD",
                                      T ~ Comparison_Group)) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = rev(c("SCZ", "BPD", "ASD", "ADHD")))) %>%
  ggplot(data=., mapping=aes(x=Balanced_Accuracy_Across_Folds,
                             y=Comparison_Group,
                             fill=sig_fill)) +
  geom_boxplot() +
  scale_fill_manual(values=c("gray85", "#9B51B4")) +
  ylab("Clinical Group") +
  xlab("Balanced Accuracy per Repeat (%)") +
  theme(legend.position = "none",
        axis.title.y = element_text(size=17),
        axis.text.y = element_text(size=15))
ggsave(glue("{plot_path}/Combo_wise_results.png"),
       width=4.5, height=5.5, units="in", dpi=300)
