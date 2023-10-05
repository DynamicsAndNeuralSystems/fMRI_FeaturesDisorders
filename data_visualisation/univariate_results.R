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
                             Group_Nickname = c("SCZ", "BPD", "ADHD", "ASD"))
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

# Load t-statistics
lm_beta_stats_catch25_whole_brain <- feather::read_feather(glue("{data_path}/univariate_catch25_lm_beta_statistics_by_brain_region.feather"))

################################################################################
# Univariate region-wise analysis
################################################################################

# Cortex-wise violin plots
lobe_performance <- univariate_p_values %>%
  filter(Analysis_Type == "Univariate_Brain_Region") %>%
  mutate(sig = p_value_Bonferroni<0.05) %>%
  rename("Brain_Region" = "group_var") %>%
  left_join(., plyr::rbind.fill(UCLA_CNP_brain_region_info, ABIDE_ASD_brain_region_info)) %>%
  mutate(Group_Nickname =  case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                     Comparison_Group == "Bipolar" ~ "BPD",
                                     T ~ Comparison_Group)) %>%
  mutate(Group_Nickname = factor(Group_Nickname, levels=c("SCZ", "BPD", "ADHD", "ASD"))) 

lobe_performance %>%
  filter(sig) %>%
  mutate(Cortex = as.factor(Cortex)) %>%
  ggplot(data=., mapping=aes(x=Group_Nickname, y=100*Balanced_Accuracy_Across_Folds,
                             fill=Cortex, color=Cortex)) +
  ylab("Balanced Accuracy (%)") +
  geom_violin(scale="width", position = position_dodge(width = 1), alpha=0.6) +
  geom_jitter(position = position_dodge(width = 1)) +
  theme(legend.position="none", 
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.title.y = element_blank())  +
  labs(fill="Cortical Lobe") +
  facet_grid(Group_Nickname ~ ., scales="free") +
  scale_fill_manual(values=c("Cingulate" = "#AECDE1",
                             "Frontal" = "#3C76AF",
                             "Insula" = "#BBDE93",
                             "Occipital" = "#549E3F",
                             "Parietal" = "#ED9F9C",
                             "Temporal" = "#F4C17B",
                             "Subcortex" = "#D0352B"),
                    na.value = "white") +
  scale_color_manual(values=c("Cingulate" = "#AECDE1",
                             "Frontal" = "#3C76AF",
                             "Insula" = "#BBDE93",
                             "Occipital" = "#549E3F",
                             "Parietal" = "#ED9F9C",
                             "Temporal" = "#F4C17B",
                             "Subcortex" = "#D0352B"),
                    na.value = "white") +
  coord_flip()
ggsave(glue("{plot_path}/Region_wise_violin_plot.svg"),
       width=3.5, height=6, units="in", dpi=300)
         
# Plot balanced accuracy in the brain
plot_balacc_in_brain <- function(significant_univariate_region_wise_results, 
                                 color_palette=c("#FFCA3E", "#FF6F50", "#D03454", "#9C2162", "#772F67"),
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
  filter(p_value_Bonferroni < 0.05) %>%
  mutate(Balanced_Accuracy_Across_Folds = 100*Balanced_Accuracy_Across_Folds)

catch25_region_wise_balacc_plot_list <- plot_balacc_in_brain(significant_univariate_region_wise_results,
                                                             color_palette=c("#FFCA3E", "#FF6F50", "#D03454", "#9C2162", "#772F67"),
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

################################################################################
#
# Figure 2C feature-based
#
################################################################################

# Actual heatmap
univariate_p_values %>%
  filter(Univariate_Feature_Set == univariate_feature_set,
         Analysis_Type == "Univariate_TS_Feature",
         p_value_Bonferroni < 0.05) %>%
  dplyr::rename("feature_name" = "group_var") %>%
  left_join(., TS_feature_info) %>%
  mutate(Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                      Comparison_Group == "Bipolar" ~ "BPD",
                                      T ~ Comparison_Group),
         Balanced_Accuracy_Across_Folds = 100*Balanced_Accuracy_Across_Folds) %>%
  mutate(Figure_name = fct_reorder(Figure_name, Balanced_Accuracy_Across_Folds, .fun=sum),
         Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD"))) %>%
  ggplot(data=., mapping=aes(x=Comparison_Group, y=Figure_name, 
                             fill=Balanced_Accuracy_Across_Folds)) +
  geom_tile()+
  geom_text(aes(label = round(Balanced_Accuracy_Across_Folds, 1))) +
  scale_fill_gradientn(colors=c(alpha("#4C7FC0", 0.3), "#4C7FC0"), 
                       na.value=NA) +
  labs(fill = "Mean Balanced Accuracy (%)") +
  xlab("Clinical Group") +
  ylab("Univariate time-series feature") +
  theme(legend.position="none")
ggsave(glue("{plot_path}/Feature_wise_results.svg"),
       width=4.5, height=5.5, units="in", dpi=300)

################################################################################
# Helper function to plot the beta coefficients for a given feature in the brain
plot_feature_in_brain <- function(study_group_df, lm_beta_df, feature_name, min_fill,
                                  max_fill, bin_seq, fill_colors) {
  
  ggseg_plot_list <- list()
  
  for (i in 1:nrow(study_group_df)) {
    dataset_ID <- study_group_df$Study[i]
    comparison_group <- study_group_df$Group_Nickname[i]
    
    # Define atlas by study
    atlas <- ifelse(dataset_ID == "UCLA_CNP", "dk", "hoCort")
    
    if (dataset_ID == "ABIDE_ASD") {
      lm_beta_stat_data <- lm_beta_df %>%
        filter(Study == dataset_ID, 
               Comparison_Group == comparison_group) %>%
        left_join(., ABIDE_ASD_brain_region_info) %>%
        distinct() 
    } else {
      lm_beta_stat_data <- lm_beta_df %>%
        filter(Study == dataset_ID, 
               Comparison_Group == comparison_group) %>%
        mutate(label = ifelse(str_detect(Brain_Region, "ctx-"),
                              gsub("-", "_", Brain_Region),
                              as.character(Brain_Region))) %>%
        mutate(label = gsub("ctx_", "", label)) %>%
        distinct()
    }
    
    # Plot T stat data in cortex
    dataset_ggseg <- plot_data_with_ggseg_discrete(dataset_ID = dataset_ID,
                                                   atlas_name=atlas,
                                                   atlas_data=get(atlas) %>% as_tibble(),
                                                   data_to_plot = lm_beta_stat_data,
                                                   fill_variable = "estimate",
                                                   fill_colors = fill_colors,
                                                   bin_seq = bin_seq,
                                                   line_color = "gray30",
                                                   na_color = "white")  +
      labs(fill="Beta")
    
    ggseg_plot_list <- list.append(ggseg_plot_list, dataset_ggseg)
    
    # Add subcortical data for UCLA CNP
    if (dataset_ID == "UCLA_CNP") {
      dataset_ggseg_subctx <- plot_data_with_ggseg_discrete(dataset_ID = dataset_ID,
                                                            atlas_name = "aseg",
                                                            atlas_data = aseg %>% as_tibble(),
                                                            data_to_plot=lm_beta_stat_data,
                                                            fill_variable = "estimate",
                                                            fill_colors = fill_colors,
                                                            bin_seq = bin_seq,
                                                            line_color = "gray30",
                                                            na_color = "white")  +
        labs(fill="Beta")
      
      # Append to list
      ggseg_plot_list <- list.append(ggseg_plot_list, dataset_ggseg_subctx)
    }
  }
  
  return(ggseg_plot_list)
}

# Plot lm beta statistics for Wang Periodicity in the brain

periodicity_wang_lm_beta_for_ggseg <- lm_beta_stats_catch25_whole_brain %>%
  filter(TS_Feature == "PD_PeriodicityWang_th0_01") %>%
  mutate(estimate = estimate) %>%
  dplyr::select(Study, Comparison_Group, TS_Feature, Brain_Region, TS_Feature, estimate, p.value)

wang_periodicity_lm_in_brain <- plot_feature_in_brain(study_group_df = study_group_df,
                                                      lm_beta_df = periodicity_wang_lm_beta_for_ggseg,
                                                      feature_name = "Wang Periodicity",
                                                      min_fill = round(min(periodicity_wang_lm_beta_for_ggseg$estimate),2),
                                                      max_fill = round(max(periodicity_wang_lm_beta_for_ggseg$estimate),2),
                                                      bin_seq = round(seq(min_fill, max_fill, by=0.2),1),
                                                      fill_colors = c("#053061", "#4393C3","white", "white", '#F2686D'))

wrap_plots(wang_periodicity_lm_in_brain, 
           ncol=2, 
           byrow=T) + 
  plot_layout(guides = "collect")
ggsave(glue("{plot_path}/Wang_Periodicity_lm_beta_stats.svg"),
       width=5, height=7, units="in", dpi=300)

################################################################################
# Plot lm beta coefficients for SD + mean in the brain

# SD
SD_lm_beta_for_ggseg <- lm_beta_stats_catch25_whole_brain %>%
  filter(TS_Feature == "DN_Spread_Std") %>%
  mutate(estimate = estimate) %>%
  dplyr::select(Study, Comparison_Group, TS_Feature, Brain_Region, TS_Feature, estimate, p.value)
SD_lm_in_brain <- plot_feature_in_brain(study_group_df = study_group_df,
                                                      lm_beta_df = SD_lm_beta_for_ggseg,
                                                      feature_name = "SD",
                                                      min_fill = round(min(SD_lm_beta_for_ggseg$estimate),2),
                                                      max_fill = round(max(SD_lm_beta_for_ggseg$estimate),2),
                                                      bin_seq = round(seq(min_fill, max_fill, by=0.2),1),
                                                      fill_colors = c("#053061", "#2072A6", "#85B9CE", "white", '#FF9373', "#DB4246"))

wrap_plots(SD_lm_in_brain, 
           ncol=2, 
           byrow=T) + 
  plot_layout(guides = "collect")
ggsave(glue("{plot_path}/SD_lm_beta_stats.svg"),
       width=5, height=7, units="in", dpi=300)

# fALFF
fALFF_lm_beta_for_ggseg <- lm_beta_stats_catch25_whole_brain %>%
  filter(TS_Feature == "fALFF") %>%
  mutate(estimate = estimate) %>%
  dplyr::select(Study, Comparison_Group, TS_Feature, Brain_Region, TS_Feature, estimate, p.value)
fALFF_lm_in_brain <- plot_feature_in_brain(study_group_df = study_group_df,
                                           lm_beta_df = fALFF_lm_beta_for_ggseg,
                                           feature_name = "SD",
                                           min_fill = round(min(fALFF_lm_beta_for_ggseg$estimate),2),
                                           max_fill = round(max(fALFF_lm_beta_for_ggseg$estimate),2),
                                           bin_seq = round(seq( round(min(fALFF_lm_beta_for_ggseg$estimate),2), 
                                                                round(max(fALFF_lm_beta_for_ggseg$estimate),2), 
                                                                by=0.2),1),
                                           fill_colors = c("#053061", "#2072A6", "#85B9CE", "white", '#FF9373'))

wrap_plots(fALFF_lm_in_brain, 
           ncol=2, 
           byrow=T) + 
  plot_layout(guides = "collect")
ggsave(glue("{plot_path}/fALFF_lm_beta_stats.svg"),
       width=5, height=7, units="in", dpi=300)


################################################################################
# Ridge plot for catch25 features' lm beta coefficients across entire brain
lm_beta_stats_catch25_whole_brain %>%
  ungroup() %>%
  left_join(., TS_feature_info) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD")))%>%
  mutate(Figure_Name = fct_reorder(Figure_Name, estimate, .fun=sd, .desc=F)) %>%
  ggplot(data=., mapping=aes(x=estimate, y=Figure_Name, fill=Comparison_Group, color=Comparison_Group)) +
  geom_density_ridges(alpha=0.6, scale=1.1, rel_min_height = 0.05) +
  xlab("Beta coefficients across\nall brain regions") +
  ylab("catch25 time-series feature") +
  scale_fill_manual(values=c("Control" = "#5BB67B", 
                             "SCZ" = "#573DC7", 
                             "BPD" = "#D5492A", 
                             "ADHD" = "#0F9EA9", 
                             "ASD" = "#C47B2F")) +
  scale_color_manual(values=c("Control" = "#5BB67B", 
                              "SCZ" = "#573DC7", 
                              "BPD" = "#D5492A", 
                              "ADHD" = "#0F9EA9", 
                              "ASD" = "#C47B2F")) +
  guides(fill = guide_legend(nrow=2, byrow=T),
         color = guide_legend(nrow=2, byrow=T)) +
  scale_y_discrete(labels = wrap_format(28)) +
  scale_x_continuous(limits=c(-1,1)) +
  theme(legend.position = "bottom",
        axis.title = element_text(size=16),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=15),
        legend.text = element_text(size=16),
        legend.title = element_blank())
ggsave(glue("{plot_path}/catch25_feature_lm_beta_across_brain.svg"),
       width=4.75, height=10, units="in", dpi=300)


################################################################################
# Figure 2D combo-based
################################################################################

# Bowtie plot comparing each brain region and feature to combo-wise performance
num_comparisons <- univariate_p_values %>%
  filter(Analysis_Type != "Univariate_Combo") %>%
  group_by(Study, Comparison_Group) %>%
  count()

univariate_p_values %>%
  filter(Analysis_Type == "Univariate_Combo") %>% 
  left_join(., num_comparisons) %>%
  expandRows("n") %>%
  group_by(Study, Comparison_Group) %>%
  mutate(group_ID = paste0(Comparison_Group, "_", row_number())) %>%
  plyr::rbind.fill(., univariate_p_values %>%
                     filter(Analysis_Type != "Univariate_Combo") %>%
                     group_by(Study, Comparison_Group) %>%
                     mutate(group_ID = paste0(Comparison_Group, "_", row_number()))) %>%
  rowwise() %>%
  mutate(Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                      Comparison_Group == "Bipolar" ~ "BPD",
                                      T ~ Comparison_Group),
         Analysis_Label = case_when(Analysis_Type == "Univariate_Combo" ~ "All\nRegions \u00D7\nFeatures",
                                    Analysis_Type == "Univariate_TS_Feature" ~ "Individual\ncatch25\nFeatures",
                                    T ~ "Individual\nBrain\nRegions")) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD")),
         Analysis_Label = factor(Analysis_Label, levels = c("Individual\ncatch25\nFeatures",
                                                            "All\nRegions \u00D7\nFeatures",
                                                            "Individual\nBrain\nRegions"))) %>%
  group_by(Study, Comparison_Group, group_ID) %>%
  mutate(Comparison_Sig = paste0(Comparison_Group, "_", p_value_Bonferroni[Analysis_Type != "Univariate_Combo"]<0.05)) %>% 
  ggplot(data=., mapping=aes(x=Analysis_Label, y=100*Balanced_Accuracy_Across_Folds, 
                             group=group_ID, color=Comparison_Sig)) +
  geom_line(alpha=0.5) +
  geom_point(size=1.25) +
  ylab("Mean Balanced Accuracy Across Repeats (%)") +
  xlab("Linear SVM Type") +
  scale_x_discrete(expand=c(0.05,0.05,0.05,0.05)) +
  scale_color_manual(values = c("SCZ_TRUE"="#573DC7", 
                                "BPD_TRUE"="#D5492A", 
                                "ADHD_TRUE"="#0F9EA9",
                                "ASD_TRUE"="#C47B2F",
                                "SCZ_FALSE"="gray70", 
                                "BPD_FALSE"="gray70", 
                                "ADHD_FALSE"="gray70",
                                "ASD_FALSE"="gray70")) +
  facet_wrap(Comparison_Group ~ ., scales="free_x", ncol=4) +
  theme(legend.position = "none",
        plot.margin = margin(1,10,1,1)) + 
  coord_flip()
ggsave(glue("{plot_path}/univariate_bowtie_balanced_accuracy.svg"),
       width=9, height=2.5, units="in", dpi=300)

################################################################################
# Compare performance with L1-regularized SVM balanced accuracy implemented with LinearSVC()
SVM_L1_balanced_accuracy_by_folds <- pyarrow_feather$read_feather(glue("{data_path}/SVM_L1_Regularized_Balanced_Accuracy.feather"))
# Aggregate balanced accuracy by repeats
SVM_L1_balanced_accuracy_by_repeats <- SVM_L1_balanced_accuracy_by_folds %>%
  group_by(Study, Comparison_Group, Analysis_Type, Repeat_Number) %>%
  summarise(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
            Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T)) 
SVM_L1_balanced_accuracy <- SVM_L1_balanced_accuracy_by_repeats %>%
  group_by(Study, Comparison_Group, Analysis_Type) %>%
  summarise(Balanced_Accuracy_Across_Repeats = mean(Balanced_Accuracy_Across_Folds, na.rm=T),
            Balanced_Accuracy_Across_Repeats_SD = sd(Balanced_Accuracy_Across_Folds, na.rm=T)) 

SVM_L1_regularized_coefficients_by_folds <- pyarrow_feather$read_feather(glue("{data_path}/SVM_L1_Regularized_Coefficients.feather")) %>%
  dplyr::rename("Feature_Name" = "Feature Name", "Repeat_Number" = "Repeat Number")
SVM_L1_regularized_coefficients_by_repeats <- SVM_L1_regularized_coefficients_by_folds %>%
  group_by(Study, Comparison_Group, Analysis_Type, Feature_Name, Repeat_Number) %>%
  summarise(Coefficient_Across_Folds = mean(Coefficient, na.rm=T),
            Coefficient_Across_Folds_SD = sd(Coefficient, na.rm=T)) 
SVM_L1_regularized_coefficients <- SVM_L1_regularized_coefficients_by_repeats %>%
  group_by(Study, Comparison_Group, Analysis_Type, Feature_Name) %>%
  summarise(Coefficient_Across_Repeats = mean(Coefficient_Across_Folds, na.rm=T),
            Coefficient_Across_Repeats_SD = sd(Coefficient_Across_Folds, na.rm=T)) 

# Find number of zero vs non-zero components per group
SVM_L1_regularized_coefficients %>%
  group_by(Study, Comparison_Group) %>%
  summarise(num_zero = sum(Coefficient_Across_Repeats == 0),
            num_nonzero = sum(Coefficient_Across_Repeats != 0),
            total = n())

# Plot performance in normal SVM vs. L1-regularized SVM
SVM_L1_balanced_accuracy_by_repeats %>%
  dplyr::select(-Balanced_Accuracy_Across_Folds_SD) %>%
  plyr::rbind.fill(univariate_balanced_accuracy_by_repeats %>% filter(Analysis_Type=="Univariate_Combo")) %>%
  left_join(study_group_df) %>%
  mutate(Group_Nickname = factor(Group_Nickname, levels=c("SCZ", "BPD", "ADHD", "ASD"))) %>%
  ggplot(data=., mapping=aes(x=Group_Nickname, y=Balanced_Accuracy_Across_Folds)) +
  geom_violin(aes(fill = Analysis_Type), position = position_dodge(width = 0.75),
              scale="width", width=0.6) +
  xlab("Comparison Group") +
  ylab("Balanced Accuracy\nby Repeat") +
  scale_fill_manual(values = c("cadetblue2", "darkcyan")) +
  geom_boxplot(aes(fill = Analysis_Type), color="black", width=0.1, position = position_dodge(width = 0.75)) +
  theme(legend.position = "none")
ggsave(glue("{plot_path}/Univariate_Combo_with_vs_without_Regularization.svg"),
       width=4, height=2.5, units="in", dpi=300)
