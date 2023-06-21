################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/univariate_results/")
TAF::mkdir(plot_path)

# python_to_use <- "~/.conda/envs/pyspi/bin/python3"
python_to_use <- "/Users/abry4213/anaconda3/envs/pyspi/bin/python3"
pairwise_feature_set <- "pyspi14"
univariate_feature_set <- "catch24"
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
library(knitr)
library(kableExtra)
library(patchwork)
library(broom)
library(colorspace)
library(see)
library(ggridges)
library(splitstackshape)
library(DescTools)
library(metR)
library(scales)
theme_set(theme_cowplot())

# Source visualisation script
source(glue("{github_dir}/data_visualisation/Manuscript_Draft_Visualisations_Helper.R"))

UCLA_CNP_brain_region_info <- read.csv("~/data/UCLA_CNP/study_metadata/UCLA_CNP_Brain_Region_info.csv")
ABIDE_ASD_brain_region_info <- read.table("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_Harvard_Oxford_cort_prob_2mm_ROI_lookup.txt", sep=";", header = T) %>%
  mutate(Brain_Region = ifelse(Index==45, "Heschl's Gyrus (includes H1 and H2)", Brain_Region))

# Load in univariate time-series feature info
TS_feature_info <- read.csv(glue("{github_dir}/data_visualisation/catch24_info.csv"))
# Load data
univariate_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_mixedsigmoid_scaler_balanced_accuracy_all_folds.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)
univariate_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_mixedsigmoid_scaler_empirical_p_values.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)
univariate_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_mixedsigmoid_scaler_null_balanced_accuracy_distributions.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)

# Load study metadata
UCLA_CNP_metadata <- pyarrow_feather$read_feather("~/data/UCLA_CNP/study_metadata/UCLA_CNP_sample_metadata.feather") 
ABIDE_ASD_metadata <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_sample_metadata.feather") 

# Load t-statistics
lm_beta_stats_catch24_whole_brain <- feather::read_feather(glue("{data_path}/univariate_catch24_lm_beta_statistics_by_brain_region.feather"))

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
                                        fill_variable = "Balanced_Accuracy_Across_Repeats",
                                        fill_colors = c("#FFCA3E", "#FF6F50", "#D03454", "#9C2162", "#772F67"),
                                        bin_seq = seq(50,75,by=5),
                                        line_color = "gray30",
                                        na_color = "white") +
    labs(fill = "Mean Balanced Accuracy (%)") +
    theme(plot.title = element_blank(),
          legend.position = "bottom") +
    guides(fill = guide_colorsteps(title.position="top", ticks=TRUE, barwidth=10))
  
  # Append to list
  ggseg_plot_list <- list.append(ggseg_plot_list, dataset_ggseg)
  
  # Add subcortical data for UCLA CNP
  if (dataset_ID == "UCLA_CNP") {
    dataset_ggseg_subctx <- plot_data_with_ggseg_discrete(dataset_ID = dataset_ID,
                                                 atlas_name = "aseg",
                                                 atlas_data = aseg %>% as_tibble(),
                                                 data_to_plot = significant_data_for_ggseg,
                                                 fill_variable = "Balanced_Accuracy_Across_Repeats",
                                                 fill_colors = c("#FFCA3E", "#FF6F50", "#D03454", "#9C2162", "#772F67"),
                                                 bin_seq = seq(50,75,by=5),
                                                 line_color = "gray30",
                                                 na_color = "white") +
      labs(fill = "Mean Balanced Accuracy (%)") +
      theme(plot.title = element_blank(),
            legend.position = "bottom") +
      guides(fill = guide_colorsteps(title.position="top", ticks=TRUE, barwidth=12))
    # Append to list
    ggseg_plot_list <- list.append(ggseg_plot_list, dataset_ggseg_subctx)
  }
}

# Combine plots 
wrap_plots(ggseg_plot_list, nrow=1, widths=c(0.25, 0.12, 0.25, 0.12, 0.25, 0.12, 0.25)) + 
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom',
        legend.text = element_text(size=14),
        legend.title = element_text(size=14)) &
  guides(fill = guide_colorsteps(title.position="top", ticks=TRUE, barwidth=12))
ggsave(glue("{plot_path}/Region_wise_results.svg"),
       width=9, height=3.5, units="in", dpi=300)

# Think about hemisphere differences in performance
# Plot asymmetry index in classification performance by diagnosis group as a supplement
hemisphere_balacc_AI <- univariate_p_values %>%
  filter(Study == "UCLA_CNP", Analysis_Type == "Univariate_Brain_Region") %>%
  mutate(Hemisphere = case_when(str_detect(group_var, "Left|lh-") ~ "Left",
                                str_detect(group_var, "Right|rh-") ~ "Right")) %>%
  mutate(Brain_Region = gsub("Left-|ctx-lh-|Right-|ctx-rh-", "", group_var)) %>%
  group_by(Comparison_Group, Brain_Region) %>%
  filter(any(p_value_Bonferroni < 0.05)) %>%
  ungroup() %>% 
  pivot_wider(id_cols = c(Comparison_Group, Brain_Region),
              names_from = Hemisphere, values_from = Balanced_Accuracy_Across_Repeats) %>%
  rowwise() %>%
  mutate(Hemisphere_AI = (Left - Right) / (Left + Right)) 

# Plot regions across conditions
hemisphere_balacc_AI %>%
  left_join(., study_group_df) %>%
  mutate(Brain_Region = fct_reorder(Brain_Region, Hemisphere_AI, 
                                    .fun="sum", .desc=TRUE),
         Group_Nickname = factor(Group_Nickname, levels=c("SCZ", "BPD", "ADHD"))) %>%
  ggplot(data=., mapping=aes(x = Group_Nickname, 
                             y = Brain_Region, 
                             fill = Hemisphere_AI)) +
  geom_tile() +
  ylab("Brain region") +
  xlab("Comparison group") +
  labs(fill = "Asymmetry Index") +
  scale_fill_continuous_divergingx(palette = "PiYG", mid=0, limits=c(-0.17, 0.17)) +
  theme(legend.position = "bottom") +
  scale_y_discrete(limits=rev) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust=0.5, barwidth=10))
ggsave(glue("{plot_path}/Hemisphere_Asymmetry_Index_Balacc.svg"),
       width=4.5, height=7, units="in", dpi=300)

AI_plot_list <- list()
for (comparison_group in c("Schizophrenia", "Bipolar", "ADHD")) {
  cortex_plot <- hemisphere_balacc_AI %>% filter(Comparison_Group == comparison_group) %>%
    mutate(label = paste0("rh_", Brain_Region)) %>%
    left_join(., dk %>% as_tibble()) %>%
    filter(!is.na(region)) %>%
    dplyr::select(-region) %>%
    ggseg(atlas = "dk", mapping=aes(fill = Hemisphere_AI),
          hemisphere = "right", color="gray30") +
    scale_fill_continuous_divergingx(palette = "PiYG", mid=0, limits=c(-0.17, 0.17),
                                     na.value="white") +
    theme_void() +
    theme(legend.position = "bottom") +
    labs(fill = "Asymmetry Index") +
    guides(fill = guide_colorbar(title.position = "top", title.hjust=0.5, barwidth=10))
  
  subcortex_plot <- hemisphere_balacc_AI %>% filter(Comparison_Group == comparison_group) %>%
    mutate(label = paste0("Right-", Brain_Region)) %>%
    left_join(., aseg %>% as_tibble()) %>%
    filter(!is.na(region)) %>%
    dplyr::select(-region) %>%
    ggplot() +
    geom_brain(atlas = aseg, mapping=aes(fill = Hemisphere_AI),
          color="gray30", side = "coronal", hemi = "right") +
    scale_fill_continuous_divergingx(palette = "PiYG", mid=0, limits=c(-0.17, 0.17),
                                     na.value="white") +
    theme_void() +
    theme(legend.position = "bottom") +
    labs(fill = "Asymmetry Index") +
    guides(fill = guide_colorbar(title.position = "top", title.hjust=0.5, barwidth=10))
  AI_plot_list <- list.append(AI_plot_list, cortex_plot, subcortex_plot)
}
wrap_plots(AI_plot_list, ncol=2, widths = c(0.6, 0.4)) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom')
ggsave(glue("{plot_path}/Hemisphere_Asymmetry_Index_Balacc_in_brain.svg"),
       width=5, height=5, units="in", dpi=300)


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
ggsave(glue("{plot_path}/Feature_wise_colorbar.svg"),
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
  mutate(Figure_name = fct_reorder(Figure_name, Balanced_Accuracy_Across_Repeats, .fun=sum),
         Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD"))) %>%
  ggplot(data=., mapping=aes(x=Comparison_Group, y=Figure_name, 
                             fill=Balanced_Accuracy_Across_Repeats)) +
  geom_tile()+
  geom_text(aes(label = round(Balanced_Accuracy_Across_Repeats, 1))) +
  scale_fill_gradientn(colors=c(alpha("#4C7FC0", 0.3), "#4C7FC0"), 
                       na.value=NA) +
  labs(fill = "Mean Balanced Accuracy (%)") +
  xlab("Clinical Group") +
  ylab("Univariate time-series feature") +
  theme(legend.position="none")
ggsave(glue("{plot_path}/Feature_wise_results.svg"),
       width=4.5, height=5.5, units="in", dpi=300)


# Focus on Wang periodicity
wang_high_TS <- read_feather("~/data/TS_feature_manuscript/Wang_Periodicity_High_TS.feather")
wang_low_TS <- read_feather("~/data/TS_feature_manuscript/Wang_Periodicity_Low_TS.feather")

wang_high_TS %>% 
  filter(timepoint<=1000) %>%
  ggplot(data=., mapping=aes(x=timepoint, y=value)) +
  xlab("Time (s)") +
  ylab("Value") +
  geom_line(color="red")
ggsave(glue("{plot_path}/Wang_Periodicity_High_TS.svg"),
       width=3, height=1.25, units="in", dpi=300)
data.frame(ACF = acf(wang_high_TS %>%
  filter(timepoint<1000) %>%
  pull(value), lag.max=75, plot=F)$acf, Lag = -1:74) %>%
  filter(Lag >= 0) %>%
  ggplot(data=., mapping=aes(x=Lag, y=ACF)) +
  ylab("AC")  +
  geom_line(color="black") +
  geom_vline(xintercept = 62, color="red")
ggsave(glue("{plot_path}/Wang_Periodicity_High_ACF.svg"),
       width=3, height=1.25, units="in", dpi=300)

wang_low_TS %>% 
  filter(timepoint<150) %>%
  ggplot(data=., mapping=aes(x=timepoint, y=value)) +
  xlab("Time (s)") +
  ylab("Value") +
  geom_line(color="blue")
ggsave(glue("{plot_path}/Wang_Periodicity_Low_TS.svg"),
       width=3, height=1.25, units="in", dpi=300)
data.frame(ACF = acf(wang_low_TS %>%
                       pull(value), lag.max=75, plot=F)$acf, Lag = -1:74) %>%
  filter(Lag >= 0) %>%
  ggplot(data=., mapping=aes(x=Lag, y=ACF)) +
  ylab("AC")  +
  geom_line(color="black") +
  geom_vline(xintercept = 4, color="blue")
ggsave(glue("{plot_path}/Wang_Periodicity_Low_ACF.svg"),
       width=3, height=1.25, units="in", dpi=300)

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
      lm_beta_stat_data <- periodicity_wang_lm_beta_for_ggseg %>%
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

periodicity_wang_lm_beta_for_ggseg <- lm_beta_stats_catch24_whole_brain %>%
  filter(TS_Feature == "PD_PeriodicityWang_th0_01") %>%
  mutate(estimate = estimate) %>%
  dplyr::select(Study, Comparison_Group, TS_Feature, Brain_Region, TS_Feature, estimate, p.value)

min_fill <- floor(min(periodicity_wang_lm_beta_for_ggseg$estimate))
max_fill <- ceiling(max(periodicity_wang_lm_beta_for_ggseg$estimate))
bin_seq <- seq(min_fill, max_fill, by=1)
fill_colors = rev(c('#B2182B',"#F4A582", "white", "white", 
                    '#92C5DE', "#4393C3", "#2166AC", "#053061"))

wang_periodicity_lm_in_brain <- plot_feature_in_brain(study_group_df = study_group_df,
                                                      lm_beta_df = periodicity_wang_lm_beta_for_ggseg,
                                                      feature_name = "Wang Periodicity",
                                                      min_fill = min_fill,
                                                      max_fill = max_fill,
                                                      bin_seq = bin_seq,
                                                      fill_colors = fill_colors)

wrap_plots(wang_periodicity_lm_in_brain, 
           ncol=2, 
           byrow=T) + 
  plot_layout(guides = "collect")
ggsave(glue("{plot_path}/Wang_Periodicity_lm_beta_stats.svg"),
       width=5, height=7, units="in", dpi=300)

################################################################################
# Plot lm beta coefficients for SD + mean in the brain

SD_lm_beta_for_ggseg <- lm_beta_stats_catch24_whole_brain %>%
  filter(TS_Feature == "DN_Spread_Std") %>%
  mutate(estimate = estimate) %>%
  dplyr::select(Study, Comparison_Group, TS_Feature, Brain_Region, TS_Feature, estimate, p.value)

min_fill <- floor(min(SD_lm_beta_for_ggseg$estimate))
max_fill <- ceiling(max(SD_lm_beta_for_ggseg$estimate))
bin_seq <- seq(min_fill, max_fill, by=1)
fill_colors = rev(c("#FDDBC7", "#D1E5F0", "#92C5DE"))

SD_lm_in_brain <- plot_feature_in_brain(study_group_df = study_group_df,
                                                      lm_beta_df = SD_lm_beta_for_ggseg,
                                                      feature_name = "SD",
                                                      min_fill = min_fill,
                                                      max_fill = max_fill,
                                                      bin_seq = bin_seq,
                                                      fill_colors = fill_colors)

wrap_plots(SD_lm_in_brain, 
           ncol=2, 
           byrow=T) + 
  plot_layout(guides = "collect")
ggsave(glue("{plot_path}/SD_lm_beta_stats.svg"),
       width=5, height=7, units="in", dpi=300)

mean_lm_beta_for_ggseg <- lm_beta_stats_catch24_whole_brain %>%
  filter(TS_Feature == "DN_Mean") %>%
  mutate(estimate = estimate) %>%
  dplyr::select(Study, Comparison_Group, TS_Feature, Brain_Region, TS_Feature, estimate, p.value)

min_fill <- floor(min(mean_lm_beta_for_ggseg$estimate))
max_fill <- ceiling(max(mean_lm_beta_for_ggseg$estimate))
bin_seq <- seq(min_fill, max_fill, by=1)
fill_colors = rev(c("#FDDBC7", "#D1E5F0", "#92C5DE"))

SD_lm_in_brain <- plot_feature_in_brain(study_group_df = study_group_df,
                                        lm_beta_df = SD_lm_beta_for_ggseg,
                                        feature_name = "SD",
                                        min_fill = min_fill,
                                        max_fill = max_fill,
                                        bin_seq = bin_seq,
                                        fill_colors = fill_colors)

wrap_plots(SD_lm_in_brain, 
           ncol=2, 
           byrow=T) + 
  plot_layout(guides = "collect")
ggsave(glue("{plot_path}/SD_lm_beta_stats.svg"),
       width=5, height=7, units="in", dpi=300)

################################################################################
# Ridge plot for catch24 features' lm beta coefficients across entire brain
lm_beta_stats_catch24_whole_brain %>%
  ungroup() %>%
  left_join(., TS_feature_info) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD")))%>%
  mutate(Figure_Name = fct_reorder(Figure_Name, estimate, .fun=sd, .desc=F)) %>%
  ggplot(data=., mapping=aes(x=estimate, y=Figure_Name, fill=Comparison_Group, color=Comparison_Group)) +
  geom_density_ridges(alpha=0.6, scale=1.1, rel_min_height = 0.05) +
  xlab("Beta coefficients across\nall brain regions") +
  ylab("catch24 time-series feature") +
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
ggsave(glue("{plot_path}/catch24_feature_lm_beta_across_brain.svg"),
       width=4.75, height=10, units="in", dpi=300)

################################################################################
# Figure 2D combo-based
################################################################################
combo_null_data <- univariate_null_distribution %>%
  filter(Analysis_Type == "Univariate_Combo") %>%
  dplyr::rename("Balanced_Accuracy_Across_Folds" = Null_Balanced_Accuracy) %>%
  mutate(Balanced_Accuracy_Across_Folds = 100*Balanced_Accuracy_Across_Folds,
         Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                      Comparison_Group == "Bipolar" ~ "BPD",
                                      T ~ Comparison_Group)) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD"))) %>%
  dplyr::select(Comparison_Group, Balanced_Accuracy_Across_Folds)

univariate_balanced_accuracy_by_repeats %>%
  filter(Univariate_Feature_Set == univariate_feature_set,
         Analysis_Type == "Univariate_Combo") %>%
  mutate(Balanced_Accuracy_Across_Folds = 100*Balanced_Accuracy_Across_Folds,
         Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                      Comparison_Group == "Bipolar" ~ "BPD",
                                      T ~ Comparison_Group)) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD"))) %>%
  ggplot(data=., mapping=aes(x=Comparison_Group,
                             y=Balanced_Accuracy_Across_Folds)) +
  scale_fill_manual(values=c("#573DC7", "#D5492A", "#0F9EA9","#C47B2F")) +
  geom_violinhalf(aes(fill=Comparison_Group)) +
  geom_boxplot(width=0.1, notch=FALSE, notchwidth = 0.4, outlier.shape = NA,
               fill=NA, color=alpha("white", 0.7),
               position = position_nudge(x=0.058), coef = 0) +
  geom_violinhalf(data = combo_null_data, fill="gray80", flip=T) +
  geom_boxplot(data = combo_null_data, 
               width=0.1, notch=FALSE, 
               notchwidth = 0.4, outlier.shape = NA,
               fill=NA, color=alpha("black", 0.7),
               position = position_nudge(x=-0.058), coef = 0) +
  geom_hline(yintercept = 50, linetype=2, alpha=0.5) +
  xlab("Clinical Group") +
  ylab("Balanced Accuracy\nper Repeat (%)") +
  theme(legend.position = "none",
        axis.title = element_text(size=17), 
        axis.text = element_text(size=15)) 
ggsave(glue("{plot_path}/Combo_wise_results.svg"),
       width=4.5, height=3, units="in", dpi=300)


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
                                    Analysis_Type == "Univariate_TS_Feature" ~ "Individual\ncatch24\nFeatures",
                                    T ~ "Individual\nBrain\nRegions")) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD")),
         Analysis_Label = factor(Analysis_Label, levels = c("Individual\nBrain\nRegions",
                                                            "All\nRegions \u00D7\nFeatures",
                                                            "Individual\ncatch24\nFeatures"))) %>%
  group_by(Study, Comparison_Group, group_ID) %>%
  mutate(Comparison_Sig = paste0(Comparison_Group, "_", p_value_Bonferroni[Analysis_Type != "Univariate_Combo"]<0.05)) %>% 
  ggplot(data=., mapping=aes(x=Analysis_Label, y=Balanced_Accuracy_Across_Repeats, 
                             group=group_ID, color=Comparison_Sig)) +
  geom_line(alpha=0.5) +
  geom_point() +
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
       width=9, height=3, units="in", dpi=300)

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
  geom_violin(aes(fill = Analysis_Type), position = position_dodge(width = 0.75)) +
  xlab("Comparison Group") +
  ylab("Balanced Accuracy by Repeat") +
  scale_fill_manual(values = c("#FFD1E1", "#F05B9D")) +
  geom_boxplot(aes(fill = Analysis_Type), color="black", width=0.1, position = position_dodge(width = 0.75)) +
  theme(legend.position = "none")
ggsave(glue("{plot_path}/Univariate_Combo_with_vs_without_Regularization.svg"),
       width=9, height=3, units="in", dpi=300)
