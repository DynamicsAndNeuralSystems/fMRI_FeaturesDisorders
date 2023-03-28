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
                             Comparison_Group = c("Schizophrenia", "Bipolar", "ADHD", "ASD"),
                             Group_Nickname = c("SCZ", "BPD", "ADHD", "ASD"))


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
library(ggseg)
library(broom)
library(colorspace)
library(see)
library(ggridges)
library(splitstackshape)
theme_set(theme_cowplot())

# Source visualisation script
source(glue("{github_dir}/data_visualisation/Manuscript_Draft_Visualisations_Helper.R"))

# Load in univariate time-series feature info
TS_feature_info <- read.csv(glue("{github_dir}/data_visualisation/catch24_info.csv"))
# Load data
univariate_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_mixedsigmoid_scaler_balanced_accuracy_all_folds.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)
univariate_SVM_coefficients <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_mixedsigmoid_scaler_SVM_coefficients.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)
univariate_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_mixedsigmoid_scaler_empirical_p_values.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)
univariate_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_mixedsigmoid_scaler_null_balanced_accuracy_distributions.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)

# Load study metadata
UCLA_CNP_metadata <- pyarrow_feather$read_feather("~/data/UCLA_CNP/study_metadata/UCLA_CNP_sample_metadata.feather") 
ABIDE_ASD_metadata <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_sample_metadata.feather") 

# Load raw feature data
UCLA_CNP_catch24 <- pyarrow_feather$read_feather("~/data/UCLA_CNP/processed_data/UCLA_CNP_AROMA_2P_GMR_catch24_filtered.feather")  %>%
  left_join(., UCLA_CNP_metadata) %>%
  mutate(label = ifelse(str_detect(Brain_Region, "ctx-"),
                        gsub("-", "_", Brain_Region),
                        as.character(Brain_Region))) %>%
  mutate(label = gsub("ctx_", "", label))
ABIDE_ASD_catch24 <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/processed_data/ABIDE_ASD_FC1000_catch24_filtered.feather")  %>%
  left_join(., ABIDE_ASD_metadata) %>%
  left_join(., ABIDE_ASD_brain_region_info) %>%
  dplyr::rename("region" = "ggseg")

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

################################################################################
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
# Plot significant results relative to their respective null distributions
null_data_for_plot <- univariate_null_distribution %>%
  filter(Analysis_Type == "Univariate_Brain_Region") %>%
  dplyr::select(Study, Comparison_Group, Null_Balanced_Accuracy) %>%
  mutate(Type = "Null",
         Null_Balanced_Accuracy = 100*Null_Balanced_Accuracy) %>%
  dplyr::rename("Balanced_Accuracy_Across_Repeats" = "Null_Balanced_Accuracy") %>%
  plyr::rbind.fill(., significant_univariate_region_wise_results %>%
                     dplyr::select(Study, Comparison_Group, Balanced_Accuracy_Across_Repeats) %>%
                     mutate(Type = "Main")) %>%
  mutate(Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                      Comparison_Group == "Bipolar" ~ "BPD",
                                      T ~ Comparison_Group)) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD")))
null_data_for_plot %>%
  ggplot(data=., mapping=aes(x=Balanced_Accuracy_Across_Repeats)) +
  geom_histogram(data = subset(null_data_for_plot, Type=="Null"),
                 fill="gray80", bins=50) +
  geom_vline(data = subset(null_data_for_plot, Type=="Main"),
             aes(xintercept = Balanced_Accuracy_Across_Repeats,
                 color = Comparison_Group),
             linewidth=0.2) +
  scale_color_manual(values=c("#573DC7", "#D5492A", "#0F9EA9", "#C47B2F")) +
  facet_wrap(Comparison_Group ~ ., ncol=2, scales="free") +
  xlab("Balanced Accuracy Across Repeats (%)") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.placement = "outside") +
  theme(legend.position="none")
ggsave(glue("{plot_path}/univariate_region_main_vs_null_balanced_acc.png"),
       width=6, height=3.5, units="in", dpi=300)

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
# Plot T statistics for Wang Periodicity in the brain
# Helper function to run t-test for given statistic
run_t_test_for_feature <- function(all_feature_values, input_feature, study, comparison_group) {
  results <- all_feature_values %>%
    filter(names==input_feature) %>%
    dplyr::select(Brain_Region, Diagnosis, values) %>%
    mutate(Diagnosis = factor(Diagnosis, levels = c(comparison_group, "Control"))) %>%
    group_by(Brain_Region) %>%
    nest() %>%
    mutate(
      fit = map(data, ~ t.test(values ~ Diagnosis, data = .x)),
      tidied = map(fit, tidy)
    ) %>% 
    unnest(tidied) %>%
    dplyr::select(-data, -fit) %>%
    arrange(p.value) %>%
    ungroup() %>%
    dplyr::select(Brain_Region, statistic) %>%
    mutate(Comparison_Group = comparison_group,
           Study = study)
  return(results)
}

# Run t-test for PD_PeriodicityWang_th0_01
periodicity_wang_Tdata_UCLA_CNP <- c("Schizophrenia", "Bipolar", "ADHD") %>%
  purrr::map_df(~ run_t_test_for_feature(all_feature_values = UCLA_CNP_catch24,
                                         comparison_group = .x,
                                         study = "UCLA_CNP",
                                         input_feature = "PD_PeriodicityWang_th0_01")) %>%
  mutate(label = ifelse(str_detect(Brain_Region, "ctx-"),
                        gsub("-", "_", Brain_Region),
                        as.character(Brain_Region))) %>%
  mutate(label = gsub("ctx_", "", label))
periodicity_wang_Tdata_ABIDE <- c("ASD") %>%
  purrr::map_df(~ run_t_test_for_feature(all_feature_values = ABIDE_ASD_catch24,
                                         comparison_group = .x,
                                         study = "ABIDE_ASD",
                                         input_feature = "PD_PeriodicityWang_th0_01")) %>%
  left_join(., ABIDE_ASD_brain_region_info) %>%
  dplyr::rename("region" = "ggseg")
periodicity_wang_Tdata_for_ggseg <- plyr::rbind.fill(periodicity_wang_Tdata_UCLA_CNP, 
                                                     periodicity_wang_Tdata_ABIDE)

ggseg_plot_list <- list()
legend_list <- list()

for (i in 1:nrow(study_group_df)) {
  dataset_ID <- study_group_df$Study[i]
  comparison_group <- study_group_df$Comparison_Group[i]
  
  # Define atlas by study
  atlas <- ifelse(dataset_ID == "UCLA_CNP", "dk", "hoCort")
  
  if (dataset_ID == "ABIDE_ASD") {
    t_stat_data <- periodicity_wang_Tdata_for_ggseg %>%
      filter(Study == dataset_ID, 
             Comparison_Group == comparison_group) %>%
      distinct() %>%
      dplyr::select(-label)
  } else {
    t_stat_data <- periodicity_wang_Tdata_for_ggseg %>%
      filter(Study == dataset_ID, 
             Comparison_Group == comparison_group) %>%
      distinct() %>%
      dplyr::select(-Index, -region)
  }
  
  # Plot T stat data in cortex
  dataset_ggseg <- plot_data_with_ggseg_diverging(dataset_ID=dataset_ID,
                                                  atlas_name=atlas,
                                                  atlas_data=get(atlas),
                                                  data_to_plot=t_stat_data,
                                                  min_fill = min(t_stat_data$statistic),
                                                  max_fill = max(t_stat_data$statistic),
                                                  fill_variable="statistic",
                                                  fill_palette="RdBu") 
  
  # Add subcortical data for UCLA CNP
  if (dataset_ID == "UCLA_CNP") {
    dataset_ggseg_subctx <- plot_data_with_ggseg_diverging(dataset_ID = dataset_ID,
                                                           atlas_name = "aseg",
                                                           atlas_data = aseg,
                                                           data_to_plot=t_stat_data,
                                                           min_fill = min(t_stat_data$statistic),
                                                           max_fill = max(t_stat_data$statistic),
                                                           fill_variable="statistic",
                                                           fill_palette="RdBu") +
      labs(fill = glue("{comparison_group}\nT-statistic"))
    
    # Extract just legend
    subctx_legend <- ggpubr::as_ggplot(ggpubr::get_legend(dataset_ggseg_subctx))
    
    # Extract just brain
    dataset_ggseg_subctx <- dataset_ggseg_subctx + 
      theme(legend.position = "none")
    
    # Append to list
    dataset_ggseg <- dataset_ggseg + theme(legend.position = "none")
    ggseg_plot_list <- list.append(ggseg_plot_list, dataset_ggseg)
    ggseg_plot_list <- list.append(ggseg_plot_list, dataset_ggseg_subctx)
    legend_list <- list.append(legend_list, subctx_legend)
  } else {
    # For ABIDE
    # Extract just legend
    ctx_legend <- ggpubr::as_ggplot(ggpubr::get_legend(dataset_ggseg + 
                                                         labs(fill = glue("{comparison_group}\nT-Statistic"))))
    
    # Append to list
    dataset_ggseg <- dataset_ggseg + theme(legend.position = "none")
    ggseg_plot_list <- list.append(ggseg_plot_list, dataset_ggseg)
    ggseg_plot_list <- list.append(ggseg_plot_list, plot_spacer())
    legend_list <- list.append(legend_list, ctx_legend)
  }
}

wrap_plots(ggseg_plot_list, 
           ncol=2, 
           byrow=T)
ggsave(glue("{plot_path}/Wang_Periodicity_T_Stats.png"),
       width=4, height=7, units="in", dpi=300)
wrap_plots(legend_list, 
           nrow=1, 
           byrow=T)
ggsave(glue("{plot_path}/Wang_Periodicity_T_Stats_legends.png"),
       width=5, height=3, units="in", dpi=300)

################################################################################
# Plot significant results relative to their respective null distributions
null_data_for_plot <- univariate_null_distribution %>%
  filter(Analysis_Type == "Univariate_TS_Feature") %>%
  dplyr::select(Study, Comparison_Group, Null_Balanced_Accuracy) %>%
  mutate(Type = "Null") %>%
  dplyr::rename("Balanced_Accuracy_Across_Repeats" = "Null_Balanced_Accuracy") %>%
  plyr::rbind.fill(., univariate_p_values %>%
                     filter(Univariate_Feature_Set == univariate_feature_set,
                            Analysis_Type == "Univariate_TS_Feature",
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
             linewidth=0.3) +
  scale_color_manual(values=c("#573DC7", "#D5492A", "#0F9EA9", "#C47B2F")) +
  facet_wrap(Comparison_Group ~ ., ncol=2, scales="free") +
  xlab("Balanced Accuracy Across Repeats (%)") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.placement = "outside") +
  theme(legend.position="none")
ggsave(glue("{plot_path}/univariate_feature_main_vs_null_balanced_acc.png"),
       width=6, height=3.5, units="in", dpi=300)

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
ggsave(glue("{plot_path}/Combo_wise_results.png"),
       width=4.5, height=3, units="in", dpi=300)

################################################################################
# Plot average magnitude of all features' coefficients per brain region
group_colors <- c("#573DC7", "#D5492A", "#0F9EA9","#C47B2F")
region_wise_avg_coefs <- univariate_SVM_coefficients %>%
  filter(Analysis_Type == "Combo") %>%
  rowwise() %>%
  mutate(Brain_Region = str_split(Feature_Name, "_", 2)[[1]][1],
         TS_Feature = str_split(Feature_Name, "_", 2)[[1]][2]) %>%
  group_by(Comparison_Group, Brain_Region) %>%
  summarise(mean_coef_magnitude = mean(abs(Coefficient), na.rm = T)) %>%
  mutate(normalised_mean_coef_magnitude = mean_coef_magnitude/max(mean_coef_magnitude)) %>%
  arrange(desc(normalised_mean_coef_magnitude)) %>%
  mutate(label = ifelse(str_detect(Brain_Region, "ctx-"),
                        gsub("-", "_", Brain_Region),
                        as.character(Brain_Region))) %>%
  mutate(label = gsub("ctx_", "", label)) %>%
  left_join(., ABIDE_ASD_brain_region_info) %>%
  dplyr::rename("region" = "ggseg")

ggseg_plot_list <- list()
legend_list <- list()

for (i in 1:nrow(study_group_df)) {
  dataset_ID <- study_group_df$Study[i]
  comparison_group <- study_group_df$Comparison_Group[i]
  group_color <- group_colors[i]
  
  # Define atlas by study
  atlas <- ifelse(dataset_ID == "UCLA_CNP", "dk", "hoCort")
  
  if (dataset_ID == "ABIDE_ASD") {
    coef_data <- region_wise_avg_coefs %>%
      filter(Comparison_Group == comparison_group) %>%
      distinct() %>%
      dplyr::select(-label)
  } else {
    coef_data <- region_wise_avg_coefs %>%
      filter(Comparison_Group == comparison_group) %>%
      distinct() %>%
      dplyr::select(-Index, -region)
  }
  
  # Plot T stat data in cortex
  dataset_ggseg <- plot_data_with_ggseg(dataset_ID=dataset_ID,
                                        atlas_name=atlas,
                                        atlas_data=get(atlas),
                                        data_to_plot=coef_data,
                                        min_fill = 0,
                                        max_fill = 1,
                                        line_color = "gray30",
                                        fill_variable="normalised_mean_coef_magnitude",
                                        fill_colors = c("white", group_color)) 
  
  # Add subcortical data for UCLA CNP
  if (dataset_ID == "UCLA_CNP") {
    dataset_ggseg_subctx <- plot_data_with_ggseg(dataset_ID = dataset_ID,
                                                 atlas_name = "aseg",
                                                 atlas_data = aseg,
                                                 data_to_plot=coef_data,
                                                 min_fill = 0,
                                                 max_fill = 1,
                                                 line_color = "gray30",
                                                 fill_variable="normalised_mean_coef_magnitude",
                                                 fill_colors = c("white", group_color)) +
      labs(fill = glue("{comparison_group}\nNormalised\nSVM Coef"))
    
    # Extract just legend
    subctx_legend <- ggpubr::as_ggplot(ggpubr::get_legend(dataset_ggseg_subctx + 
                                                            theme(legend.position = "bottom",
                                                                  legend.text = element_text(size=12)) +
                                                            labs(fill = "Normalised Absolute SVM Coefficient") +
                                                            guides(fill = guide_colorbar(title.position = "top", 
                                                                                         nrow = 1,
                                                                                         barwidth = 12, 
                                                                                         barheight = 0.75,
                                                                                         title.hjust = 0.5,
                                                                                         label.position = "bottom"))))
    
    # Extract just brain
    dataset_ggseg_subctx <- dataset_ggseg_subctx + 
      theme(legend.position = "none")
    
    # Append to list
    dataset_ggseg <- dataset_ggseg + theme(legend.position = "none")
    ggseg_plot_list <- list.append(ggseg_plot_list, dataset_ggseg)
    ggseg_plot_list <- list.append(ggseg_plot_list, dataset_ggseg_subctx)
    legend_list <- list.append(legend_list, subctx_legend)
  } else {
    # For ABIDE
    # Extract just legend
    ctx_legend <- ggpubr::as_ggplot(ggpubr::get_legend(dataset_ggseg + 
                                                         theme(legend.position = "bottom",
                                                               legend.text = element_text(size=12)) +
                                                         labs(fill = "Normalised Absolute SVM Coefficient") +
                                                         guides(fill = guide_colorbar(title.position = "top", 
                                                                                      nrow = 1,
                                                                                      barwidth = 12, 
                                                                                      barheight = 0.75,
                                                                                      title.hjust = 0.5,
                                                                                      label.position = "bottom"))))
    
    # Append to list
    dataset_ggseg <- dataset_ggseg + theme(legend.position = "none")
    ggseg_plot_list <- list.append(ggseg_plot_list, dataset_ggseg)
    ggseg_plot_list <- list.append(ggseg_plot_list, plot_spacer())
    legend_list <- list.append(legend_list, ctx_legend)
  }
}

wrap_plots(ggseg_plot_list, 
           ncol=2, 
           byrow=T)
ggsave(glue("{plot_path}/Combo_Region_Wise_SVM_Coefs.png"),
       width=4, height=7, units="in", dpi=300)
wrap_plots(legend_list, 
           nrow=4, 
           byrow=T)
ggsave(glue("{plot_path}/Combo_Region_Wise_SVM_Coefs_legend.png"),
       width=3, height=4, units="in", dpi=300)


################################################################################
# Plot average magnitude of all features' coefficients per brain region
group_colors <- c("#573DC7", "#D5492A", "#0F9EA9","#C47B2F")
feature_wise_avg_coefs <- univariate_SVM_coefficients %>%
  filter(Analysis_Type == "Combo") %>%
  rowwise() %>%
  mutate(Brain_Region = str_split(Feature_Name, "_", 2)[[1]][1],
         TS_Feature = str_split(Feature_Name, "_", 2)[[1]][2]) %>%
  group_by(Comparison_Group, TS_Feature) %>%
  summarise(mean_coef_magnitude = mean(abs(Coefficient), na.rm = T)) %>%
  mutate(normalised_mean_coef_magnitude = mean_coef_magnitude/max(mean_coef_magnitude)) %>%
  arrange(desc(normalised_mean_coef_magnitude)) 

# Annotation bar with feature type
feature_wise_avg_coefs %>%
  ungroup() %>%
  left_join(., TS_feature_info) %>%
  group_by(TS_Feature, Category) %>%
  summarise(Norm_Coef_Mean = mean(normalised_mean_coef_magnitude)) %>%
  ungroup() %>%
  mutate(TS_Feature = fct_reorder(TS_Feature, Norm_Coef_Mean),
         Category = fct_reorder(Category, Norm_Coef_Mean, .fun=mean, .desc=T)) %>%
  ggplot(data=., mapping=aes(x=0, y=TS_Feature, fill=Category)) +
  geom_tile() +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text=element_text(size=14)) +
  guides(fill = guide_legend(title.position = "top", 
                             ncol = 2,
                             byrow=T,
                             title.hjust = 0.5)) 
ggsave(glue("{plot_path}/Combo_feature_wise_SVM_coef_colorbar.png"),
       width=6, height=6, units="in", dpi=300)

# Actual heatmap
feature_wise_avg_coefs %>%
  ungroup() %>%
  left_join(., TS_feature_info) %>%
  mutate(Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                      Comparison_Group == "Bipolar" ~ "BPD",
                                      T ~ Comparison_Group)) %>%
  mutate(Figure_Name = fct_reorder(Figure_Name, normalised_mean_coef_magnitude, .fun=mean),
         Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD"))) %>%
  ggplot(data=., mapping=aes(x=Comparison_Group, y=Figure_Name, 
                             fill=normalised_mean_coef_magnitude)) +
  geom_tile()+
  geom_text(aes(label = round(normalised_mean_coef_magnitude, 1))) +
  scale_fill_gradientn(colors=c(alpha("#56C82E", 0.3), "#56C82E"), 
                       na.value=NA) +
  labs(fill = "Mean Normalised\nAbsolute SVM Coefficient") +
  xlab("Clinical Group") +
  ylab("Univariate time-series feature") +
  theme(legend.position="bottom")  +
  guides(fill = guide_colorbar(title.position = "top", 
                               nrow = 1,
                               barwidth = 14, 
                               barheight = 1,
                               title.hjust = 0.5)) 
ggsave(glue("{plot_path}/Combo_feature_wise_SVM_coef.png"),
       width=6, height=6, units="in", dpi=300)

################################################################################
# Plot significant results relative to their respective null distributions
null_data_for_plot <- univariate_null_distribution %>%
  filter(Analysis_Type == "Univariate_Combo") %>%
  dplyr::select(Study, Comparison_Group, Null_Balanced_Accuracy) %>%
  mutate(Type = "Null") %>%
  dplyr::rename("Balanced_Accuracy_Across_Repeats" = "Null_Balanced_Accuracy") %>%
  plyr::rbind.fill(., univariate_p_values %>%
                     filter(Univariate_Feature_Set == univariate_feature_set,
                            Analysis_Type == "Univariate_Combo") %>%
                     dplyr::select(Study, Comparison_Group, Balanced_Accuracy_Across_Repeats) %>%
                     mutate(Type = "Main")) %>%
  mutate(Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                      Comparison_Group == "Bipolar" ~ "BPD",
                                      T ~ Comparison_Group)) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD")))
null_data_for_plot %>%
  ggplot(data=., mapping=aes(x=100*Balanced_Accuracy_Across_Repeats)) +
  geom_histogram(data = subset(null_data_for_plot, Type=="Null"),
                 fill="gray80", bins=30) +
  geom_vline(data = subset(null_data_for_plot, Type=="Main"),
             aes(xintercept = 100*Balanced_Accuracy_Across_Repeats,
                 color = Comparison_Group),
             linewidth=1) +
  scale_color_manual(values=c("#573DC7", "#D5492A","#C47B2F", "#0F9EA9")) +
  facet_wrap(Comparison_Group ~ ., ncol=2, scales="free") +
  xlab("Balanced Accuracy Across Repeats (%)") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.placement = "outside") +
  theme(legend.position="none")
ggsave(glue("{plot_path}/univariate_combo_main_vs_null_balanced_acc.png"),
       width=4.5, height=3.5, units="in", dpi=300)

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
ggsave(glue("{plot_path}/univariate_bowtie_balanced_accuracy.png"),
       width=9, height=3, units="in", dpi=300)


################################################################################
# Summary tables for paper
################################################################################

top_coefs_overall <- univariate_SVM_coefficients %>%
  filter(Analysis_Type == "Combo") %>%
  group_by(Comparison_Group) %>%
  slice_max(order_by = abs(Coefficient), n = 10) %>%
  rowwise() %>%
  mutate(Brain_Region = str_split(Feature_Name, "_", 2)[[1]][1],
         TS_Feature = str_split(Feature_Name, "_", 2)[[1]][2])