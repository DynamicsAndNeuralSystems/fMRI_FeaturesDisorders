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
study_group_df <- data.frame(Study = rep("UCLA_CNP", 3),
                             Noise_Proc = rep("AROMA+2P+GMR", 3),
                             Comparison_Group = c("Schizophrenia", "ADHD", "Bipolar"),
                             Group_Nickname = c("SCZ", "BPD", "ADHD"))
# study_group_df <- data.frame(Study = c(rep("UCLA_CNP", 3), "ABIDE_ASD"),
#                              Noise_Proc = c(rep("AROMA+2P+GMR",3), "FC1000"),
#                              Comparison_Group = c("Schizophrenia", "ADHD", "Bipolar", "ASD"),
# Group_Nickname = c("SCZ", "BPD", "ADHD", "ASD"))

# Load brain region info
UCLA_CNP_brain_region_info <- read.csv("~/data/UCLA_CNP/study_metadata/UCLA_CNP_Brain_Region_info.csv")
ABIDE_ASD_brain_region_info <- read.csv("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv")
region_node_to_from <- read.csv("~/data/TS_feature_manuscript/node_to_from_structure.csv")

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
library(ggnewscale)
library(scales)
library(splitstackshape)
library(LaCroixColoR)
library(correctR)
library(ggpubr)
library(ggsignif)
library(poolr)
library(ggraph)
library(igraph)
theme_set(theme_cowplot())

# Source visualisation script
source(glue("{github_dir}/data_visualisation/Manuscript_Draft_Visualisations_Helper.R"))

# Load in SPI info
SPI_info <- read.csv(glue("{github_dir}/data_visualisation/SPI_info.csv"))

# Load participants included
UCLA_CNP_subjects_to_keep <- pyarrow_feather$read_feather("~/data/UCLA_CNP/processed_data/UCLA_CNP_filtered_sample_info_AROMA_2P_GMR_catch24_pyspi14.feather")
  
# Load study metadata
UCLA_CNP_metadata <- pyarrow_feather$read_feather("~/data/UCLA_CNP/study_metadata/UCLA_CNP_sample_metadata.feather") %>%
  mutate(Study = "UCLA_CNP") %>%
  filter(Sample_ID %in% UCLA_CNP_subjects_to_keep$Sample_ID)
ABIDE_ASD_metadata <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_sample_metadata.feather") %>%
  mutate(Study = "ABIDE")

# Load stats data
univariate_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_mixedsigmoid_scaler_balanced_accuracy_all_folds.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)
univariate_balanced_accuracy <- univariate_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Analysis_Type, group_var, Repeat_Number) %>%
  reframe(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
          Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T)) %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Analysis_Type, group_var) %>%
  reframe(Balanced_Accuracy_Across_Repeats = mean(Balanced_Accuracy_Across_Folds, na.rm=T),
          Balanced_Accuracy_Across_Repeats_SD = sd(Balanced_Accuracy_Across_Folds, na.rm=T))

pairwise_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_mixedsigmoid_scaler_balanced_accuracy_all_folds.feather"))
pairwise_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_mixedsigmoid_scaler_empirical_p_values.feather"))
pairwise_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_mixedsigmoid_scaler_null_balanced_accuracy_distributions.feather"))
univariate_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_mixedsigmoid_scaler_null_balanced_accuracy_distributions.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)

combo_univariate_pairwise_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combined_univariate_pairwise_mixedsigmoid_scaler_balanced_accuracy_all_folds.feather")) %>%
  mutate(Analysis_Type = "SPI_Univariate_Combo")
combo_univariate_pairwise_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combined_univariate_pairwise_mixedsigmoid_scaler_empirical_p_values.feather"))
combo_univariate_pairwise_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combined_univariate_pairwise_mixedsigmoid_scaler_null_balanced_accuracy_distributions.feather"))

pairwise_all_SPIs_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_all_SPIs_mixedsigmoid_scaler_balanced_accuracy_all_folds.feather"))

# Aggregate the main results across folds and then across repeats
# Pairwise SPI-wise
pairwise_balanced_accuracy <- pairwise_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Pairwise_Feature_Set, Analysis_Type, group_var, Repeat_Number) %>%
  reframe(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
          Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T)) %>%
  group_by(Study, Comparison_Group, Pairwise_Feature_Set, Analysis_Type, group_var) %>%
  reframe(Balanced_Accuracy_Across_Repeats = mean(Balanced_Accuracy_Across_Folds, na.rm=T),
          Balanced_Accuracy_Across_Repeats_SD = sd(Balanced_Accuracy_Across_Folds, na.rm=T))
pairwise_balanced_accuracy_by_repeats <- pairwise_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Pairwise_Feature_Set, Analysis_Type, group_var, Repeat_Number) %>%
  summarise(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
            Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T)) %>%
  left_join(., pairwise_p_values %>% dplyr::select(Study:group_var, p_value:p_value_Bonferroni))

# Pairwise SPI-wise with univariate combo data
combo_univariate_pairwise_balanced_accuracy <- combo_univariate_pairwise_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Pairwise_Feature_Set, Analysis_Type, group_var, Repeat_Number) %>%
  reframe(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
          Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T)) %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Pairwise_Feature_Set, Analysis_Type, group_var) %>%
  reframe(Balanced_Accuracy_Across_Repeats = mean(Balanced_Accuracy_Across_Folds, na.rm=T),
          Balanced_Accuracy_Across_Repeats_SD = sd(Balanced_Accuracy_Across_Folds, na.rm=T))
combo_univariate_pairwise_balanced_accuracy_by_repeats <- combo_univariate_pairwise_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Pairwise_Feature_Set, Analysis_Type, group_var, Repeat_Number) %>%
  summarise(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
            Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T)) %>%
  left_join(., combo_univariate_pairwise_p_values %>% dplyr::select(Study:group_var, p_value:p_value_Bonferroni))

# Pairwise all SPIs together
pairwise_all_SPIs_balanced_accuracy <- pairwise_all_SPIs_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Pairwise_Feature_Set, Analysis_Type, group_var, Repeat_Number) %>%
  reframe(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
          Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T)) %>%
  group_by(Study, Comparison_Group, Pairwise_Feature_Set, Analysis_Type, group_var) %>%
  reframe(Balanced_Accuracy_Across_Repeats = mean(Balanced_Accuracy_Across_Folds, na.rm=T),
          Balanced_Accuracy_Across_Repeats_SD = sd(Balanced_Accuracy_Across_Folds, na.rm=T))
pairwise_all_SPIs_balanced_accuracy_by_repeats <- pairwise_all_SPIs_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Pairwise_Feature_Set, Analysis_Type, group_var, Repeat_Number) %>%
  summarise(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
            Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T)) %>%
  left_join(., pairwise_p_values %>% dplyr::select(Study:group_var, p_value:p_value_Bonferroni))


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
  theme(legend.position="none")
ggsave(glue("{plot_path}/SPI_wise_results.png"),
       width=5, height=4.5, units="in", dpi=300)


################################################################################
# Compare SPIs with vs without univariate info
################################################################################

# Use correctR to test for difference across resamples for FTM vs catch22+FTM
run_correctR_group <- function(comparison_group, study, metadata, results_df) {
  num_subjects <- metadata %>%
    filter(Study == study, 
           Diagnosis %in% c("Control", comparison_group)) %>%
    distinct(Sample_ID) %>%
    nrow()
  training_size <- ceiling(0.9*num_subjects)
  test_size <- floor(0.1*num_subjects)
  
  data_for_correctR <- results_df %>%
    filter(Study == study, 
           Comparison_Group == comparison_group) %>%
    group_by(group_var) %>%
    filter(any(p_value_Bonferroni < 0.05)) %>%
    ungroup() %>%
    pivot_wider(id_cols = c(Repeat_Number, group_var), 
                names_from = Analysis_Type,
                values_from = Balanced_Accuracy_Across_Folds) %>%
    dplyr::rename("x" = "Pairwise_SPI", "y" = "SPI_Univariate_Combo") %>%
    group_by(group_var) %>%
    group_split()
  
  res <- data_for_correctR %>%
    purrr::map_df(~ as.data.frame(resampled_ttest(x=.x$x, 
                                    y=.x$y, 
                                    n=10, 
                                    n1=training_size, n2=test_size)) %>%
                    mutate(SPI = unique(.x$group_var))) %>%
    ungroup() %>%
    mutate(p_value_Bonferroni = p.adjust(p.value, method="bonferroni"),
           Comparison_Group = comparison_group)
  
  return(res)
  
}

metadata <- UCLA_CNP_metadata
# metadata <- plyr::rbind.fill(UCLA_CNP_metadata, ABIDE_ASD_metadata)
results_df = plyr::rbind.fill(pairwise_balanced_accuracy_by_repeats, 
                              combo_univariate_pairwise_balanced_accuracy_by_repeats)

corrected_SPI_T_res <- 1:nrow(study_group_df) %>%
  purrr::map_df(~ run_correctR_group(comparison_group = study_group_df$Comparison_Group[.x],
                                     study = study_group_df$Study[.x],
                                     metadata = metadata,
                                     results_df = results_df)) %>%
  filter(p_value_Bonferroni < 0.05)

plyr::rbind.fill(pairwise_p_values,
                 combo_univariate_pairwise_p_values) %>%
  dplyr::rename("SPI" = group_var) %>%
  semi_join(., corrected_SPI_T_res %>% dplyr::select(SPI, Comparison_Group)) %>%
  mutate(Analysis_Type = ifelse(Analysis_Type == "Pairwise_SPI", "SPI\nOnly", "SPI + Univariate\nRegion ×\nFeature"),
         sig = ifelse(p_value_Bonferroni < 0.05, "Significant", "Not significant")) %>%
  mutate(Analysis_Type = factor(Analysis_Type, levels=c("SPI\nOnly", "SPI + Univariate\nRegion ×\nFeature"))) %>%
  mutate(Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                      Comparison_Group == "Bipolar" ~ "BPD",
                                      T ~ Comparison_Group)) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD"))) %>%
  ggplot(data=., mapping=aes(x=Analysis_Type, y=Balanced_Accuracy_Across_Repeats,
                             group = SPI)) +
  geom_line(aes(color = Comparison_Group), show.legend = FALSE) +
  scale_color_manual(values=c("Control" = "#5BB67B", 
                              "SCZ" = "#573DC7", 
                              "BPD" = "#D5492A", 
                              "ADHD" = "#0F9EA9", 
                              "ASD" = "#C47B2F")) +
  new_scale_colour() +  # start a new scale
  geom_point(aes(color = sig)) +
  scale_color_manual(values = c("gray60", "#5BB67B")) +
  # scale_x_discrete(labels = wrap_format(7)) +
  xlab("Analysis Type") +
  ylab("Mean Balanced Accuracy (%)") +
  facet_wrap(Comparison_Group ~ ., ncol=2, scales="fixed") +
  scale_x_discrete(expand=c(0,0.2,0,0.2)) +
  theme(legend.position = "bottom",
        plot.margin = margin(1,30,1,1, unit="pt"),
        legend.title = element_blank())
ggsave(glue("{plot_path}/SPI_with_vs_without_univariate_spaghetti.png"),
       width=4, height=4.5, units="in", dpi=300)

################################################################################
# Distribution of SPI-wise T-statistics
################################################################################
t_stats_pyspi14_whole_brain <- feather::read_feather(glue("{data_path}/pairwise_pyspi14_t_statistics_by_region_pair.feather"))

t_stats_pyspi14_whole_brain %>%
  ungroup() %>%
  left_join(., SPI_info) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD")))%>%
  mutate(Nickname = fct_reorder(Nickname, statistic, .fun=sd)) %>%
  ggplot(data=., mapping=aes(x=statistic, y=Nickname, fill=Comparison_Group, color=Comparison_Group)) +
  geom_density_ridges(alpha=0.6, scale=1.1) +
  xlab("T-statistic across\nall brain regions") +
  ylab("pyspi14 time-series feature") +
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
  guides(fill = guide_legend(nrow=2),
         color = guide_legend(nrow=2)) +
  scale_y_discrete(labels = wrap_format(22)) +
  theme(legend.position = "bottom",
        axis.title = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        legend.text = element_text(size=16),
        legend.title = element_blank())
ggsave(glue("{plot_path}/pyspi14_feature_t_statistics_across_brain.png"),
       width=5.5, height=10, units="in", dpi=300)


################################################################################
# Condense T-statistics down to the regional level
################################################################################

# Demo brain figure
dk %>%
  as_tibble() %>%
  mutate(fillval = case_when(label == "lh_caudalmiddlefrontal" ~ "1",
                             label == "lh_bankssts" ~ "2",
                             label == "lh_lateralorbitofrontal" ~ "3", 
                             label == "lh_superiorparietal" ~ "4",
                             label == "lh_lateraloccipital" ~ "5",
                             T ~ NA_character_)) %>%
  ggseg(atlas = "dk", mapping = aes(fill = fillval),
        hemisphere="left",
        view = "lateral",
        position = "stacked", colour = "gray50") +
  scale_fill_manual(values=lacroix_palette("PassionFruit", n=5), na.value="white") +
  theme_void() +
  theme(plot.title = element_blank(),
        legend.position = "none")
ggsave(glue("{plot_path}/demo_brain_for_FC.png"), width = 3, height=2, units="in", dpi=300)

# Find regions most disrupted across all pairwise connections
pairwise_t_stats_by_region_from <- t_stats_pyspi14_whole_brain %>%
  separate(Region_Pair, c("region_from", "region_to"),
           sep="_") %>%
  group_by(region_from, SPI, Study, Comparison_Group) %>%
  summarise(mean_T_magnitude = mean(abs(statistic))) %>%
  dplyr::rename("Brain_Region" = "region_from") %>%
  mutate(Direction = "from")

pairwise_t_stats_by_region_to <- t_stats_pyspi14_whole_brain %>%
  separate(Region_Pair, c("region_from", "region_to"),
           sep="_") %>%
  group_by(region_to, SPI, Study, Comparison_Group) %>%
  summarise(mean_T_magnitude = mean(abs(statistic)))%>%
  dplyr::rename("Brain_Region" = "region_to") %>%
  mutate(Direction = "to")

# Function to plot T-statistics in the brain for a given SPI
plot_data_per_group <- function(SPI_nickname, 
                                bin_seq,
                                input_data, 
                                group_colors,
                                bin_interval) { 
  ggseg_plot_list <- list()
  
  color_range = c(plyr::round_any(min(input_data$mean_T_magnitude), bin_interval, f=floor), 
                  plyr::round_any(max(input_data$mean_T_magnitude), bin_interval, f=ceiling))
  color_seq = seq(color_range[1], color_range[2], by=bin_interval)
  
  for (i in 1:nrow(study_group_df)) {
    dataset_ID <- study_group_df$Study[i]
    comparison_group <- study_group_df$Group_Nickname[i]
    group_color <- group_colors[i]
    
    # Define atlas by study
    atlas <- ifelse(dataset_ID == "UCLA_CNP", "dk", "hoCort")
    
    T_data_to_plot <- input_data %>%
      filter(Study == dataset_ID,
             Comparison_Group == comparison_group) %>%
      select(where(function(x) any(!is.na(x)))) %>%
      arrange(desc(mean_T_magnitude)) %>%
      ungroup() 
    
    # Plot T stat data in cortex
    dataset_ggseg <- plot_data_with_ggseg_discrete(dataset_ID=dataset_ID,
                                                   atlas_name=atlas,
                                                   atlas_data=get(atlas),
                                                   data_to_plot=T_data_to_plot,
                                                   num_bins = length(color_seq) - 1,
                                                   line_color = "gray30",
                                                   fill_variable="mean_T_magnitude",
                                                   bin_seq = color_seq,
                                                   fill_colors = colorRampPalette(c("white", group_color))(length(color_seq) - 1)) 
    
    # Append to list
    ggseg_plot_list <- list.append(ggseg_plot_list, dataset_ggseg)
    
    # Add subcortical data for UCLA CNP
    if (dataset_ID == "UCLA_CNP") {
      dataset_ggseg_subctx <- plot_data_with_ggseg_discrete(dataset_ID = dataset_ID,
                                                            atlas_name = "aseg",
                                                            atlas_data = aseg,
                                                            data_to_plot=T_data_to_plot,
                                                            num_bins = length(color_seq) - 1,
                                                            line_color = "gray30",
                                                            fill_variable="mean_T_magnitude",
                                                            bin_seq = color_seq,
                                                            fill_colors = colorRampPalette(c("white", group_color))(length(color_seq) - 1))
      # Append to list
      ggseg_plot_list <- list.append(ggseg_plot_list, dataset_ggseg_subctx)
    }
  }
  return(ggseg_plot_list)
}

group_colors <- c("#573DC7", "#D5492A", "#0F9EA9")
pearson_regional_data <- pairwise_t_stats_by_region_from %>%
  filter(SPI == "cov_EmpiricalCovariance") %>%
  left_join(., UCLA_CNP_brain_region_info)  %>%
  dplyr::select(-Index) %>%
  left_join(., ABIDE_ASD_brain_region_info) %>%
  filter(Comparison_Group != "ASD")

correlation_plots <- plot_data_per_group(SPI_nickname = "Pearson R",
                                         input_data = pearson_regional_data,
                                         group_colors = group_colors,
                                         bin_interval = 0.3)

wrap_plots(correlation_plots, 
           ncol=2, 
           byrow=T)  + 
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom',
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.text = element_text(size=14),
        legend.title = element_blank()) &
  guides(fill = guide_colorsteps(title.position="top", ticks=TRUE, barwidth=12,
                                 nrow=3))
ggsave(glue("{plot_path}/Region_wise_avg_t_stat_pearson_corrs.png"),
       width=5, height=8, units="in", dpi=300)

# Then visualize Gaussian DI from each region
gaussian_DI_from_regional_data <- pairwise_t_stats_by_region_from %>%
  filter(SPI == "di_gaussian") %>%
  left_join(., UCLA_CNP_brain_region_info)  %>%
  dplyr::select(-Index) %>%
  left_join(., ABIDE_ASD_brain_region_info) %>%
  filter(Comparison_Group != "ASD")

DI_gaussian_from_plots <- plot_data_per_group(SPI_nickname = "DI Gaussian",
                                         input_data = gaussian_DI_from_regional_data,
                                         group_colors = group_colors,
                                         bin_interval = 0.4)

wrap_plots(DI_gaussian_from_plots, 
           ncol=2, 
           byrow=T)  + 
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom',
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.text = element_text(size=14),
        legend.title = element_blank()) &
  guides(fill = guide_colorsteps(title.position="top", ticks=TRUE, barwidth=12,
                                 nrow=3))
ggsave(glue("{plot_path}/Region_wise_avg_t_stat_DI_gaussian_from_corrs.png"),
       width=5, height=8, units="in", dpi=300)

# Then visualize Gaussian DI to each region
gaussian_DI_to_regional_data <- pairwise_t_stats_by_region_to %>%
  filter(SPI == "di_gaussian") %>%
  left_join(., UCLA_CNP_brain_region_info)  %>%
  dplyr::select(-Index) %>%
  left_join(., ABIDE_ASD_brain_region_info) %>%
  filter(Comparison_Group != "ASD")

DI_gaussian_to_plots <- plot_data_per_group(SPI_nickname = "DI Gaussian",
                                              input_data = gaussian_DI_to_regional_data,
                                              group_colors = group_colors,
                                              bin_interval = 0.8)

wrap_plots(DI_gaussian_to_plots, 
           ncol=2, 
           byrow=T)  + 
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom',
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.text = element_text(size=14),
        legend.title = element_blank()) &
  guides(fill = guide_colorsteps(title.position="top", ticks=TRUE, barwidth=12,
                                 nrow=3))
ggsave(glue("{plot_path}/Region_wise_avg_t_stat_DI_gaussian_to_corrs.png"),
       width=5, height=8, units="in", dpi=300)

################################################################################
# Compare regional FC T-statistics with univariate region-wise accuracy
################################################################################

plot_regional_SPI_vs_balacc <- function(input_data, 
                                        label.y.pos,
                                        ylab) {
  # Plot SPI vs univariate balanced accuracy in a scatter plot
  p <- input_data %>%
    dplyr::select(Study, Comparison_Group, mean_T_magnitude, Brain_Region) %>%
    left_join(., univariate_balanced_accuracy %>%
                dplyr::select(Study, Comparison_Group, Analysis_Type, group_var, Balanced_Accuracy_Across_Repeats) %>%
                filter(Analysis_Type == "Univariate_Brain_Region") %>%
                mutate(Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                                    Comparison_Group == "Bipolar" ~ "BPD",
                                                    T ~ Comparison_Group)) %>%
                dplyr::rename("Brain_Region" = "group_var")) %>%
    mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD"))) %>%
    ggplot(data=., mapping=aes(x=100*Balanced_Accuracy_Across_Repeats, y=mean_T_magnitude)) +
    geom_point(aes(color = Comparison_Group)) +
    facet_wrap(Comparison_Group ~ ., nrow=2) +
    stat_smooth(method="lm", 
                color="black",
                se=F) +
    stat_cor(method="spearman", 
             size = 5,
             label.y = label.y.pos,
             cor.coef.name = "rho", 
             p.accuracy = 0.01) +
    scale_y_continuous(expand=c(0,0,0.1,0)) +
    scale_color_manual(values = group_colors) +
    ylab(ylab) +
    xlab("Mean Balanced Accuracy for Brain Region") +
    theme(legend.position = "none")
  
  return(p)
}

# Pearson correlation
plot_regional_SPI_vs_balacc(input_data=pearson_regional_data, 
                            label.y.pos = 1.75,
                            ylab = "Mean Pearson T for Brain Region")
ggsave(glue("{plot_path}/Mean_Pearson_T_vs_Univariate_Balanced_Accuracy.png"),
       width=6, height=5, units="in", dpi=300)

# DI-Gaussian, from
plot_regional_SPI_vs_balacc(input_data=gaussian_DI_from_regional_data, 
                            label.y.pos = 2.6,
                            ylab = "Mean DI-Gaussian from each Brain Region")
ggsave(glue("{plot_path}/Mean_DI_Gaussian_From_vs_Univariate_Balanced_Accuracy.png"),
       width=6, height=5, units="in", dpi=300)

# DI-Gaussian, to
plot_regional_SPI_vs_balacc(input_data=gaussian_DI_to_regional_data, 
                            label.y.pos = 3.2,
                            ylab = "Mean DI-Gaussian to each Brain Region")
ggsave(glue("{plot_path}/Mean_DI_Gaussian_To_vs_Univariate_Balanced_Accuracy.png"),
       width=6, height=5, units="in", dpi=300)

# Generate table for all undirected SPIs
pairwise_t_stats_by_region_from %>%
  filter(Study == "UCLA_CNP") %>%
  left_join(., UCLA_CNP_brain_region_info) %>%
  dplyr::select(-Index) %>%
  dplyr::select(Study, Comparison_Group, SPI, mean_T_magnitude, Brain_Region) %>%
  left_join(SPI_info) %>%
  filter(Directed == "No") %>%
  left_join(., univariate_balanced_accuracy %>%
              dplyr::select(Study, Comparison_Group, Analysis_Type, group_var, Balanced_Accuracy_Across_Repeats) %>%
              filter(Analysis_Type == "Univariate_Brain_Region") %>%
              mutate(Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                                  Comparison_Group == "Bipolar" ~ "BPD",
                                                  T ~ Comparison_Group)) %>%
              dplyr::rename("Brain_Region" = "group_var")) %>%
  group_by(Comparison_Group, SPI) %>%
  nest() %>%
  mutate(
    fit = map(data, ~ cor.test(.x$Balanced_Accuracy_Across_Repeats, 
                               .x$mean_T_magnitude, method="spearman")),
    tidied = map(fit, tidy)
  ) %>% 
  unnest(tidied) %>%
  dplyr::select(-data, -fit)
  

################################################################################
# Plot T-statistics to vs. from

# Spaghetti plot to vs from for gaussian DI per brain region by group, DIRECTED
pairwise_t_stats_by_region_from %>%
  plyr::rbind.fill(pairwise_t_stats_by_region_to) %>%
  left_join(., SPI_info) %>%
  filter(Directed == "Yes") %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD")),
         SPI = gsub("_", " ", SPI)) %>%
  ggplot(data=., mapping=aes(x=Direction, y=mean_T_magnitude, group=Brain_Region, color=Comparison_Group)) +
  geom_line(linewidth=0.4, alpha=0.5) +
  facet_grid(Nickname ~ Comparison_Group, scales="free", switch="y",
             labeller = labeller(Nickname = label_wrap_gen(15))) +
  ylab("Mean T-Statistic by Brain Region") +
  xlab("Functional Connectivity Direction")  +
  scale_color_manual(values = c("SCZ"="#573DC7", 
                                "BPD"="#D5492A", 
                                "ADHD"="#0F9EA9",
                                "ASD"="#C47B2F"))  +
  theme(legend.position = "none",
        strip.placement = "outside",
        strip.text.y.left = element_text(angle=0))
ggsave(glue("{plot_path}/Directed_SPI_T_Stats_To_vs_From.png"),
       width=9, height=7.5, units="in", dpi=300)

# Spaghetti plot to vs from for gaussian DI per brain region by group, UNDIRECTED
pairwise_t_stats_by_region_from %>%
  plyr::rbind.fill(pairwise_t_stats_by_region_to) %>%
  left_join(., SPI_info) %>%
  filter(Directed == "No") %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD")),
         SPI = gsub("_", " ", SPI)) %>%
  ggplot(data=., mapping=aes(x=Direction, y=mean_T_magnitude, group=Brain_Region, color=Comparison_Group)) +
  geom_line(linewidth=0.4, alpha=0.5) +
  facet_grid(Nickname ~ Comparison_Group, scales="free", switch="y",
             labeller = labeller(Nickname = label_wrap_gen(15))) +
  ylab("Mean T-Statistic by Brain Region") +
  xlab("Functional Connectivity Direction")  +
  scale_color_manual(values = c("SCZ"="#573DC7", 
                                "BPD"="#D5492A", 
                                "ADHD"="#0F9EA9",
                                "ASD"="#C47B2F"))  +
  theme(legend.position = "none",
        strip.placement = "outside",
        strip.text.y.left = element_text(angle=0))
ggsave(glue("{plot_path}/Undirected_SPI_T_Stats_To_vs_From.png"),
       width=9, height=7.5, units="in", dpi=300)

# # Check if 'directed' features are actually bidirectionally symmetric and vice versa
# UCLA_CNP_pyspi14 <- pyarrow_feather$read_feather("~/data/UCLA_CNP/processed_data/UCLA_CNP_AROMA_2P_GMR_pyspi14_filtered.feather")  %>%
#   left_join(., UCLA_CNP_metadata) %>%
#   filter(!is.na(Diagnosis)) %>%
#   mutate(Study = "UCLA_CNP")
# 
# # Take sub-10171 phase slope index, frequency domain
# UCLA_CNP_pyspi14 %>%
#   filter(Sample_ID == "sub-10171",
#          SPI == "psi_multitaper_mean_fs-1_fmin-0_fmax-0-5") %>%
#   filter(brain_region_from == 'ctx-lh-bankssts' | brain_region_to == 'ctx-lh-bankssts') %>%
#   write.table("~/Desktop/temp.csv", row.names=F, sep=",")

################################################################################
# All SPIs in one model
################################################################################

null_data_for_plot <- pairwise_null_distribution %>%
  filter(Analysis_Type == "Pairwise_SPI") %>%
  dplyr::rename("Balanced_Accuracy_Across_Folds" = Null_Balanced_Accuracy) %>%
  mutate(Balanced_Accuracy_Across_Folds = 100*Balanced_Accuracy_Across_Folds,
         Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                      Comparison_Group == "Bipolar" ~ "BPD",
                                      T ~ Comparison_Group)) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD"))) %>%
  dplyr::select(Comparison_Group, Balanced_Accuracy_Across_Folds)

pairwise_all_SPIs_balanced_accuracy_by_repeats %>%
  mutate(Balanced_Accuracy_Across_Folds = 100*Balanced_Accuracy_Across_Folds,
         Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                      Comparison_Group == "Bipolar" ~ "BPD",
                                      T ~ Comparison_Group)) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD"))) %>%
  ggplot(data=., mapping=aes(x=Comparison_Group,
                             y=Balanced_Accuracy_Across_Folds)) +
  geom_violin(aes(fill=Comparison_Group), trim=TRUE) +
  scale_fill_manual(values=c("#573DC7", "#D5492A", "#0F9EA9","#C47B2F")) +
  scale_color_manual(values=c("#573DC7", "#D5492A", "#0F9EA9","#C47B2F")) +
  geom_boxplot(width=0.1, notch=FALSE, notchwidth = 0.4, outlier.shape = NA,
               fill=NA, color="black",
               coef = 0) +
  geom_hline(yintercept = 50, linetype=2, alpha=0.5) +
  xlab("Clinical Group") +
  ylab("Balanced Accuracy\nper Repeat (%)") +
  theme(legend.position = "none",
        axis.title = element_text(size=17), 
        axis.text = element_text(size=15)) 
ggsave(glue("{plot_path}/Pairwise_all_SPIs_results.png"),
       width=4.5, height=3, units="in", dpi=300)

################################################################################
# Hypo vs hyper connectivity
################################################################################

plot_hyper_hypo_FC_group <- function(comparison_group,
                                     study,
                                     t_stat_df,
                                     t_mag_threshold) {
  # Edges are defined as cortical lobe --> specific ROI connection
  edges <- region_node_to_from %>%
    filter(Study == study) %>% 
    distinct() %>%
    dplyr::select(-Study)
  
  # ROIs don't include the origin --> cortical lobe connection
  rois <- edges %>% filter(!(to %in% c("Cingulate", "Frontal", "Insula",
                                       "Occipital", "Parietal", "Temporal", "Subcortex")))
  
  # Create a dataframe of vertices, one line per object in the ROI cortical lobe hierarchy
  vertices = data.frame(name = unique(c(as.character(edges$from), as.character(edges$to))))
  vertices$group <- edges$from[match(vertices$name, edges$to)]
  
  # Create an igraph object
  mygraph <- graph_from_data_frame(d=edges, vertices=vertices)
  
  # connect = dataframe of pairwise correlations between cortical ROIs
  hyper_data <- t_stats_pyspi14_whole_brain %>%
    filter(Comparison_Group == comparison_group) %>%
    group_by(Comparison_Group, Region_Pair) %>%
    summarise(mean_T = mean(estimate)) %>%
    separate(Region_Pair, into=c("from", "to"), sep="_") %>%
    filter(mean_T > t_mag_threshold)
  
  hyper_connect <- hyper_data %>%
    rename("value" = "mean_T") %>%
    arrange(from, to)
  
  # mygraph = igraph object linking each cortical ROI
  # convert to a circular dendrogram-shaped ggraph object
  p_hyper <- ggraph(mygraph, layout = 'dendrogram', circular = TRUE)
  
  from <- match(hyper_connect$from, vertices$name)
  to <- match(hyper_connect$to, vertices$name)
  
  p_hyper <- p_hyper + geom_conn_bundle(data = get_con(from = from, to = to, 
                                                       value=hyper_connect$value), 
                                        tension=0.7, width=1.5,
                                        aes(color=value))  +
    scale_edge_color_gradientn(colors=c(alpha("red", 0.1), "red")) +
    labs(edge_color="Multi-metric\nT Statistic") + 
    geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour=group),   
                                       size=2) +
    scale_color_manual(values=c("Cingulate" = "#F8756D",
                                "Frontal" = "#C39A00",
                                "Occipital" = "#53B400",
                                "Parietal" = "#01BF93",
                                "Insula" = "#00B6EB",
                                "Subcortex" = "#A58AFF",
                                "Temporal" = "#FB61D7"
                                ), guide="none") +
    labs(color="Cortex") +
    theme_void() + 
    theme(plot.title=element_text(size=14, face="bold", hjust=0.5),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(10,10,10,10))
  
  # hypo
  hypo_data <- t_stats_pyspi14_whole_brain %>%
    filter(Comparison_Group == comparison_group) %>%
    group_by(Comparison_Group, Region_Pair) %>%
    summarise(mean_T = mean(estimate)) %>%
    separate(Region_Pair, into=c("from", "to"), sep="_") %>%
    filter(mean_T < -1*t_mag_threshold)
  
  hypo_connect <- hypo_data %>%
    rename("value" = "mean_T") %>%
    arrange(from, to)
  
  p_hypo <- ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
    theme_void()
  
  from <- match(hypo_connect$from, vertices$name)
  to <- match(hypo_connect$to, vertices$name)
  
  p_hypo <- p_hypo + 
    geom_conn_bundle(data = get_con(from = from, to = to, 
                                    value=hypo_connect$value), 
                     tension=0.7, width=1.5, aes(color=value))  +
    scale_edge_color_gradientn(colors=c("blue", alpha("blue", 0.2)), 
                               guide = guide_colourbar(available_aes = "edge_colour",
                                                       reverse = TRUE)) +
    labs(edge_color="Multi-metric\nT Statistic") + 
    geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour=group),
                    size=2) +
    scale_color_manual(values=c("Cingulate" = "#F8756D",
                                "Frontal" = "#C39A00",
                                "Occipital" = "#53B400",
                                "Parietal" = "#01BF93",
                                "Insula" = "#00B6EB",
                                "Subcortex" = "#A58AFF",
                                "Temporal" = "#FB61D7"
    ), guide="none") +
    labs(color="Cortex") +
    theme_void() + 
    theme(plot.title=element_text(size=14, face="bold", hjust=0.5),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(10,10,10,10))
  
  # Combine
  final_plot <- p_hyper / p_hypo
  
  return(final_plot)
  
}


# Schizophrenia
plot_hyper_hypo_FC_group(comparison_group = "SCZ",
  study = "UCLA_CNP",
  t_stat_df = t_stats_pyspi14_whole_brain,
  t_mag_threshold = 0.6)
ggsave(glue("{plot_path}/SCZ_network_plots.png"),
       width=4.25, height=6, units="in", dpi=300)

# BPD
plot_hyper_hypo_FC_group(comparison_group = "BPD",
                         study = "UCLA_CNP",
                         t_stat_df = t_stats_pyspi14_whole_brain,
                         t_mag_threshold = 0.6)
ggsave(glue("{plot_path}/BPD_network_plots.png"),
       width=4.25, height=6, units="in", dpi=300)

# ADHD
plot_hyper_hypo_FC_group(comparison_group = "ADHD",
                         study = "UCLA_CNP",
                         t_stat_df = t_stats_pyspi14_whole_brain,
                         t_mag_threshold = 0.6)
ggsave(glue("{plot_path}/ADHD_network_plots.png"),
       width=4.25, height=6, units="in", dpi=300)

t_stat_to_plot <- t_stat_df %>% 
  group_by(Comparison_Group, Region_Pair) %>%
  summarise(mean_T = mean(estimate)) %>%
  separate(Region_Pair, into=c("from", "to"), sep="_") %>%
  filter(abs(mean_T) > 0.6) %>%
  left_join(., vertices, by=c("from"="name")) %>%
  dplyr::rename("from_cortex" = "group") %>%
  left_join(., vertices, by=c("to"="name")) %>%
  dplyr::rename("to_cortex" = "group") %>%
  dplyr::select(from_cortex, to_cortex, mean_T, Comparison_Group, from, to)