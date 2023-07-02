################################################################################
# Load libraries
################################################################################

# python_to_use <- "~/.conda/envs/pyspi/bin/python3"
python_to_use <- "/Users/abry4213/anaconda3/envs/pyspi/bin/python3"
pairwise_feature_set <- "pyspi14"
univariate_feature_set <- "catch24"
study_group_df <- data.frame(Study = c(rep("UCLA_CNP", 3)),
                             Noise_Proc = c(rep("AROMA+2P+GMR",3)),
                             Comparison_Group = c("Schizophrenia", "Bipolar", "ADHD"),
                             Group_Nickname = c("SCZ", "BPD", "ADHD"))
reticulate::use_python(python_to_use)

library(reticulate)
library(tidyverse)
library(icesTAF)
library(cowplot)
library(theft)
library(glue)
library(ggseg)
theme_set(theme_cowplot())
pyarrow_feather <- import("pyarrow.feather")

################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
source(glue("{github_dir}/data_visualisation/Manuscript_Draft_Visualisations_Helper.R"))
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/univariate_results/")
icesTAF::mkdir(plot_path)

# Load in univariate time-series feature info
TS_feature_info <- read.csv(glue("{github_dir}/data_visualisation/catch24_info.csv"))
UCLA_CNP_brain_region_info <- read.csv("~/data/UCLA_CNP/study_metadata/UCLA_CNP_Brain_Region_info.csv")

# Load catch24 data for UCLA CNP
UCLA_CNP_catch24 <- pyarrow_feather$read_feather("~/data/UCLA_CNP/processed_data/UCLA_CNP_AROMA_2P_GMR_catch24_filtered.feather")

# Load study metadata
UCLA_CNP_metadata <- pyarrow_feather$read_feather("~/data/UCLA_CNP/study_metadata/UCLA_CNP_sample_metadata.feather") 


################################################################################
# Think about hemisphere differences in performance
run_correctR_group <- function(metadata, study, results_df) {
  results_across_groups <- list()
  
  for (comparison_group in unique(results_df$Comparison_Group)) {
    group_data = subset(results_df, Comparison_Group == comparison_group & Study == study)
    
    if (nrow(group_data) > 0) {
      num_subjects <- metadata %>%
        filter(Diagnosis %in% c("Control", comparison_group)) %>%
        distinct(Sample_ID) %>%
        nrow()
      training_size <- ceiling(0.9*num_subjects)
      test_size <- floor(0.1*num_subjects)
      
      data_for_correctR <- group_data %>% 
        mutate(Hemisphere = case_when(str_detect(group_var, "Left|lh-") ~ "Left",
                                      str_detect(group_var, "Right|rh-") ~ "Right")) %>%
        mutate(Brain_Region = gsub("Left-|ctx-lh-|Right-|ctx-rh-", "", group_var)) %>%
        left_join(., univariate_p_values) %>%
        group_by(Brain_Region) %>%
        filter(any(p_value_Bonferroni < 0.05)) %>%
        pivot_wider(id_cols = c(Repeat_Number, Brain_Region), 
                    names_from = Hemisphere,
                    values_from = Balanced_Accuracy_Across_Folds) %>%
        dplyr::rename("x" = "Right", "y" = "Left") %>%
        group_by(Brain_Region) %>%
        group_split()
      
      res <- data_for_correctR %>%
        purrr::map_df(~ as.data.frame(resampled_ttest(x=.x$x, 
                                                      y=.x$y, 
                                                      n=10, 
                                                      n1=training_size, n2=test_size)) %>%
                        mutate(Brain_Region = unique(.x$Brain_Region))) %>%
        ungroup() %>%
        dplyr::rename("p_value_corr"="p.value") %>%
        mutate(p_value_corr_Bonferroni = p.adjust(p_value_corr, method="bonferroni"),
               Comparison_Group = comparison_group)
      
      results_across_groups <- list.append(results_across_groups, res)
    }
  }
  results_across_groups_df <- do.call(plyr::rbind.fill, results_across_groups)
  return(results_across_groups_df)
}

hemisphere_balacc_t_res <- run_correctR_group(metadata = UCLA_CNP_metadata,
                                              study = "UCLA_CNP",
                                              results_df = subset(univariate_balanced_accuracy_AUC_by_repeats, 
                                                                  Analysis_Type == "Univariate_Brain_Region"))

# Plot regions across conditions
hemisphere_balacc_t_res %>%
  filter(p_value_corr_Bonferroni<0.05) %>%
  left_join(., study_group_df) %>%
  mutate(Brain_Region = fct_reorder(Brain_Region, statistic, 
                                    .fun="sum", .desc=TRUE),
         Group_Nickname = factor(Group_Nickname, levels=c("SCZ", "BPD", "ADHD"))) %>%
  ggplot(data=., mapping=aes(x = Group_Nickname, 
                             y = Brain_Region, 
                             fill = statistic)) +
  geom_tile() +
  ylab("Brain region") +
  xlab("Comparison group") +
  labs(fill = "T-statistic") +
  scale_fill_continuous_divergingx(palette = "PiYG", rev=T, mid=0, limits=c(-12, 12)) +
  geom_text(aes(label = round(statistic, 1))) +
  theme(legend.position = "bottom") +
  scale_y_discrete(limits=rev) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust=0.5, barwidth=10))
ggsave(glue("{plot_path}/Hemisphere_Asymmetry_T_Stat_Balacc.svg"),
       width=4.5, height=7, units="in", dpi=300)

# Plot T statistic in the brain
# Helper function to plot the beta coefficients for a given feature in the brain
plot_hemi_asym_in_brain <- function(study_group_df, 
                                    groups,
                                    hemi_data, 
                                    feature_name, min_fill,
                                    max_fill, fill_colors) {
  
  ggseg_plot_list <- list()
  
  for (comparison_group in groups) {
    
    group_hemi_data <- hemi_data %>%
      filter(Comparison_Group == comparison_group) %>%
      dplyr::rename("label"="Brain_Region") %>%
      distinct()
    
    # Plot T stat data in cortex
    dataset_ggseg <- plot_data_with_ggseg_gradient(dataset_ID = dataset_ID,
                                                   hemisphere = "right",
                                                   atlas_name = "dk",
                                                   atlas_data = dk %>% as_tibble(),
                                                   data_to_plot = group_hemi_data %>%
                                                     mutate(label = paste0("rh_", label)),
                                                   fill_variable = "statistic",
                                                   line_color = "gray30",
                                                   na_color = "white",
                                                   fill_colors = fill_colors,
                                                   min_fill = min_fill,
                                                   max_fill = max_fill)  +
      labs(fill="T")
    
    ggseg_plot_list <- list.append(ggseg_plot_list, dataset_ggseg)
    
    # Add subcortical data
    dataset_ggseg_subctx <- plot_data_with_ggseg_gradient(dataset_ID = dataset_ID,
                                                          hemisphere = "right",
                                                          atlas_name = "aseg",
                                                          atlas_data = aseg %>% as_tibble(),
                                                          data_to_plot = group_hemi_data %>%
                                                            mutate(label = paste0("Right-", label)),
                                                          fill_variable = "statistic",
                                                          fill_colors = fill_colors,
                                                          min_fill = min_fill,
                                                          max_fill = max_fill,
                                                          line_color = "gray30",
                                                          na_color = "white")  +
      labs(fill="T")
    
    # Append to list
    ggseg_plot_list <- list.append(ggseg_plot_list, dataset_ggseg_subctx)
  }
  
  return(ggseg_plot_list)
}

# Plot lm beta statistics for Hemisphere asymmetry T-statistic
min_fill <- floor(min(hemisphere_balacc_t_res$statistic))
max_fill <- ceiling(max(hemisphere_balacc_t_res$statistic))
fill_colors <- rev(RColorBrewer::brewer.pal(8, "PiYG"))

hemi_asymmetry_T_in_brain <- plot_hemi_asym_in_brain(study_group_df = study_group_df,
                                                     groups = c("Schizophrenia", "Bipolar", "ADHD"),
                                                     hemi_data = hemisphere_balacc_t_res,
                                                     feature_name = "Hemisphere T Stat",
                                                     min_fill = min_fill,
                                                     max_fill = max_fill,
                                                     fill_colors = fill_colors)

wrap_plots(hemi_asymmetry_T_in_brain, 
           ncol=2, 
           byrow=T) + 
  plot_layout(guides = "collect")  & 
  theme(legend.position = 'bottom')
ggsave(glue("{plot_path}/Hemisphere_Asymmetry_T_Stat_in_Brain.svg"),
       width=5, height=7, units="in", dpi=300)