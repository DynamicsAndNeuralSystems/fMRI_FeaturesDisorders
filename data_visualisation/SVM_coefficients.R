################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/SVM_coefficients/")
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
ggsave(glue("{plot_path}/Combo_Region_Wise_SVM_Coefs.svg"),
       width=4, height=7, units="in", dpi=300)
wrap_plots(legend_list, 
           nrow=4, 
           byrow=T)
ggsave(glue("{plot_path}/Combo_Region_Wise_SVM_Coefs_legend.svg"),
       width=3, height=4, units="in", dpi=300)


################################################################################
# Plot average magnitude of all brain regions' coefficients per catch24 feature
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
ggsave(glue("{plot_path}/Combo_feature_wise_SVM_coef_colorbar.svg"),
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
ggsave(glue("{plot_path}/Combo_feature_wise_SVM_coef.svg"),
       width=6, height=6, units="in", dpi=300)
