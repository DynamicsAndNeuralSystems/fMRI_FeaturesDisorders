################################################################################
# Load libraries
################################################################################

python_to_use <- "/Users/abry4213/opt/anaconda3/envs/pyspi/bin/python3"
reticulate::use_python(python_to_use)

library(reticulate)
library(tidyverse)
library(icesTAF)
library(cowplot)
library(theft)
library(glue)
library(ggridges)
library(scales)
library(patchwork)
theme_set(theme_cowplot())

# Import pyarrow.feather as pyarrow_feather
pyarrow_feather <- import("pyarrow.feather")

################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
plot_path <- glue("{github_dir}/plots/Manuscript_Draft/null_visualisations/")
TAF::mkdir(plot_path)

UCLA_CNP_data_path <- "~/data/UCLA_CNP/processed_data/"
ABIDE_ASD_data_path <- "~/data/ABIDE_ASD/processed_data/"

# UNIVARIATE 
# Load catch24 mixed sigmoid-transformed null balanced accuracy for each dataset
UCLA_CNP_Schizophrenia_univariate_null <- pyarrow_feather$read_feather(glue("{UCLA_CNP_data_path}/UCLA_CNP_Schizophrenia_Univariate_catch24_mixedsigmoid_scaler_SVM_null_balanced_accuracy_distributions.feather")) %>%
  mutate(Analysis_Type = case_when(str_detect(group_var, "Left|Right|ctx") ~ "Brain_Region",
                                   str_detect(group_var, "Combo") ~ "Combo",
                                   T ~ "TS_Feature"),
         Study_Group = "UCLA Schizophrenia") %>%
  filter(Analysis_Type == "Brain_Region")

UCLA_CNP_Bipolar_univariate_null <- pyarrow_feather$read_feather(glue("{UCLA_CNP_data_path}/UCLA_CNP_Bipolar_Univariate_catch24_mixedsigmoid_scaler_SVM_null_balanced_accuracy_distributions.feather")) %>%
  mutate(Analysis_Type = case_when(str_detect(group_var, "Left|Right|ctx") ~ "Brain_Region",
                                   str_detect(group_var, "Combo") ~ "Combo",
                                   T ~ "TS_Feature"),
         Study_Group = "UCLA Bipolar") %>%
  filter(Analysis_Type == "Brain_Region")

UCLA_CNP_ADHD_univariate_null <- pyarrow_feather$read_feather(glue("{UCLA_CNP_data_path}/UCLA_CNP_ADHD_Univariate_catch24_mixedsigmoid_scaler_SVM_null_balanced_accuracy_distributions.feather")) %>%
  mutate(Analysis_Type = case_when(str_detect(group_var, "Left|Right|ctx") ~ "Brain_Region",
                                   str_detect(group_var, "Combo") ~ "Combo",
                                   T ~ "TS_Feature"),
         Study_Group = "UCLA ADHD") %>%
  filter(Analysis_Type == "Brain_Region")

ABIDE_ASD_univariate_null <- pyarrow_feather$read_feather(glue("{ABIDE_ASD_data_path}/ABIDE_ASD_ASD_Univariate_catch24_mixedsigmoid_scaler_SVM_null_balanced_accuracy_distributions.feather")) %>%
  mutate(Analysis_Type = case_when(str_detect(group_var, "_") ~ "TS_Feature",
                                   str_detect(group_var, "Combo") ~ "Combo",
                                   T ~ "Brain_Region"),
         Study_Group = "ABIDE ASD") %>%
  filter(Analysis_Type == "Brain_Region")

# Combine into one null dataset
UCLA_CNP_all_univariate_null <- do.call(plyr::rbind.fill, list(UCLA_CNP_ADHD_univariate_null, 
                                            UCLA_CNP_univariate_Bipolar_null,
                                            UCLA_CNP_univariate_Schizophrenia_null))


# Plot in one ridgeline plot
UCLA_univariate_null_plot <- UCLA_CNP_all_univariate_null %>%
  mutate(group_var = str_replace_all(group_var, "-", " ")) %>%
  mutate(group_var = str_replace_all(group_var, "ctx rh ", "Right ")) %>%
  mutate(group_var = str_replace_all(group_var, "ctx lh ", "Left ")) %>%
  ggplot(data=., mapping=aes(x=Null_Balanced_Accuracy, 
                             y=group_var, 
                             height = after_stat(density),
                             fill = group_var)) +
  geom_density_ridges(stat = "binline", bins = 20, scale = 1.75, 
                      draw_baseline = FALSE, alpha=0.6) +
  # scale_x_continuous(limits = c(0.2, 0.8)) +
  xlab("Null Balanced Accuracy") +
  ylab("Brain Regions")  +
  scale_y_discrete(labels = wrap_format(35))+
  facet_wrap(. ~ Study_Group, scales="free_x", nrow=1)  +
  theme(legend.position = "none",
        axis.text.y = element_text(size=9))

ABIDE_univariate_null_plot <- ABIDE_ASD_univariate_null %>%
  ggplot(data=., mapping=aes(x=Null_Balanced_Accuracy, 
                             y=group_var, 
                             height = after_stat(density),
                             fill = group_var)) +
  geom_density_ridges(stat = "binline", bins = 20, scale = 1.75, 
                      draw_baseline = FALSE, alpha=0.6) +
  xlab("Null Balanced Accuracy") +
  # scale_x_continuous(limits = c(0.2, 0.8)) +
  facet_wrap(. ~ Study_Group, scales="free", nrow=1)  +
  theme(legend.position = "none", 
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=8)) +
  scale_y_discrete(labels = wrap_format(45))

wrap_plots(list(UCLA_univariate_null_plot, ABIDE_univariate_null_plot), widths = c(0.75, 0.25))

ggsave(glue("{plot_path}/univariate_null_distributions.png"), 
       width = 12, height = 10, units="in", dpi=300)


################################################################################
# Load data
univariate_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_mixedsigmoid_scaler_balanced_accuracy_all_folds.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)
univariate_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_mixedsigmoid_scaler_null_balanced_accuracy_distributions.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)
univariate_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_mixedsigmoid_scaler_empirical_p_values.feather"))%>%
  filter(Univariate_Feature_Set == univariate_feature_set)

pairwise_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_mixedsigmoid_scaler_balanced_accuracy_all_folds.feather"))
pairwise_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_mixedsigmoid_scaler_empirical_p_values.feather"))
pairwise_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_mixedsigmoid_scaler_null_balanced_accuracy_distributions.feather"))

combined_univariate_pairwise_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combined_univariate_pairwise_mixedsigmoid_scaler_balanced_accuracy_all_folds.feather"))
combined_univariate_pairwise_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combined_univariate_pairwise_mixedsigmoid_scaler_empirical_p_values.feather"))
combined_univariate_pairwise_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combined_univariate_pairwise_mixedsigmoid_scaler_null_balanced_accuracy_distributions.feather"))

plot_data_vs_null <- function(input_null_distribution, main_results, 
                              analysis_type, line_width, line_colors,
                              num_bins = 50) {
  null_data_for_plot <- input_null_distribution %>%
    dplyr::select(Study, Analysis_Type, Comparison_Group, Null_Balanced_Accuracy) %>%
    mutate(Type = "Null",
           Null_Balanced_Accuracy = 100*Null_Balanced_Accuracy) %>%
    dplyr::rename("Balanced_Accuracy_Across_Repeats" = "Null_Balanced_Accuracy") %>%
    plyr::rbind.fill(., main_results %>%
                       dplyr::select(Study, Comparison_Group, Analysis_Type, Balanced_Accuracy_Across_Repeats, p_value_Bonferroni) %>%
                       mutate(Type = "Main",
                              Balanced_Accuracy_Across_Repeats = 100*Balanced_Accuracy_Across_Repeats)) %>%
    filter(Analysis_Type == analysis_type) %>%
    mutate(Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                        Comparison_Group == "Bipolar" ~ "BPD",
                                        T ~ Comparison_Group),
           significance = p_value_Bonferroni < 0.05) %>%
    mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD")))
  
  p <- null_data_for_plot %>%
    ggplot(data=., mapping=aes(x=Balanced_Accuracy_Across_Repeats)) +
    geom_vline(data = subset(null_data_for_plot, Type=="Main" & !(significance)),
               aes(xintercept = Balanced_Accuracy_Across_Repeats),
               color="gray90",
               linewidth=line_width) +
    geom_histogram(data = subset(null_data_for_plot, Type=="Null"),
                   fill="gray70", bins=num_bins) +
    geom_vline(data = subset(null_data_for_plot, Type=="Main" & significance),
               aes(xintercept = Balanced_Accuracy_Across_Repeats,
                   color = Comparison_Group),
               linewidth=line_width) +
    scale_color_manual(values=line_colors) +
    facet_wrap(Comparison_Group ~ ., ncol=2, scales="free") +
    xlab("Balanced Accuracy Across Repeats (%)") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          strip.placement = "outside") +
    theme(legend.position="none")
  
  return(p)
}

# Univariate region-wise
plot_data_vs_null(input_null_distribution=univariate_null_distribution, 
                                                      main_results=univariate_p_values, 
                                                      analysis_type= "Univariate_Brain_Region", 
                                                      line_width=0.2, 
                                                      line_colors=c("#573DC7", "#D5492A", "#0F9EA9", "#C47B2F"))
ggsave(glue("{plot_path}/univariate_region_main_vs_null_balanced_acc.png"),
       width=6, height=3.5, units="in", dpi=300)

# Univariate feature-wise
plot_data_vs_null(input_null_distribution=univariate_null_distribution, 
                                                      main_results=univariate_p_values, 
                                                      analysis_type= "Univariate_TS_Feature", 
                                                      line_width=0.3, 
                                                      line_colors=c("#573DC7", "#D5492A", "#0F9EA9", "#C47B2F"))
ggsave(glue("{plot_path}/univariate_feature_main_vs_null_balanced_acc.png"),
       width=6, height=3.5, units="in", dpi=300)

# Univariate combo-wise
plot_data_vs_null(input_null_distribution=univariate_null_distribution, 
                  main_results=univariate_p_values, 
                  analysis_type= "Univariate_Combo", 
                  line_width=1, 
                  num_bins=30,
                  line_colors=c("#573DC7", "#D5492A", "#0F9EA9", "#C47B2F"))
ggsave(glue("{plot_path}/univariate_combo_main_vs_null_balanced_acc.png"),
       width=4.5, height=3.5, units="in", dpi=300)

# Pairwise region-wise
plot_data_vs_null(input_null_distribution=pairwise_null_distribution, 
                  main_results=pairwise_p_values, 
                  analysis_type= "Pairwise_SPI", 
                  line_width=0.5, 
                  line_colors=c("#573DC7", "#D5492A", "#0F9EA9", "#C47B2F"))
ggsave(glue("{plot_path}/pairwise_SPI_main_vs_null_balanced_acc.png"),
       width=6, height=3.5, units="in", dpi=300)

# Pairwise SPI + univariate combo
plot_data_vs_null(input_null_distribution=combined_univariate_pairwise_null_distribution, 
                  main_results=combined_univariate_pairwise_p_values, 
                  analysis_type= "SPI_Univariate_Combo", 
                  line_width=0.5, 
                  line_colors=c("#573DC7", "#D5492A", "#0F9EA9", "#C47B2F"))
ggsave(glue("{plot_path}/combined_univariate_pairwise_SPI_wise_main_vs_null_balanced_acc.png"),
       width=6, height=3.5, units="in", dpi=300)