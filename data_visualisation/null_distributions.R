################################################################################
# Load libraries
################################################################################

python_to_use <- "/Users/abry4213/anaconda3/envs/pyspi/bin/python3"
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

github_dir <- "~/Library/CloudStorage/OneDrive-TheUniversityofSydney(Students)/github/fMRI_FeaturesDisorders/"
plot_path <- glue("{github_dir}/plots/Manuscript_Draft/methods_supplement/null_visualisations/")
TAF::mkdir(plot_path)

data_path <- "~/data/TS_feature_manuscript/"
UCLA_CNP_data_path <- "~/data/UCLA_CNP/processed_data/"
ABIDE_ASD_data_path <- "~/data/ABIDE_ASD/processed_data/"
univariate_feature_set <- "catch25"

study_group_df <- data.frame(Study = c(rep("UCLA_CNP", 3), "ABIDE_ASD"),
                             Comparison_Group = c("Schizophrenia", "Bipolar", "ADHD", "ASD"),
                             Group_Nickname = c("SCZ", "BP", "ADHD", "ASD"))


# UNIVARIATE 
# Univariate regional nulls
univariate_nulls <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_null_balanced_accuracy_distributions.feather"))  
univariate_regional_nulls <- univariate_nulls %>% filter(Analysis_Type == "Univariate_Brain_Region") %>%
  mutate(Null_Balanced_Accuracy = 100*Null_Balanced_Accuracy) %>%
  group_by(Analysis_Type, Comparison_Group, Study, group_var) %>%
  mutate(mu = mean(Null_Balanced_Accuracy),
         sigma = sd(Null_Balanced_Accuracy))


# Plot in one ridgeline plot
UCLA_CNP_summary_curves <- univariate_regional_nulls %>%
  filter(Study=="UCLA_CNP") %>%
  left_join(., study_group_df) %>%
  group_by(group_var, Group_Nickname) %>%
  mutate(Group_Nickname = factor(Group_Nickname, levels=c("SCZ", "BP", "ADHD"))) %>%
  reframe(x = Null_Balanced_Accuracy,
          y = dnorm(Null_Balanced_Accuracy, 
                    mean = mean(Null_Balanced_Accuracy),
                    sd = sd(Null_Balanced_Accuracy)))

UCLA_CNP_univariate_null_plot <- univariate_regional_nulls %>%
  left_join(., study_group_df) %>%
  filter(Study=="UCLA_CNP") %>%
  mutate(Group_Nickname = factor(Group_Nickname, levels=c("SCZ", "BP", "ADHD"))) %>%
  ggplot(data=., mapping=aes(fill = group_var, x=Null_Balanced_Accuracy)) +
  geom_histogram(aes(y = after_stat(density)), alpha=0.7) +
  facet_grid(group_var~Group_Nickname, scales="free", switch="y") +
  xlab("Null Balanced Accuracy (%)") +
  ylab("Brain Regions") +
  geom_line(data=UCLA_CNP_summary_curves, aes(x=x, y=y), color="black") +
  geom_vline(xintercept = 50, linetype=2) +
  theme(legend.position="none",
        strip.text.y.left = element_text(angle=0, size=9, lineheight = 0.6, hjust=1),
        strip.background = element_blank(),
        strip.placement="outside",
        panel.spacing = unit(-0.25, "lines"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) 

# Calculate the normal curve density for each group
ABIDE_summary_curves <- univariate_regional_nulls %>%
  filter(Study=="ABIDE_ASD") %>%
  left_join(., study_group_df) %>%
  group_by(group_var) %>%
  reframe(x = Null_Balanced_Accuracy,
          y = dnorm(Null_Balanced_Accuracy, 
                    mean = mean(Null_Balanced_Accuracy),
                    sd = sd(Null_Balanced_Accuracy)))

ABIDE_univariate_null_plot <- univariate_regional_nulls %>%
  filter(Study=="ABIDE_ASD") %>%
  left_join(., study_group_df) %>%
  ggplot(data=., mapping=aes(fill = group_var, x=Null_Balanced_Accuracy)) +
  geom_histogram(aes(y = after_stat(density)), alpha=0.7) +
  facet_grid(group_var~Group_Nickname, scales="free", switch="y",
             labeller = labeller(group_var = label_wrap_gen(45))) +
  xlab("Null Balanced Accuracy (%)") +
  ylab("Brain Regions") +
  geom_line(data=ABIDE_summary_curves, aes(x=x, y=y, fill=group_var), color="black") +
  geom_vline(xintercept = 50, linetype=2) +
  theme(legend.position="none",
        strip.text.y.left = element_text(angle=0, size=9, lineheight = 0.6, hjust=1),
        strip.background = element_blank(),
        strip.placement="outside",
        panel.spacing = unit(-0.25, "lines"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) 

wrap_plots(list(UCLA_CNP_univariate_null_plot, ABIDE_univariate_null_plot), widths = c(0.75, 0.25))

ggsave(glue("{plot_path}/univariate_null_distributions.svg"), 
       width = 12, height = 10, units="in", dpi=300)


################################################################################
# Load data
univariate_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_balanced_accuracy_all_folds.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)
univariate_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_null_balanced_accuracy_distributions.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)
univariate_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_empirical_p_values.feather"))%>%
  filter(Univariate_Feature_Set == univariate_feature_set)

pairwise_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_balanced_accuracy_all_folds.feather"))
pairwise_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_empirical_p_values.feather"))
pairwise_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_null_balanced_accuracy_distributions.feather"))

combined_univariate_pairwise_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combined_univariate_pairwise_balanced_accuracy_all_folds.feather"))
combined_univariate_pairwise_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combined_univariate_pairwise_empirical_p_values.feather"))
combined_univariate_pairwise_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combined_univariate_pairwise_null_balanced_accuracy_distributions.feather"))

plot_data_vs_null <- function(input_null_distribution, main_results, 
                              analysis_type, line_width, line_colors,
                              num_bins = 50) {
  null_data_for_plot <- input_null_distribution %>%
    dplyr::select(Study, Analysis_Type, Comparison_Group, Null_Balanced_Accuracy) %>%
    mutate(Type = "Null",
           Null_Balanced_Accuracy = 100*Null_Balanced_Accuracy) %>%
    dplyr::rename("Balanced_Accuracy_Across_Folds" = "Null_Balanced_Accuracy") %>%
    plyr::rbind.fill(., main_results %>%
                       dplyr::select(Study, Comparison_Group, Analysis_Type, Balanced_Accuracy_Across_Folds, p_value_Bonferroni) %>%
                       mutate(Type = "Main",
                              Balanced_Accuracy_Across_Folds = 100*Balanced_Accuracy_Across_Folds)) %>%
    filter(Analysis_Type == analysis_type) %>%
    mutate(Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                        Comparison_Group == "Bipolar" ~ "BP",
                                        T ~ Comparison_Group),
           significance = p_value_Bonferroni < 0.05) %>%
    mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BP", "ADHD", "ASD")))
  
  p <- null_data_for_plot %>%
    ggplot(data=., mapping=aes(x=Balanced_Accuracy_Across_Folds)) +
    geom_vline(data = subset(null_data_for_plot, Type=="Main" & !(significance)),
               aes(xintercept = Balanced_Accuracy_Across_Folds),
               color="gray90",
               linewidth=line_width) +
    geom_histogram(data = subset(null_data_for_plot, Type=="Null"),
                   fill="gray70", bins=num_bins) +
    geom_vline(data = subset(null_data_for_plot, Type=="Main" & significance),
               aes(xintercept = Balanced_Accuracy_Across_Folds,
                   color = Comparison_Group),
               linewidth=line_width) +
    scale_color_manual(values=line_colors) +
    facet_wrap(Comparison_Group ~ ., ncol=2, scales="free") +
    xlab("Balanced Accuracy (%)") +
    theme(axis.text.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face="bold"),
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
ggsave(glue("{plot_path}/univariate_region_main_vs_null_balanced_acc.svg"),
       width=6, height=3.5, units="in", dpi=300)

# Univariate feature-wise
plot_data_vs_null(input_null_distribution=univariate_null_distribution, 
                                                      main_results=univariate_p_values, 
                                                      analysis_type= "Univariate_TS_Feature", 
                                                      line_width=0.3, 
                                                      line_colors=c("#573DC7", "#D5492A", "#0F9EA9", "#C47B2F"))
ggsave(glue("{plot_path}/univariate_feature_main_vs_null_balanced_acc.svg"),
       width=6, height=3.5, units="in", dpi=300)

# Univariate combo-wise
plot_data_vs_null(input_null_distribution=univariate_null_distribution, 
                  main_results=univariate_p_values, 
                  analysis_type= "Univariate_Combo", 
                  line_width=1, 
                  num_bins=30,
                  line_colors=c("#573DC7", "#D5492A", "#0F9EA9", "#C47B2F"))
ggsave(glue("{plot_path}/univariate_combo_main_vs_null_balanced_acc.svg"),
       width=4.5, height=3.5, units="in", dpi=300)

# Pairwise region-wise
plot_data_vs_null(input_null_distribution=pairwise_null_distribution, 
                  main_results=pairwise_p_values, 
                  analysis_type= "Pairwise_SPI", 
                  line_width=0.5, 
                  line_colors=c("#573DC7", "#D5492A", "#0F9EA9", "#C47B2F"))
ggsave(glue("{plot_path}/pairwise_SPI_main_vs_null_balanced_acc.svg"),
       width=6, height=3.5, units="in", dpi=300)

# Pairwise SPI + univariate combo
plot_data_vs_null(input_null_distribution=combined_univariate_pairwise_null_distribution, 
                  main_results=combined_univariate_pairwise_p_values, 
                  analysis_type= "SPI_Univariate_Combo", 
                  line_width=0.5, 
                  line_colors=c("#573DC7", "#D5492A", "#0F9EA9", "#C47B2F"))
ggsave(glue("{plot_path}/combined_univariate_pairwise_SPI_wise_main_vs_null_balanced_acc.svg"),
       width=6, height=3.5, units="in", dpi=300)
