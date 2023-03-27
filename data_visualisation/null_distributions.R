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