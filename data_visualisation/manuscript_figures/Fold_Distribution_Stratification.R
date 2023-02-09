################################################################################
# Load libraries
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
UCLA_CNP_data_path <- "~/data/UCLA_CNP/"
ABIDE_ASD_data_path <- "~/data/ABIDE_ASD/"
plot_path <- "~/github/fMRI_FeaturesDisorders/plots/Manuscript_Draft/FigureS7/"
TAF::mkdir(plot_path)

python_to_use <- "~/.conda/envs/pyspi/bin/python3"

reticulate::use_python(python_to_use)

library(reticulate)

# Import pyarrow.feather as pyarrow_feather
pyarrow_feather <- import("pyarrow.feather")


################################################################################
# Define study/data paths
################################################################################

library(tidyverse)
library(knitr)
library(kableExtra)
library(glue)
library(patchwork)

# Load metadata
UCLA_CNP_sample_metdata <- pyarrow_feather$read_feather(glue("{UCLA_CNP_data_path}/study_metadata/UCLA_CNP_sample_metadata.feather"))
ABIDE_ASD_sample_metdata <- pyarrow_feather$read_feather(glue("{ABIDE_ASD_data_path}/study_metadata/ABIDE_ASD_sample_metadata.feather"))

# Load fold assignments
UCLA_CNP_schizophrenia_catch22_fold_assignments <- pyarrow_feather$read_feather(glue("{UCLA_CNP_data_path}/processed_data/UCLA_CNP_Schizophrenia_Univariate_catch22_robustsigmoid_scaler_SVM_fold_assignments.feather")) %>%
  left_join(., UCLA_CNP_sample_metdata)

UCLA_CNP_ADHD_catch22_fold_assignments <- pyarrow_feather$read_feather(glue("{UCLA_CNP_data_path}/processed_data/UCLA_CNP_ADHD_Univariate_catch22_robustsigmoid_scaler_SVM_fold_assignments.feather")) %>%
  left_join(., UCLA_CNP_sample_metdata)

UCLA_CNP_bipolar_catch22_fold_assignments <- pyarrow_feather$read_feather(glue("{UCLA_CNP_data_path}/processed_data/UCLA_CNP_Bipolar_Univariate_catch22_robustsigmoid_scaler_SVM_fold_assignments.feather")) %>%
  left_join(., UCLA_CNP_sample_metdata)

ABIDE_ASD_catch22_fold_assignments <- pyarrow_feather$read_feather(glue("{ABIDE_ASD_data_path}/processed_data/ABIDE_ASD_ASD_Univariate_catch22_robustsigmoid_scaler_SVM_fold_assignments.feather")) %>%
  left_join(., ABIDE_ASD_sample_metdata)

################################################################################
# Heatmap visualisations
################################################################################

plot_heatmap_for_dataset <- function(fold_assignments_df, group_var_to_use, plot_title) {
  p <- fold_assignments_df %>%
          filter(group_var ==group_var_to_use) %>%
          mutate(Repeat = factor(Repeat),
                Fold = factor(Fold)) %>%
          ggplot(data=., mapping=aes(x = Sample_ID, y = Repeat, fill = Fold )) +
          facet_grid(. ~ Diagnosis, scales="free", space="free") +
          geom_tile() +
          xlab("Samples") +
          ylab("CV-SVM Repeat") +
          labs(fill = "k-fold") +
          ggtitle(plot_title) +
          scale_fill_viridis_d() +
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                plot.title = element_text(hjust=0.5),
                legend.position = "bottom") +
          guides(fill = guide_legend(title.position = "top", 
                                    nrow = 1,
                                      title.hjust = 0.5,
                                      label.position = "bottom")) 

}

# UCLA CNP schizophrenia vs. control ctx-lh-caudalanteriorcingulate
UCLA_scz_heatmap <- plot_heatmap_for_dataset(fold_assignments_df = UCLA_CNP_schizophrenia_catch22_fold_assignments, 
                                             group_var_to_use = "ctx-lh-caudalanteriorcingulate", 
                                             plot_title = "UCLA CNP -- Schizophrenia")


# UCLA CNP ADHD vs. control ctx-lh-caudalanteriorcingulate
UCLA_ADHD_heatmap <- plot_heatmap_for_dataset(fold_assignments_df = UCLA_CNP_ADHD_catch22_fold_assignments, 
                                             group_var_to_use = "ctx-lh-caudalanteriorcingulate", 
                                             plot_title = "UCLA CNP -- ADHD")

# UCLA CNP ADHD vs. control ctx-lh-caudalanteriorcingulate
UCLA_bipolar_heatmap <- plot_heatmap_for_dataset(fold_assignments_df = UCLA_CNP_bipolar_catch22_fold_assignments, 
                                             group_var_to_use = "ctx-lh-caudalanteriorcingulate", 
                                             plot_title = "UCLA CNP -- Bipolar")


# Use ABIDE ASD vs. control Superior Frontal Gyrus as example
ABIDE_ASD_heatmap <- plot_heatmap_for_dataset(fold_assignments_df = ABIDE_ASD_catch22_fold_assignments, 
                                             group_var_to_use = "Superior Frontal Gyrus", 
                                             plot_title = "ABIDE -- ASD")

wrap_plots(list(UCLA_scz_heatmap,
                UCLA_ADHD_heatmap,
                UCLA_bipolar_heatmap,
                ABIDE_ASD_heatmap),
           ncol = 2) + 
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')
ggsave(glue("{plot_path}/All_catch22_robustsigmoid_scaler_fold_distributions.png"),
       width = 12, height = 5, units="in", dpi=300)

################################################################################
# Calculate distributions across folds
################################################################################

# UCLA Schizophrenia vs. control
UCLA_CNP_schizophrenia_catch22_fold_assignments %>%
  mutate(num_samples = length(unique(Sample_ID))) %>%
  group_by(group_var, Repeat) %>%
  mutate(Control_Prop = sum(Diagnosis == "Control")/num_samples) %>%
  group_by(group_var, Fold, Repeat, Control_Prop) %>%
  summarise(num_subjects_in_test_fold = n(),
         num_controls_in_test_fold = sum(Diagnosis == "Control"),
         control_prop_in_fold = num_controls_in_test_fold/num_subjects_in_test_fold) %>%
  group_by(Control_Prop) %>%
  summarise(mean_control_prop_across_folds = 100*mean(control_prop_in_fold),
            min_control_prop_across_folds = 100*min(control_prop_in_fold),
            max_control_prop_across_folds = 100*max(control_prop_in_fold),
            SD_control_prop_across_folds = 100*sd(control_prop_in_fold))


# UCLA ADHD vs. control
UCLA_CNP_ADHD_catch22_fold_assignments %>%
  mutate(num_samples = length(unique(Sample_ID))) %>%
  group_by(group_var, Repeat) %>%
  mutate(Control_Prop = sum(Diagnosis == "Control")/num_samples) %>%
  group_by(group_var, Fold, Repeat, Control_Prop) %>%
  summarise(num_subjects_in_test_fold = n(),
            num_controls_in_test_fold = sum(Diagnosis == "Control"),
            control_prop_in_fold = num_controls_in_test_fold/num_subjects_in_test_fold) %>%
  group_by(Control_Prop) %>%
  summarise(mean_control_prop_across_folds = 100*mean(control_prop_in_fold),
            SD_control_prop_across_folds = 100*sd(control_prop_in_fold),
            min_control_prop_across_folds = 100*min(control_prop_in_fold),
            max_control_prop_across_folds = 100*max(control_prop_in_fold)
            )


# UCLA bipolar vs. control
UCLA_CNP_bipolar_catch22_fold_assignments %>%
  mutate(num_samples = length(unique(Sample_ID))) %>%
  group_by(group_var, Repeat) %>%
  mutate(Control_Prop = sum(Diagnosis == "Control")/num_samples) %>%
  group_by(group_var, Fold, Repeat, Control_Prop) %>%
  summarise(num_subjects_in_test_fold = n(),
            num_controls_in_test_fold = sum(Diagnosis == "Control"),
            control_prop_in_fold = num_controls_in_test_fold/num_subjects_in_test_fold) %>%
  group_by(Control_Prop) %>%
  summarise(mean_control_prop_across_folds = 100*mean(control_prop_in_fold),
            SD_control_prop_across_folds = 100*sd(control_prop_in_fold),
            min_control_prop_across_folds = 100*min(control_prop_in_fold),
            max_control_prop_across_folds = 100*max(control_prop_in_fold)
  )


# ABIDE ASD vs. control
ABIDE_ASD_catch22_fold_assignments %>%
  mutate(num_samples = length(unique(Sample_ID))) %>%
  group_by(group_var, Repeat) %>%
  mutate(Control_Prop = sum(Diagnosis == "Control")/num_samples) %>%
  group_by(group_var, Fold, Repeat, Control_Prop) %>%
  summarise(num_subjects_in_test_fold = n(),
            num_controls_in_test_fold = sum(Diagnosis == "Control"),
            control_prop_in_fold = num_controls_in_test_fold/num_subjects_in_test_fold) %>%
  group_by(Control_Prop) %>%
  summarise(mean_control_prop_across_folds = 100*mean(control_prop_in_fold),
            min_control_prop_across_folds = 100*min(control_prop_in_fold),
            max_control_prop_across_folds = 100*max(control_prop_in_fold),
            SD_control_prop_across_folds = 100*sd(control_prop_in_fold))
  