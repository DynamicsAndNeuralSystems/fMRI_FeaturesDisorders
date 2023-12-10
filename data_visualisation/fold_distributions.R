################################################################################
# Load libraries
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
UCLA_CNP_data_path <- "~/data/UCLA_CNP/"
ABIDE_ASD_data_path <- "~/data/ABIDE_ASD/"
data_path <- "~/data/TS_feature_manuscript"
plot_path <- "~/github/fMRI_FeaturesDisorders/plots/Manuscript_Draft/fold_assignment/"
TAF::mkdir(plot_path)

# python_to_use <- "~/.conda/envs/pyspi/bin/python3"
python_to_use <- "/Users/abry4213/anaconda3/envs/pyspi/bin/python3"

reticulate::use_python(python_to_use)

library(reticulate)

# Import pyarrow.feather as pyarrow_feather
pyarrow_feather <- import("pyarrow.feather")


################################################################################
# Define study/data paths
################################################################################

library(tidyverse)
library(glue)
library(patchwork)
library(cowplot)
theme_set(theme_cowplot())

# Load metadata
UCLA_CNP_sample_metadata <- pyarrow_feather$read_feather(glue("{UCLA_CNP_data_path}/study_metadata/UCLA_CNP_sample_metadata.feather"))
ABIDE_ASD_sample_metadata <- pyarrow_feather$read_feather(glue("{ABIDE_ASD_data_path}/study_metadata/ABIDE_ASD_sample_metadata.feather"))

# Load fold assignments
univariate_fold_assignments <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_fold_assignments.feather")) %>%
  filter(Univariate_Feature_Set == "catch25")

################################################################################
# Heatmap visualisations
################################################################################


# UCLA CNP schizophrenia vs. control ctx-lh-caudalanteriorcingulate
UCLA_scz_heatmap <- univariate_fold_assignments %>%
  filter(Study == "UCLA_CNP", Comparison_Group  == "Schizophrenia") %>%
  left_join(., UCLA_CNP_sample_metadata) %>%
  plot_heatmap_for_dataset(fold_assignments_df = .,
                           group_var_to_use = "ctx-lh-caudalanteriorcingulate",
                           plot_title = "UCLA CNP -- Schizophrenia")


# UCLA CNP ADHD vs. control ctx-lh-caudalanteriorcingulate
UCLA_ADHD_heatmap <- univariate_fold_assignments %>%
  filter(Study == "UCLA_CNP", Comparison_Group  == "ADHD") %>%
  left_join(., UCLA_CNP_sample_metadata) %>%
  plot_heatmap_for_dataset(fold_assignments_df = .,
                           group_var_to_use = "ctx-lh-caudalanteriorcingulate",
                           plot_title = "UCLA CNP -- ADHD")


# UCLA CNP ADHD vs. control ctx-lh-caudalanteriorcingulate
UCLA_bipolar_heatmap <- univariate_fold_assignments %>%
  filter(Study == "UCLA_CNP", Comparison_Group  == "Bipolar") %>%
  left_join(., UCLA_CNP_sample_metadata) %>%
  plot_heatmap_for_dataset(fold_assignments_df = .,
                           group_var_to_use = "ctx-lh-caudalanteriorcingulate",
                           plot_title = "UCLA CNP -- Bipolar")


# Use ABIDE ASD vs. control Superior Frontal Gyrus as example
ABIDE_ASD_heatmap <- univariate_fold_assignments %>%
  filter(Study == "ABIDE_ASD", Comparison_Group  == "ASD") %>%
  left_join(., ABIDE_ASD_sample_metadata) %>%
  plot_heatmap_for_dataset(fold_assignments_df = .,
                           group_var_to_use = "Superior Frontal Gyrus",
                           plot_title = "ABIDE -- ASD")

wrap_plots(list(UCLA_scz_heatmap,
                UCLA_bipolar_heatmap,
                UCLA_ADHD_heatmap,
                ABIDE_ASD_heatmap),
           ncol = 1) + 
  plot_layout(guides = "collect") & theme(legend.position = 'bottom',
                                          plot.title = element_text(size=12))
ggsave(glue("{plot_path}/All_catch25_robustsigmoid_scaler_fold_distributions.svg"),
       width = 8, height = 9, units="in", dpi=300)

################################################################################
# Heatmap for all 82 brain regions in the UCLA schizophrenia vs control dataset 
# to show that samples are allocated in the same way across regions within repeat 1
################################################################################

univariate_fold_assignments %>%
  filter(Study == "UCLA_CNP", Comparison_Group == "Schizophrenia",
         Analysis_Type == "Brain_Region", Repeat == 1) %>%
  mutate(Fold = factor(Fold)) %>%
  ggplot(data = ., mapping = aes(x=fct_reorder(Sample_ID, as.numeric(Fold)), y = group_var, fill = Fold)) +
  geom_raster() +
  scale_fill_viridis_d() +
  ylab("Brain Region") +
  labs(fill = "Fold # in Repeat 1") +
  ggtitle("Sample Allocations to Folds for\nUCLA CNP Schizophrenia Cohort") +
  xlab("Samples") +
  theme(legend.position="bottom",
        plot.title = element_text(hjust=0.5, size=12),
        axis.text.y = element_text(size=8),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(title.position = "top", 
                             nrow = 1,
                             title.hjust = 0.5,
                             label.position = "bottom")) 
ggsave(glue("{plot_path}/UCLA_CNP_Schizophrenia_Brain_Regions_Repeat1_Allocations.svg"),
       width = 5, height = 9, units="in", dpi=300)


################################################################################
# Calculate distributions across folds
################################################################################

univariate_fold_assignments %>%
  left_join(., plyr::rbind.fill(UCLA_CNP_sample_metadata, ABIDE_ASD_sample_metadata)) %>%
  group_by(Comparison_Group, Study) %>%
  mutate(num_samples = length(unique(Sample_ID))) %>%
  group_by(Study, Comparison_Group, group_var, Repeat) %>%
  mutate(Control_Prop = sum(Diagnosis == "Control")/num_samples) %>%
  group_by(Study, Comparison_Group, group_var, Fold, Repeat, Control_Prop) %>%
  summarise(num_subjects_in_test_fold = n(),
         num_controls_in_test_fold = sum(Diagnosis == "Control"),
         control_prop_in_fold = num_controls_in_test_fold/num_subjects_in_test_fold) %>%
  group_by(Study, Comparison_Group, Control_Prop) %>%
  summarise(mean_control_prop_across_folds = 100*mean(control_prop_in_fold),
            min_control_prop_across_folds = 100*min(control_prop_in_fold),
            max_control_prop_across_folds = 100*max(control_prop_in_fold),
            SD_control_prop_across_folds = 100*sd(control_prop_in_fold))


# UCLA ADHD vs. control
UCLA_CNP_ADHD_fold_assignments %>%
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
UCLA_CNP_bipolar_fold_assignments %>%
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
ABIDE_ASD_fold_assignments %>%
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
  