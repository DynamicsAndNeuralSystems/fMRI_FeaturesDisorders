# python_to_use <- "~/.conda/envs/pyspi/bin/python3"
python_to_use <- "/Users/abry4213/anaconda3/envs/pyspi/bin/python3"
reticulate::use_python(python_to_use)

library(tidyverse)
library(reticulate)
library(icesTAF)
library(cowplot)
library(ggpubr)
library(ggsignif)
library(feather)
library(glue)
library(broom)
library(patchwork)
theme_set(theme_cowplot())

# Import pyarrow.feather as pyarrow_feather
pyarrow_feather <- import("pyarrow.feather")

# DIY rlist::list.append
list.append <- function (.data, ...) 
{
  if (is.list(.data)) {
    c(.data, list(...))
  }
  else {
    c(.data, ..., recursive = FALSE)
  }
}

################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/Library/CloudStorage/OneDrive-TheUniversityofSydney(Students)/github/fMRI_FeaturesDisorders/"

plot_path <- paste0(github_dir, "plots/Manuscript_Draft/volumetric_analysis/")
TAF::mkdir(plot_path)
univariate_feature_set <- "catch24"
SVM_kernel <- "linear"

UCLA_CNP_data_path <- "~/data/UCLA_CNP"
data_path <- "~/data/TS_feature_manuscript"

# Read in metadata
UCLA_CNP_sample_metadata <- feather::read_feather(glue("{UCLA_CNP_data_path}/study_metadata/UCLA_CNP_sample_metadata.feather"))
UCLA_CNP_brain_region_info <- read.csv("~/data/UCLA_CNP/study_metadata/UCLA_CNP_Brain_Region_info.csv")
aparc_aseg_LUT <- read.table("~/data/neuroimaging_atlases/FreeSurferLUT.txt",
                             header=T) %>%
  dplyr::rename("ROI_Index" = "Value")

study_group_df <- data.frame(Study = rep("UCLA_CNP", 3),
                             Noise_Proc = rep("AROMA+2P+GMR",3),
                             Comparison_Group = c("Schizophrenia", "Bipolar", "ADHD"),
                             Group_Nickname = c("SCZ", "BPD", "ADHD"))

# Read in region-wise volumes
region_wise_volumes <- pyarrow_feather$read_feather(glue("{UCLA_CNP_data_path}/processed_data/aparc_aseg_BOLD_space_voxel_volumes.feather")) %>%
  left_join(., UCLA_CNP_sample_metadata) %>%
  left_join(., aparc_aseg_LUT) %>%
  left_join(., UCLA_CNP_brain_region_info) %>%
  filter(!is.na(Index))

# Fit OLS models to extract beta coefficient for regional volumes in each group relative to control
run_lm_beta_stats_for_group <- function(comparison_group, region_wise_volumes){
  res <- region_wise_volumes %>%
    filter(Diagnosis %in% c(comparison_group, "Control")) %>%
    dplyr::select(Brain_Region, Diagnosis, Num_Voxels) %>%
    mutate(Diagnosis = factor(Diagnosis, levels = c("Control", comparison_group))) %>%
    group_by(Brain_Region) %>%
    nest() %>%
    mutate(
      fit = map(data, ~ lm(Num_Voxels ~ Diagnosis, data = .x)),
      tidied = map(fit, tidy)
    ) %>% 
    unnest(tidied) %>%
    dplyr::select(-data, -fit) %>%
    ungroup() %>%
    filter(term != "(Intercept)") %>%
    mutate(Comparison_Group = comparison_group)
  
  return(res)
}

ROI_volume_beta_by_group <- 1:3 %>%
  purrr::map_df(~ run_lm_beta_stats_for_group(region_wise_volumes = region_wise_volumes,
                                          comparison_group = study_group_df$Comparison_Group[.x])) %>%
  group_by(Comparison_Group) %>%
  mutate(p_value_Bonferroni = p.adjust(p.value, method="bonferroni"))

# Load univariate classification results across all folds
univariate_balanced_accuracy_AUC_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_mixedsigmoid_scaler_balanced_accuracy_AUC_all_folds.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set, kernel==SVM_kernel)
# Compute mean + SD performance across all folds
univariate_balanced_accuracy <- univariate_balanced_accuracy_AUC_all_folds %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Analysis_Type, group_var, kernel) %>%
  reframe(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
          Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T),
          ROC_AUC_Across_Folds = mean(ROC_AUC, na.rm=T),
          ROC_AUC_Across_Folds_SD = sd(ROC_AUC, na.rm=T))
# Load p-values
univariate_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_mixedsigmoid_scaler_empirical_p_values.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set) %>%
  dplyr::select(-Balanced_Accuracy_Across_Repeats, -Balanced_Accuracy_Across_Repeats_SD, 
                -ROC_AUC_Across_Repeats, -ROC_AUC_Across_Repeats_SD) %>%
  left_join(., univariate_balanced_accuracy)


# Plot region-wise volume beta coefficients vs balanced accuracy
all_regions_plot <- ROI_volume_beta_by_group %>%
  dplyr::select(Brain_Region, Comparison_Group, estimate) %>%
  dplyr::rename("beta_coef" = "estimate") %>%
  left_join(., univariate_p_values %>% dplyr::rename("Brain_Region" = "group_var")) %>%
  ggplot(data=., mapping=aes(x=beta_coef, y=100*Balanced_Accuracy_Across_Folds, color=Comparison_Group)) +
  geom_point() +
  facet_grid(Comparison_Group ~ ., scales="free", switch="both") +
  ylab("Mean Balanced Accuracy for Region") +
  xlab("\u03b2 for Region Volume") +
  stat_cor(method="spearman", cor.coef.name="rho", color="black", label.sep = "\n") +
  stat_smooth(method="lm", color="black") +
  scale_color_manual(values=c("Schizophrenia" = "#573DC7", 
                              "Bipolar" = "#D5492A", 
                              "ADHD" = "#0F9EA9")) +
  theme(legend.position = "none")


# Re-plot only with significant regions
sig_regions_plot <- ROI_volume_beta_by_group %>%
  dplyr::select(Brain_Region, Comparison_Group, estimate) %>%
  dplyr::rename("beta_coef" = "estimate") %>%
  left_join(., univariate_p_values %>% dplyr::rename("Brain_Region" = "group_var")) %>%
  filter(p_value_Bonferroni < 0.05) %>%
  ggplot(data=., mapping=aes(x=beta_coef, y=100*Balanced_Accuracy_Across_Folds, color=Comparison_Group)) +
  geom_point() +
  facet_grid(Comparison_Group ~ ., scales="free", switch="both") +
  ylab("Mean Balanced Accuracy for Region") +
  xlab("\u03b2 for Region Volume") +
  stat_cor(method="spearman", cor.coef.name="rho", color="black", label.sep = "\n") +
  stat_smooth(method="lm", color="black") +
  scale_color_manual(values=c("Schizophrenia" = "#573DC7", 
                              "Bipolar" = "#D5492A", 
                              "ADHD" = "#0F9EA9")) +
  theme(legend.position = "none")

all_regions_plot + sig_regions_plot
ggsave(glue("{plot_path}/Volume_vs_BalAcc_res.svg"),
       width=6, height=6, units="in", dpi=300)


# Plot volume in the left vs right hemispheres by region/condition
region_wise_volumes %>%
  mutate(Hemisphere = case_when(str_detect(Brain_Region, "Left|lh-") ~ "Left",
                                str_detect(Brain_Region, "Right|rh-") ~ "Right")) %>%
  mutate(Brain_Region = gsub("Left-|ctx-lh-|Right-|ctx-rh-", "", Brain_Region)) %>%
  dplyr::select(Sample_ID, Diagnosis, Brain_Region, Hemisphere, Num_Voxels) %>%
  pivot_wider(id_cols = c(Sample_ID, Diagnosis, Brain_Region),
              names_from = Hemisphere,
              values_from = Num_Voxels) %>%
  ggplot(data=., mapping=aes(x=Left, y=Right, color=Diagnosis)) +
  geom_point() +
  theme(legend.position="none") +
  geom_abline(slope=1, intercept=0, color="black") +
  facet_wrap(Diagnosis ~ .)


# Are there any left-right volumetric differences by cohort?
brain_regions_we_used <- univariate_p_values %>%
  filter(Analysis_Type == "Univariate_Brain_Region",
         Study=="UCLA_CNP") %>%
  distinct(group_var) %>%
  mutate(Hemisphere = case_when(str_detect(group_var, "Left|lh-") ~ "Left",
                                str_detect(group_var, "Right|rh-") ~ "Right")) %>%
  mutate(Brain_Region = gsub("Left-|ctx-lh-|Right-|ctx-rh-", "", group_var)) %>%
  distinct(Brain_Region) %>%
  pull(Brain_Region)

left_right_paired_T_results <- region_wise_volumes %>%
  mutate(Hemisphere = case_when(str_detect(Brain_Region, "Left|lh-") ~ "Left",
                                str_detect(Brain_Region, "Right|rh-") ~ "Right")) %>%
  mutate(Brain_Region = gsub("Left-|ctx-lh-|Right-|ctx-rh-", "", Brain_Region)) %>%
  dplyr::select(Sample_ID, Diagnosis, Brain_Region, Hemisphere, Num_Voxels) %>%
  filter(!is.na(Hemisphere) & Brain_Region %in% brain_regions_we_used) %>%
  group_by(Diagnosis, Brain_Region) %>%
  nest() %>%
  mutate(
    test = map(data, ~ t.test(Num_Voxels ~ Hemisphere, data=.x, paired=TRUE)), # S3 list-col
    tidied = map(test, tidy)
  ) %>%
  unnest(tidied) %>%
  ungroup() %>%
  group_by(Diagnosis) %>%
  mutate(p_val_Bonferroni = p.adjust(p.value, method="bonferroni"))
