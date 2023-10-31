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
univariate_feature_set <- "catch25"
SVM_kernel <- "Linear"

UCLA_CNP_data_path <- "~/data/UCLA_CNP"
data_path <- "~/data/TS_feature_manuscript"

# Load in univariate time-series feature info
TS_feature_info <- read.csv(glue("{github_dir}/data_visualisation/catch25_info.csv"))

# Read in metadata
UCLA_CNP_sample_metadata <- feather::read_feather(glue("{UCLA_CNP_data_path}/study_metadata/UCLA_CNP_sample_metadata.feather"))
UCLA_CNP_brain_region_info <- read.csv("~/data/UCLA_CNP/study_metadata/UCLA_CNP_Brain_Region_info.csv")
aparc_aseg_LUT <- read.table("~/data/neuroimaging_atlases/FreeSurferLUT.txt",
                             header=T) %>%
  dplyr::rename("ROI_Index" = "Value")

study_group_df <- data.frame(Study = rep("UCLA_CNP", 3),
                             Noise_Proc = rep("AROMA+2P+GMR",3),
                             Comparison_Group = c("Schizophrenia", "Bipolar", "ADHD"),
                             Group_Nickname = c("SCZ", "BP", "ADHD"))

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
univariate_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_balanced_accuracy_all_folds.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set, kernel==SVM_kernel)
# Compute mean + SD performance across all folds
univariate_balanced_accuracy <- univariate_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Analysis_Type, group_var, kernel) %>%
  reframe(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
          Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T))

# Load p-values
univariate_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_empirical_p_values.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set) %>%
  dplyr::select(-Balanced_Accuracy_Across_Folds) %>%
  left_join(., univariate_balanced_accuracy)

################################################################################
# Plot region-wise volume beta coefficients vs balanced accuracy
ROI_volume_beta_by_group %>%
  dplyr::select(Brain_Region, Comparison_Group, estimate) %>%
  dplyr::rename("beta_coef" = "estimate") %>%
  left_join(., univariate_p_values %>% dplyr::rename("Brain_Region" = "group_var")) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("Schizophrenia", "Bipolar", "ADHD"))) %>%
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
  theme(legend.position = "none",
        panel.border = element_blank(),
        strip.background = element_blank(), 
        strip.placement = "outside")
ggsave(glue("{plot_path}/Volume_vs_BalAcc_res.svg"),
       width=3, height=6, units="in", dpi=300)

# Do correlation test with Bonferroni correction
ROI_volume_beta_by_group %>%
  dplyr::select(Brain_Region, Comparison_Group, estimate) %>%
  dplyr::rename("beta_coef" = "estimate") %>%
  left_join(., univariate_p_values %>% dplyr::rename("Brain_Region" = "group_var")) %>%
  group_by(Comparison_Group) %>%
  do(tidy(cor.test(.$beta_coef, .$Balanced_Accuracy_Across_Folds, method="spearman"))) %>%
  ungroup() %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("Schizophrenia", "Bipolar", "ADHD")),
         p_val_Bonferroni = p.adjust(p.value, method="bonferroni")) %>%
  arrange(Comparison_Group)

################################################################################
# Does the average volume of a region relate to its balanced accuracy?
region_wise_volumes %>%
  group_by(Diagnosis, Brain_Region) %>%
  summarise(mean_volume = mean(Num_Voxels)) %>%
  filter(Diagnosis != "Control") %>%
  dplyr::rename("Comparison_Group"="Diagnosis") %>%
  left_join(., univariate_p_values %>% dplyr::rename("Brain_Region" = "group_var")) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("Schizophrenia", "Bipolar", "ADHD"))) %>%
  ggplot(data=., mapping=aes(x=mean_volume, y=100*Balanced_Accuracy_Across_Folds, color=Comparison_Group)) +
  geom_point() +
  facet_grid(Comparison_Group ~ ., scales="free", switch="both") +
  ylab("Mean Balanced Accuracy for Region") +
  xlab("Mean Region Volume") +
  stat_cor(method="spearman", cor.coef.name="rho", color="black", label.sep = "\n",
           label.x.npc="right") +
  stat_smooth(method="lm", color="black") +
  scale_color_manual(values=c("Schizophrenia" = "#573DC7",
                              "Bipolar" = "#D5492A",
                              "ADHD" = "#0F9EA9")) +
  theme(legend.position = "none",
        panel.border = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside")
ggsave(glue("{plot_path}/Average_Volume_vs_BalAcc_res.svg"),
       width=3, height=6, units="in", dpi=300)


# Do correlation test with Bonferroni correction
region_wise_volumes %>%
  group_by(Diagnosis, Brain_Region) %>%
  summarise(mean_volume = mean(Num_Voxels)) %>%
  filter(Diagnosis != "Control") %>%
  dplyr::rename("Comparison_Group"="Diagnosis") %>%
  left_join(., univariate_p_values %>% dplyr::rename("Brain_Region" = "group_var")) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("Schizophrenia", "Bipolar", "ADHD"))) %>%
  group_by(Comparison_Group) %>%
  do(tidy(cor.test(.$mean_volume, .$Balanced_Accuracy_Across_Folds, method="spearman"))) %>%
  ungroup() %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("Schizophrenia", "Bipolar", "ADHD")),
         p_val_Bonferroni = p.adjust(p.value, method="bonferroni")) %>%
  arrange(Comparison_Group)

################################################################################
# Do any features correlate with general region size?
UCLA_CNP_catch25 <- pyarrow_feather$read_feather("~/data/UCLA_CNP/processed_data/UCLA_CNP_AROMA_2P_GMR_catch25_filtered.feather")  %>%
  left_join(., UCLA_CNP_sample_metadata)

UCLA_CNP_catch25 %>%
  left_join(., region_wise_volumes) %>%
  dplyr::select(Sample_ID, Brain_Region, names, Num_Voxels, values) %>%
  left_join(., TS_feature_info, by=c("names"="feature_name")) %>%
  filter(!is.na(Num_Voxels)) %>%
  group_by(Brain_Region, names, Figure_name) %>%
  summarise(mean_num_voxels = mean(Num_Voxels),
            mean_feature_values = mean(values)) %>%
  ggplot(data=., mapping=aes(x=mean_num_voxels, y=mean_feature_values, color=Figure_name)) +
  geom_point() + 
  ylab("Average feature value in region") +
  xlab("Average # voxels in region") +
  stat_smooth(se=FALSE, color="black", geom="line", alpha=0.7, size=1) +
  facet_wrap(Figure_name ~ ., nrow=5, scales="free_y") +
  theme(legend.position="none",
        strip.background = element_blank(),
        strip.text = element_text(face="bold"),
        axis.ticks.y = element_blank(),
        axis.text=element_blank())
ggsave(glue("{plot_path}/Volume_vs_Feature_Values.svg"),
       width=8, height=6, units="in", dpi=300)

