################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/Library/CloudStorage/OneDrive-TheUniversityofSydney(Students)/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/analysis_summary/")
TAF::mkdir(plot_path)

# python_to_use <- "~/.conda/envs/pyspi/bin/python3"
python_to_use <- "/Users/abry4213/anaconda3/envs/pyspi/bin/python3"
pairwise_feature_set <- "pyspi14"
univariate_feature_set <- "catch25"
SVM_kernel <- "Linear"
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
library(see)
theme_set(theme_cowplot())

# Source visualisation script
source(glue("{github_dir}/data_visualisation/Manuscript_Draft_Visualisations_Helper.R"))

UCLA_CNP_brain_region_info <- read.csv("~/data/UCLA_CNP/study_metadata/UCLA_CNP_Brain_Region_info.csv")
ABIDE_ASD_brain_region_info <- read.table("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_Harvard_Oxford_cort_prob_2mm_ROI_lookup.txt", sep=";", header = T) %>%
  mutate(Brain_Region = ifelse(Index==45, "Heschl's Gyrus (includes H1 and H2)", Brain_Region))

# Load in univariate time-series feature info
TS_feature_info <- read.csv(glue("{github_dir}/data_visualisation/catch25_info.csv"))


# Load study metadata
UCLA_CNP_metadata <- pyarrow_feather$read_feather("~/data/UCLA_CNP/study_metadata/UCLA_CNP_sample_metadata.feather") %>%
  mutate(Study="UCLA_CNP")
ABIDE_ASD_metadata <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_sample_metadata.feather")  %>%
  mutate(Study="ABIDE_ASD")

# Univariate results
univariate_balanced_accuracy <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_balanced_accuracy_all_folds.feather")) %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Analysis_Type, group_var, kernel) %>%
  reframe(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
          Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T))
univariate_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_empirical_p_values.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set) %>%
  dplyr::select(-Balanced_Accuracy_Across_Folds) %>%
  left_join(., univariate_balanced_accuracy)

# Pairwise results
pairwise_balanced_accuracy <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_balanced_accuracy_all_folds.feather")) %>%
  group_by(Study, Comparison_Group, Pairwise_Feature_Set, Analysis_Type, group_var) %>%
  reframe(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
          Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T))
pairwise_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_empirical_p_values.feather")) %>%
  dplyr::select(-Balanced_Accuracy_Across_Folds) %>%
  left_join(., pairwise_balanced_accuracy)

# Univariate+Pairwise Combo Results
combo_univariate_pairwise_balanced_accuracy <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combined_univariate_pairwise_balanced_accuracy_all_folds.feather")) %>%
  mutate(Analysis_Type = "SPI_Univariate_Combo") %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Pairwise_Feature_Set, Analysis_Type, group_var) %>%
  reframe(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
          Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T))
combo_univariate_pairwise_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combined_univariate_pairwise_empirical_p_values.feather")) %>%
  dplyr::select(-Balanced_Accuracy_Across_Folds) %>%
  left_join(., combo_univariate_pairwise_balanced_accuracy)

# Merge the p-value results together
all_p_values <- do.call(plyr::rbind.fill, list(univariate_p_values,
                                               pairwise_p_values,
                                               combo_univariate_pairwise_p_values)) %>%
  mutate(Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                      Comparison_Group == "Bipolar" ~ "BPD",
                                      T ~ Comparison_Group),
         Balanced_Accuracy_Across_Folds = 100*Balanced_Accuracy_Across_Folds) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD")))
  
################################################################################
# Violin plot by disorder
max_values_by_disorder <- all_p_values %>%
  group_by(Comparison_Group) %>%
  summarise(max_balacc = max(Balanced_Accuracy_Across_Folds))

all_p_values %>%
  mutate(Analysis_Type = factor(Analysis_Type, levels = c("Univariate_Brain_Region",
                                                          "Univariate_TS_Feature",
                                                          "Univariate_Combo",
                                                          "Pairwise_SPI",
                                                          "SPI_Univariate_Combo"))) %>%
  ggplot(data=., mapping=aes(x=Analysis_Type, y=Balanced_Accuracy_Across_Folds)) +
  
  geom_hline(yintercept=50, linetype=2, color="gray50") + # Null line
  geom_hline(data = max_values_by_disorder, aes(yintercept = max_balacc), linetype=2, size=0.8) + # Top performance line
  geom_violinhalf(aes(fill=Analysis_Type), scale="width", 
                  position = position_nudge(x=0.1))  +
  geom_boxplot(width=0.1, notch=FALSE, notchwidth = 0.4, outlier.shape = NA,
               fill=NA, color="black",
               position = position_nudge(x=0.15), coef = 0) +
  geom_point(aes(color = Analysis_Type), position = position_jitterdodge(dodge.width = 1,
                                                                     jitter.width = 0.5),
             size = 1.7, alpha=0.6, stroke=0) +
  facet_grid(Comparison_Group ~ ., switch="both") +
  ylab("Balanced Accuracy (%)") +
  xlab("Representation Type") +
  scale_color_manual(values = c("Univariate_Brain_Region" = "#FF003B",
                               "Univariate_TS_Feature" = "#3F84C6",
                               "Univariate_Combo" = "#C083CB",
                               "Pairwise_SPI" = "#4FCD1F",
                               "SPI_Univariate_Combo" = "#FCB33B")) +
  scale_fill_manual(values = c("Univariate_Brain_Region" = "#FF003B",
                               "Univariate_TS_Feature" = "#3F84C6",
                               "Univariate_Combo" = "#C083CB",
                               "Pairwise_SPI" = "#4FCD1F",
                               "SPI_Univariate_Combo" = "#FCB33B")) +
  theme(strip.placement="outside",
        strip.background=element_blank(),
        strip.text.y.left = element_text(angle=0, face="bold"),
        legend.position="none")
ggsave(paste0(plot_path, "Performance_across_representations.svg"),
       width = 6, height= 6, units="in", dpi=300)