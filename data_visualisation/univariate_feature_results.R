################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/Library/CloudStorage/OneDrive-TheUniversityofSydney(Students)/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/univariate_results/")
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
                             Group_Nickname = c("SCZ", "BP", "ADHD", "ASD"))
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
library(patchwork)
library(broom)
library(colorspace)
library(scales)
library(correctR)
library(splitstackshape)
library(ComplexHeatmap)
library(circlize)
theme_set(theme_cowplot())

# Source visualisation script
source(glue("{github_dir}/data_visualisation/Manuscript_Draft_Visualisations_Helper.R"))

UCLA_CNP_brain_region_info <- read.csv("~/data/UCLA_CNP/study_metadata/UCLA_CNP_Brain_Region_info.csv")
ABIDE_ASD_brain_region_info <- read.table("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_Harvard_Oxford_cort_prob_2mm_ROI_lookup.txt", sep=";", header = T) %>%
  mutate(Brain_Region = ifelse(Index==45, "Heschl's Gyrus (includes H1 and H2)", Brain_Region))

# Load in univariate time-series feature info
TS_feature_info <- read.csv(glue("{github_dir}/data_visualisation/catch25_info.csv"))

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

# Load null distribution
univariate_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_null_balanced_accuracy_distributions.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)

# Load study metadata
UCLA_CNP_metadata <- pyarrow_feather$read_feather("~/data/UCLA_CNP/study_metadata/UCLA_CNP_sample_metadata.feather") %>%
  mutate(Study="UCLA_CNP")
ABIDE_ASD_metadata <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_sample_metadata.feather")  %>%
  mutate(Study="ABIDE_ASD")


################################################################################
# Balanced accuracy heatmap
univariate_p_values %>%
  filter(Univariate_Feature_Set == univariate_feature_set,
         Analysis_Type == "Univariate_TS_Feature") %>%
  dplyr::rename("feature_name" = "group_var") %>%
  left_join(., TS_feature_info) %>%
  filter(p_value_Bonferroni<0.05) %>%
  mutate(Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                      Comparison_Group == "Bipolar" ~ "BP",
                                      T ~ Comparison_Group),
         Balanced_Accuracy_Across_Folds = 100*Balanced_Accuracy_Across_Folds) %>%
  mutate(Figure_name = fct_reorder(Figure_name, Balanced_Accuracy_Across_Folds, .fun=sum),
         Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BP", "ADHD", "ASD"))) %>%
  ggplot(data=., mapping=aes(x=Comparison_Group, y=Figure_name, 
                             fill=Balanced_Accuracy_Across_Folds)) +
  geom_tile()+
  geom_text(aes(label = round(Balanced_Accuracy_Across_Folds, 1))) +
  scale_fill_gradientn(colors=c(alpha("#4C7FC0", 0.1), "#4C7FC0"), 
                       na.value=NA) +
  labs(fill = "Mean Balanced Accuracy (%)") +
  xlab("Clinical Group") +
  ylab("Univariate catch25 feature") +
  theme(legend.position="none",
        strip.background = element_blank(),
        strip.text = element_blank())
ggsave(glue("{plot_path}/Feature_wise_results.svg"),
       width=4.5, height=6.5, units="in", dpi=300)


################################################################################
# Create adjacency matrix of feature values across all participants

top10_features <- univariate_p_values %>%
  filter(Univariate_Feature_Set == univariate_feature_set,
         Analysis_Type == "Univariate_TS_Feature") %>%
  filter(p_value_Bonferroni<0.05) %>%
  group_by(group_var) %>%
  summarise(sum_balacc = sum(Balanced_Accuracy_Across_Folds)) %>%
  slice_max(sum_balacc, n=10) %>%
  pull(group_var)

data_for_corr_heatmap <- plyr::rbind.fill(UCLA_CNP_catch25, ABIDE_ASD_catch25) %>%
  filter(names %in% top10_features) %>%
  mutate(unique_ID = paste0(Sample_ID, "__", Brain_Region), .keep="unused") %>%
  left_join(., TS_feature_info, by=c("names"="feature_name")) %>%
  dplyr::select(unique_ID, Figure_name, values) %>%
  pivot_wider(id_cols=unique_ID, names_from=Figure_name, values_from=values) %>%
  dplyr::select(-unique_ID) %>%
  cor(method="spearman", use="complete.obs") %>%
  abs()

# Convert to long for easy querying
data_for_corr_long <- data_for_corr_heatmap %>%
  as.data.frame() %>%
  rownames_to_column(var="feature1") %>%
  pivot_longer(cols=c(-feature1), names_to="feature2", values_to="spearman_corr_abs")

num_branches <- 3

ht1 <- ComplexHeatmap::Heatmap(data_for_corr_heatmap,
                               clustering_distance_rows = "spearman",
                               clustering_distance_columns = "spearman",
                               clustering_method_rows = "average",
                               clustering_method_columns = "average",
                               row_names_side = "right",
                               row_dend_side = "left", 
                               row_dend_width = unit(3, "cm"),
                               row_dend_gp = gpar(lwd=unit(3, "cm")),
                               row_split = num_branches,
                               column_split = num_branches,
                               column_gap = unit(3, "mm"),
                               row_title = NULL,
                               column_title = NULL,
                               show_row_names = TRUE,
                               show_column_names = FALSE,
                               show_column_dend = FALSE,
                               col = c("white", RColorBrewer::brewer.pal(12, "Reds")),
                               name = "Spearman rank correlation, \u03c1",
                               heatmap_legend_param = list(legend_direction = "horizontal",
                                                           legend_width = unit(5, "cm")))

svg(glue("{plot_path}/catch25_top_feature_spearman_corr.svg"),
    width=7, height=6, bg=NA)
draw(ht1, heatmap_legend_side = "bottom",
     background = "transparent")
dev.off()

################################################################################
# How are different disorders correlated in terms of classification accuracy by catch25 feature?

feature_accuracy_corr_data <- univariate_p_values %>%
  filter(Analysis_Type=="Univariate_TS_Feature") %>%
  left_join(., study_group_df) %>%
  dplyr::select(group_var, Group_Nickname, Balanced_Accuracy_Across_Folds) %>%
  pivot_wider(id_cols=group_var, names_from=Group_Nickname, values_from=Balanced_Accuracy_Across_Folds) %>%
  dplyr::select(-group_var) %>%
  cor(method="spearman")

num_branches <- 3
ht2 <- ComplexHeatmap::Heatmap(feature_accuracy_corr_data,
                               clustering_distance_rows = "spearman",
                               clustering_distance_columns = "spearman",
                               clustering_method_rows = "average",
                               clustering_method_columns = "average",
                               row_names_side = "right",
                               row_dend_side = "left", 
                               row_dend_width = unit(3, "cm"),
                               row_dend_gp = gpar(lwd=unit(3, "cm")),
                               row_split = num_branches,
                               row_gap = unit(3, "mm"),
                               column_split = num_branches,
                               column_gap = unit(3, "mm"),
                               row_title = NULL,
                               column_title = NULL,
                               show_row_names = TRUE,
                               show_column_names = FALSE,
                               show_column_dend = FALSE,
                               cell_fun = function(j, i, x, y, width, height, fill) {
                                 grid.text(sprintf("%.2f", feature_accuracy_corr_data[i, j]), x, y, gp = gpar(fontsize = 16))
                               },
                               col = colorRamp2(seq(-1, 1, length.out=9), brewer_pal("div", "RdBu", direction=-1)(9)),
                               name = "Spearman rank correlation, \u03c1",
                               heatmap_legend_param = list(legend_direction = "horizontal",
                                                           legend_width = unit(5, "cm")))

svg(glue("{plot_path}/catch25_feature_performance_correlations.svg"),
    width=5, height=4, bg=NA)
draw(ht2, heatmap_legend_side = "bottom",
     background = "transparent")
dev.off()