library(tidyverse)
python_to_use <- "/Users/abry4213/anaconda3/envs/pyspi/bin/python3"
reticulate::use_python(python_to_use)

library(cowplot)
theme_set(theme_cowplot())

library(reticulate)
library(see)
library(glue)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(FactoMineR)
library("factoextra")

# Import pyarrow.feather as pyarrow_feather
pyarrow_feather <- import("pyarrow.feather")

# Define data paths
github_dir <- "~/Library/CloudStorage/OneDrive-TheUniversityofSydney(Students)/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/methods_supplement/SVM_weighting_analysis/")
TAF::mkdir(plot_path)
data_path <- "~/data/TS_feature_manuscript"
UCLA_data_path <- "~/data/UCLA_CNP"
ABIDE_data_path <- "~/data/ABIDE_ASD"

# Load catch25 data
UCLA_CNP_sample_metadata <- feather::read_feather(glue("{UCLA_data_path}/study_metadata/UCLA_CNP_sample_metadata.feather"))
UCLA_CNP_catch25 <- pyarrow_feather$read_feather(glue("{UCLA_data_path}/processed_data/UCLA_CNP_AROMA_2P_GMR_catch25_filtered.feather"))  %>%
  left_join(., UCLA_CNP_sample_metadata)
ABIDE_sample_metadata <- feather::read_feather(glue("{ABIDE_data_path}/study_metadata/ABIDE_ASD_sample_metadata.feather"))
ABIDE_catch25 <- pyarrow_feather$read_feather(glue("{ABIDE_data_path}/processed_data/ABIDE_ASD_FC1000_catch25_filtered.feather"))  %>%
  left_join(., ABIDE_sample_metadata)

# Load SVM data with vs without inverse probability weighting
SVM_balacc_with_vs_without_inv_prob <- pyarrow_feather$read_feather(glue("{data_path}/SVM_with_vs_without_inv_prob_weighting.feather")) 

# For each brain region, compare the difference in mean balanced accuracy with versus without inverse probability weighting
SVM_balacc_with_vs_without_inv_prob %>%
  group_by(Comparison_Group, Analysis_Type, group_var, Weighting_Type) %>%
  summarise(Mean_Balanced_Accuracy = mean(Balanced_Accuracy)) %>%
  group_by(Comparison_Group, Analysis_Type, group_var) %>%
  summarise(Weighting_Diff = Mean_Balanced_Accuracy[Weighting_Type=="Balanced"] - Mean_Balanced_Accuracy[Weighting_Type=="None"]) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("Schizophrenia", "Bipolar", "ADHD"))) %>%
  ungroup() %>%
  mutate(group_var = fct_reorder(group_var, Weighting_Diff, .fun=mean)) %>%
  ggplot(data=., mapping=aes(x=Comparison_Group, y=group_var, fill=Weighting_Diff)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red", limits=c(-0.15, 0.15)) +
  ylab("Representation for SVM") +
  xlab("Disorder") +
  facet_grid(Analysis_Type ~ ., scales="free", space="free", switch="both") +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size=8),
        strip.placement="outside",
        axis.text.x = element_text(face="bold"),
        strip.text.y.left = element_text(angle=0, face="bold"),
        strip.background = element_blank())

ggsave(paste0(plot_path, "Linear_SVM_Performance_based_on_weighting_type.svg"),
       width = 7, height=10, units="in", dpi=300)


################################################################################

# Heatmap of low freq power values across the brain in ADHD vs controls
data_for_heatmap <- UCLA_CNP_catch25 %>%
  filter(Diagnosis %in% c("Control", "ADHD"), names=="SP_Summaries_welch_rect_area_5_1") %>%
  dplyr::select(Sample_ID, Diagnosis, Brain_Region, values) %>%
  pivot_wider(id_cols = c(Sample_ID, Diagnosis), names_from = Brain_Region, values_from=values)

mat_data <- data_for_heatmap %>%
  dplyr::select(-Diagnosis) %>%
  column_to_rownames(var="Sample_ID") %>%
  as.matrix()

row_ha = rowAnnotation(foo2 = runif(10), bar2 = anno_barplot(runif(10)))

ha = rowAnnotation(df = data_for_heatmap %>% pull(Diagnosis),
                   col = list(df = c("Control" = "red", "ADHD" = "green")
                   )
)

ht2 <- ComplexHeatmap::Heatmap(mat_data,
                               clustering_distance_rows = "spearman",
                               clustering_distance_columns = "spearman",
                               clustering_method_rows = "average",
                               clustering_method_columns = "average",
                               row_names_side = "right",
                               row_dend_side = "left", 
                               row_dend_width = unit(3, "cm"),
                               row_dend_gp = gpar(lwd=unit(1, "cm")),
                               row_title = "Subjects",
                               column_title = "Brain Regions",
                               show_row_names = TRUE,
                               show_column_names = FALSE,
                               show_column_dend = FALSE,
                               right_annotation = ha,
                               # col = colorRamp2(seq(0, 0.9, length.out=9), scales::brewer_pal("div", "RdBu", direction=-1)(9)),
                               name = "low_freq_power values across the brain",
                               heatmap_legend_param = list(legend_direction = "horizontal",
                                                           legend_width = unit(5, "cm")))

draw(ht2, heatmap_legend_side = "bottom",
     background = "transparent")


###########################################################################
# Running nulls with vs. without inverse probability weighting

ADHD_low_freq_power_main <-  pyarrow_feather$read_feather(glue("{data_path}/ADHD_low_freq_power_main_balacc_res.feather")) 

ADHD_low_freq_power_main_linear_mean <- ADHD_low_freq_power_main %>%
  filter(Kernel=="linear") %>%
  group_by(SVM_Weighting, group_var) %>%
  summarise(Balanced_Accuracy_Across_Folds = 100*mean(Balanced_Accuracy))

ADHD_low_freq_power_nulls <- pyarrow_feather$read_feather(glue("{data_path}/ADHD_low_freq_power_nulls.feather"))

ADHD_low_freq_power_nulls %>%
  ggplot(data=., mapping=aes(x=100*Null_Balanced_Accuracy, fill=SVM_Weighting)) +
  geom_histogram(alpha=0.4, position="identity") +
  geom_vline(data = ADHD_low_freq_power_main_linear_mean, aes(xintercept=Balanced_Accuracy_Across_Folds, 
                                                  color=SVM_Weighting),
             linewidth=1.5) + 
  scale_color_manual(values=c("Balanced" = "#169412",
                              "None" = "#B965D1")) +
  scale_fill_manual(values=c("Balanced" = "#169412",
                              "None" = "#B965D1")) +
  ylab("# Null Iterations") +
  xlab("Balanced Accuracy (%)") +
  labs(color="SVM Weighting", fill="SVM Weighting") +
  theme(legend.position = "bottom")


###########################################################################
# Linear kernel vs RBF kernel

ADHD_low_freq_power_main %>%
  mutate(Kernel_Weighting = paste0(Kernel, " kernel, ", SVM_Weighting)) %>%
  ggplot(data=., mapping=aes(x=Kernel_Weighting, y=100*Balanced_Accuracy)) +
  theme(legend.position="none") +
  geom_violinhalf(aes(fill=Kernel_Weighting), scale="width", 
                  position = position_nudge(x=0.1))  +
  geom_hline(yintercept = 50, linetype=2, color="black") +
  geom_boxplot(width=0.1, notch=FALSE, notchwidth = 0.4, outlier.shape = NA,
               fill=NA, color="white",
               position = position_nudge(x=0.16), coef = 0) +
  geom_point(aes(color = Kernel_Weighting), position = position_jitterdodge(dodge.width = 1,
                                                                            jitter.width = 0.2),
             size = 1, alpha=0.7) +
  ylab("Balanced Accuracy (%)") +
  xlab("Kernel Type + Weighting Type") +
  coord_flip()

###########################################################################
# Normalization based on training data versus full data

ADHD_low_freq_power_normalization_comparison <- pyarrow_feather$read_feather(glue("{data_path}/ADHD_low_freq_power_normalization.feather"))

ADHD_low_freq_power_normalization_comparison %>%
  mutate(Repeat_Fold = paste0(Repeat_Number,"__",Fold),
         Normalization = factor(Normalization, levels=c("train", "full"))) %>%
  ggplot(data=., mapping=aes(x=Normalization, y=100*Balanced_Accuracy, fill=Normalization)) +
  geom_violin(alpha=0.4) +
  geom_line(aes(group=Repeat_Fold), alpha=0.2) + 
  stat_summary(fun=mean, geom="crossbar", linewidth=0.7, aes(color=Normalization))  + 
  scale_color_manual(values=c("full" = "#169412",
                              "train" = "#B965D1")) +
  scale_fill_manual(values=c("full" = "#169412",
                             "train" = "#B965D1")) +
  ylab("Balanced Accuracy (%)") +
  xlab("Normalization Parameter Learning Dataset") +
  theme(legend.position="none")

# Mean balanced accuracy by model type
ADHD_low_freq_power_normalization_comparison %>%
  group_by(Normalization) %>%
  summarise(mean_balacc=100*mean(Balanced_Accuracy))




###########################################################################

# PCA based on low_freq_power feature
# Load PCA scores
ADHD_low_freq_power_PCA_scores <- pyarrow_feather$read_feather(glue("{data_path}/ADHD_low_freq_power_pca_scores.feather"))

ADHD_low_freq_power_PCA_scores %>%
  ggplot(data=., mapping=aes(x=PC1, y=PC2,color=Diagnosis)) +
  geom_vline(linetype=2, xintercept=0) +
  geom_hline(linetype=2, yintercept=0)+
  ggtitle("First two components from low_freq_power PCA") +
  geom_point() +
  theme(plot.title = element_text(hjust=0.5))
 
# Load classification results
ADHD_low_freq_power_PCA_classification_res <- pyarrow_feather$read_feather(glue("{data_path}/ADHD_low_freq_power_pca_classification_res.feather"))

ADHD_low_freq_power_PCA_classification_res %>%
  mutate(Repeat_Fold = paste0(Repeat_Number,"__",Fold)) %>%
  ggplot(data=., mapping=aes(x=SVM_Weighting, y=100*Balanced_Accuracy, fill=SVM_Weighting)) +
  geom_violin(alpha=0.4) +
  geom_line(aes(group=Repeat_Fold), alpha=0.2) + 
  stat_summary(fun=mean, geom="crossbar", linewidth=0.7, aes(color=SVM_Weighting))  + 
  scale_color_manual(values=c("Balanced" = "#169412",
                              "None" = "#B965D1")) +
  scale_fill_manual(values=c("Balanced" = "#169412",
                             "None" = "#B965D1")) +
  ylab("Balanced Accuracy (%)") +
  xlab("SVM Weighting Type") +
  theme(legend.position="none")

# Take the means
ADHD_low_freq_power_PCA_classification_res %>%
  group_by(SVM_Weighting) %>%
  summarise(mean=mean(Balanced_Accuracy))