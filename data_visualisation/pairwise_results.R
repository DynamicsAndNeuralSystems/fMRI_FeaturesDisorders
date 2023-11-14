################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/Library/CloudStorage/OneDrive-TheUniversityofSydney(Students)/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/pairwise_results/")
TAF::mkdir(plot_path)

# python_to_use <- "~/.conda/envs/pyspi/bin/python3"
python_to_use <- "/Users/abry4213/anaconda3/envs/pyspi/bin/python3"
pairwise_feature_set <- "pyspi14"
univariate_feature_set <- "catch25"
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
library(patchwork)
library(ggseg)
library(broom)
library(colorspace)
library(scales)
library(correctR)
library(circlize)
library(ComplexHeatmap)
theme_set(theme_cowplot())

# Source visualisation script
source(glue("{github_dir}/data_visualisation/Manuscript_Draft_Visualisations_Helper.R"))

# Load in SPI info
SPI_info <- read.csv(glue("{github_dir}/data_visualisation/SPI_info.csv"))

# Load brain region info
UCLA_CNP_brain_region_info <- read.csv("~/data/UCLA_CNP/study_metadata/UCLA_CNP_Brain_Region_info.csv")
ABIDE_ASD_brain_region_info <- read.table("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_Harvard_Oxford_cort_prob_2mm_ROI_lookup.txt", sep=";", header = T) %>%
  mutate(Brain_Region = ifelse(Index==45, "Heschl's Gyrus (includes H1 and H2)", Brain_Region))

# Load participants included
UCLA_CNP_subjects_to_keep <- pyarrow_feather$read_feather("~/data/UCLA_CNP/processed_data/UCLA_CNP_filtered_sample_info_AROMA_2P_GMR_catch25_pyspi14.feather")
  
# Load study metadata
UCLA_CNP_metadata <- pyarrow_feather$read_feather("~/data/UCLA_CNP/study_metadata/UCLA_CNP_sample_metadata.feather") %>%
  mutate(Study = "UCLA_CNP") %>%
  filter(Sample_ID %in% UCLA_CNP_subjects_to_keep$Sample_ID)
ABIDE_ASD_metadata <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_sample_metadata.feather") %>%
  mutate(Study = "ABIDE_ASD")

# Load stats data
pairwise_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_balanced_accuracy_all_folds.feather"))
# Compute mean + SD performance across all folds
pairwise_balanced_accuracy <-pairwise_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Pairwise_Feature_Set, Analysis_Type, group_var) %>%
  reframe(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
          Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T))

pairwise_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_empirical_p_values.feather")) %>%
  dplyr::select(-Balanced_Accuracy_Across_Folds) %>%
  left_join(., pairwise_balanced_accuracy)

################################################################################
# Balanced accuracy heatmap

pairwise_p_values %>%
  filter(Pairwise_Feature_Set == pairwise_feature_set,
         Analysis_Type == "Pairwise_SPI",
         p_value_HolmBonferroni < 0.05) %>%
  dplyr::rename("pyspi_name" = "group_var") %>%
  left_join(., SPI_info) %>%
  mutate(Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                      Comparison_Group == "Bipolar" ~ "BP",
                                      T ~ Comparison_Group),
         Balanced_Accuracy_Across_Folds = 100*Balanced_Accuracy_Across_Folds) %>%
  mutate(Figure_name = fct_reorder(Figure_name, Balanced_Accuracy_Across_Folds, .fun=mean),
         Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BP", "ADHD", "ASD"))) %>%
  ggplot(data=., mapping=aes(x=Comparison_Group, y=Figure_name, 
                             fill=Balanced_Accuracy_Across_Folds)) +
  geom_tile()+
  geom_text(aes(label = round(Balanced_Accuracy_Across_Folds, 1))) +
  scale_fill_gradientn(colors=c(alpha("#AC77BD", 0.3), "#AC77BD"), 
                       na.value=NA)  + 
  scale_y_discrete(labels = wrap_format(28)) +
  labs(fill = "Mean Balanced Accuracy (%)") +
  xlab("Clinical Group") +
  ylab("Pairwise SPI") +
  theme(legend.position="none")
ggsave(glue("{plot_path}/SPI_wise_results.svg"),
       width=4.25, height=4.5, units="in", dpi=300)


################################################################################
# Correlogram by SPI for 10 SPIs with significant balanced accuracies
if (!file.exists(glue("{data_path}/SPI_performance_correlation_across_disorders.feather"))) {
  UCLA_CNP_pyspi14 <- pyarrow_feather$read_feather("~/data/UCLA_CNP/processed_data/UCLA_CNP_AROMA_2P_GMR_pyspi14_filtered.feather")  %>%
    left_join(., UCLA_CNP_metadata) %>%
    filter(!is.na(Diagnosis))
  ABIDE_ASD_pyspi14 <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/processed_data/ABIDE_ASD_FC1000_pyspi14_filtered.feather")  %>%
    left_join(., ABIDE_ASD_metadata) %>%
    filter(!is.na(Diagnosis))
  merged_pyspi14 <- plyr::rbind.fill(UCLA_CNP_pyspi14, ABIDE_ASD_pyspi14) %>%
    mutate(unique_ID = paste0(Sample_ID, "__", brain_region_from, "__", brain_region_to), .keep="unused") 
  
  rm(UCLA_CNP_pyspi14)
  rm(ABIDE_ASD_pyspi14)
  
  # Create adjacency matrix of feature values across all participants
  unique_groups <- unique(merged_pyspi14$SPI)
  group_combinations <- combn(unique_groups, 2, simplify = FALSE)
  
  cor_list <- list()
  
  for (i in 1:length(group_combinations)) {
    group_pair <- sort(group_combinations[[i]])
    
    SPI1 <- group_pair[1]
    SPI2 <- group_pair[2]
    
    SPI1_data <- merged_pyspi14 %>% filter(SPI==SPI1) %>% pull(value)
    SPI2_data <- merged_pyspi14 %>% filter(SPI==SPI2) %>% pull(value)
    
    correlation_df <- data.frame(SPI1 = SPI1,
                                 SPI2 = SPI2,
                                 SPI1_nickname = SPI1_nickname,
                                 SPI2_nickname = SPI2_nickname,
                                 spearman_corr_abs = abs(cor(SPI1_data, SPI2_data, method = "spearman",  use="complete.obs")))
    
    cor_list <- list.append(cor_list, correlation_df)
    
    rm(SPI1_data)
    rm(SPI2_data)
    gc()
  }
  
  data_for_corr_long <- do.call(plyr::rbind.fill, cor_list)
  
  # Write to feather file
  pyarrow_feather$write_feather(data_for_corr_long, glue("{data_path}/SPI_performance_correlation_across_disorders.feather"))
} else {
  data_for_corr_long <- pyarrow_feather$read_feather(glue("{data_path}/SPI_performance_correlation_across_disorders.feather"))
}

sig_SPIs <- pairwise_p_values %>%
  filter(Pairwise_Feature_Set == pairwise_feature_set,
         Analysis_Type == "Pairwise_SPI",
         p_value_HolmBonferroni < 0.05) %>%
  distinct(group_var) %>%
  dplyr::rename("pyspi_name" = "group_var") %>%
  left_join(., SPI_info) %>%
  pull(Figure_name)

# Convert to long for easy querying
data_for_corr_heatmap <- data_for_corr_long %>%
  plyr::rbind.fill(., data_for_corr_long %>% rename("SPI1_nickname"="SPI2_nickname", "SPI2_nickname"="SPI1_nickname")) %>%
  pivot_wider(id_cols=SPI1_nickname, names_from=SPI2_nickname, values_from=spearman_corr_abs) %>%
  column_to_rownames(var="SPI1_nickname") %>%
  select(all_of(sig_SPIs)) %>%
  dplyr::filter(row.names(.) %in% sig_SPIs) %>%
  as.matrix()
data_for_corr_heatmap[is.na(data_for_corr_heatmap)] <- 1

# Rearrange rows and columns
data_for_corr_heatmap <- data_for_corr_heatmap[sig_SPIs, sig_SPIs]

num_branches <- 6

ht1 <- ComplexHeatmap::Heatmap(data_for_corr_heatmap,
                               clustering_distance_rows = "spearman",
                               clustering_distance_columns = "spearman",
                               clustering_method_rows = "average",
                               clustering_method_columns = "average",
                               row_names_side = "right",
                               row_dend_side = "left", 
                               row_dend_width = unit(2, "cm"),
                               row_dend_gp = gpar(lwd=unit(2, "cm")),
                               row_split = num_branches,
                               row_gap = unit(3, "mm"),
                               column_split = num_branches,
                               column_gap = unit(3, "mm"),
                               row_title = NULL,
                               column_title = NULL,
                               show_row_names = TRUE,
                               row_names_gp = grid::gpar(fontsize = 26),
                               show_column_names = FALSE,
                               show_column_dend = FALSE,
                               border = TRUE,
                               col = c("white", RColorBrewer::brewer.pal(12, "Reds")),
                               name = "Spearman rank correlation, \u03c1",
                               heatmap_legend_param = list(legend_direction = "horizontal",
                                                           legend_width = unit(5, "cm")))

svg(glue("{plot_path}/pyspi14_feature_spearman_corr.svg"),
    width=8.25, height=6, bg=NA)
draw(ht1, heatmap_legend_side = "bottom",
     background = "transparent")
dev.off()



################################################################################
# How are different disorders correlated in terms of classification accuracy by pyspi14 SPI?

pairwise_p_values %>%
  left_join(., study_group_df) %>%
  dplyr::select(group_var, Group_Nickname, Balanced_Accuracy_Across_Folds) %>%
  pivot_wider(id_cols=group_var, names_from=Group_Nickname, values_from=Balanced_Accuracy_Across_Folds) %>%
  dplyr::select(-group_var) %>%
  cor(method="spearman") %>%
  as.data.frame() %>%
  rownames_to_column(var="group1") %>%
  pivot_longer(cols=c(-group1), names_to="group2", values_to="Spearman_Corr") %>%
  mutate(group1 = factor(group1, levels = c("SCZ", "BP", "ADHD", "ASD")),
         group2 = factor(group2, levels = rev(c("SCZ", "BP", "ADHD", "ASD")))) %>%
  ggplot(data=., mapping=aes(x=group1, y=group2, fill=Spearman_Corr)) + 
  geom_tile() +
  geom_text(aes(label=round(Spearman_Corr,2))) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  scale_fill_distiller(palette="RdBu", limits=c(-1,1)) +
  labs(fill="Spearman rank correlation, \u03c1") +
  theme(legend.position = "bottom",
        axis.title = element_blank())  +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5, barwidth=10))
ggsave(glue("{plot_path}/pyspi14_feature_performance_correlations.svg"),
       width=3, height=3.5, units="in", dpi=300)


################################################################################
# How are different disorders correlated in terms of classification accuracy by catch25 feature?

SPI_accuracy_corr_data <- pairwise_p_values %>%
  left_join(., study_group_df) %>%
  dplyr::select(group_var, Group_Nickname, Balanced_Accuracy_Across_Folds) %>%
  pivot_wider(id_cols=group_var, names_from=Group_Nickname, values_from=Balanced_Accuracy_Across_Folds) %>%
  dplyr::select(-group_var) %>%
  cor(method="spearman")

num_branches <- 2
ht2 <- ComplexHeatmap::Heatmap(SPI_accuracy_corr_data,
                               clustering_distance_rows = "spearman",
                               clustering_distance_columns = "spearman",
                               clustering_method_rows = "average",
                               clustering_method_columns = "average",
                               row_names_side = "right",
                               row_dend_side = "left", 
                               row_dend_width = unit(2, "cm"),
                               row_dend_gp = gpar(lwd=unit(2, "cm")),
                               row_split = num_branches,
                               row_gap = unit(3, "mm"),
                               column_split = num_branches,
                               column_gap = unit(3, "mm"),
                               row_title = NULL,
                               column_title = NULL,
                               show_row_names = TRUE,
                               show_column_names = FALSE,
                               show_column_dend = FALSE,
                               border = TRUE,
                               cell_fun = function(j, i, x, y, width, height, fill) {
                                 grid.text(sprintf("%.2f", SPI_accuracy_corr_data[i, j]), x, y, gp = gpar(fontsize = 16))
                               },
                               col = colorRamp2(seq(-1, 1, length.out=9), brewer_pal("div", "RdBu", direction=-1)(9)),
                               name = "Spearman rank correlation, \u03c1",
                               heatmap_legend_param = list(legend_direction = "horizontal",
                                                           legend_width = unit(5, "cm")))

svg(glue("{plot_path}/pyspi14_SPI_performance_correlations.svg"),
    width=5, height=4, bg=NA)
draw(ht2, heatmap_legend_side = "bottom",
     background = "transparent")
dev.off()

################################################################################
# How do directed SPI values differ to vs. from each brain region across cohorts?
if (!(file.exists(glue("{data_path}/pyspi14_directed_to_vs_from_regions_by_group.feather")))) {
  UCLA_CNP_pyspi14 <- pyarrow_feather$read_feather("~/data/UCLA_CNP/processed_data/UCLA_CNP_AROMA_2P_GMR_pyspi14_filtered.feather")  %>%
    left_join(., UCLA_CNP_metadata) %>%
    filter(!is.na(Diagnosis))
  ABIDE_ASD_pyspi14 <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/processed_data/ABIDE_ASD_FC1000_pyspi14_filtered.feather")  %>%
    left_join(., ABIDE_ASD_metadata) %>%
    filter(!is.na(Diagnosis))
  merged_pyspi14 <- plyr::rbind.fill(UCLA_CNP_pyspi14, ABIDE_ASD_pyspi14) 
  
  rm(UCLA_CNP_pyspi14)
  rm(ABIDE_ASD_pyspi14)
  
  # Iterate over each DIRECTED SPI
  for (directed_SPI in SPI_info %>% filter(Directed=="Yes") %>% pull(pyspi_name)) {
    directed_SPI_data <- subset(merged_pyspi14, SPI==directed_SPI)
    
    # For each brain region, compute the 
  }
  
  cor_list <- list()
}

pyspi14_lm_beta_coefficients %>%
  separate(Region_Pair, into=c("From", "To"), sep="_") %>%
  dplyr::select(Comparison_Group, SPI, From, To, estimate) %>%
  pivot_longer(cols=c(From, To), names_to="Direction", values_to="Brain_Region") %>%
  group_by(Comparison_Group, SPI, Brain_Region, Direction) %>%
  summarise(mean_estimate = mean(abs(estimate)),
            sd_estimate = sd(abs(estimate))) %>%
  left_join(., SPI_info, by=c("SPI"="pyspi_name")) %>%
  filter(Directed=="Yes") %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BP", "ADHD", "ASD"))) %>%
  ggplot(data=., mapping=aes(x=Direction, y=mean_estimate, group=Brain_Region, color=Comparison_Group)) +
  geom_line(alpha=0.8) +
  facet_grid(Figure_name ~ Comparison_Group, switch="y", space="free") +
  scale_color_manual(values = c("SCZ"="#573DC7", 
                                "BP"="#D5492A", 
                                "ADHD"="#0F9EA9",
                                "ASD"="#C47B2F")) +
  ylab("Mean |\u03b2| per Region") +
  scale_x_discrete(expand=c(0.1,0.1)) +
  theme(legend.position = "none",
        axis.text.y = element_text(size=10),
        strip.text = element_text(angle=0, face="bold"),
        strip.text.y.left = element_text(angle=0, face="bold"),
        strip.placement = "outside",
        strip.background = element_blank())
ggsave(glue("{plot_path}/pyspi14_SPI_beta_coefficients_directed_SPIs.svg"),
       width=7, height=6, units="in", dpi=300)

################################################################################
# How do lm beta coefficients differ TO vs. FROM each region for the directed SPIs?

pyspi14_lm_beta_coefficients %>%
  separate(Region_Pair, into=c("From", "To"), sep="_") %>%
  dplyr::select(Comparison_Group, SPI, From, To, estimate) %>%
  pivot_longer(cols=c(From, To), names_to="Direction", values_to="Brain_Region") %>%
  group_by(Comparison_Group, SPI, Brain_Region, Direction) %>%
  summarise(mean_estimate = mean(abs(estimate)),
            sd_estimate = sd(abs(estimate))) %>%
  left_join(., SPI_info, by=c("SPI"="pyspi_name")) %>%
  filter(Directed=="Yes") %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BP", "ADHD", "ASD"))) %>%
  ggplot(data=., mapping=aes(x=Direction, y=mean_estimate, group=Brain_Region, color=Comparison_Group)) +
  geom_line(alpha=0.8) +
  facet_grid(Figure_name ~ Comparison_Group, switch="y", space="free") +
  scale_color_manual(values = c("SCZ"="#573DC7", 
                                "BP"="#D5492A", 
                                "ADHD"="#0F9EA9",
                                "ASD"="#C47B2F")) +
  ylab("Mean |\u03b2| per Region") +
  scale_x_discrete(expand=c(0.1,0.1)) +
  theme(legend.position = "none",
        axis.text.y = element_text(size=10),
        strip.text = element_text(angle=0, face="bold"),
        strip.text.y.left = element_text(angle=0, face="bold"),
        strip.placement = "outside",
        strip.background = element_blank())
ggsave(glue("{plot_path}/pyspi14_SPI_beta_coefficients_directed_SPIs.svg"),
       width=7, height=6, units="in", dpi=300)
