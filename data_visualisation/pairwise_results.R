################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/pairwise_results/")
TAF::mkdir(plot_path)

# python_to_use <- "~/.conda/envs/pyspi/bin/python3"
python_to_use <- "/Users/abry4213/anaconda3/envs/pyspi/bin/python3"
pairwise_feature_set <- "pyspi14"
univariate_feature_set <- "catch24"
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
library(ggseg)
library(ggsegHO)
library(knitr)
library(kableExtra)
library(patchwork)
library(ggseg)
library(broom)
library(colorspace)
library(see)
library(ggridges)
library(ggnewscale)
library(scales)
library(splitstackshape)
library(LaCroixColoR)
library(correctR)
library(ggpubr)
library(ggsignif)
library(poolr)
library(gtools)
theme_set(theme_cowplot())

# Source visualisation script
source(glue("{github_dir}/data_visualisation/Manuscript_Draft_Visualisations_Helper.R"))

# Load in SPI info
SPI_info <- read.csv(glue("{github_dir}/data_visualisation/SPI_info.csv"))

# Load brain region info
UCLA_CNP_brain_region_info <- read.csv("~/data/UCLA_CNP/study_metadata/UCLA_CNP_Brain_Region_info.csv")
ABIDE_ASD_brain_region_info <- read.table("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_Harvard_Oxford_cort_prob_2mm_ROI_lookup.txt", sep=";", header = T) %>%
  mutate(Brain_Region = ifelse(Index==45, "Heschl's Gyrus (includes H1 and H2)", Brain_Region))
region_node_to_from <- read.csv("~/data/TS_feature_manuscript/node_to_from_structure.csv")

# Map to the coordinates in MNI152 space
aparc_aseg_coords <- read.csv("/Users/abry4213/data/TS_feature_manuscript/aparc_aseg_MNI_coords.csv", header=F) %>%
  mutate(Index = row_number())
HO_coords <- read.csv("/Users/abry4213/data/TS_feature_manuscript/HO_MNI_coords.csv", header=F) %>%
  mutate(Index = row_number())

# Load participants included
UCLA_CNP_subjects_to_keep <- pyarrow_feather$read_feather("~/data/UCLA_CNP/processed_data/UCLA_CNP_filtered_sample_info_AROMA_2P_GMR_catch24_pyspi14.feather")
  
# Load study metadata
UCLA_CNP_metadata <- pyarrow_feather$read_feather("~/data/UCLA_CNP/study_metadata/UCLA_CNP_sample_metadata.feather") %>%
  mutate(Study = "UCLA_CNP") %>%
  filter(Sample_ID %in% UCLA_CNP_subjects_to_keep$Sample_ID)
ABIDE_ASD_metadata <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_sample_metadata.feather") %>%
  mutate(Study = "ABIDE_ASD")

# Load stats data
univariate_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_mixedsigmoid_scaler_empirical_p_values.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)

pairwise_balanced_accuracy_AUC_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_mixedsigmoid_scaler_balanced_accuracy_AUC_all_folds.feather"))
pairwise_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_mixedsigmoid_scaler_empirical_p_values.feather"))
pairwise_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_mixedsigmoid_scaler_null_balanced_accuracy_distributions.feather"))
univariate_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_mixedsigmoid_scaler_null_balanced_accuracy_distributions.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set)

combo_univariate_pairwise_balanced_accuracy_AUC_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combined_univariate_pairwise_mixedsigmoid_scaler_balanced_accuracy_AUC_all_folds.feather")) %>%
  mutate(Analysis_Type = "SPI_Univariate_Combo")
combo_univariate_pairwise_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combined_univariate_pairwise_mixedsigmoid_scaler_empirical_p_values.feather"))
combo_univariate_pairwise_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combined_univariate_pairwise_mixedsigmoid_scaler_null_balanced_accuracy_distributions.feather"))

# Load TPR/FPR data
pairwise_TPR_FPR <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_mixedsigmoid_scaler_ROC_TPR_FPR.feather")) %>%
  filter(Pairwise_Feature_Set == pairwise_feature_set) 

# Load lm beta statistics
lm_beta_stats_pyspi14_whole_brain <- pyarrow_feather$read_feather(glue("{data_path}/pairwise_pyspi14_lm_beta_statistics_by_region_pair.feather"))

################################################################################
# SPI-wise SVM results
################################################################################

# Plot ROC of top-performing features
top_SPIs_to_find_AUC <- pairwise_p_values %>%
  filter(Analysis_Type == "Pairwise_SPI") %>%
  group_by(Comparison_Group, Study) %>%
  slice_max(n=1, order_by=Balanced_Accuracy_Across_Repeats) %>%
  mutate(Group_Nickname =  case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                     Comparison_Group == "Bipolar" ~ "BPD",
                                     T ~ Comparison_Group)) %>%
  dplyr::select(Study, Comparison_Group, Group_Nickname, group_var, ROC_AUC_Across_Repeats)

pairwise_TPR_FPR  %>%
  mutate(Group_Nickname =  case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                     Comparison_Group == "Bipolar" ~ "BPD",
                                     T ~ Comparison_Group)) %>%
  semi_join(top_SPIs_to_find_AUC)  %>%
  ggplot(data=.) +
  geom_abline(slope=1, linetype=2) +
  geom_smooth(se=T, aes(color=Comparison_Group, x=fpr,y=tpr)) +
  xlab("FPR") +
  ylab("TPR") +
  coord_equal() +
  geom_text(data = top_SPIs_to_find_AUC,
            aes(label=glue("{Group_Nickname}: {round(ROC_AUC_Across_Repeats, 2)}"),
                color=Comparison_Group),
            x = 1, y=c(0.1, 0.2, 0.3, 0.4), 
            size=4.5, hjust=1) +
  theme(legend.position = "none") +
  scale_color_manual(values=c("Control" = "#5BB67B", 
                              "Schizophrenia" = "#573DC7", 
                              "Bipolar" = "#D5492A", 
                              "ADHD" = "#0F9EA9", 
                              "ASD" = "#C47B2F"))
ggsave(glue("{plot_path}/pairwise_top_SPI_ROC.svg"),
       width=3, height=3, units="in", dpi=300)

# Annotation bar with SPI type
pairwise_p_values %>%
  filter(Pairwise_Feature_Set == pairwise_feature_set,
         Analysis_Type == "Pairwise_SPI",
         p_value_Bonferroni < 0.05) %>%
  dplyr::rename("SPI" = "group_var") %>%
  left_join(., SPI_info) %>%
  mutate(Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                      Comparison_Group == "Bipolar" ~ "BPD",
                                      T ~ Comparison_Group),
         Balanced_Accuracy_Across_Repeats = 100*Balanced_Accuracy_Across_Repeats) %>%
  mutate(Nickname = fct_reorder(Nickname, Balanced_Accuracy_Across_Repeats, .fun=mean),
         Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD"))) %>%
  group_by(SPI, Category) %>%
  summarise(Balacc_Sum = sum(Balanced_Accuracy_Across_Repeats)) %>%
  ungroup() %>%
  mutate(SPI = fct_reorder(SPI, Balacc_Sum),
         Category = fct_reorder(Category, Balacc_Sum, .fun=sum, .desc=T)) %>%
  ggplot(data=., mapping=aes(x=0, y=SPI, fill=Category)) +
  geom_tile() +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text=element_text(size=14)) +
  guides(fill = guide_legend(title.position = "top", 
                             ncol = 2,
                             byrow=T,
                             title.hjust = 0.5)) 
ggsave(glue("{plot_path}/SPI_wise_colorbar.svg"),
       width=6, height=6, units="in", dpi=300)

# Actual heatmap
pairwise_p_values %>%
  filter(Pairwise_Feature_Set == pairwise_feature_set,
         Analysis_Type == "Pairwise_SPI",
         p_value_Bonferroni < 0.05) %>%
  dplyr::rename("SPI" = "group_var") %>%
  left_join(., SPI_info) %>%
  mutate(Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                      Comparison_Group == "Bipolar" ~ "BPD",
                                      T ~ Comparison_Group),
         Balanced_Accuracy_Across_Repeats = 100*Balanced_Accuracy_Across_Repeats) %>%
  mutate(Nickname = fct_reorder(Nickname, Balanced_Accuracy_Across_Repeats, .fun=mean),
         Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD"))) %>%
  ggplot(data=., mapping=aes(x=Comparison_Group, y=Nickname, 
                             fill=Balanced_Accuracy_Across_Repeats)) +
  geom_tile()+
  geom_text(aes(label = round(Balanced_Accuracy_Across_Repeats, 1))) +
  scale_fill_gradientn(colors=c(alpha("#AC77BD", 0.3), "#AC77BD"), 
                       na.value=NA)  + 
  scale_y_discrete(labels = wrap_format(28)) +
  labs(fill = "Mean Balanced Accuracy (%)") +
  xlab("Clinical Group") +
  ylab("Pairwise SPI") +
  theme(legend.position="none")
ggsave(glue("{plot_path}/SPI_wise_results.svg"),
       width=5, height=4.5, units="in", dpi=300)

################################################################################
# Condense lm beta statistics down to the regional level
################################################################################

# Demo brain figure
dk %>%
  as_tibble() %>%
  mutate(fillval = case_when(label == "lh_caudalmiddlefrontal" ~ "1",
                             label == "lh_bankssts" ~ "2",
                             label == "lh_lateralorbitofrontal" ~ "3", 
                             label == "lh_superiorparietal" ~ "4",
                             label == "lh_lateraloccipital" ~ "5",
                             T ~ NA_character_)) %>%
  ggseg(atlas = "dk", mapping = aes(fill = fillval),
        hemisphere="left",
        view = "lateral",
        position = "stacked", colour = "gray50") +
  scale_fill_manual(values=lacroix_palette("PassionFruit", n=5), na.value="white") +
  theme_void() +
  theme(plot.title = element_blank(),
        legend.position = "none")
ggsave(glue("{plot_path}/demo_brain_for_FC.svg"), width = 3, height=2, units="in", dpi=300)

# Find regions most disrupted across all pairwise connections
pairwise_lm_beta_stats_by_region_from <- lm_beta_stats_pyspi14_whole_brain %>%
  separate(Region_Pair, c("region_from", "region_to"),
           sep="_") %>%
  group_by(region_from, SPI, Study, Comparison_Group) %>%
  summarise(mean_lm_beta_magnitude = mean(abs(estimate))) %>%
  dplyr::rename("Brain_Region" = "region_from") %>%
  mutate(Direction = "from")

pairwise_lm_beta_stats_by_region_to <- lm_beta_stats_pyspi14_whole_brain %>%
  separate(Region_Pair, c("region_from", "region_to"),
           sep="_") %>%
  group_by(region_to, SPI, Study, Comparison_Group) %>%
  summarise(mean_lm_beta_magnitude = mean(abs(estimate)))%>%
  dplyr::rename("Brain_Region" = "region_to") %>%
  mutate(Direction = "to")

# Helper function to plot the beta coefficients for a given feature in the brain
plot_feature_in_brain <- function(study_group_df, lm_beta_df, SPI_name, min_fill,
                                  max_fill, bin_seq, fill_colors) {
  
  ggseg_plot_list <- list()
  
  for (i in 1:nrow(study_group_df)) {
    dataset_ID <- study_group_df$Study[i]
    comparison_group <- study_group_df$Group_Nickname[i]
    
    # Define atlas by study
    atlas <- ifelse(dataset_ID == "UCLA_CNP", "dk", "hoCort")
    
    if (dataset_ID == "ABIDE_ASD") {
      lm_beta_stat_data <- lm_beta_df %>%
        filter(Study == dataset_ID, 
               Comparison_Group == comparison_group) %>%
        left_join(., ABIDE_ASD_brain_region_info) %>%
        distinct() 
    } else {
      lm_beta_stat_data <- lm_beta_df %>%
        filter(Study == dataset_ID, 
               Comparison_Group == comparison_group) %>%
        mutate(label = ifelse(str_detect(Brain_Region, "ctx-"),
                              gsub("-", "_", Brain_Region),
                              as.character(Brain_Region))) %>%
        mutate(label = gsub("ctx_", "", label)) %>%
        distinct()
    }
    
    # Plot T stat data in cortex
    dataset_ggseg <- plot_data_with_ggseg_discrete(dataset_ID = dataset_ID,
                                                   atlas_name=atlas,
                                                   atlas_data=get(atlas) %>% as_tibble(),
                                                   data_to_plot = lm_beta_stat_data,
                                                   fill_variable = "mean_lm_beta_magnitude",
                                                   fill_colors = fill_colors,
                                                   bin_seq = bin_seq,
                                                   line_color = "gray30",
                                                   na_color = "white")  +
      labs(fill="Beta")
    
    ggseg_plot_list <- list.append(ggseg_plot_list, dataset_ggseg)
    
    # Add subcortical data for UCLA CNP
    if (dataset_ID == "UCLA_CNP") {
      dataset_ggseg_subctx <- plot_data_with_ggseg_discrete(dataset_ID = dataset_ID,
                                                            atlas_name = "aseg",
                                                            atlas_data = aseg %>% as_tibble(),
                                                            data_to_plot=lm_beta_stat_data,
                                                            fill_variable = "mean_lm_beta_magnitude",
                                                            fill_colors = fill_colors,
                                                            bin_seq = bin_seq,
                                                            line_color = "gray30",
                                                            na_color = "white")  +
        labs(fill="Beta")
      
      # Append to list
      ggseg_plot_list <- list.append(ggseg_plot_list, dataset_ggseg_subctx)
    }
  }
  
  return(ggseg_plot_list)
}

# Plot lm beta statistics for average Pearson correlation in the brain
pearson_fc_lm_beta_for_ggseg <- pairwise_lm_beta_stats_by_region_from %>%
  filter(SPI == "cov_EmpiricalCovariance") %>%
  dplyr::select(Study, Comparison_Group, SPI, Brain_Region, mean_lm_beta_magnitude)

min_fill <- 0
max_fill <- round(max(pearson_fc_lm_beta_for_ggseg$mean_lm_beta_magnitude), 1)
bin_seq <- seq(min_fill, max_fill, length.out=7)
fill_colors <- c("white", "white", colorRampPalette(c("white", "#BB8EC8", "#884c9c"))(4))

pearson_fc_lm_in_brain <- plot_feature_in_brain(study_group_df = study_group_df,
                                                      lm_beta_df = pearson_fc_lm_beta_for_ggseg,
                                                      SPI_name = "cov_EmpiricalCovariance",
                                                      min_fill = min_fill,
                                                      max_fill = max_fill,
                                                      bin_seq = bin_seq,
                                                      fill_colors = fill_colors)

wrap_plots(pearson_fc_lm_in_brain, 
           ncol=2, 
           byrow=T) + 
  plot_layout(guides = "collect")
ggsave(glue("{plot_path}/Region_wise_avg_lm_beta_stat_pearson_corrs.svg"),
       width=5, height=7, units="in", dpi=300)

################################################################################
# Plot Pearson correlation mean value vs. univariate region balanced accuracy
pearson_balacc_data <- pearson_fc_lm_beta_for_ggseg %>% 
  dplyr::select(-SPI) %>%
  left_join(., univariate_p_values %>% 
              filter(Analysis_Type=="Univariate_Brain_Region")
            %>% dplyr::rename("Brain_Region" = "group_var") %>%
              mutate(Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                                  Comparison_Group == "Bipolar" ~ "BPD",
                                                  T ~ Comparison_Group))) %>%
  mutate(univar_sig = ifelse(p_value_Bonferroni < 0.05, "Significant", "Not Significant"),
         Comparison_Group = factor(Comparison_Group, levels=c("SCZ", "BPD", "ADHD", "ASD"))) 

pearson_balacc_data %>%
  ggplot(data=., mapping=aes(x=mean_lm_beta_magnitude, y=Balanced_Accuracy_Across_Repeats)) +
  geom_point(aes(color = univar_sig)) +
  scale_color_manual(values=c("Significant"="#884c9c", "Not Significant" = "gray70")) +
  facet_wrap(Comparison_Group ~ ., scales="free", nrow=1) +
  theme(legend.position="bottom",
        axis.title = element_text(size=14),
        axis.text = element_text(size=10)) +
  ylab("Univariate Mean\nBalanced Accuracy") +
  xlab("Mean Pearson FC |Beta| from OLS Model") +
  labs(color = "Univariate Classification") +
  guides(color = guide_legend(nrow=1)) + 
  scale_x_continuous(breaks = seq(min(floor(pearson_balacc_data$mean_lm_beta_magnitude / 0.05) * 0.05), 
                                  max(ceiling(pearson_balacc_data$mean_lm_beta_magnitude / 0.05) * 0.05), 
                                  by = 0.05))
ggsave(glue("{plot_path}/Region_wise_avg_Pearson_vs_BalAcc.svg"),
       width=8, height=2.75, units="in", dpi=300)

################################################################################
# Hypo vs hyper connectivity
################################################################################

# Prep data for charticulate
# https://charticulator.com/gallery/les_miserables_circular.html

# attributes: UCLA CNP
UCLA_CNP_attributes <- region_node_to_from %>%
  filter(Study=="UCLA_CNP") %>%
  filter(from != "Origin") %>%
  dplyr::select("to", "from") %>%
  dplyr::rename("id" = "to",
                "cortex" = "from") %>%
  arrange(cortex)

# ggseg vis for lobes
UCLA_CNP_attributes %>%
  mutate(label = ifelse(str_detect(id, "ctx-"),
                        gsub("-", "_", id),
                        as.character(id))) %>%
  mutate(label = gsub("ctx_", "", label)) %>%
  left_join(., dk %>% as_tibble()) %>%
  dplyr::select(-label) %>%
  ggseg(atlas = dk, mapping = aes(fill = cortex),
        position = "stacked", colour = "gray20", hemisphere="left") +
  scale_fill_manual(values=c("Cingulate" = "#AECDE1",
                             "Frontal" = "#3C76AF",
                             "Insula" = "#BBDE93",
                             "Occipital" = "#549E3F",
                             "Parietal" = "#ED9F9C",
                             "Temporal" = "#F4C17B",
                             "Subcortex" = "#D0352B"),
                    na.value = "white") +
  theme_void() +
  theme(legend.position = "none")
ggsave(glue("{plot_path}/cortical_lobes.svg"),
       width=4, height=2, units="in", dpi=300)
# Subcortical lobes
UCLA_CNP_attributes %>%
  mutate(label = ifelse(str_detect(id, "ctx-"),
                        gsub("-", "_", id),
                        as.character(id))) %>%
  mutate(label = gsub("ctx_", "", label)) %>%
  left_join(., aseg %>% as_tibble()) %>%
  dplyr::select(-label) %>%
  filter(!is.na(atlas)) %>%
  ggplot() +
  geom_brain(atlas = aseg, mapping = aes(fill=cortex), 
             side = "coronal", colour = "gray20")  +
  scale_fill_manual(values=c("Cingulate" = "#AECDE1",
                             "Frontal" = "#3C76AF",
                             "Insula" = "#BBDE93",
                             "Occipital" = "#549E3F",
                             "Parietal" = "#ED9F9C",
                             "Temporal" = "#F4C17B",
                             "Subcortex" = "#D0352B"),
                    na.value = "white") +
  theme_void() +
  theme(legend.position = "none")
ggsave(glue("{plot_path}/subcortical_lobes.svg"),
       width=2, height=1.5, units="in", dpi=300)

# Write all pairs to a pandas dataframe
write_all_connected_pairs <- function(study, 
                                      lm_beta_stats,
                                      comparison_group) {
  all_pairs <- lm_beta_stats %>%
    filter(Study==study, Comparison_Group==comparison_group) %>%
    group_by(Region_Pair) %>%
    summarise(mean_beta = mean(estimate)) %>%
    separate(Region_Pair, into=c("from", "to"), sep="_") %>%
    dplyr::rename("source_id" = "from",
                  "target_id" = "to",
                  "value" = "mean_beta") %>%
    arrange(value)
  
  write.table(all_pairs, glue("{data_path}/{study}_{comparison_group}_all_FC_pairs.csv"), 
              row.names=F, col.names = T, sep=",")
}

# Calculate adjacency matrix for UCLA CNP SCZ
write_top_hyper_hypo_matrices <- function(study,
                                          comparison_group,
                                          lm_beta_stats,
                                          brain_region_lookup,
                                          MNI_coords_df,
                                          node_color_lookup) {
  hyper_adj_matrix <- lm_beta_stats %>%
    filter(Study==study, Comparison_Group==comparison_group) %>%
    group_by(Region_Pair) %>%
    summarise(mean_beta = mean(estimate)) %>%
    separate(Region_Pair, into=c("from", "to"), sep="_") %>%
    dplyr::rename("source_id" = "from",
                  "target_id" = "to",
                  "value" = "mean_beta") %>%
    slice_max(order_by=value, n=5) %>%
    left_join(., brain_region_lookup, by=c("source_id" = "Brain_Region")) %>%
    dplyr::rename("Index_from" = "Index") %>%
    left_join(., brain_region_lookup, by=c("target_id" = "Brain_Region")) %>%
    dplyr::rename("Index_to" = "Index") %>%
    dplyr::select(Index_from, Index_to, value) %>%
    mutate(Index_from = as.character(Index_from), Index_to = as.character(Index_to)) %>%
    arrange(Index_from, Index_to) %>%
    pivot_wider(id_cols="Index_from", names_from=Index_to, values_from=value) %>%
    column_to_rownames(var="Index_from") %>%
    as.matrix()

  hypo_adj_matrix <- lm_beta_stats %>%
    filter(Study==study, Comparison_Group==comparison_group) %>%
    group_by(Region_Pair) %>%
    summarise(mean_beta = mean(estimate)) %>%
    separate(Region_Pair, into=c("from", "to"), sep="_") %>%
    dplyr::rename("source_id" = "from",
                  "target_id" = "to",
                  "value" = "mean_beta") %>%
    slice_min(order_by=value, n=5) %>%
    left_join(., brain_region_lookup, by=c("source_id" = "Brain_Region")) %>%
    dplyr::rename("Index_from" = "Index") %>%
    left_join(., brain_region_lookup, by=c("target_id" = "Brain_Region")) %>%
    dplyr::rename("Index_to" = "Index") %>%
    dplyr::select(Index_from, Index_to, value) %>%
    mutate(Index_from = as.character(Index_from), Index_to = as.character(Index_to)) %>%
    arrange(Index_from, Index_to) %>%
    pivot_wider(id_cols="Index_from", names_from=Index_to, values_from=value) %>%
    column_to_rownames(var="Index_from") %>%
    as.matrix()
  
  # Create a mapping of row and column labels
  labels_hyper <- unique(c(rownames(hyper_adj_matrix), colnames(hyper_adj_matrix)))
  labels_hypo <- unique(c(rownames(hypo_adj_matrix), colnames(hypo_adj_matrix)))
  
  # Create a 10x10 matrix with NA values
  adjacency_10x10_hyper <- matrix(NA, nrow = length(unique(labels_hyper)), ncol = length(unique(labels_hyper)))
  adjacency_10x10_hypo <-matrix(NA, nrow = length(unique(labels_hypo)), ncol = length(unique(labels_hypo)))
  
  # Print the resulting 10x10 matrix with row and column labels
  colnames(adjacency_10x10_hyper) <- rownames(adjacency_10x10_hyper) <- labels_hyper
  colnames(adjacency_10x10_hypo) <- rownames(adjacency_10x10_hypo) <- labels_hypo
  
  # Fill in with real values
  for (i in 1:nrow(permutations(n=length(labels_hyper), r=2, v=labels_hyper))) {
    index = permutations(n=length(labels_hyper), r=2, v=labels_hyper)[i,]
    if (index[1] %in% rownames(hyper_adj_matrix) & index[2] %in% colnames(hyper_adj_matrix)) {
      mat_value = hyper_adj_matrix[index[1], index[2]]
      adjacency_10x10_hyper[index[1], index[2]] <- mat_value
    }
  }
  for (i in 1:nrow(permutations(n=length(labels_hypo), r=2, v=labels_hypo))) {
    index = permutations(n=length(labels_hypo), r=2, v=labels_hypo)[i,]
    if (index[1] %in% rownames(hypo_adj_matrix) & index[2] %in% colnames(hypo_adj_matrix)) {
      mat_value = hypo_adj_matrix[index[1], index[2]]
      adjacency_10x10_hypo[index[1], index[2]] <- mat_value
    }
  }
  
  # Write to output tables
  write.table(adjacency_10x10_hyper, file=glue("/Users/abry4213/data/TS_feature_manuscript/{study}_{comparison_group}_hyper_FC.csv"), 
              sep=",", quote=F, na="NA", row.names=F, col.names=F)
  write.table(adjacency_10x10_hypo, file=glue("/Users/abry4213/data/TS_feature_manuscript/{study}_{comparison_group}_hypo_FC.csv"), 
              sep=",", quote=F, na="NA", row.names=F, col.names=F)
  
  # Filter coords and save to tables
  MNI_coords_df %>%
    filter(Index %in% labels_hyper) %>%
    mutate(Index = factor(Index, levels = labels_hyper)) %>%
    arrange(Index) %>%
    dplyr::select(-Index) %>%
    write.table(., file=glue("/Users/abry4213/data/TS_feature_manuscript/{study}_{comparison_group}_hyper_FC_coords.csv"), 
                sep=",", quote=F, na="NA", row.names=F, col.names=F)
  MNI_coords_df %>%
    filter(Index %in% labels_hypo) %>%
    mutate(Index = factor(Index, levels = labels_hypo)) %>%
    arrange(Index) %>%
    dplyr::select(-Index) %>%
    write.table(., file=glue("/Users/abry4213/data/TS_feature_manuscript/{study}_{comparison_group}_hypo_FC_coords.csv"), 
                sep=",", quote=F, na="NA", row.names=F, col.names=F)
  
  # Extract node colors
  hyper_nodes <- brain_region_lookup %>% 
    filter(Index %in% as.numeric(labels_hyper)) %>%
    mutate(Index = factor(Index, levels=labels_hyper)) %>%
    arrange(Index) %>%
    left_join(., node_color_lookup) %>%
    pull(color)
  hypo_nodes <- brain_region_lookup %>% 
    filter(Index %in% as.numeric(labels_hypo)) %>%
    mutate(Index = factor(Index, levels=labels_hypo)) %>%
    arrange(Index) %>%
    left_join(., node_color_lookup) %>%
    pull(color)
  data.frame(Hyper_Nodes = hyper_nodes, Hypo_Nodes = hypo_nodes) %>%
    write.table(., file=glue("/Users/abry4213/data/TS_feature_manuscript/{study}_{comparison_group}_node_colors.csv"), 
                sep=",", quote=F, na="NA", row.names=F, col.names=T)
}

node_color_lookup <- data.frame(Cortex=c("Cingulate", "Frontal", "Insula", "Occipital", "Parietal", "Temporal", "Subcortex"),
                                color = c("#AECDE1","#3C76AF","#BBDE93","#549E3F","#ED9F9C","#F4C17B","#D0352B"))
   

# SCZ
write_all_connected_pairs(study = "UCLA_CNP",
                          lm_beta_stats = lm_beta_stats_pyspi14_whole_brain,
                          comparison_group = "SCZ")
write_top_hyper_hypo_matrices(study = "UCLA_CNP",
                              lm_beta_stats = lm_beta_stats_pyspi14_whole_brain,
                              comparison_group = "SCZ",
                              brain_region_lookup = UCLA_CNP_brain_region_info,
                              MNI_coords_df = aparc_aseg_coords,
                              node_color_lookup=node_color_lookup)

# BPD
write_all_connected_pairs(study = "UCLA_CNP",
                          lm_beta_stats = lm_beta_stats_pyspi14_whole_brain,
                          comparison_group = "BPD")
write_top_hyper_hypo_matrices(study = "UCLA_CNP",
                              lm_beta_stats = lm_beta_stats_pyspi14_whole_brain,
                              comparison_group = "BPD",
                              brain_region_lookup = UCLA_CNP_brain_region_info,
                              MNI_coords_df = aparc_aseg_coords,
                              node_color_lookup=node_color_lookup)

# ADHD
write_all_connected_pairs(study = "UCLA_CNP",
                          lm_beta_stats = lm_beta_stats_pyspi14_whole_brain,
                          comparison_group = "ADHD")
write_top_hyper_hypo_matrices(study = "UCLA_CNP",
                              lm_beta_stats = lm_beta_stats_pyspi14_whole_brain,
                              comparison_group = "ADHD",
                              brain_region_lookup = UCLA_CNP_brain_region_info,
                              MNI_coords_df = aparc_aseg_coords,
                              node_color_lookup=node_color_lookup)


# ASD
write_all_connected_pairs(study = "ABIDE_ASD",
                          lm_beta_stats = lm_beta_stats_pyspi14_whole_brain,
                          comparison_group = "ASD")
write_top_hyper_hypo_matrices(study = "ABIDE_ASD",
                              lm_beta_stats = lm_beta_stats_pyspi14_whole_brain,
                              comparison_group = "ASD",
                              brain_region_lookup = ABIDE_ASD_brain_region_info,
                              MNI_coords_df = HO_coords,
                              node_color_lookup=node_color_lookup)

################################################################################
# Plot composite mean value per region vs. univariate region balanced accuracy
# Plot from vs to each region separately
lm_beta_stats_pyspi14_whole_brain %>%
  group_by(Study, Comparison_Group, Region_Pair) %>%
  summarise(composite_beta = mean(abs(estimate))) %>%
  separate(Region_Pair, into=c("from", "to"), sep="_") %>%
  pivot_longer(cols=c(from, to), names_to="Direction", values_to="Brain_Region") %>%
  dplyr::select(Study, Comparison_Group, composite_beta, Direction, Brain_Region) %>%
  group_by(Study, Comparison_Group, Direction, Brain_Region) %>%
  summarise(mean_composite_beta = mean(composite_beta)) %>%
  left_join(., univariate_p_values %>% 
              filter(Analysis_Type=="Univariate_Brain_Region")
            %>% dplyr::rename("Brain_Region" = "group_var") %>%
              mutate(Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                                  Comparison_Group == "Bipolar" ~ "BPD",
                                                  T ~ Comparison_Group))) %>%
  mutate(univar_sig = ifelse(p_value_Bonferroni < 0.05, "Significant", "Not Significant")) %>%
  ggplot(data=., mapping=aes(x=mean_composite_beta, y=Balanced_Accuracy_Across_Repeats)) +
  geom_point(aes(color=univar_sig)) +
  scale_color_manual(values=c("Significant"="#884c9c", "Not Significant" = "gray70")) +
  facet_wrap(Comparison_Group ~ Direction, dir="v", nrow=2, scales="free_y") +
  theme(legend.position="bottom",
        axis.title = element_text(size=14),
        axis.text = element_text(size=10)) +
  ylab("Univariate Mean Balanced Accuracy") +
  xlab("Mean Composite Beta from OLS Model") +
  labs(color = "Univariate Classification")
ggsave(glue("{plot_path}/Region_wise_composite_mean_beta_vs_BalAcc.svg"),
       width=7, height=4.5, units="in", dpi=300)


################################################################################
# Distribution of SPI-wise lm beta statistics
################################################################################
lm_beta_stats_pyspi14_whole_brain <- feather::read_feather(glue("{data_path}/pairwise_pyspi14_lm_beta_statistics_by_region_pair.feather"))

lm_beta_stats_pyspi14_whole_brain %>%
  ungroup() %>%
  left_join(., SPI_info) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD")))%>%
  mutate(Nickname = fct_reorder(Nickname, estimate, .fun=sd)) %>%
  ggplot(data=., mapping=aes(x=estimate, y=Nickname, fill=Comparison_Group, color=Comparison_Group)) +
  geom_density_ridges(alpha=0.6, scale=1.1, rel_min_height = 0.01) +
  scale_x_continuous(limits=c(-0.5,0.5)) +
  xlab("Beta coefficient across\nall brain regions") +
  ylab("pyspi14 time-series feature") +
  scale_fill_manual(values=c("Control" = "#5BB67B", 
                             "SCZ" = "#573DC7", 
                             "BPD" = "#D5492A", 
                             "ADHD" = "#0F9EA9", 
                             "ASD" = "#C47B2F")) +
  scale_color_manual(values=c("Control" = "#5BB67B", 
                              "SCZ" = "#573DC7", 
                              "BPD" = "#D5492A", 
                              "ADHD" = "#0F9EA9", 
                              "ASD" = "#C47B2F")) +
  guides(fill = guide_legend(nrow=2, byrow=T),
         color = guide_legend(nrow=2, byrow=T)) +
  scale_y_discrete(labels = wrap_format(20)) +
  theme(legend.position = "bottom",
        axis.title = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        legend.text = element_text(size=16),
        legend.title = element_blank())
ggsave(glue("{plot_path}/pyspi14_feature_lm_beta_statistics_across_brain.svg"),
       width=5, height=10, units="in", dpi=300)

