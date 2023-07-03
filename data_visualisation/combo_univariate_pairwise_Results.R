################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/combo_univariate_pairwise_results/")
TAF::mkdir(plot_path)

# python_to_use <- "~/.conda/envs/pyspi/bin/python3"
python_to_use <- "/Users/abry4213/anaconda3/envs/pyspi/bin/python3"
univariate_feature_set <- "catch24"
pairwise_feature_set <- "pyspi14"
data_path <- "~/data/TS_feature_manuscript"
UCLA_CNP_sample_metadata <- "~/data/UCLA_CNP/study_metadata/UCLA_CNP_sample_metadata.feather"
ABIDE_ASD_sample_metadata <- "~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_sample_metadata.feather"
study_group_df <- data.frame(Study = c(rep("UCLA_CNP", 3), "ABIDE_ASD"),
                             Num_Samples = c(166, 157, 167, 1150), 
                             Noise_Proc = c(rep("AROMA+2P+GMR",3), "FC1000"),
                             Comparison_Group = c("Schizophrenia", "ADHD", "Bipolar", "ASD"))
univariate_feature_sets <- c("catch22", "catch2", "catch24")

reticulate::use_python(python_to_use)

library(reticulate)

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
# Load libraries
################################################################################
library(feather)
library(tidyverse)
library(glue)
library(icesTAF)
library(cowplot)
library(correctR)
library(scales)
library(ggnewscale)
theme_set(theme_cowplot())

# Source visualisation script
source(glue("{github_dir}/data_visualisation/Manuscript_Draft_Visualisations_Helper.R"))

# Load in SPI info
SPI_info <- read.csv(glue("{github_dir}/data_visualisation/SPI_info.csv"))

# Load brain region info
UCLA_CNP_brain_region_info <- read.csv("~/data/UCLA_CNP/study_metadata/UCLA_CNP_Brain_Region_info.csv")
ABIDE_ASD_brain_region_info <- read.table("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_Harvard_Oxford_cort_prob_2mm_ROI_lookup.txt", sep=";", header = T) %>%
  mutate(Brain_Region = ifelse(Index==45, "Heschl's Gyrus (includes H1 and H2)", Brain_Region))
region_node_to_from <- read.csv("~/data/TS_feature_manuscript/node_to_from_structure.csv") %>%
  mutate(Study = ifelse(Study == "ABIDE", "ABIDE_ASD", "ABIDE"))

# Load participants included
UCLA_CNP_subjects_to_keep <- pyarrow_feather$read_feather("~/data/UCLA_CNP/processed_data/UCLA_CNP_filtered_sample_info_AROMA_2P_GMR_catch24_pyspi14.feather")

# Load study metadata
UCLA_CNP_metadata <- pyarrow_feather$read_feather("~/data/UCLA_CNP/study_metadata/UCLA_CNP_sample_metadata.feather") %>%
  mutate(Study = "UCLA_CNP") %>%
  filter(Sample_ID %in% UCLA_CNP_subjects_to_keep$Sample_ID)
ABIDE_ASD_metadata <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_sample_metadata.feather") %>%
  mutate(Study = "ABIDE_ASD")

# Load stats data
pairwise_balanced_accuracy_AUC_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_mixedsigmoid_scaler_balanced_accuracy_AUC_all_folds.feather"))
pairwise_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_mixedsigmoid_scaler_empirical_p_values.feather"))
pairwise_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_mixedsigmoid_scaler_null_balanced_accuracy_distributions.feather"))

combo_univariate_pairwise_balanced_accuracy_AUC_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combined_univariate_pairwise_mixedsigmoid_scaler_balanced_accuracy_AUC_all_folds.feather")) %>%
  mutate(Analysis_Type = "SPI_Univariate_Combo")
combo_univariate_pairwise_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combined_univariate_pairwise_mixedsigmoid_scaler_empirical_p_values.feather"))
combo_univariate_pairwise_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combined_univariate_pairwise_mixedsigmoid_scaler_null_balanced_accuracy_distributions.feather"))


# Take performance by repeats
pairwise_balanced_accuracy_AUC_by_repeats <- pairwise_balanced_accuracy_AUC_all_folds %>%
  group_by(Study, Comparison_Group, Analysis_Type, group_var, Repeat_Number) %>%
  summarise(Balanced_Accuracy_Across_Repeats = mean(Balanced_Accuracy, na.rm=T))
combo_univariate_pairwise_balanced_accuracy_AUC_by_repeats <- combo_univariate_pairwise_balanced_accuracy_AUC_all_folds %>%
  group_by(Study, Comparison_Group, Analysis_Type, group_var, Repeat_Number) %>%
  summarise(Balanced_Accuracy_Across_Repeats = mean(Balanced_Accuracy, na.rm=T))

# Load TPR/FPR data
pairwise_TPR_FPR <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_mixedsigmoid_scaler_ROC_TPR_FPR.feather")) %>%
  filter(Pairwise_Feature_Set == pairwise_feature_set) 
combo_univariate_pairwise_TPR_FPR <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combined_univariate_pairwise_mixedsigmoid_scaler_ROC_TPR_FPR.feather")) %>%
  filter(Pairwise_Feature_Set == pairwise_feature_set,
         Univariate_Feature_Set == univariate_feature_set) 


################################################################################
# Actual heatmap
################################################################################

combo_univariate_pairwise_p_values %>%
  filter(Pairwise_Feature_Set == pairwise_feature_set,
         Analysis_Type == "SPI_Univariate_Combo",
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
  scale_fill_gradientn(colors=c(alpha("darkgoldenrod2", 0.3), "darkgoldenrod2"), 
                       na.value=NA)  + 
  scale_y_discrete(labels = wrap_format(28)) +
  labs(fill = "Mean Balanced Accuracy (%)") +
  xlab("Clinical Group") +
  ylab("Pairwise SPI with All Univariate Info") +
  theme(legend.position="none",
        axis.text.y = element_text(size=10))
ggsave(glue("{plot_path}/Combo_univariate_pairwise_SPI_wise_heatmap.svg"),
       width=4.5, height=4.5, units="in", dpi=300)

################################################################################
# Compare SPIs with vs without univariate info
################################################################################

# Use correctR to test for difference across resamples for FTM vs catch22+FTM
run_correctR_group <- function(comparison_group, study, metadata, results_df) {
  num_subjects <- metadata %>%
    filter(Study == study, 
           Diagnosis %in% c("Control", comparison_group)) %>%
    distinct(Sample_ID) %>%
    nrow()
  training_size <- ceiling(0.9*num_subjects)
  test_size <- floor(0.1*num_subjects)
  
  data_for_correctR <- results_df %>%
    filter(Study == study, 
           Comparison_Group == comparison_group) %>%
    group_by(group_var) %>%
    filter(any(p_value_Bonferroni < 0.05)) %>%
    ungroup() %>%
    pivot_wider(id_cols = c(Repeat_Number, group_var), 
                names_from = Analysis_Type,
                values_from = Balanced_Accuracy_Across_Repeats) %>%
    dplyr::rename("x" = "Pairwise_SPI", "y" = "SPI_Univariate_Combo") %>%
    group_by(group_var) %>%
    group_split()
  
  res <- data_for_correctR %>%
    purrr::map_df(~ as.data.frame(resampled_ttest(x=.x$x, 
                                                  y=.x$y, 
                                                  n=10, 
                                                  n1=training_size, n2=test_size)) %>%
                    mutate(SPI = unique(.x$group_var))) %>%
    ungroup() %>%
    dplyr::rename("p_value_corr"="p.value") %>%
    mutate(p_value_corr_Bonferroni = p.adjust(p_value_corr, method="bonferroni"),
           Comparison_Group = comparison_group)
  
  return(res)
  
}

# metadata <- plyr::rbind.fill(UCLA_CNP_metadata, ABIDE_ASD_metadata)
results_df = plyr::rbind.fill(pairwise_balanced_accuracy_AUC_by_repeats %>% left_join(pairwise_p_values %>%
                                                                                        dplyr::select(Study:group_var, p_value_Bonferroni)), 
                              combo_univariate_pairwise_balanced_accuracy_AUC_by_repeats %>% left_join(combo_univariate_pairwise_p_values %>%
                                                                                                         dplyr::select(Study:group_var, p_value_Bonferroni)))

corrected_SPI_T_res <- 1:nrow(study_group_df) %>%
  purrr::map_df(~ run_correctR_group(comparison_group = study_group_df$Comparison_Group[.x],
                                     study = study_group_df$Study[.x],
                                     metadata = plyr::rbind.fill(UCLA_CNP_metadata, ABIDE_ASD_metadata),
                                     results_df = results_df)) 

plyr::rbind.fill(pairwise_p_values,
                 combo_univariate_pairwise_p_values) %>%
  dplyr::rename("SPI" = group_var) %>%
  semi_join(., corrected_SPI_T_res %>% dplyr::select(SPI, Comparison_Group)) %>%
  left_join(., corrected_SPI_T_res) %>%
  mutate(Analysis_Type = ifelse(Analysis_Type == "Pairwise_SPI", "SPI\nOnly", "SPI + Univariate\nRegion ×\nFeature"),
         individually_significant = ifelse(p_value_Bonferroni < 0.05, "Significant", "Not significant"),
         significant_diff_with_univariate = ifelse(p_value_corr_Bonferroni < 0.05, "Sig Diff", "No Sig Diff"),
         Analysis_Type = factor(Analysis_Type, levels=c("SPI\nOnly", "SPI + Univariate\nRegion ×\nFeature")),
         Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                      Comparison_Group == "Bipolar" ~ "BPD",
                                      T ~ Comparison_Group)) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD"))) %>%
  ggplot(data=., mapping=aes(x=Analysis_Type, y=100*Balanced_Accuracy_Across_Repeats,
                             group = SPI)) +
  geom_line(aes(color = Comparison_Group, 
                alpha = significant_diff_with_univariate), show.legend = FALSE) +
  scale_alpha_manual(values=c("Sig Diff" = 1, "No Sig Diff" = 0.2)) +
  scale_color_manual(values=c("Control" = "#5BB67B", 
                              "SCZ" = "#573DC7", 
                              "BPD" = "#D5492A", 
                              "ADHD" = "#0F9EA9", 
                              "ASD" = "#C47B2F")) +
  facet_wrap(Comparison_Group ~ ., ncol=1, scales="fixed", strip.position = "right") +
  new_scale_colour() +  # start a new scale
  geom_point(aes(color = individually_significant)) +
  scale_color_manual(values = c("gray60", "#5BB67B")) +
  # scale_x_discrete(labels = wrap_format(7)) +
  xlab("Analysis Type") +
  ylab("Mean Balanced Accuracy (%)") +
  scale_y_continuous(breaks=c(45, 55, 65)) +
  scale_x_discrete(expand=c(0,0.2,0,0.2)) +
  theme(legend.position = "bottom",
        strip.placement = "outside",
        strip.text.y.right = element_text(angle=0),
        plot.margin = margin(1,30,1,1, unit="pt"),
        legend.title = element_blank())
ggsave(glue("{plot_path}/SPI_with_vs_without_univariate_spaghetti.svg"),
       width=3.5, height=5, units="in", dpi=300)

################################################################################
# Plot ROC of top-performing features
top_SPIs_to_find_AUC <- combo_univariate_pairwise_p_values %>%
  filter(Analysis_Type == "SPI_Univariate_Combo") %>%
  group_by(Comparison_Group, Study) %>%
  slice_max(n=1, order_by=Balanced_Accuracy_Across_Repeats) %>%
  mutate(Group_Nickname =  case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                     Comparison_Group == "Bipolar" ~ "BPD",
                                     T ~ Comparison_Group)) %>%
  dplyr::select(Study, Comparison_Group, Group_Nickname, group_var, ROC_AUC_Across_Repeats)

combo_univariate_pairwise_TPR_FPR  %>%
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
ggsave(glue("{plot_path}/combined_univariate_pairwise_top_SPI_ROC.svg"),
       width=3, height=3, units="in", dpi=300)