################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/Library/CloudStorage/OneDrive-TheUniversityofSydney(Students)/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/combo_univariate_pairwise_results/")
TAF::mkdir(plot_path)

# python_to_use <- "~/.conda/envs/pyspi/bin/python3"
python_to_use <- "/Users/abry4213/anaconda3/envs/pyspi/bin/python3"
univariate_feature_set <- "catch25"
pairwise_feature_set <- "pyspi14"
data_path <- "~/data/TS_feature_manuscript"
UCLA_CNP_sample_metadata <- "~/data/UCLA_CNP/study_metadata/UCLA_CNP_sample_metadata.feather"
ABIDE_ASD_sample_metadata <- "~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_sample_metadata.feather"
study_group_df <- data.frame(Study = c(rep("UCLA_CNP", 3), "ABIDE_ASD"),
                             Num_Samples = c(166, 157, 167, 1150), 
                             Noise_Proc = c(rep("AROMA+2P+GMR",3), "FC1000"),
                             Comparison_Group = c("Schizophrenia", "ADHD", "Bipolar", "ASD"))
univariate_feature_sets <- c("catch22", "catch2", "catch25")

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

# Load participants included
UCLA_CNP_subjects_to_keep <- pyarrow_feather$read_feather("~/data/UCLA_CNP/processed_data/UCLA_CNP_filtered_sample_info_AROMA_2P_GMR_catch25_pyspi14.feather")

# Load study metadata
UCLA_CNP_metadata <- pyarrow_feather$read_feather("~/data/UCLA_CNP/study_metadata/UCLA_CNP_sample_metadata.feather") %>%
  mutate(Study = "UCLA_CNP") %>%
  filter(Sample_ID %in% UCLA_CNP_subjects_to_keep$Sample_ID)
ABIDE_ASD_metadata <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_sample_metadata.feather") %>%
  mutate(Study = "ABIDE_ASD")

# Load pairwise stats data
pairwise_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_balanced_accuracy_all_folds.feather"))
# Compute mean + SD performance across all folds
pairwise_balanced_accuracy <- pairwise_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Pairwise_Feature_Set, Analysis_Type, group_var) %>%
  reframe(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
          Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T))
pairwise_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_empirical_p_values.feather")) %>%
  dplyr::select(-Balanced_Accuracy_Across_Folds) %>%
  left_join(., pairwise_balanced_accuracy)
pairwise_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_pairwise_null_balanced_accuracy_distributions.feather"))

# Load combined pairwise + univariate data
combo_univariate_pairwise_balanced_accuracy_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combined_univariate_pairwise_balanced_accuracy_all_folds.feather")) %>%
  mutate(Analysis_Type = "SPI_Univariate_Combo")
# Compute mean + SD performance across all folds
combo_univariate_pairwise_balanced_accuracy <- combo_univariate_pairwise_balanced_accuracy_all_folds %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Pairwise_Feature_Set, Analysis_Type, group_var) %>%
  reframe(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
          Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T))
combo_univariate_pairwise_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combined_univariate_pairwise_empirical_p_values.feather")) %>%
  dplyr::select(-Balanced_Accuracy_Across_Folds) %>%
  left_join(., combo_univariate_pairwise_balanced_accuracy)
combo_univariate_pairwise_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_combined_univariate_pairwise_null_balanced_accuracy_distributions.feather"))

################################################################################
# Actual heatmap
################################################################################

combo_univariate_pairwise_p_values %>%
  filter(Pairwise_Feature_Set == pairwise_feature_set,
         Analysis_Type == "SPI_Univariate_Combo",
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
  scale_fill_gradientn(colors=c(alpha("darkgoldenrod2", 0.3), "darkgoldenrod2"), 
                       na.value=NA)  + 
  scale_y_discrete(labels = wrap_format(28)) +
  labs(fill = "Mean Balanced Accuracy (%)") +
  xlab("Disorder") +
  ylab("FC + Brain-wide local dynamics") +
  theme(legend.position="none",
        axis.text.y = element_text(size=10))
ggsave(glue("{plot_path}/Combo_univariate_pairwise_SPI_wise_heatmap.svg"),
       width=4, height=4.5, units="in", dpi=300)

################################################################################
# Compare SPIs with vs without univariate info
################################################################################

repkfold_ttest <- function(data, n1, n2, k, r){
  
  # Arg checks
  
  '%ni%' <- Negate('%in%')
  
  if("model" %ni% colnames(data) || "values" %ni% colnames(data) || "k" %ni% colnames(data) || "r" %ni% colnames(data)){
    stop("data should contain at least four columns called 'model', 'values', 'k', and 'r'.")
  }
  
  if(!is.numeric(data$values) || !is.numeric(data$k) || !is.numeric(data$r)){
    stop("data should be a data.frame with only numerical values in columns 'values', 'k', and 'r'.")
  }
  
  if(!is.numeric(n1) || !is.numeric(n2) || !is.numeric(k) || !is.numeric(r) ||
     length(n1) != 1 || length(n2) != 1 || length(k) != 1 || length(r) != 1){
    stop("n1, n2, k, and r should all be integer scalars.")
  }
  
  if(length(unique(data$model)) != 2){
    stop("Column 'model' in data should only have two unique labels (one for each model to compare).")
  }
  
  # Calculations
  
  d <- c()
  
  for(i in 1:k){
    for(j in 1:r){
      x <- data[data$k == i, ]
      x <- x[x$r == j, ]
      d <- c(d, x[x$model == unique(x$model)[1], c("values")] - x[x$model == unique(x$model)[2], c("values")]) # Differences
    }
  }
  
  # Catch for when there is zero difference(s) between the models
  
  if (sum(unlist(d)) == 0) {
    tmp <- data.frame(statistic = 0, p.value = 1)
  } else{
    
    statistic <- mean(unlist(d), na.rm = TRUE) / sqrt(stats::var(unlist(d), na.rm = TRUE) * ((1/(k * r)) + (n2/n1))) # Calculate t-statistic
    df <- n1 + n2 - 2
    
    if(statistic < 0){
      p.value <- stats::pt(statistic, (k * r) - 1) # p-value for left tail
    } else{
      p.value <- stats::pt(statistic, (k * r) - 1, lower.tail = FALSE) # p-value for right tail
    }
    
    tmp <- data.frame(statistic = statistic, p.value = p.value, df = df)
  }
  
  return(tmp)
}

# Use correctR to test for difference across resamples for FTM vs catch22+FTM
run_correctR_group <- function(comparison_group, study, metadata, results_df) {
  # Find number of subjects for the specified comparison group
  num_subjects <- metadata %>%
    filter(Study == study, 
           Diagnosis %in% c("Control", comparison_group)) %>%
    distinct(Sample_ID) %>%
    nrow()
  
  # Compute the training and test fold sizes for 10-fold CV
  training_size <- ceiling(0.9*num_subjects)
  test_size <- floor(0.1*num_subjects)
  
  # Prep the resulting balanced accuracies with vs without univariate data
  data_for_correctR <- results_df %>%
    filter(Study == study, 
           Comparison_Group == comparison_group) %>%
    group_by(group_var) %>%
    filter(any(p_value_HolmBonferroni < 0.05)) %>%
    ungroup() %>%
    dplyr::rename("model" = "Analysis_Type",
                  "SPI" = "group_var",
                  "k" = "Fold",
                  "r" = "Repeat_Number",
                  "values" = "Balanced_Accuracy") %>%
    dplyr::select(model, SPI, k, r, values) %>%
    dplyr::mutate(r = r + 1) %>%
    group_by(SPI) %>%
    group_split()

  res <- data_for_correctR %>%
    purrr::map_df(~ as.data.frame(repkfold_ttest(data = .x %>% dplyr::select(-SPI), 
                                                 n1 = training_size,
                                                 n2 = test_size,
                                                 k = 10,
                                                 r = 10)) %>%
                    mutate(SPI = unique(.x$SPI))) %>%
    ungroup() %>%
    dplyr::rename("p_value_corr"="p.value") %>%
    mutate(p_value_corr_HolmBonferroni = p.adjust(p_value_corr, method="bonferroni"),
           Comparison_Group = comparison_group)
  
  return(res)
  
}

results_df = plyr::rbind.fill(pairwise_balanced_accuracy_all_folds %>% left_join(pairwise_p_values %>%
                                                                                       dplyr::select(Study:group_var, p_value_HolmBonferroni)), 
                              combo_univariate_pairwise_balanced_accuracy_all_folds %>% left_join(combo_univariate_pairwise_p_values %>%
                                                                                                        dplyr::select(Study:group_var, p_value_HolmBonferroni)))

corrected_SPI_T_res <- 1:nrow(study_group_df) %>%
  purrr::map_df(~ run_correctR_group(comparison_group = study_group_df$Comparison_Group[.x],
                                     study = study_group_df$Study[.x],
                                     metadata = plyr::rbind.fill(UCLA_CNP_metadata, ABIDE_ASD_metadata),
                                     results_df = results_df)) %>%
  left_join(., SPI_info, by=c("SPI"="pyspi_name"))

plyr::rbind.fill(pairwise_p_values,
                 combo_univariate_pairwise_p_values) %>%
  dplyr::rename("SPI" = group_var) %>%
  semi_join(., corrected_SPI_T_res %>% dplyr::select(SPI, Comparison_Group)) %>%
  left_join(., corrected_SPI_T_res) %>%
  mutate(Analysis_Type = ifelse(Analysis_Type == "Pairwise_SPI", "FC", "FC + Local\nDynamics"),
         individually_significant = ifelse(p_value_HolmBonferroni < 0.05, "pcorr < 0.05", "Not significant"),
         significant_diff_with_univariate = ifelse(p_value_corr_HolmBonferroni < 0.05, "Sig Diff", "No Sig Diff"),
         Analysis_Type = factor(Analysis_Type, levels=c("FC", "FC + Local\nDynamics")),
         Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                      Comparison_Group == "Bipolar" ~ "BP",
                                      T ~ Comparison_Group)) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BP", "ADHD", "ASD"))) %>%
  ggplot(data=., mapping=aes(x=Analysis_Type, y=100*Balanced_Accuracy_Across_Folds,
                             group = SPI)) +
  geom_line(aes(color = Comparison_Group, 
                alpha = significant_diff_with_univariate), show.legend = FALSE) +
  scale_alpha_manual(values=c("Sig Diff" = 1, "No Sig Diff" = 0.2)) +
  scale_color_manual(values=c("Control" = "#5BB67B", 
                              "SCZ" = "#573DC7", 
                              "BP" = "#D5492A", 
                              "ADHD" = "#0F9EA9", 
                              "ASD" = "#C47B2F")) +
  facet_wrap(Comparison_Group ~ ., ncol=1, scales="fixed", strip.position = "left") +
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
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle=0, face="bold"),
        plot.margin = margin(1,30,1,1, unit="pt"),
        legend.title = element_blank())
ggsave(glue("{plot_path}/SPI_with_vs_without_univariate_spaghetti.svg"),
       width=3.5, height=5, units="in", dpi=300)

# Save stats to a CSV file as a supplementary table
corrected_SPI_T_res %>%
  dplyr::select(Figure_name, Comparison_Group, statistic, df, p_value_corr_HolmBonferroni) %>%
  dplyr::rename("SPI" = "Figure_name",
                "Disorder" = "Comparison_Group",
                "T-statistic" = "statistic",
                "Holm-Bonferroni corrected P" = "p_value_corr_HolmBonferroni") %>%
  write.csv(., glue("{github_dir}/plots/Manuscript_Draft/tables/combined_pairwise_univariate_corr_T_test_res.csv"),
            row.names = F)