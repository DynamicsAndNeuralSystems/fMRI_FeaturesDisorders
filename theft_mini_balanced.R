github_dir <- "D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders/"
data_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/"
rdata_path <- paste0(data_path, "Rdata/")
set.seed(127)

#### load libraries
library(theft)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

#### sample subjects
subjects <- readRDS(paste0(rdata_path, "Filtered_subject_info_catch22.Rds"))

subjects_balanced <- subjects %>%
  group_by(group) %>%
  sample_n(47)

saveRDS(subjects_balanced, paste0(rdata_path, "Filtered_subject_info_catch22_small_balanced.Rds"))

#### load catch22 data for these subjects using AROMA+2P+GMR

feature_matrix_full <- readRDS(paste0(rdata_path, "UCLA_AROMA_2P_GMR_catch22_filtered.Rds"))

feature_matrix_small_balanced <- subset(feature_matrix_full,
                                        Subject_ID %in% subjects_balanced$Subject_ID)
saveRDS(feature_matrix_small_balanced, 
        file = paste0(rdata_path, "UCLA_AROMA_2P_GMR_catch22_filtered_small_balanced.Rds"))

################################################################################
# ROI-wise
################################################################################
null_model_res_list <- list()
for (brain_region in unique(feature_matrix_small_balanced$Brain_Region)) {
  # Full dataset
  ROI_feature_matrix <- feature_matrix_full %>%
    dplyr::filter(Brain_Region == brain_region)
  
  # Small balanced dataset
  ROI_feature_matrix_small <- feature_matrix_small_balanced %>%
    dplyr::filter(Brain_Region == brain_region)
  
  # Model-free shuffles
  ## Full
  ROI_wise_theft_model_free_shuffles_full <- fit_multi_feature_classifier(ROI_feature_matrix, 
                                                                          id_var = "Subject_ID", 
                                                                          group_var = "group", 
                                                                          by_set = TRUE, 
                                                                          test_method = "svmLinear",
                                                                          use_balanced_accuracy = TRUE,
                                                                          use_k_fold = TRUE,
                                                                          num_folds = 10,
                                                                          use_empirical_null = FALSE,
                                                                          null_testing_method = "model free shuffles",
                                                                          p_value_method = "empirical",
                                                                          num_permutations = 10,
                                                                          seed = 127)$RawClassificationResults %>%
    mutate(category = ifelse(method == "catch22", "main", "model free shuffle null"),
           dataset = "full")
  
  ## Small balanced
  ROI_wise_theft_model_free_shuffles_small <- fit_multi_feature_classifier(ROI_feature_matrix_small, 
                                                                           id_var = "Subject_ID", 
                                                                           group_var = "group", 
                                                                           by_set = TRUE, 
                                                                           test_method = "svmLinear",
                                                                           use_balanced_accuracy = TRUE,
                                                                           use_k_fold = TRUE,
                                                                           num_folds = 10,
                                                                           use_empirical_null = FALSE,
                                                                           null_testing_method = "model free shuffles",
                                                                           p_value_method = "empirical",
                                                                           num_permutations = 10,
                                                                           seed = 127)$RawClassificationResults %>%
    mutate(category = ifelse(method == "catch22", "main", "model free shuffle null"),
           dataset = "small balanced")
  
  # Null model fits
  ## Full
  ROI_wise_theft_null_model_fit_full <- fit_multi_feature_classifier(ROI_feature_matrix, 
                                                                id_var = "Subject_ID", 
                                                                group_var = "group", 
                                                                by_set = TRUE, 
                                                                test_method = "svmLinear",
                                                                use_balanced_accuracy = TRUE,
                                                                use_k_fold = TRUE,
                                                                num_folds = 10,
                                                                use_empirical_null =  TRUE,
                                                                null_testing_method = "null model fits",
                                                                p_value_method = "empirical",
                                                                num_permutations = 10,
                                                                seed = 127)$RawClassificationResults %>%
    filter(category != "Main") %>%
    mutate(category = "null model fits",
           dataset = "full")
  
  ## Small balanced
  ROI_wise_theft_null_model_fit_small <- fit_multi_feature_classifier(ROI_feature_matrix_small, 
                                                                      id_var = "Subject_ID", 
                                                                      group_var = "group", 
                                                                      by_set = TRUE, 
                                                                      test_method = "svmLinear",
                                                                      use_balanced_accuracy = TRUE,
                                                                      use_k_fold = TRUE,
                                                                      num_folds = 10,
                                                                      use_empirical_null =  TRUE,
                                                                      null_testing_method = "null model fits",
                                                                      p_value_method = "empirical",
                                                                      num_permutations = 10,
                                                                      seed = 127)$RawClassificationResults %>%
    filter(category != "Main") %>%
    mutate(category = "null model fits",
           dataset = "small balanced")
  
  ROI_null_res <- plyr::rbind.fill(ROI_wise_theft_model_free_shuffles_full,
                                   ROI_wise_theft_model_free_shuffles_small,
                                   ROI_wise_theft_null_model_fit_full,
                                   ROI_wise_theft_null_model_fit_small) %>%
    mutate(Brain_Region = brain_region)
  
  # Append to list
  null_model_res_list <- rlist::list.append(null_model_res_list,
                                            ROI_null_res)
}

ROI_wise_null_model_res <- do.call(plyr::rbind.fill, null_model_res_list)
saveRDS(ROI_wise_null_model_res, file = paste0(rdata_path, "UCLA_AROMA_2P_GMR_catch22_theft_null_models.Rds"))


ROI_wise_null_model_res %>%
  filter(str_detect(category, "null")) %>%
  ggplot(data = ., mapping = aes(x = balanced_accuracy)) +
  geom_histogram(mapping = aes(fill = category, y=0.5*..density..), 
                 alpha = 0.6, position="identity",
                 bins = 50) +
  xlab("10-fold CV Balanced Accuracy") +
  ylab("Scaled Density") +
  labs(fill = "Null Type") +
  ggtitle("Distribution of Null Balanced Accuracy Values from theft\nfor ROI-wise catch22 feature linear SVM") +
  facet_grid(dataset ~ ., scales = "free", switch = "both") +
  theme(plot.title = element_text(hjust=0.5),
        strip.placement = "outside")

################################################################################
# Combo-wise
################################################################################