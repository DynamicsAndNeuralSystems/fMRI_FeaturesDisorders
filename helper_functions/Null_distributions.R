library(caret)
library(dplyr)
library(rlist)
library(knitr)
library(kableExtra)

#-------------------------------------------------------------------------------
# Model-free shuffles
#-------------------------------------------------------------------------------

# Helper function to calculate accuracy + balanced accuracy per shuffle
calc_acc_bacc_for_shuffle <- function(x) {
  
  x <- factor(x, levels = unique(x))
  y <- sample(x, replace=F)
  
  # Calculate raw accuracy
  acc <- sum(x == y)/length(x)
  
  # Calculate balanced accuracy
  bacc <- caret::confusionMatrix(x, y)$byClass[["Balanced Accuracy"]]
  bacc <- ifelse(is.na(bacc), 0.5, bacc)
  
  output <- data.frame(accuracy = acc,
                       balanced_accuracy = bacc)
  
  return(output)
}

# Run given number of model free class label shuffles and calculate accuracy
run_model_free_n_shuffles <- function(rdata_path,
                                      feature_set = "catch22",
                                      num_shuffles = 1000000) {
  
  input_groups <- readRDS(paste0(rdata_path, "Filtered_subject_info_",
                                 feature_set, ".Rds")) %>%
    pull(group)
  
  # Run model-free shuffles
  null_res <- 1:num_shuffles %>%
    purrr::map_df(~ calc_acc_bacc_for_shuffle(input_groups)) %>%
    dplyr::mutate(Type = "null")
  
  return(null_res)
  
}

#-------------------------------------------------------------------------------
# Run given number of model-based permutations and calculate in-sample accuracy
#-------------------------------------------------------------------------------

run_null_model_n_permutations <- function(rdata_path,
                                          noise_procs = c("AROMA+2P", 
                                                          "AROMA+2P+GMR", 
                                                          "AROMA+2P+DiCER"),
                                          feature_set = "catch22",
                                          grouping_var = "Brain_Region",
                                          svm_feature_var = "Feature",
                                          test_package = "e1071",
                                          svm_kernel = "linear",
                                          num_permutations = 100,
                                          use_inv_prob_weighting = FALSE,
                                          use_SMOTE = FALSE) {
  
  nullOuts <- 1:num_permutations %>%
    purrr::map_df( ~ run_univariate_cv_svm_by_input_var(rdata_path = rdata_path,
                                                        svm_kernel = svm_kernel,
                                                        feature_set = feature_set,
                                                        test_package = test_package,
                                                        grouping_var = grouping_var,
                                                        svm_feature_var = svm_feature_var,
                                                        noise_procs = noise_procs,
                                                        use_inv_prob_weighting = use_inv_prob_weighting,
                                                        use_SMOTE = use_SMOTE,
                                                        shuffle_labels = T))
  
  
  null_res <- nullOuts %>%
    dplyr::mutate(Type = "null")
  
  return(null_res)
  
}

# Pairwise
run_null_model_n_permutations_pairwise <- function(pairwise_data,
                                                   noise_proc = "AROMA+2P+GMR",
                                                   feature_set = "pyspi_19",
                                                   test_package = "e1071",
                                                   svm_kernel = "linear",
                                                   grouping_var = "SPI",
                                                   svm_feature_var = "region_pair",
                                                   return_all_fold_metrics = FALSE,
                                                   SPI_directionality,
                                                   num_permutations = 50,
                                                   use_inv_prob_weighting = FALSE,
                                                   use_SMOTE = FALSE) {
  
  nullOuts <- 1:num_permutations %>%
    purrr::map_df( ~ run_pairwise_cv_svm_by_input_var(pairwise_data = pairwise_data,
                                                      SPI_directionality = SPI_directionality,
                                                      svm_kernel = svm_kernel,
                                                      grouping_var = grouping_var,
                                                      svm_feature_var = svm_feature_var,
                                                      test_package = test_package,
                                                      noise_proc = noise_proc,
                                                      return_all_fold_metrics = return_all_fold_metrics,
                                                      use_inv_prob_weighting = use_inv_prob_weighting,
                                                      use_SMOTE = use_SMOTE,
                                                      shuffle_labels = TRUE))
  
  
  null_res <- nullOuts %>%
    dplyr::mutate(Type = "null")
  
  return(null_res)
  
}

# Univariate + Pairwise Combined
# Pairwise
run_null_model_n_permutations_univariate_pairwise_combo <- function(univariate_data,
                                                                    univariate_feature_set = "catch22",
                                                                    pairwise_data,
                                                                    pairwise_feature_set = "pyspi_19",
                                                                    noise_proc = "AROMA+2P+GMR",
                                                                    test_package = "e1071",
                                                                    svm_kernel = "linear",
                                                                    grouping_var = "SPI",
                                                                    svm_feature_var = "region_pair",
                                                                    return_all_fold_metrics = FALSE,
                                                                    SPI_directionality,
                                                                    num_permutations = 40,
                                                                    use_inv_prob_weighting = FALSE,
                                                                    use_SMOTE = FALSE) {
  
  nullOuts <- 1:num_permutations %>%
    purrr::map_df( ~ run_combined_uni_pairwise_cv_svm_by_input_var(univariate_data = univariate_data,
                                                                   univariate_feature_set = univariate_feature_set,
                                                                   pairwise_data = pairwise_data,
                                                                   pairwise_feature_set = pairwise_feature_set,
                                                                   SPI_directionality = SPI_directionality,
                                                                   svm_kernel = svm_kernel,
                                                                   test_package = test_package,
                                                                   noise_proc = noise_proc,
                                                                   return_all_fold_metrics = return_all_fold_metrics,
                                                                   use_inv_prob_weighting = use_inv_prob_weighting,
                                                                   use_SMOTE = use_SMOTE,
                                                                   shuffle_labels = TRUE))
  
  
  null_res <- nullOuts %>%
    dplyr::mutate(Type = "null")
  
  return(null_res)
  
}

#-------------------------------------------------------------------------------
# Helper function to calculate empirical p-values based on null distribution
#-------------------------------------------------------------------------------

calc_empirical_nulls <- function(class_res,
                                 null_data,
                                 feature_set,
                                 is_main_data_averaged = TRUE,
                                 is_null_data_averaged = TRUE,
                                 grouping_var = "Brain_Region") {
  merged_list <- list()
  
  if (!("Sample_Type" %in% colnames(null_data))) {
    null_in <- null_data %>%
      mutate(Sample_Type = "In-sample")
    null_out <- null_data %>%
      mutate(Sample_Type = "Out-of-sample")
    null_data <- plyr::rbind.fill(null_in, null_out)
  }
  
  if (!("grouping_var" %in% colnames(class_res))) {
    class_res <- class_res %>%
      dplyr::rename("grouping_var" = grouping_var)
  }
  
  if (!is_main_data_averaged) {
    class_res <- class_res %>%
      group_by(Sample_Type, grouping_var, Noise_Proc, use_inv_prob_weighting, use_SMOTE) %>%
      summarise(accuracy_avg = mean(accuracy, na.rm=T),
                accuracy_SD = sd(accuracy, na.rm=T),
                balanced_accuracy_avg = mean(balanced_accuracy, na.rm=T),
                balanced_accuracy_SD = sd(balanced_accuracy, na.rm=T)) %>%
      dplyr::rename("accuracy"="accuracy_avg",
                    "balanced_accuracy"="balanced_accuracy_avg")
  } 

  main_res <- class_res %>%
    dplyr::select(grouping_var, Sample_Type, Noise_Proc, accuracy, balanced_accuracy) %>%
    mutate(Type = "main") %>%
    pivot_longer(cols=c(accuracy, balanced_accuracy),
                 names_to="Metric",
                 values_to="Value") %>%
    pivot_wider(id_cols = c(grouping_var, Noise_Proc, 
                            Type, Sample_Type),
                names_from = Metric, values_from = Value)
  
  for (group_var in unique(class_res$grouping_var)) {
    
    if ("grouping_var" %in% colnames(null_data)) {
      group_null <- null_data %>%
        dplyr::filter(grouping_var == group_var)
    } else {
      group_null <- null_data %>% 
        dplyr::mutate(grouping_var = group_var)
    }
    
    group_main <- subset(main_res, grouping_var == group_var)
    group_null$Type <- "null"
    
    # If null dataset is specific to each noise-processing method
    if ("Noise_Proc" %in% colnames(group_null)) {
      p_value_res <- plyr::rbind.fill(group_main,
                                      group_null) %>%
        group_by(grouping_var, Noise_Proc, Sample_Type) %>%
        dplyr::summarise(main_accuracy = unique(accuracy[Type=="main"]),
                         main_balanced_accuracy = unique(balanced_accuracy[Type=="main"]),
                         acc_p = 1 - stats::ecdf(accuracy[Type=="null"])(main_accuracy),
                         bal_acc_p = 1 - stats::ecdf(balanced_accuracy[Type=="null"])(main_balanced_accuracy),
        ) %>%
        ungroup() %>%
        distinct() %>%
        dplyr::rename("accuracy" = "main_accuracy",
                      "balanced_accuracy" = "main_balanced_accuracy") 
    } else {
      null_accuracy <- group_null$accuracy
      null_balanced_accuracy <- group_null$balanced_accuracy
      p_value_res <- group_main %>%
        group_by(grouping_var, Noise_Proc, Sample_Type) %>%
        dplyr::summarise(main_accuracy = accuracy,
                         main_balanced_accuracy = balanced_accuracy,
                         acc_p = 1 - stats::ecdf(null_accuracy)(main_accuracy),
                         bal_acc_p = 1 - stats::ecdf(null_balanced_accuracy)(main_balanced_accuracy),
        ) %>%
        ungroup() %>%
        distinct() %>%
        dplyr::rename("accuracy" = "main_accuracy",
                      "balanced_accuracy" = "main_balanced_accuracy") 
    }
    
    merged_list <- rlist::list.append(merged_list, p_value_res)
  }
  main_p_values <- do.call(plyr::rbind.fill, merged_list) %>%
    ungroup() %>%
    mutate(Noise_Proc = factor(Noise_Proc, levels = c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER"))) %>%
    group_by(Noise_Proc) %>%
    mutate(acc_p_adj = p.adjust(acc_p, method="BH"),
           bal_acc_p_adj = p.adjust(bal_acc_p, method="BH"),
           feature_set = feature_set) %>%
    arrange(Noise_Proc)
  
  return(main_p_values)
}

#-------------------------------------------------------------------------------
# Truncate p-values to N significant digits in scientific notation
#-------------------------------------------------------------------------------

truncate_p_values <- function(pvalue_df, N=3) {
  truncated_pvalue_df <- pvalue_df %>%
    mutate(acc_p = scales::scientific(acc_p, digits = N),
           bal_acc_p = scales::scientific(bal_acc_p, digits = N),
           acc_p_adj = scales::scientific(acc_p_adj, digits = N),
           bal_acc_p_adj = scales::scientific(bal_acc_p_adj, digits = N))
}

#-------------------------------------------------------------------------------
# Print kable with the number of significant features per noise-processing method
#-------------------------------------------------------------------------------

summarise_num_sig_features <- function(pvalue_df) {
  pvalue_df %>%
    mutate(Noise_Proc = factor(Noise_Proc, levels = noise_procs)) %>%
    group_by(Noise_Proc, Sample_Type) %>%
    summarise(num_sig_acc = sum(acc_p < 0.05),
              num_sig_acc_fdr = sum(acc_p_adj < 0.05),
              num_sig_bacc = sum(bal_acc_p < 0.05),
              num_sig_bacc_fdr = sum(bal_acc_p_adj < 0.05)) %>%
    kable(.) %>%
    kable_styling(full_width = F)
}