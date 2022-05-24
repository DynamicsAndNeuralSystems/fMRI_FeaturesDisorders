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
                                      noise_procs = c("AROMA+2P", 
                                                      "AROMA+2P+GMR", 
                                                      "AROMA+2P+DiCER"),
                                      num_shuffles = 1000000) {
  
  set.seed(123)
  
  # Initialize empty list
  null_list <- list()
  
  # Iterate over each noise-processing method
  for (noise_proc in noise_procs) {
    # Clean up names
    noise_label <- gsub("\\+", "_", noise_proc)
    
    # Load catch22 data for current noise processing method
    feature_matrix <- readRDS(paste0(rdata_path, sprintf("UCLA_%s_catch22_zscored.Rds", 
                                                         noise_label)))   
    
    input_groups <- feature_matrix %>% 
      distinct(Subject_ID, group) %>%
      mutate(group = factor(group, levels = c("Schz", "Control"))) %>%
      pull(group)
    
    # Run model-free shuffles
    output <- 1:num_shuffles %>%
      purrr::map_df(~ calc_acc_bacc_for_shuffle(input_groups)) %>%
      mutate(Noise_Proc = noise_proc)
    
    null_list <- rlist::list.append(null_list, output)
  }
  null_res <- do.call(plyr::rbind.fill, null_list) %>%
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
                                          grouping_var = "Brain_Region",
                                          svm_feature_var = "Feature",
                                          test_package = "e1071",
                                          svm_kernel = "linear",
                                          num_permutations = 100,
                                          cross_validate = FALSE,
                                          use_inv_prob_weighting = FALSE,
                                          use_SMOTE = FALSE) {
  
  # Initialize empty list
  null_list <- list()
  
  if (cross_validate) {
    nullOuts <- 1:num_permutations %>%
      purrr::map_df( ~ run_cv_svm_by_input_var(rdata_path = rdata_path,
                                               svm_kernel = svm_kernel,
                                               test_package = test_package,
                                               grouping_var = grouping_var,
                                               svm_feature_var = svm_feature_var,
                                               noise_procs = noise_procs,
                                               use_inv_prob_weighting = use_inv_prob_weighting,
                                               use_SMOTE = use_SMOTE,
                                               shuffle_labels = T))
    
  } else {
    nullOuts <- 1:num_permutations %>%
      purrr::map_df( ~ run_in_sample_svm_by_input_var(rdata_path = rdata_path,
                                                      svm_kernel = svm_kernel,
                                                      test_package = test_package,
                                                      grouping_var = grouping_var,
                                                      svm_feature_var = svm_feature_var,
                                                      noise_procs = noise_procs,
                                                      use_inv_prob_weighting = use_inv_prob_weighting,
                                                      use_SMOTE = use_SMOTE,
                                                      shuffle_labels = T))
  }
  null_res <- nullOuts %>%
    dplyr::mutate(Type = "null")
  
  return(null_res)
  
}

#-------------------------------------------------------------------------------
# Helper function to calculate empirical p-values based on null distribution
#-------------------------------------------------------------------------------

calc_empirical_nulls <- function(class_res,
                                 null_data,
                                 grouping_var = "Brain_Region") {
  merged_list <- list()
  main_res <- class_res %>%
    dplyr::select(grouping_var, Noise_Proc, accuracy, balanced_accuracy) %>%
    mutate(Type = "main")
  
  for (group_var in unique(class_res$grouping_var)) {
    
    group_null <- null_data %>% mutate(grouping_var = group_var)
    group_main <- subset(main_res, grouping_var == group_var)
    
    group_merged <- plyr::rbind.fill(group_main,
                                     group_null) %>%
      group_by(grouping_var, Noise_Proc) %>%
      dplyr::summarise(main_accuracy = unique(accuracy[Type=="main"]),
                       main_balanced_accuracy = unique(balanced_accuracy[Type=="main"]),
                       num_acc_greater = sum(main_accuracy > accuracy[Type=="null"],
                                             na.rm=T),
                       num_bal_acc_greater = sum(main_balanced_accuracy > 
                                                   balanced_accuracy[Type=="null"],
                                                 na.rm=T),
                       total_null = sum(Type == "null"),
                       acc_p = 1 - (num_acc_greater/total_null),
                       bal_acc_p = 1 - (num_bal_acc_greater/total_null)) %>%
      ungroup() %>%
      distinct() %>%
      dplyr::rename("accuracy" = "main_accuracy",
                    "balanced_accuracy" = "main_balanced_accuracy")
    
    merged_list <- rlist::list.append(merged_list, group_merged)
  }
  main_p_values <- do.call(plyr::rbind.fill, merged_list) %>%
    ungroup() %>%
    group_by(Noise_Proc) %>%
    mutate(acc_p_adj = p.adjust(acc_p, method="BH"),
           bal_acc_p_adj = p.adjust(bal_acc_p, method="BH"))
  
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
    group_by(Noise_Proc) %>%
    summarise(num_sig_acc = sum(acc_p < 0.05),
              num_sig_acc_fdr = sum(acc_p_adj < 0.05),
              num_sig_bacc = sum(bal_acc_p < 0.05),
              num_sig_bacc_fdr = sum(bal_acc_p_adj < 0.05)) %>%
    kable(.) %>%
    kable_styling(full_width = F)
}