#------------------------------------
# This script sets out to produce a
# function for reading in matlab time
# series files into R
#------------------------------------

#--------------------------------------
# Author: Annie Bryant, 23 March 2022
#--------------------------------------

require(rlist)
library(tidyverse)
library(theft)
library(caret)
library(broom)

#-------------------------------------------------------------------------------
# Calculate balanced accuracy in caret
#-------------------------------------------------------------------------------
calculate_balanced_accuracy <- function(data, lev = NULL, model = NULL) {
  
  # Calculate balanced accuracy from confusion matrix as the average of class recalls as per https://arxiv.org/pdf/2008.05756.pdf
  
  cm <- as.matrix(caret::confusionMatrix(data$pred, data$obs)$table)
  
  recall <- 1:nrow(cm) %>%
    purrr::map(~ calculate_recall(cm, x = .x)) %>%
    unlist()
  
  balanced_accuracy <- sum(recall) / length(recall)
  
  # Calculate accuracy
  
  accuracy <- sum(diag(cm)) / sum(cm)
  
  # Return results
  
  out <- c(accuracy, balanced_accuracy)
  names(out) <- c("Accuracy", "Balanced_Accuracy")
  return(out)
}

calculate_in_sample_balanced_accuracy <- function(cm) {
  recall <- 1:nrow(cm) %>%
    purrr::map(~ calculate_recall(cm, x = .x)) %>%
    unlist()
  
  balanced_accuracy <- sum(recall) / length(recall)
  
  # Calculate accuracy
  
  accuracy <- sum(diag(cm)) / sum(cm)
  
  # Return results
  
  out <- c(accuracy, balanced_accuracy)
  names(out) <- c("Accuracy", "Balanced_Accuracy")
  return(out)
}


#-------------------------------------------------------------------------------
# Run theft's multivariate classifier
#-------------------------------------------------------------------------------
run_theft_multivar_classifier_ROI_feature <- function(rdata_path,
                                                      test_method = "svmLinear",
                                                      noise_proc = "AROMA+2P",
                                                      use_balanced_accuracy = FALSE) {
  
  # Clean up names
  test_label <- gsub("-|\\.", "_", test_method)
  noise_label <- gsub("\\+", "_", noise_proc)
  
  # Load corresponding feature matrix
  feature_matrix <- readRDS(paste0(rdata_path, sprintf("UCLA_%s_catch22.Rds", 
                                                       noise_label)))
  
  # Prepare feature matrix for theft such that each feature + ROI combo
  # is treated as a feature
  feature_matrix_prepped <- feature_matrix %>%
    tidyr::unite("names", c("names", "Brain_Region"), remove=T)
  
  # Run theft's multivariable classifier
  multi_classifier_outputs <- fit_multi_feature_classifier(feature_matrix_prepped,
                                                           id_var = "Subject_ID",
                                                           group_var = "group",
                                                           by_set = F,
                                                           test_method = test_method,
                                                           use_balanced_accuracy = TRUE,
                                                           use_k_fold = T,
                                                           num_folds = 10,
                                                           use_empirical_null = T,
                                                           null_testing_method = "model free shuffles",
                                                           p_value_method = "gaussian",
                                                           num_permutations = 10000)
  
  # Compile results
  res_df <- multi_classifier_outputs$TestStatistics %>%
    dplyr::mutate(Noise_Proc = noise_proc,
                  Test_Method = test_method,
                  Norm_Method = norm_method)
    
  # Return dataframe
  return(res_df)
}

#-------------------------------------------------------------------------------
# Run simple e1071 in-sample multi-feature SVM
#-------------------------------------------------------------------------------
run_in_sample_e1071_SVM_by_all_feature_region <- function(rdata_path,
                                                          svm_kernel = "linear",
                                                          noise_procs = c("AROMA+2P", 
                                                                          "AROMA+2P+GMR", 
                                                                          "AROMA+2P+DiCER")) {
  noise_proc_res_list <- list()
  for (noise_proc in noise_procs) {
    # Clean up names
    noise_label <- gsub("\\+", "_", noise_proc)
    
    # Load catch22 data for current noise processing method
    feature_matrix <- readRDS(paste0(rdata_path, sprintf("UCLA_%s_catch22.Rds", 
                                                         noise_label)))      
    # Get control/schz proportions
    sample_props <- feature_matrix %>%
      dplyr::group_by(Subject_ID, Brain_Region) %>%
      dplyr::filter(!any(is.na(values))) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(Subject_ID, group) %>%
      dplyr::summarise(control_prop = sum(group=="Control") / n(),
                       schz_prop = sum(group=="Schz")/n())
    
    # Convert to sample weights based on inverse of probability
    sample_wts <- list("Control" = 1/sample_props$control_prop,
                       "Schz" = 1/sample_props$schz_prop)
    
      
      # Combine feature + ROI and convert data to wide format
      data_for_svm <- feature_matrix %>%
        dplyr::group_by(Subject_ID, Brain_Region) %>%
        dplyr::filter(!any(is.na(values))) %>%
        dplyr::ungroup() %>%
        tidyr::unite("names", c(names, Brain_Region), sep="_") %>%
        dplyr::select(Subject_ID, group, names, values) %>%
        tidyr::pivot_wider(id_cols = c(Subject_ID, group),
                           names_from = names,
                           values_from 
                           = values) %>%
        dplyr::select(-Subject_ID)
      
      # Run SVM with supplied kernel type
      svmModel <- e1071::svm(factor(group) ~ .,
                             kernel = svm_kernel,
                             data = data_for_svm, 
                             class.weights = sample_wts)
      
      # Generate in-sample predictions based on SVM model
      pred <- stats::predict(svmModel)
      data_for_svm$group <- factor(data_for_svm$group, levels = levels(pred))
      
      # Calculate accuracy and balanced accuracy
      accuracy <- sum(pred == data_for_svm$group)/length(pred)
      cm <- as.matrix(caret::confusionMatrix(pred, data_for_svm$group)$table)
      balanced_accuracy <- calculate_in_sample_balanced_accuracy(cm = cm)
      
      # Compile results into a dataframe
      class_res <- data.frame(Accuracy = accuracy,
                                  Balanced_Accuracy = balanced_accuracy,
                                  Noise_Proc = noise_proc)
      
      # Append results to list
      noise_proc_res_list <- rlist::list.append(noise_proc_res_list,
                                                class_res)
    }
  # Combine results into a dataframe
  res_df <- do.call(plyr::rbind.fill, noise_proc_res_list)
  
  # Return dataframe
  return(res_df)
}

#-------------------------------------------------------------------------------
# Run linear multi-feature SVM in caret with e1071 per brain region and feature
#-------------------------------------------------------------------------------
run_caret_multi_SVM_by_region_feature_inv_prop <- function(rdata_path,
                                                   use_inv_prob_weighting = TRUE,
                                                   noise_procs = c("AROMA+2P", 
                                                                   "AROMA+2P+GMR", 
                                                                   "AROMA+2P+DiCER")) {
  
  noise_proc_res_list <- list()
  for (noise_proc in noise_procs) {
    # Clean up names
    noise_label <- gsub("\\+", "_", noise_proc)
    
    # Load catch22 data for current noise processing method
    feature_matrix <- readRDS(paste0(rdata_path, sprintf("UCLA_%s_catch22.Rds", 
                                                         noise_label)))   
    # Get control/schz proportions
    sample_props <- feature_matrix %>%
      dplyr::group_by(Subject_ID, Brain_Region) %>%
      dplyr::filter(!any(is.na(values))) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(Subject_ID, group) %>%
      dplyr::summarise(control_prop = sum(group=="Control") / n(),
                       schz_prop = sum(group=="Schz")/n())
    
      # Prep catch22 feature data for SVM
      data_for_svm <- feature_matrix %>%
        dplyr::rename("id" = "Subject_ID") %>%
        tidyr::unite("names", c(names, Brain_Region), sep="_") %>%
        dplyr::select(id, group, names, values) %>%
        tidyr::pivot_wider(id_cols = c("id", "group"), names_from = "names", values_from = "values") %>%
        dplyr::select_if(~sum(!is.na(.)) > 0) %>%
        dplyr::select(where(~dplyr::n_distinct(.) > 1))  %>%
        tidyr::drop_na() %>%
        janitor::clean_names() %>%
        tidyr::pivot_longer(cols = 3:ncol(.), names_to = "names", values_to = "values") %>%
        dplyr::mutate(method = gsub("_.*", "\\1", names)) %>%
        dplyr::mutate(group = make.names(group),
                      group = as.factor(group))  %>%
        tidyr::pivot_wider(id_cols = c("id", "group"), names_from = "names", values_from = "values") %>%
        dplyr::select(-c(id))
      
      # Calculate model weights as inverse probability
      model_weights <- ifelse(data_for_svm$group == "Control", 
                              1/sample_props$control_prop,
                              1/sample_props$schz_prop)
      
      # Train SVM model
      fitControl <- caret::trainControl(method = "cv",
                                        number = 10,
                                        savePredictions = "all",
                                        summaryFunction = calculate_balanced_accuracy,
                                        classProbs = TRUE)
      
      # Run e1071 SVM with caret
      mod <- caret::train(group ~ .,
                          data = data_for_svm,
                          method = "svmLinear",
                          trControl = fitControl,
                          metric = "Balanced_Accuracy",
                          maximize = T,
                          weights = model_weights,
                          preProcess = c("center", "scale", "nzv"))
      
      # Generate predictions
      preds <- mod$pred
      data_for_svm$group <- factor(data_for_svm$group, levels = levels(preds$pred))
      cm <- list(caret::confusionMatrix(preds$pred, data_for_svm$group)$table)
      
      # Get accuracy + balanced accuracy results
      mod_res <- mod$results %>%
        mutate(Test_Method = "linear_SVM",
               Noise_Proc = noise_proc,
               Confusion_Matrix = I(cm))
      
      # Append results to list
      noise_proc_res_list <- rlist::list.append(noise_proc_res_list,
                                                mod_res)
    
  }
  # Combine results into a dataframe
  class_res_df <- do.call(plyr::rbind.fill, noise_proc_res_list)
  
  # Return dataframe
  return(class_res_df)
}
