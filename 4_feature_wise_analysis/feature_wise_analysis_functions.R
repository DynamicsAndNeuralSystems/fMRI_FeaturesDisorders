#------------------------------------
# This script sets out to produce a
# function for reading in matlab time
# series files into R
#------------------------------------

#--------------------------------------
# Author: Annie Bryant, 14 April 2022
#--------------------------------------

require(rlist)
library(tidyverse)
library(theft)
library(caret)
library(broom)
library(kernlab)
library(yardstick)


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
# Run simple in-sample multi-feature kernlab linear SVM by brain region
#-------------------------------------------------------------------------------
run_in_sample_ksvm_by_region <- function(rdata_path,
                                         svm_kernel = "linear",
                                         noise_procs = c("AROMA+2P", 
                                                         "AROMA+2P+GMR", 
                                                         "AROMA+2P+DiCER"),
                                         use_inv_prob_weighting = FALSE,
                                         upsample_minority = FALSE) {
  
  ROI_wise_class_res_list <- list()
  for (noise_proc in noise_procs) {
    # Clean up names
    noise_label <- gsub("\\+", "_", noise_proc)
    
    # Load catch22 data for current noise processing method
    feature_matrix <- readRDS(paste0(rdata_path, sprintf("UCLA_%s_catch22.Rds", 
                                                         noise_label)))      
    
    
    if (use_inv_prob_weighting) {
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
    }
    
    if (upsample_minority) {
      ctrl_subjects <- subset(feature_matrix, group=="Control") %>% 
        distinct(Subject_ID) %>% 
        pull(Subject_ID)
      
      schz_subjects_upsampled <- subset(feature_matrix, group=="Schz") %>%
        distinct(Subject_ID) %>%
        pull(Subject_ID) %>%
        sample(., length(ctrl_subjects), replace=T)
      
      upsampled_data <- data.frame(Subject_ID = c(ctrl_subjects, 
                                                  schz_subjects_upsampled)) %>%
        mutate(Unique_ID = make.unique(Subject_ID))
    }
    
    # Iterate over each brain region (ROI) for e1071 SVM
    for (this_ROI in unique(feature_matrix$Brain_Region)) {
      
      if (upsample_minority) {
        # Subset region data and convert to wide format
        data_for_svm <- upsampled_data %>%
          left_join(., feature_matrix) %>%
          filter(Brain_Region == this_ROI) %>%
          dplyr::group_by(Unique_ID, Brain_Region) %>%
          dplyr::filter(!any(is.na(values))) %>%
          dplyr::ungroup() %>%
          dplyr::select(Unique_ID, group, names, values) %>%
          tidyr::pivot_wider(id_cols = c(Unique_ID, group),
                             names_from = names,
                             values_from 
                             = values) %>%
          dplyr::select(-Unique_ID)
      } else {
        data_for_svm <- subset(feature_matrix, Brain_Region==this_ROI) %>%
          dplyr::group_by(Subject_ID, Brain_Region) %>%
          dplyr::filter(!any(is.na(values))) %>%
          dplyr::ungroup() %>%
          dplyr::select(Subject_ID, group, names, values) %>%
          tidyr::pivot_wider(id_cols = c(Subject_ID, group),
                             names_from = names,
                             values_from 
                             = values) %>%
          dplyr::select(-Subject_ID)
      }
      
      # Run SVM with supplied kernel type
      svmModel <- kernlab::ksvm(factor(group) ~ .,
                                type = "C-svc",
                                kernel = "vanilladot",
                                data = data_for_svm,
                                prob.model=F)
      
      # Generate in-sample predictions based on SVM model
      pred <- predict(svmModel, data_for_svm)
      data_for_svm$group <- factor(data_for_svm$group, levels = levels(pred))
      
      # Calculate accuracy and balanced accuracy
      accuracy <- sum(pred == data_for_svm$group)/length(pred)
      cm <- as.matrix(caret::confusionMatrix(pred, data_for_svm$group)$table)
      balanced_accuracy <- calculate_in_sample_balanced_accuracy(cm = cm)
      
      # Compile results into a dataframe
      region_df_res <- data.frame(Brain_Region = this_ROI,
                                  Accuracy = accuracy,
                                  Balanced_Accuracy = balanced_accuracy,
                                  Noise_Proc = noise_proc)
      
      # Append results to list
      ROI_wise_class_res_list <- rlist::list.append(ROI_wise_class_res_list,
                                                    region_df_res)
    }
  }
  # Combine results from all regions into a dataframe
  ROI_wise_class_res_df <- do.call(plyr::rbind.fill, ROI_wise_class_res_list)
  
  # Return dataframe
  return(ROI_wise_class_res_df)
}

#-------------------------------------------------------------------------------
# Run theft's multivariate classifier for a given catch22 feature
#-------------------------------------------------------------------------------
run_theft_multivar_classifier_by_feature <- function(rdata_path,
                                                     norm_method = "z-score",
                                                     test_method = "svmLinear",
                                                     noise_proc = "AROMA+2P",
                                                     use_balanced_accuracy = FALSE) {
  
  # Clean up names
  test_label <- gsub("-|\\.", "_", test_method)
  noise_label <- gsub("\\+", "_", noise_proc)
  
  # Load corresponding feature matrix
  feature_matrix <- readRDS(paste0(rdata_path, sprintf("UCLA_%s_catch22.Rds", 
                                                       noise_label)))
  
  # Iterate over each catch22 feature in the feature matrix
  feature_test_statistics_list <- list()
  for (this_feature in unique(feature_matrix$names)) {
    # Subset to given feature
    feature_matrix_feature <- subset(feature_matrix, names==this_feature) %>%
      dplyr::select(-names) %>%
      dplyr::rename("names" = "Brain_Region")
    
    
    cat("\nRunning theft multivariable classifier for", this_feature, noise_proc, "\n")
    
    multi_classifier_outputs <- fit_multi_feature_classifier(feature_matrix_feature,
                                                             id_var = "Subject_ID",
                                                             group_var = "group",
                                                             by_set = F,
                                                             test_method = test_method,
                                                             use_balanced_accuracy = use_balanced_accuracy,
                                                             use_k_fold = T,
                                                             num_folds = 10,
                                                             use_empirical_null = T,
                                                             null_testing_method = "model free shuffles",
                                                             p_value_method = "gaussian",
                                                             num_permutations = 10000)
    
    # Compile results
    test_stat <- multi_classifier_outputs$TestStatistics %>%
      mutate(catch22_Feature = this_feature,
             Noise_Proc = noise_proc,
             Test_Method = test_method,
             Norm_Method = norm_method)
    
    # Append results to list
    feature_test_statistics_list <- rlist::list.append(feature_test_statistics_list, test_stat)
  }
  
  # Combine list of results into a dataframe
  feature_test_statistics <- do.call(plyr::rbind.fill, feature_test_statistics_list)
  
  # Return dataframe
  return(feature_test_statistics)
}

#-------------------------------------------------------------------------------
# Plot distribution of classifier accuracies across features 
# along with control vs schz proportions
#-------------------------------------------------------------------------------
plot_class_acc_w_props <- function(class_res,
                                   cv = TRUE,
                                   rdata_path,
                                   noise_procs = c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER")) {
  
  # Calculate the proportion of control subjects after omitting NAs for each method
  ctrl_prop_list <- list()
  for (noise_proc in noise_procs) {
    # Clean up names
    noise_label <- gsub("\\+", "_", noise_proc)
    
    # Load corresponding feature matrix
    feature_matrix <- readRDS(paste0(rdata_path, sprintf("UCLA_%s_catch22.Rds", 
                                                         noise_label)))
    
    # Calculate proportion of controls after dropping NA
    ctrl_proportion <- feature_matrix %>%
      group_by(Subject_ID, Brain_Region) %>%
      filter(!any(is.na(values))) %>%
      ungroup() %>%
      distinct(Subject_ID, group) %>%
      summarise(ctrl_prop = sum(group=="Control") / n()) %>%
      mutate(Noise_Proc = noise_proc)
    
    ctrl_prop_list <- rlist::list.append(ctrl_prop_list, ctrl_proportion)
  }
  ctrl_prop <- do.call(plyr::rbind.fill, ctrl_prop_list)
  
  # Define x-axis label
  xlab <- ifelse(cv, "Mean statistic over 10-fold CV", "In-sample statistic")
  
  # Plot accuracy + balanced accuracy in histograms
  # Control subject proportion is highlighted for accuracy, 0.5 is highlighted for balanced accuracy
  class_res %>%
    dplyr::select(catch22_Feature, Noise_Proc, Accuracy, Balanced_Accuracy) %>%
    left_join(., ctrl_prop) %>%
    pivot_longer(cols=c(Accuracy, Balanced_Accuracy),
                 names_to = "Metric",
                 values_to = "Value") %>%
    mutate(Metric = gsub("_", " ", Metric),
           ctrl_prop = ifelse(Metric == "Balanced Accuracy", 0.5, ctrl_prop)) %>%
    mutate(Noise_Proc = factor(Noise_Proc, levels = noise_procs)) %>%
    ggplot(data=., mapping=aes(x=Value)) +
    geom_histogram(fill="lightsteelblue", bins=50) +
    geom_vline(aes(xintercept = ctrl_prop), linetype=2, color="gray30") +
    facet_grid(Noise_Proc ~ Metric, scales="free", switch="y") +
    xlab(xlab) +
    ylab("Number of features") +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle=0))
}

#-------------------------------------------------------------------------------
# Run simple e1071 in-sample SVM by catch22 feature 
# with inverse probability weighting
#-------------------------------------------------------------------------------
run_in_sample_e1071_SVM_by_feature <- function(rdata_path,
                                    svm_kernel = "linear",
                                    noise_procs = c("AROMA+2P", 
                                                    "AROMA+2P+GMR", 
                                                    "AROMA+2P+DiCER")) {
  
  feature_wise_class_res_list <- list()
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
    
    # Iterate over each catch22 feature for e1071 SVM
    for (this_feature in unique(feature_matrix$names)) {
      
      # Subset feature data and convert to wide format
      feature_data <- subset(feature_matrix, names==this_feature) %>%
        dplyr::group_by(Subject_ID, Brain_Region) %>%
        dplyr::filter(!any(is.na(values))) %>%
        dplyr::ungroup() %>%
        dplyr::select(Subject_ID, group, Brain_Region, values) %>%
        tidyr::pivot_wider(id_cols = c(Subject_ID, group),
                           names_from = Brain_Region,
                           values_from 
                           = values) %>%
        dplyr::select(-Subject_ID)
      
      # Run SVM with supplied kernel type
      svmModel <- e1071::svm(factor(group) ~ .,
                             kernel = svm_kernel,
                             data = feature_data, 
                             class.weights = sample_wts)
      
      # Generate in-sample predictions based on SVM model
      pred <- stats::predict(svmModel)
      feature_data$group <- factor(feature_data$group, levels = levels(pred))
      
      # Calculate accuracy and balanced accuracy
      accuracy <- sum(pred == feature_data$group)/length(pred)
      cm <- as.matrix(caret::confusionMatrix(pred, feature_data$group)$table)
      balanced_accuracy <- calculate_in_sample_balanced_accuracy(cm = cm)
      
      # Compile results into a dataframe
      region_df_res <- data.frame(catch22_Feature = this_feature,
                                  Accuracy = accuracy,
                                  Balanced_Accuracy = balanced_accuracy,
                                  Noise_Proc = noise_proc)
      
      # Append results to list
      feature_wise_class_res_list <- rlist::list.append(feature_wise_class_res_list,
                                                    region_df_res)
    }
  }
  # Combine results from all regions into a dataframe
  feature_wise_class_res_df <- do.call(plyr::rbind.fill, feature_wise_class_res_list)
  
  # Return dataframe
  return(feature_wise_class_res_df)
}


#-------------------------------------------------------------------------------
# Run linear multi-ROI SVM in caret with e1071 per catch22 feature
#-------------------------------------------------------------------------------
run_caret_multi_SVM_by_feature_inv_prop <- function(rdata_path,
                                          use_inv_prob_weighting = TRUE,
                                          noise_procs = c("AROMA+2P", 
                                                          "AROMA+2P+GMR", 
                                                          "AROMA+2P+DiCER")) {
  
  feature_wise_class_res_list <- list()
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
    
    for (this_feature in unique(feature_matrix$names)) {
      cat("\nNow running linear SVM for", this_feature, noise_label, "\n")
      
      feature_matrix_feature <- subset(feature_matrix, names==this_feature)
      
      # Prep catch22 feature data for SVM
      data_for_svm <- feature_matrix_feature %>%
        dplyr::rename("id" = "Subject_ID") %>%
        dplyr::select(id, group, Brain_Region, values) %>%
        tidyr::pivot_wider(id_cols = c("id", "group"), names_from = "Brain_Region", values_from = "values") %>%
        dplyr::select_if(~sum(!is.na(.)) > 0) %>%
        dplyr::select(where(~dplyr::n_distinct(.) > 1))  %>%
        tidyr::drop_na() %>%
        janitor::clean_names() %>%
        tidyr::pivot_longer(cols = 3:ncol(.), names_to = "Brain_Region", values_to = "values") %>%
        dplyr::mutate(method = gsub("_.*", "\\1", Brain_Region)) %>%
        dplyr::mutate(group = make.names(group),
                      group = as.factor(group))  %>%
        tidyr::pivot_wider(id_cols = c("id", "group"), names_from = "Brain_Region", values_from = "values") %>%
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
      ROI_res <- mod$results %>%
        mutate(catch22_Feature = this_feature,
               Test_Method = "linear_SVM",
               Noise_Proc = noise_proc,
               Confusion_Matrix = I(cm))
      
      # Append results to list
      feature_wise_class_res_list <- rlist::list.append(feature_wise_class_res_list,
                                                    ROI_res)
    }
    
  }
  # Combine results from all regions into a dataframe
  feature_wise_class_res_df <- do.call(plyr::rbind.fill, feature_wise_class_res_list)
  
  # Return dataframe
  return(feature_wise_class_res_df)
}