#------------------------------------
# Functions to run linear SVM
# Options to use cross-validation with
# Inverse probability weighting or SMOTE
#------------------------------------

#--------------------------------------
# Author: Annie Bryant, 16 May 2022
#--------------------------------------

require(rlist)
library(tidyverse)
library(e1071)
library(kernlab)
library(cowplot)
library(caret)
library(smotefamily)
theme_set(theme_cowplot())

#-------------------------------------------------------------------------------
# In-sample SVM with option to use inverse probability weighting or SMOTE
#-------------------------------------------------------------------------------
in_sample_linear_SVM <- function(input_data,
                                 svm_kernel = "linear",
                                 sample_wts = list("Control" = 1,
                                                   "Schz" = 1),
                                 shuffle_labels = FALSE) {
  
  # Shuffle labels if specified
  if (shuffle_labels) {
    input_data <- transform(input_data, group = sample(group, replace = FALSE))
  }
  
  # Run linear SVM
  svmModel <- e1071::svm(factor(group) ~ .,
                         kernel = svm_kernel,
                         cost = 1,
                         data = input_data,
                         class.weights = sample_wts)
  
  # Generate in-sample predictions based on SVM model
  pred <- predict(svmModel, input_data)
  input_data$group <- factor(input_data$group, levels = levels(pred))
  
  # Calculate accuracy and balanced accuracy
  accuracy <- sum(pred == input_data$group)/length(pred)
  balanced_accuracy <- caret::confusionMatrix(data=pred, 
                                              reference=input_data$group)$byClass[["Balanced Accuracy"]]
  balanced_accuracy <- ifelse(is.na(balanced_accuracy), 0.5, balanced_accuracy)
  
  # Compile results into a dataframe
  df_res <- data.frame(accuracy = accuracy,
                       balanced_accuracy = balanced_accuracy)
  
  # Return resulting dataframe
  return(df_res)
  
}

#-------------------------------------------------------------------------------
# k-fold CV SVM with option to use inverse probability weighting or SMOTE
#-------------------------------------------------------------------------------
k_fold_CV_linear_SVM <- function(input_data,
                                 k = 10,
                                 svm_kernel = "linear",
                                 sample_wts = list("Control" = 1,
                                                   "Schz" = 1),
                                 use_SMOTE = FALSE,
                                 shuffle_labels = FALSE) {
  
  # Shuffle labels if specified
  if (shuffle_labels) {
    input_data <- transform(input_data, group = sample(group, replace = FALSE))
  }
  
  # Specify that group is a factor so that createFolds creates stratified folds
  input_data$group <- factor(input_data$group)
  
  # Create train/test data folds
  flds <- caret::createFolds(input_data$group, k = k, list = TRUE, returnTrain = FALSE)
  
  # Initialize list for performance metrics
  accuracy_list <- list()
  balanced_accuracy_list <- list()
  
  # Iterate over folds 1 through k
  for (i in 1:k) {
    
    # Define test and train data
    test_i <- flds[[i]]
    train_i <- setdiff(1:nrow(input_data), test_i)
    
    train_data <- input_data[train_i, ]
    test_data <- input_data[test_i, ]
    
    # Apply SMOTE to training data only if indicated
    if (use_SMOTE) {
      train_data <- smotefamily::SMOTE(train_data[,-1], train_data$group, K = 5)$data %>%
        dplyr::rename("group" = "class")
    }
    
    # Run linear SVM on fold
    svmModel <- e1071::svm(factor(group) ~ .,
                           kernel = svm_kernel,
                           cost = 1,
                           data = train_data,
                           class.weights = sample_wts)
    
    # Generate out-of-sample predictions based on SVM model
    pred <- predict(svmModel, test_data)
    test_data$group <- factor(test_data$group, levels = levels(pred))
    
    # Calculate accuracy and balanced accuracy
    accuracy <- sum(pred == test_data$group)/length(pred)
    balanced_accuracy <- caret::confusionMatrix(reference=test_data$group, 
                                                data=pred)$byClass[["Balanced Accuracy"]]
    
    # Append results to list for given fold
    accuracy_list[[i]] <- accuracy
    balanced_accuracy_list[[i]] <- balanced_accuracy
  } 
  
  # Calculate summary statistics for accuracy and balanced accuracy
  accuracy_avg <- mean(unlist(accuracy_list), na.rm=T)
  accuracy_sd <- sd(unlist(accuracy_list), na.rm=T)
  balanced_accuracy_avg <- mean(unlist(balanced_accuracy_list), na.rm=T)
  balanced_accuracy_sd <- sd(unlist(balanced_accuracy_list), na.rm=T)
  
  # Compile results into a dataframe
  df_res <- data.frame(accuracy = accuracy_avg,
                       accuracy_SD = accuracy_sd,
                       balanced_accuracy = balanced_accuracy_avg,
                       balanced_accuracy_SD = balanced_accuracy_sd)
}

#-------------------------------------------------------------------------------
# Run simple in-sample multi-feature linear SVM by given grouping var
#-------------------------------------------------------------------------------
run_in_sample_svm_by_input_var <- function(rdata_path,
                                           svm_kernel = "linear",
                                           test_package = "e1071",
                                           grouping_var = "Brain_Region",
                                           svm_feature_var = "Feature",
                                           noise_procs = c("AROMA+2P", 
                                                           "AROMA+2P+GMR", 
                                                           "AROMA+2P+DiCER"),
                                           use_inv_prob_weighting = FALSE,
                                           use_SMOTE = FALSE,
                                           shuffle_labels = FALSE) {
  
  class_res_list <- list()
  for (noise_proc in noise_procs) {
    # Clean up names
    noise_label <- gsub("\\+", "_", noise_proc)
    
    # Load catch22 data for current noise processing method
    feature_matrix <- readRDS(paste0(rdata_path, sprintf("UCLA_%s_catch22_zscored.Rds", 
                                                         noise_label)))      
    
    if (svm_feature_var == "Feature") {
      svm_feature_var_name = "names"
      grouping_var_name = "Brain_Region"
      grouping_var_vector <- unique(feature_matrix$Brain_Region)
      
    } else if (svm_feature_var == "Brain_Region") {
      svm_feature_var_name = svm_feature_var
      grouping_var_name = "names"
      grouping_var_vector <- unique(feature_matrix$names)
    } else {
      svm_feature_var_name = "Combo"
      grouping_var_name = "Group_Var"
      
      feature_matrix <- feature_matrix %>%
        unite("Combo", c("Brain_Region", "names"), sep="_", remove=FALSE)
      
      grouping_var_vector <- c("All")
    }
    
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
    } else {
      sample_wts <- list("Control" = 1,
                         "Schz" = 1)
    }
    
    for (group_var in grouping_var_vector) {
      if (grouping_var == "Combo") {
        data_for_SVM <- feature_matrix %>%
          dplyr::group_by(Subject_ID) %>%
          dplyr::filter(!any(is.na(values))) %>%
          dplyr::ungroup() %>%
          dplyr::select(Subject_ID, group, Combo, values) %>%
          tidyr::pivot_wider(id_cols = c(Subject_ID, group),
                             names_from = Combo,
                             values_from 
                             = values) %>%
          dplyr::select(-Subject_ID)
        
      } else {
        # Otherwise iterate over each separate group
        data_for_SVM <- subset(feature_matrix, get(grouping_var_name) == group_var) %>%
          dplyr::group_by(Subject_ID, get(grouping_var_name)) %>%
          dplyr::filter(!any(is.na(values))) %>%
          dplyr::ungroup() %>%
          dplyr::select(Subject_ID, group, svm_feature_var_name, values) %>%
          tidyr::pivot_wider(id_cols = c(Subject_ID, group),
                             names_from = svm_feature_var_name,
                             values_from 
                             = values) %>%
          dplyr::select(-Subject_ID)
      }
      
      # Apply SMOTE if indicated
      if (use_SMOTE) {
        data_for_SVM <- smotefamily::SMOTE(data_for_SVM[,-1], data_for_SVM$group, K = 5)$data %>%
          dplyr::rename("group" = "class")
      }
      
      # Pass data_for_SVM to in_sample_linear_SVM
      SVM_results <- in_sample_linear_SVM(input_data = data_for_SVM,
                                          svm_kernel = svm_kernel,
                                          sample_wts = sample_wts,
                                          shuffle_labels = shuffle_labels) %>%
        mutate(grouping_var = group_var,
               Noise_Proc = noise_proc)
      
      # Append results to list
      class_res_list <- rlist::list.append(class_res_list, SVM_results)
    }
  }
  
  # Combine results from all regions into a dataframe
  class_res_df <- do.call(plyr::rbind.fill, class_res_list)
  
  # Return dataframe
  return(class_res_df)
}

#-------------------------------------------------------------------------------
# Run 10-fold cross-validated multi-feature linear SVM by given grouping var
#-------------------------------------------------------------------------------

run_cv_svm_by_input_var <- function(rdata_path,
                                    svm_kernel = "linear",
                                    test_package = "e1071",
                                    grouping_var = "Brain_Region",
                                    svm_feature_var = "Feature",
                                    noise_procs = c("AROMA+2P", 
                                                    "AROMA+2P+GMR", 
                                                    "AROMA+2P+DiCER"),
                                    use_inv_prob_weighting = FALSE,
                                    use_SMOTE = FALSE,
                                    shuffle_labels = FALSE) {
  
  class_res_list <- list()
  for (noise_proc in noise_procs) {
    # Clean up names
    noise_label <- gsub("\\+", "_", noise_proc)
    
    # Load catch22 z-scored data for current noise processing method
    feature_matrix <- readRDS(paste0(rdata_path, sprintf("UCLA_%s_catch22_zscored.Rds", 
                                                         noise_label)))      
    
    if (svm_feature_var == "Feature") {
      svm_feature_var_name = "names"
      grouping_var_name = "Brain_Region"
      grouping_var_vector <- unique(feature_matrix$Brain_Region)
      
    } else if (svm_feature_var == "Brain_Region") {
      svm_feature_var_name = svm_feature_var
      grouping_var_name = "names"
      grouping_var_vector <- unique(feature_matrix$names)
    } else {
      svm_feature_var_name = "Combo"
      grouping_var_name = "Group_Var"
      
      feature_matrix <- feature_matrix %>%
        unite("Combo", c("Brain_Region", "names"), sep="_", remove=F)
      
      grouping_var_vector <- c("All")
    }
    
    # Define sample weights
    # Default is 1 and 1 if use_inv_prob_weighting is not included
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
    } else {
      sample_wts <- list("Control" = 1,
                         "Schz" = 1)
    }
    
    # Reshape data from long to wide for SVM
    for (group_var in grouping_var_vector) {
      if (grouping_var == "Combo") {
        data_for_SVM <- feature_matrix %>%
          dplyr::group_by(Subject_ID) %>%
          dplyr::filter(!any(is.na(values))) %>%
          dplyr::ungroup() %>%
          dplyr::select(Subject_ID, group, Combo, values) %>%
          tidyr::pivot_wider(id_cols = c(Subject_ID, group),
                             names_from = Combo,
                             values_from 
                             = values) %>%
          dplyr::select(-Subject_ID)
        
      } else {
        # Otherwise iterate over each separate group
        data_for_SVM <- subset(feature_matrix, get(grouping_var_name) == group_var) %>%
          dplyr::group_by(Subject_ID, get(grouping_var_name)) %>%
          dplyr::filter(!any(is.na(values))) %>%
          dplyr::ungroup() %>%
          dplyr::select(Subject_ID, group, svm_feature_var_name, values) %>%
          tidyr::pivot_wider(id_cols = c(Subject_ID, group),
                             names_from = svm_feature_var_name,
                             values_from 
                             = values) %>%
          dplyr::select(-Subject_ID)
      }
      
      # Pass data_for_SVM to in_sample_linear_SVM
      SVM_results <- k_fold_CV_linear_SVM(input_data = data_for_SVM,
                                          k = 10,
                                          svm_kernel = svm_kernel,
                                          sample_wts = sample_wts,
                                          use_SMOTE = use_SMOTE,
                                          shuffle_labels = shuffle_labels) %>%
        mutate(grouping_var = group_var,
               Noise_Proc = noise_proc)
      
      # Append results to list
      class_res_list <- rlist::list.append(class_res_list,
                                           SVM_results)
    }
  }
  
  # Combine results from all regions into a dataframe
  class_res_df <- do.call(plyr::rbind.fill, class_res_list)
  
  # Return dataframe
  return(class_res_df)
}






