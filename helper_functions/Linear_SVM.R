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
                                 use_SMOTE = FALSE,
                                 shuffle_labels = FALSE) {
  
  # Shuffle labels if specified
  if (shuffle_labels) {
    input_data <- transform(input_data, group = sample(group, replace = FALSE))
  }
  
  # Apply SMOTE to training data only if indicated
  if (use_SMOTE) {
    input_data <- smotefamily::SMOTE(input_data[,-1], input_data$group, K = 5)$data %>%
      dplyr::rename("group" = "class")
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
                                 c_values = c(1),
                                 svm_kernel = "linear",
                                 sample_wts = list("Control" = 1,
                                                   "Schz" = 1),
                                 use_SMOTE = FALSE,
                                 shuffle_labels = FALSE,
                                 return_all_fold_metrics = FALSE) {
  
  # Shuffle labels if specified
  if (shuffle_labels) {
    input_data <- transform(input_data, group = sample(group, replace = FALSE))
  }
  
  # Specify that group is a factor so that createFolds creates stratified folds
  input_data$group <- factor(input_data$group)
  
  # Create train/test data folds
  flds <- caret::createFolds(input_data$group, k = k, list = TRUE, returnTrain = FALSE)
  
  # Initialize list for performance metrics
  c_value_list <- list()
  in_accuracy_list <- list()
  in_balanced_accuracy_list <- list()
  out_accuracy_list <- list()
  out_balanced_accuracy_list <- list()
  
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
    
    # Iterate over each value of c
    for (c in c_values) {
      svmModel <- e1071::svm(factor(group) ~ .,
                             kernel = svm_kernel,
                             cost = c,
                             data = train_data,
                             class.weights = sample_wts)
      
      # Generate in-sample predictions based on SVM model
      in_sample_pred <- predict(svmModel, train_data)
      train_data$group <- factor(train_data$group, levels = levels(in_sample_pred))
      
      # Calculate accuracy and balanced accuracy
      in_accuracy <- sum(in_sample_pred == train_data$group)/length(in_sample_pred)
      in_balanced_accuracy <- caret::confusionMatrix(reference=train_data$group, 
                                                     data=in_sample_pred)$byClass[["Balanced Accuracy"]]
      
      # Generate out-of-sample predictions based on SVM model
      out_sample_pred <- predict(svmModel, test_data)
      test_data$group <- factor(test_data$group, levels = levels(out_sample_pred))
      
      # Calculate accuracy and balanced accuracy
      out_accuracy <- sum(out_sample_pred == test_data$group)/length(out_sample_pred)
      out_balanced_accuracy <- caret::confusionMatrix(reference=test_data$group, 
                                                      data=out_sample_pred)$byClass[["Balanced Accuracy"]]
      
      # Append results to list for given fold
      c_value_list <- rlist::list.append(c_value_list, c)
      in_accuracy_list <- rlist::list.append(in_accuracy_list, in_accuracy)
      in_balanced_accuracy_list <- rlist::list.append(in_balanced_accuracy_list, in_balanced_accuracy)
      out_accuracy_list <- rlist::list.append(out_accuracy_list, out_accuracy)
      out_balanced_accuracy_list <- rlist::list.append(out_balanced_accuracy_list, out_balanced_accuracy)
    }
    
  } 
  
  df_res_in_sample <- data.frame(Null_Iteration = rep(1:k, length(c_values)),
                                 c_value = unlist(c_value_list),
                                 accuracy = unlist(in_accuracy_list),
                                 balanced_accuracy = unlist(in_balanced_accuracy_list),
                                 Sample_Type = "In-sample")
  df_res_out_sample <- data.frame(Null_Iteration = rep(1:k, length(c_values)),
                                  c_value = unlist(c_value_list),
                                  accuracy = unlist(out_accuracy_list),
                                  balanced_accuracy = unlist(out_balanced_accuracy_list),
                                  Sample_Type = "Out-of-sample")
  
  df_res <- plyr::rbind.fill(df_res_in_sample, df_res_out_sample)
  
  if (return_all_fold_metrics) {
    return(df_res)
  } else {
    df_res_avg <- df_res %>%
      group_by(Sample_Type, c_value) %>%
      summarise(accuracy_avg = mean(accuracy, na.rm=T),
                balanced_accuracy_avg = mean(balanced_accuracy, na.rm=T),
                accuracy_SD = sd(accuracy, na.rm=T),
                balanced_accuracy_SD = sd(balanced_accuracy, na.rm=T)) %>%
      dplyr::rename("accuracy" = "accuracy_avg",
                    "balanced_accuracy" = "balanced_accuracy_avg")
    
    return(df_res_avg)
  }
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
# Run simple pairwise PYSPI in-sample multi-feature linear SVM by given grouping var
#-------------------------------------------------------------------------------
run_pairwise_SVM_by_SPI <- function(pairwise_data,
                                    svm_kernel = "linear",
                                    test_package = "e1071",
                                    noise_proc = "AROMA+2P",
                                    cross_validate = FALSE,
                                    return_all_fold_metrics = FALSE,
                                    use_inv_prob_weighting = FALSE,
                                    use_SMOTE = FALSE,
                                    shuffle_labels = FALSE) {
  
  # Initialize results list for SPI-wise in-sample linear SVM
  class_res_list <- list()
  
  # Iterate over each PYSPI SPI
  for (this_SPI in unique(all_pyspi_data$SPI)) {
    # Identify whether current SPI is directed or undirected
    directionality <- SPI_directionality %>% filter(SPI == this_SPI) %>% pull(Direction)
    
    # Reshape SPI data to contain brain region names
    SPI_data <- subset(all_pyspi_data, SPI == this_SPI)  %>%
      mutate(comparison = row_number()) %>%
      pivot_longer(cols = c(brain_region_1,
                            brain_region_2),
                   names_to = "Region_Number",
                   values_to = "Index") %>%
      left_join(ROI_index) %>%
      dplyr::select(-Index) %>%
      pivot_wider(id_cols = c("Subject_ID", "group", "SPI", "value", "comparison"),
                  names_from = "Region_Number",
                  values_from = "ROI") %>%
      dplyr::select(-comparison) 
    
    # Combine brain regions into new pairwise column depending on directionality
    if (directionality == "Undirected") {
      SPI_data <- SPI_data %>%
        mutate(region_pair = ifelse(brain_region_1 < brain_region_2,
                                    paste0(brain_region_1, "_", brain_region_2),
                                    paste0(brain_region_2, "_", brain_region_1)))
    } else if (directionality == "Directed") {
      SPI_data <- SPI_data %>%
        mutate(region_pair = paste0(brain_region_1, "_", brain_region_2))
    }
    
    # Pivot data from long to wide for SVM
    data_for_SVM <- SPI_data %>%
      dplyr::select(Subject_ID, group, region_pair, value) %>%
      mutate(value = round(value, 8)) %>%
      distinct() %>%
      pivot_wider(id_cols = c(Subject_ID, group),
                  names_from = region_pair, 
                  values_from = value) %>%
      dplyr::select(-Subject_ID) %>%
      drop_na()
    
    # Get sample weights if inverse probability weighting flag is applied
    if (use_inv_prob_weighting) {
      # Get control/schz proportions
      sample_props <- data_for_SVM %>%
        ungroup() %>%
        dplyr::summarise(control_prop = sum(group=="CONTROL") / n(),
                         schz_prop = sum(group=="SCHZ")/n())
      
      # Convert to sample weights based on inverse of probability
      sample_wts <- list("CONTROL" = 1/sample_props$control_prop,
                         "SCHZ" = 1/sample_props$schz_prop)
    } else {
      sample_wts <- list("CONTROL" = 1,
                         "SCHZ" = 1)
    }
    
    # Pass data_for_SVM to in_sample_linear_SVM
    if (nrow(data_for_SVM ) > 0) {
      
      cat("\nNow running linear SVM for", this_SPI, "\n")
      if (cross_validate) {
        SVM_results <- k_fold_CV_linear_SVM(input_data = data_for_SVM,
                                            k = 10,
                                            svm_kernel = svm_kernel,
                                            sample_wts = sample_wts,
                                            use_SMOTE = use_SMOTE,
                                            return_all_fold_metrics = return_all_fold_metrics,
                                            shuffle_labels = shuffle_labels)
      } else {
        SVM_results <- in_sample_linear_SVM(input_data = data_for_SVM,
                                            svm_kernel = svm_kernel,
                                            sample_wts = sample_wts,
                                            use_SMOTE = use_SMOTE,
                                            shuffle_labels = shuffle_labels)
      }
      SVM_results <- SVM_results %>%
        mutate(SPI = this_SPI,
               Noise_Proc = noise_proc,
               use_inv_prob_weighting = use_inv_prob_weighting,
               use_SMOTE = use_SMOTE)
      
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
# Run linear SVM by iterating over number of principal components (PCs)
#-------------------------------------------------------------------------------
run_SVM_from_PCA <- function(PCA_res,
                             group_vector,
                             c_values = c(1),
                             interval = 1,
                             cross_validate = FALSE,
                             use_inv_prob_weighting = FALSE,
                             use_SMOTE = FALSE,
                             return_all_fold_metrics = FALSE) {
  total_n_PCs <- length(PCA_res$sdev)
  group <- group_vector
  
  # Initialize empty list
  PCA_SVM_res_list <- list()
  
  # Start from 2 if using SMOTE
  starting_i <- ifelse(use_SMOTE, 2, 1)
  
  # Increasingly iterate over each PCs
  for (i in seq(starting_i, total_n_PCs, by = interval)) {
    svm_for_pc <- as.data.frame(cbind(group, PCA_res$x[, 1:i])) %>%
        mutate_at(vars(contains("V")), as.numeric) %>%
        mutate_at(vars(starts_with("PC")), as.numeric) 
    
    if (use_inv_prob_weighting) {
      sample_props <- svm_for_pc %>%
        dplyr::summarise(control_prop = sum(group=="Control") / n(),
                         schz_prop = sum(group=="Schz")/n())
      
      # Convert to sample weights based on inverse of probability
      sample_wts <- list("Control" = 1/sample_props$control_prop,
                         "Schz" = 1/sample_props$schz_prop)
    } else {
      sample_wts <- list("Control" = 1, "Schz" = 1)
    }
    
    # Compile results into a dataframe
    df_res <- k_fold_CV_linear_SVM(input_data = svm_for_pc,
                                   k = 10,
                                   c_values = c_values,
                                   svm_kernel = "linear",
                                   sample_wts = sample_wts,
                                   use_SMOTE = FALSE,
                                   shuffle_labels = FALSE,
                                   return_all_fold_metrics = return_all_fold_metrics)
  
    df_res$Num_PCs <- i
    
    # Append results to list
    PCA_SVM_res_list <- rlist::list.append(PCA_SVM_res_list, df_res)
  }
  
  PCA_SVM_res <- do.call(plyr::rbind.fill, PCA_SVM_res_list)
  return(PCA_SVM_res)
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

#-------------------------------------------------------------------------------
# Run pairwise PYSPI CV linear SVM by SPI
#-------------------------------------------------------------------------------
run_in_sample_pairwise_svm_by_SPI <- function(pairwise_data,
                                              svm_kernel = "linear",
                                              test_package = "e1071",
                                              noise_proc = "AROMA+2P",
                                              use_inv_prob_weighting = FALSE,
                                              use_SMOTE = FALSE,
                                              shuffle_labels = FALSE) {
  
  # Initialize results list for SPI-wise in-sample linear SVM
  class_res_list <- list()
  
  # Iterate over each PYSPI SPI
  for (this_SPI in unique(all_pyspi_data$SPI)) {
    # Identify whether current SPI is directed or undirected
    directionality <- SPI_directionality %>% filter(SPI == this_SPI) %>% pull(Direction)
    
    # Reshape SPI data to contain brain region names
    SPI_data <- subset(all_pyspi_data, SPI == this_SPI)  %>%
      mutate(comparison = row_number()) %>%
      pivot_longer(cols = c(brain_region_1,
                            brain_region_2),
                   names_to = "Region_Number",
                   values_to = "Index") %>%
      left_join(ROI_index) %>%
      dplyr::select(-Index) %>%
      pivot_wider(id_cols = c("Subject_ID", "group", "SPI", "value", "comparison"),
                  names_from = "Region_Number",
                  values_from = "ROI") %>%
      dplyr::select(-comparison) 
    
    # Combine brain regions into new pairwise column depending on directionality
    if (directionality == "Undirected") {
      SPI_data <- SPI_data %>%
        mutate(region_pair = ifelse(brain_region_1 < brain_region_2,
                                    paste0(brain_region_1, "_", brain_region_2),
                                    paste0(brain_region_2, "_", brain_region_1)))
    } else if (directionality == "Directed") {
      SPI_data <- SPI_data %>%
        mutate(region_pair = paste0(brain_region_1, "_", brain_region_2))
    }
    
    # Pivot data from long to wide for SVM
    data_for_SVM <- SPI_data %>%
      dplyr::select(Subject_ID, group, region_pair, value) %>%
      distinct() %>%
      pivot_wider(id_cols = c(Subject_ID, group),
                  names_from = region_pair, 
                  values_from = value) %>%
      dplyr::select(-Subject_ID) %>%
      drop_na()
    
    # Get sample weights if inverse probability weighting flag is applied
    if (use_inv_prob_weighting) {
      # Get control/schz proportions
      sample_props <- data_for_SVM %>%
        ungroup() %>%
        dplyr::summarise(control_prop = sum(group=="CONTROL") / n(),
                         schz_prop = sum(group=="SCHZ")/n())
      
      # Convert to sample weights based on inverse of probability
      sample_wts <- list("CONTROL" = 1/sample_props$control_prop,
                         "SCHZ" = 1/sample_props$schz_prop)
    } else {
      sample_wts <- list("CONTROL" = 1,
                         "SCHZ" = 1)
    }
    
    # Apply SMOTE if indicated
    if (use_SMOTE) {
      data_for_SVM <- smotefamily::SMOTE(data_for_SVM[,-1], data_for_SVM$group, K = 5)$data %>%
        dplyr::rename("group" = "class")
    }
    
    # Pass data_for_SVM to in_sample_linear_SVM
    if (nrow(data_for_SVM ) > 0) {
      SVM_results <- k(input_data = data_for_SVM,
                                          svm_kernel = svm_kernel,
                                          sample_wts = sample_wts,
                                          shuffle_labels = shuffle_labels) %>%
        mutate(SPI = this_SPI,
               Noise_Proc = noise_proc,
               use_inv_prob_weighting = use_inv_prob_weighting,
               use_SMOTE = use_SMOTE)
      
      # Append results to list
      class_res_list <- rlist::list.append(class_res_list, SVM_results)
    }
  }
  
  # Combine results from all regions into a dataframe
  class_res_df <- do.call(plyr::rbind.fill, class_res_list)
  
  # Return dataframe
  return(class_res_df)
}
