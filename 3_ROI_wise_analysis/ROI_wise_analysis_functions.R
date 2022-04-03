#------------------------------------
# This script sets out to produce a
# function for reading in matlab time
# series files into R
#------------------------------------

#--------------------------------------
# Author: Annie Bryant, 29 March 2022
#--------------------------------------

require(rlist)
library(tidyverse)
library(theft)
library(caret)
library(broom)

#-------------------------------------------------------------------------------
# Plot distribution of classifier accuracies across ROIs 
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
    dplyr::select(Brain_Region, Noise_Proc, Accuracy, Balanced_Accuracy) %>%
    left_join(., ctrl_prop) %>%
    pivot_longer(cols=c(Accuracy, Balanced_Accuracy),
                 names_to = "Metric",
                 values_to = "Value") %>%
    mutate(Metric = gsub("_", " ", Metric),
           ctrl_prop = ifelse(Metric == "Balanced Accuracy", 0.5, ctrl_prop)) %>%
    ggplot(data=., mapping=aes(x=Value)) +
    geom_histogram(fill="lightsteelblue", bins=50) +
    geom_vline(aes(xintercept = ctrl_prop), linetype=2, color="gray30") +
    facet_grid(Noise_Proc ~ Metric, scales="free", switch="y") +
    xlab(xlab) +
    ylab("Number of ROIs") +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle=0))
}

#-------------------------------------------------------------------------------
# Calculate balanced accuracy in caret
#-------------------------------------------------------------------------------

calculate_balanced_accuracy <- function(data, lev = NULL, model = NULL) {
  # calculate accuracy
  accuracy <- sum(data$pred == data$obs)/length(data$obs)
  
  # Calculate balanced accuracy
  data_cm <- as.data.frame(t(caret::confusionMatrix(data$pred, data$obs)$byClass))
  balanced_accuracy <- mean(data_cm$`Balanced Accuracy`, na.rm=T)
  
  out <- c(accuracy, balanced_accuracy)
  names(out) <- c("Accuracy", "Balanced_Accuracy")
  return(out)
}

#-------------------------------------------------------------------------------
# Run theft's multivariate classifier for a given brain region
#-------------------------------------------------------------------------------
run_theft_multivar_classifier <- function(rdata_path,
                                          test_method = "svmLinear",
                                          noise_proc = "AROMA+2P",
                                          use_balanced_accuracy = FALSE) {
  
  # Clean up names
  test_label <- gsub("-|\\.", "_", test_method)
  noise_label <- gsub("\\+", "_", noise_proc)
  
  # Load corresponding feature matrix
  feature_matrix <- readRDS(paste0(rdata_path, sprintf("UCLA_%s_catch22.Rds", 
                                                       noise_label)))
  
  # Iterate over each brain region in the feature matrix
  region_test_statistics_list <- list()
  for (this_ROI in unique(feature_matrix$Brain_Region)) {
    
    # Subset to given ROI
    feature_matrix_ROI <- subset(feature_matrix, Brain_Region==this_ROI)
    
    cat("\nRunning theft multivariable classifier for", this_ROI, noise_proc, "\n")
    
    multi_classifier_outputs <- fit_multivariable_classifier(feature_matrix_ROI,
                                                             id_var = "Subject_ID",
                                                             group_var = "group",
                                                             by_set = F,
                                                             test_method = test_method,
                                                             use_balanced_accuracy = use_balanced_accuracy,
                                                             use_k_fold = T,
                                                             num_folds = 10,
                                                             use_empirical_null = T,
                                                             null_testing_method = "null model fits",
                                                             p_value_method = "gaussian",
                                                             num_permutations = 10)
    
    # Compile results
    test_stat <- multi_classifier_outputs$TestStatistics %>%
      mutate(Brain_Region = this_ROI,
             Noise_Proc = noise_proc,
             Test_Method = test_method,
             Norm_Method = norm_method)
    
    # Append results to list
    region_test_statistics_list <- rlist::list.append(region_test_statistics_list, test_stat)
  }
  
  # Combine list of results into a dataframe
  region_test_statistics <- do.call(plyr::rbind.fill, region_test_statistics_list)
  
  # Return dataframe
  return(region_test_statistics)
}

#-------------------------------------------------------------------------------
# Run simple e1071 in-sample SVM by brain region
#-------------------------------------------------------------------------------
run_e1071_SVM_by_region <- function(rdata_path,
                                    svm_kernel = "linear",
                                    noise_procs = c("AROMA+2P", 
                                                    "AROMA+2P+GMR", 
                                                    "AROMA+2P+DiCER")) {
  
  ROI_wise_class_res_list <- list()
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
    
    # Iterate over each brain region (ROI) for e1071 SVM
    for (this_ROI in unique(feature_matrix$Brain_Region)) {
      
      # Subset region data and convert to wide format
      region_data <- subset(feature_matrix, Brain_Region==this_ROI) %>%
        dplyr::group_by(Subject_ID, Brain_Region) %>%
        dplyr::filter(!any(is.na(values))) %>%
        dplyr::ungroup() %>%
        dplyr::select(Subject_ID, group, names, values) %>%
        tidyr::pivot_wider(id_cols = c(Subject_ID, group),
                           names_from = names,
                           values_from 
                           = values) %>%
        dplyr::select(-Subject_ID)
      
      # Run SVM with supplied kernel type
      svmModel <- e1071::svm(factor(group) ~ .,
                             kernel = svm_kernel,
                             data = region_data, 
                             class.weights = sample_wts)
      
      # Generate in-sample predictions based on SVM model
      pred <- stats::predict(svmModel)
      region_data$group <- factor(region_data$group, levels = levels(pred))
      
      # Calculate accuracy and balanced accuracy
      accuracy <- sum(pred == region_data$group)/length(pred)
      balanced_accuracy <- caret::confusionMatrix(pred, region_data$group)$byClass[["Balanced Accuracy"]]
      
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
# Run linear SVM in caret with e1071 per brain region
#-------------------------------------------------------------------------------
run_caret_e1071_SVM_by_region <- function(rdata_path,
                                          use_inv_prob_weighting = TRUE,
                                          noise_procs = c("AROMA+2P", 
                                                          "AROMA+2P+GMR", 
                                                          "AROMA+2P+DiCER")) {
  
  ROI_wise_class_res_list <- list()
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
    
    for (this_ROI in unique(feature_matrix$Brain_Region)) {
      cat("\nNow running linear SVM for", this_ROI, noise_label, "\n")
      
      feature_matrix_ROI <- subset(feature_matrix, Brain_Region==this_ROI)
      
      # Prep catch22 feature data for SVM
      data_for_svm <- feature_matrix_ROI %>%
        dplyr::rename("id" = "Subject_ID") %>%
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
      
      # Create tune grid
      e1071_grid <- expand.grid(cost = 1)
      
      # Train SVM model
      fitControl <- caret::trainControl(method = "cv",
                                        number = 10,
                                        savePredictions = "all",
                                        summaryFunction = calculate_balanced_accuracy,
                                        classProbs = TRUE)
      
      # Run e1071 SVM with caret
      mod <- caret::train(group ~ .,
                          data = data_for_svm,
                          method = "svmLinear2",
                          trControl = fitControl,
                          metric = "Balanced_Accuracy",
                          maximize = T,
                          weights = model_weights,
                          tuneGrid = e1071_grid,
                          preProcess = c("center", "scale", "nzv"))
      
      # Generate predictions
      preds <- mod$pred
      data_for_svm$group <- factor(data_for_svm$group, levels = levels(preds$pred))
      cm <- list(caret::confusionMatrix(preds$pred, data_for_svm$group)$table)
      
      # Get accuracy + balanced accuracy results
      ROI_res <- mod$results %>%
        mutate(Brain_Region = this_ROI,
               Test_Method = "linear_SVM",
               Noise_Proc = noise_proc,
               Confusion_Matrix = I(cm))
      
      # Append results to list
      ROI_wise_class_res_list <- rlist::list.append(ROI_wise_class_res_list,
                                                    ROI_res)
    }
    
  }
  # Combine results from all regions into a dataframe
  ROI_wise_class_res_df <- do.call(plyr::rbind.fill, ROI_wise_class_res_list)
  
  # Return dataframe
  return(ROI_wise_class_res_df)
}

#-------------------------------------------------------------------------------
# Run linear SVM in caret with e1071 per brain region with n=1000 shuffles
#-------------------------------------------------------------------------------
run_caret_e1071_SVM_by_region_reshuffle <- function(rdata_path,
                                                    num_shuffles = 1000,
                                                    use_inv_prob_weighting = TRUE,
                                                    noise_procs = c("AROMA+2P", 
                                                                    "AROMA+2P+GMR", 
                                                                    "AROMA+2P+DiCER")) {
  
  ROI_wise_class_res_list <- list()
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
    
    for (this_ROI in unique(feature_matrix$Brain_Region)) {
      cat("\nNow running linear SVM for", this_ROI, noise_label, "\n")
      
      feature_matrix_ROI <- subset(feature_matrix, Brain_Region==this_ROI)
      
      # Prep catch22 feature data for SVM
      data_for_svm <- feature_matrix_ROI %>%
        dplyr::rename("id" = "Subject_ID") %>%
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
      
      for (i in 1:num_shuffles) {
        data_for_svm_reshuffle <- data_for_svm
        data_for_svm_reshuffle$group_reshuffled <- sample(data_for_svm_reshuffle$group)
        data_for_svm_reshuffle <- data_for_svm_reshuffle %>%
          dplyr::select(-group)
        
        # Calculate model weights as inverse probability
        model_weights <- ifelse(data_for_svm_reshuffle$group_reshuffled == "Control", 
                                1/sample_props$control_prop,
                                1/sample_props$schz_prop)
        
        # Create tune grid
        e1071_grid <- expand.grid(cost = 1)
        
        # Train SVM model
        fitControl <- caret::trainControl(method = "cv",
                                          number = 10,
                                          savePredictions = "all",
                                          summaryFunction = calculate_balanced_accuracy,
                                          classProbs = TRUE)
        
        # Run e1071 SVM with caret
        mod <- caret::train(group_reshuffled ~ .,
                            data = data_for_svm_reshuffle,
                            method = "svmLinear2",
                            trControl = fitControl,
                            metric = "Balanced_Accuracy",
                            maximize = T,
                            weights = model_weights,
                            tuneGrid = e1071_grid,
                            preProcess = c("center", "scale", "nzv"))
        
        # Generate predictions
        preds <- mod$pred
        data_for_svm_reshuffle$group_reshuffled <- factor(data_for_svm_reshuffle$group_reshuffled, 
                                                          levels = levels(preds$pred))
        cm <- list(caret::confusionMatrix(preds$pred, data_for_svm_reshuffle$group_reshuffled)$table)
        
        # Get accuracy + balanced accuracy results
        ROI_res <- mod$results %>%
          mutate(Brain_Region = this_ROI,
                 Test_Method = "linear_SVM",
                 Null_Iteration = i,
                 Noise_Proc = noise_proc,
                 Confusion_Matrix = I(cm))
        
        # Append results to list
        ROI_wise_class_res_list <- rlist::list.append(ROI_wise_class_res_list,
                                                      ROI_res)
      }
    }
  }
  
  
  # Combine results from all regions into a dataframe
  ROI_wise_class_res_df <- do.call(plyr::rbind.fill, ROI_wise_class_res_list)
  
  # Return dataframe
  return(ROI_wise_class_res_df)
}
