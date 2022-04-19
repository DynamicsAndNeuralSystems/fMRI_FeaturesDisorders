#------------------------------------
# This script sets out to produce a
# function for reading in matlab time
# series files into R
#------------------------------------

#--------------------------------------
# Author: Annie Bryant, 15 April 2022
#--------------------------------------

require(rlist)
library(tidyverse)
library(theft)
library(caret)
library(broom)
library(kernlab)
library(yardstick)

#-------------------------------------------------------------------------------
# Plot distribution of classifier accuracies across features 
# along with control vs schz proportions
#-------------------------------------------------------------------------------
plot_class_acc_w_props <- function(class_res,
                                   cv = FALSE,
                                   rdata_path,
                                   noise_procs = c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER"),
                                   ylab = "Number of Feature + ROI Combos") {
  
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
      group_by(Subject_ID) %>%
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
    left_join(., ctrl_prop) %>%
    pivot_longer(cols=c(accuracy, balanced_accuracy),
                 names_to = "Metric",
                 values_to = "Value") %>%
    mutate(Metric = gsub("_", " ", Metric)) %>%
    mutate(Noise_Proc = factor(Noise_Proc, levels = noise_procs)) %>%
    ggplot(data=., mapping=aes(x=1, y=Value)) +
    geom_bar(stat="identity",fill="lightsteelblue", bins=50) +
    geom_hline(aes(yintercept = ctrl_prop), linetype=2, color="gray30") +
    facet_grid(Metric ~ Noise_Proc, scales="free", switch="y") +
    xlab(xlab) +
    ylab(ylab) +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle=0))
}

calculate_in_sample_balanced_accuracy <- function(cm) {
  recall <- 1:nrow(cm) %>%
    purrr::map(~ calculate_recall(cm, x = .x)) %>%
    unlist()
  
  balanced_accuracy <- sum(recall) / length(recall)
  
  # Return results
  
  out <- c(balanced_accuracy)
  names(out) <- c("Balanced_Accuracy")
  return(out)
}

#-------------------------------------------------------------------------------
# Run simple in-sample multi-region kernlab linear SVM by catch22 feature + ROI
#-------------------------------------------------------------------------------
run_in_sample_ksvm_by_feature_and_ROI <- function(rdata_path,
                                                  svm_kernel = "linear",
                                                  noise_procs = c("AROMA+2P", 
                                                                  "AROMA+2P+GMR", 
                                                                  "AROMA+2P+DiCER"),
                                                  use_inv_prob_weighting = FALSE,
                                                  upsample_minority = FALSE,
                                                  downsample_majority = FALSE) {
  
  feature_ROI_wise_class_res_list <- list()
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
    
    # Case when we upsample the minority class
    if (upsample_minority) {
      ctrl_subjects <- subset(feature_matrix, group=="Control") %>% 
        distinct(Subject_ID) %>% 
        pull(Subject_ID)
      
      schz_subjects_upsampled <- subset(feature_matrix, group=="Schz") %>%
        distinct(Subject_ID) %>%
        pull(Subject_ID) %>%
        sample(., length(ctrl_subjects), replace=T)
      
      resampled_data <- data.frame(Subject_ID = c(ctrl_subjects, 
                                                  schz_subjects_upsampled)) %>%
        mutate(Unique_ID = make.unique(Subject_ID))
    }
    
    # Case when we downsample the majority class
    if (downsample_majority) {
      schz_subjects <- subset(feature_matrix, group=="Schz") %>%
        distinct(Subject_ID) %>%
        pull(Subject_ID)
      
      ctrl_subjects_downsampled <- subset(feature_matrix, group=="Control") %>% 
        distinct(Subject_ID) %>% 
        pull(Subject_ID) %>%
        sample(., length(schz_subjects), replace=F)
      
      resampled_data <- data.frame(Subject_ID = c(schz_subjects, 
                                                  ctrl_subjects_downsampled)) %>%
        mutate(Unique_ID = make.unique(Subject_ID))
    }
    
    if (upsample_minority | downsample_majority) {
      # Subset region data and convert to wide format
      data_for_svm <- resampled_data %>%
        left_join(., feature_matrix) %>%
        dplyr::group_by(Unique_ID) %>%
        dplyr::filter(!any(is.na(values))) %>%
        dplyr::ungroup() %>%
        tidyr::unite("SVM_Feature", names, Brain_Region, sep="_") %>%
        dplyr::select(Unique_ID, group, SVM_Feature, values) %>%
        tidyr::pivot_wider(id_cols = c(Unique_ID, group),
                           names_from = SVM_Feature,
                           values_from 
                           = values) %>%
        dplyr::select(-Unique_ID)
    } else {
      data_for_svm <- feature_matrix %>%
        dplyr::group_by(Subject_ID) %>%
        dplyr::filter(!any(is.na(values))) %>%
        dplyr::ungroup() %>%
        tidyr::unite("SVM_Feature", names, Brain_Region, sep="_") %>%
        dplyr::select(Subject_ID, group, SVM_Feature, values) %>%
        tidyr::pivot_wider(id_cols = c(Subject_ID, group),
                           names_from = SVM_Feature,
                           values_from 
                           = values) %>%
        dplyr::select(-Subject_ID)
    }
    
    if (use_inv_prob_weighting) {
      # Run SVM with supplied kernel type
      svmModel <- kernlab::ksvm(factor(group) ~ .,
                                type = "C-svc",
                                kernel = "vanilladot",
                                data = data_for_svm,
                                class.weights = sample_wts,
                                prob.model=F)
    } else {
      # Run SVM with supplied kernel type
      svmModel <- kernlab::ksvm(factor(group) ~ .,
                                type = "C-svc",
                                kernel = "vanilladot",
                                data = data_for_svm,
                                prob.model=F)
    }

    
    # Generate in-sample predictions based on SVM model
    pred <- predict(svmModel, data_for_svm)
    data_for_svm$group <- factor(data_for_svm$group, levels = levels(pred))
    
    # Calculate accuracy and balanced accuracy
    accuracy <- sum(pred == data_for_svm$group)/length(pred)
    cm <- as.matrix(caret::confusionMatrix(pred, data_for_svm$group)$table)
    balanced_accuracy <- calculate_in_sample_balanced_accuracy(cm = cm)
    
    # Compile results into a dataframe
    feature_ROI_df_res <- data.frame(Accuracy = accuracy,
                                 Balanced_Accuracy = balanced_accuracy,
                                 Noise_Proc = noise_proc)
    
    # Append results to list
    feature_ROI_wise_class_res_list <- rlist::list.append(feature_ROI_wise_class_res_list,
                                                          feature_ROI_df_res)
  }
  # Combine results from all regions into a dataframe
  feature_ROI_wise_class_res_df <- do.call(plyr::rbind.fill, feature_ROI_wise_class_res_list)
  
  # Return dataframe
  return(feature_ROI_wise_class_res_df)
}

#-------------------------------------------------------------------------------
# Run linear multi-ROI SVM in caret with kernlab per catch22 feature + ROI
#-------------------------------------------------------------------------------
run_caret_multi_SVM_by_feature_and_ROI <- function(rdata_path,
                                                   use_inv_prob_weighting = FALSE,
                                                   upsample_minority = FALSE,
                                                   downsample_majority = FALSE,
                                                   noise_procs = c("AROMA+2P", 
                                                                   "AROMA+2P+GMR", 
                                                                   "AROMA+2P+DiCER")) {
  
  feature_ROI_wise_class_res_list <- list()
  for (noise_proc in noise_procs) {
    # Clean up names
    noise_label <- gsub("\\+", "_", noise_proc)
    
    # Load catch22 data for current noise processing method
    feature_matrix <- readRDS(paste0(rdata_path, sprintf("UCLA_%s_catch22.Rds", 
                                                         noise_label)))   %>%
      dplyr::group_by(Subject_ID, Brain_Region) %>%
      dplyr::filter(!any(is.na(values))) %>%
      dplyr::ungroup()
    
    if (use_inv_prob_weighting) { 
      # Get control/schz proportions
      subj_groups <- feature_matrix %>%
        dplyr::distinct(Subject_ID, group) 
      
      sample_props <- subj_groups %>%
        dplyr::summarise(control_prop = sum(group=="Control"),
                         schz_prop = sum(group=="Schz"))
      
      model_weights <- ifelse(subj_groups$group=="Control", 
                              (1/sample_props$control_prop)*0.5,
                              (1/sample_props$schz_prop)*0.5)
    }
    
    # Case when we upsample the minority class
    if (upsample_minority) {
      ctrl_subjects <- subset(feature_matrix, group=="Control") %>% 
        distinct(Subject_ID) %>% 
        pull(Subject_ID)
      
      schz_subjects_upsampled <- subset(feature_matrix, group=="Schz") %>%
        distinct(Subject_ID) %>%
        pull(Subject_ID) %>%
        sample(., length(ctrl_subjects), replace=T)
      
      resampled_data <- data.frame(Subject_ID = c(ctrl_subjects, 
                                                  schz_subjects_upsampled)) %>%
        mutate(Unique_ID = make.unique(Subject_ID))
    }
    
    # Case when we downsample the majority class
    if (downsample_majority) {
      schz_subjects <- subset(feature_matrix, group=="Schz") %>%
        distinct(Subject_ID) %>%
        pull(Subject_ID)
      
      ctrl_subjects_downsampled <- subset(feature_matrix, group=="Control") %>% 
        distinct(Subject_ID) %>% 
        pull(Subject_ID) %>%
        sample(., length(schz_subjects), replace=F)
      
      resampled_data <- data.frame(Subject_ID = c(schz_subjects, 
                                                  ctrl_subjects_downsampled)) %>%
        mutate(Unique_ID = make.unique(Subject_ID))
    }
    
    if (upsample_minority | downsample_majority) {
      # Subset region data and convert to wide format
      data_for_svm <- resampled_data %>%
        left_join(., feature_matrix) %>%
        dplyr::group_by(Unique_ID) %>%
        dplyr::filter(!any(is.na(values))) %>%
        dplyr::ungroup() %>%
        tidyr::unite("SVM_Feature", names, Brain_Region, sep="_") %>%
        dplyr::select(Unique_ID, group, SVM_Feature, values) %>%
        tidyr::pivot_wider(id_cols = c(Unique_ID, group),
                           names_from = SVM_Feature,
                           values_from 
                           = values) %>%
        dplyr::select(-Unique_ID)
    } else {
      data_for_svm <- feature_matrix %>%
        dplyr::group_by(Subject_ID) %>%
        dplyr::filter(!any(is.na(values))) %>%
        dplyr::ungroup() %>%
        tidyr::unite("SVM_Feature", names, Brain_Region, sep="_") %>%
        dplyr::select(Subject_ID, group, SVM_Feature, values) %>%
        tidyr::pivot_wider(id_cols = c(Subject_ID, group),
                           names_from = SVM_Feature,
                           values_from 
                           = values) %>%
        dplyr::select(-Subject_ID)
    }
    
    # Train SVM model
    fitControl <- caret::trainControl(method = "cv",
                                      number = 10,
                                      savePredictions = "all",
                                      summaryFunction = calculate_balanced_accuracy,
                                      classProbs = TRUE)
    
    if (use_inv_prob_weighting) {
      
      # Run kernlab linear SVM with caret
      mod <- caret::train(group ~ .,
                          data = data_for_svm,
                          method = "svmLinear",
                          trControl = fitControl,
                          metric = "Balanced_Accuracy",
                          maximize = T,
                          weights = model_weights,
                          preProcess = c("center", "scale", "nzv"))
    } else {
      mod <- caret::train(group ~ .,
                          data = data_for_svm,
                          method = "svmLinear",
                          trControl = fitControl,
                          metric = "Balanced_Accuracy",
                          maximize = T,
                          preProcess = c("center", "scale", "nzv"))
    }
    
    # Generate predictions
    preds <- mod$pred
    data_for_svm$group <- factor(data_for_svm$group, levels = levels(preds$pred))
    cm <- list(caret::confusionMatrix(preds$pred, data_for_svm$group)$table)
    
    # Get accuracy + balanced accuracy results
    feature_ROI_res <- mod$results %>%
      mutate(Test_Method = "linear_SVM",
             Noise_Proc = noise_proc,
             Confusion_Matrix = I(cm))
    
    # Append results to list
    feature_ROI_wise_class_res_list <- rlist::list.append(feature_ROI_wise_class_res_list,
                                                          feature_ROI_res)
  }
  
  
  # Combine results from all regions into a dataframe
  feature_ROI_wise_class_res_df <- do.call(plyr::rbind.fill, feature_ROI_wise_class_res_list)
  
  # Return dataframe
  return(feature_ROI_wise_class_res_df)
}

#-------------------------------------------------------------------------------
# Calculate empirical p-values per feature+ROI combo based on null distribution
#-------------------------------------------------------------------------------

calc_empirical_nulls <- function(class_res,
                                 null_data) {
  
  main_p_vals <- class_res %>%
    dplyr::select(Noise_Proc, accuracy, balanced_accuracy) %>%
    
    # Set main accuracy name
    mutate(Type = "main") %>%
    
    # Merge with null data
    plyr::rbind.fill(., null_data) %>%
    
    # Group by noise-processing method
    group_by(Noise_Proc) %>%
    
    # Calculate accuracy + balanced accuracy for real result, 
    # Along with empirical p-values based on rank relative to null distribution
    dplyr::summarise(main_accuracy = unique(accuracy[Type=="main"]),
                     main_balanced_accuracy = unique(balanced_accuracy[Type=="main"]),
                     acc_p = 1 - (sum(main_accuracy > accuracy[Type=="null"], 
                                      na.rm=T)/n()),
                     bal_acc_p = 1 - (sum(main_balanced_accuracy >
                                            balanced_accuracy[Type=="null"])/n())) %>%
    ungroup() %>%
    distinct() %>%
    dplyr::rename("accuracy" = "main_accuracy",
                  "balanced_accuracy" = "main_balanced_accuracy") %>%
    mutate(acc_p_adj = stats::p.adjust(acc_p, method="BH"),
           bal_acc_p_adj = stats::p.adjust(bal_acc_p, method="BH"))
  
  return(main_p_vals)
}