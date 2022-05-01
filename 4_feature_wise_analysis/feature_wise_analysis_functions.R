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
# Plot distribution of classifier accuracies across features 
# along with control vs schz proportions
#-------------------------------------------------------------------------------
plot_class_acc_w_props <- function(class_res,
                                   cv = FALSE,
                                   rdata_path,
                                   noise_procs = c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER"),
                                   ylab = "Number of catch22 features") {
  
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
    dplyr::select(catch22_Feature, Noise_Proc, accuracy, balanced_accuracy) %>%
    left_join(., ctrl_prop) %>%
    pivot_longer(cols=c(accuracy, balanced_accuracy),
                 names_to = "Metric",
                 values_to = "Value") %>%
    mutate(Metric = gsub("_", " ", Metric),
           ctrl_prop = ifelse(Metric == "balanced accuracy", NA_real_, ctrl_prop)) %>%
    mutate(Noise_Proc = factor(Noise_Proc, levels = noise_procs)) %>%
    ggplot(data=., mapping=aes(x=Value)) +
    geom_histogram(fill="lightsteelblue", bins=50) +
    geom_vline(aes(xintercept = ctrl_prop), linetype=2, color="gray30") +
    facet_grid(Noise_Proc ~ Metric, scales="free", switch="y") +
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
# Run simple in-sample multi-region linear SVM by catch22 feature
#-------------------------------------------------------------------------------
run_in_sample_svm_by_feature <- function(rdata_path,
                                          svm_kernel = "linear",
                                          test_package = "e1071",
                                          noise_procs = c("AROMA+2P", 
                                                          "AROMA+2P+GMR", 
                                                          "AROMA+2P+DiCER"),
                                          use_inv_prob_weighting = FALSE,
                                          upsample_minority = FALSE,
                                          downsample_majority = FALSE) {
  
  feature_wise_class_res_list <- list()
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
    } else {
      sample_wts <- list("Control" = 1,
                         "Schz" = 1)
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
    
    # Iterate over each brain region (ROI) for e1071 SVM
    for (this_feature in unique(feature_matrix$names)) {
      cat("\nNow running", test_package, "SVM for", this_feature, noise_proc, "\n")
      
      if (upsample_minority | downsample_majority) {
        # Subset region data and convert to wide format
        data_for_svm <- resampled_data %>%
          left_join(., feature_matrix) %>%
          filter(names == this_feature) %>%
          dplyr::group_by(Unique_ID, names) %>%
          dplyr::filter(!any(is.na(values))) %>%
          dplyr::ungroup() %>%
          dplyr::select(Unique_ID, group, Brain_Region, values) %>%
          tidyr::pivot_wider(id_cols = c(Unique_ID, group),
                             names_from = Brain_Region,
                             values_from 
                             = values) %>%
          dplyr::select(-Unique_ID)
      } else {
        data_for_svm <- feature_matrix %>%
          filter(names == this_feature) %>%
          dplyr::group_by(Subject_ID, names) %>%
          dplyr::filter(!any(is.na(values))) %>%
          dplyr::ungroup() %>%
          dplyr::select(Subject_ID, group, Brain_Region, values) %>%
          tidyr::pivot_wider(id_cols = c(Subject_ID, group),
                             names_from = Brain_Region,
                             values_from 
                             = values) %>%
          dplyr::select(-Subject_ID)
      }
      
      if (test_package == "e1071") {
        # Run SVM with supplied kernel type
        svmModel <- e1071::svm(factor(group) ~ .,
                               kernel = svm_kernel,
                               cost = 1,
                               data = data_for_svm,
                               class.weights = sample_wts)
      } else if (test_package=="kernlab") {
        # Run SVM with supplied kernel type
        svmModel <- kernlab::ksvm(factor(group) ~ .,
                                  kernel = svm_kernel,
                                  type = "C-svc",
                                  C = 1,
                                  prob.model = F,
                                  data = data_for_svm,
                                  class.weights = sample_wts)
      }

      
      # Generate in-sample predictions based on SVM model
      pred <- predict(svmModel, data_for_svm)
      data_for_svm$group <- factor(data_for_svm$group, levels = levels(pred))
      
      # Calculate accuracy and balanced accuracy
      accuracy <- sum(pred == data_for_svm$group)/length(pred)
      cm <- as.matrix(caret::confusionMatrix(pred, data_for_svm$group)$table)
      balanced_accuracy <- calculate_in_sample_balanced_accuracy(cm = cm)
      
      # Compile results into a dataframe
      feature_df_res <- data.frame(catch22_Feature = this_feature,
                                  Accuracy = accuracy,
                                  Balanced_Accuracy = balanced_accuracy,
                                  Noise_Proc = noise_proc,
                                  SVM_Package = test_package)
      
      # Append results to list
      feature_wise_class_res_list <- rlist::list.append(feature_wise_class_res_list,
                                                        feature_df_res)
    }
  }
  # Combine results from all regions into a dataframe
  feature_wise_class_res_df <- do.call(plyr::rbind.fill, feature_wise_class_res_list)
  
  # Return dataframe
  return(feature_wise_class_res_df)
}


#-------------------------------------------------------------------------------
# Run linear multi-ROI SVM in caret with kernlab per catch22 feature
#-------------------------------------------------------------------------------
run_caret_multi_SVM_by_feature <- function(rdata_path,
                                          use_inv_prob_weighting = FALSE,
                                          upsample_minority = FALSE,
                                          downsample_majority = FALSE,
                                          test_package = "e1071",
                                          noise_procs = c("AROMA+2P", 
                                                          "AROMA+2P+GMR", 
                                                          "AROMA+2P+DiCER")) {
  
  feature_wise_class_res_list <- list()
  for (noise_proc in noise_procs) {
    cat("\nNow on", noise_proc, "\n")
    # Clean up names
    noise_label <- gsub("\\+", "_", noise_proc)
    
    # Load catch22 data for current noise processing method
    feature_matrix <- readRDS(paste0(rdata_path, sprintf("UCLA_%s_catch22.Rds", 
                                                         noise_label)))   %>%
      dplyr::group_by(Subject_ID, Brain_Region) %>%
      dplyr::filter(!any(is.na(values))) %>%
      dplyr::ungroup()
    
    subj_groups <- feature_matrix %>%
      dplyr::group_by(Subject_ID, Brain_Region) %>%
      dplyr::filter(!any(is.na(values))) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(Subject_ID, group) 
    
    if (use_inv_prob_weighting) {
      # Get control/schz proportions
      sample_props <- subj_groups %>%
        dplyr::summarise(control_prop = sum(group=="Control") / n(),
                         schz_prop = sum(group=="Schz")/n())
      
      # Convert to sample weights based on inverse of probability
      model_wts <- ifelse(subj_groups$group == "Control", 
                           1/sample_props$control_prop,
                           1/sample_props$schz_prop)
      
      # Convert to sample weights based on inverse of probability
      sample_wts <- list("Control" = 1/sample_props$control_prop,
                         "Schz" = 1/sample_props$schz_prop)
    } else {
      model_wts <- rep(1, nrow(subj_groups))
      
      sample_wts <- list("Control" = 1, "Schz" = 1)
    }
    
    for (this_feature in unique(feature_matrix$names)) {
      cat("\nNow running caret", test_package, "SVM for", this_feature, noise_proc, "\n")
      
      feature_matrix_feature <- subset(feature_matrix, names==this_feature)
      
      data_for_svm <- feature_matrix_feature %>%
        dplyr::group_by(Subject_ID, names) %>%
        dplyr::filter(!any(is.na(values))) %>%
        dplyr::ungroup() %>%
        dplyr::select(Subject_ID, group, Brain_Region, values) %>%
        tidyr::pivot_wider(id_cols = c(Subject_ID, group),
                           names_from = Brain_Region,
                           values_from = values) %>%
        dplyr::select(-Subject_ID)
      
      flds <- createFolds(data_for_svm$group, k = 10, list = TRUE, returnTrain = FALSE)
      
      accuracy_list <- list()
      balanced_accuracy_list <- list()
      
      for (i in 1:length(flds)) {
        test_i <- flds[[i]]
        train_i <- setdiff(1:nrow(data_for_svm), test_i)
        
        test_data <- data_for_svm[test_i, ]
        train_data <- data_for_svm[train_i, ]
        
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
        cm <- as.matrix(caret::confusionMatrix(pred, test_data$group)$table)
        balanced_accuracy <- calculate_in_sample_balanced_accuracy(cm = cm)
        
        accuracy_list[[i]] <- accuracy
        balanced_accuracy_list[[i]] <- balanced_accuracy
      } 
      
      accuracy_avg <- mean(unlist(accuracy_list), na.rm=T)
      accuracy_sd <- sd(unlist(accuracy_list), na.rm=T)
      balanced_accuracy_avg <- mean(unlist(balanced_accuracy_list), na.rm=T)
      balanced_accuracy_sd <- sd(unlist(balanced_accuracy_list), na.rm=T)
      
      feature_res <- data.frame(catch22_Feature = this_feature,
                            Accuracy = accuracy_avg,
                            AccuracySD = accuracy_sd,
                            Balanced_Accuracy = balanced_accuracy_avg,
                            Balanced_AccuracySD = balanced_accuracy_sd,
                            Test_Package = test_package, 
                            Noise_Proc = noise_proc)
      
      
      # # Train SVM model
      # fitControl <- caret::trainControl(method = "cv",
      #                                   number = 10,
      #                                   savePredictions = "all",
      #                                   summaryFunction = calculate_balanced_accuracy,
      #                                   classProbs = TRUE)
      # 
      # # Define tuning grid
      # if (test_package == "e1071") {
      #   test_method <- "svmLinearWeights"
      #   grid <- expand.grid(cost = 1)
      # } else {
      #   test_method <- "svmLinear"
      #   grid <- expand.grid(C = 1)
      # }
      # 
      # 
      # mod <- caret::train(group ~ .,
      #                     data = data_for_svm,
      #                     method = test_method,
      #                     trControl = fitControl,
      #                     metric = "Balanced_Accuracy",
      #                     maximize = T,
      #                     weights = model_wts,
      #                     tuneGrid = grid,
      #                     preProcess = c("center", "scale", "nzv"))
      # 
      # 
      # # Generate predictions
      # preds <- mod$pred
      # data_for_svm$group <- factor(data_for_svm$group, levels = levels(preds$pred))
      # 
      # # Get accuracy + balanced accuracy results
      # feature_res <- mod$results %>%
      #   mutate(catch22_Feature = this_feature,
      #          Noise_Proc = noise_proc,
      #          Test_Package = test_package)
      
      # Append results to list
      feature_wise_class_res_list <- rlist::list.append(feature_wise_class_res_list,
                                                    feature_res)
    }
    
  }
  # Combine results from all regions into a dataframe
  feature_wise_class_res_df <- do.call(plyr::rbind.fill, feature_wise_class_res_list)
  
  # Return dataframe
  return(feature_wise_class_res_df)
}


#-------------------------------------------------------------------------------
# Calculate empirical p-values per feature based on null distribution
#-------------------------------------------------------------------------------

calc_empirical_nulls <- function(class_res,
                                 null_data) {
  main_res <- class_res %>%
    dplyr::select(catch22_Feature, Noise_Proc, accuracy, balanced_accuracy) %>%
    mutate(Type = "main")
  
  merged_list <- list()
  
  for (this_feature in unique(class_res$catch22_Feature)) {
    feature_null <- model_free_shuffle_null_upsampled %>% mutate(catch22_Feature = this_feature)
    feature_main <- subset(main_res, catch22_Feature == this_feature)
    
    feature_merged <- plyr::rbind.fill(feature_main,
                                      feature_null) %>%
      group_by(catch22_Feature, Noise_Proc) %>%
      dplyr::summarise(main_accuracy = unique(accuracy[Type=="main"]),
                       main_balanced_accuracy = unique(balanced_accuracy[Type=="main"]),
                       acc_p = 1 - (sum(main_accuracy > accuracy[Type=="null"], 
                                        na.rm=T)/n()),
                       bal_acc_p = 1 - (sum(main_balanced_accuracy >
                                              balanced_accuracy[Type=="null"])/n())) %>%
      ungroup() %>%
      distinct() %>%
      dplyr::rename("accuracy" = "main_accuracy",
                    "balanced_accuracy" = "main_balanced_accuracy") 
    
    
    merged_list <- rlist::list.append(merged_list, feature_merged)
  }
  main_p_values <- do.call(plyr::rbind.fill, merged_list) %>%
    ungroup() %>%
    mutate(acc_p_adj = p.adjust(acc_p, method="BH"),
           bal_acc_p_adj = p.adjust(bal_acc_p, method="BH")) 
  
  return(main_p_values)
}