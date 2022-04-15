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
library(foreach)
library(doParallel)
library(yardstick)

#-------------------------------------------------------------------------------
# Plot distribution of classifier accuracies across ROIs 
# along with control vs schz proportions
#-------------------------------------------------------------------------------
plot_class_acc_w_props <- function(class_res,
                                   cv = TRUE,
                                   rdata_path,
                                   noise_procs = c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER"),
                                   ylab = "Number of ROIs") {
  
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
    dplyr::select(Brain_Region, Noise_Proc, accuracy, balanced_accuracy) %>%
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
# Run model-free shuffles for a given distribution of schizophrenia vs controls
#-------------------------------------------------------------------------------

# Helper function to calculate accuracy + balanced accuracy per shuffle
calc_acc_bacc_for_shuffle <- function(x) {
  
  x <- factor(x, levels = unique(x))
  y <- sample(x, replace=F)
  
  # Calculate raw accuracy
  acc <- yardstick::accuracy_vec(truth = x,
                                 estimate = y)
  
  # Calculate balanced accuracy
  bacc <- yardstick::bal_accuracy_vec(truth = x,
                                      estimate = y)
  
  output <- data.frame(accuracy = acc,
                       balanced_accuracy = bacc)
  
  return(output)
}

# Run given number of model free class label shuffles and calculate accuracy
run_model_free_n_shuffles <- function(rdata_path,
                                      noise_procs = c("AROMA+2P", 
                                                      "AROMA+2P+GMR", 
                                                      "AROMA+2P+DiCER"),
                                      num_shuffles = 1000000,
                                      use_upsampling = FALSE) {
  
  set.seed(123)
  
  # Initialize empty list
  null_list <- list()
  
  # Iterate over each noise-processing method
  for (noise_proc in noise_procs) {
    # Clean up names
    noise_label <- gsub("\\+", "_", noise_proc)
    
    # Load catch22 data for current noise processing method
    feature_matrix <- readRDS(paste0(rdata_path, sprintf("UCLA_%s_catch22.Rds", 
                                                         noise_label)))   
    
    # Pull out the vector of diagnoses
    if (use_upsampling) {
      ctrl_subjects <- subset(feature_matrix, group=="Control") %>% 
        distinct(Subject_ID) %>% 
        pull(Subject_ID)
      
      schz_subjects_upsampled <- subset(feature_matrix, group=="Schz") %>%
        distinct(Subject_ID) %>%
        pull(Subject_ID) %>%
        sample(., length(ctrl_subjects), replace=T)
      
      input_groups <- data.frame(Subject_ID = c(ctrl_subjects, 
                                                  schz_subjects_upsampled)) %>%
        left_join(., feature_matrix %>% distinct(Subject_ID, group)) %>%
        pull(group)
    } else {
      input_groups <- feature_matrix %>% 
        distinct(Subject_ID, group) %>%
        mutate(group = factor(group, levels = c("Schz", "Control"))) %>%
        pull(group)
    }
    
    # Run model-free shuffles
    output <- 1:num_shuffles %>%
      purrr::map_df(~ calc_acc_bacc_for_shuffle(input_groups)) %>%
      mutate(Noise_Proc = noise_proc)
    
    null_list <- rlist::list.append(null_list, output)
  }
  region_wise_null_res <- do.call(plyr::rbind.fill, null_list) %>%
    dplyr::mutate(Type = "null")
  
  return(region_wise_null_res)
  
}

# Helper function to calculate empirical p-values based on null distribution
calc_empirical_nulls <- function(class_res,
                                 null_data,
                                 by_feature = F) {
  merged_list <- list()
  if (by_feature) {
    main_res <- class_res %>%
      dplyr::select(Brain_Region, Feature, Noise_Proc, accuracy, balanced_accuracy) %>%
      mutate(Type = "main")
  } else {
    main_res <- class_res %>%
      dplyr::select(Brain_Region, Noise_Proc, accuracy, balanced_accuracy) %>%
      mutate(Type = "main")
  }
  
  for (brain_region in unique(class_res$Brain_Region)) {
    region_null <- model_free_shuffle_null_upsampled %>% mutate(Brain_Region = brain_region)
    region_main <- subset(main_res, Brain_Region == brain_region)
    
    if (by_feature) {
      feature_list <- list()
      for (feature in unique(region_main$Feature)) {
        feature_null <- region_null %>%
          mutate(Feature = feature)
        feature_main <- subset(region_main, Feature == feature)
        feature_merged <- plyr::rbind.fill(feature_main,
                                          feature_null) %>%
          group_by(Brain_Region, Noise_Proc, Feature) %>%
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
        feature_list <- rlist::list.append(feature_list, feature_merged)
      }
      region_merged <- do.call(plyr::rbind.fill, feature_list)
      
    } else {
      region_merged <- plyr::rbind.fill(region_main,
                                        region_null) %>%
        group_by(Brain_Region, Noise_Proc) %>%
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
    }
    
    merged_list <- rlist::list.append(merged_list, region_merged)
  }
  main_p_values <- do.call(plyr::rbind.fill, merged_list) %>%
    ungroup() %>%
    mutate(acc_p_adj = p.adjust(acc_p, method="BH"),
           bal_acc_p_adj = p.adjust(bal_acc_p, method="BH")) 
  
  return(main_p_values)
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
                                         upsample_minority = FALSE,
                                         downsample_majority = FALSE) {
  
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
    for (this_ROI in unique(feature_matrix$Brain_Region)) {
      
      if (upsample_minority | downsample_majority) {
        # Subset region data and convert to wide format
        data_for_svm <- resampled_data %>%
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
    
    tryCatch({
      multi_classifier_outputs <- fit_multi_feature_classifier(feature_matrix_ROI,
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
      test_stat <- multi_classifier_outputs$TestStatistics %>%
        mutate(Brain_Region = this_ROI,
               Noise_Proc = noise_proc,
               Test_Method = test_method,
               Norm_Method = norm_method)
      
      # Append results to list
      region_test_statistics_list <- rlist::list.append(region_test_statistics_list, test_stat)
    }, error = function(e) {
      cat("\ntheft did not work for", this_ROI, "\n")
    })
  }
  
  # Combine list of results into a dataframe
  region_test_statistics <- do.call(plyr::rbind.fill, region_test_statistics_list)
  
  # Return dataframe
  return(region_test_statistics)
}


#-------------------------------------------------------------------------------
# Run linear multi-feature SVM in caret with kernlab ksvm per brain region
#-------------------------------------------------------------------------------
run_caret_multi_SVM_by_region <- function(rdata_path,
                                          use_inv_prob_weighting = FALSE,
                                          upsample_minority = FALSE,
                                          downsample_majority = FALSE,
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
    
    if (use_inv_prob_weighting) { 
      # Get control/schz proportions
      subj_groups <- feature_matrix %>%
        dplyr::group_by(Subject_ID, Brain_Region) %>%
        dplyr::filter(!any(is.na(values))) %>%
        dplyr::ungroup() %>%
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
    
    for (this_ROI in unique(feature_matrix$Brain_Region)) {
      cat("\nNow running linear SVM for", this_ROI, noise_label, "\n")
      
      if (upsample_minority | downsample_majority) {
        # Subset region data and convert to wide format
        data_for_svm <- resampled_data %>%
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
        data_for_svm <- feature_matrix %>%
          filter(Brain_Region == this_ROI) %>%
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
# Run linear single-feature SVM in-sample with kernlab per brain region
#-------------------------------------------------------------------------------
run_in_sample_uni_feature_linear_SVM_by_region <- function(rdata_path,
                                                          svm_kernel = "vanilladot",
                                                          noise_procs = c("AROMA+2P", 
                                                                          "AROMA+2P+GMR", 
                                                                          "AROMA+2P+DiCER"),
                                                          use_inv_prob_weighting = FALSE,
                                                          use_resampling = FALSE) {
  
  res_list <- list()
  for (noise_proc in noise_procs) {
    # Clean up names
    noise_label <- gsub("\\+", "_", noise_proc)
    
    # Load catch22 data for current noise processing method
    feature_matrix <- readRDS(paste0(rdata_path, sprintf("UCLA_%s_catch22.Rds", 
                                                         noise_label)))     
    
    if (use_resampling) {
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
    } else if (use_inv_prob_weighting) {
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
    
    # Iterate over each brain region (ROI) for e1071 SVM
    for (this_ROI in unique(feature_matrix$Brain_Region)) {
      
      if (use_resampling) {
        # Subset region data and convert to wide format
        region_data <- upsampled_data %>%
          left_join(., feature_matrix) %>%
          filter(Brain_Region == this_ROI)
      } else {
        region_data <- feature_matrix %>%
          filter(Brain_Region == this_ROI)
      }
      
      for (feature in unique(region_data$names)) {
        # Subset region data and convert to wide format
        data_for_svm <- region_data %>%
          filter(names == feature) %>%
          dplyr::group_by(Unique_ID, Subject_ID, Brain_Region) %>%
          dplyr::filter(!any(is.na(values))) %>%
          dplyr::ungroup() %>%
          dplyr::select(Unique_ID, Subject_ID, group, names, values) %>%
          tidyr::pivot_wider(id_cols = c(Unique_ID, Subject_ID, group),
                             names_from = names,
                             values_from 
                             = values) %>%
          dplyr::select(-Subject_ID, -Unique_ID)
        
        # Run SVM with supplied kernel type
        if (use_inv_prob_weighting) {
          svmModel <- kernlab::ksvm(factor(group) ~ .,
                                    type = "C-svc",
                                    kernel = svm_kernel,
                                    data = data_for_svm,
                                    prob.model=F, 
                                    class.weights = sample_wts)
        } else {
          svmModel <- kernlab::ksvm(factor(group) ~ .,
                                    type = "C-svc",
                                    kernel = svm_kernel,
                                    data = data_for_svm,
                                    prob.model=F)
        }

        
        # Generate in-sample predictions based on SVM model
        pred <- predict(svmModel)
        data_for_svm$group <- factor(data_for_svm$group, levels = levels(pred))
        
        # Calculate accuracy and balanced accuracy
        accuracy <- sum(pred == data_for_svm$group)/length(pred)
        cm <- as.matrix(caret::confusionMatrix(pred, data_for_svm$group)$table)
        balanced_accuracy <- calculate_in_sample_balanced_accuracy(cm = cm)
        
        # Compile results into a dataframe
        region_df_res <- data.frame(Brain_Region = this_ROI,
                                    Feature = feature,
                                    accuracy = accuracy,
                                    balanced_accuracy = balanced_accuracy,
                                    Noise_Proc = noise_proc)
        
        # Append results to list
        res_list <- rlist::list.append(res_list, region_df_res)
      }
      
    }
  }
  # Combine results from all regions into a dataframe
  res_df <- do.call(plyr::rbind.fill, res_list)
  
  # Return dataframe
  return(res_df)
}

#-------------------------------------------------------------------------------
# Run linear single-feature SVM in caret with kernlab per brain region
#-------------------------------------------------------------------------------
run_CV_uni_feature_linear_SVM_by_region <- function(rdata_path,
                                                    use_inv_prob_weighting = FALSE,
                                                    use_resampling = FALSE,
                                                    noise_procs = c("AROMA+2P", 
                                                                    "AROMA+2P+GMR", 
                                                                    "AROMA+2P+DiCER")) {
  
  res_list <- list()
  for (noise_proc in noise_procs) {
    # Clean up names
    noise_label <- gsub("\\+", "_", noise_proc)
    
    # Load catch22 data for current noise processing method
    feature_matrix <- readRDS(paste0(rdata_path, sprintf("UCLA_%s_catch22.Rds", 
                                                         noise_label)))   
    if (use_resampling) {
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
    } else if (use_inv_prob_weighting) {
      # Get control/schz proportions
      subj_groups <- feature_matrix %>%
        dplyr::group_by(Subject_ID, Brain_Region) %>%
        dplyr::filter(!any(is.na(values))) %>%
        dplyr::ungroup() %>%
        dplyr::distinct(Subject_ID, group) 
      
      sample_props <- subj_groups %>%
        dplyr::summarise(control_prop = sum(group=="Control"),
                         schz_prop = sum(group=="Schz"))
      
      model_weights <- ifelse(subj_groups$group=="Control", 
                              (1/sample_props$control_prop)*0.5,
                              (1/sample_props$schz_prop)*0.5)
    }
    
    for (this_ROI in unique(feature_matrix$Brain_Region)) {
      cat("\nNow running linear SVM for", this_ROI, noise_label, "\n")
      
      if (use_resampling) {
        # Subset region data and convert to wide format
        region_data <- upsampled_data %>%
          left_join(., feature_matrix) %>%
          filter(Brain_Region == this_ROI)
      } else {
        region_data <- feature_matrix %>%
          filter(Brain_Region == this_ROI)
      }
      
      # Train SVM model
      fitControl <- caret::trainControl(method = "cv",
                                        number = 10,
                                        savePredictions = "all",
                                        summaryFunction = calculate_balanced_accuracy,
                                        classProbs = TRUE)
      
      # Iterate over each feature
      for (feature in unique(region_data$names)) {
        
        # Prep catch22 feature data for SVM
        # Subset region data and convert to wide format
        data_for_svm <- region_data %>%
          filter(names == feature) %>%
          dplyr::group_by(Unique_ID, Subject_ID, Brain_Region) %>%
          dplyr::filter(!any(is.na(values))) %>%
          dplyr::ungroup() %>%
          dplyr::select(Unique_ID, Subject_ID, group, names, values) %>%
          tidyr::pivot_wider(id_cols = c(Unique_ID, Subject_ID, group),
                             names_from = names,
                             values_from 
                             = values) %>%
          dplyr::select(-Subject_ID, -Unique_ID)

        tryCatch({
          # Run SVM with caret
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
          ROI_feature_res <- mod$results %>%
            mutate(Brain_Region = this_ROI,
                   Test_Method = "linear_SVM",
                   Feature = feature,
                   Noise_Proc = noise_proc,
                   Confusion_Matrix = I(cm))
          
          # Append results to list
          res_list <- rlist::list.append(res_list, ROI_feature_res)
        }, error = function(e) {
          cat("\nWarning:", this_ROI, feature, noise_proc, "did not work.\n")
        })
        
      }
    }
    
  }
  # Combine results from all regions into a dataframe
  res <- do.call(plyr::rbind.fill, res_list)
  
  # Return dataframe
  return(res)
}