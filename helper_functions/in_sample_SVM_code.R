#------------------------------------
# This script sets out to produce a
# function for reading in matlab time
# series files into R
#------------------------------------

#--------------------------------------
# Author: Annie Bryant, 24 April 2022
#--------------------------------------

require(rlist)
library(tidyverse)
library(e1071)
library(kernlab)
library(cowplot)
library(caret)
theme_set(theme_cowplot())

#-------------------------------------------------------------------------------
# Plot distribution of classifier accuracies across given feature variable
# along with control vs schz proportions
#-------------------------------------------------------------------------------
plot_class_acc_w_props <- function(class_res,
                                   rdata_path,
                                   group_var = NULL,
                                   noise_procs = c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER"),
                                   cv = FALSE,
                                   ylab = "Number of ROIs") {
  
  # Calculate the proportion of control subjects after omitting NAs for each method
  ctrl_prop_list <- list()
  for (noise_proc in noise_procs) {
    # Clean up names
    noise_label <- gsub("\\+", "_", noise_proc)
    
    # Load corresponding feature matrix
    feature_matrix <- readRDS(paste0(rdata_path, sprintf("UCLA_%s_catch22.Rds", 
                                                         noise_label)))
    
    group_var_name <- ifelse(group_var == "Feature", "names", group_var)
    
    # Calculate proportion of controls after dropping NA
    if (!is.null(group_var)) {
      ctrl_proportion <- feature_matrix %>%
        group_by(Subject_ID, get(group_var_name)) %>%
        filter(!any(is.na(values))) %>%
        ungroup() %>%
        distinct(Subject_ID, group) %>%
        summarise(ctrl_prop = sum(group=="Control") / n()) %>%
        mutate(Noise_Proc = noise_proc)
    } else {
      ctrl_proportion <- feature_matrix %>%
        group_by(Subject_ID) %>%
        filter(!any(is.na(values))) %>%
        ungroup() %>%
        distinct(Subject_ID, group) %>%
        summarise(ctrl_prop = sum(group=="Control") / n()) %>%
        mutate(Noise_Proc = noise_proc)
    }
    
    
    ctrl_prop_list <- rlist::list.append(ctrl_prop_list, ctrl_proportion)
  }
  ctrl_prop <- do.call(plyr::rbind.fill, ctrl_prop_list)
  
  # Define x-axis label
  xlab <- ifelse(cv, "Mean statistic over 10-fold CV", "In-sample statistic")
  
  # Plot accuracy + balanced accuracy in histograms
  # Control subject proportion is highlighted for accuracy, 0.5 is highlighted for balanced accuracy
  if (!is.null(group_var)) {
    p_data <- class_res %>%
      dplyr::select(grouping_var, Noise_Proc, accuracy, balanced_accuracy)
  } else {
    p_data <- class_res %>%
      dplyr::select(Noise_Proc, accuracy, balanced_accuracy)
  }
  
  p_data %>%
    left_join(., ctrl_prop) %>%
    pivot_longer(cols=c(accuracy, balanced_accuracy),
                 names_to = "Metric",
                 values_to = "Value") %>%
    mutate(Metric = gsub("_", " ", Metric),
           ctrl_prop = ifelse(Metric == "balanced accuracy", 0.5, ctrl_prop)) %>%
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
                                           use_inv_prob_weighting = FALSE) {
  
  class_res_list <- list()
  for (noise_proc in noise_procs) {
    # Clean up names
    noise_label <- gsub("\\+", "_", noise_proc)
    
    # Load catch22 data for current noise processing method
    feature_matrix <- readRDS(paste0(rdata_path, sprintf("UCLA_%s_catch22.Rds", 
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
      cat("\nNow running linear SVM with", test_package, "for", group_var, "\n")
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
      
      
      # Run linear SVM with the supplied kernel type
      if (test_package == "e1071") {
        svmModel <- e1071::svm(factor(group) ~ .,
                               kernel = svm_kernel,
                               cost = 1,
                               data = data_for_SVM,
                               class.weights = sample_wts)
      } else if (test_package=="kernlab") {
        svmModel <- kernlab::ksvm(factor(group) ~ .,
                                  kernel = svm_kernel,
                                  type = "C-svc",
                                  C = 1,
                                  prob.model = F,
                                  data = data_for_SVM,
                                  class.weights = sample_wts)
      }
      
      
      # Generate in-sample predictions based on SVM model
      pred <- predict(svmModel, data_for_SVM)
      data_for_SVM$group <- factor(data_for_SVM$group, levels = levels(pred))
      
      # Calculate accuracy and balanced accuracy
      accuracy <- sum(pred == data_for_SVM$group)/length(pred)
      balanced_accuracy <- caret::confusionMatrix(data=pred, 
                                                  reference=data_for_SVM$group)$byClass[["Balanced Accuracy"]]
      balanced_accuracy <- ifelse(is.na(balanced_accuracy), 0.5, balanced_accuracy)
      
      # Compile results into a dataframe
      df_res <- data.frame(grouping_var = group_var,
                           Accuracy = accuracy,
                           Balanced_Accuracy = balanced_accuracy,
                           Noise_Proc = noise_proc,
                           SVM_Package = test_package)
      
      
      # Append results to list
      class_res_list <- rlist::list.append(class_res_list,
                                           df_res)
      
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
                                           seed = 127) {
  
  class_res_list <- list()
  for (noise_proc in noise_procs) {
    # Clean up names
    noise_label <- gsub("\\+", "_", noise_proc)
    
    # Load catch22 data for current noise processing method
    feature_matrix <- readRDS(paste0(rdata_path, sprintf("UCLA_%s_catch22.Rds", 
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
      cat("\nNow running linear SVM with", test_package, "for", group_var, "\n")
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
      
      set.seed(seed)
      flds <- createFolds(data_for_SVM$group, k = 10, list = TRUE, returnTrain = FALSE)
      
      accuracy_list <- list()
      balanced_accuracy_list <- list()
      
      
      for (i in 1:length(flds)) {
        test_i <- flds[[i]]
        train_i <- setdiff(1:nrow(data_for_SVM), test_i)
        
        test_data <- data_for_SVM[test_i, ]
        train_data <- data_for_SVM[train_i, ]
        
        if (test_package == "e1071") {
          # Run SVM with supplied kernel type
          svmModel <- e1071::svm(factor(group) ~ .,
                                 kernel = svm_kernel,
                                 cost = 1,
                                 data = train_data,
                                 class.weights = sample_wts)
        } else if (test_package=="kernlab") {
          # Run SVM with supplied kernel type
          svmModel <- kernlab::ksvm(factor(group) ~ .,
                                    kernel = svm_kernel,
                                    type = "C-svc",
                                    C = 1,
                                    prob.model = F,
                                    data = train_data,
                                    class.weights = sample_wts)
        }
        
        # Generate out-of-sample predictions based on SVM model
        pred <- predict(svmModel, test_data)
        test_data$group <- factor(test_data$group, levels = levels(pred))
        
        # Calculate accuracy and balanced accuracy
        accuracy <- sum(pred == test_data$group)/length(pred)
        balanced_accuracy <- caret::confusionMatrix(reference=test_data$group, 
                                                    data=pred)$byClass[["Balanced Accuracy"]]
        
        accuracy_list[[i]] <- accuracy
        balanced_accuracy_list[[i]] <- balanced_accuracy
      } 
      
      accuracy_avg <- mean(unlist(accuracy_list), na.rm=T)
      accuracy_sd <- sd(unlist(accuracy_list), na.rm=T)
      balanced_accuracy_avg <- mean(unlist(balanced_accuracy_list), na.rm=T)
      balanced_accuracy_sd <- sd(unlist(balanced_accuracy_list), na.rm=T)
      
      # Compile results into a dataframe
      df_res <- data.frame(grouping_var = group_var,
                           Accuracy = accuracy_avg,
                           AccuracySD = accuracy_sd,
                           Balanced_Accuracy = balanced_accuracy_avg,
                           Balanced_AccuracySD = balanced_accuracy_sd,
                           Noise_Proc = noise_proc,
                           SVM_Package = test_package)
      
      
      # Append results to list
      class_res_list <- rlist::list.append(class_res_list,
                                           df_res)
    }
  }
  
  # Combine results from all regions into a dataframe
  class_res_df <- do.call(plyr::rbind.fill, class_res_list)
  
  # Return dataframe
  return(class_res_df)
}

#-------------------------------------------------------------------------------
# Run model-free shuffles for a given distribution of schizophrenia vs controls
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
    feature_matrix <- readRDS(paste0(rdata_path, sprintf("UCLA_%s_catch22.Rds", 
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

# Helper function to calculate empirical p-values based on null distribution
calc_empirical_nulls <- function(class_res,
                                 null_data,
                                 grouping_var = "Brain_Region") {
  merged_list <- list()
  main_res <- class_res %>%
    dplyr::select(grouping_var, Noise_Proc, accuracy, balanced_accuracy) %>%
    mutate(Type = "main")
  
  for (group_var in unique(class_res$grouping_var)) {
    group_null <- model_free_shuffle_null_res %>% mutate(grouping_var = group_var)
    group_main <- subset(main_res, grouping_var == group_var)

    group_merged <- plyr::rbind.fill(group_main,
                                      group_null) %>%
      group_by(grouping_var, Noise_Proc) %>%
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

    merged_list <- rlist::list.append(merged_list, group_merged)
  }
  main_p_values <- do.call(plyr::rbind.fill, merged_list) %>%
    ungroup() %>%
    mutate(acc_p_adj = p.adjust(acc_p, method="BH"),
           bal_acc_p_adj = p.adjust(bal_acc_p, method="BH"))

  return(main_p_values)
}


