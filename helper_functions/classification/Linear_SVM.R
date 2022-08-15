#------------------------------------
# Functions to run linear SVM
# Options to use cross-validation with
# Inverse probability weighting
#------------------------------------

#--------------------------------------
# Author: Annie Bryant, 16 May 2022
#--------------------------------------

library(tidyverse)
library(e1071)
library(kernlab)
library(cowplot)
library(caret)
theme_set(theme_cowplot())

# DIY rlist::list.append
list.append <- function (.data, ...) 
{
  if (is.list(.data)) {
    c(.data, list(...))
  }
  else {
    c(.data, ..., recursive = FALSE)
  }
}


#-------------------------------------------------------------------------------
# k-fold CV SVM with option to use inverse probability weighting
#-------------------------------------------------------------------------------
k_fold_CV_linear_SVM <- function(input_data,
                                 flds = NULL,
                                 k = k,
                                 c = 1,
                                 svm_kernel = "linear",
                                 sample_wts = list("Control" = 1,
                                                   "Schz" = 1),
                                 out_of_sample_only = TRUE,
                                 shuffle_labels = FALSE) {
  
  # Shuffle labels if specified
  if (shuffle_labels) {
    input_data <- transform(input_data, Diagnosis = sample(Diagnosis, replace = FALSE))
  }
  
  # Specify that Diagnosis is a factor so that createFolds creates stratified folds
  input_data$Diagnosis <- factor(input_data$Diagnosis)
  
  # Use pre-specified folds if passed in,
  # Otherwise create train/test data folds
  if (is.null(flds)) {
    flds <- caret::createFolds(input_data$Diagnosis, k = k, list = TRUE, returnTrain = FALSE)
  }
  
  # Create dataframe to store subject IDs and whether or not they were properly classified
  subject_classification_list <- list()
  
  
  # Iterate over folds 1 through k
  for (i in 1:k) {
    
    # Define test and train data
    test_i <- flds[[i]]
    train_i <- setdiff(1:nrow(input_data), test_i)
    
    # Training data
    train_data <- input_data[train_i, ] %>%
      filter(!is.na(Diagnosis))
    train_subjects <- train_data$Sample_ID
    
    train_data <- train_data %>% dplyr::select(-Sample_ID) 
    
    # Testing data
    test_data <- input_data[test_i, ] %>%
      filter(!is.na(Diagnosis))
    test_subjects <- test_data$Sample_ID
    
    test_data <- test_data %>% dplyr::select(-Sample_ID) 
    
    # Run linear SVM on fold
    
    svmModel <- e1071::svm(factor(Diagnosis) ~ .,
                           kernel = svm_kernel,
                           cost = c,
                           data = train_data,
                           class.weights = sample_wts)
    
    if (!out_of_sample_only) {
      # Generate in-sample predictions based on SVM model
      in_sample_pred <- predict(svmModel, train_data)
      train_data$Diagnosis <- factor(train_data$Diagnosis, levels = levels(in_sample_pred))
      
      # Create dataframe containing subject ID and whether out-of-sample prediction was correct
      in_fold_predictions_by_subject <- data.frame(Sample_ID = train_subjects,
                                                   Sample_Type = "In-sample",
                                                   fold_number = i,
                                                   Actual_Diagnosis = train_data$Diagnosis,
                                                   Predicted_Diagnosis = in_sample_pred) %>%
        mutate(Prediction_Correct = Actual_Diagnosis == Predicted_Diagnosis)
      subject_classification_list <- list.append(subject_classification_list,
                                            in_fold_predictions_by_subject)
    }

    # Generate out-of-sample predictions based on SVM model
    out_sample_pred <- predict(svmModel, test_data)
    test_data$Diagnosis <- factor(test_data$Diagnosis, levels = levels(out_sample_pred))
    
    # Create dataframe containing subject ID and whether out-of-sample prediction was correct
    out_fold_predictions_by_subject <- data.frame(Sample_ID = test_subjects,
                                                  Sample_Type = "Out-of-sample",
                                                  fold_number = i,
                                                  Actual_Diagnosis = test_data$Diagnosis,
                                                  Predicted_Diagnosis = out_sample_pred) %>%
      mutate(Prediction_Correct = Actual_Diagnosis == Predicted_Diagnosis)
    
    subject_classification_list <- list.append(subject_classification_list,
                                          out_fold_predictions_by_subject)
    
  } 
  
  # Compile classification res
  classification_res <- do.call(plyr::rbind.fill, subject_classification_list)
  
  # Subset to just out-of-sample data if requested
  if (out_of_sample_only) {
    classification_res <- classification_res %>%
      filter(Sample_Type == "Out-of-sample")
  }
  
  # Return results 
  return(classification_res)
}

#-------------------------------------------------------------------------------
# Run 10-fold cross-validated multi-feature linear SVM by given grouping var
#-------------------------------------------------------------------------------

run_univariate_cv_svm_by_input_var <- function(data_path,
                                               dataset_ID,
                                               svm_kernel = "linear",
                                               univariate_feature_set = "catch22",
                                               pairwise_feature_set = "pyspi_19",
                                               grouping_var = "Brain_Region",
                                               svm_feature_var = "Feature",
                                               noise_procs = c("AROMA+2P", 
                                                               "AROMA+2P+GMR", 
                                                               "AROMA+2P+DiCER"),
                                               flds = NULL,
                                               num_k_folds = 10,
                                               out_of_sample_only = TRUE,
                                               use_inv_prob_weighting = FALSE,
                                               shuffle_labels = FALSE) {
  
  rdata_path <- paste0(data_path, "Rdata/")
  
  # Get diagnosis proportions
  sample_groups <- readRDS(paste0(data_path, sprintf("%s_samples_with_univariate_%s_and_pairwise_%s.Rds",
                                                     dataset_ID,
                                                     univariate_feature_set,
                                                     pairwise_feature_set))) %>%
    dplyr::select(Sample_ID, Diagnosis)
  
  # Define sample weights
  # Default is 1 and 1 if use_inv_prob_weighting is not included
  if (use_inv_prob_weighting) {
    
    # Convert to sample weights based on inverse of probability
    sample_wts <- as.list(1/prop.table(table(sample_groups$Diagnosis)))
    
  } else {
    sample_wts <- as.list(rep(1, length(unique(sample_groups$Diagnosis))))
    names(sample_wts) = unique(sample_groups$Diagnosis)
  }
  
  class_res_list <- list()
  for (noise_proc in noise_procs) {
    # Clean up names
    noise_label <- gsub("\\+", "_", noise_proc)
    
    # Load z-scored feature data for current noise processing method and
    # time-series feature set
    feature_matrix <- readRDS(paste0(rdata_path, sprintf("%s_%s_filtered_zscored.Rds", 
                                                         dataset_ID,
                                                         univariate_feature_set))) %>%
      dplyr::filter(Noise_Proc == noise_proc)
    
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
    
    
    # Reshape data from long to wide for SVM
    for (group_var in grouping_var_vector) {
      if (grouping_var == "Combo") {
        data_for_SVM <- feature_matrix %>%
          left_join(., sample_groups) %>%
          dplyr::select(Sample_ID, Diagnosis, Combo, values) %>%
          distinct(.keep_all = T) %>%
          tidyr::pivot_wider(id_cols = c(Sample_ID, Diagnosis),
                             names_from = Combo,
                             values_from 
                             = values) %>%
          # Drop columns that are all NA/NAN
          dplyr::select(where(function(x) any(!is.na(x)))) %>%
          # Drop rows with NA for one or more column
          drop_na()
        
      } else {
        # Otherwise iterate over each separate group
        data_for_SVM <- subset(feature_matrix, get(grouping_var_name) == group_var) %>%
          dplyr::ungroup() %>%
          left_join(., sample_groups) %>%
          dplyr::select(Sample_ID, Diagnosis, svm_feature_var_name, values) %>%
          distinct(.keep_all = T) %>%
          tidyr::pivot_wider(id_cols = c(Sample_ID, Diagnosis),
                             names_from = svm_feature_var_name,
                             values_from 
                             = values) %>%
          # Drop columns that are all NA/NAN
          dplyr::select(where(function(x) any(!is.na(x)))) %>%
          # Drop rows with NA for one or more column
          drop_na()
      }
      
      # Pass data_for_SVM to in_sample_linear_SVM
      SVM_results <- k_fold_CV_linear_SVM(input_data = data_for_SVM,
                                          flds = flds,
                                          k = num_k_folds,
                                          svm_kernel = svm_kernel,
                                          sample_wts = sample_wts,
                                          shuffle_labels = shuffle_labels,
                                          out_of_sample_only = out_of_sample_only) %>%
        dplyr::mutate(grouping_var = group_var,
                      feature_set = univariate_feature_set,
                      Noise_Proc = noise_proc)
      
      # Append results to list
      class_res_list <- list.append(class_res_list, SVM_results)
    }
  }
  
  # Combine results from all regions into a dataframe
  class_res_df <- do.call(plyr::rbind.fill, class_res_list)
  
  # Return dataframe
  return(class_res_df)
}

#-------------------------------------------------------------------------------
# Run pairwise PYSPI multi-feature linear SVM by input feature
#-------------------------------------------------------------------------------
run_pairwise_cv_svm_by_input_var <- function(pairwise_data,
                                             data_path,
                                             SPI_directionality,
                                             svm_kernel = "linear",
                                             grouping_var = "SPI",
                                             svm_feature_var = "region_pair",
                                             noise_proc = "AROMA+2P+GMR",
                                             num_k_folds = 10,
                                             flds = NULL,
                                             out_of_sample_only = TRUE,
                                             use_inv_prob_weighting = FALSE,
                                             shuffle_labels = FALSE) {
  
  
  rdata_path <- paste0(data_path, "Rdata/")
  
  # Get diagnosis proportions
  sample_groups <- readRDS(paste0(data_path, sprintf("%s_samples_with_univariate_%s_and_pairwise_%s.Rds",
                                                     dataset_ID,
                                                     univariate_feature_set,
                                                     pairwise_feature_set))) %>%
    dplyr::select(Sample_ID, Diagnosis)
  
  # Define sample weights
  # Default is 1 and 1 if use_inv_prob_weighting is not included
  if (use_inv_prob_weighting) {
    
    # Convert to sample weights based on inverse of probability
    sample_wts <- as.list(1/prop.table(table(sample_groups$Diagnosis)))
    
  } else {
    sample_wts <- as.list(rep(1, length(unique(sample_groups$Diagnosis))))
    names(sample_wts) = unique(sample_groups$Diagnosis)
  }
  
  # Initialize results list for SVM
  class_res_list <- list()
  
  # Combine region pair names
  if (svm_feature_var == "region_pair") {
    svm_feature_var_name = svm_feature_var
    grouping_var_name = "SPI"
    grouping_var_vector <- unique(pairwise_data$SPI)
    
    # Filter by directionality
    pairwise_data <- pairwise_data %>%
      dplyr::rename("group_SPI" = "SPI") %>%
      group_by(group_SPI) %>%
      mutate(Direction = SPI_directionality %>% 
               filter(SPI == unique(group_SPI)) %>% 
               distinct(Direction) %>%
               pull(Direction)) %>%
      dplyr::rename("SPI" = "group_SPI") %>%
      mutate(region_pair = case_when(Direction == "Undirected" ~ ifelse(brain_region_1 < brain_region_2,
                                                                        paste0(brain_region_1, "_", brain_region_2),
                                                                        paste0(brain_region_2, "_", brain_region_1)),
                                     Direction == "Directed" ~ paste0(brain_region_1, "_", brain_region_2))) %>%
      dplyr::select(-brain_region_1, -brain_region_2)  %>%
      distinct(Sample_ID, SPI, region_pair, .keep_all = T)
    
  } else if (svm_feature_var == "SPI") {
    
    # Don't want to filter by directionality
    pairwise_data <- pairwise_data %>%
      rowwise() %>%
      tidyr::unite("region_pair", c(brain_region_1, brain_region_2), sep="_") %>%
      distinct(Sample_ID, SPI, region_pair, .keep_all = T)
    
    svm_feature_var_name = svm_feature_var
    grouping_var_name = "region_pair"
    grouping_var_vector <- unique(pairwise_data$region_pair)
    
  } else {
    svm_feature_var_name = "Combo"
    grouping_var_name = "Group_Var"
    
    # Filter by directionality
    pyspi_data <- pyspi_data %>%
      # Special cases
      filter(SPI != "sgc_nonparametric_mean_fs-1_fmin-0_fmax-0-5",
             !(Sample_ID == "sub-10171" & SPI == "di_gaussian")) %>%
      rowwise() %>%
      tidyr::unite("region_pair", c(brain_region_1, brain_region_2), sep="_") %>%
      distinct(Sample_ID, SPI, region_pair, .keep_all = T) %>%
      group_by(SPI, region_pair) %>%
      filter(!all(is.na(value))) %>%
      dplyr::select(where(function(x) any(!is.na(x)))) %>%
      unite("Combo", c("region_pair", "SPI"), sep="_", remove=F)
    
    grouping_var_vector <- c("All")
    
  }
  
  # Reshape data from long to wide for SVM
  for (group_var in unique(grouping_var_vector)) {
    if (grouping_var == "Combo") {
      data_for_SVM <- pairwise_data %>%
        # Impute missing data with the mean
        group_by(Diagnosis, Combo) %>%
        mutate(value = ifelse(is.na(value), mean(value, na.rm=T), value)) %>%
        dplyr::select(Sample_ID, Diagnosis, Combo, value) %>%
        tidyr::pivot_wider(id_cols = c(Sample_ID, Diagnosis),
                           names_from = Combo,
                           values_from 
                           = value) %>%
        # Drop columns that are all NA/NAN
        dplyr::select(where(function(x) any(!is.na(x)))) %>%
        # Drop rows with NA for one or more column
        drop_na()
      
    } else {
      # Otherwise iterate over each separate group
      data_for_SVM <- subset(pairwise_data, get(grouping_var_name) == group_var) %>%
        dplyr::ungroup() %>%
        dplyr::select(Sample_ID, Diagnosis, svm_feature_var_name, value) %>%
        tidyr::pivot_wider(id_cols = c(Sample_ID, Diagnosis),
                           names_from = svm_feature_var_name,
                           values_from 
                           = value) %>%
        # Drop columns that are all NA/NAN
        dplyr::select(where(function(x) any(!is.na(x)))) %>%
        # Drop rows with NA for one or more column
        drop_na()
    }
    
    if (nrow(data_for_SVM) > 0) {
      # Run k-fold linear SVM
      SVM_results <- k_fold_CV_linear_SVM(input_data = data_for_SVM,
                                          flds = flds,
                                          k = num_k_folds,
                                          svm_kernel = svm_kernel,
                                          sample_wts = sample_wts,
                                          shuffle_labels = shuffle_labels,
                                          out_of_sample_only = out_of_sample_only) %>%
        dplyr::mutate(grouping_var = group_var,
                      pairwise_feature_set = pairwise_feature_set,
                      use_inv_prob_weighting = use_inv_prob_weighting,
                      Noise_Proc = noise_proc)
      
      # Append results to list
      class_res_list <- list.append(class_res_list, SVM_results)
    } else {
      cat("\nNo observations available for", group_var, "after filtering.\n")
    }
  }
  
  # Combine results from all regions into a dataframe
  class_res_df <- do.call(plyr::rbind.fill, class_res_list)
  
  # Return dataframe
  return(class_res_df)
}

#-------------------------------------------------------------------------------
# Run combined univariate theft plus pairwise PYSPI 
# multi-feature linear SVM by input feature
#-------------------------------------------------------------------------------
run_combined_uni_pairwise_cv_svm_by_input_var <- function(univariate_data,
                                                          univariate_feature_set,
                                                          pairwise_data,
                                                          pairwise_feature_set,
                                                          SPI_directionality,
                                                          num_k_folds = 10,
                                                          flds = NULL,
                                                          svm_kernel = "linear",
                                                          test_package = "e1071",
                                                          noise_proc = "AROMA+2P+GMR",
                                                          return_all_fold_metrics = FALSE,
                                                          use_inv_prob_weighting = FALSE,
                                                          shuffle_labels = FALSE) {
  
  # Initialize results list for SVM
  class_res_list <- list()
  
  # Combine region pair names
  pairwise_data <- pairwise_data %>%
    left_join(., SPI_directionality) %>%
    rowwise() %>%
    mutate(region_pair = case_when(Direction == "Undirected" ~ ifelse(brain_region_1 < brain_region_2,
                                                                      paste0(brain_region_1, "_", brain_region_2),
                                                                      paste0(brain_region_2, "_", brain_region_1)),
                                   Direction == "Directed" ~ paste0(brain_region_1, "_", brain_region_2))) %>%
    dplyr::select(-brain_region_1, -brain_region_2)  %>%
    distinct(Sample_ID, SPI, region_pair, .keep_all = T)
  
  # Merge ROI plus theft feature for univariate
  univariate_combo <- univariate_data %>%
    tidyr::unite("Unique_ID", c("names", "Brain_Region"), sep="_") %>%
    dplyr::select(Sample_ID, group, Unique_ID, values)
  
  # Initialize list for each SPI
  class_res_list <- list()
  
  # Split pairwise data by SPI
  for (this_SPI in unique(pairwise_data$SPI)) {
    # Merge region-pair plus SPI data for pairwise
    pairwise_combo <- pairwise_data %>%
      filter(SPI == this_SPI) %>%
      tidyr::unite("Unique_ID", c("SPI", "region_pair"), sep="_") %>%
      dplyr::select(Sample_ID, group, Unique_ID, value) %>%
      dplyr::rename("values"="value")
    
    # Combine univariate + pairwise data for SVM
    combined_data_for_SVM <- plyr::rbind.fill(univariate_combo, pairwise_combo) %>%
      tidyr::pivot_wider(id_cols = c(Sample_ID, group),
                         names_from = Unique_ID, 
                         values_from = values) %>%
      # Drop columns that are all NA/NAN
      dplyr::select(where(function(x) any(!is.na(x)))) %>%
      # Drop rows with NA for one or more column
      drop_na()
    
    # Define sample weights
    # Default is 1 and 1 if use_inv_prob_weighting is not included
    if (use_inv_prob_weighting) {
      # Get control/schz proportions
      sample_wts <- as.list(1/prop.table(table(combined_data_for_SVM$group)))
    } else {
      sample_wts <- list("Control" = 1, "Schz" = 1)
    }
    
    if (nrow(combined_data_for_SVM) > 0) {
      # Run k-fold linear SVM
      SVM_results <- k_fold_CV_linear_SVM(input_data = combined_data_for_SVM,
                                          k = num_k_folds,
                                          flds = flds,
                                          svm_kernel = svm_kernel,
                                          sample_wts = sample_wts,
                                          shuffle_labels = shuffle_labels,
                                          return_all_fold_metrics = return_all_fold_metrics) %>%
        dplyr::mutate(SPI = this_SPI,
                      univariate_feature_set = univariate_feature_set,
                      pairwise_feature_set = pairwise_feature_set,
                      Noise_Proc = noise_proc,
                      use_inv_prob_weighting = use_inv_prob_weighting)
      
      # Append results to list
      class_res_list <- list.append(class_res_list, SVM_results)
    } else {
      cat("\nNo observations available for", this_SPI, "after filtering.\n")
      
    }
  }
  # Combine results from all regions into a dataframe
  class_res_df <- do.call(plyr::rbind.fill, class_res_list)
  
  # Return dataframe
  return(class_res_df)
}




