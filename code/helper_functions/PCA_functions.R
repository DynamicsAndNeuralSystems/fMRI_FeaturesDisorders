#------------------------------------
# This script compiles movement data across the UCLA subjects
#------------------------------------

#--------------------------------------
# Author: Annie Bryant, 5 June 2022
#--------------------------------------

library(FactoMineR)
library(factoextra)

#-------------------------------------------------------------------------------
# Function to read in univariate TS feature data and return subjects with NA 
# values. For each noise-processing method, the number of TS features with all 
# NA for all brain regions are given.
#-------------------------------------------------------------------------------

run_PCA_by_group_var <- function(feature_matrix,
                                 grouping_variable = "Brain_Region",
                                 feature_var = "names") {
  
  if (grouping_variable == "Combo") {
    grouping_var_vector <- c("All")
    feature_matrix <- feature_matrix %>%
      tidyr::unite("feature_var", c(names, Brain_Region), sep="_") %>%
      mutate(grouping_variable = "Combo")
  } else{
    grouping_var_vector <- feature_matrix %>%
      dplyr::pull(get(grouping_variable)) %>%
      unique()
    
    feature_matrix <- feature_matrix %>%
      dplyr::rename("grouping_variable" = grouping_variable,
                    "feature_var" = feature_var)
  }
  
  PCA_res_list <- list()
  for (this_group in unique(feature_matrix$grouping_variable)) {
    data_for_PCA <- subset(feature_matrix, 
                           grouping_variable == this_group) %>%
      pivot_wider(id_cols = c(Subject_ID, group),
                  names_from = feature_var,
                  values_from = values) %>%
      drop_na()
    
    data_for_PCA_mat <- data_for_PCA %>%
      dplyr::select(-Subject_ID, -group) %>%
      as.matrix()
    
    PCA_res <- prcomp(data_for_PCA_mat, center = TRUE, scale. = TRUE)
    PCA_res_list[[this_group]] <- PCA_res
  }
  return(PCA_res_list)
}

run_PCA_for_uni_pairwise_combo <- function(univariate_data,
                                           pairwise_data) {
  
  # Merge ROI plus theft feature for univariate
  univariate_combo <- univariate_data %>%
    tidyr::unite("Unique_ID", c("names", "Brain_Region"), sep="_") %>%
    dplyr::select(Subject_ID, group, Unique_ID, values)
  
  # Combine region pair names
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
    distinct(Subject_ID, SPI, region_pair, .keep_all = T)
  
  PCA_res_list <- list()
  for (this_SPI in unique(pairwise_data$SPI)) {
    # Merge region-pair plus SPI data for pairwise
    pairwise_combo <- pairwise_data %>%
      filter(SPI == this_SPI) %>%
      tidyr::unite("Unique_ID", c("SPI", "region_pair"), sep="_") %>%
      dplyr::select(Subject_ID, group, Unique_ID, value) %>%
      dplyr::rename("values"="value")
    
    # Combine univariate + pairwise data for PCA
    combined_data_for_PCA <- plyr::rbind.fill(univariate_combo, pairwise_combo) %>%
      tidyr::pivot_wider(id_cols = c(Subject_ID, group),
                         names_from = Unique_ID, 
                         values_from = values) %>%
      # Drop columns that are all NA/NAN
      dplyr::select(where(function(x) any(!is.na(x)))) %>%
      # Drop rows with NA for one or more column
      drop_na()
    
    if (nrow(combined_data_for_PCA > 0)) {
      data_for_PCA_mat <- combined_data_for_PCA %>%
        dplyr::select(-group, -Subject_ID) %>%
        as.matrix()
      
      PCA_res <- prcomp(data_for_PCA_mat, center = F, scale. = F)
      PCA_res_list[[this_SPI]] <- PCA_res
    } else {
      cat("\nNo available data for", this_SPI, "\n")
    }
    
  }
  return(PCA_res_list)
}

#-------------------------------------------------------------------------------
# Function to run linear SVM with increasing number of PCs
#-------------------------------------------------------------------------------

run_SVM_from_PCA <- function(list_of_PCA_res,
                             subject_dx_list,
                             c_values = c(1),
                             interval = 1,
                             use_inv_prob_weighting = FALSE,
                             return_all_fold_metrics = FALSE) {
  
  # Initialize results list
  PCA_SVM_res_list <- list()
  
  # Get names of groups to iterate over
  grouping_var_vector <- names(list_of_PCA_res)
  
  # Iterate over each grouping variable
  for (group in grouping_var_vector) {
    PCA_res <- list_of_PCA_res[[group]]
    total_n_PCs <- length(PCA_res$sdev)
    
    # Start from 1
    starting_i <- 1
    
    # Increasingly iterate over each PCs
    for (i in seq(starting_i, total_n_PCs, by = interval)) {
      svm_for_pc <- as.data.frame(cbind(group = subject_dx_list, PCA_res$x[, 1:i])) %>%
        mutate_at(vars(contains("V")), as.numeric) %>%
        mutate_at(vars(starts_with("PC")), as.numeric) 
      
      if (use_inv_prob_weighting) {
        sample_wts <- 1/prop.table(table(subject_dx_list))
      } else {
        sample_wts <- list("Control" = 1, "Schz" = 1)
      }
      
      # Compile results into a dataframe
      df_res <- k_fold_CV_linear_SVM(input_data = svm_for_pc,
                                     k = 10,
                                     c_values = c_values,
                                     svm_kernel = "linear",
                                     sample_wts = sample_wts,
                                     shuffle_labels = FALSE,
                                     return_all_fold_metrics = return_all_fold_metrics)
      
      df_res$Num_PCs <- i
      df_res$grouping_var <- group
      
      # Append results to list
      PCA_SVM_res_list <- rlist::list.append(PCA_SVM_res_list, df_res)
    }
  }
  PCA_SVM_res <- do.call(plyr::rbind.fill, PCA_SVM_res_list)
  return(PCA_SVM_res)
}