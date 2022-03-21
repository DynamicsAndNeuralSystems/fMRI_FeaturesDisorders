#------------------------------------
# This script sets out to produce a
# function for reading in matlab time
# series files into R
#------------------------------------

#--------------------------------------
# Author: Annie Bryant, 15 March 2022
#--------------------------------------

library(tidyverse)
library(theft)
library(caret)


#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Univariate classification analysis for one individual region

region_by_region_analysis <- function(ROI, feature_mat_region, test_method="svmLinear", display_figures=F, return_restable=T) {
  
  # Find the top ten features for given brain region
  classifier_outputs <- compute_top_features(feature_mat_region, 
                                             id_var = "Subject_ID", 
                                             group_var = "group",
                                             num_features = 22, 
                                             normalise_violin_plots = FALSE,
                                             method = "RobustSigmoid",
                                             cor_method = "pearson",
                                             test_method = test_method,
                                             use_k_fold = TRUE,
                                             num_folds = 10,
                                             use_empirical_null =  TRUE,
                                             null_testing_method = "model free shuffles",
                                             p_value_method = "empirical",
                                             num_permutations = 10,
                                             pool_empirical_null = FALSE)
  
  if (display_figures) {
    # Plot the dimensionality-reduced time-series data along PC1 and PC2
    print(plot_low_dimension(feature_mat_region, 
                             is_normalised = FALSE, 
                             id_var = "Subject_ID", 
                             group_var = "group", 
                             method = "z-score", 
                             low_dim_method = "PCA", 
                             plot = TRUE,
                             show_covariance = TRUE))
    
    
    # Plot the violin plots for the top 10 featIures
    print(classifier_outputs$ViolinPlots)
    
  } 
  
  if (return_restable) {
    return(classifier_outputs$ResultsTable)
  }
}

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Run univariate classification across all ROIs
run_region_by_region_analysis <- function(feature_matrix, test_method, rdata_path, noise_proc) {
  
  test_label <- gsub("-", "_", test_method)
  noise_label <- gsub("\\+", "_", noise_proc)
  region_wise_univ_class_res_list <- list()
  
  for (region in unique(feature_matrix$Brain_Region)) {
    feature_mat_region <- subset(feature_matrix, Brain_Region==region)
    
    tryCatch({
      region_res <- region_by_region_analysis(ROI=region,
                                              test_method=test_method,
                                              feature_mat_region=feature_mat_region,
                                              display_figures = F,
                                              return_restable = T)
      region_res$Brain_Region <- region
      region_wise_univ_class_res_list[[region]] <- region_res
    }, error = function(e) {
      cat("\nError with ROI", region, ":")
      message(e)
    })

  }
  region_wise_univ_class_res <- do.call(plyr::rbind.fill, region_wise_univ_class_res_list)
  saveRDS(region_wise_univ_class_res, paste0(rdata_path, "UCLA_", 
                                             noise_label, "_catch22_ROIwise_",  
                                             test_label, ".Rds"))
}