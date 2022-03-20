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

region_by_region_analysis <- function(ROI, feature_matrix, display_figures=F, return_restable=T) {
  feature_mat_region <- subset(feature_matrix, Brain_Region==ROI)
  
  # Find the top ten features for given brain region
  classifier_outputs <- compute_top_features(feature_mat_region, 
                                             id_var = "Subject_ID", 
                                             group_var = "group",
                                             num_features = 22, 
                                             normalise_violin_plots = FALSE,
                                             method = "RobustSigmoid",
                                             cor_method = "pearson",
                                             test_method = "svmLinear",
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
    
    
    # Plot the violin plots for the top 10 features
    print(classifier_outputs$ViolinPlots)
    
  } 
  
  if (return_restable) {
    return(classifier_outputs$ResultsTable)
  }
}