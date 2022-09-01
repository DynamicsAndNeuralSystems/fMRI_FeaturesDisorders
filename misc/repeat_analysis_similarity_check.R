library(tidyverse)

data_path <- "/headnode1/abry4213/data/UCLA_Schizophrenia/"

###############################################################
# catch22 

# Load in data for each of 5 catch22 runs after filtering and z-scoring
catch22_list <- list()
for (run_number in 1:5){
  run_number_catch22 <- readRDS(sprintf("%s/processed_data_run%s/Rdata/UCLA_Schizophrenia_catch22_filtered_zscored.Rds",
                                        data_path, run_number))
  catch22_list[[run_number]] <- run_number_catch22
}

# Check whether all five dataframes in terms are equal
all(apply(combn(length(catch22_list), 2), 2, function(x)
  all.equal(catch22_list[[x[1]]], catch22_list[[x[2]]])))
# TRUE woo

###############################################################
# Univariate linear SVM

univariate_ROI_wise_SVM_list <- list()
for (run_number in 1:5){
  run_number_ROI_SVM_balacc <- readRDS(sprintf("%s/processed_data_run%s/Rdata/ROI_wise_CV_linear_SVM_catch22_inv_prob_balacc.Rds",
                                        data_path, run_number)) %>%
    mutate(run_number = run_number)
  univariate_ROI_wise_SVM_list[[run_number]] <- run_number_ROI_SVM_balacc
}
univariate_ROI_wise_SVM_df <- do.call(plyr::rbind.fill, univariate_ROI_wise_SVM_list)

# Pull out left caudal anterior cingulate cortex as an example
univariate_ROI_wise_SVM_df %>%
  filter(grouping_var=="ctx-lh-caudalanteriorcingulate") %>%
  arrange(Noise_Proc)