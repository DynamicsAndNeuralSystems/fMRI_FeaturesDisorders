library(tidyverse)

#param dataframe
param_df <- data.frame(Round = 1:4,
                       Sampling = c("in_sample", "in_sample", "CV", "CV"),
                       Weighting = c("none", "inv_prob", "inv_prob", "SMOTE"))


rdata_path <- "/home/abry4213/data/scz/UCLA/Rdata/"
seed <- 127

set.seed(seed)

# Compare all three noise processing methods
noise_procs = c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER")

# Use e1071 SVM with a linear kernel
test_package = "e1071"
svm_kernel = "linear"

# Run 1,000 permutations
num_permutations = 1000

# Iterate over param_df
for (i in 1:nrow(param_df)) {
  sampling = param_df$Sampling[i]
  weighting = param_df$Weighting[i]
  file_label = ifelse(weighting == "none", sampling, paste0(sampling, "_", weighting))
  
  use_inv_prob_weighting = ifelse(weighting == "inv_prob", TRUE, FALSE)
  use_SMOTE = ifelse(weighting == "SMOTE", TRUE, FALSE)
  cross_validate = ifelse(sampling == "CV", TRUE, FALSE)
  
  # Cross-validated, inverse probability reweighting
  if (!file.exists(paste0(rdata_path, 
                          sprintf("ROI_wise_model_permutation_null_%s.Rds",
                                  file_label)))) {
    
    model_permutation_null_CV_inv_prob <- run_null_model_n_permutations(rdata_path,
                                                                        noise_procs = noise_procs,
                                                                        grouping_var = "Brain_Region",
                                                                        svm_feature_var = "Feature",
                                                                        num_permutations = num_permutations,
                                                                        seed = seed,
                                                                        use_inv_prob_weighting = use_inv_prob_weighting,
                                                                        cross_validate = cross_validate,
                                                                        use_SMOTE = use_SMOTE)
    
    saveRDS(model_permutation_null_CV_inv_prob, file=paste0(rdata_path, 
                                                            sprintf("ROI_wise_model_permutation_null_%s.Rds",
                                                                    file_label)))
  }
}