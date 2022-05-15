source("SVM_functions.R")
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

# Cross-validated, SMOTE
if (!file.exists(paste0(rdata_path, "ROI_wise_model_permutation_null_CV_SMOTE.Rds"))) {
  model_permutation_null_CV_SMOTE <- run_null_model_n_permutations(rdata_path,
                                                                   noise_procs = noise_procs,
                                                                   grouping_var = "Brain_Region",
                                                                   svm_feature_var = "Feature",
                                                                   num_permutations = num_permutations,
                                                                   seed = seed,
                                                                   use_inv_prob_weighting = FALSE,
                                                                   cross_validate = TRUE,
                                                                   use_SMOTE = TRUE)
  
  saveRDS(model_permutation_null_CV_SMOTE, file=paste0(rdata_path, "ROI_wise_model_permutation_null_CV_SMOTE.Rds"))
} else {
  model_permutation_null_CV_SMOTE <- readRDS(paste0(rdata_path, "ROI_wise_model_permutation_null_CV_SMOTE.Rds"))
}