
rdata_path <- "/home/abry4213/data/scz/UCLA/Rdata/"
seed <- 127

set.seed(seed)

# Compare all three noise processing methods
noise_procs = c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER")

# Use e1071 SVM with a linear kernel
test_package = "e1071"
kernel = "linear"

if (!file.exists(paste0(rdata_path, "ROI_wise_model_permutation_null_in_sample.Rds"))) {
  model_permutation_null_in_sample <- run_null_model_n_permutations(rdata_path,
                                                                    noise_procs = noise_procs,
                                                                    grouping_var = "Brain_Region",
                                                                    svm_feature_var = "Feature",
                                                                    num_permutations = 1000,
                                                                    seed = 127,
                                                                    use_inv_prob_weighting = FALSE,
                                                                    cross_validate = TRUE)
  
  saveRDS(model_permutation_null_in_sample, file=paste0(rdata_path, "ROI_wise_model_permutation_null_in_sample.Rds"))
} else {
  model_permutation_null_in_sample <- readRDS(paste0(rdata_path, "ROI_wise_model_permutation_null_in_sample.Rds"))
}