# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--rdata_path", default="/project/hctsa/annie/data/scz/UCLA/Rdata/")
parser$add_argument("--output_data_dir", default="/project/hctsa/annie/data/scz/UCLA/Rdata/Univariate_catch22_inv_prob_null_model_fits/")
parser$add_argument("--github_dir", default="/project/hctsa/annie/github/fMRI_FeaturesDisorders/")
parser$add_argument("--grouping_var", default="Feature")
parser$add_argument("--num_k_folds", default=10)
parser$add_argument("--null_iter_number", default=1)
parser$add_argument("--num_perms_for_iter", default=1)
parser$add_argument("--feature_set", default="catch22")
parser$add_argument("--svm_kernel", default="linear")
parser$add_argument("--svm_feature_var", default="Brain_Region")
parser$add_argument("--test_package", default="e1071")
parser$add_argument("--return_all_fold_metrics", action='store_true', default=FALSE)
parser$add_argument("--weighting_name", default="unweighted")
parser$add_argument("--use_inv_prob_weighting", action='store_true', default=FALSE)
parser$add_argument("--use_SMOTE", action='store_true', default=FALSE)

# Parse input arguments
args <- parser$parse_args()
rdata_path <- args$rdata_path
output_data_dir <- args$output_data_dir
github_dir <- args$github_dir
grouping_var <- args$grouping_var
num_k_folds <- as.numeric(args$num_k_folds)
null_iter_number <- args$null_iter_number
num_perms_for_iter <- args$num_perms_for_iter
feature_set <- args$feature_set
svm_kernel <- args$svm_kernel
grouping_var <- args$grouping_var
svm_feature_var <- args$svm_feature_var
test_package <- args$test_package
return_all_fold_metrics <- args$return_all_fold_metrics
weighting_name <- args$weighting_name
use_inv_prob_weighting <- args$use_inv_prob_weighting
use_SMOTE <- args$use_SMOTE

# Source linear SVM functions
source(paste0(github_dir, "helper_functions/Linear_SVM.R"))

icesTAF::mkdir(output_data_dir)

cat("\nNumber of k-folds:", num_k_folds, "\n")
cat("\nNum permutations per iteration:", num_perms_for_iter, "\n")
cat("\nData type:", typeof(num_perms_for_iter), "\n")
# Run null iteration
null_out <- 1:num_perms_for_iter  %>%
  purrr::map_df( ~ run_univariate_cv_svm_by_input_var(rdata_path = rdata_path,
                                                      svm_kernel = svm_kernel,
                                                      feature_set = feature_set,
                                                      test_package = test_package,
                                                      grouping_var = grouping_var,
                                                      svm_feature_var = svm_feature_var,
                                                      noise_procs = c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER"),
                                                      num_k_folds = num_k_folds,
                                                      use_inv_prob_weighting = use_inv_prob_weighting,
                                                      use_SMOTE = use_SMOTE,
                                                      shuffle_labels = T) %>%
                   # Keep track of which null iteration this is
                   mutate(Null_Iter_Number = .x + (.x * (as.numeric(null_iter_number) - 1))))

# Save null results to RDS
saveRDS(null_out, file=sprintf("%s/%s_wise_%s_%s_null_model_fit_iter_%s.Rds",
                               output_data_dir, grouping_var, feature_set, 
                               weighting_name, null_iter_number))