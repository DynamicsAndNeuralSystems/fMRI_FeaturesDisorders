# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--pairwise_data_file", default="/project/hctsa/annie/github/")
parser$add_argument("--SPI_directionality_file", default="/project/hctsa/annie/github/")
parser$add_argument("--rdata_path", default="/project/hctsa/annie/data/scz/UCLA/Rdata/")
parser$add_argument("--null_iter_number", default=1)
parser$add_argument("--feature_set", default="pyspi_19")
parser$add_argument("--svm_kernel", default="linear")
parser$add_argument("--grouping_var", default="region_pair")
parser$add_argument("--svm_feature_var", default="SPI")
parser$add_argument("--test_package", default="e1071")
parser$add_argument("--noise_proc", default="AROMA+2P+GMR")
parser$add_argument("--return_all_fold_metrics", action='store_true', default=TRUE)
parser$add_argument("--use_inv_prob_weighting", action='store_true', default=TRUE)
parser$add_argument("--use_SMOTE", action='store_true', default=FALSE)

# Parse input arguments
args <- parser$parse_args()
pairwise_data_file <- args$pairwise_data_file
SPI_directionality_file <- args$SPI_directionality_file
rdata_path <- args$rdata_path
null_iter_number <- args$null_iter_number
feature_set <- args$feature_set
svm_kernel <- args$svm_kernel
grouping_var <- args$grouping_var
svm_feature_var <- args$svm_feature_var
test_package <- args$test_package
noise_proc <- args$noise_proc
return_all_fold_metrics <- args$return_all_fold_metrics
use_inv_prob_weighting <- args$use_inv_prob_weighting
use_SMOTE <- args$use_SMOTE

# Load data
pairwise_data <- readRDS(pairwise_data_file)
SPI_directionality <- read.csv(SPI_directionality_file)

# Source linear SVM functions
source(paste0(github_dir, "helper_functions/Linear_SVM.R"))

# Define output directory
if (use_inv_prob_weighting) {
  output_dir <- paste0(rdata_path, sprintf("Pairwise_%s_inv_prob_null_model_fits/",
                                           feature_set))
} else {
  output_dir <- paste0(rdata_path, sprintf("Pairwise_%s_unweighted_null_model_fits/",
                                           feature_set))
}

icesTAF::mkdir(output_dir)

# Run null iteration
null_out <- run_pairwise_cv_svm_by_input_var(pairwise_data = pairwise_data,
                                             SPI_directionality = SPI_directionality,
                                             svm_kernel = svm_kernel,
                                             grouping_var = grouping_var,
                                             svm_feature_var = svm_feature_var,
                                             test_package = test_package,
                                             noise_proc = noise_proc,
                                             return_all_fold_metrics = return_all_fold_metrics,
                                             use_inv_prob_weighting = use_inv_prob_weighting,
                                             use_SMOTE = use_SMOTE,
                                             shuffle_labels = TRUE)

# Save null results to RDS
if (use_inv_prob_weighting) {
  saveRDS(null_out, file=sprintf("%s/Pairwise_%s_inv_prob_null_model_fit_iter_%s.Rds",
                                 output_dir, feature_set, null_iter_number))
} else {
  saveRDS(null_out, file=sprintf("%s/Pairwise_%s_unweighted_null_model_fit_iter_%s.Rds",
                                 output_dir, feature_set, null_iter_number))
}
