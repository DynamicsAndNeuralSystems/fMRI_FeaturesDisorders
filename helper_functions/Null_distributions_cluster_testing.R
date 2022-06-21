# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

pairwise_data_file="/project/hctsa/annie/data/scz/UCLA/Rdata/pyspi_SPI_pairwise_CV_linear_SVM_pyspi_19_inv_prob_small.Rds"
SPI_directionality_file="/project/hctsa/annie/github/pairwise_analysis/SPI_Direction_Info.csv"
rdata_path="/project/hctsa/annie/data/scz/UCLA/Rdata/"
github_dir="/project/hctsa/annie/github/fMRI_FeaturesDisorders/"
null_iter_number=1
feature_set="pyspi_19"
svm_kernel="linear"
grouping_var="region_pair"
svm_feature_var="SPI"
test_package="e1071"
noise_proc="AROMA+2P+GMR"
return_all_fold_metrics=TRUE
use_inv_prob_weighting=TRUE
use_SMOTE=FALSE

# Load data
cat("\nNow reading in pyspi RDS file\n")
pairwise_data <- readRDS(pairwise_data_file)
cat("\nNow reading in SPI directionality file\n")
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

cat("\nHead of pairwise data:\n")
head(pairwise_data)
cat("\nHead of SPI directionality data:\n")
head(SPI_directionality)
# 
# # Run null iteration
# null_out <- run_pairwise_cv_svm_by_input_var(pairwise_data = pairwise_data,
#                                              SPI_directionality = SPI_directionality,
#                                              svm_kernel = svm_kernel,
#                                              grouping_var = grouping_var,
#                                              svm_feature_var = svm_feature_var,
#                                              test_package = test_package,
#                                              noise_proc = noise_proc,
#                                              return_all_fold_metrics = return_all_fold_metrics,
#                                              use_inv_prob_weighting = use_inv_prob_weighting,
#                                              use_SMOTE = use_SMOTE,
#                                              shuffle_labels = TRUE)
# 
# # Save null results to RDS
# if (use_inv_prob_weighting) {
#   saveRDS(null_out, file=sprintf("%s/Pairwise_%s_inv_prob_null_model_fit_iter_%s.Rds",
#                                  output_dir, feature_set, null_iter_number))
# } else {
#   saveRDS(null_out, file=sprintf("%s/Pairwise_%s_unweighted_null_model_fit_iter_%s.Rds",
#                                  output_dir, feature_set, null_iter_number))
# }
