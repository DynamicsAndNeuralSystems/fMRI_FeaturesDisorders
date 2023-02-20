from sklearnex import patch_sklearn
patch_sklearn()
import sys

# append the sys.path list
sys.path.insert(0, 'helper_functions/classification/')

from core_classification_functions import *
import argparse

# Command-line arguments to parse
parser = argparse.ArgumentParser(description='Process inputs for pairwise data preparation.')
parser.add_argument('--data_path', default="/headnode1/abry4213/data/UCLA_CNP/", dest='data_path')
parser.add_argument('--metadata_file', default="UCLA_CNP_sample_metadata.feather", dest='metadata_file')
parser.add_argument('--SPI_directionality_file', default="/headnode1/abry4213/github/fMRI_FeaturesDisorders/classification_analysis/SPI_Direction_Info.csv", dest='SPI_directionality_file')
parser.add_argument('--comparison_group', default="Schizophrenia", dest='comparison_group')
parser.add_argument('--univariate_feature_set', default='catch22', dest='univariate_feature_set')
parser.add_argument('--pairwise_feature_set', default='pyspi14', dest='pairwise_feature_set')
parser.add_argument('--pairwise_feature_file', default="/headnode1/abry4213/data/UCLA_CNP/processed_data/UCLA_CNP_AROMA_2P_GMR_pyspi14_filtered_zscored.feather", dest='pairwise_feature_file')
parser.add_argument('--noise_proc', dest='noise_proc')
parser.add_argument('--scaling_type', default="robustsigmoid", dest='scaling_type')
parser.add_argument('--num_folds', default=10, dest='num_folds')
parser.add_argument('--num_repeats', default=10, dest='num_repeats')
parser.add_argument('--num_jobs', default=8, dest='num_jobs')
parser.add_argument('--dataset_ID', default="UCLA_CNP", dest='dataset_ID')

# Parse arguments
args = parser.parse_args()
data_path = args.data_path
metadata_file = args.metadata_file
SPI_directionality_file = args.SPI_directionality_file
comparison_group = args.comparison_group
univariate_feature_set = args.univariate_feature_set
pairwise_feature_set = args.pairwise_feature_set
pairwise_feature_file = args.pairwise_feature_file
scaling_type = args.scaling_type
num_folds = args.num_folds
noise_proc = args.noise_proc
num_repeats = args.num_repeats
num_jobs = args.num_jobs
dataset_ID = args.dataset_ID

# dataset_ID = "UCLA_CNP"
# data_path = "/Users/abry4213/data/UCLA_CNP/"
# metadata_file = "UCLA_CNP_sample_metadata.feather"
# SPI_directionality_file = "/Users/abry4213/github/fMRI_FeaturesDisorders/classification_analysis/SPI_Direction_Info.csv"
# comparison_group = "Schizophrenia"
# univariate_feature_set = "catch22"
# pairwise_feature_set = "pyspi14"
# pairwise_feature_file ="/Users/abry4213/data/UCLA_CNP/processed_data/UCLA_CNP_AROMA_2P_GMR_pyspi14_filtered.feather"
# noise_proc = "AROMA+2P+GMR"
# scaling_type = "robustsigmoid"
# num_repeats = 10
# num_jobs = 1

# Run the pairwise main SVM
run_pairwise_SVM_by_SPI(pairwise_feature_file=pairwise_feature_file,
                 SPI_directionality_file = SPI_directionality_file,
                 univariate_feature_set=univariate_feature_set,
                 pairwise_feature_set=pairwise_feature_set,
                       noise_proc = noise_proc,
                       dataset_ID=dataset_ID,
                       metadata_file=metadata_file,
                       comparison_to_control_group=comparison_group,
                       pydata_path=data_path + "processed_data/",
                       data_path=data_path,
                       scaling_type = scaling_type,
                       num_repeats = int(num_repeats),
                       num_jobs = int(num_jobs))
