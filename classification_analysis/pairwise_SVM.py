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
parser.add_argument('--num_null_iters', default=1000, dest='num_null_iters')
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
noise_proc = args.noise_proc
num_null_iters = args.num_null_iters
dataset_ID = args.dataset_ID

# dataset_ID = "UCLA_CNP"
# data_path = "/headnode1/abry4213/data/UCLA_CNP/"
# metadata_file = "UCLA_CNP_sample_metadata.feather"
# SPI_directionality_file = "/headnode1/abry4213/github/fMRI_FeaturesDisorders/classification_analysis/SPI_Direction_Info.csv"
# comparison_group = "Schizophrenia"
# univariate_feature_set = "catch22"
# pairwise_feature_set = "pyspi14"
# pairwise_feature_file ="/headnode1/abry4213/data/UCLA_CNP/processed_data/UCLA_CNP_AROMA_2P_GMR_pyspi14_filtered_zscored.feather"
# noise_proc = "AROMA+2P+GMR"
# num_null_iters = 2

###############################################################################
# Main analysis
###############################################################################

run_pairwise_SVM(pairwise_feature_file=pairwise_feature_file,
                 SPI_directionality_file = SPI_directionality_file,
                 univariate_feature_set=univariate_feature_set,
                 pairwise_feature_set=pairwise_feature_set,
                       noise_proc = noise_proc,
                       dataset_ID=dataset_ID,
                       metadata_file=metadata_file,
                       comparison_to_control_group=comparison_group,
                       pydata_path=data_path + "processed_data/")

