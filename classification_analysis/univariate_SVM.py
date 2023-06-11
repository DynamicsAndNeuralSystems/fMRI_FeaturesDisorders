from sklearnex import patch_sklearn
patch_sklearn()
import sys

# append the sys.path list
sys.path.insert(0, 'helper_functions/classification/')
# sys.path.insert(0, 'github/fMRI_FeaturesDisorders/helper_functions/classification/')

from core_classification_functions import *
import argparse

# Command-line arguments to parse
parser = argparse.ArgumentParser(description='Process inputs for pairwise data preparation.')
parser.add_argument('--data_path', default="/headnode1/abry4213/data/UCLA_CNP/", dest='data_path')
parser.add_argument('--metadata_file', default="UCLA_CNP_sample_metadata.feather", dest='metadata_file')
parser.add_argument('--comparison_group', default="Schizophrenia", dest='comparison_group')
parser.add_argument('--univariate_feature_set', default='catch22', dest='univariate_feature_set')
parser.add_argument('--pairwise_feature_set', default='pyspi14', dest='pairwise_feature_set')
parser.add_argument('--univariate_feature_file', default="/headnode1/abry4213/data/UCLA_CNP/processed_data/UCLA_CNP_AROMA_2P_GMR_catch22_filtered_zscored.feather", dest='univariate_feature_file')
parser.add_argument('--noise_proc', dest='noise_proc')
parser.add_argument('--scaling_type', default="mixedsigmoid", dest='scaling_type')
parser.add_argument('--num_folds', default=10, dest='num_folds')
parser.add_argument('--num_repeats', default=10, dest='num_repeats')
parser.add_argument('--num_jobs', default=10, dest='num_jobs')
parser.add_argument('--dataset_ID', default="UCLA_CNP", dest='dataset_ID')

# Parse arguments
args = parser.parse_args()
data_path = args.data_path
metadata_file = args.metadata_file
comparison_group = args.comparison_group
univariate_feature_set = args.univariate_feature_set
pairwise_feature_set = args.pairwise_feature_set
univariate_feature_file = args.univariate_feature_file
noise_proc = args.noise_proc
scaling_type = args.scaling_type
num_folds = args.num_folds
num_repeats = args.num_repeats
num_jobs = args.num_jobs
dataset_ID = args.dataset_ID

# dataset_ID = "UCLA_CNP"
# data_path = "/headnode1/abry4213/data/UCLA_CNP/"
# metadata_file = "UCLA_CNP_sample_metadata.feather"
# comparison_group = "Schizophrenia"
# univariate_feature_set = "catch24"
# pairwise_feature_set = "pyspi14"
# univariate_feature_file ="/headnode1/abry4213/data/UCLA_CNP/processed_data/UCLA_CNP_AROMA_2P_GMR_catch24_filtered.feather"
# noise_proc = "AROMA+2P+GMR"
# num_folds = 10
# num_null_iters = 2
# num_repeats = 2
# num_jobs = 1
# scaling_type = "mixedsigmoid"

# Run the univariate main SVM
run_univariate_SVM(dataset_ID=dataset_ID,
                    data_path=data_path,
                    pydata_path=data_path + "processed_data/",
                    metadata_file=metadata_file,
                    comparison_to_control_group=comparison_group,
                    univariate_feature_set=univariate_feature_set, 
                    pairwise_feature_set=pairwise_feature_set,
                    univariate_feature_file=univariate_feature_file,
                    noise_proc=noise_proc,
                    scaling_type = scaling_type,
                    num_folds = int(num_folds),
                    num_repeats = int(num_repeats),
                    num_jobs = int(num_jobs))
