from sklearnex import patch_sklearn
patch_sklearn()
import sys
import scipy.io
import pandas as pd
import os

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
parser.add_argument('--scaling_type', default="robust", dest='scaling_type')
parser.add_argument('--run_nulls', action='store_true')
parser.add_argument('--num_folds', default=10, dest='num_folds')
parser.add_argument('--num_null_iters', default=1000, dest='num_null_iters')
parser.add_argument('--num_repeats', default=10, dest='num_repeats')
parser.add_argument('--num_jobs', default=8, dest='num_jobs')
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
run_nulls = args.run_nulls
scaling_type = args.scaling_type
num_folds = args.num_folds
num_null_iters = args.num_null_iters
num_repeats = args.num_repeats
num_jobs = args.num_jobs
dataset_ID = args.dataset_ID

# dataset_ID = "UCLA_CNP"
# data_path = "/headnode1/abry4213/data/UCLA_CNP/"
# metadata_file = "UCLA_CNP_sample_metadata.feather"
# comparison_to_control_group = "Schizophrenia"
# univariate_feature_set = "catch22"
# pairwise_feature_set = "pyspi14"
# univariate_feature_file ="/headnode1/abry4213/data/UCLA_CNP/processed_data/UCLA_CNP_AROMA_2P_GMR_catch22_filtered.feather"
# noise_proc = "AROMA+2P+GMR"
# num_repeats = 10
# num_folds = 5
# run_nulls = False
# num_jobs = 1
# scaling_type = "robustsigmoid"

###########################################################################
# Function to run SVMs per brain region per dataset
###########################################################################



###########################################################################
# Load data
###########################################################################

noise_label = noise_proc.replace("+", "_")

# Load metadata
metadata = pd.read_feather(data_path + "study_metadata/" + metadata_file)

# Load in data containing subjects with both univariate and pairwise data available
samples_to_keep = pd.read_feather(f"{data_path}/processed_data/{dataset_ID}_filtered_sample_info_{noise_label}_{univariate_feature_set}_{pairwise_feature_set}.feather").Sample_ID.tolist()                                                                   

# Univariate feature data
univariate_feature_data = pd.read_feather(univariate_feature_file).merge(metadata, on='Sample_ID', how='left').drop(["Age", "Sex"],
                                                                        axis = 1)


# Filter univariate data by samples with both univariate and pairwise
# Filter by samples with univariate data available as well
univariate_feature_data = univariate_feature_data[univariate_feature_data.Sample_ID.isin(samples_to_keep)] 
brain_region_list = univariate_feature_data.Brain_Region.unique().tolist()

# Load full movement data
movement_data = pd.DataFrame(scipy.io.loadmat(f"{data_path}/movement_data/{dataset_ID}_all_FD.mat")["FD_mat"],
                             columns = ["Sample_ID", "Jenkinson", "Power", "VanDijk"])

movement_data_long = (movement_data
                      .explode(["Sample_ID"])
                      .explode(["Jenkinson", "Power", "VanDijk"])
                      .explode(["Jenkinson", "Power", "VanDijk"])
                      .query("Sample_ID.isin(@samples_to_keep)"))
movement_data_long['Frame_Number'] = movement_data_long.groupby(['Sample_ID']).cumcount()+1

# Load mean movement data
movement_data_mean = pd.read_csv(f"{data_path}/movement_data/{dataset_ID}_mFD.txt", 
                            header=None, names=["Sample_ID", "Jenkinson", "Power", "VanDijk"])
                            
# Cast sample ID as string type
movement_data_mean = movement_data_mean.astype({'Sample_ID':'string'})

# Filter by samples to keep
movement_data_mean = movement_data_mean.query("Sample_ID.isin(@samples_to_keep)")


###########################################################################
# Check if output feather file exists
###########################################################################

if os.path.isfile(f"{data_path}/processed_data/{dataset_ID}_{comparison_group}_Univariate_{univariate_feature_set}_{scaling_type}_scaler_movement_SVM_balanced_accuracy.feather"):
    exit()



###########################################################################
# Keeping all subjects -- 5-fold CV with 10 repeats
###########################################################################                                                                       

no_threshold_balanced_accuracy = run_univariate_SVM_with_movement(univariate_feature_data=univariate_feature_data,
                       threshold_type="No Movement Threshold",
                       comparison_to_control_group = comparison_group,
                       num_folds = int(num_folds),
                       scaling_type=scaling_type,
                       num_jobs = int(num_jobs),
                       num_repeats = int(num_repeats))


###########################################################################
# Lenient thresholding -- 5-fold CV with 10 repeats
###########################################################################

movement_data_lenient = movement_data_mean.query("Power < 0.55")
univariate_feature_data_lenient = univariate_feature_data[univariate_feature_data.Sample_ID.isin(movement_data_lenient.Sample_ID)]   


lenient_threshold_balanced_accuracy = run_univariate_SVM_with_movement(univariate_feature_data=univariate_feature_data_lenient,
                       threshold_type="Lenient Movement Threshold",
                       comparison_to_control_group = comparison_group,
                       num_folds = int(num_folds),
                       scaling_type=scaling_type,
                       num_jobs = int(num_jobs),
                       num_repeats = int(num_repeats))

###########################################################################
# Stringent thresholding -- 5-fold CV with 10 repeats
###########################################################################

movement_data_stringent = (movement_data_long
                           .groupby(["Sample_ID"])
                           .filter(lambda x: (x["Power"].max() < 5) & 
                                   (x["Power"].mean() < 0.25) &
                                   (((sum(x["Power"] > 0.2)) / len(x)) < 0.2))
                           .groupby(["Sample_ID"])["Power"]
                           .mean()
    ).reset_index()

univariate_feature_data_stringent = univariate_feature_data[univariate_feature_data.Sample_ID.isin(movement_data_stringent.Sample_ID)]   

stringent_threshold_balanced_accuracy = run_univariate_SVM_with_movement(univariate_feature_data=univariate_feature_data_stringent,
                       threshold_type="Stringent Movement Threshold",
                       comparison_to_control_group = comparison_group,
                       num_folds = int(num_folds),
                       scaling_type=scaling_type,
                       num_jobs = int(num_jobs),
                       num_repeats = int(num_repeats))


###########################################################################
# Combine and save results
###########################################################################
all_threshold_balanced_accuracy = pd.concat([no_threshold_balanced_accuracy,
                                             lenient_threshold_balanced_accuracy,
                                             stringent_threshold_balanced_accuracy]).reset_index()

# Write to a feather file
all_threshold_balanced_accuracy.to_feather(f"{data_path}/processed_data/{dataset_ID}_{comparison_group}_Univariate_{univariate_feature_set}_{scaling_type}_scaler_movement_SVM_balanced_accuracy.feather")


