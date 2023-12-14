# Load libraries
import numpy as np
import pandas as pd
from scipy import stats
import argparse
import os

# Command-line arguments to parse
parser = argparse.ArgumentParser(description='Process inputs for pairwise data preparation.')
parser.add_argument('--github_dir', default="/headnode1/abry4213/github/", dest='github_dir')
parser.add_argument('--data_path', default="/headnode1/abry4213/data/UCLA_CNP/", dest='data_path')
parser.add_argument('--noise_proc', dest='noise_proc')
parser.add_argument('--dataset_ID', default="UCLA_CNP", dest='dataset_ID')

# Parse arguments
args = parser.parse_args()
data_path = args.data_path
noise_proc = args.noise_proc
dataset_ID = args.dataset_ID
github_dir = args.github_dir
fmri_github_dir = github_dir + "fMRI_FeaturesDisorders/"

# github_dir = "/headnode1/abry4213/github/"
# fmri_github_dir="/headnode1/abry4213/github/fMRI_FeaturesDisorders/"
# data_path="/headnode1/abry4213/data/UCLA_CNP/"
# noise_proc="AROMA+2P+GMR"
# dataset_ID="UCLA_CNP"

# Iterate over subjects
noise_label = noise_proc.replace("+", "_")
# Define raw time-series data file
raw_TS_file_dir = data_path + "raw_data/time_series_files/" + noise_label + "/"
try:
    os.makedirs(data_path + "raw_data/numpy_files/" + noise_label, 
                exist_ok = True)
except:
    pass
for TS_file in os.listdir(raw_TS_file_dir):
    try:
        sample_ID = TS_file.replace("_TS.csv", "")
        TS_data = pd.read_csv(raw_TS_file_dir + sample_ID + "_TS.csv", header=None)
        # Convert to numpy array
        TS_array = TS_data.to_numpy()
        data_norm = np.apply_along_axis(stats.zscore, 0, TS_array)
        data_norm = np.transpose(data_norm)
    
        # Save numpy array to a numpy binary file
        np.save(f"{data_path}/raw_data/numpy_files/{noise_label}/{sample_ID}.npy",
                data_norm)
    except:
        pass