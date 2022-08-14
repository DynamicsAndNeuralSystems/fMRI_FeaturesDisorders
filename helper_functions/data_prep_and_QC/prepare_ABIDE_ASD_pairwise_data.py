# Load libraries
import numpy as np
import pandas as pd
import scipy.io
from scipy import stats
import argparse
import os
import os.path

# Command-line arguments to parse
parser = argparse.ArgumentParser(description='Process inputs for pairwise data preparation.')
parser.add_argument('--github_dir', default="/project/hctsa/annie/github/", dest='github_dir')
parser.add_argument('--data_path', default="/project/hctsa/annie/data/ABIDE_ASD/", dest='data_path')
parser.add_argument('--input_mat_file', default="", nargs="?", dest='input_mat_file')
parser.add_argument('--subject_csv', default="participants.csv", dest='subject_csv')
parser.add_argument('--noise_procs', default=["FC1000"], nargs='*', action='append', dest='noise_procs')
parser.add_argument("--brain_region_lookup", default="Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv", dest="brain_region_lookup")
parser.add_argument("--parcellation_name", default="harvard_oxford_cort_prob_2mm", dest="parcellation_name")
parser.add_argument('--dataset_ID', default="ABIDE_ASD", dest='dataset_ID')

# Parse arguments
args = parser.parse_args()
data_path = args.data_path
input_mat_file = args.input_mat_file
subject_csv = args.subject_csv
noise_procs = args.noise_procs
parcellation_name = args.parcellation_name
brain_region_lookup = args.brain_region_lookup
dataset_ID = args.dataset_ID
github_dir = args.github_dir

fmri_github_dir = github_dir + "fMRI_FeaturesDisorders/"

# fmri_github_dir="/media/sf_Shared_Folder/github/fMRI_FeaturesDisorders/"
# data_path="/media/sf_Shared_Folder/PhD_work/data/ABIDE_ASD/"
# input_mat_file=""
# subject_csv='participants.csv'
# noise_procs=["FC1000"]
# parcellation_name = "harvard_oxford_cort_prob_2mm"
# brain_region_lookup = "Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv"
# dataset_ID="ABIDE_ASD"

# Input data
input_data_path = data_path + parcellation_name + "/"
# Define output data directory
pydata_path = data_path + "pydata/"
try:
    os.mkdir(pydata_path)
except OSError as error:
    pass

noise_label = noise_procs[0].replace("+", "_")

# Iterate over subjects and save data to numpy .npy files
for sample_ID in os.listdir(input_data_path):
    sample_CSV_file = data_path + parcellation_name + "/" + sample_ID + "/run_1/" + sample_ID + "_task-Rest_confounds.csv"
    sample_TS_data = pd.read_csv(sample_CSV_file, header=None).to_numpy()
    
    num_times, num_regions = sample_TS_data.shape
    print(f"Sample {sample_ID} {noise_label}: {num_regions} regions, {num_times} time points\n")

    # Z-score the time-series data
    data_norm = np.apply_along_axis(stats.zscore, 0, sample_TS_data)
    data_norm = np.transpose(data_norm)
    
    # Save region by timepoint data to its own .npy file
    if not os.path.exists(f"{pydata_path}/{noise_label}/{sample_ID}.npy"):
        print(f"Saving .npy data for {sample_ID}")
        np.save(f"{pydata_path}/{noise_label}/{sample_ID}.npy", data_norm)

