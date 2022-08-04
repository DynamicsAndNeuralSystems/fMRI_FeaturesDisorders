# Load libraries
import numpy as np
import pandas as pd
import scipy.io
from scipy import stats
import argparse
import os

# Command-line arguments to parse
parser = argparse.ArgumentParser(description='Process inputs for pairwise data preparation.')
parser.add_argument('--github_dir', default="/project/hctsa/annie/github/fMRI_FeaturesDisorders/", dest='github_dir')
parser.add_argument('--data_path', default="/project/hctsa/annie/data/UCLA_Schizophrenia/", dest='data_path')
parser.add_argument('--input_mat_file', default="", nargs="?", dest='input_mat_file')
parser.add_argument('--subject_csv', default="participants.csv", dest='subject_csv')
parser.add_argument('--noise_procs', default=["AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER"], nargs='*', action='append', dest='noise_procs')
parser.add_argument('--parcellation_name', default="harvard_oxford_cort_prob_2mm", dest='parcellation_name', nargs='?')
parser.add_argument('--brain_region_lookup', default="Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv", dest='brain_region_lookup', nargs='?')
parser.add_argument('--dataset_ID', default="UCLA_Schizophrenia", dest='dataset_ID')

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

# github_dir="/media/sf_Shared_Folder/github/fMRI_FeaturesDisorders/"
# data_path="/media/sf_Shared_Folder/PhD_work/data/UCLA_Schizophrenia/"
# input_mat_file="new/UCLA_time_series_four_groups.mat"
# subject_csv='participants.csv'
# noise_procs=["AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER"]
# dataset_ID="UCLA_Schizophrenia"

# Define output data directory
pydata_path = data_path + "pydata/"
try:
    os.mkdir(pydata_path)
except OSError as error:
    pass

# Helper function
def flatten(t):
    return [item for sublist in t for item in sublist]

# Load subject data
sample_info = pd.read_csv(data_path + subject_csv, sep=",")

# Load TS data
fMRI_data = scipy.io.loadmat(data_path+input_mat_file)

# Load preprocessing data
noise_proc_info = flatten([x.tolist() for x in fMRI_data["noiseOptions"].tolist()[0]])

# Brain region info
Brain_Region_info = fMRI_data["StructNames"].tolist()
Brain_Region_info = [Brain_Region.replace(" ", "") for Brain_Region in Brain_Region_info]
Brain_Region_df = pd.DataFrame(data={"Brain_Region": Brain_Region_info, "Index": list(range(1, len(Brain_Region_info) + 1))})

# Write brain region index + name to a CSV
Brain_Region_df.to_csv(pydata_path + "Brain_Region_info.csv", sep=',', index=False)

# Load sample IDs
sample_IDs = flatten(flatten(fMRI_data["subject_list"].tolist()))

# Read in time-series data and split by noise processing method (4th dimension)
TS_full = fMRI_data["time_series"]
TS_data = np.split(TS_full, 3, axis=3)

# The following noise-processing methods are stored in the following slots of TS_data:
# AROMA+2P = 0
# AROMA+2P+GMR = 1
# AROMA+2P+DiCER = 2
noise_proc_indices = {0: "AROMA+2P", 1:"AROMA+2P+GMR", 2:"AROMA+2P+DiCER"}

def array_to_npy(i, noise_proc, split_data, sample_IDs, pydata_path):
    # Subset data to ith sample
    data_array = split_data[i]
    # Subset data to just 2 dimensions: brain region by timepoint
    array_2d = data_array[:,:,0]
    # Find sample ID
    sample_ID = sample_IDs[i]
    num_times, num_regions = array_2d.shape
    print(f"Sample {sample_ID} {noise_proc}: {num_regions} regions, {num_times} time points\n")
    # Z-score the time-series data
    data_norm = np.apply_along_axis(stats.zscore, 0, array_2d)
    data_norm = np.transpose(data_norm)
    # Rename noise processing method to have underscores
    noise_label = noise_proc.replace("+", "_")
    # Save region by timepoint data to its own .npy file
    np.save(f"{pydata_path}/{noise_label}/{sample_ID}.npy", data_norm)

for index in noise_proc_indices:
    noise_proc = noise_proc_indices[index]
    noise_label = noise_proc.replace("+", "_")
    
    # Make output directory
    try:
        os.mkdir(pydata_path + noise_label)
    except OSError as error:
        pass

    # Subset time-series from this noise-processing method
    noise_proc_time_series = TS_data[index][:,:,:,0]

    # Split along axis=2
    TS_data_split = np.split(noise_proc_time_series, 260, axis=2)

    # Save a .npy file for each sample + noise-processing method
    [array_to_npy(i=i, noise_proc = noise_proc, split_data = TS_data_split, 
                  sample_IDs = sample_IDs, pydata_path = pydata_path) 
                  for i in list(range(len(TS_data_split)))]

