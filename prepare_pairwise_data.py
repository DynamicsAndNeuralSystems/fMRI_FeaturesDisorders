# Load libraries
import numpy as np
import pandas as pd
import scipy.io
from scipy import stats
import argparse
import os

# Command-line arguments to parse
parser = argparse.ArgumentParser(description='Process inputs for pairwise data preparation.')
parser.add_argument('--project_path', default="/project/hctsa/annie/github/fMRI_FeaturesDisorders/", dest='project_path')
parser.add_argument('--github_dir', default="/project/hctsa/annie/", dest='github_dir')
parser.add_argument('--data_path', default="/project/hctsa/annie/data/UCLA_Schizophrenia/", dest='data_path')
parser.add_argument('--input_mat_file', default="", nargs="?", dest='input_mat_file')
parser.add_argument('--subject_csv', default="participants.csv", dest='subject_csv')
parser.add_argument('--pairwise_feature_set', default="pyspi19", dest='pairwise_feature_set')
parser.add_argument('--parcellation_name', default="harvard_oxford_cort_prob_2mm", dest='parcellation_name', nargs='?')
parser.add_argument('--brain_region_lookup', default="Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv", dest='brain_region_lookup', nargs='?')
parser.add_argument('--noise_procs', default=["AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER"], nargs='*', action='append', dest='noise_procs')
parser.add_argument('--dataset_ID', default="UCLA_Schizophrenia", dest='dataset_ID')

# Parse arguments
args = parser.parse_args()
project_path = args.project_path
github_dir = args.github_dir
data_path = args.data_path
input_mat_file = args.input_mat_file
subject_csv = args.subject_csv
pairwise_feature_set = args.pairwise_feature_set
parcellation_name = args.parcellation_name
brain_region_lookup = args.brain_region_lookup
noise_procs = args.noise_procs
dataset_ID = args.dataset_ID

# pairwise_feature_set = "pyspi19"
# subject_csv = "participants.csv"
# project_path = "/media/sf_Shared_Folder/github/"
# github_dir = "/media/sf_Shared_Folder/github/fMRI_FeaturesDisorders/"

# UCLA schizophrenia
# data_path = "/media/sf_Shared_Folder/PhD_work/data/UCLA_Schizophrenia/"
# dataset_ID = "UCLA_Schizophrenia"
# input_mat_file = "new/UCLA_time_series_four_groups.mat"
# noise_procs = ["AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER"]
# parcellation_name = "aparc+aseg"
# brain_region_lookup = ""

# ABIDE ASD
# data_path = "/media/sf_Shared_Folder/PhD_work/data/ABIDE_ASD/"
# dataset_ID = "ABIDE_ASD"
# input_mat_file = ""
# noise_procs = ["FC1000"]
# parcellation_name = "harvard_oxford_cort_prob_2mm"
# brain_region_lookup = "Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv"

# Prep noise-procs for command line
print("Noise procs:")
print(noise_procs)
noise_procs_cl = ' '.join(noise_procs)

cmd_to_execute=f"python3 {github_dir}/data_prep_and_QC/prepare_{dataset_ID}_pairwise_data.py --github_dir {github_dir} --data_path {data_path} --input_mat_file {input_mat_file} --subject_csv {subject_csv} --noise_procs {noise_procs_cl} --brain_region_lookup {brain_region_lookup} --parcellation_name {parcellation_name} --dataset_ID {dataset_ID}"

os.system(cmd_to_execute)  

