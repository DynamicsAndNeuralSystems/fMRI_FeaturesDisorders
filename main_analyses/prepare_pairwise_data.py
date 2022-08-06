# Load libraries
import numpy as np
import pandas as pd
import scipy.io
from scipy import stats
import argparse
import os

# Command-line arguments to parse
parser = argparse.ArgumentParser(description='Process inputs for pairwise data preparation.')
parser.add_argument('--project_path', default="/project/hctsa/annie/", dest='project_path')
parser.add_argument('--github_dir', default="/project/hctsa/annie/github/", dest='github_dir')
parser.add_argument('--data_path', default="/project/hctsa/annie/data/UCLA_Schizophrenia/", dest='data_path')
parser.add_argument('--input_mat_file', default="", nargs="?", dest='input_mat_file')
parser.add_argument('--subject_csv', default="participants.csv", dest='subject_csv')
parser.add_argument('--pairwise_feature_set', default="pyspi19", dest='pairwise_feature_set')
parser.add_argument('--parcellation_name', default="harvard_oxford_cort_prob_2mm", dest='parcellation_name', nargs='?')
parser.add_argument('--brain_region_lookup', default="Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv", dest='brain_region_lookup', nargs='?')
parser.add_argument('--noise_procs', default=["AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER"], nargs='*', dest='noise_procs')
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

# Define pyspi script directory
pyspi_script_dir = github_dir + "pyspi-distribute/"

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
noise_procs_cl = ' '.join(noise_procs)

# Prepare pyspi data
prepare_pyspi_cmd=f"python3 {github_dir}/fMRI_FeaturesDisorders/data_prep_and_QC/prepare_{dataset_ID}_pairwise_data.py --github_dir {github_dir} --data_path {data_path} --input_mat_file {input_mat_file} --subject_csv {subject_csv} --noise_procs {noise_procs_cl} --brain_region_lookup {brain_region_lookup} --parcellation_name {parcellation_name} --dataset_ID {dataset_ID}"
os.system(prepare_pyspi_cmd)  

# Run pyspi
for noise_proc in noise_procs:
    noise_label = noise_proc.replace("+", "_")
    run_pyspi_cmd=f"python $pyspi_script_dir/distribute_jobs.py --data_dir {project_path}/data/{dataset_ID}/pydata/{noise_label}/ --compute_file {pyspi_script_dir}/pyspi_compute.py --template_pbs_file {pyspi_script_dir}/template.pbs --pyspi_config {config_file} --sample_yaml {project_path}/data/{dataset_ID}/pydata/{noise_label}/sample.yaml --pbs_notify a  --email abry4213@uni.sydney.edu.au --walltime_hrs 2 --cpu 2 --mem 8 --table_only"
    os.system(run_pyspi_cmd)

# Write calc.table pkl files to CSVs
pkl_to_csv_cmd=f"python3 {github_dir}/fMRI_FeaturesDisorders/data_prep_and_QC/pyspi_pickle_to_csv --github_dir {github_dir} --data_path {data_path} --noise_procs {noise_procs_cl} --brain_region_lookup {brain_region_lookup} --parcellation_name {parcellation_name} --dataset_ID {dataset_ID}"
os.system(pkl_to_csv_cmd)  