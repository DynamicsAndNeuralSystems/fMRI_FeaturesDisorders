# Load modules
import numpy as np
import pandas as pd
import dill
import argparse
import os
import pyarrow.feather as feather

# Command-line arguments to parse
parser = argparse.ArgumentParser(description='Process inputs for pairwise data preparation.')
parser.add_argument('--data_path', default="/headnode1/abry4213/data/UCLA_CNP_ABIDE_ASD/", dest='data_path')
parser.add_argument('--dataset_ID', default="UCLA_CNP", dest='dataset_ID')
parser.add_argument('--pairwise_feature_set', default="pyspi14", dest='pairwise_feature_set')
parser.add_argument('--univariate_feature_set', default="catch22", dest='univariate_feature_set')
parser.add_argument('--brain_region_lookup', default="UCLA_CNP_Brain_Region_Lookup.feather", dest='brain_region_lookup')
parser.add_argument('--noise_proc', dest='noise_proc')

# Parse arguments
args = parser.parse_args()
data_path = args.data_path
noise_proc = args.noise_proc
dataset_ID = args.dataset_ID
pairwise_feature_set = args.pairwise_feature_set
univariate_feature_set = args.univariate_feature_set
brain_region_lookup = args.brain_region_lookup

# data_path = "/headnode1/abry4213/data/UCLA_CNP/"
# dataset_ID = "UCLA_CNP"
# noise_proc="AROMA+2P+GMR"
# pairwise_feature_set = "pyspi14"
# univariate_feature_set = "catch22"
# brain_region_lookup = "UCLA_CNP_Brain_Region_Lookup.feather"

proc_data_path = data_path + "processed_data/"

        
def intersection_univariate_pairwise(proc_data_path, dataset_ID, noise_proc, univariate_feature_set, pairwise_feature_set):
    # Set noise label for file paths
    noise_label = noise_proc.replace("+", "_")

    # Load in data on samples with univariate feature data
    univariate_data_to_keep = pd.read_feather(f"{proc_data_path}/{dataset_ID}_filtered_sample_info_{noise_label}_{univariate_feature_set}.feather")
    
    # Load in data on samples with pairwise feature data
    pairwise_data_to_keep = pd.read_feather(f"{proc_data_path}/{dataset_ID}_filtered_sample_info_{noise_label}_{pairwise_feature_set}.feather")
    feather.write_feather(pairwise_data_to_keep,
                          f"{proc_data_path}/{dataset_ID}_filtered_sample_info_{noise_label}_{pairwise_feature_set}.feather",
                          version=1)
    
    # Merge the two datasets
    merged_sample_info = pd.merge(univariate_data_to_keep, pairwise_data_to_keep, how="inner")
    feather.write_feather(merged_sample_info,
                          f"{data_path}/processed_data/{dataset_ID}_filtered_sample_info_{noise_label}_{univariate_feature_set}_{pairwise_feature_set}.feather",
                          version=1)
    
intersection_univariate_pairwise(proc_data_path = proc_data_path, 
                  dataset_ID = dataset_ID, 
                  noise_proc = noise_proc, 
                  univariate_feature_set = univariate_feature_set, 
                  pairwise_feature_set = pairwise_feature_set)
