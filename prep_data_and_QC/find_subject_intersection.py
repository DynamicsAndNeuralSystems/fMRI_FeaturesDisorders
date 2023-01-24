
import argparse
import pandas as pd
import numpy as np
import random

# Command-line arguments to parse
parser = argparse.ArgumentParser(description='Process inputs for pairwise data preparation.')
parser.add_argument('--dataset_ID', default="UCLA_CNP", dest='dataset_ID')
parser.add_argument('--data_path', default="/headnode1/abry4213/data/UCLA_CNP/", dest='data_path')
parser.add_argument('--metadata_file', default="UCLA_CNP_sample_metadata.feather", dest='metadata_file')
parser.add_argument('--univariate_feature_set', default='catch22', dest='univariate_feature_set')
parser.add_argument('--pairwise_feature_set', default='pyspi14', dest='pairwise_feature_set')
parser.add_argument('--noise_proc', dest='noise_proc')

# Parse arguments
args = parser.parse_args()
data_path = args.data_path
dataset_ID = args.dataset_ID
metadata_file = args.metadata_file
univariate_feature_set = args.univariate_feature_set
pairwise_feature_set = args.pairwise_feature_set
noise_proc = args.noise_proc

noise_label = noise_proc.replace("+", "_")

# dataset_ID = "UCLA_CNP"
# data_path = "/headnode1/abry4213/data/UCLA_CNP/"
# metadata_file = "UCLA_CNP_sample_metadata.feather"
# univariate_feature_set = "catch22"
# pairwise_feature_set = "pyspi14"
# noise_proc = "AROMA+2P+GMR"

# Load info on subjects with univariate data
univariate_sample_info = pd.read_feather(f"{data_path}/processed_data/{dataset_ID}_filtered_sample_info_{noise_label}_{univariate_feature_set}.feather")

# Load info on subjects with pairwise data
