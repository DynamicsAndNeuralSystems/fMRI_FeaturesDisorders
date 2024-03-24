# Load modules
import pandas as pd
import argparse
import os
import sys
import numpy as np

# add path to classification analysis functions
sys.path.insert(0, './')
from core_classification_functions import *
current_path = os.getcwd()

# Command-line arguments to parse
parser = argparse.ArgumentParser(description='Process inputs for pairwise data preparation.')
parser.add_argument('--dataset_ID', dest='dataset_ID')
parser.add_argument('--disorder', dest='disorder')
parser.add_argument('--data_path', dest='data_path')
parser.add_argument('--main_analysis_type', dest='main_analysis_type', default="Univariate_catch25")
parser.add_argument('--classifier_type', dest='classifier_type', default='Linear_SVM_sklearn')

# Select classification function based on value of analysis_type
args = parser.parse_args()
dataset_ID = args.dataset_ID
disorder = args.disorder
main_analysis_type = args.main_analysis_type
classifier_type = args.classifier_type
data_path = args.data_path

# Concatenate null results
output_null_file = f"{data_path}/{dataset_ID}_{disorder}_{main_analysis_type}_{classifier_type}_null_balanced_accuracy_distributions.feather"

if not os.path.isfile(output_null_file):
    null_results = combine_null_results(data_path, dataset_ID, disorder)
    null_results.to_feather(output_null_file)