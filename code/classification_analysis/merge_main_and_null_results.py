# Load modules
import pandas as pd
import argparse
import os
import sys
import numpy as np
import re

# add path to classification analysis functions
sys.path.insert(0, './')
from core_classification_functions import *
current_path = os.getcwd()

# Command-line arguments to parse
parser = argparse.ArgumentParser(description='Process inputs for pairwise data preparation.')
parser.add_argument('--dataset_ID', dest='dataset_ID')
parser.add_argument('--disorder', dest='disorder')
parser.add_argument('--data_path', dest='data_path')
parser.add_argument('--classifier_type', dest='classifier_type', default='Linear_SVM_sklearn')

# Select classification function based on value of analysis_type
args = parser.parse_args()
dataset_ID = args.dataset_ID
disorder = args.disorder
classifier_type = args.classifier_type
data_path = args.data_path

analysis_type_kwds = {"univariate": "Brain_Region|catch25_feature|Univariate_Combo",
                      "pairwise": f"{disorder}_pyspi14_SPI",
                      "combined_univariate_pairwise": "combined_univariate_catch25_and_pyspi14"}

for main_analysis_type in ["univariate"]:
# for main_analysis_type in ["univariate", "pairwise", "combined_univariate_pairwise"]:

    # Concatenate main results
    output_main_file = f"{data_path}/balanced_accuracy/{dataset_ID}_{disorder}_{main_analysis_type}_{classifier_type}_balanced_accuracy_all_folds.feather"

    # Find files to combine that match the corresponding list of patterns in analysis_type_kwds
    this_analysis_kwds = analysis_type_kwds[main_analysis_type]
    files = [f for f in os.listdir(f"{data_path}/balanced_accuracy/{dataset_ID}_{disorder}/") if os.path.isfile(os.path.join(f"{data_path}/balanced_accuracy/{dataset_ID}_{disorder}/", f))]
    main_files_for_this_analysis_type = [f"{data_path}/balanced_accuracy/{dataset_ID}_{disorder}/{f}" for f in files if re.search(this_analysis_kwds, f)]

    main_results = combine_main_results(files_to_merge=main_files_for_this_analysis_type, dataset_ID=dataset_ID, disorder=disorder, average_across_folds=False)
    main_results.to_feather(output_main_file)

    output_main_file_averaged = f"{data_path}/balanced_accuracy/{dataset_ID}_{disorder}_{main_analysis_type}_{classifier_type}_balanced_accuracy.feather"
    main_results_averaged = (main_results
                                .groupby(["group_var", "Classifier_Type", "Analysis_Type", "Disorder"], as_index=False)["Balanced_Accuracy"]
                                .agg(["mean", "std"])
                                .reset_index()
                                .rename(columns={"mean": "Balanced_Accuracy", "std": "Balanced_Accuracy_SD"}))
    main_results_averaged.to_feather(output_main_file_averaged)

    # Concatenate null results
    output_null_file = f"{data_path}/null_results/{dataset_ID}_{disorder}_{main_analysis_type}_{classifier_type}_null_balanced_accuracy_distributions.feather"

    files = [f for f in os.listdir(f"{data_path}/null_results/{dataset_ID}_{disorder}/") if os.path.isfile(os.path.join(f"{data_path}/null_results/{dataset_ID}_{disorder}/", f))]
    null_files_for_this_analysis_type = [f"{data_path}/null_results/{dataset_ID}_{disorder}/{f}" for f in files if re.search(this_analysis_kwds, f)]

    null_results = combine_null_results(files_to_merge=null_files_for_this_analysis_type, dataset_ID=dataset_ID, disorder=disorder, num_null_iters=1000)
    null_results.to_feather(output_null_file)