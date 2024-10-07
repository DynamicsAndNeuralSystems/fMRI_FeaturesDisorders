import os
import sys
import numpy as np
import nibabel as nib
import pandas as pd
from sklearn import svm
from sklearn.pipeline import Pipeline
import os.path
from sklearn.model_selection import RepeatedStratifiedKFold, cross_validate, permutation_test_score, StratifiedKFold, GridSearchCV, RepeatedKFold
from sklearn.ensemble import RandomForestClassifier,GradientBoostingClassifier
from sklearn.metrics import balanced_accuracy_score,make_scorer
from copy import deepcopy
import argparse
from scipy.stats import pearsonr

from core_classification_functions import *
current_path = os.getcwd()
from mixed_sigmoid_normalisation import MixedSigmoidScaler
data_path="/headnode1/abry4213/data/TS_feature_manuscript/"

# Load metadata
UCLA_CNP_metadata = pd.read_feather(f"{data_path}/input_data/UCLA_CNP_sample_metadata_filtered.feather")
ABIDE_metadata = pd.read_feather(f"{data_path}/input_data/ABIDE_sample_metadata_filtered.feather").assign(Study="ABIDE")

# Load time series
# Load filtered subject list
UCLA_CNP_filtered_subjects = pd.read_feather(f"{data_path}/time_series_features/UCLA_CNP_filtered_sample_info_catch25_pyspi14.feather")
# Load time-series for UCLA CNP
UCLA_CNP_fMRI_TS = pd.read_feather(f"{data_path}/input_data/UCLA_CNP_fMRI_TS.feather").query("Sample_ID in @UCLA_CNP_filtered_subjects.Sample_ID")

# ABIDE
ABIDE_filtered_subjects = pd.read_feather(f"{data_path}/time_series_features/ABIDE_filtered_sample_info_catch25_pyspi14.feather")
ABIDE_fMRI_TS = pd.read_feather(f"{data_path}/input_data/ABIDE_fMRI_TS.feather").query("Sample_ID in @ABIDE_filtered_subjects.Sample_ID")

# Compute brain-wide correlations
region_corr_res_list = []

# UCLA CNP
if not os.path.isfile(f"{data_path}/time_series_features/All_glocal_coupling.feather"):
    for subject_ID in UCLA_CNP_fMRI_TS['Sample_ID'].unique().tolist():
        this_subject_data = UCLA_CNP_fMRI_TS.query("Sample_ID == @subject_ID")
        this_subject_metadata = UCLA_CNP_metadata.query("Sample_ID == @subject_ID")

        for this_region in this_subject_data['Brain_Region'].unique():
            this_region_data = this_subject_data.query("Brain_Region == @this_region")
            this_rest_of_brain_data = this_subject_data.query("Brain_Region != @this_region").groupby('timepoint')['values'].mean().reset_index()

            # Compute pearson correlation between the region and the rest of the brain
            this_region_corr = pearsonr(this_region_data['values'], this_rest_of_brain_data['values'])

            this_region_res = (pd.DataFrame({"Sample_ID": subject_ID,
                                            "Brain_Region": this_region,
                                            "Pearson_Correlation": this_region_corr[0]}, index=[0])
                                            .merge(this_subject_metadata, on="Sample_ID")
                                            .assign(Study = "UCLA_CNP"))
            region_corr_res_list.append(this_region_res)

    # ABIDE
    for subject_ID in ABIDE_fMRI_TS['Sample_ID'].unique().tolist():
        this_subject_data = ABIDE_fMRI_TS.query("Sample_ID == @subject_ID")
        this_subject_metadata = ABIDE_metadata.query("Sample_ID == @subject_ID")

        for this_region in this_subject_data['Brain_Region'].unique():
            this_region_data = this_subject_data.query("Brain_Region == @this_region")
            this_rest_of_brain_data = this_subject_data.query("Brain_Region != @this_region").groupby('timepoint')['values'].mean().reset_index()

            # Compute pearson correlation between the region and the rest of the brain
            this_region_corr = pearsonr(this_region_data['values'], this_rest_of_brain_data['values'])

            this_region_res = (pd.DataFrame({"Sample_ID": subject_ID,
                                            "Brain_Region": this_region,
                                            "Pearson_Correlation": this_region_corr[0]}, index=[0])
                                            .merge(this_subject_metadata, on="Sample_ID")
                                            .assign(Study = "ABIDE"))
            region_corr_res_list.append(this_region_res)

    region_corr_res = pd.concat(region_corr_res_list).reset_index(drop=True)
    region_corr_res.to_feather(f"{data_path}/time_series_features/All_glocal_coupling.feather")
