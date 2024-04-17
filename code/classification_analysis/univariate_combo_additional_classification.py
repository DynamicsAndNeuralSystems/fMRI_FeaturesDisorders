import os
import numpy as np
import nibabel as nib
import pandas as pd
import sys
import argparse

# Command-line arguments to parse
parser = argparse.ArgumentParser()
parser.add_argument('--num_jobs', dest='num_jobs', default=4)
args = parser.parse_args()
num_jobs = int(args.num_jobs)

# Uncomment the next two lines if you're running this on an HPC and want to speed up scikit-learn
# from sklearnex import patch_sklearn
# patch_sklearn()

# add path to classification analysis functions
sys.path.insert(0, 'code/classification_analysis/')
from core_classification_functions import *
current_path = os.getcwd()

data_path = "/headnode1/abry4213/data/TS_feature_manuscript"

# Load participants included
UCLA_CNP_subjects_to_keep = pd.read_feather(f"{data_path}/time_series_features/UCLA_CNP_filtered_sample_info_catch25_pyspi14.feather")
ABIDE_subjects_to_keep = pd.read_feather(f"{data_path}/time_series_features/ABIDE_filtered_sample_info_catch25_pyspi14.feather")
merged_subjects_to_keep = pd.concat([UCLA_CNP_subjects_to_keep, ABIDE_subjects_to_keep], axis=0).Sample_ID.tolist()

# Load metadata
UCLA_CNP_metadata = (pd.read_feather(f"{data_path}/input_data/UCLA_CNP_sample_metadata.feather")
                        .assign(Study = "UCLA_CNP")
                        .query("Sample_ID in @UCLA_CNP_subjects_to_keep.Sample_ID"))
ABIDE_metadata = (pd.read_feather(f"{data_path}/input_data/ABIDE_sample_metadata.feather")
                        .assign(Study = "ABIDE")
                        .query("Sample_ID in @ABIDE_subjects_to_keep.Sample_ID"))

# Load head movement 
UCLA_CNP_head_mvmt = (pd.read_table(f"{data_path}/movement_data/UCLA_CNP_Mean_FD_Power.txt", sep=',')
                      .assign(Study = "UCLA_CNP")
                        .query("Sample_ID in @UCLA_CNP_subjects_to_keep.Sample_ID"))
ABIDE_head_mvmt = (pd.read_table(f"{data_path}/movement_data/ABIDE_Mean_FD_Power.txt", sep=',', dtype={'Sample_ID': str,
                                                                                              'Mean_FD_Power': float})
                   .assign(Study = "ABIDE")
                        .query("Sample_ID in @ABIDE_subjects_to_keep.Sample_ID"))

# Merge metadata + head movement
merged_metadata = pd.concat([UCLA_CNP_metadata, ABIDE_metadata], axis=0).merge(pd.concat([UCLA_CNP_head_mvmt, ABIDE_head_mvmt], axis=0))

# Study/disorder lookup table
study_disorder_lookup = {'SCZ': 'UCLA_CNP', 
                          'BP': 'UCLA_CNP', 
                          'ADHD': 'UCLA_CNP', 
                          'ASD': 'ABIDE'}

# Load first 25 principal components for univariate region times feature combination matrices
univariate_combo_first25_PCs = pd.read_feather(f"{data_path}/time_series_features/univariate_combo_first25_PCs.feather")

# Load univariate time-series feature data for the two datasets
UCLA_CNP_univariate_features = pd.read_feather(f"{data_path}/time_series_features/UCLA_CNP_catch25_filtered.feather").merge(UCLA_CNP_metadata)
ABIDE_univariate_features = pd.read_feather(f"{data_path}/time_series_features/ABIDE_catch25_filtered.feather").merge(ABIDE_metadata)
merged_univariate_features = pd.concat([UCLA_CNP_univariate_features, ABIDE_univariate_features], axis=0)

# Prepare combination data for each disorder
study = "UCLA_CNP" 
disorder = "SCZ" 
disorder_catch25_data = merged_univariate_features.query("Study==@study & Diagnosis in ['Control', @disorder]")

# Lastly, combine feature names and region names into one matrix
disorder_catch25_data["Combo_Feature"] = disorder_catch25_data.Brain_Region + "_" + disorder_catch25_data.names
combo_data = disorder_catch25_data.drop(["Brain_Region", "names"], axis=1)

# Pivot from long to wide
combo_data_wide = combo_data.pivot(index=["Sample_ID", "Diagnosis"],
                                columns = "Combo_Feature",
                                values = "values")

if not os.path.isfile(f"{data_path}/classification_results/robustness_analysis/All_univariate_combo_first_25_PCs_classification_res.feather"):
    main_PC_model = svm.SVC(kernel="linear", C=1, class_weight="balanced")
    main_PC_pipe = Pipeline([('scaler', MixedSigmoidScaler(unit_variance=True)),
                        ('model', main_PC_model)])
    num_folds = 10
    num_repeats = 10
    RepeatedStratifiedKFold_splitter = RepeatedStratifiedKFold(n_splits=num_folds, n_repeats=num_repeats, random_state=127)
    num_null_iters = 0

    # Define scorers
    scorers = [make_scorer(balanced_accuracy_score)]
    scoring_names = ["Balanced_Accuracy"]

    all_univariate_combo_first_25_PCs_classification_res_list = []

    # Iterate over the four disorders
    for disorder in study_disorder_lookup.keys():
        study = study_disorder_lookup[disorder]

        class_labels = np.load(f"{data_path}/input_data/{study}_{disorder}_class_labels.npy", allow_pickle=True).tolist()
        sample_IDs = np.load(f"{data_path}/input_data/{study}_{disorder}_sample_IDs.npy", allow_pickle=True).tolist()
        num_folds=10
        num_repeats=10

        # Start with SCZ 
        disorder_univariate_combo_first25_PCs = univariate_combo_first25_PCs.query("Study==@study & Disorder == @disorder & Sample_ID in @merged_subjects_to_keep")
        disorder_univariate_combo_first25_PCs_features_only = disorder_univariate_combo_first25_PCs.drop(["Sample_ID", "Diagnosis", "Disorder", "Study", "index"], axis=1)

        main_classification_res, _, _ = run_k_fold_classifier_for_feature(feature_data = disorder_univariate_combo_first25_PCs_features_only, 
                                                                                                pipe = main_PC_pipe,
                                                                                                CV_splitter = RepeatedStratifiedKFold_splitter,
                                                                                                class_labels=class_labels,
                                                                                                sample_IDs = sample_IDs,
                                                                                                scorers=scorers,
                                                                                                scoring_names=scoring_names,
                                                                                                num_null_iters=num_null_iters,
                                                                                                num_folds = 10,
                                                                                                num_repeats = 10,
                                                                                                num_jobs = num_jobs)

        # Assign key details to dataframes
        main_classification_res["Study"] = study
        main_classification_res["Disorder"] = disorder
        main_classification_res["Analysis_Type"] = "Univariate_Combo_25_PCs"
        main_classification_res["group_var"] = "Univariate_Combo_25_PCs"
        main_classification_res["Classifier_Type"] = "Linear_SVM_sklearn"

        # Append results to dataframe 
        all_univariate_combo_first_25_PCs_classification_res_list.append(main_classification_res)

    # Concatenate results
    all_univariate_combo_first_25_PCs_classification_res = pd.concat(all_univariate_combo_first_25_PCs_classification_res_list, axis=0).reset_index()

    # Save to feather file
    all_univariate_combo_first_25_PCs_classification_res.to_feather(f"{data_path}/classification_results/robustness_analysis/All_univariate_combo_first_25_PCs_classification_res.feather")

if not os.path.isfile(f"{data_path}/classification_results/robustness_analysis/All_univariate_combo_L1_regularized_linear_SVM_classification_res.feather"):
    # L1 (LASSO) regression comparison
    main_LASSO_model = svm.LinearSVC(penalty="l1", loss="squared_hinge", dual=False, C=1, class_weight="balanced")
    main_LASSO_pipe = Pipeline([('scaler', MixedSigmoidScaler(unit_variance=True)),
                        ('model', main_LASSO_model)])
    num_folds = 10
    num_repeats = 10
    RepeatedStratifiedKFold_splitter = RepeatedStratifiedKFold(n_splits=num_folds, n_repeats=num_repeats, random_state=127)
    num_null_iters = 0

    # Define scorers
    scorers = [make_scorer(balanced_accuracy_score)]
    scoring_names = ["Balanced_Accuracy"]

    all_univariate_combo_LASSO_classification_res_list = []

    # Iterate over the four disorders
    for disorder in study_disorder_lookup.keys():
        study = study_disorder_lookup[disorder]

        class_labels = np.load(f"{data_path}/input_data/{study}_{disorder}_class_labels.npy", allow_pickle=True).tolist()
        sample_IDs = np.load(f"{data_path}/input_data/{study}_{disorder}_sample_IDs.npy", allow_pickle=True).tolist()
        num_folds=10
        num_repeats=10

        # Subset combo data to given disorder 
        disorder_catch25_data = merged_univariate_features.query("Study==@study & Diagnosis in ['Control', @disorder]")

        # Lastly, combine feature names and region names into one matrix
        disorder_catch25_data["Combo_Feature"] = disorder_catch25_data.Brain_Region + "_" + disorder_catch25_data.names
        combo_data = disorder_catch25_data.drop(["Brain_Region", "names"], axis=1)

        # Pivot from long to wide
        combo_data_wide = combo_data.pivot(index=["Sample_ID", "Diagnosis"],
                                        columns = "Combo_Feature",
                                        values = "values")

        main_classification_res, _, _ = run_k_fold_classifier_for_feature(feature_data = combo_data_wide, 
                                                                                                pipe = main_LASSO_pipe,
                                                                                                CV_splitter = RepeatedStratifiedKFold_splitter,
                                                                                                class_labels=class_labels,
                                                                                                sample_IDs = sample_IDs,
                                                                                                scorers=scorers,
                                                                                                scoring_names=scoring_names,
                                                                                                num_null_iters=num_null_iters,
                                                                                                num_folds = 10,
                                                                                                num_repeats = 10,
                                                                                                num_jobs = num_jobs)

        # Assign key details to dataframes
        main_classification_res["Study"] = study
        main_classification_res["Disorder"] = disorder
        main_classification_res["Analysis_Type"] = "Univariate_Combo_LASSO"
        main_classification_res["group_var"] = "Univariate_Combo_LASSO"
        main_classification_res["Classifier_Type"] = "Linear_SVM_L1_Regularized"

        # Append results to dataframe 
        all_univariate_combo_LASSO_classification_res_list.append(main_classification_res)

    # Concatenate results
    all_univariate_combo_LASSO_classification_res = pd.concat(all_univariate_combo_LASSO_classification_res_list, axis=0).reset_index()

    # Save to feather file
    all_univariate_combo_LASSO_classification_res.to_feather(f"{data_path}/classification_results/robustness_analysis/All_univariate_combo_L1_regularized_linear_SVM_classification_res.feather")