from sklearnex import patch_sklearn
patch_sklearn()

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

# Command-line arguments to parse
parser = argparse.ArgumentParser()
parser.add_argument('--num_jobs', dest='num_jobs', default=4)
args = parser.parse_args()
num_jobs = int(args.num_jobs)

# add path to classification analysis functions
sys.path.insert(0, './')
from core_classification_functions import *
current_path = os.getcwd()
from mixed_sigmoid_normalisation import MixedSigmoidScaler

data_path="/headnode1/abry4213/data/TS_feature_manuscript/"

# Load participants included
UCLA_CNP_subjects_to_keep = pd.read_feather(f"{data_path}/time_series_features/UCLA_CNP_filtered_sample_info_catch25_pyspi14.feather")
ABIDE_subjects_to_keep = pd.read_feather(f"{data_path}/time_series_features/ABIDE_filtered_sample_info_catch25_pyspi14.feather")

# Load metadata
UCLA_CNP_metadata = (pd.read_feather(f"{data_path}/input_data/UCLA_CNP_sample_metadata.feather")
                        .assign(Study = "UCLA_CNP")
                        .query("Sample_ID in @UCLA_CNP_subjects_to_keep.Sample_ID"))
ABIDE_metadata = (pd.read_feather(f"{data_path}/input_data/ABIDE_sample_metadata.feather")
                        .assign(Study = "ABIDE")
                        .query("Sample_ID in @ABIDE_subjects_to_keep.Sample_ID"))

# Load head movement 
UCLA_CNP_head_mvmt = (pd.read_table(f'{data_path}/movement_data/UCLA_CNP_Mean_FD_Power.txt', sep=',')
                      .assign(Study = "UCLA_CNP")
                        .query("Sample_ID in @UCLA_CNP_subjects_to_keep.Sample_ID"))
ABIDE_head_mvmt = (pd.read_table(f'{data_path}/movement_data/ABIDE_Mean_FD_Power.txt', sep=',', dtype={'Sample_ID': str,
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


# Analysis 1: Predicting diagnosis based on confound variables -- age, sex, and/or head motion
if not os.path.isfile(f"{data_path}/classification_results/confound_analysis/confounds_predicting_dx_balanced_accuracy_results.feather"):
    confounds_balanced_accuracy_list = []

    # Iterate over the four disorders
    for disorder in study_disorder_lookup.keys():
        study = study_disorder_lookup[disorder]

        class_labels = np.load(f"{data_path}/input_data/{study}_{disorder}_class_labels.npy")
        num_folds=10
        num_repeats=10

        age_feature = merged_metadata.query("Diagnosis in ['Control', @disorder] & Study == @study").Age.values.reshape(-1,1) 
        sex_feature = merged_metadata.query("Diagnosis in ['Control', @disorder] & Study == @study").Sex.values.reshape(-1,1) 
        head_mvmt_feature = merged_metadata.query("Diagnosis in ['Control', @disorder] & Study == @study").Mean_FD_Power.values.reshape(-1,1) 

        # Convert M to 0 and F to 1
        sex_feature = np.where(sex_feature == "M", 0, 1)

        # One dataset that includes age, sex, and head motion in one model
        all_confound_features = np.concatenate([age_feature, sex_feature, head_mvmt_feature], axis=1)

        SVM_model = svm.SVC(kernel='linear', class_weight='balanced', C=1)
        pipeline = Pipeline([('scaler', MixedSigmoidScaler(unit_variance=True)), 
                                ('model', SVM_model)])
        RepeatedStratifiedKFold_splitter = RepeatedStratifiedKFold(n_splits=num_folds, n_repeats=num_repeats, random_state=127) 

        analysis_type_results = []
        for analysis_type in ["Age", "Sex", "Mean_FD_Power", "All_Confounds"]:
            if analysis_type == "Age":
                feature_data = age_feature
            elif analysis_type == "Sex": 
                feature_data = sex_feature
            elif analysis_type == "Mean_FD_Power": 
                feature_data = head_mvmt_feature
            else: 
                feature_data = all_confound_features

            # Find balanced accuracy for dataset 
            confound_balanced_accuracy = cross_validate(pipeline, feature_data, class_labels, 
                                                        cv=RepeatedStratifiedKFold_splitter, scoring="balanced_accuracy", n_jobs=num_jobs,
                                                        return_estimator=False)['test_score']
            
            # Create dataframe
            confound_balanced_accuracy_df = pd.DataFrame({"Study" : study, 
                                                            "Disorder": disorder,
                                                            "Analysis_Type": analysis_type.replace("_", " "),
                                                            "Balanced_Accuracy": confound_balanced_accuracy})
            
            # Assign folds and repeats 
            confound_balanced_accuracy_df["Fold"] = confound_balanced_accuracy_df.index % num_folds
            confound_balanced_accuracy_df["Repeat"] = confound_balanced_accuracy_df.index // num_repeats

            # Append results to list 
            confounds_balanced_accuracy_list.append(confound_balanced_accuracy_df)

    # Concatenate results
    confounds_balanced_accuracy_results_all_folds = pd.concat(confounds_balanced_accuracy_list, axis=0)

    # Take average across folds per disorder 
    confounds_balanced_accuracy_results = (confounds_balanced_accuracy_results_all_folds
                                    .groupby(["Study", "Disorder", "Analysis_Type"], as_index=False)['Balanced_Accuracy']
                                    .agg(['mean', 'std'])
                                    .reset_index()
                                    .rename(columns={"mean": "Balanced_Accuracy", "std": "Balanced_Accuracy_SD"}))

    # Save to feather file
    confounds_balanced_accuracy_results_all_folds.reset_index().to_feather(f"{data_path}/classification_results/confound_analysis/confounds_predicting_dx_balanced_accuracy_results_all_folds.feather")
    confounds_balanced_accuracy_results.reset_index().to_feather(f"{data_path}/classification_results/confound_analysis/confounds_predicting_dx_balanced_accuracy_results.feather")

# Analysis 2: Predicting confound variables based on time-series features
if not os.path.isfile(f"{data_path}/classification_results/confound_analysis/time_series_features_predicting_confounds_r2_results.feather"):
    confounds_prediction_list = []

    # Iterate over the four disorders
    for disorder in study_disorder_lookup.keys():
        study = study_disorder_lookup[disorder]

        class_labels = np.load(f"../../data/input_data/{study}_{disorder}_class_labels.npy")
        num_folds=10
        num_repeats=10

        age_feature = merged_metadata.query("Diagnosis in ['Control', @disorder] & Study == @study").Age.values
        sex_feature = merged_metadata.query("Diagnosis in ['Control', @disorder] & Study == @study").Sex.values
        head_mvmt_feature = merged_metadata.query("Diagnosis in ['Control', @disorder] & Study == @study").Mean_FD_Power.values

        # Convert M to 0 and F to 1
        sex_feature = np.where(sex_feature == "M", 0, 1)

        # Apply support vector regression to predict age 
        RepeatedKFold_splitter = RepeatedKFold(n_splits=num_folds, n_repeats=num_repeats, random_state=127) 
        SVR_model = svm.SVR(kernel='linear', C=1)
        regression_pipeline = Pipeline([('scaler', MixedSigmoidScaler(unit_variance=True)), 
                                ('model', SVR_model)])
        SVM_model = svm.SVC(kernel='linear', class_weight='balanced', C=1)
        classification_pipeline = Pipeline([('scaler', MixedSigmoidScaler(unit_variance=True)), 
                                ('model', SVM_model)])

        # Define all models to test
        disorder_univariate_models = pd.read_table(f"/headnode1/abry4213/data/TS_feature_manuscript/time_series_features/processed_numpy_files/{study}_{disorder}_univariate_models.txt",header=None)
        disorder_pairwise_models = pd.read_table(f"/headnode1/abry4213/data/TS_feature_manuscript/time_series_features/processed_numpy_files/{study}_{disorder}_pairwise_models.txt",header=None)
        disorder_combined_univariate_pairwise_models  = pd.read_table(f"/headnode1/abry4213/data/TS_feature_manuscript/time_series_features/processed_numpy_files/{study}_{disorder}_combined_univariate_pairwise_models.txt",header=None)
        disorder_all_models = pd.concat([disorder_univariate_models, disorder_pairwise_models, disorder_combined_univariate_pairwise_models])
        disorder_all_models.columns = ["Model_Name"]

        for confound_variable in ["Age", "Sex", "Mean_FD_Power"]:
            if confound_variable == "Age":
                confound_data = age_feature
                pipeline = deepcopy(regression_pipeline)
            elif confound_variable == "Sex": 
                confound_data = sex_feature
                pipeline = deepcopy(classification_pipeline)
            elif confound_variable == "Mean_FD_Power": 
                confound_data = head_mvmt_feature
                pipeline = deepcopy(regression_pipeline)

            # Iterate over univariate models 
            for model_name in disorder_all_models["Model_Name"].tolist()[0:2]: 
            # Define analysis type
                if "ROI" in model_name:
                    Analysis_Type = "Brain_Region"
                elif "combo_catch25_features_all_regions" in model_name:
                    Analysis_Type = "Univariate_Combo"
                elif "combined_univariate_catch25_and_pyspi14" in model_name:
                    Analysis_Type = "SPI_Combo"
                elif "catch25_feature" in model_name:
                    Analysis_Type = "catch25_feature"
                else:
                    Analysis_Type = "pyspi14_SPI"

                if Analysis_Type=="Brain_Region":
                    grouping_var = model_name.split("_ROI_")[1]
                elif Analysis_Type=="Univariate_Combo":
                    grouping_var = "Combo"
                elif Analysis_Type == "SPI_Combo":
                    grouping_var = model_name.split("combined_univariate_catch25_and_pyspi14_SPI_")[1]
                elif Analysis_Type == "catch25_feature":
                    grouping_var = model_name.split("_catch25_feature_")[1]
                else:
                    grouping_var = model_name.split("_pyspi14_SPI_")[1]

                model_data = np.load(f"/headnode1/abry4213/data/TS_feature_manuscript/time_series_features/processed_numpy_files/{model_name}.npy")

                # Find balanced accuracy for dataset 
                model_confound_r2 = cross_validate(pipeline, model_data, confound_data, 
                                                            cv=RepeatedKFold_splitter, 
                                                            n_jobs=num_jobs, scoring='r2',
                                                            return_estimator=False)['test_score']

                # Create dataframe
                model_confound_r2_df = pd.DataFrame({"Study" : study, 
                                                        "Disorder": disorder,
                                                        "Confound_Variable": confound_variable,
                                                        "Analysis_Type": Analysis_Type,
                                                        "group_var": grouping_var,
                                                        "r2": np.mean(model_confound_r2),
                                                        "r2_SD": np.std(model_confound_r2)}, index=[0])
                # Append results to list 
                confounds_prediction_list.append(model_confound_r2_df)

    # Concatenate results
    time_series_features_predicting_confounds_r2_results = pd.concat(confounds_prediction_list, axis=0)

    # Save to feather file
    time_series_features_predicting_confounds_r2_results.reset_index().to_feather(f"{data_path}/classification_results/confound_analysis/time_series_features_predicting_confounds_r2_results.feather")