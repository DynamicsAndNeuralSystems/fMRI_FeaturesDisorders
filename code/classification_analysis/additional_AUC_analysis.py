from sklearnex import patch_sklearn
patch_sklearn()

import os
import sys
import numpy as np
import pandas as pd
import os.path
import argparse

# Command-line arguments to parse
parser = argparse.ArgumentParser(description='Process inputs for pairwise data preparation.')
parser.add_argument('--num_jobs', dest='num_jobs', default=4)

# Parse input arguments
args = parser.parse_args()
num_jobs = int(args.num_jobs)

data_path="/headnode1/abry4213/data/TS_feature_manuscript/"
sys.path.insert(0, './')
from core_classification_functions import *
from mixed_sigmoid_normalisation import MixedSigmoidScaler
from sklearn.metrics import roc_auc_score

# Load data
UCLA_CNP_subjects_to_keep = pd.read_feather(f"{data_path}/time_series_features/UCLA_CNP_filtered_sample_info_catch25_pyspi14.feather")
ABIDE_subjects_to_keep = pd.read_feather(f"{data_path}/time_series_features/ABIDE_filtered_sample_info_catch25_pyspi14.feather")

# Load metadata
UCLA_CNP_metadata = (pd.read_feather(f"{data_path}/input_data/UCLA_CNP_sample_metadata.feather")
                        .assign(Study = "UCLA_CNP")
                        .query("Sample_ID in @UCLA_CNP_subjects_to_keep.Sample_ID"))
ABIDE_metadata = (pd.read_feather(f"{data_path}/input_data/ABIDE_sample_metadata.feather")
                        .assign(Study = "ABIDE")
                        .query("Sample_ID in @ABIDE_subjects_to_keep.Sample_ID"))

# Load univariate time-series feature info
univariate_feature_info = pd.read_csv(f"{data_path}/feature_info/univariate_feature_info.csv")
pairwise_feature_info = pd.read_csv(f"{data_path}/feature_info/pairwise_feature_info.csv")

# Define parameters that you can change
univariate_feature_set = "catch25"
pairwise_feature_set = "pyspi14"

# Load SPI directionality information
SPI_directionality_data = pd.read_csv("SPI_Direction_Info.csv")
SPI_directionality_dict = dict(SPI_directionality_data.values)

# Load univariate time-series feature data for the two datasets
UCLA_CNP_univariate_features = pd.read_feather(f"{data_path}/time_series_features/UCLA_CNP_catch25_filtered.feather")
ABIDE_univariate_features = pd.read_feather(f"{data_path}/time_series_features/ABIDE_catch25_filtered.feather")

# Load pyspi14 data for UCLA CNP and ABIDE
UCLA_CNP_pyspi14 = pd.read_feather(f"{data_path}/time_series_features/UCLA_CNP_pyspi14_filtered.feather")
ABIDE_pyspi14 = pd.read_feather(f"{data_path}/time_series_features/ABIDE_pyspi14_filtered.feather")

model = svm.SVC(kernel="linear", C=1, class_weight="balanced", probability=True)
pipe = Pipeline([('scaler', MixedSigmoidScaler(unit_variance=True)),
                ('model', model)])
classifier_type = "Linear_SVM_sklearn"

# Define scorers
scorers = [make_scorer(roc_auc_score, needs_proba=True)]
scoring_names = ["AUC"]

# Classification parameters
num_folds = 10
num_repeats = 10
num_null_iters = 0
RepeatedStratifiedKFold_splitter = RepeatedStratifiedKFold(n_splits=num_folds, n_repeats=num_repeats, random_state=127)

# Start list 
auc_results_list = []

study_lookup_df = pd.DataFrame({"Disorder": ["SCZ", "BP", 'ADHD', "ASD"], "Study": ["UCLA_CNP", "UCLA_CNP", "UCLA_CNP", "ABIDE"]})

# Iterate over disorder and study in study_lookup_df
for disorder, dataset_ID in study_lookup_df.values:
    # Class labels and sample IDs
    class_labels = np.load(f"{data_path}/input_data/{dataset_ID}_{disorder}_class_labels.npy")
    sample_IDs = np.load(f"{data_path}/input_data/{dataset_ID}_{disorder}_sample_IDs.npy")

    # Iterate over model types
    for model_type in ['univariate', 'pairwise', 'combined_univariate_pairwise']:
        model_list = pd.read_table(f"{data_path}/time_series_features/processed_numpy_files/{dataset_ID}_{disorder}_{model_type}_models.txt", header=None)[0].tolist()

        # Iterate over models
        for model_name in model_list:
            print(f"Running {model_name}")
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

            # Find grouping_var
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

            output_file = f"{data_path}/classification_results/auc/{model_name}_AUC_results.feather"

            # Run classification if file doesn't exist already
            if not os.path.isfile(output_file):
                feature_data = np.load(f"{data_path}/time_series_features/processed_numpy_files/{model_name}.npy")
                # # Define splits
                # # Find splits
                # splits = list(RepeatedStratifiedKFold_splitter.split(feature_data, class_labels))
                # # Convert splits to a dataframe
                # splits_df = pd.DataFrame(splits, columns = ["Train", "Test"])
                # # Assign the fold and repeat numbers
                # splits_df["Fold"] = splits_df.index % num_folds
                # splits_df["Repeat"] = splits_df.index // num_repeats

                # fold_res_list = []

                # # Iterate over each row of fold_splits dataframe
                # for i, row in splits_df.iterrows():
                #     fold_num = row["Fold"]
                #     repeat_num = row["Repeat"]
                #     train_indices = row["Train"]
                #     test_indices = row["Test"]

                #     train_data = feature_data[train_indices]
                #     test_data = feature_data[test_indices]

                #     train_labels = class_labels[train_indices]
                #     test_labels = class_labels[test_indices]

                #     # Fit pipe to train_data
                #     loop_pipe = deepcopy(pipe)
                #     loop_pipe.fit(train_data, train_labels)
                #     test_preds = loop_pipe.predict(test_data)
                #     test_preds_prob = loop_pipe.predict_proba(test_data)

                #     # Figure out which column to keep: 
                #     if len(test_preds[test_preds==0]) == 0:
                #         # Find whichever column of test_preds_prob has the lower mean
                #         prob_col = np.argmin([np.mean(test_preds_prob[:,0]), np.mean(test_preds_prob[:,1])])
                #     elif len(test_preds[test_preds==1]) == 0:
                #         # Find whichever column of test_preds_prob has the higher mean
                #         prob_col = np.argmax([np.mean(test_preds_prob[:,0]), np.mean(test_preds_prob[:,1])])
                #     elif np.mean(test_preds_prob[test_preds==0]) < 0.5: 
                #         prob_col = 1
                #     else: 
                #         prob_col = 0
                #     test_preds_prob_data = test_preds_prob[:,prob_col]

                #     # Compute AUC
                #     auc = roc_auc_score(test_labels, test_preds_prob_data)

                #     this_loop_res = pd.DataFrame({"Fold": fold_num, "Repeat": repeat_num, "AUC": auc}, index=[0])
                #     fold_res_list.append(this_loop_res)

                # main_classification_res = pd.concat(fold_res_list)

                main_classification_res, _, _ = run_k_fold_classifier_for_feature(feature_data = feature_data, 
                                                                                        pipe = pipe,
                                                                                        CV_splitter = RepeatedStratifiedKFold_splitter,
                                                                                        class_labels=class_labels,
                                                                                        sample_IDs = sample_IDs,
                                                                                        scorers=scorers,
                                                                                        scoring_names=scoring_names,
                                                                                        num_null_iters=num_null_iters,
                                                                                        num_folds = num_folds,
                                                                                        num_repeats = num_repeats,
                                                                                        num_jobs = num_jobs)

                # Assign key details to dataframes
                main_classification_res["Disorder"] = disorder
                main_classification_res["Study"] = dataset_ID
                main_classification_res["Analysis_Type"] = Analysis_Type
                main_classification_res["group_var"] = grouping_var
                main_classification_res["Classifier_Type"] = classifier_type

                # Save to feather file
                main_classification_res.reset_index().to_feather(output_file)

            else:
                main_classification_res = pd.read_feather(output_file)

            # Append to list
            auc_results_list.append(main_classification_res)


# Combine all results
all_auc_results = pd.concat(auc_results_list).reset_index(drop=True)

# Save results
all_auc_results.to_feather(f"{data_path}/classification_results/All_AUC_results.feather")