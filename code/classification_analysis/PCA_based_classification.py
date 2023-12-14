import sys
sys.path.append("/headnode1/abry4213/github/fMRI_FeaturesDisorders/helper_functions/classification/")
from mixed_sigmoid_normalisation import MixedSigmoidScaler

import pandas as pd
from sklearn import svm
from sklearn import metrics
from sklearn.pipeline import Pipeline
import os.path
from sklearn.model_selection import StratifiedKFold, cross_val_predict, cross_validate, permutation_test_score
import numpy as np


data_path_base = "/Users/abry4213/data/"
univariate_feature_set = "catch25"
pairwise_feature_set = "pyspi14"
num_repeats = 10
num_folds = 10

study_group_df = pd.DataFrame({"Study": ["UCLA_CNP", "UCLA_CNP", "UCLA_CNP", "ABIDE_ASD"],
                             "Noise_Proc": ["AROMA+2P+GMR", "AROMA+2P+GMR", "AROMA+2P+GMR", "FC1000"],
                             "Comparison_Group": ["Schizophrenia", "Bipolar", "ADHD", "ASD"],
                             "Group_Nickname": ["SCZ", "BPD", "ADHD", "ASD"]})

# Load PCA data
PCA_scores_df = pd.read_feather(f"{data_path_base}/TS_feature_manuscript/Univariate_combo_first25_PCs.feather")

first25_PCs_results_list = []

for i in range(0, len(study_group_df.index)):
    dataset_ID = study_group_df.Study[i]
    noise_proc = study_group_df.Noise_Proc[i]
    comparison_group = study_group_df.Comparison_Group[i]
    group_nickname = study_group_df.Group_Nickname[i]

    noise_label = noise_proc.replace("+", "_")
    data_path = f"{data_path_base}/{dataset_ID}/"

    # Merge brain region + univariate feature name
    feature_data = PCA_scores_df.query("Diagnosis in ['Control', @comparison_group] & Comparison_Group == @comparison_group & Study == @dataset_ID") 
    
    # Extract only the combo feature data
    features_only = feature_data.drop(["Sample_ID", "Diagnosis", "Comparison_Group", "Study"],axis=1).to_numpy()
    
    # Extract sample ID and diagnosis as lists
    class_labels = [int(i==comparison_group) for i in feature_data["Diagnosis"].tolist()]
    sample_IDs = feature_data["Sample_ID"].tolist()

    # Define pipeline
    pipe = Pipeline([('scaler', MixedSigmoidScaler(unit_variance=True)), 
                     ('SVM', svm.SVC(kernel = "linear", C = 1, class_weight = "balanced"))])

    
    # Run 10 repeats of 10-fold CV SVM to measure balanced accuracy
    test_metrics_by_fold_list = []
        
    for r in range(num_repeats):
        # Initialise sample and class dataframe for repeat
        skf = StratifiedKFold(n_splits=num_folds, shuffle=True, random_state=r)

        # Fit SVM with 10-fold cross-validation for balanced accuracy
        cv_results_balacc = cross_validate(pipe, features_only, class_labels, 
                                    cv=skf, scoring=["balanced_accuracy"],
                                    return_estimator=False)
        
        # Extract balanced accuracy and ROC AUC by fold
        test_metrics_by_fold_df = pd.DataFrame({"Analysis_Type": "Univariate_Combo_25_PCs",
                                                "Study": dataset_ID,
                                                "Comparison_Group": comparison_group,
                                        "Fold": [*range(1, num_folds + 1, 1)],
                                        "Repeat_Number": r,
                                        "Balanced_Accuracy": cv_results_balacc["test_balanced_accuracy"]})
        test_metrics_by_fold_list.append(test_metrics_by_fold_df)

        
    test_metrics_by_fold = pd.concat(test_metrics_by_fold_list)
    first25_PCs_results_list.append(test_metrics_by_fold)

# Combine regularized SVM results across groups
first25_PCs_results = pd.concat(first25_PCs_results_list).reset_index()

# Write results to feather files
first25_PCs_results.to_feather(f"{data_path_base}/TS_feature_manuscript/Univariate_combo_first25_PCs_balanced_accuracy.feather")