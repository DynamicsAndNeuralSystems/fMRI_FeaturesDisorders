import os.path
import dill
import copy
from dplython import (DplyFrame, 
                      X, 
                      diamonds, 
                      dfilter,
                      select,
                      rename,
                      sift, 
                      sample_n,
                      sample_frac, 
                      head, 
                      arrange,
                      mutate,
                      nrow,
                      group_by,
                      summarize, 
                      DelayFunction) 
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import svm
from sklearn.model_selection import StratifiedKFold, cross_val_score, cross_val_predict, cross_validate

# Define paths specific to the UCLA CNP and ABIDE ASD datasets
univariate_feature_set = "catch22"
pairwise_feature_set = "pyspi14"
github_dir = "/Users/abry4213/github/fMRI_FeaturesDisorders/"
dataset_ID = "UCLA_CNP"
data_path = "/Users/abry4213/data/UCLA_CNP/"
rdata_path = data_path + "processed_data/Rdata/"
pydata_path = data_path + "processed_data/pydata/"
noise_proc = "AROMA+2P+GMR"
metadata_file = "UCLA_CNP_sample_metadata.pkl"

# Load metadata
with open(data_path + "study_metadata/" + metadata_file, "rb") as f:
    metadata = dill.load(f)

# Univariate feature data
univariate_feature_file = rdata_path + "UCLA_CNP_AROMA_2P_GMR_catch22_filtered_zscored.feather"
univariate_feature_data = pd.read_feather(univariate_feature_file).merge(metadata, on='Sample_ID', how='left').drop(["Age", "Sex", "Study"],
                                                                           axis = 1)

# Run SVM per ROI
fold_assignments_per_ROI_list = []
SVM_coefficients_per_ROI_list = []
balanced_accuracy_per_ROI_list = []
CV_sample_predictions_per_ROI_list = []

# Split it up by first brain reigon
for ROI in univariate_feature_data.Brain_Region.unique().tolist():
    
    # Subset data to ROI
    region_data = univariate_feature_data.query("Brain_Region == @ROI & Diagnosis in ['Control', 'Schizophrenia']").drop(["Brain_Region", "Noise_Proc",
                                                                              "method"], axis=1)
    
    # Pivot from long to wide
    region_data_wide = region_data.pivot(index=['Sample_ID', 'Diagnosis'], columns='names', values='values')
    
    # Extract name of features
    feature_list = region_data_wide.columns.tolist()
    
    # Extract sample ID and diagnosis as lists
    index_data = region_data_wide.index.to_frame().reset_index(drop=True)
    sample_ID = index_data["Sample_ID"].tolist()
    class_labels = index_data["Diagnosis"].tolist()
    
    # Extract only the feature data
    features_only = region_data_wide.reset_index(drop=True).to_numpy()
    
    # Fit the SVM
    SVM_model = svm.SVC(kernel="linear", C=1, shrinking=False, class_weight="balanced")
    
    # Define lists for: 
    # (1) fold assignments by repeat,
    # (2) SVM coefficients,
    # (3) balanced accuracy by fold/repeat,
    # (4) individual sample predictions per repeat
    fold_assignments_list = []
    SVM_coefficients_list = []
    balanced_accuracy_list = []
    CV_sample_predictions_list = []
    
    # Get 10-repeat 10-fold balanced accuracy
    for i in range(10):
        skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=i)
        repeat_labels = copy.deepcopy(index_data)
        
        # Find splits
        splits = list(skf.split(features_only, class_labels))
        fold_split_list = []
        
        # Fit SVM with 10-fold cross-validation
        cv_results = cross_validate(SVM_model, features_only, class_labels, 
                                    cv=skf, scoring="balanced_accuracy",
                                    return_estimator=True)
    
        
        # Extract fold number and feature coefficients
        coef_list = [pd.DataFrame(svm.coef_.T) for svm in cv_results['estimator']]
        for f in range(10):
            # Split for fold number
            test_indices = pd.DataFrame(splits[f][1], columns=["Sample_Index"])
            test_indices["Fold"] = f+1
            test_indices["Repeat"] = i+1
            fold_split_list.append(test_indices)
            
            # Coefficients
            coef_list[f].reset_index(inplace=True)
            coef_list[f] = coef_list[f].rename(columns = {'index':'Feature_Number',
                                                          0: "Coefficient"})
            coef_list[f]["Feature Name"] = feature_list
            coef_list[f]["Fold"] = f+1
            coef_list[f]["Repeat Number"] = i+1
        # Combine lists into dataframes
        fold_splits = pd.concat(fold_split_list)
        fold_assignments_list.append(fold_splits)
        coef_df = pd.concat(coef_list)
        SVM_coefficients_list.append(coef_df)
        
        # Extract balanced accuracy by fold
        balanced_accuracy_by_fold_df = pd.DataFrame(cv_results["test_score"],
                                                    columns=["Balanced_Accuracy"])
        balanced_accuracy_by_fold_df["Fold"] = [*range(1, 11, 1)]
        balanced_accuracy_by_fold_df["Repeat_Number"] = i
        balanced_accuracy_list.append(balanced_accuracy_by_fold_df)
        
        # Generate CV predictions across folds and save
        CV_pred = cross_val_predict(SVM_model, features_only, class_labels, cv=skf)
        repeat_labels["CV_Predicted_Diagnosis"] = CV_pred
        repeat_labels["Repeat_Number"] = i+1
        CV_sample_predictions_list.append(repeat_labels)
        
    # Concatenate results and save per ROI
    fold_assignments = pd.concat(fold_assignments_list)
    fold_assignments["Brain_Region"] = ROI
    fold_assignments_per_ROI_list.append(fold_assignments)
    
    SVM_coefficients = pd.concat(SVM_coefficients_list)
    SVM_coefficients["Brain_Region"] = ROI
    SVM_coefficients_per_ROI_list.append(SVM_coefficients)
    
    balanced_accuracy = pd.concat(balanced_accuracy_list)
    balanced_accuracy["Brain_Region"] = ROI
    balanced_accuracy_per_ROI_list.append(balanced_accuracy)
    
    CV_sample_predictions = pd.concat(CV_sample_predictions_list)
    CV_sample_predictions["Brain_Region"] = ROI
    CV_sample_predictions_per_ROI_list.append(CV_sample_predictions)
    
# Save results
fold_assignments_per_ROI = pd.concat(fold_assignments_per_ROI_list).reset_index()
fold_assignments_per_ROI.to_feather(pydata_path + "UCLA_CNP_Univariate_ROI_wise_SVM_fold_assigments.feather")

SVM_coefficients_per_ROI = pd.concat(SVM_coefficients_per_ROI_list).reset_index()
SVM_coefficients_per_ROI.to_feather(pydata_path + "UCLA_CNP_Univariate_ROI_wise_SVM_fold_SVM_coefficients.feather")

balanced_accuracy_per_ROI = pd.concat(balanced_accuracy_per_ROI_list).reset_index()
balanced_accuracy_per_ROI.to_feather(pydata_path + "UCLA_CNP_Univariate_ROI_wise_SVM_balanced_accuracy.feather")

CV_sample_predictions_per_ROI = pd.concat(CV_sample_predictions_per_ROI_list).reset_index()
CV_sample_predictions_per_ROI.to_feather(pydata_path + "UCLA_CNP_Univariate_ROI_wise_SVM_sample_predictions.feather")

    


###############################################################################
# Pairwise feature data
pairwise_feature_file = rdata_path + "UCLA_CNP_AROMA_2P_GMR_pyspi14_filtered_zscored.feather"
pairwise_feature_data = pd.read_feather(pairwise_feature_file)