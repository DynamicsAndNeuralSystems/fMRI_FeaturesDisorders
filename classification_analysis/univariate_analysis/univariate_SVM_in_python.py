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
from sklearn.model_selection import StratifiedKFold, cross_val_score, cross_val_predict

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

# Merge with metadata

ROI = 'ctx-lh-bankssts'


region_SVM_list = {}
CV_predictions = {}
balanced_accuracy_res = {}

# Split it up by first brain reigon
for ROI in univariate_feature_data.Brain_Region.unique().tolist():
    
    # Subset data to ROI
    region_data = univariate_feature_data.query("Brain_Region == @ROI & Diagnosis in ['Control', 'Schizophrenia']").drop(["Brain_Region", "Noise_Proc",
                                                                              "method"], axis=1)
    
    # Pivot from long to wide
    region_data_wide = region_data.pivot(index=['Sample_ID', 'Diagnosis'], columns='names', values='values')
    
    # Extract sample ID and diagnosis as lists
    index_data = region_data_wide.index.to_frame().reset_index(drop=True)
    sample_ID = index_data["Sample_ID"].tolist()
    class_labels = index_data["Diagnosis"].tolist()
    
    # Extract only the feature data
    features_only = region_data_wide.reset_index(drop=True).to_numpy()
    
    # Fit the SVM
    SVM_model = svm.SVC(kernel="linear", C=1, shrinking=False, class_weight="balanced")
    
    # Save SVM to list
    region_SVM_list[ROI] = SVM_model
    
    CV_repeat_preds = []
    balanced_accuracy_list = []
    
    # Get 10-repeat 10-fold balanced accuracy
    for i in range(1,11):
        skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=i+10)
        repeat_labels = copy.deepcopy(index_data)
        
        # Find balanced accuracy
        balanced_accuracy_by_fold = cross_val_score(SVM_model, features_only, 
                                                    class_labels, cv=skf, 
                                                    scoring="balanced_accuracy")
        balanced_accuracy_by_fold_df = pd.DataFrame(balanced_accuracy_by_fold, columns=["Balanced_Accuracy"])
        balanced_accuracy_by_fold_df["Fold"] = [*range(1, 11, 1)]
        balanced_accuracy_by_fold_df["Repeat_Number"] = i
        balanced_accuracy_list.append(balanced_accuracy_by_fold_df)
        
        # Generate CV predictions across folds and save
        CV_pred = cross_val_predict(SVM_model, features_only, class_labels, cv=skf)
        repeat_labels["CV_Predicted_Diagnosis"] = CV_pred
        repeat_labels["Repeat_Number"] = i
        # add fold number column?
        CV_repeat_preds.append(repeat_labels)

    # Combine predictions per fold/repeat
    CV_predictions_region = pd.concat(CV_repeat_preds)
    CV_predictions[ROI] = CV_predictions_region
    # Combine balanced accuracy per fold/repeat
    balanced_accuracy_region = pd.concat(balanced_accuracy_list)
    balanced_accuracy_res[ROI] = balanced_accuracy_region

###############################################################################
# Pairwise feature data
pairwise_feature_file = rdata_path + "UCLA_CNP_AROMA_2P_GMR_pyspi14_filtered_zscored.feather"
pairwise_feature_data = pd.read_feather(pairwise_feature_file)