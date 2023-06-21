import pandas as pd
from sklearn import svm
from sklearn import metrics
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, RobustScaler
import os.path
from sklearn.model_selection import StratifiedKFold, cross_val_predict, cross_validate, permutation_test_score
import numpy as np

# dataset_ID = "UCLA_CNP"
# data_path = "/headnode1/abry4213/data/UCLA_CNP/"
# metadata_file = "UCLA_CNP_sample_metadata.feather"
# comparison_group = "Schizophrenia"
# univariate_feature_set = "catch24"
# pairwise_feature_set = "pyspi14"
# univariate_feature_file = f"{data_path}/processed_data/UCLA_CNP_AROMA_2P_GMR_catch24_filtered.feather"
# noise_proc = "AROMA+2P+GMR"
# num_folds = 10
# num_null_iters = 2
# num_repeats = 2
# num_jobs = 1
# scaling_type = "mixedsigmoid"

data_path_base = "/headnode1/abry4213/data/"
univariate_feature_set = "catch24"
pairwise_feature_set = "pyspi14"
num_repeats = 10
num_folds = 10

study_group_df = pd.DataFrame({"Study": ["UCLA_CNP", "UCLA_CNP", "UCLA_CNP", "ABIDE_ASD"],
                             "Noise_Proc": ["AROMA+2P+GMR", "AROMA+2P+GMR", "AROMA+2P+GMR", "FC1000"],
                             "Comparison_Group": ["Schizophrenia", "Bipolar", "ADHD", "ASD"],
                             "Group_Nickname": ["SCZ", "BPD", "ADHD", "ASD"]})


regularized_SVM_results_list = []
regularized_SVM_coefs_list = []

for i in range(0, len(study_group_df.index)):
    dataset_ID = study_group_df.Study[i]
    noise_proc = study_group_df.Noise_Proc[i]
    comparison_group = study_group_df.Comparison_Group[i]
    group_nickname = study_group_df.Group_Nickname[i]

    noise_label = noise_proc.replace("+", "_")
    data_path = f"{data_path_base}/{dataset_ID}/"
    metadata_file = f"{dataset_ID}_sample_metadata.feather"
    univariate_feature_file = f"{data_path}/processed_data/{dataset_ID}_{noise_label}_{univariate_feature_set}_filtered.feather"

    # Load metadata
    metadata = pd.read_feather(data_path + "study_metadata/" + metadata_file)

    # Load in data containing subjects with both univariate and pairwise data available
    samples_to_keep = pd.read_feather(f"{data_path}/processed_data/{dataset_ID}_filtered_sample_info_{noise_label}_{univariate_feature_set}_{pairwise_feature_set}.feather")                                                                           

    # Univariate feature data
    univariate_feature_data = pd.read_feather(univariate_feature_file).merge(metadata, on='Sample_ID', how='left').drop(["Age", "Sex"],
                                                                            axis = 1)

    # Filter univariate data by samples with both univariate and pairwise
    # Filter by samples with univariate data available as well
    univariate_feature_data = univariate_feature_data[univariate_feature_data.Sample_ID.isin(samples_to_keep.Sample_ID)]   


    # Merge brain region + univariate feature name
    combo_data = (univariate_feature_data
                .query("Diagnosis in ['Control', @comparison_group]")
                .drop(["method", "Noise_Proc"], axis=1))
    combo_data["Combo_Feature"] = combo_data.Brain_Region + "_" + combo_data.names
    combo_data = combo_data.drop(["Brain_Region", "names"], axis=1)
    
    # Pivot from long to wide
    combo_data_wide = combo_data.pivot(index=["Sample_ID", "Diagnosis"],
                                    columns = "Combo_Feature",
                                    values = "values")
    
    # Extract name of combo features
    combo_features = combo_data_wide.columns.tolist()
    
    # Extract only the combo feature data
    features_only = combo_data_wide.reset_index(drop=True).to_numpy()
    
    # Extract sample ID and diagnosis as lists
    index_data = combo_data_wide.index.to_frame().reset_index(drop=True)
    class_labels = [int(i==comparison_group) for i in index_data["Diagnosis"].tolist()]
    sample_IDs = index_data["Sample_ID"].tolist()


    # Define pipeline
    pipe = Pipeline([('scaler', MixedSigmoidScaler(unit_variance=True)), 
                     ('SVM', svm.LinearSVC(penalty="l1", C = 1, dual=False, class_weight = "balanced"))])
    
    
    # Run 10 repeats of 10-fold CV SVM to measure balanced accuracy
    test_metrics_by_fold_list = []
    SVM_coefficients_list = []
        
    for r in range(num_repeats):
        # Initialise sample and class dataframe for repeat
        skf = StratifiedKFold(n_splits=num_folds, shuffle=True, random_state=r)

        # Fit SVM with 10-fold cross-validation for balanced accuracy
        cv_results_balacc = cross_validate(pipe, features_only, class_labels, 
                                    cv=skf, scoring=["balanced_accuracy"],
                                    return_estimator=True)
        
        # Extract balanced accuracy and ROC AUC by fold
        test_metrics_by_fold_df = pd.DataFrame({"Analysis_Type": "Univariate_Combo_Regularized",
                                                "Study": dataset_ID,
                                                "Comparison_Group": comparison_group,
                                        "Fold": [*range(1, num_folds + 1, 1)],
                                        "Repeat_Number": r,
                                        "Balanced_Accuracy": cv_results_balacc["test_balanced_accuracy"]})
        test_metrics_by_fold_list.append(test_metrics_by_fold_df)
        
        # Extract SVM coefficients
        coef_list = []
        for f in range(num_folds):
            # Coefficients
            SVM_for_fold = cv_results_balacc["estimator"][f]["SVM"]
            SVM_coefs = pd.DataFrame(SVM_for_fold.coef_.T)
            SVM_coefs.reset_index(inplace=True)
            SVM_coefs = SVM_coefs.rename(columns = {'index':'Feature_Number',
                                                          0: "Coefficient"})
            SVM_coefs["Feature Name"] = combo_features
            SVM_coefs["Fold"] = f+1
            SVM_coefs["Repeat Number"] = r+1
            coef_list.append(SVM_coefs)
        
        # Add sample ID name based on index
        coef_df = pd.concat(coef_list)
        coef_df["Analysis_Type"] = "Univariate_Combo_Regularized"
        coef_df["Study"] = dataset_ID
        coef_df["Comparison_Group"] = comparison_group
        SVM_coefficients_list.append(coef_df)

        
    test_metrics_by_fold = pd.concat(test_metrics_by_fold_list)
    regularized_SVM_results_list.append(test_metrics_by_fold)
    
    SVM_coefficients = pd.concat(SVM_coefficients_list)
    regularized_SVM_coefs_list.append(SVM_coefficients)

# Combine regularized SVM results across groups
regularized_SVM_results = pd.concat(regularized_SVM_results_list).reset_index()
regularized_SVM_coefficients = pd.concat(regularized_SVM_coefs_list).reset_index()

# Write results to feather files
regularized_SVM_results.to_feather(f"{data_path_base}/TS_feature_manuscript/SVM_L1_Regularized_Balanced_Accuracy.feather")
regularized_SVM_coefficients.to_feather(f"{data_path_base}/TS_feature_manuscript/SVM_L1_Regularized_Coefficients.feather")