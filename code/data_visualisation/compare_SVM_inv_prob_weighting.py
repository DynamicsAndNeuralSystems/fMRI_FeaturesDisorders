import pandas as pd
import scipy.io 
import sys
import os
import pyarrow.feather
from sklearn import svm
from sklearn.pipeline import Pipeline
import os.path
from sklearn.model_selection import StratifiedKFold, cross_val_predict, cross_validate, permutation_test_score
import numpy as np
import random
from sklearn.decomposition import PCA

comparison_to_control_group = "ADHD"
analysis_type = "Univariate_TS_Feature"
classifier_Type = "Linear_SVM"
num_folds = 10
num_repeats = 10

# Load metadata
UCLA_CNP_metadata = pd.read_feather("/Users/abry4213/data/UCLA_CNP/study_metadata/UCLA_CNP_sample_metadata.feather")
ABIDE_metadata = pd.read_feather("/Users/abry4213/data/ABIDE_ASD/study_metadata/ABIDE_ASD_sample_metadata.feather")

# Load in data containing subjects with both univariate and pairwise data available
UCLA_CNP_samples_to_keep = pd.read_feather("/Users/abry4213/data/UCLA_CNP/processed_data/UCLA_CNP_filtered_sample_info_AROMA_2P_GMR_final.feather") 
ABIDE_samples_to_keep = pd.read_feather("/Users/abry4213/data/ABIDE_ASD/processed_data/ABIDE_ASD_filtered_sample_info_FC1000_final.feather")       
samples_to_keep = pd.concat([UCLA_CNP_samples_to_keep, ABIDE_samples_to_keep], axis=0).Sample_ID.tolist()                                                  

# Load in catch25 data
UCLA_CNP_univariate_feature_data = pd.read_feather("/Users/abry4213/data/UCLA_CNP/processed_data/UCLA_CNP_AROMA_2P_GMR_catch25_filtered.feather").merge(UCLA_CNP_metadata, on='Sample_ID', how='left').drop(["Age", "Sex"],
                                                                        axis = 1).assign(Study = "UCLA_CNP")
ABIDE_univariate_feature_data = pd.read_feather("/Users/abry4213/data/ABIDE_ASD/processed_data/ABIDE_ASD_FC1000_catch25_filtered.feather").merge(ABIDE_metadata, on='Sample_ID', how='left').drop(["Age", "Sex"],
                                                                        axis = 1).assign(Study = "ABIDE_ASD")


# Filter univariate data by samples with both univariate and pairwise
# Filter by samples with univariate data available as well
univariate_feature_data = pd.concat([UCLA_CNP_univariate_feature_data, ABIDE_univariate_feature_data], axis=0)
univariate_feature_data = univariate_feature_data[univariate_feature_data.Sample_ID.isin(samples_to_keep)]                                                                           


# Define helper function to compare linear SVM with versus without inverse probability weighting
def compare_SVM_weightings(univariate_feature_data, study, comparison_to_control_group, grouping_var_name, grouping_var_type, num_repeats=10, num_folds=10):
    
    # Subset data to TS feature
    if grouping_var_type == "TS_Feature":
        TS_feature_data = univariate_feature_data.query("names == @grouping_var_name & Study == @study & Diagnosis in ['Control', @comparison_to_control_group]").drop(["names", "Noise_Proc",
                                                                                "method"], axis=1)
        # Pivot from long to wide
        feature_data_wide = TS_feature_data.pivot(index=['Sample_ID', 'Diagnosis'], columns='Brain_Region', values='values')
        
    else:
        TS_feature_data = univariate_feature_data.query("Brain_Region == @grouping_var_name & Study == @study & Diagnosis in ['Control', @comparison_to_control_group]").drop(["Brain_Region", "Noise_Proc",
                                                                                "method"], axis=1)
        # Pivot from long to wide
        feature_data_wide = TS_feature_data.pivot(index=['Sample_ID', 'Diagnosis'], columns='names', values='values')
        
    # Check that there are two groups for Diagnosis
    if len(TS_feature_data.Diagnosis.unique()) < 2:
        return

    # Extract sample ID and diagnosis as lists
    index_data = feature_data_wide.index.to_frame().reset_index(drop=True)
    class_labels = [int(i==comparison_to_control_group) for i in index_data["Diagnosis"].tolist()]
    sample_IDs = index_data["Sample_ID"].tolist()

    # Extract only the feature data
    features_only = feature_data_wide.reset_index(drop=True).to_numpy()


    # Define pipeline WITH vs WITHOUT balanced accuracy
    pipe = Pipeline([('scaler', MixedSigmoidScaler(unit_variance=True)), 
                     ('SVM', svm.SVC(kernel = "linear", C = 1))])
    
    pipe_bal = Pipeline([('scaler', MixedSigmoidScaler(unit_variance=True)), 
                     ('SVM', svm.SVC(kernel = "linear", C = 1, class_weight = "balanced"))])
    
    test_metrics_by_fold_list = []
    
    for i in range(num_repeats):
        # Initialise sample and class dataframe for repeat
        skf = StratifiedKFold(n_splits=num_folds, shuffle=True, random_state=i)
        
        # Find splits
        splits = list(skf.split(features_only, class_labels))
        fold_split_list = []
        
        # Fit SVM with 10-fold cross-validation for balanced accuracy
        cv_results_balacc = cross_validate(pipe, features_only, class_labels, 
                                    cv=skf, scoring=["accuracy", "balanced_accuracy"],
                                    n_jobs = 1,
                                    return_estimator=True)
        
        cv_results_balacc_invprob = cross_validate(pipe_bal, features_only, class_labels, 
                                    cv=skf, scoring=["accuracy", "balanced_accuracy"],
                                    n_jobs = 1,
                                    return_estimator=True)
        
        
        # Extract balanced accuracy by fold
        
        # One for no reweighitng
        test_metrics_by_fold_df = pd.DataFrame({"group_var": grouping_var_name,
                                                "Analysis_Type": grouping_var_type,
                                                "Weighting_Type": "None",
                                                "Fold": [*range(1, num_folds + 1, 1)],
                                                "Repeat_Number": i,
                                                "Balanced_Accuracy": cv_results_balacc["test_balanced_accuracy"],
                                                "Accuracy": cv_results_balacc["test_accuracy"]})
        
        # One with inverse probability reweiging
        test_metrics_by_fold_df_invprob = pd.DataFrame({"group_var": grouping_var_name,
                                                "Analysis_Type": grouping_var_type,
                                                "Weighting_Type": "Balanced",
                                                "Fold": [*range(1, num_folds + 1, 1)],
                                                "Repeat_Number": i,
                                                "Balanced_Accuracy": cv_results_balacc_invprob["test_balanced_accuracy"],
                                                "Accuracy": cv_results_balacc_invprob["test_accuracy"]})
        
        
        test_metrics_by_fold_list.append(test_metrics_by_fold_df)
        test_metrics_by_fold_list.append(test_metrics_by_fold_df_invprob)
    
    
    test_metrics_by_fold = pd.concat(test_metrics_by_fold_list)
    test_metrics_by_fold["group_var"] = grouping_var_name
    
    return test_metrics_by_fold


# Iterate over each brain region and time-series feature across all three UCLA CNP disorders
study_group_df = pd.DataFrame({"study": ["UCLA_CNP","UCLA_CNP","UCLA_CNP","ABIDE_ASD"],
                               "disorder": ["Schizophrenia", "Bipolar", "ADHD", "ASD"]})

balacc_res_list = []

# UCLA CNP disorders
for index, row in study_group_df.iterrows():
    study = row["study"]
    disorder = row["disorder"]
    
    for brain_region in univariate_feature_data.Brain_Region.unique().tolist():
        try:
            region_balacc = compare_SVM_weightings(univariate_feature_data = univariate_feature_data, 
                                               study = study,
                                                         comparison_to_control_group = disorder, 
                                                         grouping_var_type = "Brain_Region",
                                                         grouping_var_name = brain_region)
            region_balacc["Comparison_Group"] = disorder
            balacc_res_list.append(region_balacc)
        except:
            pass
        
    # Iterate over catch25 features
    for catch25_feature in univariate_feature_data.names.unique().tolist():
        try:
            feature_balacc = compare_SVM_weightings(univariate_feature_data = univariate_feature_data, 
                                                   study = study,
                                                             comparison_to_control_group = disorder, 
                                                             grouping_var_type = "TS_Feature",
                                                             grouping_var_name = catch25_feature)
            feature_balacc["Comparison_Group"] = disorder
            balacc_res_list.append(feature_balacc)
        except:
            pass
        

# Concatenate results and save to feather file
balacc_res = pd.concat(balacc_res_list).reset_index()
balacc_res.to_feather("/Users/abry4213/data/TS_feature_manuscript/SVM_with_vs_without_inv_prob_weighting.feather")


###############################################################################
# Also re-run the analysis and null dist for ADHD low_freq_power


def run_k_fold_SVM_for_feature_weighting(feature_data, 
                               grouping_var_name,
                               analysis_type,
                               sample_IDs,
                               class_labels,
                               inv_prob_weighting=True,
                               kernel="linear",
                               num_folds = 10,
                               num_repeats = 10,
                               num_jobs = 10):
        
    
    print(f"Creating {classifier_type} pipeline.")

    # Define the pipeline
    if inv_prob_weighting:
        if kernel=="linear":
            pipe = Pipeline([('scaler', MixedSigmoidScaler(unit_variance=True)), 
                             ('SVM', svm.SVC(kernel = "linear", C = 1, class_weight = "balanced"))])
        else:
            pipe = Pipeline([('scaler', MixedSigmoidScaler(unit_variance=True)), 
                             ('SVM', svm.SVC(kernel = "rbf", C = 1, class_weight = "balanced"))])
            
    else: 
        if kernel=="linear":
            pipe = Pipeline([('scaler', MixedSigmoidScaler(unit_variance=True)),
                             ('SVM', svm.SVC(kernel = "linear", C = 1))])
        else:
            pipe = Pipeline([('scaler', MixedSigmoidScaler(unit_variance=True)), 
                             ('SVM', svm.SVC(kernel = "rbf", C = 1))])

    # Define lists for: 
    # (1) balanced accuracy by fold/repeat
    test_metrics_by_fold_list = []
    
    for i in range(num_repeats):
        # Initialise sample and class dataframe for repeat
        skf = StratifiedKFold(n_splits=num_folds, shuffle=True, random_state=i)
        
        # Fit SVM with 10-fold cross-validation for balanced accuracy
        cv_results_balacc = cross_validate(pipe, feature_data, class_labels, 
                                    cv=skf, scoring=["balanced_accuracy"],
                                    n_jobs = int(num_jobs),
                                    return_estimator=True)
        

        # Extract balanced accuracy and ROC AUC by fold
        test_metrics_by_fold_df = pd.DataFrame({"group_var": grouping_var_name,
                                                "Analysis_Type": analysis_type,
                                        "Fold": [*range(1, num_folds + 1, 1)],
                                        "Repeat_Number": i,
                                        "Balanced_Accuracy": cv_results_balacc["test_balanced_accuracy"]})
        test_metrics_by_fold_list.append(test_metrics_by_fold_df)

    test_metrics_by_fold = pd.concat(test_metrics_by_fold_list).reset_index()
    
    
    test_metrics_by_fold["Kernel"] = kernel
    if inv_prob_weighting:
        test_metrics_by_fold["SVM_Weighting"] = "Balanced"
    else:
        test_metrics_by_fold["SVM_Weighting"] = "None"
    return (test_metrics_by_fold)


comparison_to_control_group = "ADHD"
classifier_type = "Linear_SVM"
num_jobs = 1
TS_feature = "SP_Summaries_welch_rect_area_5_1"
TS_feature_data = univariate_feature_data.query("names == @TS_feature & Diagnosis in ['Control', @comparison_to_control_group]").drop(["names", "Noise_Proc",
                                                                        "method"], axis=1)
# Pivot from long to wide
feature_data_wide = TS_feature_data.pivot(index=['Sample_ID', 'Diagnosis'], columns='Brain_Region', values='values')

# Extract sample ID and diagnosis as lists
index_data = feature_data_wide.index.to_frame().reset_index(drop=True)
class_labels = [int(i==comparison_to_control_group) for i in index_data["Diagnosis"].tolist()]
sample_IDs = index_data["Sample_ID"].tolist()

# Extract only the feature data
features_only = feature_data_wide.reset_index(drop=True).to_numpy()

# Main analysis with linear kernel
ADHD_low_freq_linear_balacc_by_fold_inv_prob = run_k_fold_SVM_for_feature_weighting(feature_data = features_only, 
    grouping_var_name = TS_feature,
    analysis_type = "Univariate_TS_Feature",
    sample_IDs = sample_IDs,
    class_labels = class_labels,
    inv_prob_weighting=True,
    kernel="linear",
    num_folds = num_folds,
    num_jobs = num_jobs,
    num_repeats = num_repeats)

ADHD_low_freq_linear_balacc_by_fold_no_reweighting = run_k_fold_SVM_for_feature_weighting(feature_data = features_only, 
    grouping_var_name = TS_feature,
    analysis_type = "Univariate_TS_Feature",
    sample_IDs = sample_IDs,
    class_labels = class_labels,
    inv_prob_weighting=False,
    kernel="linear",
    num_folds = num_folds,
    num_jobs = num_jobs,
    num_repeats = num_repeats)

# Main analysis with RBF kernel
ADHD_low_freq_rbf_balacc_by_fold_inv_prob = run_k_fold_SVM_for_feature_weighting(feature_data = features_only, 
    grouping_var_name = TS_feature,
    analysis_type = "Univariate_TS_Feature",
    sample_IDs = sample_IDs,
    class_labels = class_labels,
    inv_prob_weighting=True,
    kernel="rbf",
    num_folds = num_folds,
    num_jobs = num_jobs,
    num_repeats = num_repeats)
ADHD_low_freq_rbf_balacc_by_fold_no_reweighting = run_k_fold_SVM_for_feature_weighting(feature_data = features_only, 
    grouping_var_name = TS_feature,
    analysis_type = "Univariate_TS_Feature",
    sample_IDs = sample_IDs,
    class_labels = class_labels,
    inv_prob_weighting=False,
    kernel="rbf",
    num_folds = num_folds,
    num_jobs = num_jobs,
    num_repeats = num_repeats)


# Concatenate and save results
ADHD_main_balacc_res = pd.concat([ADHD_low_freq_linear_balacc_by_fold_inv_prob,
                                  ADHD_low_freq_linear_balacc_by_fold_no_reweighting,
                                  ADHD_low_freq_rbf_balacc_by_fold_inv_prob,
                                  ADHD_low_freq_rbf_balacc_by_fold_no_reweighting]).reset_index()
ADHD_main_balacc_res.to_feather("/Users/abry4213/data/TS_feature_manuscript/ADHD_low_freq_power_main_balacc_res.feather")

##############################################################################
# Null analyses
null_balacc_res_list = []
for null_iter in range(0, 1000):
    iter_class_labels = [int(i==comparison_to_control_group) for i in index_data["Diagnosis"].tolist()]
    iter_class_labels_shuffled = iter_class_labels.copy()
    random.shuffle(iter_class_labels_shuffled)
    
    # Inverse probability weighting
    ADHD_low_freq_power_null_inv_prob = run_k_fold_SVM_for_feature_weighting(feature_data = features_only, 
        grouping_var_name = TS_feature,
        analysis_type = "Univariate_TS_Feature",
        sample_IDs = sample_IDs,
        class_labels = iter_class_labels_shuffled,
        kernel = "linear",
        inv_prob_weighting=True,
        num_folds = num_folds,
        num_jobs = num_jobs,
        num_repeats = num_repeats)
    
    ADHD_low_freq_power_null_inv_prob["Null_Iteration"] = null_iter
    
    # Take the mean balanced accuracy for this iteration
    ADHD_low_freq_power_null_inv_prob_averaged = (ADHD_low_freq_power_null_inv_prob
                              .groupby("Repeat_Number")["Balanced_Accuracy"]
                              .mean()
                              .to_frame()
                              .reset_index())["Balanced_Accuracy"].mean()
    
    ADHD_low_freq_power_null_inv_prob_res = pd.DataFrame(data = [[comparison_to_control_group, TS_feature, ADHD_low_freq_power_null_inv_prob_averaged, null_iter, "Balanced"]],
                                                columns = ["Comparison_Group", "grouping_var", "Null_Balanced_Accuracy", "Null_Iteration", "SVM_Weighting"])
        
    
    null_balacc_res_list.append(ADHD_low_freq_power_null_inv_prob_res)
    
    # No reweighting
    ADHD_low_freq_power_null_no_reweighting = run_k_fold_SVM_for_feature_weighting(feature_data = features_only, 
        grouping_var_name = TS_feature,
        analysis_type = "Univariate_TS_Feature",
        sample_IDs = sample_IDs,
        class_labels = iter_class_labels_shuffled,
        kernel = "linear",
        inv_prob_weighting=False,
        num_folds = num_folds,
        num_jobs = num_jobs,
        num_repeats = num_repeats)
    
    ADHD_low_freq_power_null_no_reweighting["Null_Iteration"] = null_iter
    
    # Take the mean balanced accuracy for this iteration
    ADHD_low_freq_power_null_no_reweighting_averaged = (ADHD_low_freq_power_null_no_reweighting
                              .groupby("Repeat_Number")["Balanced_Accuracy"]
                              .mean()
                              .to_frame()
                              .reset_index())["Balanced_Accuracy"].mean()
    
    ADHD_low_freq_power_null_no_reweighting_res = pd.DataFrame(data = [[comparison_to_control_group, TS_feature, ADHD_low_freq_power_null_no_reweighting_averaged, null_iter, "None"]],
                                                columns = ["Comparison_Group", "grouping_var", "Null_Balanced_Accuracy", "Null_Iteration", "SVM_Weighting"])
        
        
    null_balacc_res_list.append(ADHD_low_freq_power_null_no_reweighting_res)
    
ADHD_low_freq_power_null_res = pd.concat(null_balacc_res_list).reset_index()
ADHD_low_freq_power_null_res.to_feather("/Users/abry4213/data/TS_feature_manuscript/ADHD_low_freq_power_nulls.feather")


###############################################################################
# Normalization for both train/test vs. train only


def run_k_fold_SVM_for_normalization_type(feature_data, 
                               grouping_var_name,
                               analysis_type,
                               sample_IDs,
                               class_labels,
                               inv_prob_weighting=True,
                               normalization_type="train",
                               num_folds = 10,
                               num_repeats = 10,
                               num_jobs = 10):
        
    
    print(f"Creating {classifier_type} pipeline.")

    # Define the pipeline
    if inv_prob_weighting:
        if normalization_type=="train":
            print("Running linear SVM with inverse probability weighting and training normalization")
            pipe = Pipeline([('scaler', MixedSigmoidScaler(unit_variance=True)), 
                             ('SVM', svm.SVC(kernel = "linear", C = 1, class_weight = "balanced"))])
            feature_data_to_use = feature_data
        else:
            print("Running linear SVM with inverse probability weighting and full normalization")
            feature_data_to_use = feature_data
            transformer = MixedSigmoidScaler(unit_variance=True).fit(feature_data_to_use)
            feature_data_to_use = transformer.transform(feature_data_to_use)
            pipe = Pipeline([('SVM', svm.SVC(kernel = "linear", C = 1, class_weight = "balanced"))])
            
    else: 
        if normalization_type=="train":
            print("Running linear SVM with training normalization")
            pipe = Pipeline([('scaler', MixedSigmoidScaler(unit_variance=True)),
                             ('SVM', svm.SVC(kernel = "linear", C = 1))])
            feature_data_to_use = feature_data
        else:
            print("Running linear SVM with full normalization")
            feature_data_to_use = feature_data
            transformer = MixedSigmoidScaler(unit_variance=True).fit(feature_data_to_use)
            feature_data_to_use = transformer.transform(feature_data_to_use)
            pipe = Pipeline([('SVM', svm.SVC(kernel = "linear", C = 1))])

    # Define lists for: 
    # (1) balanced accuracy by fold/repeat
    test_metrics_by_fold_list = []
    
    for i in range(num_repeats):
        # Initialise sample and class dataframe for repeat
        skf = StratifiedKFold(n_splits=num_folds, shuffle=True, random_state=i)
        
        # Fit SVM with 10-fold cross-validation for balanced accuracy
        cv_results_balacc = cross_validate(pipe, feature_data_to_use, class_labels, 
                                    cv=skf, scoring=["balanced_accuracy"],
                                    n_jobs = int(num_jobs),
                                    return_estimator=True)
        

        # Extract balanced accuracy and ROC AUC by fold
        test_metrics_by_fold_df = pd.DataFrame({"group_var": grouping_var_name,
                                                "Analysis_Type": analysis_type,
                                        "Fold": [*range(1, num_folds + 1, 1)],
                                        "Repeat_Number": i,
                                        "Balanced_Accuracy": cv_results_balacc["test_balanced_accuracy"]})
        test_metrics_by_fold_list.append(test_metrics_by_fold_df)

    test_metrics_by_fold = pd.concat(test_metrics_by_fold_list).reset_index()
    
    
    test_metrics_by_fold["Normalization"] = normalization_type
    if inv_prob_weighting:
        test_metrics_by_fold["SVM_Weighting"] = "Balanced"
    else:
        test_metrics_by_fold["SVM_Weighting"] = "None"
    return (test_metrics_by_fold, feature_data_to_use)

comparison_to_control_group = "ADHD"
classifier_type = "Linear_SVM"
num_jobs = 1
TS_feature = "SP_Summaries_welch_rect_area_5_1"
TS_feature_data = univariate_feature_data.query("names == @TS_feature & Diagnosis in ['Control', @comparison_to_control_group]").drop(["names", "Noise_Proc",
                                                                        "method"], axis=1)
# Pivot from long to wide
feature_data_wide = TS_feature_data.pivot(index=['Sample_ID', 'Diagnosis'], columns='Brain_Region', values='values')

# Extract sample ID and diagnosis as lists
index_data = feature_data_wide.index.to_frame().reset_index(drop=True)
class_labels = [int(i==comparison_to_control_group) for i in index_data["Diagnosis"].tolist()]
sample_IDs = index_data["Sample_ID"].tolist()

# Extract only the feature data
features_only = feature_data_wide.reset_index(drop=True).to_numpy()

(ADHD_low_freq_power_train_normalized, ADsaHD_low_freq_power_train_normalized_input) = run_k_fold_SVM_for_normalization_type(feature_data = features_only, 
    grouping_var_name = TS_feature,
    analysis_type = "Univariate_TS_Feature",
    sample_IDs = sample_IDs,
    class_labels = class_labels,
    normalization_type="train",
    inv_prob_weighting=True,
    num_folds = num_folds,
    num_jobs = num_jobs,
    num_repeats = num_repeats)

(ADHD_low_freq_power_full_normalized, ADHD_low_freq_power_full_normalized_input) = run_k_fold_SVM_for_normalization_type(feature_data = features_only, 
    grouping_var_name = TS_feature,
    analysis_type = "Univariate_TS_Feature",
    sample_IDs = sample_IDs,
    class_labels = class_labels,
    normalization_type="full",
    inv_prob_weighting=True,
    num_folds = num_folds,
    num_jobs = num_jobs,
    num_repeats = num_repeats)

ADHD_low_freq_power_normalization_comparison_res = pd.concat([ADHD_low_freq_power_train_normalized, ADHD_low_freq_power_full_normalized]).reset_index()
ADHD_low_freq_power_normalization_comparison_res.to_feather("/Users/abry4213/data/TS_feature_manuscript/ADHD_low_freq_power_normalization.feather")


###############################################################################
# First two PCs as classifier input

# Using the same features_only dataset as above

# Apply MixedSigmoid scaling
features_only_scaled = MixedSigmoidScaler().fit_transform(features_only)

# Run PCA
ADHD_low_freq_power_pca = PCA(n_components=2)
ADHD_low_freq_power_pca_scores = ADHD_low_freq_power_pca.fit_transform(features_only_scaled)
# Save scores to feather
ADHD_low_freq_power_pca_scores_df = pd.DataFrame(ADHD_low_freq_power_pca_scores, columns=["PC1", "PC2"])
ADHD_low_freq_power_pca_scores_df = pd.concat([ADHD_low_freq_power_pca_scores_df, index_data], axis=1).reset_index()
ADHD_low_freq_power_pca_scores_df.to_feather("/Users/abry4213/data/TS_feature_manuscript/ADHD_low_freq_power_pca_scores.feather")


# Run classifier based on these first two components

ADHD_low_freq_power_pca_inv_prob_res = run_k_fold_SVM_for_feature_weighting(feature_data=ADHD_low_freq_power_pca_scores, 
                               grouping_var_name="PC_Scores",
                               analysis_type="Univariate_TS_Feature",
                               sample_IDs = sample_IDs,
                               class_labels = class_labels,
                               inv_prob_weighting=True,
                               kernel="linear",
                               num_folds = 10,
                               num_repeats = 10,
                               num_jobs = 10)

ADHD_low_freq_power_pca_no_reweighting_res = run_k_fold_SVM_for_feature_weighting(feature_data=ADHD_low_freq_power_pca_scores, 
                               grouping_var_name="PC_Scores",
                               analysis_type="Univariate_TS_Feature",
                               sample_IDs = sample_IDs,
                               class_labels = class_labels,
                               inv_prob_weighting=False,
                               kernel="linear",
                               num_folds = 10,
                               num_repeats = 10,
                               num_jobs = 10)

# Concatenate and save to feather file
ADHD_low_freq_power_pca_classification_res = pd.concat([ADHD_low_freq_power_pca_inv_prob_res, ADHD_low_freq_power_pca_no_reweighting_res]).reset_index()
ADHD_low_freq_power_pca_classification_res.to_feather("/Users/abry4213/data/TS_feature_manuscript/ADHD_low_freq_power_pca_classification_res.feather")
