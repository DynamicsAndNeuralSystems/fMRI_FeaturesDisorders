import argparse
import pandas as pd
import numpy as np
import random
import pyarrow.feather as feather
from sklearn import svm
from sklearn.model_selection import StratifiedKFold, cross_val_predict, cross_validate, cross_val_score

# Command-line arguments to parse
parser = argparse.ArgumentParser(description='Process inputs for pairwise data preparation.')
parser.add_argument('--data_path', default="/headnode1/abry4213/data/UCLA_CNP/", dest='data_path')
parser.add_argument('--metadata_file', default="UCLA_CNP_sample_metadata.feather", dest='metadata_file')
parser.add_argument('--SPI_directionality_file', default="/headnode1/abry4213/github/fMRI_FeaturesDisorders/classification_analysis/pairwise_analysis/SPI_Direction_Info.csv", dest='SPI_directionality_file')
parser.add_argument('--comparison_group', default="Schizophrenia", dest='comparison_group')
parser.add_argument('--univariate_feature_set', default='catch22', dest='univariate_feature_set')
parser.add_argument('--pairwise_feature_set', default='pyspi14', dest='pairwise_feature_set')
parser.add_argument('--pairwise_feature_file', default="/headnode1/abry4213/data/UCLA_CNP/processed_data/UCLA_CNP_AROMA_2P_GMR_pyspi14_filtered_zscored.feather", dest='pairwise_feature_file')
parser.add_argument('--noise_proc', dest='noise_proc')
parser.add_argument('--num_null_iters', default=1000, dest='num_null_iters')
parser.add_argument('--dataset_ID', default="UCLA_CNP", dest='dataset_ID')

# Parse arguments
args = parser.parse_args()
data_path = args.data_path
metadata_file = args.metadata_file
SPI_directionality_file = args.SPI_directionality_file
comparison_group = args.comparison_group
univariate_feature_set = args.univariate_feature_set
pairwise_feature_set = args.pairwise_feature_set
pairwise_feature_file = args.pairwise_feature_file
noise_proc = args.noise_proc
num_null_iters = args.num_null_iters
dataset_ID = args.dataset_ID

# dataset_ID = "UCLA_CNP"
# data_path = "/headnode1/abry4213/data/UCLA_CNP/"
# metadata_file = "UCLA_CNP_sample_metadata.feather"
# SPI_directionality_file = "/headnode1/abry4213/github/fMRI_FeaturesDisorders/classification_analysis/pairwise_analysis/SPI_Direction_Info.csv"
# comparison_group = "Schizophrenia"
# univariate_feature_set = "catch22"
# pairwise_feature_set = "pyspi14"
# pairwise_feature_file ="/headnode1/abry4213/data/UCLA_CNP/processed_data/UCLA_CNP_AROMA_2P_GMR_pyspi14_filtered_zscored.feather"
# noise_proc = "AROMA+2P+GMR"
# num_null_iters = 2

###############################################################################
# Function definitions
###############################################################################

def run_k_fold_SVM_for_feature(feature_data, 
                               feature_list,
                               grouping_var_name,
                               scoring_method,
                               sample_and_class_df,
                               class_labels):
        
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
            
            # Find splits
            splits = list(skf.split(feature_data, class_labels))
            fold_split_list = []
            
            # Fit SVM with 10-fold cross-validation
            cv_results = cross_validate(SVM_model, feature_data, class_labels, 
                                        cv=skf, scoring=scoring_method,
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
            CV_pred = cross_val_predict(SVM_model, feature_data, class_labels, cv=skf)
            sample_and_class_df["CV_Predicted_Diagnosis"] = CV_pred
            sample_and_class_df["Repeat_Number"] = i+1
            CV_sample_predictions_list.append(sample_and_class_df)
            
        # Concatenate results and save per ROI
        fold_assignments = pd.concat(fold_assignments_list)
        fold_assignments["group_var"] = grouping_var_name
        # fold_assignments_per_ROI_list.append(fold_assignments)
        
        SVM_coefficients = pd.concat(SVM_coefficients_list)
        SVM_coefficients["group_var"] = grouping_var_name
        # SVM_coefficients_per_ROI_list.append(SVM_coefficients)
        
        balanced_accuracy = pd.concat(balanced_accuracy_list)
        balanced_accuracy["group_var"] = grouping_var_name
        # balanced_accuracy_per_ROI_list.append(balanced_accuracy)
        
        CV_sample_predictions = pd.concat(CV_sample_predictions_list)
        CV_sample_predictions["group_var"] = grouping_var_name
        
        return (fold_assignments,
                SVM_coefficients,
                balanced_accuracy,
                CV_sample_predictions)

def run_nulls(feature_data, 
              num_iters,
              grouping_type,
                               grouping_var_name,
                               scoring_method,
                               class_labels):
        
        # Fit the SVM
        SVM_model = svm.SVC(kernel="linear", C=1, shrinking=False, class_weight="balanced")
        
        # Define list for null balanced accuracy
        null_balanced_accuracy_list = []
        
        # Shuffle labels and run SVM per null iteration
        for i in range(int(num_iters)):
            
            # Shuffle class labels
            class_labels_shuffled = random.sample(class_labels, len(class_labels))
            
            # Initialize CV splitter
            skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=i)
            
            # Fit SVM with 10-fold cross-validation
            cv_results_null_iter = np.nanmean(cross_val_score(SVM_model,
                                                   feature_data,
                                                   class_labels_shuffled,
                                                   cv = skf,
                                                   scoring = scoring_method))
        
            null_balanced_accuracy_list.append(cv_results_null_iter)
            
        null_balanced_accuracy = pd.DataFrame(null_balanced_accuracy_list, columns=["Balanced_Accuracy"])
        null_balanced_accuracy["group_var"] = grouping_var_name
        null_balanced_accuracy["Grouping_Type"] = grouping_type
        null_balanced_accuracy["Null_Iteration"] = range(int(num_iters))
        
        return null_balanced_accuracy
    
def run_pairwise_SVM(pairwise_feature_file,
                     SPI_directionality_file,
                       univariate_feature_set, 
                       pairwise_feature_set,
                       dataset_ID,
                       metadata_file,
                       comparison_to_control_group,
                       pydata_path,
                       noise_proc,
                       num_null_iters):

    # Define noise label
    noise_label = noise_proc.replace("+", "_")
    
    # Read in directionality data
    SPI_directionality_data = pd.read_csv(SPI_directionality_file)
    SPI_directionality_dict = dict(SPI_directionality_data.values)

    # Load metadata
    metadata = pd.read_feather(data_path + "study_metadata/" + metadata_file)

    # Load in data containing subjects with both univariate and pairwise data available
    samples_to_keep = pd.read_feather(f"{pydata_path}/{dataset_ID}_filtered_sample_info_{noise_label}_{univariate_feature_set}_{pairwise_feature_set}.feather")                                                                           
    
    # Pairwise feature data
    pairwise_feature_data = pd.read_feather(pairwise_feature_file).merge(metadata, on='Sample_ID', how='left').drop(["Age", "Sex"],
                                                                               axis = 1)

    # Filter univariate data by samples with both univariate and pairwise
    # Filter by samples with univariate data available as well
    pairwise_feature_data = pairwise_feature_data[pairwise_feature_data.Sample_ID.isin(samples_to_keep.Sample_ID)]                                                                           

    # Initialise lists for results
    fold_assignments_list = []
    SVM_coefficients_list = []
    balanced_accuracy_list = []
    CV_sample_predictions_list = []
    null_balanced_accuracy_list = []
    
    ###########################################################################
    # SPI-wise

    for this_SPI in pairwise_feature_data.SPI.unique().tolist():
        
        # Subset data to SPI
        SPI_data = pairwise_feature_data.query("SPI == @this_SPI & Diagnosis in ['Control', @comparison_to_control_group]").drop(["SPI"], axis=1)
        
        # Find directionality of SPI
        SPI_directionality = SPI_directionality_dict[this_SPI]
        
        # Merge brain regions according to directionality
        if SPI_directionality == "Directed":
            SPI_data["region_pair"] = SPI_data.brain_region_from + "_" + SPI_data.brain_region_to
            SPI_data = SPI_data.drop(["brain_region_from", "brain_region_to"], axis=1)
        else:
            SPI_data_sorted = [sorted(pair) for pair in SPI_data[["brain_region_from", "brain_region_to"]].values.tolist()]
            SPI_data['region_pair'] = ['_'.join(string) for string in SPI_data_sorted]
            SPI_data = (SPI_data
                        .drop(["brain_region_from", "brain_region_to"], axis=1)
                        .drop_duplicates(ignore_index=True,
                                                                   subset=['Sample_ID', 'region_pair'])
                        )
        
        # Pivot from long to wide
        SPI_data_wide = SPI_data.pivot(index=['Sample_ID', 'Diagnosis'], columns='region_pair', values='value')
        
        # Extract name of features
        feature_list = SPI_data_wide.columns.tolist()
        
        # Extract sample ID and diagnosis as lists
        index_data = SPI_data_wide.index.to_frame().reset_index(drop=True)
        class_labels = index_data["Diagnosis"].tolist()
        
        # Impute any NaN with column mean
        SPI_data_imputed = SPI_data_wide.fillna(SPI_data_wide.mean())
        
        # Extract only the feature data
        features_only = SPI_data_imputed.reset_index(drop=True).to_numpy()
        
        # Run main SVM
        (fold_assignments, SVM_coefficients, balanced_accuracy, CV_sample_predictions) = run_k_fold_SVM_for_feature(feature_data = features_only, 
                                       feature_list = feature_list,
                                       grouping_var_name = this_SPI,
                                       scoring_method = "balanced_accuracy",
                                       sample_and_class_df = index_data,
                                       class_labels = class_labels)
        
        # Run empirical null model permutations
        null_balacc = run_nulls(feature_data=features_only, 
                                num_iters = num_null_iters,
                                grouping_var_name=this_SPI,
                                grouping_type = "SPI",
                                scoring_method = "balanced_accuracy",
                                class_labels = class_labels)
        
        # Save to list of dataframes
        fold_assignments_list.append(fold_assignments)
        SVM_coefficients_list.append(SVM_coefficients)
        balanced_accuracy_list.append(balanced_accuracy)
        CV_sample_predictions_list.append(CV_sample_predictions)
        null_balanced_accuracy_list.append(null_balacc)
    
    ###########################################################################
    # Merge + save results
    fold_assignments_res = pd.concat(fold_assignments_list).reset_index()
    SVM_coefficients_res = pd.concat(SVM_coefficients_list).reset_index()
    balanced_accuracy_res = pd.concat(balanced_accuracy_list).reset_index()
    CV_sample_predictions_res = pd.concat(CV_sample_predictions_list).reset_index()
    null_balanced_accuracy_res = pd.concat(null_balanced_accuracy_list).reset_index()
    
    # Add comparison group info
    fold_assignments_res["Comparison_Group"] = comparison_to_control_group
    SVM_coefficients_res["Comparison_Group"] = comparison_to_control_group
    balanced_accuracy_res["Comparison_Group"] = comparison_to_control_group
    CV_sample_predictions_res["Comparison_Group"] = comparison_to_control_group
    null_balanced_accuracy_res["Comparison Group"] = comparison_to_control_group
        
    # Save results
    feather.write_feather(fold_assignments_res, pydata_path + f"{dataset_ID}_{comparison_to_control_group}_Pairwise_{pairwise_feature_set}_SVM_fold_assigments.feather", version=1)
    feather.write_feather(SVM_coefficients_res, pydata_path + f"{dataset_ID}_{comparison_to_control_group}_Pairwise_{pairwise_feature_set}_SVM_fold_SVM_coefficients.feather", version=1)
    feather.write_feather(balanced_accuracy_res, pydata_path + f"{dataset_ID}_{comparison_to_control_group}_Pairwise_{pairwise_feature_set}_SVM_balanced_accuracy.feather", version=1)
    feather.write_feather(CV_sample_predictions_res, pydata_path + f"{dataset_ID}_{comparison_to_control_group}_Pairwise_{pairwise_feature_set}_SVM_sample_predictions.feather", version=1)
    feather.write_feather(null_balanced_accuracy_res, pydata_path + f"{dataset_ID}_{comparison_to_control_group}_Pairwise_{pairwise_feature_set}_SVM_null_balanced_accuracy.feather", version=1)

                               
###############################################################################
# Main analysis
###############################################################################

run_pairwise_SVM(pairwise_feature_file=pairwise_feature_file,
                 SPI_directionality_file = SPI_directionality_file,
                 univariate_feature_set=univariate_feature_set,
                 pairwise_feature_set=pairwise_feature_set,
                       num_null_iters=num_null_iters,
                       noise_proc = noise_proc,
                       dataset_ID=dataset_ID,
                       metadata_file=metadata_file,
                       comparison_to_control_group=comparison_group,
                       pydata_path=data_path + "processed_data/")

