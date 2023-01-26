
import pandas as pd
import numpy as np
import random
import pyarrow.feather as feather
from sklearn import svm
import os.path
from sklearn.model_selection import StratifiedKFold, cross_val_predict, cross_validate, cross_val_score, permutation_test_score

def run_k_fold_SVM_for_feature(feature_data, 
                                feature_list,
                                grouping_var_name,
                                scoring_method,
                                sample_and_class_df,
                                class_labels,
                                num_null_iters):
        
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
        
        SVM_coefficients = pd.concat(SVM_coefficients_list)
        SVM_coefficients["group_var"] = grouping_var_name
        
        balanced_accuracy = pd.concat(balanced_accuracy_list)
        balanced_accuracy["group_var"] = grouping_var_name
        
        CV_sample_predictions = pd.concat(CV_sample_predictions_list)
        CV_sample_predictions["group_var"] = grouping_var_name
        
        # Run 10-fold CV null model permutations with num_null_iters permutations
        null_model_balanced_accuracy_list = permutation_test_score(SVM_model,
                                                feature_data,
                                                class_labels,
                                                cv=10,
                                                random_state=27,
                                                n_jobs = 16,
                                                n_permutations = num_null_iters,
                                                scoring=scoring_method)[1]
        null_model_balanced_accuracy = pd.DataFrame(null_model_balanced_accuracy_list,
                                                    columns = ["Null_Balanced_Accuracy"])
        null_model_balanced_accuracy["group_var"] = grouping_var_name
        null_model_balanced_accuracy["Null_Iter"] = [*range(1, num_null_iters + 1, 1)]
        
        return (fold_assignments,
                SVM_coefficients,
                balanced_accuracy,
                CV_sample_predictions,
                null_model_balanced_accuracy)

def run_univariate_SVM(univariate_feature_file,
                       univariate_feature_set, 
                       pairwise_feature_set,
                       dataset_ID,
                       metadata_file,
                       comparison_to_control_group,
                       data_path,
                       pydata_path,
                       noise_proc,
                       num_null_iters,
                       overwrite=False):

    # Check if file already exists or overwrite flag is set
    if not os.path.isfile(f"{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_SVM_balanced_accuracy.feather") or overwrite:
        noise_label = noise_proc.replace("+", "_")

        # Load metadata
        metadata = pd.read_feather(data_path + "study_metadata/" + metadata_file)

        # Load in data containing subjects with both univariate and pairwise data available
        samples_to_keep = pd.read_feather(f"{pydata_path}/{dataset_ID}_filtered_sample_info_{noise_label}_{univariate_feature_set}_{pairwise_feature_set}.feather")                                                                           
        
        # Univariate feature data
        univariate_feature_data = pd.read_feather(univariate_feature_file).merge(metadata, on='Sample_ID', how='left').drop(["Age", "Sex"],
                                                                                axis = 1)

        # Filter univariate data by samples with both univariate and pairwise
        # Filter by samples with univariate data available as well
        univariate_feature_data = univariate_feature_data[univariate_feature_data.Sample_ID.isin(samples_to_keep.Sample_ID)]                                                                           

        # Initialise lists for results
        fold_assignments_list = []
        SVM_coefficients_list = []
        balanced_accuracy_list = []
        CV_sample_predictions_list = []
        null_balanced_accuracy_list = []
        
        ###########################################################################
        # Region-wise

        for ROI in univariate_feature_data.Brain_Region.unique().tolist():
            
            # Subset data to ROI
            region_data = univariate_feature_data.query("Brain_Region == @ROI & Diagnosis in ['Control', @comparison_to_control_group]").drop(["Brain_Region", "Noise_Proc",
                                                                                    "method"], axis=1)
            
            # Pivot from long to wide
            region_data_wide = region_data.pivot(index=['Sample_ID', 'Diagnosis'], columns='names', values='values')
            
            # Extract name of features
            feature_list = region_data_wide.columns.tolist()
            
            # Extract sample ID and diagnosis as lists
            index_data = region_data_wide.index.to_frame().reset_index(drop=True)
            class_labels = index_data["Diagnosis"].tolist()
            
            # Extract only the feature data
            features_only = region_data_wide.reset_index(drop=True).to_numpy()
            
            # Run SVM
            (fold_assignments, SVM_coefficients, balanced_accuracy, CV_sample_predictions, null_balacc_dist) = run_k_fold_SVM_for_feature(feature_data = features_only, 
                                        feature_list = feature_list,
                                        grouping_var_name = ROI,
                                        scoring_method = "balanced_accuracy",
                                        sample_and_class_df = index_data,
                                        class_labels = class_labels,
                                        num_null_iters = num_null_iters)
            
            # Save to list of dataframes
            fold_assignments_list.append(fold_assignments)
            SVM_coefficients_list.append(SVM_coefficients)
            balanced_accuracy_list.append(balanced_accuracy)
            CV_sample_predictions_list.append(CV_sample_predictions)
            null_balanced_accuracy_list.append(null_balacc_dist)
            
        ###########################################################################
        # TS Feature-wise
        for TS_feature in univariate_feature_data.names.unique().tolist():
            
            # Subset data to TS feature
            TS_feature_data = univariate_feature_data.query("names == @TS_feature & Diagnosis in ['Control', @comparison_to_control_group]").drop(["names", "Noise_Proc",
                                                                                    "method"], axis=1)
            
            # Pivot from long to wide
            TS_feature_data_wide = TS_feature_data.pivot(index=['Sample_ID', 'Diagnosis'], columns='Brain_Region', values='values')
            
            # Extract name of features
            region_list = TS_feature_data_wide.columns.tolist()
            
            # Extract sample ID and diagnosis as lists
            index_data = TS_feature_data_wide.index.to_frame().reset_index(drop=True)
            class_labels = index_data["Diagnosis"].tolist()
            
            # Extract only the feature data
            features_only = TS_feature_data_wide.reset_index(drop=True).to_numpy()
            
            (fold_assignments, SVM_coefficients, balanced_accuracy, CV_sample_predictions, null_balacc_dist) = run_k_fold_SVM_for_feature(feature_data = features_only, 
                                        feature_list = region_list,
                                        grouping_var_name = TS_feature,
                                        scoring_method = "balanced_accuracy",
                                        sample_and_class_df = index_data,
                                        class_labels = class_labels,
                                        num_null_iters = num_null_iters)
            
            # Save to list of dataframes
            fold_assignments_list.append(fold_assignments)
            SVM_coefficients_list.append(SVM_coefficients)
            balanced_accuracy_list.append(balanced_accuracy)
            CV_sample_predictions_list.append(CV_sample_predictions)
            null_balanced_accuracy_list.append(null_balacc_dist)
            
        ###########################################################################
        # Combo-wise
        
        # Merge brain region + univariate feature name
        combo_data = (univariate_feature_data
                    .query("Diagnosis in ['Control', @comparison_to_control_group]")
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
        class_labels = index_data["Diagnosis"].tolist()
        
        (fold_assignments, SVM_coefficients, balanced_accuracy, CV_sample_predictions, null_balacc_dist) = run_k_fold_SVM_for_feature(feature_data = features_only, 
                                    feature_list = combo_features,
                                    grouping_var_name = "Combo",
                                    scoring_method = "balanced_accuracy",
                                    sample_and_class_df = index_data,
                                    class_labels = class_labels,
                                    num_null_iters = num_null_iters)
        
        # Save to list of dataframes
        fold_assignments_list.append(fold_assignments)
        SVM_coefficients_list.append(SVM_coefficients)
        balanced_accuracy_list.append(balanced_accuracy)
        CV_sample_predictions_list.append(CV_sample_predictions)
        null_balanced_accuracy_list.append(null_balacc_dist)
        
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
        null_balanced_accuracy_res["Comparison_Group"] = comparison_to_control_group
            
        # Save results
        fold_assignments_res.to_feather(pydata_path + f"{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_SVM_fold_assignments.feather")
        SVM_coefficients_res.to_feather(pydata_path + f"{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_SVM_fold_SVM_coefficients.feather")
        balanced_accuracy_res.to_feather(pydata_path + f"{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_SVM_balanced_accuracy.feather")
        CV_sample_predictions_res.to_feather(pydata_path + f"{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_SVM_sample_predictions.feather")
        null_balanced_accuracy_res.to_feather(pydata_path + f"{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_SVM_null_balanced_accuracy_distributions.feather")
   
def run_pairwise_SVM(pairwise_feature_file,
                     SPI_directionality_file,
                       univariate_feature_set, 
                       pairwise_feature_set,
                       dataset_ID,
                       metadata_file,
                       comparison_to_control_group,
                       pydata_path,
                       data_path,
                       noise_proc,
                       num_null_iters,
                       overwrite=False):
    

    # Check if file already exists or overwrite flag is set
    if not os.path.isfile(pydata_path + f"{dataset_ID}_{comparison_to_control_group}_Pairwise_{pairwise_feature_set}_SVM_balanced_accuracy.feather") or overwrite:

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
            (fold_assignments, SVM_coefficients, balanced_accuracy, CV_sample_predictions, null_balacc_dist) = run_k_fold_SVM_for_feature(feature_data = features_only, 
                                        feature_list = feature_list,
                                        grouping_var_name = this_SPI,
                                        scoring_method = "balanced_accuracy",
                                        sample_and_class_df = index_data,
                                        class_labels = class_labels,
                                        num_null_iters = num_null_iters)

            # Save to list of dataframes
            fold_assignments_list.append(fold_assignments)
            SVM_coefficients_list.append(SVM_coefficients)
            balanced_accuracy_list.append(balanced_accuracy)
            CV_sample_predictions_list.append(CV_sample_predictions)
            null_balanced_accuracy_list.append(null_balacc_dist)
        
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
        null_balanced_accuracy_res["Comparison_Group"] = comparison_to_control_group
            
        # Save results
        fold_assignments_res.to_feather(pydata_path + f"{dataset_ID}_{comparison_to_control_group}_Pairwise_{pairwise_feature_set}_SVM_fold_assignments.feather")
        SVM_coefficients_res.to_feather(pydata_path + f"{dataset_ID}_{comparison_to_control_group}_Pairwise_{pairwise_feature_set}_SVM_fold_SVM_coefficients.feather")
        balanced_accuracy_res.to_feather(pydata_path + f"{dataset_ID}_{comparison_to_control_group}_Pairwise_{pairwise_feature_set}_SVM_balanced_accuracy.feather")
        CV_sample_predictions_res.to_feather(pydata_path + f"{dataset_ID}_{comparison_to_control_group}_Pairwise_{pairwise_feature_set}_SVM_sample_predictions.feather")
        null_balanced_accuracy_res.to_feather(pydata_path + f"{dataset_ID}_{comparison_to_control_group}_Pairwise_{pairwise_feature_set}_SVM_null_balanced_accuracy_distributions.feather")
    
    
def run_combined_SVM(univariate_feature_file,
        pairwise_feature_file,
                     SPI_directionality_file,
                       univariate_feature_set, 
                       pairwise_feature_set,
                       dataset_ID,
                       metadata_file,
                       comparison_to_control_group,
                       pydata_path,
                       data_path,
                       noise_proc,
                       num_null_iters,
                       overwrite=False):

    # Check if file already exists or overwrite flag is set
    if not os.path.isfile(pydata_path + f"{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_Pairwise_{pairwise_feature_set}_SVM_balanced_accuracy.feather") or overwrite:
        # Define noise label
        noise_label = noise_proc.replace("+", "_")
        
        # Read in directionality data
        SPI_directionality_data = pd.read_csv(SPI_directionality_file)
        SPI_directionality_dict = dict(SPI_directionality_data.values)

        # Load metadata
        metadata = pd.read_feather(data_path + "study_metadata/" + metadata_file)

        # Load in data containing subjects with both univariate and pairwise data available
        samples_to_keep = pd.read_feather(f"{pydata_path}/{dataset_ID}_filtered_sample_info_{noise_label}_{univariate_feature_set}_{pairwise_feature_set}.feather")                                                                           
        
        # Load in univariate feature data
        univariate_feature_data = pd.read_feather(univariate_feature_file).merge(metadata, on='Sample_ID', how='left').drop(["Age", "Sex"],
                                                                                axis = 1)
        # Load in pairwise feature data
        pairwise_feature_data = pd.read_feather(pairwise_feature_file).merge(metadata, on='Sample_ID', how='left').drop(["Age", "Sex"],
                                                                                axis = 1)

        # Filter data by samples with both univariate and pairwise
        # Filter by samples with univariate data available as well
        univariate_feature_data = univariate_feature_data[univariate_feature_data.Sample_ID.isin(samples_to_keep.Sample_ID)]
        pairwise_feature_data = pairwise_feature_data[pairwise_feature_data.Sample_ID.isin(samples_to_keep.Sample_ID)]                                                                           

        # Initialise lists for results
        fold_assignments_list = []
        SVM_coefficients_list = []
        balanced_accuracy_list = []
        CV_sample_predictions_list = []
        null_balanced_accuracy_list = []
        
        ###########################################################################
        # Prepare univariate combo data
        univariate_combo_data = (univariate_feature_data
                    .query("Diagnosis in ['Control', @comparison_to_control_group]")
                    .drop(["method", "Noise_Proc"], axis=1))
        univariate_combo_data["Combo_Feature"] = univariate_combo_data.Brain_Region + "_" + univariate_combo_data.names
        univariate_combo_data = univariate_combo_data.drop(["Brain_Region", "names"], axis=1)
        
        # Pivot from long to wide
        univariate_combo_data_wide = univariate_combo_data.pivot(index=["Sample_ID", "Diagnosis"],
                                        columns = "Combo_Feature",
                                        values = "values")
        
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
            
            # Merge widened SPI data with widened univariate combo data
            SPI_combo_data_wide = pd.merge(univariate_combo_data_wide,
                                        SPI_data_wide,
                                        how = "inner",
                                        on = ["Sample_ID", "Diagnosis"])
            
            # Extract name of features
            feature_list = SPI_combo_data_wide.columns.tolist()
            
            # Extract sample ID and diagnosis as lists
            index_data = SPI_combo_data_wide.index.to_frame().reset_index(drop=True)
            class_labels = index_data["Diagnosis"].tolist()
            
            # Impute any NaN with column mean
            SPI_combo_data_imputed = SPI_combo_data_wide.fillna(SPI_combo_data_wide.mean())
            
            # Extract only the feature data
            features_only = SPI_combo_data_imputed.reset_index(drop=True).to_numpy()
            
            # Run main SVM
            (fold_assignments, SVM_coefficients, balanced_accuracy, CV_sample_predictions, null_balacc_dist) = run_k_fold_SVM_for_feature(feature_data = features_only, 
                                        feature_list = feature_list,
                                        grouping_var_name = this_SPI,
                                        scoring_method = "balanced_accuracy",
                                        sample_and_class_df = index_data,
                                        class_labels = class_labels,
                                        num_null_iters = num_null_iters)
            
            # Save to list of dataframes
            fold_assignments_list.append(fold_assignments)
            SVM_coefficients_list.append(SVM_coefficients)
            balanced_accuracy_list.append(balanced_accuracy)
            CV_sample_predictions_list.append(CV_sample_predictions)
            null_balanced_accuracy_list.append(null_balacc_dist)
        
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
        null_balanced_accuracy_res["Comparison_Group"] = comparison_to_control_group
            
        # Save results
        fold_assignments_res.to_feather(pydata_path + f"{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_Pairwise_{pairwise_feature_set}_SVM_fold_assignments.feather")
        SVM_coefficients_res.to_feather(pydata_path + f"{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_Pairwise_{pairwise_feature_set}_SVM_fold_SVM_coefficients.feather")
        balanced_accuracy_res.to_feather(pydata_path + f"{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_Pairwise_{pairwise_feature_set}_SVM_balanced_accuracy.feather")
        CV_sample_predictions_res.to_feather(pydata_path + f"{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_Pairwise_{pairwise_feature_set}_SVM_sample_predictions.feather")
        null_balanced_accuracy_res.to_feather(pydata_path + f"{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_Pairwise_{pairwise_feature_set}_SVM_null_balanced_accuracy_distributions.feather")
      