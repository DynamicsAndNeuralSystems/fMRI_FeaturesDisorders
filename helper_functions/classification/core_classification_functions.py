
import pandas as pd
from sklearn import svm
from sklearn import metrics
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, RobustScaler
import os.path
from sklearn.model_selection import StratifiedKFold, cross_val_predict, cross_validate, permutation_test_score
import numpy as np
from mixed_sigmoid_normalisation import MixedSigmoidScaler


def run_k_fold_SVM_for_feature(feature_data, 
                               grouping_var_name,
                               analysis_type,
                               sample_IDs,
                               class_labels,
                               scaling_type,
                               kernel = "linear",
                               num_folds = 10,
                               num_repeats = 10,
                               num_jobs = 10):
        
    
    print(f"Creating pipeline with {scaling_type} scaler.")

    # Define the pipeline -- one for binary predictions, one for probability predictions for ROC
    if scaling_type == "mixedsigmoid":
        pipe = Pipeline([('scaler', MixedSigmoidScaler(unit_variance=True)), 
                         ('SVM', svm.SVC(kernel = kernel, C = 1, 
                         class_weight = "balanced"))])
        prob_pipe = Pipeline([('scaler',  MixedSigmoidScaler(unit_variance=True)), 
                              ('SVM', svm.SVC(kernel=kernel, C=1, 
                             class_weight="balanced", probability=True))])
    else: 
        pipe = Pipeline([('scaler', StandardScaler()), 
                         ('SVM', svm.SVC(kernel=kernel, C=1, 
                                     class_weight="balanced"))])
        prob_pipe = Pipeline([('scaler',  StandardScaler()), 
                              ('SVM', svm.SVC(kernel=kernel, C=1, 
                                 class_weight="balanced", probability=True))])
    
    # Define lists for: 
    # (1) balanced accuracy and AUC by fold/repeat,
    # (2) fold assignments by repeat,
    # (3) TPR/FPR by fold/repeat
    # (4) The proportion of folds for which a subject was predicted correctly
    test_metrics_by_fold_list = []
    fold_assignments_list = []
    TPR_FPR_list = []
    CV_sample_predictions_list = []
    
    for i in range(num_repeats):
        # Initialise sample and class dataframe for repeat
        skf = StratifiedKFold(n_splits=num_folds, shuffle=True, random_state=i)
        
        # Find splits
        splits = list(skf.split(feature_data, class_labels))
        fold_split_list = []
        
        # Fit SVM with 10-fold cross-validation for balanced accuracy
        cv_results_balacc = cross_validate(pipe, feature_data, class_labels, 
                                    cv=skf, scoring=["balanced_accuracy"],
                                    n_jobs = int(num_jobs),
                                    return_estimator=True)
        
        # Define a probabilistic SVM for ROC AUC
        CV_prob = cross_val_predict(prob_pipe, feature_data, class_labels, cv=skf, method="predict_proba")
        ROC_AUC_by_fold = []
        
        # Iterate over folds to save participant assignments and TPR/FPR
        for f in range(num_folds):
            # Split for fold number
            test_indices = splits[f][1]
            test_sample_IDs = [sample_IDs[index] for index in test_indices]
            test_indices_df = pd.DataFrame(test_indices, columns=["Sample_Index"])
            test_indices_df["Sample_ID"] = test_sample_IDs
            test_indices_df["Fold"] = f+1
            test_indices_df["Repeat"] = i+1        
            fold_split_list.append(test_indices_df)

            # TPR and FPR to plot ROC
            y_true_fold = [class_labels[index] for index in test_indices]
            y_pred_fold = np.asarray(CV_prob[test_indices,1])
            roc_auc_for_fold = metrics.roc_auc_score(y_true_fold, y_pred_fold)
            ROC_AUC_by_fold.append(roc_auc_for_fold)
        
        # Combine lists into dataframes
        fold_splits = pd.concat(fold_split_list)
        fold_splits["Analysis_Type"] = analysis_type
        fold_assignments_list.append(fold_splits)

        # FPR/TPR for ROC
        y_true = np.asarray([int(u==1) for u in class_labels])
        y_pred = np.asarray(CV_prob[:,1])
        fpr, tpr, _ = metrics.roc_curve(y_true, y_pred, drop_intermediate=False)
        TPR_FPR_data = pd.DataFrame({"fpr": fpr,
                                 "tpr": tpr,
                                 "Repeat_Number": i+1})
        
        # Append repeat results
        TPR_FPR_list.append(TPR_FPR_data)

        # Extract balanced accuracy and ROC AUC by fold
        test_metrics_by_fold_df = pd.DataFrame({"group_var": grouping_var_name,
                                                "Analysis_Type": analysis_type,
                                        "Fold": [*range(1, num_folds + 1, 1)],
                                        "Repeat_Number": i,
                                        "Balanced_Accuracy": cv_results_balacc["test_balanced_accuracy"],
                                        "ROC_AUC": ROC_AUC_by_fold})
        test_metrics_by_fold_list.append(test_metrics_by_fold_df)

        # Find the proportion of folds for which a subject was predicted correctly as a DataFrame
        CV_pred = cross_val_predict(pipe, feature_data, class_labels, cv=skf)
        CV_pred_df = pd.DataFrame({"Sample_ID": sample_IDs,
                                                         "True_Diagnosis": class_labels})
        CV_pred_df["CV_Predicted_Diagnosis"] = CV_pred
        CV_pred_df["Repeat_Number"] = i+1
        CV_sample_predictions_list.append(CV_pred_df)
    
    test_metrics_by_fold = pd.concat(test_metrics_by_fold_list)
    test_metrics_by_fold["group_var"] = grouping_var_name

    fold_assignments = pd.concat(fold_assignments_list)
    fold_assignments["group_var"] = grouping_var_name

    TPR_FPR = pd.concat(TPR_FPR_list)
    TPR_FPR["group_var"] = grouping_var_name

    # Find the proportion of repeats for which each subject was correctly predicted as a DataFrame
    CV_sample_predictions = pd.concat(CV_sample_predictions_list)
    CV_sample_predictions["group_var"] = grouping_var_name

    # Group by 'Sample_ID' and calculate the proportion of matches
    prop_correct_predictions = CV_sample_predictions.groupby('Sample_ID').apply(lambda x: pd.Series({
        'True_Diagnosis': x['True_Diagnosis'].iloc[0],
        'Proportion_Match': (x['True_Diagnosis'] == x['CV_Predicted_Diagnosis']).mean()
    })).reset_index()
    prop_correct_predictions["Analysis_Type"] = analysis_type

    return (test_metrics_by_fold,
            TPR_FPR,
            fold_assignments,
            prop_correct_predictions)

def run_univariate_SVM(univariate_feature_file,
                       univariate_feature_set, 
                       pairwise_feature_set,
                       dataset_ID,
                       metadata_file,
                       comparison_to_control_group,
                       data_path,
                       pydata_path,
                       noise_proc,
                       kernel="linear",
                       num_folds = 10,
                       scaling_type="mixedsigmoid",
                       num_jobs = 8,
                       num_repeats = 10,
                       overwrite=False):

    # Check if file already exists or overwrite flag is set
    if not os.path.isfile(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_{scaling_type}_scaler_{kernel}_SVM_balanced_accuracy_AUC.feather") or overwrite:
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
        test_metrics_list = []
        TPR_FPR_list = []
        prop_correct_predictions_list = prop_correct_predictions
        
        ###########################################################################
        # Region-wise

        for ROI in univariate_feature_data.Brain_Region.unique().tolist():
            
            # Subset data to ROI
            region_data = univariate_feature_data.query("Brain_Region == @ROI & Diagnosis in ['Control', @comparison_to_control_group]").drop(["Brain_Region", "Noise_Proc",
                                                                                    "method"], axis=1)
            
            # Pivot from long to wide
            region_data_wide = region_data.pivot(index=['Sample_ID', 'Diagnosis'], columns='names', values='values')
            
            # Extract sample ID and diagnosis as lists
            index_data = region_data_wide.index.to_frame().reset_index(drop=True)
            class_labels = [int(i==comparison_to_control_group) for i in index_data["Diagnosis"].tolist()]
            sample_IDs = index_data["Sample_ID"].tolist()

            # Extract only the feature data
            features_only = region_data_wide.reset_index(drop=True).to_numpy()
            
            # Run SVM
            (test_metrics_by_fold, TPR_FPR, fold_assignments, prop_correct_predictions) = run_k_fold_SVM_for_feature(feature_data = features_only, 
                grouping_var_name = ROI,
                analysis_type = "Univariate_Brain_Region",
                sample_IDs = sample_IDs,
                class_labels = class_labels,
                scaling_type = scaling_type,
                kernel = kernel,
                num_folds = num_folds,
                num_jobs = num_jobs,
                num_repeats = num_repeats)
            
            # Save to list of dataframes
            test_metrics_list.append(test_metrics_by_fold)
            TPR_FPR_list.append(TPR_FPR)
            fold_assignments_list.append(fold_assignments)
            prop_correct_predictions_list.append(prop_correct_predictions)
            
        ###########################################################################
        # TS Feature-wise
        for TS_feature in univariate_feature_data.names.unique().tolist():
            
            # Subset data to TS feature
            TS_feature_data = univariate_feature_data.query("names == @TS_feature & Diagnosis in ['Control', @comparison_to_control_group]").drop(["names", "Noise_Proc",
                                                                                    "method"], axis=1)
            
            # Pivot from long to wide
            TS_feature_data_wide = TS_feature_data.pivot(index=['Sample_ID', 'Diagnosis'], columns='Brain_Region', values='values')
            
            # Extract sample ID and diagnosis as lists
            index_data = TS_feature_data_wide.index.to_frame().reset_index(drop=True)
            class_labels = [int(i==comparison_to_control_group) for i in index_data["Diagnosis"].tolist()]
            sample_IDs = index_data["Sample_ID"].tolist()
            
            # Extract only the feature data
            features_only = TS_feature_data_wide.reset_index(drop=True).to_numpy()
            
            # Run SVM
            (test_metrics_by_fold, TPR_FPR, fold_assignments, prop_correct_predictions) = run_k_fold_SVM_for_feature(feature_data = features_only, 
                grouping_var_name = TS_feature,
                analysis_type = "Univariate_TS_Feature",
                sample_IDs = sample_IDs,
                class_labels = class_labels,
                scaling_type = scaling_type,
                kernel = kernel,
                num_folds = num_folds,
                num_jobs = num_jobs,
                num_repeats = num_repeats)
            
            # Save to list of dataframes
            test_metrics_list.append(test_metrics_by_fold)
            TPR_FPR_list.append(TPR_FPR)
            fold_assignments_list.append(fold_assignments)
            prop_correct_predictions_list.append(prop_correct_predictions)
            
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
        
        # Extract only the combo feature data
        features_only = combo_data_wide.reset_index(drop=True).to_numpy()
        
        # Extract sample ID and diagnosis as lists
        index_data = combo_data_wide.index.to_frame().reset_index(drop=True)
        class_labels = [int(i==comparison_to_control_group) for i in index_data["Diagnosis"].tolist()]
        sample_IDs = index_data["Sample_ID"].tolist()
        
        # Run SVM
        (test_metrics_by_fold, TPR_FPR, fold_assignments, prop_correct_predictions) = run_k_fold_SVM_for_feature(feature_data = features_only, 
                    sample_IDs = sample_IDs,
                    class_labels = class_labels,
                    scaling_type = scaling_type,
                    kernel = kernel,
                    num_folds = num_folds,
                    num_jobs = num_jobs,
                    num_repeats = num_repeats)
        
        # Save to list of dataframes
        test_metrics_list.append(test_metrics_by_fold)
        TPR_FPR_list.append(TPR_FPR)
        fold_assignments_list.append(fold_assignments)
        prop_correct_predictions_list.append(prop_correct_predictions)
        
        ###########################################################################
        # Merge + save results
        test_metrics_res = pd.concat(test_metrics_list).reset_index()
        TPR_FPR_res = pd.concat(TPR_FPR_list).reset_index()
        fold_assignments_res = pd.concat(fold_assignments_list).reset_index()
        prop_correct_predictions_res = pd.concat(prop_correct_predictions_list).reset_index()
        
        # Add comparison group info and normalisation method info
        test_metrics_res["Comparison_Group"] = comparison_to_control_group
        test_metrics_res["Scaling_Type"] = scaling_type
        test_metrics_res["Kernel"] = kernel

        TPR_FPR_res["Comparison_Group"] = comparison_to_control_group
        TPR_FPR_res["Scaling_Type"] = scaling_type
        TPR_FPR_res["Kernel"] = kernel

        fold_assignments_res["Comparison_Group"] = comparison_to_control_group
        fold_assignments_res["Scaling_Type"] = scaling_type
        fold_assignments_res["Kernel"] = kernel

        prop_correct_predictions_res["Comparison_Group"] = comparison_to_control_group
        prop_correct_predictions_res["Scaling_Type"] = scaling_type
        prop_correct_predictions_res["Kernel"] = kernel

        test_metrics_res.to_feather(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_{scaling_type}_scaler_{kernel}_SVM_balanced_accuracy_AUC.feather")
        TPR_FPR_res.to_feather(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_{scaling_type}_scaler_{kernel}_SVM_ROC_TPR_FPR.feather")
        fold_assignments_res.to_feather(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_{scaling_type}_scaler_{kernel}_SVM_fold_assignments.feather")
        prop_correct_predictions_res.to_feather(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_{scaling_type}_scaler_{kernel}_SVM_prop_correct_pred.feather")
        
def run_pairwise_SVM_by_SPI(pairwise_feature_file,
                     SPI_directionality_file,
                       univariate_feature_set, 
                       pairwise_feature_set,
                       dataset_ID,
                       metadata_file,
                       comparison_to_control_group,
                       pydata_path,
                       data_path,
                       noise_proc,
                       kernel = "linear",
                       num_folds = 10,
                       scaling_type="mixedsigmoid",
                       num_jobs = 10,
                       num_repeats = 10,
                       overwrite=False):
    

    # Check if file already exists or overwrite flag is set
    if not os.path.isfile(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Pairwise_{pairwise_feature_set}_{scaling_type}_scaler_{kernel}_SVM_balanced_accuracy_AUC.feather") or overwrite:

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
        test_metrics_list = []
        TPR_FPR_list = []
        fold_assignments_list = []
        prop_correct_predictions_list = []
        
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
            
            # Extract sample ID and diagnosis as lists
            index_data = SPI_data_wide.index.to_frame().reset_index(drop=True)
            class_labels = [int(i==comparison_to_control_group) for i in index_data["Diagnosis"].tolist()]
            sample_IDs = index_data["Sample_ID"].tolist()
            
            # Impute any NaN with column mean
            SPI_data_imputed = SPI_data_wide.fillna(SPI_data_wide.mean())
            
            # Extract only the feature data
            features_only = SPI_data_imputed.reset_index(drop=True).to_numpy()
            
            # Run main SVM
            (test_metrics_by_fold, TPR_FPR, fold_assignments, prop_correct_predictions) = run_k_fold_SVM_for_feature(feature_data = features_only,
                                        grouping_var_name = this_SPI,
                                        analysis_type = "Pairwise_SPI",
                                        sample_IDs = sample_IDs,
                                        kernel=kernel,
                                        class_labels = class_labels,
                                        scaling_type = scaling_type,
                                        num_folds = num_folds,
                                        num_jobs = num_jobs,
                                        num_repeats = num_repeats)

            # Save to list of dataframes
            test_metrics_list.append(test_metrics_by_fold)
            TPR_FPR_list.append(TPR_FPR)
            fold_assignments_list.append(fold_assignments)
            prop_correct_predictions_list.append(prop_correct_predictions)
            
        ###########################################################################
        # Merge + save results
        test_metrics_res = pd.concat(test_metrics_list).reset_index()
        TPR_FPR_res = pd.concat(TPR_FPR_list).reset_index()
        fold_assignments_res = pd.concat(fold_assignments_list).reset_index()
        prop_correct_predictions_res = pd.concat(prop_correct_predictions_list).reset_index()
        
        # Add comparison group info and normalisation method info
        test_metrics_res["Comparison_Group"] = comparison_to_control_group
        test_metrics_res["Scaling_Type"] = scaling_type
        test_metrics_res["Kernel"] = kernel

        TPR_FPR_res["Comparison_Group"] = comparison_to_control_group
        TPR_FPR_res["Scaling_Type"] = scaling_type
        TPR_FPR_res["Kernel"] = kernel

        fold_assignments_res["Comparison_Group"] = comparison_to_control_group
        fold_assignments_res["Scaling_Type"] = scaling_type
        fold_assignments_res["Kernel"] = kernel

        prop_correct_predictions_res["Comparison_Group"] = comparison_to_control_group
        prop_correct_predictions_res["Scaling_Type"] = scaling_type
        prop_correct_predictions_res["Kernel"] = kernel
            
        # Save results
        fold_assignments_res.to_feather(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Pairwise_{pairwise_feature_set}_{scaling_type}_scaler_{kernel}_SVM_fold_assignments.feather")
        test_metrics_res.to_feather(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Pairwise_{pairwise_feature_set}_{scaling_type}_scaler_{kernel}_SVM_balanced_accuracy_AUC.feather")
        TPR_FPR_res.to_feather(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Pairwise_{pairwise_feature_set}_{scaling_type}_scaler_{kernel}_SVM_ROC_TPR_FPR.feather")
        prop_correct_predictions_res.to_feather(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Pairwise_{pairwise_feature_set}_{scaling_type}_scaler_{kernel}_SVM_prop_correct_pred.feather")

def run_combined_uni_pairwise_SVM_by_SPI(univariate_feature_file,
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
                        kernel = "linear",
                       scaling_type = "mixedsigmoid",
                       num_jobs = 10,
                       num_repeats = 10,
                       num_folds=10,
                       overwrite=False):

    # Check if file already exists or overwrite flag is set
    if not os.path.isfile(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_Pairwise_{pairwise_feature_set}_{scaling_type}_scaler_{kernel}_SVM_balanced_accuracy_AUC.feather") or overwrite:
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
        test_metrics_list = []
        TPR_FPR_list = []
        fold_assignments_list = []
        prop_correct_predictions_list = []
        
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
            
            # Extract sample ID and diagnosis as lists
            index_data = SPI_combo_data_wide.index.to_frame().reset_index(drop=True)
            class_labels = [int(i==comparison_to_control_group) for i in index_data["Diagnosis"].tolist()]
            sample_IDs = index_data["Sample_ID"].tolist()
            
            # Impute any NaN with column mean
            SPI_combo_data_imputed = SPI_combo_data_wide.fillna(SPI_combo_data_wide.mean())
            
            # Extract only the feature data
            features_only = SPI_combo_data_imputed.reset_index(drop=True).to_numpy()
            
            # Run main SVM
            (test_metrics_by_fold, TPR_FPR, fold_assignments, prop_correct_predictions) = run_k_fold_SVM_for_feature(feature_data = features_only, 
                                        grouping_var_name = this_SPI,
                                        analysis_type = "Univariate_Pairwise_Combo",
                                        sample_IDs = sample_IDs,
                                        class_labels = class_labels,
                                        kernel=kernel,
                                        scaling_type = scaling_type,
                                        num_jobs = num_jobs,
                                        num_folds = num_folds,
                                        num_repeats = num_repeats)
            
            # Save to list of dataframes
            test_metrics_list.append(test_metrics_by_fold)
            TPR_FPR_list.append(TPR_FPR)
            fold_assignments_list.append(fold_assignments)
            prop_correct_predictions_list.append(prop_correct_predictions)
        
        ###########################################################################
        # Merge + save results
        test_metrics_res = pd.concat(test_metrics_list).reset_index()
        TPR_FPR_res = pd.concat(TPR_FPR_list).reset_index()
        fold_assignments_res = pd.concat(fold_assignments_list).reset_index()
        prop_correct_predictions_res = pd.concat(prop_correct_predictions).reset_index()
        
        # Add comparison group info and normalisation method info
        test_metrics_res["Comparison_Group"] = comparison_to_control_group
        test_metrics_res["Scaling_Type"] = scaling_type
        test_metrics_res["Kernel"] = kernel

        TPR_FPR_res["Comparison_Group"] = comparison_to_control_group
        TPR_FPR_res["Scaling_Type"] = scaling_type
        TPR_FPR_res["Kernel"] = kernel

        fold_assignments_res["Comparison_Group"] = comparison_to_control_group
        fold_assignments_res["Scaling_Type"] = scaling_type
        fold_assignments_res["Kernel"] = kernel

        prop_correct_predictions_res["Comparison_Group"] = comparison_to_control_group
        prop_correct_predictions_res["Scaling_Type"] = scaling_type
        prop_correct_predictions_res["Kernel"] = kernel
            
        # Save results
        test_metrics_res.to_feather(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_Pairwise_{pairwise_feature_set}_{scaling_type}_scaler_{kernel}_balanced_accuracy_AUC.feather")
        TPR_FPR_res.to_feather(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_Pairwise_{pairwise_feature_set}_{scaling_type}_scaler_{kernel}_SVM_ROC_TPR_FPR.feather")
        fold_assignments_res.to_feather(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_Pairwise_{pairwise_feature_set}_{scaling_type}_scaler_{kernel}_SVM_fold_assignments.feather")
        prop_correct_predictions_res.to_feather(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_Pairwise_{pairwise_feature_set}_{scaling_type}_scaler_{kernel}_SVM_prop_correct_pred.feather")

###################### NULLS ########################

def run_nulls_for_feature(feature_data, 
                               output_file_base,
                               sample_and_class_df,
                               scaling_type,
                               num_folds = 10,
                               num_null_iters = 1000,
                               num_repeats = 10,
                               num_jobs = 10):

    print(f"Creating pipeline with {scaling_type} scaler.")

    # Define the pipeline
    if scaling_type == "robust":
        pipe = Pipeline([('scaler', RobustScaler()), 
                         ('SVM', svm.SVC(kernel="linear", C=1, shrinking=False, 
                                     class_weight="balanced"))])
    elif scaling_type == "mixedsigmoid":
        pipe = Pipeline([('scaler', MixedSigmoidScaler(unit_variance=True)), 
                         ('SVM', svm.SVC(kernel = "linear", C = 1, shrinking = False, 
                         class_weight = "balanced"))])
    else: 
        pipe = Pipeline([('scaler', StandardScaler()), 
                         ('SVM', svm.SVC(kernel="linear", C=1, shrinking=False, 
                                     class_weight="balanced"))])
    
    # Run per number of null iterations
    for j in range(num_null_iters):
        
        # Check if null iter results already exist
        if not os.path.isfile(f"{output_file_base}_null_iter_{str(j)}.feather"):
            
            # print(f"Now running null iter str({j}) for {output_file_base}")
        
            # Create list to store null balanced accuracy results across repeats for given iter
            null_iter_list = []
            
            # For each iteration of nulls, run 10-repeat 10-fold CV
            # Take the mean value across all folds (which is what the permutation_test_score outputs),
            # Then take the mean value across the 10 repeats
            
            # USING MANUAL NULLS
            for n in range(num_repeats):
                # print(f"Repeat number {str(n)}")
                # Initialise sample and class dataframe for repeat
                sample_and_class_df_for_repeat = sample_and_class_df.copy(deep=True)
                
                # Shuffle class labels
                sample_and_class_df_for_repeat["Reshuffled_Label"] = sample_and_class_df_for_repeat.Diagnosis.sample(frac=1).values
                reshuffled_labels = [int(i != "Control") for i in sample_and_class_df_for_repeat.Reshuffled_Label.tolist()]
                
                # Split into ten folds
                skf = StratifiedKFold(n_splits=num_folds, shuffle=True)
                
                # Fit SVM with 10-fold cross-validation
                cv_results = cross_validate(pipe, 
                                            feature_data, 
                                            reshuffled_labels, 
                                            cv=skf, 
                                            scoring="balanced_accuracy",
                                            n_jobs = int(num_jobs),
                                            return_estimator=False)
    
                # Extract balanced accuracy by fold
                null_balanced_accuracy_by_fold_df = pd.DataFrame(cv_results["test_score"],
                                                            columns=["Null_Balanced_Accuracy"])
                null_balanced_accuracy_by_fold_df["Null_Iter"] = j + 1
                null_balanced_accuracy_by_fold_df["Fold"] = [*range(1, num_folds + 1, 1)]
                null_balanced_accuracy_by_fold_df["Repeat_Number"] = n + 1
                
                # Append repeat results to list for the given iter
                null_iter_list.append(null_balanced_accuracy_by_fold_df)
                
            # Combine all results across repeats for the given iter
            null_iter_res = pd.concat(null_iter_list)
            
            # Save to feather
            null_iter_res = null_iter_res.reset_index()
            null_iter_res.to_feather(f"{output_file_base}_null_iter_{str(j)}.feather")

def combine_nulls_for_feature(grouping_var_name,
                               analysis_type,
                               output_file_base,
                               num_null_iters = 1000):

    # Initialise list
    null_balacc_list = []
        
    # Load null data per iter
    for j in range(num_null_iters):
        null_iter_res = pd.read_feather(f"{output_file_base}_null_iter_{str(j)}.feather")


        # Take average balanced accuracy first across folds, then across repeats
        null_iter_res_averaged = (null_iter_res
                                  .groupby("Repeat_Number")["Null_Balanced_Accuracy"]
                                  .mean()
                                  .to_frame()
                                  .reset_index())["Null_Balanced_Accuracy"].mean()
        
        # Add as dataframe
        null_iter_res_averaged_df = pd.DataFrame(data = [[j, null_iter_res_averaged]],
                                                 columns = ["Null_Iter_Number", "Null_Balanced_Accuracy"])
        null_balacc_list.append(null_iter_res_averaged_df)
    
    # Concatenate results across null iters
    null_balacc_res = pd.concat(null_balacc_list)
    null_balacc_res["Analysis_Type"] = analysis_type
    null_balacc_res["group_var"] = grouping_var_name

    # Return the results
    return(null_balacc_res)

def run_univariate_nulls(univariate_feature_file,
                       univariate_feature_set, 
                       pairwise_feature_set,
                       dataset_ID,
                       metadata_file,
                       comparison_to_control_group,
                       data_path,
                       pydata_path,
                       noise_proc,
                       num_null_iters = 1000,
                       num_folds = 10,
                       scaling_type="mixedsigmoid",
                       num_jobs = 10,
                       num_repeats = 10,
                       overwrite=False):
    
    # Try making output directory
    if not os.path.isdir(f"{pydata_path}/null_results"):
        os.mkdir(f"{pydata_path}/null_results")

    # Check if file already exists or overwrite flag is set
    if not os.path.isfile(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_{scaling_type}_scaler_SVM_null_balanced_accuracy_distributions.feather") or overwrite:
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
        null_balanced_accuracy_list = []
        
        ###########################################################################
        # Region-wise

        for ROI in univariate_feature_data.Brain_Region.unique().tolist():
            
            # Subset data to ROI
            region_data = univariate_feature_data.query("Brain_Region == @ROI & Diagnosis in ['Control', @comparison_to_control_group]").drop(["Brain_Region", "Noise_Proc",
                                                                                    "method"], axis=1)
            
            # Pivot from long to wide
            region_data_wide = region_data.pivot(index=['Sample_ID', 'Diagnosis'], columns='names', values='values')
            
            # Extract sample ID and diagnosis as lists
            index_data = region_data_wide.index.to_frame().reset_index(drop=True)
            
            # Extract only the feature data
            features_only = region_data_wide.reset_index(drop=True).to_numpy()
            
            # Run nulls
            run_nulls_for_feature(feature_data = features_only, 
                                    output_file_base = f"{pydata_path}/null_results/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_{scaling_type}_scaler_{ROI}_SVM",
                                    sample_and_class_df = index_data,
                                    scaling_type = scaling_type,
                                    num_null_iters = num_null_iters,
                                    num_folds = num_folds,
                                    num_jobs = num_jobs,
                                    num_repeats = num_repeats)

            # Combine nulls
            null_balanced_accuracies = combine_nulls_for_feature(grouping_var_name = ROI, 
                                                                 analysis_type="Univariate_Brain_Region", 
                                                                 output_file_base = f"{pydata_path}/null_results/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_{scaling_type}_scaler_{ROI}_SVM", 
                                                                 num_null_iters=num_null_iters)

            # Save to list of dataframes
            null_balanced_accuracy_list.append(null_balanced_accuracies)
            
        ###########################################################################
        # TS Feature-wise
        for TS_feature in univariate_feature_data.names.unique().tolist():
            
            # Subset data to TS feature
            TS_feature_data = univariate_feature_data.query("names == @TS_feature & Diagnosis in ['Control', @comparison_to_control_group]").drop(["names", "Noise_Proc",
                                                                                    "method"], axis=1)
            
            # Pivot from long to wide
            TS_feature_data_wide = TS_feature_data.pivot(index=['Sample_ID', 'Diagnosis'], columns='Brain_Region', values='values')
            
            # Extract sample ID and diagnosis as lists
            index_data = TS_feature_data_wide.index.to_frame().reset_index(drop=True)
            
            # Extract only the feature data
            features_only = TS_feature_data_wide.reset_index(drop=True).to_numpy()
            
            # Run nulls
            run_nulls_for_feature(feature_data = features_only, 
                                    output_file_base = f"{pydata_path}/null_results/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_{scaling_type}_scaler_{TS_feature}_SVM",
                                    sample_and_class_df = index_data,
                                    scaling_type = scaling_type,
                                    num_null_iters = num_null_iters,
                                    num_folds = num_folds,
                                    num_jobs = num_jobs,
                                    num_repeats = num_repeats)

            # Combine nulls
            null_balanced_accuracies = combine_nulls_for_feature(grouping_var_name = TS_feature, 
                                                                 analysis_type="Univariate_TS_Feature", 
                                                                 output_file_base = f"{pydata_path}/null_results/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_{scaling_type}_scaler_{TS_feature}_SVM", 
                                                                 num_null_iters=num_null_iters)
            
            # Save to list of dataframes
            null_balanced_accuracy_list.append(null_balanced_accuracies)
            
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
        
        # Extract only the combo feature data
        features_only = combo_data_wide.reset_index(drop=True).to_numpy()
        
        # Extract sample ID and diagnosis as lists
        index_data = combo_data_wide.index.to_frame().reset_index(drop=True)
        
        # Run nulls
        run_nulls_for_feature(feature_data = features_only, 
                                output_file_base = f"{pydata_path}/null_results/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_{scaling_type}_scaler_Combo_SVM",
                                sample_and_class_df = index_data,
                                scaling_type = scaling_type,
                                num_null_iters = num_null_iters,
                                num_folds = num_folds,
                                num_jobs = num_jobs,
                                num_repeats = num_repeats)

        # Combine nulls
        null_balanced_accuracies = combine_nulls_for_feature(grouping_var_name = "Combo", 
                                                                analysis_type="Univariate_Combo", 
                                                                output_file_base = f"{pydata_path}/null_results/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_{scaling_type}_scaler_Combo_SVM", 
                                                                num_null_iters=num_null_iters)
        
        # Save to list of dataframes
        null_balanced_accuracy_list.append(null_balanced_accuracies)
        
        ###########################################################################
        # Merge + save results
        null_balanced_accuracy_res = pd.concat(null_balanced_accuracy_list).reset_index()
        
        # Add comparison group info and normalisation method info
        null_balanced_accuracy_res["Comparison_Group"] = comparison_to_control_group
        null_balanced_accuracy_res["Scaling_Type"] = scaling_type
        
        # Save results
        null_balanced_accuracy_res.to_feather(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_{scaling_type}_scaler_SVM_null_balanced_accuracy_distributions.feather")
        
def run_pairwise_nulls_by_SPI(pairwise_feature_file,
                     SPI_directionality_file,
                       univariate_feature_set, 
                       pairwise_feature_set,
                       dataset_ID,
                       metadata_file,
                       comparison_to_control_group,
                       pydata_path,
                       data_path,
                       noise_proc,
                       num_folds = 10,
                       scaling_type="mixedsigmoid",
                       num_null_iters = 1000,
                       num_jobs = 10,
                       num_repeats = 10,
                       overwrite=False):
    
    # Try making output directory
    if not os.path.isdir(f"{pydata_path}/null_results"):
        os.mkdir(f"{pydata_path}/null_results")

    # Check if file already exists or overwrite flag is set
    if not os.path.isfile(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Pairwise_{pairwise_feature_set}_{scaling_type}_scaler_SVM_null_balanced_accuracy_distributions.feather") or overwrite:
        
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
            
            # Extract sample ID and diagnosis as lists
            index_data = SPI_data_wide.index.to_frame().reset_index(drop=True)
            class_labels = index_data["Diagnosis"].tolist()
            
            # Impute any NaN with column mean
            SPI_data_imputed = SPI_data_wide.fillna(SPI_data_wide.mean())
            
            # Extract only the feature data
            features_only = SPI_data_imputed.reset_index(drop=True).to_numpy()
            
            # Run nulls
            run_nulls_for_feature(feature_data = features_only, 
                                        output_file_base = f"{pydata_path}/null_results/{dataset_ID}_{comparison_to_control_group}_Pairwise_{pairwise_feature_set}_{scaling_type}_scaler_{this_SPI}_SVM",
                                        sample_and_class_df = index_data,
                                        scaling_type = scaling_type,
                                        num_null_iters = num_null_iters,
                                        num_folds = num_folds,
                                        num_jobs = num_jobs,
                                        num_repeats = num_repeats)

            # Combine nulls
            null_balanced_accuracies = combine_nulls_for_feature(grouping_var_name = this_SPI, 
                                                                    analysis_type="Pairwise_SPI", 
                                                                    output_file_base = f"{pydata_path}/null_results/{dataset_ID}_{comparison_to_control_group}_Pairwise_{pairwise_feature_set}_{scaling_type}_scaler_{this_SPI}_SVM", 
                                                                    num_null_iters=num_null_iters)
            
            # Save to list of dataframes
            null_balanced_accuracy_list.append(null_balanced_accuracies)
            
        ###########################################################################
        # Merge + save results
        null_balanced_accuracy_res = pd.concat(null_balanced_accuracy_list).reset_index()
        
        # Add comparison group info and normalisation method info
        null_balanced_accuracy_res["Comparison_Group"] = comparison_to_control_group
        null_balanced_accuracy_res["Scaling_Type"] = scaling_type
        
        # Save results
        null_balanced_accuracy_res.to_feather(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Pairwise_{pairwise_feature_set}_{scaling_type}_scaler_SVM_null_balanced_accuracy_distributions.feather")
        
def run_combined_uni_pairwise_nulls_by_SPI(univariate_feature_file,
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
                       scaling_type = "mixedsigmoid",
                       num_jobs = 10,
                       num_repeats = 10,
                       num_folds=10,
                       num_null_iters=1000,
                       overwrite=False):
    
    # Try making output directory
    if not os.path.isdir(f"{pydata_path}/null_results"):
        os.mkdir(f"{pydata_path}/null_results")

    # Check if file already exists or overwrite flag is set
    if not os.path.isfile(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_Pairwise_{pairwise_feature_set}_{scaling_type}_scaler_SVM_null_balanced_accuracy_distributions.feather") or overwrite:
        
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

        # # Print head of univariate_feature_data
        # print("Head of univariate_feature_data:")
        # print(univariate_feature_data.head())

        # Load in pairwise feature data
        pairwise_feature_data = pd.read_feather(pairwise_feature_file).merge(metadata, on='Sample_ID', how='left').drop(["Age", "Sex"],
                                                                                axis = 1)

        # Filter data by samples with both univariate and pairwise
        # Filter by samples with univariate data available as well
        univariate_feature_data = univariate_feature_data[univariate_feature_data.Sample_ID.isin(samples_to_keep.Sample_ID)]
        pairwise_feature_data = pairwise_feature_data[pairwise_feature_data.Sample_ID.isin(samples_to_keep.Sample_ID)]                                                                           

        # Initialise lists for results
        null_balanced_accuracy_list = []

        ###########################################################################
        # Prepare univariate combo data
        univariate_combo_data = (univariate_feature_data
                    .query("Diagnosis in ['Control', @comparison_to_control_group]"))

        univariate_combo_data["Combo_Feature"] = univariate_combo_data.Brain_Region + "_" + univariate_combo_data.names
        univariate_combo_data = univariate_combo_data.drop(["Brain_Region", "names"], axis=1)
        
        # Pivot from long to wide
        univariate_combo_data_wide = univariate_combo_data.pivot(index=["Sample_ID", "Diagnosis"],
                                        columns = "Combo_Feature",
                                        values = "values")



        
        ###########################################################################
        # SPI-wise
        # for this_SPI in ["sgc_nonparametric_mean_fs-1_fmin-0_fmax-0-5"]:
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
            
            # Extract sample ID and diagnosis as lists
            index_data = SPI_combo_data_wide.index.to_frame().reset_index(drop=True)
            
            # Impute any NaN with column mean
            SPI_combo_data_imputed = SPI_combo_data_wide.fillna(SPI_combo_data_wide.mean())
            
            # Extract only the feature data
            features_only = SPI_combo_data_imputed.reset_index(drop=True).to_numpy()
            
            # Run nulls
            run_nulls_for_feature(feature_data = features_only, 
                                        output_file_base = f"{pydata_path}/null_results/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_Pairwise_{pairwise_feature_set}_{scaling_type}_scaler_{this_SPI}_SVM",
                                        sample_and_class_df = index_data,
                                        scaling_type = scaling_type,
                                        num_null_iters = num_null_iters,
                                        num_folds = num_folds,
                                        num_jobs = num_jobs,
                                        num_repeats = num_repeats)

            # Combine nulls
            null_balanced_accuracies = combine_nulls_for_feature(grouping_var_name = this_SPI, 
                                                                    analysis_type="Univariate_Pairwise_Combo", 
                                                                    output_file_base = f"{pydata_path}/null_results/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_Pairwise_{pairwise_feature_set}_{scaling_type}_scaler_{this_SPI}_SVM", 
                                                                    num_null_iters=num_null_iters)
            
            # Save to list of dataframes
            null_balanced_accuracy_list.append(null_balanced_accuracies)
            
        ###########################################################################
        # Merge + save results
        null_balanced_accuracy_res = pd.concat(null_balanced_accuracy_list).reset_index()
        
        # Add comparison group info and normalisation method info
        null_balanced_accuracy_res["Comparison_Group"] = comparison_to_control_group
        null_balanced_accuracy_res["Scaling_Type"] = scaling_type
        
        # Save results
        null_balanced_accuracy_res.to_feather(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_Pairwise_{pairwise_feature_set}_{scaling_type}_scaler_SVM_null_balanced_accuracy_distributions.feather")
        