
import pandas as pd
from sklearn import svm
from sklearn.pipeline import Pipeline
import os.path
from sklearn.model_selection import RepeatedStratifiedKFold, cross_val_predict, cross_validate, permutation_test_score
import numpy as np
from sklearn.ensemble import RandomForestClassifier,GradientBoostingClassifier
from mixed_sigmoid_normalisation import MixedSigmoidScaler
from sklearn.metrics import balanced_accuracy_score,make_scorer

def run_k_fold_classifier_for_feature(feature_data, 
                                pipe,
                                CV_splitter,
                               grouping_var_name,
                               analysis_type, 
                               sample_IDs,
                               class_labels,
                               scorers,
                               scoring_names,
                               random_num_for_perm = 127,
                               num_null_iters = 1000,
                               num_folds = 10,
                               num_repeats = 10,
                               num_jobs = 10):
    
    # Fit the main and null results
    main_res_by_fold, null_res = permutation_test_score(pipe,
                                                feature_data, 
                                                class_labels,
                                                cv=CV_splitter,
                                                scoring=scorers,
                                                n_permutations=num_null_iters,
                                                n_jobs=num_jobs,
                                                random_state=random_num_for_perm,
                                                return_all_folds = True)
    
    # Unpack the main and null results from the dictionary
    if isinstance(main_res_by_fold, dict):
        main_res_by_fold = pd.DataFrame(main_res_by_fold)

        # Get the original column names dynamically
        original_column_names = main_res_by_fold.columns.tolist()

        # Rename the columns
        new_columns = dict(zip(original_column_names, scoring_names))
        main_res_by_fold = main_res_by_fold.rename(columns=new_columns)

    if num_null_iters > 0:
        if isinstance(null_res, dict):
            null_res = pd.DataFrame(null_res)
            # Check number of columns in case AUC results were split across jobs
            if len(null_res.columns) > 2:
                null_res['coalesce'] = null_res.bfill(axis=1).iloc[:, 0]
                null_res = pd.DataFrame({"Balanced_Accuracy": null_res.coalesce})

            # Assign null iteration number
            null_res["Null_Iteration"] = range(1, num_null_iters+1)
        
        # Confirm string column names
        null_res.columns = ["Balanced_Accuracy", "Null_Iteration"]

    else:
        null_res = pd.DataFrame(columns=["Balanced_Accuracy", "Null_Iteration"])

    # Find splits
    splits = list(CV_splitter.split(feature_data, class_labels))

    # Convert splits to a dataframe
    splits_df = pd.DataFrame(splits, columns = ["Train", "Test"])

    # Map the indices in splits_df to the corresponding sample_IDs
    splits_df["Train"] = splits_df["Train"].map(lambda x: [sample_IDs[i] for i in x])
    splits_df["Test"] = splits_df["Test"].map(lambda x: [sample_IDs[i] for i in x])

    # Assign the fold and repeat numbers
    splits_df["Fold"] = splits_df.index % num_folds
    splits_df["Repeat"] = splits_df.index // num_repeats
    splits_df["group_var"] = grouping_var_name
    splits_df["Analysis_Type"] = analysis_type

    main_res_by_fold["Fold"] = main_res_by_fold.index % num_folds
    main_res_by_fold["Repeat"] = main_res_by_fold.index // num_repeats
    main_res_by_fold["group_var"] = grouping_var_name
    main_res_by_fold["Analysis_Type"] = analysis_type

    if num_null_iters > 0:
        return (main_res_by_fold, splits_df, null_res)
    else:
        return (main_res_by_fold, splits_df)

def combine_null_results(data_path, dataset_ID, disorder):
    disorder_path = f"{data_path}/{dataset_ID}_{disorder}"
    disorder_nulls_list = []

    # Iterate over the files in data_path
    for feather_file in os.listdir(disorder_path):
        this_null_res = pd.read_feather(f"{disorder_path}/{feather_file}")
        starting_point = int(feather_file.split("from_")[1].replace(".feather", "")) - 1
        # Fix the null_iteration number by adding starting_point to the values
        this_null_res["Null_Iteration"] = this_null_res["Null_Iteration"] + starting_point
        disorder_nulls_list.append(this_null_res)

    disorder_nulls = pd.concat(disorder_nulls_list, axis=0).reset_index(drop=True)

    # Fix the analysis type
    def fix_analysis(analysis):
        if "combo_catch25" in analysis:
            return 'Univariate_Combo'
        elif "ROI" in analysis:
            return 'Brain_Region'
        elif "catch25_feature" in analysis:
            return 'catch25_feature'
        else:
            return "Other"
        

    disorder_nulls = pd.concat(disorder_nulls_list, axis=0).reset_index(drop=True)
    disorder_nulls["Analysis_Type"] = disorder_nulls["Feature_Name"].map(lambda x: fix_analysis(x))

    # Return the nulls
    return disorder_nulls