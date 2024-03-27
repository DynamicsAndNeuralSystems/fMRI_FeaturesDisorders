
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

    main_res_by_fold["Fold"] = main_res_by_fold.index % num_folds
    main_res_by_fold["Repeat"] = main_res_by_fold.index // num_repeats

    return (main_res_by_fold, splits_df, null_res)

def combine_main_results(files_to_merge, dataset_ID, disorder, average_across_folds=True):
    disorder_main_res_list = []
        
    # Iterate over the files in data_path
    for feather_file in files_to_merge:
        this_main_res = pd.read_feather(feather_file)
        this_main_res["Study"] = dataset_ID
        this_main_res["Disorder"] = disorder
        disorder_main_res_list.append(this_main_res)

    disorder_main_res = pd.concat(disorder_main_res_list, axis=0).reset_index(drop=True)

    # Average across folds/repeats if requested
    if average_across_folds:
        disorder_main_res = (disorder_main_res
                             .groupby(["group_var", "Analysis_Type"], as_index=False)['Balanced_Accuracy']
                             .agg(["mean", "std"])
                             .reset_index()
                             .rename(columns={"mean": "Balanced_Accuracy", "std": "Balanced_Accuracy_SD"}))
        
    # Return the nulls
    return disorder_main_res


def combine_null_results(files_to_merge, dataset_ID, disorder, num_null_iters=1000):
    disorder_nulls_list = []
        
    # Iterate over the files in data_path
    for feather_file in files_to_merge:
        this_null_res = pd.read_feather(feather_file)
        this_null_res["Null_Iteration"] = np.arange(1, num_null_iters+1)
        this_null_res["Study"] = dataset_ID
        this_null_res["Disorder"] = disorder
        disorder_nulls_list.append(this_null_res)

    disorder_nulls = pd.concat(disorder_nulls_list, axis=0).reset_index(drop=True)

    # Return the nulls
    return disorder_nulls
