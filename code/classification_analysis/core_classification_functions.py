
import pandas as pd
from sklearn import svm
from sklearn.pipeline import Pipeline
import os.path
from sklearn.model_selection import RepeatedStratifiedKFold, cross_validate, permutation_test_score, StratifiedKFold, GridSearchCV, RandomizedSearchCV
import numpy as np
from sklearn.ensemble import RandomForestClassifier,GradientBoostingClassifier
from mixed_sigmoid_normalisation import MixedSigmoidScaler
from sklearn.metrics import balanced_accuracy_score,make_scorer
from copy import deepcopy

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


#######################################################################################################
# Robustness analysis functions
def fit_nested_CV(X, y, inner_cv, main_cv, pipe, num_folds=10, num_repeats=10, 
                  scoring="balanced_accuracy", num_jobs=2, 
                  param_grid={"model__C": [0.0001, 0.001, 0.01, 0.1, 1, 10, 100],
                              "model__class_weight": [None, "balanced"]}):

    # Define train/test splits from the passed cv argument
    repeat_folds = list(main_cv.split(X, y))
    training_folds, testing_folds = zip(*repeat_folds)

    # Reshape folds into a dataframe
    fold_repeat_structure = pd.DataFrame({"Fold_Number": np.arange(num_folds).tolist()*num_repeats,
                         "Repeat_Number": [i for i in range(num_repeats) for _ in range(num_folds)],
                         "Train_Folds": training_folds,
                         "Test_Folds": testing_folds})
    
    # Start a dataframe to compile the results
    fold_wise_results_list = []

    # Iterate over each row of the dataframe
    for index, row in fold_repeat_structure.iterrows():
        fold_number = row["Fold_Number"]
        repeat_number = row["Repeat_Number"]
        train_folds = row["Train_Folds"]
        test_folds = row["Test_Folds"]
        
        # Subset the training and test data
        training_data = X[train_folds,:]
        testing_data = X[test_folds,:]
        
        training_labels = np.array(y)[train_folds]
        testing_labels = np.array(y)[test_folds]
        
        # Define an inner 10-fold cross-validation for the grid search
        search = RandomizedSearchCV(pipe, 
                                    param_grid, 
                                    scoring=scoring, 
                                    cv=inner_cv, 
                                    refit=True, 
                                    n_jobs=num_jobs)
        
        # Execute the search
        search_result = search.fit(training_data, training_labels)
        
        # Extract the best-performing model fit on the whole training set
        best_model = search_result.best_estimator_
        
        # Predict the class values with best_model
        predicted_labels = best_model.predict(testing_data)
        
        # Compute the balanced accuracy
        fold_balacc = balanced_accuracy_score(testing_labels, predicted_labels)
    
        # Append to list
        fold_wise_results_list.append(fold_balacc)

    nested_CV_balacc_mean = np.mean(fold_wise_results_list)

    return nested_CV_balacc_mean

def robustness_analysis(X, y, model_dict, inner_cv, main_cv, num_folds=10, num_repeats=10, base_model_name="Linear_SVM_sklearn", scoring="balanced_accuracy", num_jobs=2):
    base_model = model_dict[base_model_name]

    base_pipe = Pipeline([('scaler', MixedSigmoidScaler(unit_variance=True)), 
                          ('model', base_model)])
    
    #######################################################################################################
    # Compare train vs. test balanced accuracy

    # Fit classifier
    training_classification_res = cross_validate(deepcopy(base_pipe), X, y,
                                                     cv=main_cv, scoring=scoring, n_jobs=num_jobs, 
                                                     return_estimator=False, return_train_score=True)
    test_balacc = training_classification_res["test_score"].mean()
    train_balacc = training_classification_res["train_score"].mean()

    # Save results to list
    training_balacc_df = pd.DataFrame({"Classifier_Type": base_model_name,
                                        "Train_Balanced_Accuracy": test_balacc,
                                        "Test_Balanced_Accuracy": train_balacc}, index=[0])

    #######################################################################################################
    # Compare across classifier types
    this_classifier_type_df_list = []
    for classifier_type in model_dict.keys():
        model = model_dict[classifier_type]
        
        pipe = Pipeline([('scaler', MixedSigmoidScaler(unit_variance=True)), 
                          ('model', model)])
        
        # Fit classifier
        this_classifier_res = cross_validate(pipe, X, y, cv=main_cv, scoring=scoring, n_jobs=num_jobs, 
                                            return_estimator=False, return_train_score=False)
        this_classifier_test_balacc = this_classifier_res["test_score"].mean()

        # Save results to list
        this_classifier_type_df = pd.DataFrame({"Classifier_Type": classifier_type,
                                            "Balanced_Accuracy": this_classifier_test_balacc}, index=[0])
        this_classifier_type_df_list.append(this_classifier_type_df)

    classifier_type_df = pd.concat(this_classifier_type_df_list, axis=0).reset_index(drop=True)

    #######################################################################################################
    
    # Compare regular vs. nested cross-validation for hyperparameter tuning
    inner_cv = StratifiedKFold(n_splits=num_folds, shuffle=True, random_state=127)

    nested_CV_res = fit_nested_CV(X, y, inner_cv, main_cv, deepcopy(base_pipe), num_folds=num_folds, num_repeats=num_repeats, scoring=scoring, num_jobs=num_jobs)
    nested_CV_df = pd.DataFrame({"Classifier_Type": f"{base_model_name}_nested_CV",
                                 "Balanced_Accuracy": nested_CV_res}, index=[0])
    
    # Return the three sets of results as a tuple
    return (training_balacc_df, classifier_type_df, nested_CV_df)