
import pandas as pd
from sklearn import svm
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, RobustScaler
import os.path
from sklearn.model_selection import StratifiedKFold, cross_val_predict, cross_validate, permutation_test_score
from sklearn.base import BaseEstimator, TransformerMixin
import numpy as np
from sklearn.utils.validation import (
    check_is_fitted,
    check_random_state,
    _check_sample_weight,
    FLOAT_DTYPES,
)
import math

class RobustSigmoidScaler(TransformerMixin, BaseEstimator):
    """Scale features using statistics that are robust to outliers.
    This Scaler removes the median and scales the data according to
    the quantile range (defaults to IQR: Interquartile Range).
    The IQR is the range between the 1st quartile (25th quantile)
    and the 3rd quartile (75th quantile).
    Centering and scaling happen independently on each feature by
    computing the relevant statistics on the samples in the training
    set. Median and interquartile range are then stored to be used on
    later data using the :meth:`transform` method.
    Standardization of a dataset is a common requirement for many
    machine learning estimators. Typically this is done by removing the mean
    and scaling to unit variance. However, outliers can often influence the
    sample mean / variance in a negative way. In such cases, the median and
    the interquartile range often give better results.
    .. versionadded:: 0.17
    Read more in the :ref:`User Guide <preprocessing_scaler>`.
    Parameters
    ----------
    with_centering : bool, default=True
        If `True`, center the data before scaling.
        This will cause :meth:`transform` to raise an exception when attempted
        on sparse matrices, because centering them entails building a dense
        matrix which in common use cases is likely to be too large to fit in
        memory.
    with_scaling : bool, default=True
        If `True`, scale the data to interquartile range.
    quantile_range : tuple (q_min, q_max), 0.0 < q_min < q_max < 100.0, \
        default=(25.0, 75.0)
        Quantile range used to calculate `scale_`. By default this is equal to
        the IQR, i.e., `q_min` is the first quantile and `q_max` is the third
        quantile.
        .. versionadded:: 0.18
    copy : bool, default=True
        If `False`, try to avoid a copy and do inplace scaling instead.
        This is not guaranteed to always work inplace; e.g. if the data is
        not a NumPy array or scipy.sparse CSR matrix, a copy may still be
        returned.
    unit_variance : bool, default=False
        If `True`, scale data so that normally distributed features have a
        variance of 1. In general, if the difference between the x-values of
        `q_max` and `q_min` for a standard normal distribution is greater
        than 1, the dataset will be scaled down. If less than 1, the dataset
        will be scaled up.
        .. versionadded:: 0.24
    Attributes
    ----------
    center_ : array of floats
        The median value for each feature in the training set.
    scale_ : array of floats
        The (scaled) interquartile range for each feature in the training set.
        .. versionadded:: 0.17
           *scale_* attribute.
    n_features_in_ : int
        Number of features seen during :term:`fit`.
        .. versionadded:: 0.24
    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.
        .. versionadded:: 1.0
    See Also
    --------
    robust_scale : Equivalent function without the estimator API.
    sklearn.decomposition.PCA : Further removes the linear correlation across
        features with 'whiten=True'.
    Notes
    -----
    For a comparison of the different scalers, transformers, and normalizers,
    see :ref:`examples/preprocessing/plot_all_scaling.py
    <sphx_glr_auto_examples_preprocessing_plot_all_scaling.py>`.
    https://en.wikipedia.org/wiki/Median
    https://en.wikipedia.org/wiki/Interquartile_range
    Examples
    --------
    >>> from sklearn.preprocessing import RobustScaler
    >>> X = [[ 1., -2.,  2.],
    ...      [ -2.,  1.,  3.],
    ...      [ 4.,  1., -2.]]
    >>> transformer = RobustScaler().fit(X)
    >>> transformer
    RobustScaler()
    >>> transformer.transform(X)
    array([[ 0. , -2. ,  0. ],
           [-1. ,  0. ,  0.4],
           [ 1. ,  0. , -1.6]])
    """

    _parameter_constraints: dict = {
         "quantile_range": [tuple],
         "copy": ["boolean"],
        "unit_variance": ["boolean"]
    }

    def __init__(
        self,
        *,
        quantile_range=(25.0, 75.0),
        copy = True,
        unit_variance=False
    ):
        self.quantile_range = quantile_range
        self.copy = copy
        self.unit_variance = unit_variance

    def fit(self, X, y=None):
        """Compute the median and quantiles to be used for scaling.
        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The data used to compute the median and quantiles
            used for later scaling along the features axis.
        y : Ignored
            Not used, present here for API consistency by convention.
        Returns
        -------
        self : object
            Fitted scaler.
        """
        self._validate_params()

        # at fit, convert sparse matrices to csc for optimized computation of
        # the quantiles
        X = self._validate_data(
            X,
            accept_sparse="csc",
            dtype=FLOAT_DTYPES,
            force_all_finite="allow-nan",
        )

        q_min, q_max = self.quantile_range
        if not 0 <= q_min <= q_max <= 100:
            raise ValueError("Invalid quantile range: %s" % str(self.quantile_range))
            
        self.IQR_ = np.nanpercentile(X, self.quantile_range, axis=0)
        self.IQR_diffs_ = self.IQR_[1,:] - self.IQR_[0,:]
        self.median_ = np.nanmedian(X, axis=0)

        return self


    def transform(self, X):
        """Center and scale the data.
        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The data used to scale along the specified axis.
        Returns
        -------
        X_tr : {ndarray, sparse matrix} of shape (n_samples, n_features)
            Transformed array.
        """
        check_is_fitted(self)
        X = self._validate_data(
            X,
            accept_sparse=("csr", "csc"),
            copy=self.copy,
            dtype=FLOAT_DTYPES,
            reset=False,
            force_all_finite="allow-nan",
        )

        X_transformed = np.zeros(X.shape)

        for i in range(X_transformed.shape[1]):
            X_i = 1/(1 + math.e**(- ( (X[:,i] - self.median_[i]) / (self.IQR_diffs_[i] / 1.35) ) ))
            
            if self.unit_variance:
                X_i_min = np.nanmin(X_i)
                X_i_max = np.nanmax(X_i)
                X_i_diff = np.ptp(X_i)
                X_transformed[:,i] = (X_i - X_i_min) / X_i_diff
            else: 
                X_transformed[:,i] = X_i

        return(X_transformed)


def run_k_fold_SVM_for_feature(feature_data, 
                               feature_list,
                               grouping_var_name,
                               scoring_method,
                               sample_and_class_df,
                               scaling_type,
                               class_labels,
                               num_null_iters,
                               num_jobs = 8):
        
    # Define the pipeline
    if scaling_type == "robust":
        pipe = Pipeline([('scaler', RobustScaler()), 
                         ('SVM', svm.SVC(kernel="linear", C=1, shrinking=False, 
                                     class_weight="balanced"))])
    elif scaling_type == "robustsigmoid":
        pipe = Pipeline([('scaler', RobustSigmoidScaler()), 
                         ('SVM', svm.SVC(kernel = "linear", C = 1, shrinking = False, 
                         class_weight = "balanced"))])
    else: 
        pipe = Pipeline([('scaler', StandardScaler()), 
                         ('SVM', svm.SVC(kernel="linear", C=1, shrinking=False, 
                                     class_weight="balanced"))])
    
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
        # Initialise sample and class dataframe for repeat
        sample_and_class_df_for_repeat = sample_and_class_df.copy(deep=True)
        skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=i)
        
        # Find splits
        splits = list(skf.split(feature_data, class_labels))
        fold_split_list = []
        
        # Fit SVM with 10-fold cross-validation
        cv_results = cross_validate(pipe, feature_data, class_labels, 
                                    cv=skf, scoring=scoring_method,
                                    n_jobs = int(num_jobs),
                                    return_estimator=True)
    
        
        # Extract fold number and feature coefficients
        coef_list = []
        
        for f in range(10):
            # Split for fold number
            test_indices = splits[f][1]
            test_indices_df = pd.DataFrame(test_indices, columns=["Sample_Index"])
            test_indices_df["Sample_ID"] = [sample_and_class_df_for_repeat.Sample_ID.tolist()[index] for index in test_indices]
            test_indices_df["Fold"] = f+1
            test_indices_df["Repeat"] = i+1        
            fold_split_list.append(test_indices_df)
            
            # Coefficients
            SVM_for_fold = cv_results["estimator"][f]["SVM"]
            SVM_coefs = pd.DataFrame(SVM_for_fold.coef_.T)
            SVM_coefs.reset_index(inplace=True)
            SVM_coefs = SVM_coefs.rename(columns = {'index':'Feature_Number',
                                                          0: "Coefficient"})
            SVM_coefs["Feature Name"] = feature_list
            SVM_coefs["Fold"] = f+1
            SVM_coefs["Repeat Number"] = i+1
            coef_list.append(SVM_coefs)
        
        # Combine lists into dataframes
        fold_splits = pd.concat(fold_split_list)
        # Add sample ID name based on index

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
        CV_pred = cross_val_predict(pipe, feature_data, class_labels, cv=skf)
        sample_and_class_df_for_repeat["CV_Predicted_Diagnosis"] = CV_pred
        sample_and_class_df_for_repeat["Repeat_Number"] = i+1
        CV_sample_predictions_list.append(sample_and_class_df_for_repeat)
        
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
    null_model_balanced_accuracy_list = permutation_test_score(pipe,
                                            feature_data,
                                            class_labels,
                                            cv=10,
                                            random_state=27,
                                            n_jobs = int(num_jobs),
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
                       scaling_type="robust",
                       num_jobs = 8,
                       overwrite=False):

    # Check if file already exists or overwrite flag is set
    if not os.path.isfile(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_{scaling_type}_scaler_SVM_balanced_accuracy.feather") or overwrite:
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
                                        scaling_type = scaling_type,
                                        num_null_iters = num_null_iters,
                                        num_jobs = num_jobs)
            
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
                                        scaling_type = scaling_type,
                                        num_null_iters = num_null_iters,
                                        num_jobs = num_jobs)
            
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
                                    scaling_type = scaling_type,
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
        fold_assignments_res.to_feather(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_{scaling_type}_scaler_SVM_fold_assignments.feather")
        SVM_coefficients_res.to_feather(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_{scaling_type}_scaler_SVM_fold_SVM_coefficients.feather")
        balanced_accuracy_res.to_feather(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_{scaling_type}_scaler_SVM_balanced_accuracy.feather")
        CV_sample_predictions_res.to_feather(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_{scaling_type}_scaler_SVM_sample_predictions.feather")
        null_balanced_accuracy_res.to_feather(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_{scaling_type}_scaler_SVM_null_balanced_accuracy_distributions.feather")
   
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
                       scaling_type="robust",
                       num_jobs = 8,
                       overwrite=False):
    

    # Check if file already exists or overwrite flag is set
    if not os.path.isfile(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Pairwise_{pairwise_feature_set}_{scaling_type}_scaler_SVM_balanced_accuracy.feather") or overwrite:

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
                                        scaling_type = scaling_type,
                                        num_null_iters = num_null_iters,
                                        num_jobs = num_jobs)

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
        fold_assignments_res.to_feather(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Pairwise_{pairwise_feature_set}_{scaling_type}_scaler_SVM_fold_assignments.feather")
        SVM_coefficients_res.to_feather(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Pairwise_{pairwise_feature_set}_{scaling_type}_scaler_SVM_fold_SVM_coefficients.feather")
        balanced_accuracy_res.to_feather(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Pairwise_{pairwise_feature_set}_{scaling_type}_scaler_SVM_balanced_accuracy.feather")
        CV_sample_predictions_res.to_feather(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Pairwise_{pairwise_feature_set}_{scaling_type}_scaler_SVM_sample_predictions.feather")
        null_balanced_accuracy_res.to_feather(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Pairwise_{pairwise_feature_set}_{scaling_type}_scaler_SVM_null_balanced_accuracy_distributions.feather")
    
    
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
                       num_null_iters,
                       scaling_type = "robust",
                       num_jobs = 8,
                       overwrite=False):

    # Check if file already exists or overwrite flag is set
    if not os.path.isfile(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_Pairwise_{pairwise_feature_set}_{scaling_type}_scaler_SVM_balanced_accuracy.feather") or overwrite:
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
                                        scaling_type = scaling_type,
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
        fold_assignments_res.to_feather(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_Pairwise_{pairwise_feature_set}_{scaling_type}_scaler_SVM_fold_assignments.feather")
        SVM_coefficients_res.to_feather(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_Pairwise_{pairwise_feature_set}_{scaling_type}_scaler_SVM_fold_SVM_coefficients.feather")
        balanced_accuracy_res.to_feather(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_Pairwise_{pairwise_feature_set}_{scaling_type}_scaler_SVM_balanced_accuracy.feather")
        CV_sample_predictions_res.to_feather(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_Pairwise_{pairwise_feature_set}_{scaling_type}_scaler_SVM_sample_predictions.feather")
        null_balanced_accuracy_res.to_feather(f"{pydata_path}/{dataset_ID}_{comparison_to_control_group}_Univariate_{univariate_feature_set}_Pairwise_{pairwise_feature_set}_{scaling_type}_scaler_SVM_null_balanced_accuracy_distributions.feather")
      