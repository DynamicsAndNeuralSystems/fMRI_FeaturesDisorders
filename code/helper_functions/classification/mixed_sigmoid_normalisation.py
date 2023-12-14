from sklearn.utils.validation import (
    check_is_fitted,
    check_random_state,
    _check_sample_weight,
    FLOAT_DTYPES,
)
from sklearn.preprocessing import StandardScaler, RobustScaler
import pandas as pd
import math
from sklearn.base import BaseEstimator, TransformerMixin
import numpy as np

class MixedSigmoidScaler(TransformerMixin, BaseEstimator):
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
        # self._validate_params()

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
        self.mean_ = np.nanmean(X, axis=0)
        self.sd_ = np.nanstd(X, axis=0)
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

        # Iterate over each feature column
        for i in range(X_transformed.shape[1]):
            # catch for if the IQR difference is zero, in which case we apply standard logistic transformation using mean and SD
            if self.IQR_diffs_[i] == 0:
                scale = 1
                X_i = 1/(1 + math.e**(- ( (X[:,i] - self.mean_[i]) / self.sd_[i] ) ))
            # otherwise, set the scaling factor to the IQR/1.35
            else: 
                IQR = self.IQR_diffs_[i]
                # Center by removing the median, standarize by dividing by scale value
                X_i = 1/(1 + math.e**(- ( (X[:,i] - self.median_[i]) / IQR/1.35 ) ))
            
            # Scale to unit variance [0,1] if specified
            if self.unit_variance:
                # Find the min value in this feature
                X_i_min = np.nanmin(X_i)
                # Find the range of values in this feature
                X_i_diff = np.ptp(X_i)

                # If the range is 0, return X_i_diff since the scaling above would already set all values to 0.5
                if X_i_diff == 0:
                    X_transformed[:,i] = X_i
                # Otherwise, subtract the minimum and divide by the range
                else:
                    X_transformed[:,i] = (X_i - X_i_min) / X_i_diff
            else: 
                X_transformed[:,i] = X_i

        return(X_transformed)

