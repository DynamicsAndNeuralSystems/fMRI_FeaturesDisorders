# https://scikit-learn.org/stable/auto_examples/compose/plot_digits_pipe.html

import numpy as np
import pandas as pd

from sklearn import datasets
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

# Define a pipeline to search for the best combination of PCA truncation
# and classifier regularization.
pca = PCA()
# Define a Standard Scaler to normalize inputs
scaler = StandardScaler()

# Define one pipeline for linear SVM and one for RBF
pipe_SVM_linear = Pipeline(steps=[("scaler", scaler), ("pca", pca), ("SVM", SVC(kernel="linear"))])
pipe_SVM_radial = Pipeline(steps=[("scaler", scaler), ("pca", pca), ("SVM", SVC(kernel="rbf"))])

X_digits, y_digits = datasets.load_digits(return_X_y=True)
# Parameters of pipelines can be set using '__' separated parameter names:
param_grid = {
    "pca__n_components": np.arange(5, 65, 5)
}

# Search through PC components in-sample
search_in_sample_linear = GridSearchCV(pipe_SVM_linear, param_grid, cv=[(slice(None), slice(None))], scoring="balanced_accuracy")
search_in_sample_linear.fit(X_digits, y_digits)
search_in_sample_radial = GridSearchCV(pipe_SVM_radial, param_grid, cv=[(slice(None), slice(None))], scoring="balanced_accuracy")
search_in_sample_radial.fit(X_digits, y_digits)


# Separately search through PC components out-of-sample
search_out_of_sample_linear = GridSearchCV(pipe_SVM_linear, param_grid, cv=10, scoring="balanced_accuracy")
search_out_of_sample_linear.fit(X_digits, y_digits)
search_out_of_sample_radial = GridSearchCV(pipe_SVM_radial, param_grid, cv=10, scoring="balanced_accuracy")
search_out_of_sample_radial.fit(X_digits, y_digits)


# For each number of components, find the best classifier results
in_sample_PC_performance_results_linear = (pd.DataFrame(search_in_sample_linear.cv_results_)[["param_pca__n_components",
                                                                               "mean_test_score"]]
                                    .assign(SVM_Type = "Linear")
                                    .rename(columns={'param_pca__n_components': 'Num_PCs',
                                                     'mean_test_score': 'In_Balanced_Accuracy'})
                                    )

in_sample_PC_performance_results_radial = (pd.DataFrame(search_in_sample_radial.cv_results_)[["param_pca__n_components",
                                                                               "mean_test_score"]]
                                    .assign(SVM_Type = "RBF")
                                    .rename(columns={'param_pca__n_components': 'Num_PCs',
                                                     'mean_test_score': 'In_Balanced_Accuracy'})
                                    )

# TODO: also add a 10x repeat outer loop to this and store the average balanced accuracy at each # of PCs
out_of_sample_PC_performance_results_linear = (pd.DataFrame(search_out_of_sample_linear.cv_results_)[["param_pca__n_components",
                                                                               "mean_test_score"]]
                                    .assign(SVM_Type = "Linear")
                                    .rename(columns={'param_pca__n_components': 'Num_PCs',
                                                     'mean_test_score': 'Out_Balanced_Accuracy'})
                                    )

out_of_sample_PC_performance_results_radial = (pd.DataFrame(search_out_of_sample_radial.cv_results_)[["param_pca__n_components",
                                                                               "mean_test_score"]]
                                    .assign(SVM_Type = "RBF")
                                    .rename(columns={'param_pca__n_components': 'Num_PCs',
                                                     'mean_test_score': 'Out_Balanced_Accuracy'})
                                    )
