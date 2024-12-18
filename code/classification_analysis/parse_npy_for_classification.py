from sklearnex import patch_sklearn
patch_sklearn()

# Load modules
import pandas as pd
import argparse
import os
import sys
import numpy as np

# add path to classification analysis functions
sys.path.insert(0, './')
from core_classification_functions import *
current_path = os.getcwd()
from mixed_sigmoid_normalisation import MixedSigmoidScaler

# Command-line arguments to parse
parser = argparse.ArgumentParser(description='Process inputs for pairwise data preparation.')
parser.add_argument('--dataset_ID', dest='dataset_ID')
parser.add_argument('--disorder', dest='disorder')
parser.add_argument('--input_data_file', dest='input_data_file')
parser.add_argument('--sample_IDs_file', dest='sample_IDs_file')
parser.add_argument('--class_labels_file', dest='class_labels_file')
parser.add_argument('--classifier_type', dest='classifier_type', default='Linear_SVM_sklearn')
parser.add_argument('--output_data_path', dest='output_data_path')
parser.add_argument('--feature_name', dest='feature_name')
parser.add_argument('--num_folds', dest='num_folds', default=10)
parser.add_argument('--num_repeats', dest='num_repeats', default=10)
parser.add_argument('--num_jobs', dest='num_jobs', default=4)
parser.add_argument('--num_null_iters', dest='num_null_iters', default=1)

# Parse input arguments
args = parser.parse_args()
dataset_ID = args.dataset_ID
disorder = args.disorder
input_data_file = args.input_data_file
sample_IDs_file = args.sample_IDs_file
class_labels_file = args.class_labels_file
classifier_type = args.classifier_type
output_data_path = args.output_data_path
classifier_type = args.classifier_type
num_folds = int(args.num_folds)
num_repeats = int(args.num_repeats)
num_jobs = int(args.num_jobs)
num_null_iters = int(args.num_null_iters)

# Read in input data
feature_data = np.load(input_data_file)

# Define analysis type
if "ROI" in input_data_file:
    Analysis_Type = "Brain_Region"
elif "combo_catch25_features_all_regions" in input_data_file:
    Analysis_Type = "Univariate_Combo"
elif "combined_univariate_catch25_and_pyspi14" in input_data_file:
    Analysis_Type = "SPI_Combo"
elif "catch25_feature" in input_data_file:
    Analysis_Type = "catch25_feature"
else:
    Analysis_Type = "pyspi14_SPI"

# Find grouping_var
model_name = input_data_file.split(".")[0]

if Analysis_Type=="Brain_Region":
    grouping_var = model_name.split("_ROI_")[1]
elif Analysis_Type=="Univariate_Combo":
    grouping_var = "Combo"
elif Analysis_Type == "SPI_Combo":
    grouping_var = model_name.split("combined_univariate_catch25_and_pyspi14_SPI_")[1]
elif Analysis_Type == "catch25_feature":
    grouping_var = model_name.split("_catch25_feature_")[1]
else:
    grouping_var = model_name.split("_pyspi14_SPI_")[1]

# Read in class labels
class_labels = np.load(class_labels_file, allow_pickle=True).tolist()
sample_IDs = np.load(sample_IDs_file, allow_pickle=True).tolist()
    
# Define main output data file for this feature
main_output_file_base = f"{dataset_ID}_{disorder}_{Analysis_Type}_{grouping_var}_{classifier_type}_{num_repeats}_repeats_{num_folds}_folds_CV"

os.makedirs(output_data_path, exist_ok=True)
os.makedirs(f"{output_data_path}/balanced_accuracy/{dataset_ID}_{disorder}", exist_ok=True)
os.makedirs(f"{output_data_path}/null_results/{dataset_ID}_{disorder}", exist_ok=True)
os.makedirs(f"{output_data_path}/robustness_analysis/{dataset_ID}_{disorder}", exist_ok=True)


# Define classification pipelines
if classifier_type=="Linear_SVM_sklearn":
    model = svm.SVC(kernel="linear", C=1, class_weight="balanced")
    model_noreg = svm.SVC(kernel="linear")
elif classifier_type=="RBF_SVM_sklearn":
    model = svm.SVC(kernel="rbf", C=1, class_weight="balanced")
    model_noreg = svm.SVC(kernel="rbf")
elif classifier_type=="RandomForest":
    model = model_noreg = RandomForestClassifier(n_estimators=100, class_weight="balanced", n_jobs=1)
elif classifier_type=="GradientBoost":
    model = model_noreg = GradientBoostingClassifier(n_estimators=100, learning_rate=1.0, max_depth=1)

pipe = Pipeline([('scaler', MixedSigmoidScaler(unit_variance=True)),
                ('model', model)])
pipe_noreg = Pipeline([('scaler', MixedSigmoidScaler(unit_variance=True)),
                        ('model', model_noreg)])

############### Main analysis with default C-value ###############
    
# # Main analysis
# if not os.path.isfile(f"{output_data_path}/null_results/{dataset_ID}_{disorder}/{main_output_file_base}_nulls.feather"):
#     # Define RepeatedStratifiedKFold splitter
#     RepeatedStratifiedKFold_splitter = RepeatedStratifiedKFold(n_splits=num_folds, n_repeats=num_repeats, random_state=127)

#     # Define classification pipelines
#     if classifier_type=="Linear_SVM_sklearn":
#         model = svm.SVC(kernel="linear", C=1, class_weight="balanced")
#         model_noreg = svm.SVC(kernel="linear")
#     elif classifier_type=="RBF_SVM_sklearn":
#         model = svm.SVC(kernel="rbf", C=1, class_weight="balanced")
#         model_noreg = svm.SVC(kernel="rbf")
#     elif classifier_type=="RandomForest":
#         model = model_noreg = RandomForestClassifier(n_estimators=100, class_weight="balanced", n_jobs=1)
#     elif classifier_type=="GradientBoost":
#         model = model_noreg = GradientBoostingClassifier(n_estimators=100, learning_rate=1.0, max_depth=1)

#     pipe = Pipeline([('scaler', MixedSigmoidScaler(unit_variance=True)),
#                     ('model', model)])
#     pipe_noreg = Pipeline([('scaler', MixedSigmoidScaler(unit_variance=True)),
#                             ('model', model_noreg)])
        
#     # Define scorers
#     scorers = [make_scorer(balanced_accuracy_score)]
#     scoring_names = ["Balanced_Accuracy"]

#     print(f"Now running {classifier_type} for {disorder} {grouping_var}")
#     main_classification_res, splits_df, null_classification_res = run_k_fold_classifier_for_feature(feature_data = feature_data, 
#                                                                                         pipe = pipe,
#                                                                                         CV_splitter = RepeatedStratifiedKFold_splitter,
#                                                                                         class_labels=class_labels,
#                                                                                         sample_IDs = sample_IDs,
#                                                                                         scorers=scorers,
#                                                                                         scoring_names=scoring_names,
#                                                                                         num_null_iters=num_null_iters,
#                                                                                         num_folds = num_folds,
#                                                                                         num_repeats = num_repeats,
#                                                                                         num_jobs = num_jobs)

#     # Assign key details to dataframes
#     main_classification_res["Disorder"] = disorder
#     main_classification_res["Study"] = dataset_ID
#     main_classification_res["Analysis_Type"] = Analysis_Type
#     main_classification_res["group_var"] = grouping_var
#     main_classification_res["Classifier_Type"] = classifier_type

#     # Save results to feather file
#     main_classification_res.to_feather(f"{output_data_path}/balanced_accuracy/{dataset_ID}_{disorder}/{main_output_file_base}_main_res.feather")

#     # If nulls were requested, save those too
#     if num_null_iters > 0:
#         null_classification_res["Analysis_Type"] = Analysis_Type
#         null_classification_res["Disorder"] = disorder
#         null_classification_res["Study"] = dataset_ID
#         null_classification_res["Classifier_Type"] = classifier_type
#         null_classification_res["group_var"] = grouping_var
#         null_classification_res.to_feather(f"{output_data_path}/null_results/{dataset_ID}_{disorder}/{main_output_file_base}_nulls.feather")


# # Robustness analysis
# robustness_output_file_base = f"{dataset_ID}_{disorder}_{Analysis_Type}_{grouping_var}"
# print(f"Robustness output file: {robustness_output_file_base}")
# if not os.path.isfile(f"{output_data_path}/robustness_analysis/{dataset_ID}_{disorder}/{robustness_output_file_base}_training_balacc_df.feather"):
#     print(f"Now running robustness analysis for {disorder} {grouping_var}")

#     # Define CV splitters
#     main_cv = RepeatedStratifiedKFold(n_splits=num_folds, n_repeats=num_repeats, random_state=127)
#     inner_cv = StratifiedKFold(n_splits=num_folds, shuffle=True, random_state=127)

#     # Define the model dict
#     model_dict = {"Linear_SVM_sklearn": svm.SVC(kernel="linear", class_weight="balanced", C=1),
#                    "Linear_SVM_libsvm": svm.LinearSVC(C=1, dual=False, penalty='l1', class_weight='balanced'),
#                    "RBF_SVM_sklearn": svm.SVC(kernel="rbf", class_weight="balanced", C=1),
#                    "RandomForest": RandomForestClassifier(n_estimators=100, class_weight="balanced"),
#                    "GradientBoost": GradientBoostingClassifier(n_estimators=100, learning_rate=1.0, max_depth=1)}
    
#     # Define class weight and C parameter grid for tuning
#     param_grid={"model__C": [0.0001, 0.001, 0.01, 0.1, 1, 10, 100],
#                 "model__class_weight": [None, "balanced"]}

#     # Run the robustness analysis
#     (training_balacc_df, classifier_type_df, nested_CV_df) = robustness_analysis(feature_data, class_labels, model_dict, inner_cv, main_cv, 
#                                                                                  num_folds=num_folds, num_repeats=num_repeats, 
#                                                                                  base_model_name="Linear_SVM_sklearn", 
#                                                                                  scoring="balanced_accuracy", num_jobs=num_jobs)

#     # Assign Analysis_Type, group_var, Disorder, and Dataset columns
#     for df in [training_balacc_df, classifier_type_df, nested_CV_df]:
#         df["Analysis_Type"] = Analysis_Type
#         df["group_var"] = grouping_var
#         df["Disorder"] = disorder
#         df["Dataset"] = dataset_ID

#     # Save results to feather files
#     training_balacc_df.to_feather(f"{output_data_path}/robustness_analysis/{dataset_ID}_{disorder}/{robustness_output_file_base}_training_balacc_df.feather")
#     classifier_type_df.to_feather(f"{output_data_path}/robustness_analysis/{dataset_ID}_{disorder}/{robustness_output_file_base}_classifier_type_df.feather")
#     nested_CV_df.to_feather(f"{output_data_path}/robustness_analysis/{dataset_ID}_{disorder}/{robustness_output_file_base}_nested_CV_df.feather")



# First 10 principal components analysis
first10_PC_output_file_base = f"{dataset_ID}_{disorder}_{Analysis_Type}_{grouping_var}"
print(f"First 10 PCs output file: {first10_PC_output_file_base}")
if not os.path.isfile(f"{output_data_path}/robustness_analysis/{dataset_ID}_{disorder}/{first10_PC_output_file_base}_first10_PCs_main_classification_res.feather"):
    
    # Define RepeatedStratifiedKFold splitter
    RepeatedStratifiedKFold_splitter = RepeatedStratifiedKFold(n_splits=num_folds, n_repeats=num_repeats, random_state=127)

    # Run the robustness analysis
    X_10_PCs_classification_res_df = first10_PC_analysis(X=feature_data, 
                                                         y=class_labels, 
                                                         pipe=pipe, 
                                                         CV_splitter=RepeatedStratifiedKFold_splitter,
                                                         scoring="balanced_accuracy",
                                                         num_folds = num_folds,
                                                         num_repeats = num_repeats,
                                                         num_jobs = num_jobs)

    # Assign Analysis_Type, group_var, Disorder, and Dataset columns
    X_10_PCs_classification_res_df = (X_10_PCs_classification_res_df.assign(Analysis_Type=Analysis_Type,
                                                                            group_var=grouping_var,
                                                                            Disorder=disorder,
                                                                            Dataset=dataset_ID))
    
    # Save results to feather files
    X_10_PCs_classification_res_df.to_feather(f"{output_data_path}/robustness_analysis/{dataset_ID}_{disorder}/{first10_PC_output_file_base}_first10_PCs_main_classification_res.feather")