import os
import numpy as np
import nibabel as nib
import pandas as pd
import sys
import argparse

# Command-line arguments to parse
parser = argparse.ArgumentParser()
parser.add_argument('--num_jobs', dest='num_jobs', default=4)
args = parser.parse_args()
num_jobs = int(args.num_jobs)

# add path to classification analysis functions
sys.path.insert(0, 'code/classification_analysis/')
from core_classification_functions import *
current_path = os.getcwd()

data_path = "/headnode1/abry4213/data/TS_feature_manuscript"

# Load participants included
UCLA_CNP_subjects_to_keep = pd.read_feather(f"{data_path}/time_series_features/UCLA_CNP_filtered_sample_info_catch25_pyspi14.feather")
ABIDE_subjects_to_keep = pd.read_feather(f"{data_path}/time_series_features/ABIDE_filtered_sample_info_catch25_pyspi14.feather")
merged_subjects_to_keep = pd.concat([UCLA_CNP_subjects_to_keep, ABIDE_subjects_to_keep], axis=0).Sample_ID.tolist()

# Load metadata
UCLA_CNP_metadata = (pd.read_feather(f"{data_path}/input_data/UCLA_CNP_sample_metadata.feather")
                        .assign(Study = "UCLA_CNP")
                        .query("Sample_ID in @UCLA_CNP_subjects_to_keep.Sample_ID"))
ABIDE_metadata = (pd.read_feather(f"{data_path}/input_data/ABIDE_sample_metadata.feather")
                        .assign(Study = "ABIDE")
                        .query("Sample_ID in @ABIDE_subjects_to_keep.Sample_ID"))

# Load head movement 
UCLA_CNP_head_mvmt = (pd.read_table(f"{data_path}/movement_data/UCLA_CNP_Mean_FD_Power.txt", sep=',')
                      .assign(Study = "UCLA_CNP")
                        .query("Sample_ID in @UCLA_CNP_subjects_to_keep.Sample_ID"))
ABIDE_head_mvmt = (pd.read_table(f"{data_path}/movement_data/ABIDE_Mean_FD_Power.txt", sep=',', dtype={'Sample_ID': str,
                                                                                              'Mean_FD_Power': float})
                   .assign(Study = "ABIDE")
                        .query("Sample_ID in @ABIDE_subjects_to_keep.Sample_ID"))

# Merge metadata + head movement
merged_metadata = pd.concat([UCLA_CNP_metadata, ABIDE_metadata], axis=0).merge(pd.concat([UCLA_CNP_head_mvmt, ABIDE_head_mvmt], axis=0))

# Study/disorder lookup table
study_disorder_lookup = {'SCZ': 'UCLA_CNP', 
                          'BP': 'UCLA_CNP', 
                          'ADHD': 'UCLA_CNP', 
                          'ASD': 'ABIDE'}

# We'll work with the top two sites: site 20 (N=106) and site 5 (N=98).
ABIDE_site20_indices = ABIDE_metadata.Sample_ID[ABIDE_metadata.Site=='20'].index.tolist()
ABIDE_site20 = ABIDE_metadata.query("Site=='20'")
ABIDE_site20_participants = ABIDE_site20.Sample_ID.tolist()
ABIDE_site20_class_labels = [int(i=='ASD') for i in ABIDE_site20["Diagnosis"].tolist()]

ABIDE_site5_indices = ABIDE_metadata.Sample_ID[ABIDE_metadata.Site=='5'].index.tolist()
ABIDE_site5 = ABIDE_metadata.query("Site=='5'")
ABIDE_site5_participants = ABIDE_site5.Sample_ID.tolist()
ABIDE_site5_class_labels = [int(i=='ASD') for i in ABIDE_site5["Diagnosis"].tolist()]

# Define parameters for classification
num_folds = 10
num_repeats = 10
num_null_iters = 0
RepeatedStratifiedKFold_splitter = RepeatedStratifiedKFold(n_splits=num_folds, n_repeats=num_repeats, random_state=127)
model = svm.SVC(kernel="linear", C=1, class_weight="balanced")
pipe = Pipeline([('scaler', MixedSigmoidScaler(unit_variance=True)),
                    ('model', model)])
classifier_type = "Linear_SVM_sklearn"

# Define scorers
scorers = [make_scorer(balanced_accuracy_score)]
scoring_names = ["Balanced_Accuracy"]

# Define all models to test
study = "ABIDE" 
disorder = "ASD" 

ASD_univariate_models = pd.read_table(f"/headnode1/abry4213/data/TS_feature_manuscript/time_series_features/processed_numpy_files/{study}_{disorder}_univariate_models.txt",header=None)
ASD_pairwise_models = pd.read_table(f"/headnode1/abry4213/data/TS_feature_manuscript/time_series_features/processed_numpy_files/{study}_{disorder}_pairwise_models.txt",header=None)
ASD_combined_univariate_pairwise_models  = pd.read_table(f"/headnode1/abry4213/data/TS_feature_manuscript/time_series_features/processed_numpy_files/{study}_{disorder}_combined_univariate_pairwise_models.txt",header=None)
ASD_all_models = pd.concat([ASD_univariate_models, ASD_pairwise_models, ASD_combined_univariate_pairwise_models])
ASD_all_models.columns = ["Model_Name"]

def call_models_for_ABIDE_site(ASD_all_models, site_number, site_indices, class_labels, site_sample_IDs):
    for model_name in ASD_all_models["Model_Name"].tolist(): 
        if not os.path.isfile(f"{data_path}/classification_results/robustness_analysis/ABIDE_ASD/{model_name}_Site{site_number}_main_classification_res.feather"):
            # Define analysis type
            if "ROI" in model_name:
                Analysis_Type = "Brain_Region"
            elif "combo_catch25_features_all_regions" in model_name:
                Analysis_Type = "Univariate_Combo"
            elif "combined_univariate_catch25_and_pyspi14" in model_name:
                Analysis_Type = "SPI_Combo"
            elif "catch25_feature" in model_name:
                Analysis_Type = "catch25_feature"
            else:
                Analysis_Type = "pyspi14_SPI"

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

            model_data = np.load(f"/headnode1/abry4213/data/TS_feature_manuscript/time_series_features/processed_numpy_files/{model_name}.npy")[site_indices,:]

            # Find balanced accuracy for dataset 
            try:
                main_classification_res, _, _ = run_k_fold_classifier_for_feature(feature_data = model_data, 
                                                                                                    pipe = pipe,
                                                                                                    CV_splitter = RepeatedStratifiedKFold_splitter,
                                                                                                    class_labels=class_labels,
                                                                                                    sample_IDs = site_sample_IDs,
                                                                                                    scorers=scorers,
                                                                                                    scoring_names=scoring_names,
                                                                                                    num_null_iters=num_null_iters,
                                                                                                    num_folds = num_folds,
                                                                                                    num_repeats = num_repeats,
                                                                                                    num_jobs = num_jobs)
                
                # Assign key details to dataframes
                main_classification_res["Disorder"] = disorder
                main_classification_res["Study"] = study
                main_classification_res["Analysis_Type"] = Analysis_Type
                main_classification_res["group_var"] = grouping_var
                main_classification_res["Classifier_Type"] = classifier_type
                main_classification_res["Site_Number"] = site_number
                main_classification_res["Fold"] = main_classification_res.index % num_folds
                main_classification_res["Repeat"] = main_classification_res.index // num_folds

                # Save results
                main_classification_res.to_feather(f"{data_path}/classification_results/robustness_analysis/ABIDE_ASD/{model_name}_Site{site_number}_main_classification_res.feather")
            except:
                print(f"Error with {model_name} for site {site_number}")

# Call on site 20
call_models_for_ABIDE_site(ASD_all_models=ASD_all_models, site_number=20, site_indices=ABIDE_site20_indices, class_labels=ABIDE_site20_class_labels, site_sample_IDs=ABIDE_site20_participants)

# Call on site 5
call_models_for_ABIDE_site(ASD_all_models=ASD_all_models, site_number=5, site_indices=ABIDE_site5_indices, class_labels=ABIDE_site5_class_labels, site_sample_IDs=ABIDE_site5_participants)