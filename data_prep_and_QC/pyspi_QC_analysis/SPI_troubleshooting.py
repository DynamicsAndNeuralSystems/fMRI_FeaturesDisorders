# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import os
import dill
from pyspi.calculator import Calculator
from pyspi.data import Data
from copy import deepcopy
import csv

# github path
github_path = "/headnode1/abry4213/github/fMRI_FeaturesDisorders/"
# Define data path
SCZ_data_path = "/headnode1/abry4213/data/UCLA_Schizophrenia/raw_data/numpy_files/AROMA_2P_GMR/" 
ASD_data_path = "/headnode1/abry4213/data/ABIDE_ASD/raw_data/numpy_files/FC1000/" 

# config file
config_file = github_path + "data_prep_and_QC/pyspi_QC_analysis/pyspi_di_gaussian_config.yaml"

############################# di_gaussian redo #################################

############################ UCLA Schizophrenia ################################

# subject_list = ["sub-10206", "sub-10292", "sub-10438", "sub-10624",
#                 "sub-10674", "sub-11019", "sub-11066", "sub-11122",
#                 "sub-11131", "sub-11143", "sub-50015", "sub-50016",
#                 "sub-50022"]

# for subject in subject_list:
#     if not os.path.exists(SCZ_data_path + subject + "/calc_di_gaussian.pkl"):
#         subject_data = np.load(SCZ_data_path + subject + ".npy")
    
#         calc = Calculator(dataset=subject_data, configfile=config_file)
#         calc.compute()
    
#         calc_res = calc.table
    
#         # Save calc results to a pickle file
#         with open(SCZ_data_path + subject + "/calc_di_gaussian.pkl", 'wb') as f:
#             dill.dump(calc_res, f)

################################ ABIDE ASD ####################################

# # ABIDE ASD di_gaussian only, round 1
# subjects_to_run = pd.read_csv(github_path + "data_prep_and_QC/pyspi_QC_analysis/ABIDE_ASD_di_gaussian_NaN_subjects.csv")["x"].tolist()
# subjects_to_run = [str(s) for s in subjects_to_run]

# for subject in subjects_to_run:
#         if not os.path.exists(ASD_data_path + subject + "/calc_di_gaussian.pkl"):
#             print("Now running di_gaussian for " + subject + "\n")
#             subject_data = np.load(ASD_data_path + subject + ".npy")
        
#             calc = Calculator(dataset=subject_data, configfile=config_file)
#             calc.compute()
        
#             calc_res = calc.table
        
#             # Save calc results to a pickle file
#             with open(ASD_data_path + subject + "/calc_di_gaussian.pkl", 'wb') as f:
#                 dill.dump(calc_res, f)
                
# # ABIDE ASD di_gaussian only, round 2
# subjects_to_run_v2 = pd.read_csv(github_path + "data_prep_and_QC/pyspi_QC_analysis/ABIDE_ASD_di_gaussian_NaN_subjects_v2.csv")["x"].tolist()
# subjects_to_run_v2 = [str(s) for s in subjects_to_run_v2]

# for subject in subjects_to_run_v2:
#         if not os.path.exists(ASD_data_path + subject + "/calc_di_gaussian_v2.pkl"):
#             print("Now running di_gaussian round 2 for " + subject + "\n")
#             subject_data = np.load(ASD_data_path + subject + ".npy")
        
#             calc = Calculator(dataset=subject_data, configfile=config_file)
#             calc.compute()
        
#             calc_res = calc.table
        
#             # Save calc results to a pickle file
#             with open(ASD_data_path + subject + "/calc_di_gaussian_v2.pkl", 'wb') as f:
#                 dill.dump(calc_res, f)

################### Run di_gaussian the same brain region pair 1000x ####################

# Load sub-10159 AROMA+2P+GMR left bankssts and left entorhinal cortex time-series
# prepared in SPI_troubleshooting.R
SCZ_pydata_path = "/headnode1/abry4213/data/UCLA_Schizophrenia/raw_data/pydata/" 
subject_data = pd.read_csv(SCZ_pydata_path + "sub-10159_lh_bankssts_lh_entorhinal.csv",
                           header = None).to_numpy()

# Initialize Calculator object to copy in each iteration
# using custom di_gaussian config file
basecalc = Calculator(configfile=config_file)

# Initialize a list to store the di_gaussian estimate for the given
# brain pair across iterations
di_gauss_res = []

for i in range(1,1001):
    # Load Data object using subject's data
    data = Data(data=subject_data,dim_order="ps",name="sub-10159",
                normalise=True, procnames=["ctx-lh-bankssts",
                                           "ctx-lh-entorhinal"])

    # Create a deepcopy of the original basecalc
    calc = deepcopy(basecalc)

    # Load the data object
    calc.load_dataset(data)

    # Compute di_gaussian SPI
    calc.compute()
    
    # Extract results and append to di_gauss_res list
    calc_res = calc.table
    lh_bankssts_to_lh_entorhinal = calc_res.iloc[0,1]
    di_gauss_res.append(lh_bankssts_to_lh_entorhinal)
    
# Write resulting list of 1,000 di_gaussian values to a CSV
with open(SCZ_pydata_path + "sub-10159_lh_bankssts_lh_entorhinal_di_gaussian.csv", "w") as f:
    writer = csv.writer(f)
    for val in di_gauss_res:
        writer.writerow([val])