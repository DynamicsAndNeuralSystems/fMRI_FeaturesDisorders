# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import os
import dill
from pyspi.calculator import Calculator
from pyspi.data import Data

# github path
github_path = "/headnode1/abry4213/github/fMRI_FeaturesDisorders/"
# Define data path
SCZ_data_path = "/headnode1/abry4213/data/UCLA_Schizophrenia/raw_data/numpy_files/AROMA_2P_GMR/" 
ASD_data_path = "/headnode1/abry4213/data/ABIDE_ASD/raw_data/numpy_files/FC1000/" 

############################# di_gaussian redo #################################
# config file
config_file = github_path + "data_prep_and_QC/pyspi_QC_analysis/pyspi_di_gaussian_config.yaml"

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

