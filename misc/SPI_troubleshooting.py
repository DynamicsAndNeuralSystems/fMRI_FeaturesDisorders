# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import numpy as np
from pyspi.calculator import Calculator
from pyspi.data import Data


# github path
github_path = "/Users/abry4213/github/fMRI_FeaturesDisorders/misc/"

# Define data path
SCZ_data_path = "/Users/abry4213/data/UCLA_Schizophrenia/raw_data/numpy_files/AROMA_2P_GMR/" 

subject = "sub-10527"
subject_data = np.load(SCZ_data_path + subject + ".npy")

# config file
config_file = github_path + "pyspi_sgc_config.yaml"

calc = Calculator(dataset=subject_data, configfile=config_file)

calc_res = calc.table
