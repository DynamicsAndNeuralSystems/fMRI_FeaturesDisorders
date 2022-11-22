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
github_path = "/Users/abry4213/github/fMRI_FeaturesDisorders/data_visualisation/quality_control/"

# Define data path
SCZ_data_path = "/Users/abry4213/data/UCLA_Schizophrenia/raw_data/numpy_files/AROMA_2P_GMR/" 

subject = "sub-10171"
subject_data = np.load(SCZ_data_path + subject + ".npy")

# config file
config_file = github_path + "pyspi_di_gaussian_config.yaml"

calc = Calculator(dataset=subject_data, configfile=config_file)
calc.compute()

calc_res = calc.table


# Iterate over each SPI
calc_res.columns = calc_res.columns.to_flat_index()

# Convert index to column
calc_res.reset_index(level=0, inplace=True)

 # Rename index as first brain region
calc_res = calc_res.rename(columns={"index": "brain_region_1"})

# Pivot data from wide to long
SPI_res_long = pd.melt(calc_res, id_vars="brain_region_1")
SPI_res_long['SPI'], SPI_res_long['brain_region_2'] = SPI_res_long.variable.str

 # Remove variable column
SPI_res_long = SPI_res_long.drop("variable", 1)

SPI_res_long.to_csv("temp.csv")