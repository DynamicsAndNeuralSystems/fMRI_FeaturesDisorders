# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import os
import dill
from pyspi.calculator import Calculator
from pyspi.data import Data
from copy import deepcopy
import csv
import random

# Set the seed
random.seed(127)

# github path
github_path = "/headnode1/abry4213/github/fMRI_FeaturesDisorders/"
# Define data path
SCZ_pydata_path = "/headnode1/abry4213/data/UCLA_Schizophrenia/raw_data/pydata/" 
ASD_pydata_path = "/headnode1/abry4213/data/ABIDE_ASD/raw_data/pydata/" 

# config file
config_file = github_path + "data_prep_and_QC/pyspi_QC_analysis/pyspi_di_gaussian_config.yaml"

################### Run all SPIs for the same brain region pair 100x ####################
def process_SPIs_100x(sample_ID, pydata_path, subject_data_file):
    # Read in subject's time-series data and convert to a numpy array
    subject_data = pd.read_csv(pydata_path + subject_data_file + ".csv", header=None).to_numpy()

    # Initialize Calculator object to copy in each iteration
    basecalc = Calculator()

    # Load Data object using subject's data
    data = Data(data=subject_data, dim_order="ps", name=sample_ID, normalise=True)

    # Initialize a list to store the di_gaussian estimate for the given
    # brain pair across iterations
    full_pyspi_res = []

    for i in range(1,101):

        # Create a deepcopy of the original basecalc
        calc = deepcopy(basecalc)

        # Load the data object
        calc.load_dataset(data)

        # Compute all SPIs
        calc.compute()
        
        # Extract results and append to full_pyspi_res list
        calc_res = calc.table
        calc_res_filt = calc_res.filter(items = ["proc-0"], axis = 0)
        
        # Convert index to column
        calc_res_filt.columns = calc_res_filt.columns.to_flat_index()
        calc_res_filt.reset_index(level=0, inplace=True)
        
        # Rename index as first brain region
        calc_res_filt = calc_res_filt.rename(columns={"index": "brain_region_from"})
        
        # Pivot data from wide to long
        calc_res_filt_long = pd.melt(calc_res_filt, id_vars="brain_region_from")
        calc_res_filt_long["SPI"], calc_res_filt_long["brain_region_to"] = calc_res_filt_long.variable.str
        
        # Remove variable column
        calc_res_filt_long = calc_res_filt_long.drop("variable", 1)
        
        # Filter by region from = 0 and region to = 1
        calc_res_filt_long = calc_res_filt_long.loc[(calc_res_filt_long["brain_region_from"] == "proc-0") & (calc_res_filt_long["brain_region_to"] == "proc-1")]

        # Add iteration number
        calc_res_filt_long["Iteration"] = i
        full_pyspi_res.append(calc_res_filt_long)
        
    full_pyspi_merged = pd.concat(full_pyspi_res)
    full_pyspi_merged.value = np.real(full_pyspi_merged.value)

    # Write resulting dataframe to a CSV
    full_pyspi_merged.to_csv(ASD_pydata_path + subject_data_file + "_all_SPIs.csv")


# UCLA Schizophrenia with AROMA+2P+GMR

# # sub-10159, left bankssts --> left entorhinal cortex
# process_SPIs_100x(sample_ID="sub-10159", 
#                   pydata_path=SCZ_pydata_path, 
#                   subject_data_file="sub-10159_lh_bankssts_lh_entorhinal")
# # sub-10527, left rostral anterior cingulate --> left caudal middle frontal
# process_SPIs_100x(sample_ID="sub-10527",
#                   pydata_path=SCZ_pydata_path,
#                   subject_data_file="sub-10527_lh_rostralanteriorcingulate_lh_caudalmiddlefrontal")


# ABIDE ASD with FC1000 noise-processing
process_SPIs_100x(sample_ID="10021451277603445196",
                  pydata_path=ASD_pydata_path,
                  subject_data_file="subject_10021451277603445196_precentral_angular")


