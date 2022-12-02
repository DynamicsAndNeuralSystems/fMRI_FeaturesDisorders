from pyspi.statistics.infotheory import directed_info
import pandas as pd
import numpy as np
import os
import dill
from pyspi.calculator import Calculator
from pyspi.data import Data
from copy import deepcopy
import csv
import random
from pyspi.statistics.infotheory import directed_info

# github path
github_path = "/headnode1/abry4213/github/fMRI_FeaturesDisorders/"
# Define data path
SCZ_pydata_path = "/headnode1/abry4213/data/UCLA_Schizophrenia/raw_data/pydata/" 
ASD_pydata_path = "/headnode1/abry4213/data/ABIDE_ASD/raw_data/pydata/" 

subject_data_file = "sub-10159_lh_bankssts_lh_entorhinal"
pydata_path = SCZ_pydata_path
subject_data = pd.read_csv(pydata_path + subject_data_file + ".csv", header=None).to_numpy()

####################### direct directed_info measurement ########################

x = subject_data[0]
y = subject_data[1]

di_gauss_vec = list()

# Run directed info
for i in range(1,101):
    test = directed_info().bivariate(subject_data[0], subject_data[1])
    di_gauss_vec.append(test)
# This outputs different results
    
# set the seed once outside the loop
random.seed(27)
di_gauss_one_seed = list()
for i in range(1,101):
    test = directed_info().bivariate(subject_data[0], subject_data[1])
    di_gauss_one_seed.append(test)
# This outputs the same result all 100 loops
    
# set the seed once per loop
di_gauss_one_seed_per_loop = list()
for i in range(1,101):
    random.seed(127)
    test = directed_info().bivariate(subject_data[0], subject_data[1])
    di_gauss_one_seed_per_loop.append(test)
# This outputs different results


####################### using pyspi calculator ########################

# config file
config_file = github_path + "data_prep_and_QC/pyspi_QC_analysis/pyspi_di_gaussian_config.yaml"
    
def process_SPIs_100x(pydata_object, config_file, set_seed = False):

    # Initialize Calculator object to copy in each iteration
    basecalc = Calculator(configfile=config_file)

    # Initialize a list to store the di_gaussian estimate for the given
    # brain pair across iterations
    full_pyspi_res = []

    for i in range(1,101):
        # Set seed within loop if indicated
        if set_seed:
            random.seed(127)

        # Create a deepcopy of the original basecalc
        calc = deepcopy(basecalc)

        # Load the data object
        calc.load_dataset(pydata_object)

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

    return(full_pyspi_merged)

# Load Data object using subject's data
data = Data(data=subject_data, dim_order="ps", name="sub-10159", normalise=True)

# Without setting the seed
di_gauss_pyspi_no_seed = process_SPIs_100x(pydata_object = data, 
    config_file = config_file, 
    set_seed = False)
di_gauss_pyspi_no_seed.to_csv(pydata_path + subject_data_file + "_di_gaussian_no_seed.csv")

# Set the seed once outside the loop
random.seed(27)
di_gauss_pyspi_one_seed = process_SPIs_100x(pydata_object = data, 
    config_file = config_file, 
    set_seed = False)

# set the seed once per loop
di_gauss_pyspi_seed_per_loop = process_SPIs_100x(pydata_object = data, 
    config_file = config_file, 
    set_seed = True)

