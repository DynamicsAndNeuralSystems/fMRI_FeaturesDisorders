# Parse command-line arguments
import argparse
import os

parser = argparse.ArgumentParser(description="SPI troubleshooting script.")

parser.add_argument('--sample_ID', dest='sample_ID',
                    help='Sample ID')
parser.add_argument('--TS_file', dest='TS_file',
                    help="File path to CSV file containing the two brain regions' time-series")       
parser.add_argument('--data_out_file', dest='data_out_file',
                    help="File path to the output CSV file containing SPI results")  
parser.add_argument('--set_global_seed', dest='set_global_seed',
                    default = False, action="store_true",
                    help="Flag to set a global seed, outside of the for loop")  
parser.add_argument('--set_for_loop_seed', dest='set_for_loop_seed',
                    default = False, action="store_true",
                    help="Flag to set a seed within the for loop")

# Parse the arguments
args = parser.parse_args()
sample_ID = args.sample_ID
TS_file = args.TS_file
data_out_file = args.data_out_file
set_global_seed = args.set_global_seed
set_for_loop_seed = args.set_for_loop_seed

# Import other needed modules
import pandas as pd
import numpy as np
from pyspi.calculator import Calculator
from pyspi.data import Data
from copy import deepcopy
import random

# Set global seed to 27 if flag was indicated
if set_global_seed:
    random.seed(27)

################### Run all SPIs for the same brain region pair 100x ####################
def process_SPIs_100x(sample_ID, TS_file, data_out_file, set_loop_seed = False):

    # sub-10527, left rostral anterior cingulate --> left caudal middle frontal
    # Read in subject's time-series data and convert to a numpy array
    subject_data = pd.read_csv(TS_file, header=None).to_numpy()

    calc = Calculator(Data(subject_data, dim_order="ps", normalise=True))
    calc.compute()

    # Initialize Calculator object to copy in each iteration
    basecalc = Calculator()

    # Load Data object using subject's data
    data = Data(data=subject_data, dim_order="ps", name=sample_ID, normalise=True)

    # Initialize a list to store the di_gaussian estimate for the given
    # brain pair across iterations
    full_pyspi_res = []

    for i in range(1,101):
        # Set seed within loop if indicated
        if set_loop_seed:
            random.seed(127)

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
    full_pyspi_merged.to_csv(data_out_file)



############ Run all SPIs for the same brain region pair 100x ############

process_SPIs_100x(sample_ID=sample_ID,
                  TS_file = TS_file,
                  data_out_file = data_out_file,
                  set_loop_seed = set_for_loop_seed)

# # Seeded
# random.seed(127)
# process_SPIs_100x(sample_ID="sub-10527",
#                   subject_data = subject_data,
#                   pydata_path=SCZ_pydata_path,
#                   data_out_file="sub-10527_lh_rostralanteriorcingulate_lh_caudalmiddlefrontal_seeded")