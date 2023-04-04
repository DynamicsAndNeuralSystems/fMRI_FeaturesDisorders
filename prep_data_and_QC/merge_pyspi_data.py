# Load modules
import numpy as np
import pandas as pd
import dill
import argparse
import os
import pyarrow.feather as feather

# Command-line arguments to parse
parser = argparse.ArgumentParser(description='Process inputs for pairwise data preparation.')
parser.add_argument('--data_path', default="/headnode1/abry4213/data/UCLA_CNP_ABIDE_ASD/", dest='data_path')
parser.add_argument('--dataset_ID', default="UCLA_CNP", dest='dataset_ID')
parser.add_argument('--pkl_file', default="calc14_pyspi.pkl", dest='pkl_file')
parser.add_argument('--pairwise_feature_set', default="pyspi14", dest='pairwise_feature_set')
parser.add_argument('--brain_region_lookup', default="UCLA_CNP_Brain_Region_Lookup.feather", dest='brain_region_lookup')
parser.add_argument('--noise_proc', dest='noise_proc')

# Parse arguments
args = parser.parse_args()
data_path = args.data_path
noise_proc = args.noise_proc
dataset_ID = args.dataset_ID
pkl_file = args.pkl_file
pairwise_feature_set = args.pairwise_feature_set
brain_region_lookup = args.brain_region_lookup

# data_path = "/headnode1/abry4213/data/UCLA_CNP/"
# dataset_ID = "UCLA_CNP"
# pkl_file = "calc_pyspi14.pkl"
# noise_proc="AROMA+2P+GMR"
# pairwise_feature_set = "pyspi14"
# brain_region_lookup = "UCLA_CNP_Brain_Region_Lookup.feather"

proc_data_path = data_path + "processed_data/"
pkl_data_path = data_path + "raw_data/numpy_files/"

def merge_calcs_into_df(proc_data_path, 
                        pkl_data_path, 
                        brain_region_lookup,
                        pairwise_feature_set, 
                        noise_proc):
    
    # Set noise label for file paths
    noise_label = noise_proc.replace("+", "_")
    
    # Check if feather data file already exists
    if not os.path.isfile(f"{proc_data_path}/{dataset_ID}_{noise_label}_{pairwise_feature_set}_cov_EmpiricalCovariance_filtered.feather"):
        # Where individual pickle data files will be read in from
        input_np_data_path = pkl_data_path + noise_label + "/"
    
        # Read in ROI index data
        ROI_lookup = pd.read_feather(brain_region_lookup)
        
        # Find subjects
        samples = os.listdir(input_np_data_path)
        
        # Initialise list for each subject's pyspi14 data
        samples_with_pyspi_data = []
        sample_data_list = []
        for sample in samples:
            try:
                # Load in subject's pyspi14 data and filter
                with open(input_np_data_path + sample + "/" + pkl_file, "rb") as f:
                    sample_data = dill.load(f)
                
                # Exclude data points between the same brain region
                sample_data = sample_data.query("brain_region_from != brain_region_to")
                
                # Add Sample_ID column
                sample_data["Sample_ID"] = sample
                
                # Append to list of pyspi14 data
                samples_with_pyspi_data.append(sample)
                sample_data_list.append(sample_data)
            except:
                print("Error for " + sample)

        # Save a list of subjects who have pyspi14 data to a feather file with a column name
        pd.DataFrame(samples_with_pyspi_data, columns=["Sample_ID"]).to_feather(f"{proc_data_path}/{dataset_ID}_filtered_sample_info_{noise_label}_{pairwise_feature_set}.feather")
        
        # Switch brain region indices for region names to/from
        full_pyspi_res = pd.concat(sample_data_list).reset_index()
        full_pyspi_res["brain_region_from"] = full_pyspi_res["brain_region_from"].str.replace("proc-", "").astype(int) + 1
        full_pyspi_res["brain_region_to"] = full_pyspi_res["brain_region_to"].str.replace("proc-", "").astype(int) + 1
        
        full_pyspi_res_ROI = ((pd.merge(full_pyspi_res, 
                                      ROI_lookup,
                                      how = "left",
                                      left_on = "brain_region_from",
                                      right_on = "Index")
                              .drop(["index", "brain_region_from", "Index"], axis=1)
                              .rename(columns={"Brain_Region": "brain_region_from"})
                              ).merge(ROI_lookup,
                                      how="left",
                                      left_on = "brain_region_to",
                                      right_on = "Index")
                                      .drop(["brain_region_to", "Index"], axis=1)
                                      .rename(columns={"Brain_Region": "brain_region_to"})
                                      )
        # Get a list of unique values in the SPI column
        SPI_list = full_pyspi_res_ROI["SPI"].unique().tolist()

        # Split into individual dataframes per SPI
        full_pyspi_res_ROI_split = [full_pyspi_res_ROI.query(f"SPI == '{i}'") for i in SPI_list]

        # Save individual SPI dataframes into feather files with the SPI in the filename
        for i in range(1, 15):
            SPI_name = SPI_list[i-1]
            feather.write_feather(full_pyspi_res_ROI_split[i-1], f"{proc_data_path}/{dataset_ID}_{noise_label}_{pairwise_feature_set}_{SPI_name}_filtered.feather", version=1)
        

def filter_pyspi_data(proc_data_path,
                      dataset_ID,
                      noise_proc,
                      pairwise_feature_set):
    
    # Set noise label for file paths
    noise_label = noise_proc.replace("+", "_")
    
    # Check if feather data file already exists
    if not os.path.isfile(f"{proc_data_path}/{dataset_ID}_{noise_label}_{pairwise_feature_set}_filtered.feather"):
        raw_pyspi_res = pd.read_feather(f"{proc_data_path}/{dataset_ID}_{noise_label}_{pairwise_feature_set}.feather")
        
        # Merge data into region pairs
        raw_pyspi_res["Region_Pair"] = raw_pyspi_res.brain_region_from + "_" + raw_pyspi_res.brain_region_to
        
        # Find number of unique region pairs
        num_region_pairs = raw_pyspi_res.Region_Pair.nunique()
        
        # Find number of unique SPIs
        num_SPIs = raw_pyspi_res.SPI.nunique()    
        
        # Find all NA data
        all_NA_data = (raw_pyspi_res.loc[pd.isnull(raw_pyspi_res.value)]
                       .groupby(["Sample_ID", "SPI"])
                       .count())
        
        samples_to_drop = []
        
        # Check if there are any samples that yielded all NaN
        for index, row in all_NA_data.iterrows():
            sample = index[0]
            sample_num_NaN_region_pairs = row["Region_Pair"]
            
            # If the number of NaN region pairs for this subject/SPI equals the total
            # number of unique region pairs, drop that sample ID
            if sample_num_NaN_region_pairs == num_region_pairs:
                samples_to_drop.append(sample)
        
        # Drop any samples retained in samples_to_drop
        filtered_pyspi_res = (raw_pyspi_res[~raw_pyspi_res.Sample_ID.isin(samples_to_drop)]
                              .drop(["Region_Pair"], axis=1))


        # Save filtered pyspi results to a feather file
        feather.write_feather(filtered_pyspi_res, f"{proc_data_path}/{dataset_ID}_{noise_label}_{pairwise_feature_set}_filtered.feather", version=1)
       
        
merge_calcs_into_df(proc_data_path = proc_data_path,
                    pkl_data_path = pkl_data_path,
                    brain_region_lookup = data_path + "study_metadata/" + brain_region_lookup,
                    pairwise_feature_set = pairwise_feature_set,
                    noise_proc = noise_proc)
# filter_pyspi_data(proc_data_path = proc_data_path,
#                   dataset_ID = dataset_ID,
#                   noise_proc = noise_proc,
#                   pairwise_feature_set = pairwise_feature_set)


