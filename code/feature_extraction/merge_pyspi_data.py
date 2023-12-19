# Load modules
import pandas as pd
import dill
import argparse
import os
import pyarrow.feather as feather

# Command-line arguments to parse
parser = argparse.ArgumentParser(description='Process inputs for pairwise data preparation.')
parser.add_argument('--data_path', default="data/", dest='data_path')
parser.add_argument('--dataset_ID', default="UCLA_CNP", dest='dataset_ID')
parser.add_argument('--pkl_file', default="calc.pkl", dest='pkl_file')
parser.add_argument('--pairwise_feature_set', default="pyspi14", dest='pairwise_feature_set')
parser.add_argument('--brain_region_lookup', default="UCLA_CNP_Brain_Region_Lookup.feather", dest='brain_region_lookup')

# Parse arguments
args = parser.parse_args()
data_path = args.data_path
dataset_ID = args.dataset_ID
pkl_file = args.pkl_file
pairwise_feature_set = args.pairwise_feature_set
brain_region_lookup = args.brain_region_lookup

# Define output data paths
output_data_path = data_path + 'time_series_features/'
pkl_data_path = f"{data_path}/time_series_features/{dataset_ID}/"

def merge_calcs_into_df(output_data_path, 
                        pkl_data_path, 
                        dataset_ID,
                        brain_region_lookup,
                        pairwise_feature_set):
    
    # Check if feather data file already exists
    if not os.path.isfile(f"{output_data_path}/{dataset_ID}_{pairwise_feature_set}.feather"):
    
        # Read in ROI index data
        ROI_lookup = pd.read_feather(brain_region_lookup)
        
        # Find subjects
        samples = os.listdir(pkl_data_path)
        
        # Initialise list for each subject's pyspi14 data
        sample_data_list = []
        for sample in samples:
            try:
                # Load in subject's pyspi14 data and filter
                with open(pkl_data_path + sample + "/" + pkl_file, "rb") as f:
                    sample_data = dill.load(f)
                
                # Exclude data points between the same brain region
                sample_data = sample_data.query("brain_region_from != brain_region_to")
                
                # Add Sample_ID column
                sample_data["Sample_ID"] = sample
                
                # Append to list of pyspi14 data
                sample_data_list.append(sample_data)
            except:
                print("Error for " + sample)
        
        
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
        
        # Save merged res into a feather file
        feather.write_feather(full_pyspi_res_ROI, f"{output_data_path}/{dataset_ID}_{pairwise_feature_set}.feather", version=1)

def filter_pyspi_data(output_data_path,
                      dataset_ID,
                      pairwise_feature_set):
    
    
    # Check if feather data file already exists
    if not os.path.isfile(f"{output_data_path}/{dataset_ID}_{pairwise_feature_set}_filtered.feather"):
        raw_pyspi_res = pd.read_feather(f"{output_data_path}/{dataset_ID}_{pairwise_feature_set}.feather")
        
        # Merge data into region pairs
        raw_pyspi_res["Region_Pair"] = raw_pyspi_res.brain_region_from + "_" + raw_pyspi_res.brain_region_to
        
        # Find number of unique region pairs
        num_region_pairs = raw_pyspi_res.Region_Pair.nunique()
        
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
        feather.write_feather(filtered_pyspi_res, f"{output_data_path}/{dataset_ID}_{pairwise_feature_set}_filtered.feather", version=1)
             
merge_calcs_into_df(output_data_path = output_data_path,
                    pkl_data_path = pkl_data_path,
                    dataset_ID = dataset_ID,
                    brain_region_lookup = brain_region_lookup,
                    pairwise_feature_set = pairwise_feature_set)
filter_pyspi_data(output_data_path = output_data_path,
                  dataset_ID = dataset_ID,
                  pairwise_feature_set = pairwise_feature_set)


