import pandas as pd
import scipy.io 
import sys
import os

input_feather_file_base = sys.argv[1]
# input_feather_file_base =  "/headnode1/abry4213/data/UCLA_CNP/raw_data/UCLA_CNP_AROMA_2P_GMR_fMRI_TS"

if not os.path.isfile(input_feather_file_base+".mat"):
    loaded_data = pd.read_feather(input_feather_file_base+".feather")
    # Iterate over all subjects and brain regions
    unique_subjects = loaded_data.Sample_ID.unique().tolist()
    unique_regions = loaded_data.Brain_Region.unique().tolist()
    
    # Initialize list to store dataframes
    list_of_subjects_and_regions = []
    
    for sample_ID in unique_subjects:
        for brain_region in unique_regions:
            sample_data = loaded_data.query("Sample_ID == @sample_ID & Brain_Region == @brain_region")
            sample_TS = sample_data["values"].tolist()
            sample_dict = {'Sample_ID': sample_ID,
                            'Brain_Region': brain_region,
                            'Time_Series': [sample_TS]}
            
            sample_df = pd.DataFrame(data = sample_dict)
            
            list_of_subjects_and_regions.append(sample_df)
            
    # Concatenate list of dataframes
    subjects_and_regions = pd.concat(list_of_subjects_and_regions, axis=0)
    
    # Save to a .mat file
    scipy.io.savemat(input_feather_file_base+".mat",
                     {'feather_to_mat_struct': subjects_and_regions.to_dict('list')})

    
