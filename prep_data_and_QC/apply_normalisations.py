
from sklearnex import patch_sklearn
patch_sklearn()
from multiprocessing import Pool, Process

from numpy.random import randint
import typing
import pandas as pd
import numpy as np
import math
import os
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.preprocessing import StandardScaler
import sys
from pyarrow import feather
sys.path.append("/headnode1/abry4213/github/fMRI_FeaturesDisorders/helper_functions/classification/")
from robust_sigmoid_normalisation import RobustSigmoidScaler

# Define data paths
UCLA_CNP_data_path = "/headnode1/abry4213/data/UCLA_CNP/processed_data/"
ABIDE_ASD_data_path = "/headnode1/abry4213/data/ABIDE_ASD/processed_data/"

# Define normalisation function
def apply_transform_by_region(input_data, transform_type, output_file):
    
    if not os.path.isfile(output_file):
        # Initialise list for transformed data by region
        region_transformed_list = []
        
        for brain_region in input_data.Brain_Region.unique():
            region_data = input_data.query("Brain_Region == @brain_region")
            
            # Pivot from long to wide
            region_data_wide = region_data.pivot(index = "Sample_ID", columns = "names", values = "values")    
            
            # Apply given transformation
            if transform_type == "z-score":
                transformer = StandardScaler().fit(region_data_wide)
                
            else:
                transformer = RobustSigmoidScaler(unit_variance=True).fit(region_data_wide)
                
            region_data_transformed = pd.DataFrame(transformer.transform(region_data_wide),
                                                columns = region_data_wide.columns)
            
            # Add column for sample ID
            region_data_transformed["Sample_ID"] = region_data_wide.index
            
            # Pivot back to long
            region_data_transformed_long = region_data_transformed.melt(id_vars = "Sample_ID", 
                                                                        var_name = "names",
                                                                        value_name = "values")
            
        
            # Add relevant info
            region_data_transformed_long["Brain_Region"] = brain_region
            region_data_transformed_long["Normalisation"] = transform_type
            
            # Append to list
            region_transformed_list.append(region_data_transformed_long)
        
        # Concatenate transformed data into one dataframe
        region_transformed_data = pd.concat(region_transformed_list)
        
        # Save transformed data
        region_transformed_data = region_transformed_data.reset_index()
        region_transformed_data.to_feather(output_file)

####################### Z-score #######################

if __name__ == '__main__':

    # Load all needed data
    UCLA_CNP_catch22_data = feather.read_feather(f"{UCLA_CNP_data_path}/UCLA_CNP_AROMA_2P_GMR_catch22_filtered.feather")

    ABIDE_ASD_catch22_data = feather.read_feather(f"{ABIDE_ASD_data_path}/ABIDE_ASD_FC1000_catch22_filtered.feather")

    UCLA_CNP_pyspi14_data = feather.read_feather(f"{UCLA_CNP_data_path}/UCLA_CNP_AROMA_2P_GMR_pyspi14_filtered.feather")
    UCLA_CNP_pyspi14_data.rename(columns={"SPI": "names", "value": "values"}, inplace=True)
    UCLA_CNP_pyspi14_data["Brain_Region"] = UCLA_CNP_pyspi14_data["brain_region_from"].astype(str) + '_' + UCLA_CNP_pyspi14_data["brain_region_to"].astype(str)

    ABIDE_ASD_pyspi14_data = feather.read_feather(f"{ABIDE_ASD_data_path}/ABIDE_ASD_FC1000_pyspi14_filtered.feather")
    ABIDE_ASD_pyspi14_data.rename(columns={"SPI": "names", "value": "values"}, inplace=True)
    ABIDE_ASD_pyspi14_data["Brain_Region"] = ABIDE_ASD_pyspi14_data["brain_region_from"].astype(str) + '_' + ABIDE_ASD_pyspi14_data["brain_region_to"].astype(str)

    # Define z-score processes
    UCLA_CNP_catch22_zscore_p = Process(target=apply_transform_by_region, args=(UCLA_CNP_catch22_data, "z-score", f"{UCLA_CNP_data_path}/UCLA_CNP_AROMA_2P_GMR_catch22_filtered_zscored.feather"))
    ABIDE_ASD_catch22_zscore_p = Process(target=apply_transform_by_region, args=(ABIDE_ASD_catch22_data, "z-score", f"{ABIDE_ASD_data_path}/ABIDE_ASD_FC1000_catch22_filtered_zscored.feather"))
    UCLA_CNP_pyspi14_zscore_p = Process(target=apply_transform_by_region, args=(UCLA_CNP_pyspi14_data, "z-score", f"{UCLA_CNP_data_path}/UCLA_CNP_AROMA_2P_GMR_pyspi14_filtered_zscored.feather"))
    ABIDE_ASD_pyspi14_zscore_p = Process(target=apply_transform_by_region, args=(ABIDE_ASD_pyspi14_data, "z-score", f"{ABIDE_ASD_data_path}/ABIDE_ASD_FC1000_pyspi14_filtered_zscored.feather"))

    # Define robust sigmoid processes
    UCLA_CNP_catch22_RobustSigmoid_p = Process(target=apply_transform_by_region, args=(UCLA_CNP_catch22_data, "RobustSigmoid", f"{UCLA_CNP_data_path}/UCLA_CNP_AROMA_2P_GMR_catch22_filtered_RobustSigmoid.feather"))
    ABIDE_ASD_catch22_RobustSigmoid_p = Process(target=apply_transform_by_region, args=(ABIDE_ASD_catch22_data, "RobustSigmoid", f"{ABIDE_ASD_data_path}/ABIDE_ASD_FC1000_catch22_filtered_RobustSigmoid.feather"))
    UCLA_CNP_pyspi14_RobustSigmoid_p = Process(target=apply_transform_by_region, args=(UCLA_CNP_pyspi14_data, "RobustSigmoid", f"{UCLA_CNP_data_path}/UCLA_CNP_AROMA_2P_GMR_pyspi14_filtered_RobustSigmoid.feather"))
    ABIDE_ASD_pyspi14_RobustSigmoid_p = Process(target=apply_transform_by_region, args=(ABIDE_ASD_pyspi14_data, "RobustSigmoid", f"{ABIDE_ASD_data_path}/ABIDE_ASD_FC1000_pyspi14_filtered_RobustSigmoid.feather"))

    # Start the processes
    UCLA_CNP_catch22_zscore_p.start()
    ABIDE_ASD_catch22_zscore_p.start()
    UCLA_CNP_pyspi14_zscore_p.start()
    ABIDE_ASD_pyspi14_zscore_p.start()
    UCLA_CNP_catch22_RobustSigmoid_p.start()
    ABIDE_ASD_catch22_RobustSigmoid_p.start()
    UCLA_CNP_pyspi14_RobustSigmoid_p.start()
    ABIDE_ASD_pyspi14_RobustSigmoid_p.start()

    # Wait for the processes to finish
    UCLA_CNP_catch22_zscore_p.join()
    ABIDE_ASD_catch22_zscore_p.join()
    UCLA_CNP_pyspi14_zscore_p.join()
    ABIDE_ASD_pyspi14_zscore_p.join()
    UCLA_CNP_catch22_RobustSigmoid_p.join()
    ABIDE_ASD_catch22_RobustSigmoid_p.join()
    UCLA_CNP_pyspi14_RobustSigmoid_p.join()
    ABIDE_ASD_pyspi14_RobustSigmoid_p.join()