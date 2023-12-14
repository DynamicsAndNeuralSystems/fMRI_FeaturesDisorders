
from sklearnex import patch_sklearn
patch_sklearn()
from multiprocessing import Process

import pandas as pd
import os
from sklearn.preprocessing import StandardScaler
import sys
from pyarrow import feather
sys.path.append("/headnode1/abry4213/github/fMRI_FeaturesDisorders/helper_functions/classification/")
from mixed_sigmoid_normalisation import MixedSigmoidScaler

# Define data paths
UCLA_CNP_data_path = "/headnode1/abry4213/data/UCLA_CNP/processed_data/"
ABIDE_ASD_data_path = "/headnode1/abry4213/data/ABIDE_ASD/processed_data/"
final_data_path = "/headnode1/abry4213/data/TS_feature_manuscript/"

# Helper function to bin data
def bin_feature_values(input_df):

    results_list = []
    for feature in input_df.names.unique().tolist():
        df_feature = input_df.query("names==@feature")
        df_feature['bin'] = pd.cut(df_feature['values'], bins=100)
        
        df_feature_binned = (df_feature
                             .groupby(["names", "Normalisation", "bin"])
                             .agg({"bin": "count"})
                             .rename(columns={"bin": "count"})
                             .reset_index()
                             )
        
        results_list.append(df_feature_binned)
        
    all_binned_res = pd.concat(results_list, axis=0).reset_index()
    return all_binned_res

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
                transformer = MixedSigmoidScaler(unit_variance=True).fit(region_data_wide)
                
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
        
        # Bin data for raw values
        region_transformed_data_counts = bin_feature_values(region_transformed_data)

        # Fix bug with parentheses in feather file
        region_transformed_data_counts['bin'] = region_transformed_data_counts['bin'].astype(str).str.replace('(', '[')
        region_transformed_data_counts.to_feather(output_file)

####################### Z-score #######################

    
if __name__ == '__main__':

    # Load all needed data
    UCLA_CNP_catch25_data = feather.read_feather(f"{UCLA_CNP_data_path}/UCLA_CNP_AROMA_2P_GMR_catch25_filtered.feather")

    ABIDE_ASD_catch25_data = feather.read_feather(f"{ABIDE_ASD_data_path}/ABIDE_ASD_FC1000_catch25_filtered.feather")

    UCLA_CNP_pyspi14_data = feather.read_feather(f"{UCLA_CNP_data_path}/UCLA_CNP_AROMA_2P_GMR_pyspi14_filtered.feather")
    UCLA_CNP_pyspi14_data.rename(columns={"SPI": "names", "value": "values"}, inplace=True)
    UCLA_CNP_pyspi14_data["Brain_Region"] = UCLA_CNP_pyspi14_data["brain_region_from"].astype(str) + '_' + UCLA_CNP_pyspi14_data["brain_region_to"].astype(str)

    ABIDE_ASD_pyspi14_data = feather.read_feather(f"{ABIDE_ASD_data_path}/ABIDE_ASD_FC1000_pyspi14_filtered.feather")
    ABIDE_ASD_pyspi14_data.rename(columns={"SPI": "names", "value": "values"}, inplace=True)
    ABIDE_ASD_pyspi14_data["Brain_Region"] = ABIDE_ASD_pyspi14_data["brain_region_from"].astype(str) + '_' + ABIDE_ASD_pyspi14_data["brain_region_to"].astype(str)

    # Define z-score processes
    UCLA_CNP_catch25_zscore_p = Process(target=apply_transform_by_region, args=(UCLA_CNP_catch25_data, "z-score", f"{UCLA_CNP_data_path}/UCLA_CNP_AROMA_2P_GMR_catch25_filtered_zscored_counts.feather"))
    ABIDE_ASD_catch25_zscore_p = Process(target=apply_transform_by_region, args=(ABIDE_ASD_catch25_data, "z-score", f"{ABIDE_ASD_data_path}/ABIDE_ASD_FC1000_catch25_filtered_zscored_counts.feather"))
    UCLA_CNP_pyspi14_zscore_p = Process(target=apply_transform_by_region, args=(UCLA_CNP_pyspi14_data, "z-score", f"{UCLA_CNP_data_path}/UCLA_CNP_AROMA_2P_GMR_pyspi14_filtered_zscored_counts.feather"))
    ABIDE_ASD_pyspi14_zscore_p = Process(target=apply_transform_by_region, args=(ABIDE_ASD_pyspi14_data, "z-score", f"{ABIDE_ASD_data_path}/ABIDE_ASD_FC1000_pyspi14_filtered_zscored_counts.feather"))

    # Define mixed sigmoid processes
    UCLA_CNP_catch25_MixedSigmoid_p = Process(target=apply_transform_by_region, args=(UCLA_CNP_catch25_data, "MixedSigmoid", f"{UCLA_CNP_data_path}/UCLA_CNP_AROMA_2P_GMR_catch25_filtered_MixedSigmoid_counts.feather"))
    ABIDE_ASD_catch25_MixedSigmoid_p = Process(target=apply_transform_by_region, args=(ABIDE_ASD_catch25_data, "MixedSigmoid", f"{ABIDE_ASD_data_path}/ABIDE_ASD_FC1000_catch25_filtered_MixedSigmoid_counts.feather"))
    UCLA_CNP_pyspi14_MixedSigmoid_p = Process(target=apply_transform_by_region, args=(UCLA_CNP_pyspi14_data, "MixedSigmoid", f"{UCLA_CNP_data_path}/UCLA_CNP_AROMA_2P_GMR_pyspi14_filtered_MixedSigmoid_counts.feather"))
    ABIDE_ASD_pyspi14_MixedSigmoid_p = Process(target=apply_transform_by_region, args=(ABIDE_ASD_pyspi14_data, "MixedSigmoid", f"{ABIDE_ASD_data_path}/ABIDE_ASD_FC1000_pyspi14_filtered_MixedSigmoid_counts.feather"))

    # Start the processes
    UCLA_CNP_catch25_zscore_p.start()
    ABIDE_ASD_catch25_zscore_p.start()
    UCLA_CNP_pyspi14_zscore_p.start()
    ABIDE_ASD_pyspi14_zscore_p.start()
    UCLA_CNP_catch25_MixedSigmoid_p.start()
    ABIDE_ASD_catch25_MixedSigmoid_p.start()
    UCLA_CNP_pyspi14_MixedSigmoid_p.start()
    ABIDE_ASD_pyspi14_MixedSigmoid_p.start()

    # Wait for the processes to finish
    UCLA_CNP_catch25_zscore_p.join()
    ABIDE_ASD_catch25_zscore_p.join()
    UCLA_CNP_pyspi14_zscore_p.join()
    ABIDE_ASD_pyspi14_zscore_p.join()
    UCLA_CNP_catch25_MixedSigmoid_p.join()
    ABIDE_ASD_catch25_MixedSigmoid_p.join()
    UCLA_CNP_pyspi14_MixedSigmoid_p.join()
    ABIDE_ASD_pyspi14_MixedSigmoid_p.join()

    # Bin data for raw values
    UCLA_CNP_catch25_raw_binned = bin_feature_values(UCLA_CNP_catch25_data.assign(Normalisation = "Raw_Values"))
    ABIDE_ASD_catch25_raw_binned = bin_feature_values(ABIDE_ASD_catch25_data.assign(Normalisation = "Raw_Values"))
    UCLA_CNP_pyspi14_raw_binned = bin_feature_values(UCLA_CNP_pyspi14_data.assign(Normalisation = "Raw_Values"))
    ABIDE_ASD_pyspi14_raw_binned = bin_feature_values(ABIDE_ASD_pyspi14_data.assign(Normalisation = "Raw_Values"))

    # Merge the results into one dataframe
    if not os.path.isfile(f"{final_data_path}/all_normalisations_counts.feather"):
        UCLA_CNP_catch25_data_counts = bin_feature_values(UCLA_CNP_catch25_data.assign(Normalisation = "Raw_Values"))
        ABIDE_ASD_catch25_data_counts = bin_feature_values(ABIDE_ASD_catch25_data.assign(Normalisation = "Raw_Values"))
        UCLA_CNP_pyspi14_data_counts = bin_feature_values(UCLA_CNP_pyspi14_data.assign(Normalisation = "Raw_Values"))
        ABIDE_ASD_pyspi14_data_counts = bin_feature_values(ABIDE_ASD_pyspi14_data.assign(Normalisation = "Raw_Values"))

        UCLA_CNP_catch25_data_z = feather.read_feather(f"{UCLA_CNP_data_path}/UCLA_CNP_AROMA_2P_GMR_catch25_filtered_zscored_counts.feather")
        ABIDE_ASD_catch25_data_z = feather.read_feather(f"{ABIDE_ASD_data_path}/ABIDE_ASD_FC1000_catch25_filtered_zscored_counts.feather")
        UCLA_CNP_pyspi14_data_z = feather.read_feather(f"{UCLA_CNP_data_path}/UCLA_CNP_AROMA_2P_GMR_pyspi14_filtered_zscored_counts.feather")
        ABIDE_ASD_pyspi14_data_z = feather.read_feather(f"{ABIDE_ASD_data_path}/ABIDE_ASD_FC1000_pyspi14_filtered_zscored_counts.feather")

        UCLA_CNP_catch25_data_MS = feather.read_feather(f"{UCLA_CNP_data_path}/UCLA_CNP_AROMA_2P_GMR_catch25_filtered_MixedSigmoid_counts.feather")
        ABIDE_ASD_catch25_data_MS = feather.read_feather(f"{ABIDE_ASD_data_path}/ABIDE_ASD_FC1000_catch25_filtered_MixedSigmoid_counts.feather")
        UCLA_CNP_pyspi14_data_MS = feather.read_feather(f"{UCLA_CNP_data_path}/UCLA_CNP_AROMA_2P_GMR_pyspi14_filtered_MixedSigmoid_counts.feather")
        ABIDE_ASD_pyspi14_data_MS = feather.read_feather(f"{ABIDE_ASD_data_path}/ABIDE_ASD_FC1000_pyspi14_filtered_MixedSigmoid_counts.feather")

        # Concatenate all results
        all_results = pd.concat([UCLA_CNP_catch25_data_counts, ABIDE_ASD_catch25_data_counts, UCLA_CNP_pyspi14_data_counts, ABIDE_ASD_pyspi14_data_counts,
                                UCLA_CNP_catch25_data_z, ABIDE_ASD_catch25_data_z, UCLA_CNP_pyspi14_data_z, ABIDE_ASD_pyspi14_data_z,
                                UCLA_CNP_catch25_data_MS, ABIDE_ASD_catch25_data_MS, UCLA_CNP_pyspi14_data_MS, ABIDE_ASD_pyspi14_data_MS]).reset_index()
        all_results.to_feather(f"{final_data_path}/all_normalisations_counts.feather")