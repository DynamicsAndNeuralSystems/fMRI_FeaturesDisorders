import pandas as pd
import numpy as np
import os

# Amyloid centiloid data
data_path="/headnode1/abry4213/data/TS_feature_manuscript/"

disorder_dataset_df = pd.DataFrame({"Dataset_ID": ["UCLA_CNP", "UCLA_CNP", "UCLA_CNP", "ABIDE"], "Disorder": ["SCZ", "BP", "ADHD", "ASD"]})

# Iterate over rows of the dataframe
for index, row in disorder_dataset_df.iterrows():
    dataset_ID = row["Dataset_ID"]
    disorder = row["Disorder"]

    # Load sample metadata
    sample_metadata = pd.read_feather(f"{data_path}/input_data/{dataset_ID}_sample_metadata.feather")

    # Load catch25 results
    catch25_data = pd.read_feather(f"{data_path}/time_series_features/{dataset_ID}_catch25_filtered.feather")

    # Load sample metadata
    sample_metadata = pd.read_feather(f"{data_path}/input_data/{dataset_ID}_sample_metadata.feather")

    # Load catch25 results
    catch25_data = (pd.read_feather(f"{data_path}/time_series_features/{dataset_ID}_catch25_filtered.feather")
                    .merge(sample_metadata)
                    .query("Diagnosis in ['Control', @disorder]"))
        
    # First, iterate over each brain region
    for ROI in catch25_data.Brain_Region.unique().tolist():
        # Subset data to ROI
        ROI_data = catch25_data.query("Brain_Region == @ROI").drop(["Brain_Region", "Noise_Proc"], axis=1)

        # Replace spaces in ROI with underscores
        ROI = ROI.replace(" ", "_").replace(",_", "_")
                
        # Pivot from long to wide
        ROI_data_wide = ROI_data.pivot(index=['Sample_ID', 'Diagnosis'], columns='names', values='values')
        
        # Extract sample ID and diagnosis as lists
        index_data = ROI_data_wide.index.to_frame().reset_index(drop=True)
        class_labels = [int(i==disorder) for i in index_data["Diagnosis"].tolist()]
        sample_IDs = index_data["Sample_ID"].tolist()

        # Extract only the feature data
        regional_features_only = ROI_data_wide.reset_index(drop=True).to_numpy()

        # Save to a .npy file
        np.save(f"{data_path}/time_series_features/processed_numpy_files/{dataset_ID}_{disorder}_ROI_{ROI}.npy", regional_features_only)

        # Append text file with the models
        with open(f"{data_path}/time_series_features/processed_numpy_files/{dataset_ID}_{disorder}_models.txt", "a") as f:
            f.write(f"{dataset_ID}_{disorder}_ROI_{ROI}\n")

    # Then, iterate over each catch25 feature
    for catch25_feature in catch25_data.names.unique().tolist():
        # Subset data to ROI
        catch25_feature_data = catch25_data.query("names == @catch25_feature").drop(["names", "Noise_Proc"], axis=1)
                
        # Pivot from long to wide
        catch25_feature_data_wide = catch25_feature_data.pivot(index=['Sample_ID', 'Diagnosis'], columns='Brain_Region', values='values')
        
        # Extract sample ID and diagnosis as lists
        index_data = catch25_feature_data_wide.index.to_frame().reset_index(drop=True)
        class_labels = [int(i==disorder) for i in index_data["Diagnosis"].tolist()]
        sample_IDs = index_data["Sample_ID"].tolist()

        # Extract only the feature data
        catch25_features_only = catch25_feature_data_wide.reset_index(drop=True).to_numpy()

        # Save to a .npy file
        np.save(f"{data_path}/time_series_features/processed_numpy_files/{dataset_ID}_{disorder}_catch25_feature_{catch25_feature}.npy", catch25_features_only)

        # Append text file with the models
        with open(f"{data_path}/time_series_features/processed_numpy_files/{dataset_ID}_{disorder}_models.txt", "a") as f:
            f.write(f"{dataset_ID}_{disorder}_catch25_feature_{catch25_feature}\n")

    # Lastly, combine feature names and region names into one matrix
    catch25_data["Combo_Feature"] = catch25_data.Brain_Region + "_" + catch25_data.names
    combo_data = catch25_data.drop(["Brain_Region", "names"], axis=1)

    # Pivot from long to wide
    combo_data_wide = combo_data.pivot(index=["Sample_ID", "Diagnosis"],
                                    columns = "Combo_Feature",
                                    values = "values")

    # Extract only the combo feature data
    combo_features_only = combo_data_wide.reset_index(drop=True).to_numpy()

    # Save the combo feature data
    np.save(f"{data_path}/time_series_features/processed_numpy_files/{dataset_ID}_{disorder}_combo_catch25_features_all_regions.npy", combo_features_only)

    # Append text file with the models
    with open(f"{data_path}/time_series_features/processed_numpy_files/{dataset_ID}_{disorder}_models.txt", "a") as f:
        f.write(f"{dataset_ID}_{disorder}_combo_catch25_features_all_regions\n")