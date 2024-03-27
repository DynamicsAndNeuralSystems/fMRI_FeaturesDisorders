import pandas as pd
import numpy as np
import os

# Define data path
data_path="/headnode1/abry4213/data/TS_feature_manuscript/"

# Load SPI directionality information
SPI_directionality_data = pd.read_csv("/headnode1/abry4213/github/fMRI_FeaturesDisorders/code/classification_analysis/SPI_Direction_Info.csv")
SPI_directionality_dict = dict(SPI_directionality_data.values)

# dataframe that maps dataset ID to disorder
disorder_dataset_df = pd.DataFrame({"Dataset_ID": ["UCLA_CNP", "UCLA_CNP", "UCLA_CNP", "ABIDE"], "Disorder": ["SCZ", "BP", "ADHD", "ASD"]})

# Iterate over rows of the dataframe
for index, row in disorder_dataset_df.iterrows():
    dataset_ID = row["Dataset_ID"]
    disorder = row["Disorder"]

    univariate_output_model_file = f"{data_path}/time_series_features/processed_numpy_files/{dataset_ID}_{disorder}_univariate_models.txt"
    pairwise_output_model_file = f"{data_path}/time_series_features/processed_numpy_files/{dataset_ID}_{disorder}_pairwise_models.txt"
    univariate_pairwise_combo_output_model_file = f"{data_path}/time_series_features/processed_numpy_files/{dataset_ID}_{disorder}_combined_univariate_pairwise_models.txt"

    open(univariate_output_model_file, 'w').close()
    open(pairwise_output_model_file, 'w').close()
    open(univariate_pairwise_combo_output_model_file, 'w').close()

    # Load sample metadata
    sample_metadata = pd.read_feather(f"{data_path}/input_data/{dataset_ID}_sample_metadata.feather")
    # Load catch25 results
    catch25_data = pd.read_feather(f"{data_path}/time_series_features/{dataset_ID}_catch25_filtered.feather")

    # Load sample metadata
    sample_metadata = pd.read_feather(f"{data_path}/input_data/{dataset_ID}_sample_metadata_filtered.feather")

    # Load catch25 results
    catch25_data = (pd.read_feather(f"{data_path}/time_series_features/{dataset_ID}_catch25_filtered.feather")
                    .merge(sample_metadata)
                    .query("Diagnosis in ['Control', @disorder] & Sample_ID in @sample_metadata.Sample_ID"))
        
    # First, iterate over each brain region
    for ROI in catch25_data.Brain_Region.unique().tolist():
        if not os.path.isfile(f"{data_path}/time_series_features/processed_numpy_files/{dataset_ID}_{disorder}_ROI_{ROI}.npy"):
            # Subset data to ROI
            ROI_data = catch25_data.query("Brain_Region == @ROI").drop(["Brain_Region", "Noise_Proc"], axis=1)

            # Replace spaces in ROI with underscores
            ROI = ROI.replace(" ", "_").replace(",_", "_")
                    
            # Pivot from long to wide
            ROI_data_wide = ROI_data.pivot(index=['Sample_ID', 'Diagnosis'], columns='names', values='values')

            # Extract only the feature data
            regional_features_only = ROI_data_wide.reset_index(drop=True).to_numpy()

            # Save to a .npy file
            np.save(f"{data_path}/time_series_features/processed_numpy_files/{dataset_ID}_{disorder}_ROI_{ROI}.npy", regional_features_only)

        # Append text file with the models
        with open(univariate_output_model_file, "a") as f:
            f.write(f"{dataset_ID}_{disorder}_ROI_{ROI}\n")

    # Then, iterate over each catch25 feature
    for catch25_feature in catch25_data.names.unique().tolist():
        if not os.path.isfile(f"{data_path}/time_series_features/processed_numpy_files/{dataset_ID}_{disorder}_catch25_feature_{catch25_feature}.npy"):
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
        with open(univariate_output_model_file, "a") as f:
            f.write(f"{dataset_ID}_{disorder}_catch25_feature_{catch25_feature}\n")

    # Lastly, combine feature names and region names into one matrix
    catch25_data["Combo_Feature"] = catch25_data.Brain_Region + "_" + catch25_data.names
    combo_data = catch25_data.drop(["Brain_Region", "names"], axis=1)

    # Pivot from long to wide
    combo_data_wide = combo_data.pivot(index=["Sample_ID", "Diagnosis"],
                                    columns = "Combo_Feature",
                                    values = "values")
    
    # Save class_labels and sample_IDs if they aren't already files
    if not os.path.isfile(f"{data_path}/input_data/{dataset_ID}_{disorder}_class_labels.npy"):
        index_data = combo_data_wide.index.to_frame().reset_index(drop=True)
        class_labels = np.array([int(i==disorder) for i in index_data["Diagnosis"].tolist()])
        sample_IDs = np.array(index_data["Sample_ID"].tolist())
        np.save(f"{data_path}/input_data/{dataset_ID}_{disorder}_class_labels.npy", class_labels)
        np.save(f"{data_path}/input_data/{dataset_ID}_{disorder}_sample_IDs.npy", sample_IDs)

    # Extract only the combo feature data
    combo_features_only = combo_data_wide.reset_index(drop=True).to_numpy()

    # Save the combo feature data
    np.save(f"{data_path}/time_series_features/processed_numpy_files/{dataset_ID}_{disorder}_combo_catch25_features_all_regions.npy", combo_features_only)

    # Append text file with the models
    with open(univariate_output_model_file, "a") as f:
        f.write(f"{dataset_ID}_{disorder}_combo_catch25_features_all_regions\n")

    #########################################################################################################
    # Load pyspi14 results
    pyspi14_data = (pd.read_feather(f"{data_path}/time_series_features/{dataset_ID}_pyspi14_filtered.feather")
                    .merge(sample_metadata)
                    .query("Diagnosis in ['Control', @disorder] & Sample_ID in @sample_metadata.Sample_ID"))
    
    # Iterate over each SPI
    for pyspi14_SPI in pyspi14_data.SPI.unique().tolist():

        # Save one numpy file for just the SPI
        if not os.path.isfile(f"{data_path}/time_series_features/processed_numpy_files/{dataset_ID}_{disorder}_pyspi14_SPI_{pyspi14_SPI}.npy"):
            # Subset data to ROI
            pyspi14_SPI_data = pyspi14_data.query("SPI == @pyspi14_SPI").drop(["SPI"], axis=1)

            # Find directionality of SPI
            SPI_directionality = SPI_directionality_dict[pyspi14_SPI]

            # Merge brain regions according to directionality
            if SPI_directionality == "Directed":
                pyspi14_SPI_data["region_pair"] = pyspi14_SPI_data.brain_region_from + "_" + pyspi14_SPI_data.brain_region_to
                pyspi14_SPI_data = pyspi14_SPI_data.drop(["brain_region_from", "brain_region_to"], axis=1)
            else:
                pyspi14_SPI_data_sorted = [sorted(pair) for pair in pyspi14_SPI_data[["brain_region_from", "brain_region_to"]].values.tolist()]
                pyspi14_SPI_data['region_pair'] = ['_'.join(string) for string in pyspi14_SPI_data_sorted]
                pyspi14_SPI_data = (pyspi14_SPI_data
                            .drop(["brain_region_from", "brain_region_to"], axis=1)
                            .drop_duplicates(ignore_index=True, subset=['Sample_ID', 'region_pair'])
                            )
                    
            # Pivot from long to wide
            pyspi14_SPI_data_wide = pyspi14_SPI_data.pivot(index=['Sample_ID', 'Diagnosis'], columns='region_pair', values='value')
            
            # Impute any NaN with column mean
            pyspi14_SPI_data_imputed = pyspi14_SPI_data_wide.fillna(pyspi14_SPI_data_wide.mean())

            # Extract only the feature data
            pyspi14_SPIs_only = pyspi14_SPI_data_imputed.reset_index(drop=True).to_numpy()

            # Save to a .npy file
            np.save(f"{data_path}/time_series_features/processed_numpy_files/{dataset_ID}_{disorder}_pyspi14_SPI_{pyspi14_SPI}.npy", pyspi14_SPIs_only)

        # Append text file with the models
        with open(pairwise_output_model_file, "a") as f:
            f.write(f"{dataset_ID}_{disorder}_pyspi14_SPI_{pyspi14_SPI}\n")

        # Another numpy file for the concatenated SPI data + all univariate data
        if not os.path.isfile(f"{data_path}/time_series_features/processed_numpy_files/{dataset_ID}_{disorder}_combined_univariate_catch25_and_pyspi14_SPI_{pyspi14_SPI}.npy"):
            pyspi14_SPI_and_univariate_data_combined = np.concatenate([pyspi14_SPIs_only, combo_features_only], axis=1)
            
            # Save to a .npy file
            np.save(f"{data_path}/time_series_features/processed_numpy_files/{dataset_ID}_{disorder}_combined_univariate_catch25_and_pyspi14_SPI_{pyspi14_SPI}.npy", pyspi14_SPI_and_univariate_data_combined)

        # Append text file with the models
        with open(univariate_pairwise_combo_output_model_file, "a") as f:
            f.write(f"{dataset_ID}_{disorder}_combined_univariate_catch25_and_pyspi14_SPI_{pyspi14_SPI}\n")
