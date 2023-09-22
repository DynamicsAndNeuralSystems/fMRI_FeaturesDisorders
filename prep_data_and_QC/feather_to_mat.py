import pandas as pd
import scipy.io 
import sys
import os

input_file_base = sys.argv[1]
output_file_format = sys.argv[2]

if len(sys.argv) == 6:
    dataset_ID = sys.argv[3]
    noise_proc = sys.argv[4]
    subjects_to_keep_file = sys.argv[5]

# Check if output format is .feather or .mat
if output_file_format == "feather":
    if not os.path.isfile(input_file_base+"_filtered.feather"):
        loaded_data = scipy.io.loadmat(input_file_base+".mat")

        fALFF_data = pd.DataFrame.from_dict(loaded_data['fALFF_data'])
        
        # Unlist columns
        fALFF_data = fALFF_data.applymap(lambda x: x[0])
        
        # Unlist columns 2 and 3 again
        fALFF_data[2] = fALFF_data[2].apply(lambda x: x[0])
        
        # Rename columns
        fALFF_data.columns = ['Sample_ID', 'Brain_Region', 'fALFF']

        # Reshape from wide to long
        fALFF_data = pd.melt(fALFF_data, id_vars=['Sample_ID', 'Brain_Region'], var_name='names', value_name='values')
        
        # Add columns
        fALFF_data = (fALFF_data
                           .assign(Noise_Proc = noise_proc,
                                   method = "fALFF",
                                   Feature_Set = "fALFF"))
        
        # Write output to feather file
        fALFF_data.to_feather(input_file_base+".feather")

elif output_file_format == "mat":
    if not os.path.isfile(input_file_base+".mat"):
        loaded_data = pd.read_feather(input_file_base+".feather")
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
        scipy.io.savemat(input_file_base+".mat",
                        {'feather_to_mat_struct': subjects_and_regions.to_dict('list')})

        
