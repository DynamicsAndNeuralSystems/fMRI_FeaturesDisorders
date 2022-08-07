import pandas as pd
import dill
import argparse
import os

# Command-line arguments to parse
parser = argparse.ArgumentParser(description='Process inputs for pairwise data preparation.')
parser.add_argument('--github_dir', default="/project/hctsa/annie/github/", dest='github_dir')
parser.add_argument('--data_path', default="/project/hctsa/annie/data/UCLA_Schizophrenia/", dest='data_path')
parser.add_argument('--noise_procs', default=["AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER"], nargs='*', action='append', dest='noise_procs')
parser.add_argument('--parcellation_name', default="harvard_oxford_cort_prob_2mm", dest='parcellation_name', nargs='?')
parser.add_argument('--brain_region_lookup', default="Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv", dest='brain_region_lookup', nargs='?')
parser.add_argument('--dataset_ID', default="UCLA_Schizophrenia", dest='dataset_ID')
parser.add_argument('--overwrite', default=False, action="store_true")

# Parse arguments
args = parser.parse_args()
github_dir = args.github_dir
data_path = args.data_path
noise_procs = args.noise_procs
parcellation_name = args.parcellation_name
brain_region_lookup = args.brain_region_lookup
dataset_ID = args.dataset_ID
overwrite = args.overwrite

fmri_github_dir = github_dir + "fMRI_FeaturesDisorders/"

# github_dir = "/media/sf_Shared_Folder/github/"
# fmri_github_dir="/media/sf_Shared_Folder/github/fMRI_FeaturesDisorders/"
# data_path="/media/sf_Shared_Folder/PhD_work/data/ABIDE_ASD/"
# noise_procs=["FC1000"]
# dataset_ID="ABIDE_ASD"


def pkl_to_csv(pkl_file, output_csv):
    
    with open(pkl_file,'rb') as f:
        SPI_res = dill.load(f)
    
    # Iterate over each SPI
    SPI_res.columns = SPI_res.columns.to_flat_index()
      
    # Convert index to column
    SPI_res.reset_index(level=0, inplace=True)
      
    # Rename index as first brain region
    SPI_res = SPI_res.rename(columns={"index": "brain_region_1"})
      
    # Pivot data from wide to long
    SPI_res_long = pd.melt(SPI_res, id_vars="brain_region_1")
    SPI_res_long['SPI'], SPI_res_long['brain_region_2'] = SPI_res_long.variable.str
    
    # Remove variable column
    SPI_res_long = SPI_res_long.drop("variable", 1)
    
    # Write to CSV
    SPI_res_long.to_csv(output_csv, index=None)
    
for noise_proc in noise_procs:
    # Rename noise processing method to have underscores
    noise_label = noise_proc.replace("+", "_")
    
    # Get list of subjects
    input_data_path = data_path + "pydata/" + noise_label + "/"
    subject_list = [filename for filename in os.listdir(input_data_path) if 
                    os.path.isdir(os.path.join(input_data_path,filename))]
    
    # Iterate over each subject
    for subject in subject_list:
        subject_pkl = f"{input_data_path}/{subject}/calc.pkl"
        subject_csv = f"{input_data_path}/{subject}/calc.csv"
        
        # Write calc.table from pkl file to a CSV
        if not os.path.isFile(subject_csv) or overwrite:
            pkl_to_csv(pkl_file = subject_pkl,
                       output_csv = subject_csv)